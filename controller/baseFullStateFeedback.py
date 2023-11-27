from parameters import config
from parameters import simulation_parameters as SIM
from tools.toString import arrToStr, numToStr

import numpy as np
import control, os, platform, subprocess, datetime
import scipy.io

class BaseFullStateFeedBack:
    name, omega_ns, zetas = None, None, None
    has_integrator = False
    Ts = SIM.ts_simulation
    
    integrator = np.zeros(3) # integrator
    error_d1 = np.zeros(3) # error signal delayed by 1 sample

    def __init__(self, name:str) -> None:
        self.name = name
        self.cwd  = os.path.dirname(__file__)
        self.file_stem = os.path.join(self.cwd, name)
        print(2*" "+self.name+":")

    def generateGains(self, compute_gains=True, add_integrator=False):
        self.t_r = max(2.2/self.omega_ns)
        self.t_s = max(4/(self.zetas[i]*self.omega_ns[i]) for i in range(len(self.zetas)))
        print("    Max rise time:   ", self.t_r)
        print("    Max setting time:", self.t_s)

        if compute_gains:
            if add_integrator:
                print("    Generating Integrator State Space")
                A = np.vstack([
                    np.hstack([ self.A, np.zeros((self.A.shape[0], self.C.shape[0]))]),
                    np.hstack([-self.C, np.zeros((self.C.shape[0], self.C.shape[0]))])
                ])

                B = np.vstack(
                    [self.B, np.zeros((self.C.shape[0], self.B.shape[1]))]
                )
                print("    Checking Controllability with Integrator of:", self.name)
                print("      Computing Controllability Matrix")
                CC = control.ctrb(A, B)
                CC_rank = np.linalg.matrix_rank(CC)
                if CC_rank != A.shape[0]:
                    raise ValueError(f"State space is not controllable\n  Controlability rank: {CC_rank}\n  # A rows: {A.shape[0]}")
                print("        done")

            else:
                A = self.A.copy()
                B = self.B.copy()
                print("    Checking Controllability of:", self.name)
                if self.CC is None:
                    print("      Computing Controllability Matrix")
                    self.CC = control.ctrb(self.A, self.B)
                CC_rank = np.linalg.matrix_rank(self.CC)
                if CC_rank != self.A.shape[0]:
                    raise ValueError(f"State space is not controllable\n  Controlability rank: {CC_rank}\n  # A rows: {self.A.shape[0]}")
                print("        done")

            if len(self.poles) != len(A): raise ValueError("number of poles must be the same as the number of rows of A, # rows of A: " + str(len(A)) + "  # poles provided: " + str(len(self.poles)))
            print("    Generating Matlab File: ", end="")
            A_str    = arrToStr(    A, joins=[";\n     "],     starts=["[", ""], ends=["]", ""])
            B_str    = arrToStr(    B, joins=[";\n     "],     starts=["[", ""], ends=["]", ""])
            pole_str = arrToStr(self.poles, joins=[";\n         "], starts=["[", ""], ends=["]", ""])
            
            with open(self.file_stem+".m", 'w') as f:
                f.write(f"A = {A_str};\n")
                f.write(f"B = {B_str};\n")
                f.write(f"poles = {pole_str};\n")
                f.write("[K, prec] = place(A, B, poles);\n")
                f.write(f"save('{self.file_stem}.mat', 'K', 'prec');\n")
                f.write(f"exit;\n")
            print("done")

            print("    Computing Gains:", end=" ")
            mat_file = self.file_stem+".mat"
            if os.path.exists(mat_file):
                # Get the last modification time in seconds since the epoch
                modification_time = os.path.getmtime(mat_file)
                t = datetime.datetime.fromtimestamp(modification_time)
            else:
                t = datetime.datetime.fromtimestamp(0)

            if "linux" in platform.system().lower() and compute_gains:
                os.system(f'{config.matlab} -nodisplay -nosplash -nodesktop -batch "{self.name}" -sd "{self.cwd}" > /dev/null')
            else:
                subprocess.call(f'"{config.matlab}" -nosplash -nodesktop -batch "{self.name}" -sd "{self.cwd}"', shell=True)
            
            if os.path.exists(mat_file):
                # Get the last modification time in seconds since the epoch
                modification_time = os.path.getmtime(mat_file)
                t_new = datetime.datetime.fromtimestamp(modification_time)
                if t_new == t:
                    raise Exception("Did not generate gains, Matlab did not run properly")
            else:
                raise Exception("Did not generate gains, Matlab did not run properly")
            print("done")

        print(f"    Loading gain matrix:", end=" ")
        mat_data = scipy.io.loadmat(f"{self.file_stem}.mat")
        print("Done")
        K    = mat_data["K"]
        prec = mat_data["prec"].item(0)
        print(f"       Gains have a decimal precision of: {prec}")


        if add_integrator:
            print("    Separating Integrator and PD gains")
            self.K = K[:,:len(self.A)]
            self.k_I = K[:,len(self.A):]
            self.integrator = np.zeros(self.C.shape[0])
            self.error_d1 = np.zeros(self.C.shape[0])
        else:
            self.K = K
        self.has_integrator = add_integrator
        
        # try:
        #     self.K_r = -1*np.linalg.inv(self.C_r @ (np.linalg.inv(self.A - self.B @ self.K) @ self.B))
        # except Exception as e:
        #     print("    Can not compute K_r, returning zeros")
        #     self.K_r = np.zeros((self.C_r.shape[0], self.C_r.shape[0]))
    
    def computePoles(self):
        self.poles = []
        for i in range(len(self.A)//2):
            zeta    = self.zetas[i]
            omega_n = self.omega_ns[i]
            self.poles.append(-zeta*omega_n + omega_n*np.sqrt(1 - zeta**2)*1j)
            self.poles.append(-zeta*omega_n - omega_n*np.sqrt(1 - zeta**2)*1j)

    def update(self, _x, _x_r):
        # x_vars = "p_n,p_e,p_d,u,v,w,e_0,e_3,q,r".split(",") # order matters
        error   = self.global_x_r_to_local_x_r @ (_x -     _x_r)
        x_tilde = self.global_x_to_local_x     @ (_x - self.x_e)

        if self.has_integrator:
            self.integrateError(error)
            u = self.local_u_to_global_u @ (-1*self.K @ x_tilde - self.k_I @ self.integrator)
        else:
            u = self.local_u_to_global_u @ (-1*self.K @ x_tilde)
        return self.u_e + u
    
    def eigenVectorAnalysis(self):
        poles, vecs = np.linalg.eig(self.A - self.B@self.K)
        print("    Eigenvector Analysis:")
        for i in range(len(poles)):
            if np.iscomplex(poles[i]):
                print("      pole:", numToStr(poles[i], ":.6f"))
            print("      vector:", arrToStr(vecs[:,i], joins="\n"+15*" ", num_format=":.6f", shift_for_minus=True, pad=" "))

    def integrateError(self, error):
        self.integrator = self.integrator + (self.Ts/2.0)*(error + self.error_d1)
        self.error_d1 = error