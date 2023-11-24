from controller.stateSpace import LandingStateSpace, FlipStateSpaceCP, DescentStateSpaceCP, BaseStateSpace
from tools.rotations import Euler2Quaternion, Quaternion2Euler
from parameters.baseClass import Base
from tools.toString import arrToStr, numToStr
import numpy as np
import control, os, platform, warnings, json, subprocess
import scipy.io

warnings.filterwarnings("error")

class BaseFullStateFeedBack:
    def __init__(self) -> None:
        super().__init__()
        

    def generateGains(self, A:np.ndarray, B:np.ndarray, C_r:np.ndarray, poles:list[np.complex64], name:str, compute_gains=True, CC=None):
        cwd  = os.path.dirname(__file__)
        data = {}
        with open(os.path.join(os.path.dirname(cwd), "data", "user_config.json")) as f:
            data = json.loads(f.read())
        self.matlab = data["matlab"]
        file_stem = os.path.join(cwd, name)
        
        if len(poles) != len(A): raise ValueError("number of poles must be the same as the number of rows of A")

        if compute_gains:
            print("    Checking Controllability of:", name)
            if CC is None:
                print("      Computing Controllability Matrix")
                CC = control.ctrb(A, B)
            CC_rank = np.linalg.matrix_rank(CC)
            if CC_rank != A.shape[0]:
                print()
                raise ValueError(f"State space is not controllable\n  Controlability rank: {CC_rank}\n  # A rows: {A.shape[0]}")
            print("        done")

            print("    Generating Matlab File: ", end="")
            A_str    = arrToStr(    A, joins=[";\n     "], starts=["[", ""], ends=["]", ""])
            B_str    = arrToStr(    B, joins=[";\n     "], starts=["[", ""], ends=["]", ""])
            pole_str = arrToStr(poles, joins=[";\n         "], starts=["[", ""], ends=["]", ""])
            
            with open(file_stem+".m", 'w') as f:
                f.write(f"A = {A_str};\n")
                f.write(f"B = {B_str};\n")
                f.write(f"poles = {pole_str};\n")
                f.write("[K, prec] = place(A, B, poles);\n")
                f.write(f"save('{file_stem}.mat', 'K', 'prec');\n")
                f.write(f"exit;\n")
            print("done")

            print("    Computing Gains:", end=" ")
            if "linux" in platform.system().lower() and compute_gains:
                os.system(f'{self.matlab} -nodisplay -nosplash -nodesktop -batch "{name}" -sd "{cwd}" > /dev/null')
            else:
                subprocess.call(f'& "{self.matlab}" -nodisplay -nosplash -nodesktop -batch "{name}" -sd "{cwd}"', shell=True)
            print("done")
        else:
            if compute_gains:
                print("    Can not compute gains without linux (do not know matlab is installed)")

        print(f"    Loading gain matrix:", end=" ")
        mat_data = scipy.io.loadmat(f"{file_stem}.mat")
        print("Done")
        K    = mat_data["K"]
        prec = mat_data["prec"].item(0)
        print(f"       Gains have a decimal precision of: {prec}")

        try:
            K_r = -1*np.linalg.inv(C_r @ (np.linalg.inv(A - B @ K) @ B))
        except Exception as e:
            print("    Can not compute K_r, returning zeros")
            K_r = np.zeros((C_r.shape[0], C_r.shape[0]))

        return K, K_r
    
    def update(self, _x, _x_r):
        # x_vars = "p_n,p_e,p_d,u,v,w,e_0,e_3,q,r".split(",") # order matters
        _x_tilde   = _x - self.x_e
        x_r        = self.global_x_r_to_local_x_r @ _x_r
        x_tilde    = self.global_x_to_local_x @ _x_tilde
        y_r_tilde  = x_r - self.y_re
        
        u = self.local_u_to_global_u @ (self.K @ x_tilde + self.K_r @ y_r_tilde)
        return self.u_e - u
    
    def eigenVectorAnalysis(self):
        poles, vecs = np.linalg.eig(self.A - self.B@self.K)
        print("  Eigenvector Analysis:")
        for i in range(len(poles)):
            if np.iscomplex(poles[i]):
                print("    pole:", numToStr(poles[i], ":.6f"))
            print("    vector:", arrToStr(vecs[:,i], joins="\n"+13*" ", num_format=":.6f"))


class FullStateFeedBack(Base):
    def __init__(self, compute_gains=False) -> None:
        super().__init__()
        print("\nFullStateFeedBack:")
        self.landing    = Landing(compute_gains=compute_gains)
        self.descentCP  = DescentCP(compute_gains=False)
        self.flipCP     = FlipCP(compute_gains=False)
        self.state = 2 # 0 = landing, 1 = flip, 2 = descent

        phi_e = 0
        theta_e = 0
        psi_e = 0
        self.x_r = np.array([
            0, # p_n
            0, # p_e
            0, # p_d
            0, # u
            0, # v
            0, # w
            Euler2Quaternion(phi_e, theta_e, psi_e).item(0), # e_0
            Euler2Quaternion(phi_e, theta_e, psi_e).item(1), # e_1
            Euler2Quaternion(phi_e, theta_e, psi_e).item(2), # e_2
            Euler2Quaternion(phi_e, theta_e, psi_e).item(3), # e_3
            0, # p
            0, # q
            0 # r
        ], dtype=float)

    def update(self, x):
        last_state = self.state
        if self.state == 2:
            if x.item(2) > -3000: # determining from position
                self.state = 1
                phi_e = 0
                theta_e = np.deg2rad(90)
                psi_e = 0
                self.x_r[6] = Euler2Quaternion(phi_e, theta_e, psi_e).item(0) # e_0
                self.x_r[7] = Euler2Quaternion(phi_e, theta_e, psi_e).item(1) # e_1
                self.x_r[8] = Euler2Quaternion(phi_e, theta_e, psi_e).item(2) # e_2
                self.x_r[9] = Euler2Quaternion(phi_e, theta_e, psi_e).item(3) # e_3

            else:
                phi_e = 0
                theta_e = 0
                psi_e = 0
                self.x_r[6] = Euler2Quaternion(phi_e, theta_e, psi_e).item(0) # e_0
                self.x_r[7] = Euler2Quaternion(phi_e, theta_e, psi_e).item(1) # e_1
                self.x_r[8] = Euler2Quaternion(phi_e, theta_e, psi_e).item(2) # e_2
                self.x_r[9] = Euler2Quaternion(phi_e, theta_e, psi_e).item(3) # e_3
        
        elif self.state == 1:
            phi_e = 0
            theta_e = 90
            psi_e = 0

            e = x[6:10]
            theta = np.rad2deg(Quaternion2Euler(e)[1])
            if abs(theta_e - theta) < 10:
                theta_e = np.deg2rad(theta_e)
                self.state = 0
                self.x_r[6] = Euler2Quaternion(phi_e, theta_e, psi_e).item(0) # e_0
                self.x_r[7] = Euler2Quaternion(phi_e, theta_e, psi_e).item(1) # e_1
                self.x_r[8] = Euler2Quaternion(phi_e, theta_e, psi_e).item(2) # e_2
                self.x_r[9] = Euler2Quaternion(phi_e, theta_e, psi_e).item(3) # e_3
    
        if self.state != last_state:
            print("\nChanged state from", last_state, "to", self.state, end="\n\n")

        if self.state == 0 or self.state == 1:
            u = self.landing.update(x, self.x_r)
        # elif self.state == 1:
        #     u = self.flipCP.update(x, self.x_r)
        elif self.state == 2:
            u = self.descentCP.update(x, self.x_r)

        F_E                   = u[0:3]
        F_cp_port_canard      = u[3:6]
        F_cp_starboard_canard = u[6:9]

        F_cp_port_fin         = u[9:12]
        F_cp_starboard_fin    = u[12:15]

        F_cp_canard           = F_cp_port_canard + F_cp_starboard_canard
        F_cp_fin              = F_cp_port_fin + F_cp_starboard_fin

        F_cp = F_cp_fin + F_cp_canard

        tau_cp_port_canard      = np.cross(     self.r_cp_port_canard,      F_cp_port_canard)
        tau_cp_starboard_canard = np.cross(self.r_cp_starboard_canard, F_cp_starboard_canard)
        tau_cp_port_fin         = np.cross(        self.r_cp_port_fin,         F_cp_port_fin)
        tau_cp_starboard_fin    = np.cross(   self.r_cp_starboard_fin,    F_cp_starboard_fin)

        tau = tau_cp_port_canard + tau_cp_port_fin + tau_cp_starboard_canard + tau_cp_starboard_fin

        return F_E, F_cp, tau, self.x_r


class DescentCP(DescentStateSpaceCP, BaseFullStateFeedBack):
    def __init__(self, compute_gains=False) -> None:
        print("  DescentCP:")
        super().__init__()

        self.zetas    = 0.707*np.ones(len(self.A)//2)
        self.omega_ns = 0.1*np.array([3.0, 1.5, 2.5, 1.25, 1.0])

        self.t_r = max(2.2/self.omega_ns)
        self.t_s = max(4/(self.zetas[i]*self.omega_ns[i]) for i in range(len(self.zetas)))
        print("    Max rise time:   ", self.t_r)
        print("    Max setting time:", self.t_s)


        poles = []
        for i in range(len(self.A)//2):
            zeta    = self.zetas[i]
            omega_n = self.omega_ns[i]
            poles.append(-zeta*omega_n + omega_n*np.sqrt(1 - zeta**2)*1j)
            poles.append(-zeta*omega_n - omega_n*np.sqrt(1 - zeta**2)*1j)

        poles.append(-0.1)

        self.K, self.K_r = self.generateGains(self.A, self.B, self.C_r, poles, "descentCP", compute_gains=compute_gains, CC=self.CC)
        self.eigenVectorAnalysis()


class FlipCP(FlipStateSpaceCP, BaseFullStateFeedBack):
    def __init__(self, compute_gains=False) -> None:
        print("  FlipCP:")
        super().__init__()

        self.zetas    =  0.707*np.ones(len(self.A)//2)
        self.omega_ns = 0.1*np.array([3.0, 1.5, 2.5, 1.25, 1, 1.15])
        self.t_r = max(2.2/self.omega_ns)
        self.t_s = max(4/(self.zetas[i]*self.omega_ns[i]) for i in range(len(self.zetas)))
        print("    Max rise time:   ", self.t_r)
        print("    Max setting time:", self.t_s)

        poles = []
        for i in range(len(self.A)//2):
            zeta    = self.zetas[i]
            omega_n = self.omega_ns[i]
            poles.append(-zeta*omega_n + omega_n*np.sqrt(1 - zeta**2)*1j)
            poles.append(-zeta*omega_n - omega_n*np.sqrt(1 - zeta**2)*1j)

        self.K, self.K_r = self.generateGains(self.A, self.B, self.C_r, poles, "flipCP", compute_gains=compute_gains, CC=self.CC)
        self.eigenVectorAnalysis()


class Landing(LandingStateSpace, BaseFullStateFeedBack):
    def __init__(self, compute_gains=False) -> None:
        print("  Landing:")
        super().__init__()

        self.zetas    =  0.707*np.ones(len(self.A)//2)
        self.omega_ns = 0.05*np.array([3.0, 1.5, 2.5, 1.25, 1.75])
        self.t_r = max(2.2/self.omega_ns)
        self.t_s = max(4/(self.zetas[i]*self.omega_ns[i]) for i in range(len(self.zetas)))
        print("    Max rise time:   ", self.t_r)
        print("    Max setting time:", self.t_s)

        poles = []
        for i in range(len(self.A)//2):
            zeta    = self.zetas[i]
            omega_n = self.omega_ns[i]
            poles.append(-zeta*omega_n + omega_n*np.sqrt(1 - zeta**2)*1j)
            poles.append(-zeta*omega_n - omega_n*np.sqrt(1 - zeta**2)*1j)

        self.K, self.K_r = self.generateGains(self.A, self.B, self.C_r, poles, "landing", compute_gains=compute_gains, CC=self.CC)
        self.eigenVectorAnalysis()