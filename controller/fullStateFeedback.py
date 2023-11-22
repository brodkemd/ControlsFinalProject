from controller.stateSpace import LandingStateSpace, FlipStateSpaceCP, DescentStateSpaceCP
import numpy as np
import control, os, platform
import scipy.io

class BaseFullStateFeedBack:
    def __init__(self) -> None:
        pass

    def generateGains(self, A:np.ndarray, B:np.ndarray, C_r:np.ndarray, poles:list[np.complex64], name:str, compute_gains=True):
        cwd  = os.path.dirname(__file__)
        file_stem = os.path.join(cwd, name)
        
        if "linux" in platform.system().lower() and compute_gains:
            print("    Checking Controllability of:", name, end = " ")
            CC = control.ctrb(A, B)
            CC_rank = np.linalg.matrix_rank(CC)
            if CC_rank != A.shape[0]:
                raise ValueError(f"State space is not controllable\n  Controlability rank: {CC_rank}\n  # A rows: {A.shape[0]}")
            print("done")

            print("    Computing Gains:", end=" ")
            A_str = []
            for row in A:
                A_str_row = []
                for val in row:
                    A_str_row.append("{:.16f}".format(val))
                A_str.append(A_str_row)

            B_str = []
            for row in B:
                B_str_row = []
                for val in row:
                    B_str_row.append("{:.16f}".format(val))
                B_str.append(B_str_row)

            pole_str = []
            for pole in poles:
                if pole.imag < 0.0:
                    pole_str.append("{:.16f} - {:.16f}j".format(pole.real, abs(pole.imag)))
                else:
                    pole_str.append("{:.16f} + {:.16f}j".format(pole.real, pole.imag))

            for i in range(len(A_str)):
                A_str[i] = ",".join(A_str[i])
            for i in range(len(B_str)):
                B_str[i] = ",".join(B_str[i])

            A_str    = ";\n".join(A_str)
            B_str    = ";\n".join(B_str)
            pole_str = ",\n".join(pole_str)

            A_str    = f"A = [{A_str}];"
            B_str    = f"B = [{B_str}];"
            pole_str = f"poles = [{pole_str}];"
            
            with open(file_stem+".m", 'w') as f:
                f.write(A_str   +"\n")
                f.write(B_str   +"\n")
                f.write(pole_str+"\n")
                f.write("[K, prec] = place(A, B, poles);\n")
                f.write(f"save('{file_stem}.mat', 'K', 'prec');\n")
                f.write(f"exit;\n")

            os.system(f'/usr/local/MATLAB/R2021b/bin/matlab -nodisplay -nosplash -nodesktop -batch "{name}" -sd "{cwd}" > /dev/null')
            print("done")
        else:
            if compute_gains:
                print("    Can not compute gains without linux")

        print(f"    Loading gain matrix:", end=" ")
        mat_data = scipy.io.loadmat(f"{file_stem}.mat")
        print("Done")
        K    = mat_data["K"]
        prec = mat_data["prec"].item(0)
        print(f"       Gains have a decimal precision of: {prec}")
        # print(np.linalg.inv(A - B @ K) @ B)
        K_r = -1*np.linalg.inv(C_r @ (np.linalg.inv(A - B @ K) @ B))
        
        return K, K_r


class FullStateFeedBack:
    def __init__(self, compute_gains=False) -> None:
        print("\nFullStateFeedBack:")
        self.landing    = Landing(compute_gains)
        # self.descentCP  = DescentCP(compute_gains)
        # self.flipCP     = FlipCP(compute_gains)

        self.state = 0 # 0 = landing, 1 = flip, 2 = descent

    def update(self, x, x_r):
        u = self.landing.update(x, x_r)
        return u[0:3], u[3:6], u[6:9]


class DescentCP(DescentStateSpaceCP, BaseFullStateFeedBack):
    def __init__(self, compute_gains=False) -> None:
        print("  DescentCP:")
        super().__init__()

        self.zetas    =  0.707*np.ones(len(self.A)//2)
        self.omega_ns = 0.01*np.array([3.0, 1.5, 2.5, 1.25])

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

        self.K, self.K_r = self.generateGains(self.A, self.B, self.C_r, poles, "descentCP", compute_gains=compute_gains)

    def update(self, x, y_r):
        return self.u_e - self.K(x - self.x_e) + self.K_r*(y_r - self.y_re)


class FlipCP(FlipStateSpaceCP, BaseFullStateFeedBack):
    def __init__(self, compute_gains=False) -> None:
        print("  FlipCP:")
        super().__init__()

        self.zetas    =  0.707*np.ones(len(self.A)//2)
        self.omega_ns = 0.01*np.array([3.0, 1.5, 2.5, 1.25])
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
        poles.append(-1.0)

        self.K, self.K_r = self.generateGains(self.A, self.B, self.C_r, poles, "flipCP", compute_gains=compute_gains)

    def update(self, x, y_r):
        return self.u_e - self.K(x - self.x_e) + self.K_r*(y_r - self.y_re)


class Landing(LandingStateSpace, BaseFullStateFeedBack):
    def __init__(self, compute_gains=False) -> None:
        print("  Landing:")
        super().__init__()

        self.zetas    =  0.707*np.ones(len(self.A)//2)
        self.omega_ns = 0.01*np.array([3.0, 1.5, 2.5, 1.25, 1.75])
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

        self.K, self.K_r = self.generateGains(self.A, self.B, self.C_r, poles, "landing", compute_gains=compute_gains)


    def update(self, _x, _x_r):
        # x_vars = "p_n,p_e,p_d,u,v,w,e_0,e_3,q,r".split(",") # order matters
        _x_tilde   = _x - self.x_e
        x_tilde    = np.append(  _x_tilde[:7], np.array([  _x_tilde.item(9),   _x_tilde.item(11),   _x_tilde.item(12)]))
        y_r_tilde  = _x_r[:3] - self.y_re
        u = np.zeros(9)

        u[0:3] = self.K @ x_tilde + self.K_r @ y_r_tilde
        return self.u_e - u