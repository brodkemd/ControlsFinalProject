from controller.stateSpace import LandingStateSpace, FlipStateSpaceCP, DescentStateSpaceCP
import numpy as np
import control, os, platform
import scipy.io

class BaseFullStateFeedBack:
    def __init__(self) -> None:
        pass

    def generateGains(self, A:np.ndarray, B:np.ndarray, C:np.ndarray, poles:list[np.complex64], name:str, compute_gains=True):
        cwd  = os.path.dirname(__file__)
        file_stem = os.path.join(cwd, name)

        CC = control.ctrb(A, B)
        CC_rank = np.linalg.matrix_rank(CC)
        if CC_rank != A.shape[0]:
            raise ValueError(f"State space is not controllable\n  Controlability rank: {CC_rank}\n  # A rows: {A.shape[0]}")
        
        if "linux" in platform.system().lower() and compute_gains:
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
            
            print("Computing Gains:", end=" ")
            os.system(f'/usr/local/MATLAB/R2021b/bin/matlab -nodisplay -nosplash -nodesktop -batch "{name}" -sd "{cwd}" > /dev/null')
            print("done")
        else:
            if compute_gains:
                print("Can not compute gains without linux")

        print("Loading gain matrix:", end=" ")
        mat_data = scipy.io.loadmat(f"{file_stem}.mat")
        print("Done")
        K   = mat_data["K"]
        prec = mat_data["prec"].item(0)
        print(f"  Gains have a decimal precision of: {prec}")
        

        k_r = 0 # -1/(C @ (np.linalg.inv(A - B @ K) @ B))

        return K, k_r


class FullStateFeedBack:
    def __init__(self) -> None:
        self.landing = Landing()

    def update(self, x_r, x):
        pass


class Descent(DescentStateSpaceCP, BaseFullStateFeedBack):
    def __init__(self) -> None:
        super().__init__()

        self.zetas    =  0.707*np.ones(len(self.A)//2)
        self.omega_ns = [3.0, 1.5, 2.5, 1.25]

        poles = []
        for i in range(len(self.A)//2):
            zeta    = self.zetas[i]
            omega_n = self.omega_ns[i]
            poles.append(-zeta*omega_n + omega_n*np.sqrt(1 - zeta**2)*1j)
            poles.append(-zeta*omega_n - omega_n*np.sqrt(1 - zeta**2)*1j)

        self.K, self.k_r = self.generateGains(self.A, self.B, self.C, poles, "descent", compute_gains=False)

    def update(self, x_r, x):
        pass


class Flip(FlipStateSpaceCP, BaseFullStateFeedBack):
    def __init__(self) -> None:
        super().__init__()

        self.zetas    =  0.707*np.ones(len(self.A)//2)
        self.omega_ns = [3.0, 1.5, 2.5, 1.25]

        poles = []
        for i in range(len(self.A)//2):
            zeta    = self.zetas[i]
            omega_n = self.omega_ns[i]
            poles.append(-zeta*omega_n + omega_n*np.sqrt(1 - zeta**2)*1j)
            poles.append(-zeta*omega_n - omega_n*np.sqrt(1 - zeta**2)*1j)
        poles.append(-1.0)

        self.K, self.k_r = self.generateGains(self.A, self.B, self.C, poles, "flip", compute_gains=False)

    def update(self, x_r, x):
        pass


class Landing(LandingStateSpace, BaseFullStateFeedBack):
    def __init__(self) -> None:
        super().__init__()

        self.zetas    =  0.707*np.ones(len(self.A)//2)
        self.omega_ns = [3.0, 1.5, 2.5, 1.25, 1.75]

        poles = []
        for i in range(len(self.A)//2):
            zeta    = self.zetas[i]
            omega_n = self.omega_ns[i]
            poles.append(-zeta*omega_n + omega_n*np.sqrt(1 - zeta**2)*1j)
            poles.append(-zeta*omega_n - omega_n*np.sqrt(1 - zeta**2)*1j)

        self.K, self.k_r = self.generateGains(self.A, self.B, self.C, poles, "landing", compute_gains=False)

    def update(self, x_r, x):
        pass