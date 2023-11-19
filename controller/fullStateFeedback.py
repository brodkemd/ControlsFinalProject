from controller.stateSpace import LandingStateSpace, FlipStateSpaceCP, DescentStateSpaceCP
import numpy as np
import control, warnings

class BaseFullStateFeedBack:
    def __init__(self) -> None:
        pass

    def generateGains(self, A:np.ndarray, B:np.ndarray, C:np.ndarray, poles:list):
        CC = control.ctrb(A, B)
        CC_rank = np.linalg.matrix_rank(CC)
        if CC_rank != A.shape[0]:
            raise ValueError(f"State space is not controllable\n  Controlability rank: {CC_rank}\n  # A rows: {A.shape[0]}")

        warnings.filterwarnings("error")
        K = control.place(A, B, poles)
        warnings.filterwarnings("default")

        k_r = -1/(C @ (np.linalg.inv(A - B @ K) @ B))

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

        self.K, self.k_r = self.generateGains(self.A, self.B, self.C, poles)

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

        self.K, self.k_r = self.generateGains(self.A, self.B, self.C, poles)

    def update(self, x_r, x):
        pass


class Landing(LandingStateSpace, BaseFullStateFeedBack):
    def __init__(self) -> None:
        super().__init__()

        self.zetas    =  0.707*np.zeros(len(self.A)//2)
        self.omega_ns = [3.0, 1.5, 2.5, 1.25]

        poles = []
        polys = []
        for i in range(len(self.A)//2):
            zeta    = self.zetas[i]
            omega_n = self.omega_ns[i]
            polys.append([1, 2*zeta*omega_n, omega_n**2])
            poles.append(-zeta*omega_n + omega_n*np.sqrt(1 - zeta**2)*1j)
            poles.append(-zeta*omega_n - omega_n*np.sqrt(1 - zeta**2)*1j)

        #print(poles)
        #poles = np.roots(poly)
        print(poles)

        self.K, self.k_r = self.generateGains(self.A, self.B, self.C, poles)

    def update(self, x_r, x):
        pass