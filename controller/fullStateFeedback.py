from controller.stateSpace import LandingStateSpace, FlipStateSpaceCP, DescentStateSpaceCP
import numpy as np
import control


class FullStateFeedBack:
    def __init__(self) -> None:
        self.landing = Landing()

    def update(self, x_r, x):
        pass


class Landing(LandingStateSpace):
    def __init__(self) -> None:
        super().__init__()

        self.zetas    = 0.707*np.ones(len(self.A)//2)
        self.omega_ns = [-3.0, -1.5, -2.5, -1.25]

        self.CC = control.ctrb(self.A, self.B)
        print(np.linalg.matrix_rank(self.CC), self.A.shape[0])

        poles = []
        for i in range(len(self.A)//2):
            zeta    = self.zetas[i]
            omega_n = self.omega_ns[i]
            print(1 - zeta**2)
            poles.append(-zeta*omega_n + omega_n*np.sqrt(1 - zeta**2)*1j)
            poles.append(-zeta*omega_n - omega_n*np.sqrt(1 - zeta**2)*1j)
        
        print(poles)

        self.K = control.place(self.A, self.B, [1, 1, 2, 4, 6, 7, 4, 6])
        print(self.K)
        self.k_r = -1/(self.C @ (np.linalg.inv(self.A - self.B @ self.K) @ self.B))

    def update(self, x_r, x):
        pass