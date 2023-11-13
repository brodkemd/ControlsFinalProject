import numpy as np
import parameters.body_parameters as BODY


class ForcesMoments:
    def __init__(self) -> None:
        pass

    def update(self):
        pass


class CPFromAerodynamics:
    def __init__(self) -> None:
        pass

    def update(self):
        pass


class ForcesMomentsFromCP:
    def __init__(self) -> None:
        self.r_E = BODY.r_E
        self.m   = BODY.mass
        self.g   = BODY.gravity

    def update(self, _state, F_E, r_cp, F_cp):
        e_0     = _state.item(6)
        e_1     = _state.item(7)
        e_2     = _state.item(8)
        e_3     = _state.item(9)

        F_g = self.m*self.g*np.array([
            2*(e_1*e_3 - e_2*e_0),
            2*(e_2*e_3 +  e_1*e_0),
            e_3**2 + e_0**2 - e_1**2 - e_2**2
        ])

        tau = np.cross(self.r_E, F_E) + np.cross(r_cp, F_cp)
        F   = F_E + F_g + F_cp

        return np.hstack((F, tau))