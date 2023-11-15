import numpy as np
import parameters.body_parameters as BODY


class ForcesMoments:
    def __init__(self) -> None:
        # Surface fit coefficients for the force and moments
        self.Cdf = np.array[0.03605, -0.01182, -0.0001243, 0.01985, -0.0001574, 7.297E-5, 0.001095, -8.369E-6, -6.47E-7]
        self.Clf = np.array[0.01872, -0.07913, -0.002173, 0.06861, 0.004314, 0.0001143, -0.002112, -2.849E-5, -1.013E-6]
        self.Cmf = np.array[-0.02598, 0.1033, 0.0002977, -0.08253, -0.003769, -0.0002425, -0.001543, 5.302E-5, 2.275E-6]
       
        pass

    def update(self, M, AOA, B, delSurf, delT, delE):
        '''
        Inputs:
        M = Current Mach
        AOA = Current AOA [deg]
        B = Current sideslip, beta [deg]
        delSurf = Array of control surface deflections [delPortCanard, delStarCanard, delPortFin, delStarFin]
        delT = Throttle percentage of full throttle (can come from controller) [0,1]
        delE = Array of engine deflections [deflection from x-axis, rotation]

        '''
        # Defines the coefficients based on the current Mach and AOA
        Cd = 
        Cl =
        Cm =

        # Calculates the forces for the body and the control surfaces
        # body forces
        Fxb = 
        Fyb = 
        Fzb =

        #surface forces 
        Fxs = 
        Fsy = 
        Fsz = 

        #Engine forces
        FEx = 
        FEy = 
        FEz = 

        return 


class CPFromAerodynamics:
    def __init__(self) -> None:
        #Need locations from the c.g. of the canards and fins
        
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