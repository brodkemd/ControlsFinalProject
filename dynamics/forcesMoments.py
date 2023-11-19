import numpy as np
import parameters.body_parameters as BODY


class ForcesMoments:
    def __init__(self) -> None:
        #Need locations from the c.g. of the canards and fins
        self.FMFromCP = ForcesMomentsFromCP()

    def update(self):
        pass



class CPFromAerodynamics:
    
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
        Cd = self.Cdf(0) + self.Cdf(1)*M + self.Cdf(2)*AOA + self.Cdf(3)*M**2 + self.Cdf(4)*M*AOA + self.Cdf(5)*AOA**2 + self.Cdf(6)*(M**2)*AOA + self.Cdf(7)*M*(AOA**2) + self.Cdf(8)*AOA**3
        Cl = self.Clf(0) + self.Clf(1)*M + self.Clf(2)*AOA + self.Clf(3)*M**2 + self.Clf(4)*M*AOA + self.Clf(5)*AOA**2 + self.Clf(6)*(M**2)*AOA + self.Clf(7)*M*(AOA**2) + self.Clf(8)*AOA**3
        Cm = self.Cmf(0) + self.Cmf(1)*M + self.Cmf(2)*AOA + self.Cmf(3)*M**2 + self.Cmf(4)*M*AOA + self.Cmf(5)*AOA**2 + self.Cmf(6)*(M**2)*AOA + self.Cmf(7)*M*(AOA**2) + self.Cmf(8)*AOA**3

        # Calculates the forces for the body and the control surfaces
        # body forces
        Fxb = 0
        Fyb = 0
        Fzb = 0

        # surface forces 
        # Canard forces (portCanard (pc), starboardCanard (sc))
        Fpcx = 0
        Fpcy = 0
        Fpcz = 0

        Fscx = 0
        Fscy = 0
        Fscz = 0
 
        # Fin forces (portFin (pf), starboardFin (sf))
        Fpfx = 0
        Fpfy = 0
        Fpfz = 0

        Fsfx = 0
        Fsfy = 0
        Fsfz = 0
        
        # Total forces
        Fsx = 0
        Fsy = 0
        Fsz = 0

        #Engine forces
        FEx = 0
        FEy = 0
        FEz = 0

        return 


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