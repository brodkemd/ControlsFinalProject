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

    def update(self, M, AOA, B, delSurf, delT):
        '''
        Inputs:
        M = Current Mach
        AOA = Current AOA [deg]
        B = Current sideslip, beta [deg]
        delSurf = Array of control surface deflections [delPortCanard, delStarCanard, delPortFin, delStarFin]
        delT = Throttle percentage of full throttle (can come from controller) [0,1] array: [delTx, delTy, delTz]
        delE = Array of engine deflections [deflection from x-axis, rotation]

        '''
        # Defines the coefficients based on the current Mach and AOA
        Cd = self.Cdf(0) + self.Cdf(1)*M + self.Cdf(2)*AOA + self.Cdf(3)*M**2 + self.Cdf(4)*M*AOA + self.Cdf(5)*AOA**2 + self.Cdf(6)*(M**2)*AOA + self.Cdf(7)*M*(AOA**2) + self.Cdf(8)*AOA**3
        Cl = self.Clf(0) + self.Clf(1)*M + self.Clf(2)*AOA + self.Clf(3)*M**2 + self.Clf(4)*M*AOA + self.Clf(5)*AOA**2 + self.Clf(6)*(M**2)*AOA + self.Clf(7)*M*(AOA**2) + self.Clf(8)*AOA**3
        Cm = self.Cmf(0) + self.Cmf(1)*M + self.Cmf(2)*AOA + self.Cmf(3)*M**2 + self.Cmf(4)*M*AOA + self.Cmf(5)*AOA**2 + self.Cmf(6)*(M**2)*AOA + self.Cmf(7)*M*(AOA**2) + self.Cmf(8)*AOA**3
        # Calculates the current velocity based on Mach
        V = M*(BODY.gamma*BODY.R*BODY.Tavg)
        # Calculates the forces for the body and the control surfaces
        # body forces
        Fxb = (-Cd*np.cos(AOA) - Cl*np.sin(AOA))*(0.5*BODY.rhoAvg*(V**2)*BODY.Aref)*np.cos(B)
        Fyb = (-Cd*np.cos(AOA) - Cl*np.sin(AOA))*(0.5*BODY.rhoAvg*(V**2)*BODY.Aref)*np.sin(B)
        Fzb = (-Cd*np.sin(AOA) - Cl*np.cos(AOA))*(0.5*BODY.rhoAvg*(V**2)*BODY.Aref)

        # surface forces 
        # Defines the surface coefficients [delPortCanard, delStarCanard, delPortFin, delStarFin]
        Cds = 2*np.pi*delSurf 
        Cls = 1.28*np.sin(delSurf)

        # Canard forces (portCanard (pc), starboardCanard (sc))
        Fpcx = (-Cds(1)*np.cos(AOA)*np.cos(B) - Cls(1)*np.sin(AOA)*np.sin(B))*(0.5*BODY.rhoAvg*(V**2)*BODY.AsurfC)
        Fpcy = (Cds(1)*np.cos(AOA)*np.sin(B) - Cls(1)*np.sin(AOA)*np.cos(B))*(0.5*BODY.rhoAvg*(V**2)*BODY.AsurfC)
        Fpcz = (-Cd(1)*np.sin(AOA))*(0.5*BODY.rhoAvg*(V**2)*BODY.AsurfC)

        Fscx = (-Cds(1)*np.cos(AOA)*np.cos(B) + Cls(1)*np.sin(AOA)*np.sin(B))*(0.5*BODY.rhoAvg*(V**2)*BODY.AsurfC)
        Fscy = (Cds(1)*np.cos(AOA)*np.sin(B) + Cls(1)*np.sin(AOA)*np.cos(B))*(0.5*BODY.rhoAvg*(V**2)*BODY.AsurfC)
        Fscz = (-Cd(1)*np.sin(AOA))*(0.5*BODY.rhoAvg*(V**2)*BODY.AsurfC)
 
        # Fin forces (portFin (pf), starboardFin (sf))
        Fpfx = (-Cds(1)*np.cos(AOA)*np.cos(B) - Cls(1)*np.sin(AOA)*np.sin(B))*(0.5*BODY.rhoAvg*(V**2)*BODY.AsurfF)
        Fpfy = (Cds(1)*np.cos(AOA)*np.sin(B) - Cls(1)*np.sin(AOA)*np.cos(B))*(0.5*BODY.rhoAvg*(V**2)*BODY.AsurfF)
        Fpfz = (-Cd(1)*np.sin(AOA))*(0.5*BODY.rhoAvg*(V**2)*BODY.AsurfF)

        Fsfx = (-Cds(3)*np.cos(AOA)*np.cos(B) + Cls(3)*np.sin(AOA)*np.sin(B))*(0.5*BODY.rhoAvg*(V**2)*BODY.AsurfF)
        Fsfy = (Cds(3)*np.cos(AOA)*np.sin(B) + Cls(3)*np.sin(AOA)*np.cos(B))*(0.5*BODY.rhoAvg*(V**2)*BODY.AsurfF)
        Fsfz = (-Cd(3)*np.sin(AOA))*(0.5*BODY.rhoAvg*(V**2)*BODY.AsurfF)
        
        # Total forces
        Fsx = Fpcx + Fscx + Fpfx + Fsfx
        Fsy = Fpcy + Fscy + Fpfy + Fsfy
        Fsz = Fpcz + Fscz + Fpfz + Fsfz

        #Engine forces
        FEx = BODY.maxThrust*delT(0)
        FEy = BODY.maxThrust*delT(1)
        FEz = BODY.maxThrust*delT(2)

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