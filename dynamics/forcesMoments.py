import numpy as np
import parameters.body_parameters as BODY

#delta_sc: (-21755.5291261093*AOA*A_c*Mz - 909465.715980001*A_c**2*V**2*rho + 1155468.73751558*A_c*A_f*V**2*rho + 75640.0*A_c*Fy - 21755.5291261093*A_c*My + 3904.0*A_c*Mz - 96100.0*A_f*Fy - 4960.0*A_f*Mz)/(424667.928541654*A_c**2*V**2*rho - 539537.122327511*A_c*A_f*V**2*rho)
#delta_pc: (25.0*AOA*Mz - 488.0*A_c*V**2*rho + 620.0*A_f*V**2*rho + 25.0*My)/(488.0*A_c*V**2*rho - 620.0*A_f*V**2*rho)
#delta_sf: (28434476.5678249*AOA**2*A_c*Mz + 316815354.09096*AOA*A_c**2*V**2*rho - 1303051520.57752*AOA*A_c*A_f*V**2*rho - 9243940.0*AOA*A_c*Fy + 28434476.5678249*AOA*A_c*My - 5186464.0*AOA*A_c*Mz + 1144128955.39408*AOA*A_f**2*V**2*rho + 11744350.0*AOA*A_f*Fy + 6589360.0*AOA*A_f*Mz - 13517600.0*A_c*Mx + 17174000.0*A_f*Mx)/(573301703.531233*AOA*A_c*A_f*V**2*rho - 728375115.14214*AOA*A_f**2*V**2*rho)
#delta_pf: (-28434476.5678249*AOA**2*A_c*Mz - 316815354.09096*AOA*A_c**2*V**2*rho - 498028899.526717*AOA*A_c*A_f*V**2*rho + 89617540.0*AOA*A_c*Fy - 28434476.5678249*AOA*A_c*My - 83936.0*AOA*A_c*Mz + 1144128955.39408*AOA*A_f**2*V**2*rho - 113858350.0*AOA*A_f*Fy + 106640.0*AOA*A_f*Mz + 13517600.0*A_c*Mx - 17174000.0*A_f*Mx)/(573301703.531233*AOA*A_c*A_f*V**2*rho - 728375115.14214*AOA*A_f**2*V**2*rho)


class ForcesMoments:
    def __init__(self) -> None:
        #Need locations from the c.g. of the canards and fins
        self.FMFromCP = ForcesMomentsFromCP()
        self.CPFromAerodynamics = CPFromAerodynamics()

    def update(self, _state, F_E, angles):
        return self.CPFromAerodynamics.update(_state, F_E, angles)



class CPFromAerodynamics:
    
    def __init__(self) -> None:
        # Surface fit coefficients for the force and moments
        self.Cdf = np.array([0.03605, -0.01182, -0.0001243, 0.01985, -0.0001574, 7.297E-5, 0.001095, -8.369E-6, -6.47E-7])
        self.Clf = np.array([0.01872, -0.07913, -0.002173, 0.06861, 0.004314, 0.0001143, -0.002112, -2.849E-5, -1.013E-6])
        self.Cmf = np.array([-0.02598, 0.1033, 0.0002977, -0.08253, -0.003769, -0.0002425, -0.001543, 5.302E-5, 2.275E-6])
       
        pass

    def update(self, state, F_E, delSurf):
        '''
        Inputs:
        M = Current Mach
        AOA = Current AOA [deg]
        B = Current sideslip, beta [deg]
        delSurf = Array of control surface deflections [delPortCanard, delStarCanard, delPortFin, delStarFin]
        delT = Throttle percentage of full throttle (can come from controller) [0,1] array: [delTx, delTy, delTz]
        delE = Array of engine deflections [deflection from x-axis, rotation]

        Inverted 

        '''
        V = np.sqrt(state[3]**2 + state[4]**2 + state[5]**2)
        AOA = np.arctan2(state[5],state[3])
        B = np.arctan2(state[4],V)
        AOA = np.rad2deg(AOA)
        # Calculates the current MACH
        M = V/(BODY.gamma*BODY.R*BODY.Tavg)
        # Defines the coefficients based on the current Mach and AOA
        Cd = self.Cdf[0] + self.Cdf[1]*M + self.Cdf[2]*AOA + self.Cdf[3]*M**2 + self.Cdf[4]*M*AOA + self.Cdf[5]*AOA**2 + self.Cdf[6]*(M**2)*AOA + self.Cdf[7]*M*(AOA**2) + self.Cdf[8]*AOA**3
        Cl = self.Clf[0] + self.Clf[1]*M + self.Clf[2]*AOA + self.Clf[3]*M**2 + self.Clf[4]*M*AOA + self.Clf[5]*AOA**2 + self.Clf[6]*(M**2)*AOA + self.Clf[7]*M*(AOA**2) + self.Clf[8]*AOA**3
        Cm = self.Cmf[0] + self.Cmf[1]*M + self.Cmf[2]*AOA + self.Cmf[3]*M**2 + self.Cmf[4]*M*AOA + self.Cmf[5]*AOA**2 + self.Cmf[6]*(M**2)*AOA + self.Cmf[7]*M*(AOA**2) + self.Cmf[8]*AOA**3
        # Calculates the forces for the body and the control surfaces
        AOA = np.deg2rad(AOA)
        # body forces
        Fbx = (-Cd*np.cos(AOA) - Cl*np.sin(AOA))*(0.5*BODY.rhoAvg*(V**2)*BODY.Aref)*np.cos(B)
        Fby = (-Cd*np.cos(AOA) - Cl*np.sin(AOA))*(0.5*BODY.rhoAvg*(V**2)*BODY.Aref)*np.sin(B)
        Fbz = (-Cd*np.sin(AOA) - Cl*np.cos(AOA))*(0.5*BODY.rhoAvg*(V**2)*BODY.Aref)

        # surface forces 
        # Defines the surface coefficients from flat plate thy [delPortCanard, delStarCanard, delPortFin, delStarFin]
        # Cls = 2*np.pi*delSurf 
        # Cds = 1.28*np.sin(delSurf)
        Cls = 2*np.pi*(np.pi/2 + delSurf)
        Cds = 1.28*(1 + delSurf)

        rho = BODY.rhoAvg
        A_c = BODY.AsurfC
        A_f = BODY.AsurfF

        # Canard forces (portCanard (pc), starboardCanard (sc))
        Fpcx = (Cds[0]*AOA)*(0.5*rho*V**2*A_c)
        Fpcy = Cls[0]*(0.5*rho*V**2*A_c)
        Fpcz = (-Cds[0] + Cls[0]*AOA)*(0.5*rho*V**2*A_c)

        Fscx = (Cds[1]*AOA)*(0.5*rho*V**2*A_c)
        Fscy = Cls[1]*(0.5*rho*V**2*A_c)
        Fscz = (-Cds[1] + Cls[1]*AOA)*(0.5*rho*V**2*A_c)
 
        # Fin forces (portFin (pf), starboardFin (sf))
        Fpfx = (Cds[2]*AOA)*(0.5*rho*V**2*A_f)
        Fpfy = Cls[2]*(0.5*rho*V**2*A_f)
        Fpfz = (-Cds[2] + Cls[2]*AOA)*(0.5*rho*V**2*A_f)

        Fsfx = (Cds[3]*AOA)*(0.5*rho*V**2*A_f)
        Fsfy = Cls[3]*(0.5*rho*V**2*A_f)
        Fsfz = (-Cds[3] + Cls[3]*AOA)*(0.5*rho*V**2*A_f)
        
        # Total forces
        Fsx = Fpcx + Fscx + Fpfx + Fsfx
        Fsy = Fpcy + Fscy + Fpfy + Fsfy
        Fsz = Fpcz + Fscz + Fpfz + Fsfz

        #Engine forces
        # FEx = BODY.maxThrust*delT[0]
        # FEy = BODY.maxThrust*delT[1]
        # FEz = BODY.maxThrust*delT[2]

        # radii vectors of all the surfaces (from COM)
        # Port Canard [rpc] = (15.25, -6.535, 0) (x,y,z)
        # Starboard Canard [rsc] = (15.25, 6.535, 0) (x,y,z)
        # Port Fin [rpf] = (-19.375,-6.75, 0) (x,y,z)
        # Starboard Fin [rsf] = (-19.375, 6.75, 0) (x,y,z)
        rpc = np.array([15.25, 6.535, 0])
        rsc = np.array([15.25, -6.535, 0]) 
        rpf = np.array([-19.375, 6.75, 0])
        rsf = np.array([-19.375, -6.75, 0])
        # rb = np.array([0.0, 0.0, 0.0])
        re = BODY.r_E
        
        # Sets the force vectors
        Fpc = np.array([Fpcx,Fpcy,Fpcz])
        Fsc = np.array([Fscx,Fscy,Fscz])
        Fpf = np.array([Fpfx,Fpfy,Fpfz])
        Fsf = np.array([Fsfx,Fsfy,Fsfz])
        Fb = np.array([Fbx, Fby, Fbz])
        # Fe = np.array([FEx, FEy, FEz])
        self.r_E = BODY.r_E
        self.m   = BODY.mass
        self.g   = BODY.gravity
        e_0     = state.item(6)
        e_1     = state.item(7)
        e_2     = state.item(8)
        e_3     = state.item(9)
        F_g = self.m*self.g*np.array([
            2*(e_1*e_3 - e_2*e_0),
            2*(e_2*e_3 +  e_1*e_0),
            e_3**2 + e_0**2 - e_1**2 - e_2**2
        ])

        # Calculates the radius vector of the center of pressure
        PortCanard = np.cross(Fpc,np.cross(rpc,Fpc))
        StarCanard = np.cross(Fsc,np.cross(rsc,Fsc))
        PortFin = np.cross(Fpf,np.cross(rpf,Fpf))
        StarFin = np.cross(Fsf,np.cross(rsf,Fsf))
        normSqr = (np.linalg.norm(np.array([Fsx,Fsy,Fsz])))**2 
        COPrad = (PortCanard + StarCanard + PortFin + StarFin)/normSqr

        # Calculates the total forces and moments
        # Forces (just body acting at COM)
        F = Fb + F_E + F_g + Fpc + Fsc + Fpf + Fsf
        # Body moment using the Cm eqn (moments only in y-dir)
        Mb = np.array([0, Cm*(0.5*BODY.rhoAvg*(V**2)*BODY.Aref), 0])
        # Moments (surface + engine torques)
        Me = np.cross(re,F_E)
        Mpc = np.cross(rpc,Fpc)
        Msc = np.cross(rsc,Fsc)
        Mpf = np.cross(rpf,Fpf)
        Msf = np.cross(rsf,Fsf)
        M = Me + Mpc + Msc + Mpf + Msf + Mb
        
        return np.hstack((F, M))


class ForcesMomentsFromCP:
    def __init__(self) -> None:
        self.r_E = BODY.r_E
        self.m   = BODY.mass
        self.g   = BODY.gravity

    def update(self, _state, F_E, F_cp, tau):
        e_0     = _state.item(6)
        e_1     = _state.item(7)
        e_2     = _state.item(8)
        e_3     = _state.item(9)

        F_g = self.m*self.g*np.array([
            2*(e_1*e_3 - e_2*e_0),
            2*(e_2*e_3 +  e_1*e_0),
            e_3**2 + e_0**2 - e_1**2 - e_2**2
        ])

        tau = np.cross(self.r_E, F_E) + tau
        F   = F_E + F_g + F_cp

        return np.hstack((F, tau))