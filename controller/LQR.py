from controller.stateSpace import LandingStateSpace, FlipStateSpaceCP, DescentStateSpaceCP, BaseStateSpace
from controller.baseLQR import BaseFullStateFeedBack
from tools.rotations import Euler2Quaternion, Quaternion2Euler
from parameters.baseClass import Base
from tools.loadMathModule import LoadModule
from controller.deflectioncalc import deflection_calc

import numpy as np
import os


class LQR(Base):
    global_u_to_total_u = None

    def __init__(self, compute_gains=False, rigid_body=False) -> None:
        super().__init__()
        self.rigid_body = rigid_body

        print("\nLQR:")
        self.landing     = Landing(  compute_gains=True, add_integrator=False)
        self.descentCP   = DescentCP(compute_gains=False,          add_integrator=False)
        self.flipCP      = FlipCP(   compute_gains=compute_gains, add_integrator=False)
        self.toFinAngles = deflection_calc()
        print("  Loading Transforms: ", end="")
        cwd = os.path.dirname(__file__)
        LoadModule(os.path.join(cwd, "total.h"), self)
        print("done")

        self.state = 0 # 0 = landing, 1 = flip, 2 = descent

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
            if x.item(2) > -7000: # determining from position
                self.state = 1
                self.check = 0
                phi_e = 0
                theta_e = np.deg2rad(90)
                psi_e = 0
                self.x_r[6] = Euler2Quaternion(phi_e, theta_e, psi_e).item(0) # e_0
                self.x_r[7] = Euler2Quaternion(phi_e, theta_e, psi_e).item(1) # e_1
                self.x_r[8] = Euler2Quaternion(phi_e, theta_e, psi_e).item(2) # e_2
                self.x_r[9] = Euler2Quaternion(phi_e, theta_e, psi_e).item(3) # e_3

            else:
                phi_e = 0
                theta_e = np.deg2rad(0)
                psi_e = 0
                self.x_r[6] = Euler2Quaternion(phi_e, theta_e, psi_e).item(0) # e_0
                self.x_r[7] = Euler2Quaternion(phi_e, theta_e, psi_e).item(1) # e_1
                self.x_r[8] = Euler2Quaternion(phi_e, theta_e, psi_e).item(2) # e_2
                self.x_r[9] = Euler2Quaternion(phi_e, theta_e, psi_e).item(3) # e_3
        
        # elif self.state == 1:
        #     phi_e = 0
        #     theta_e = 90
        #     psi_e = 0

        #     e = x[6:10]
        #     theta = np.rad2deg(Quaternion2Euler(e)[1])
        #     if abs(theta_e - theta) < 10:
        #         theta_e = np.deg2rad(theta_e)
        #         self.state = 0
        #         self.x_r[6] = Euler2Quaternion(phi_e, theta_e, psi_e).item(0) # e_0
        #         self.x_r[7] = Euler2Quaternion(phi_e, theta_e, psi_e).item(1) # e_1
        #         self.x_r[8] = Euler2Quaternion(phi_e, theta_e, psi_e).item(2) # e_2
        #         self.x_r[9] = Euler2Quaternion(phi_e, theta_e, psi_e).item(3) # e_3
        
        elif self.state == 0:
            if x.item(2) > -250 and self.check == 0 and np.abs(x.item(3)) < 5:
                    self.x_r[0] = x.item(0)
                    self.x_r[1] = x.item(1)
                    self.check = 1
                    print('landing sequence')
    
        if self.state != last_state:
            print("\nChanged state from", last_state, "to", self.state, end="\n\n")

        if self.state == 0 or self.state == 1:
            u = self.landing.update(x, self.x_r)
        # elif self.state == 1:
        #     u = self.flipCP.update(x, self.x_r)
        elif self.state == 2:
            u = self.descentCP.update(x, self.x_r)
        
        if u[0] < 0:
            u[0] = 0
        if u[1] > self.MaxT*.33:
            u[1] = self.MaxT*0.33
        elif u[1] < -self.MaxT*0.33:
            u[1] = -self.MaxT*0.33
        if u[2] > self.MaxT*.33:
            u[2] = self.MaxT*0.33
        elif u[2] < -self.MaxT*0.33:
            u[2] = -self.MaxT*0.33
        
        Total_lat = np.sqrt(u[1]**2 + u[2]**2)
        if u[0] < Total_lat/np.tan(np.deg2rad(45)):
            u[0] = Total_lat/np.tan(np.deg2rad(45))

        u_total = self.global_u_to_total_u @ u
        F_E                   = u_total[0:3]
        F_cp_port_canard      = u_total[3:6]
        F_cp_starboard_canard = u_total[6:9]

        F_cp_port_fin         = u_total[9:12]
        F_cp_starboard_fin    = u_total[12:15]

        F_cp_canard           = F_cp_port_canard + F_cp_starboard_canard
        F_cp_fin              = F_cp_port_fin + F_cp_starboard_fin

        F_cp = F_cp_fin + F_cp_canard

        tau_cp_port_canard      = np.cross(     self.r_cp_port_canard,      F_cp_port_canard)
        tau_cp_starboard_canard = np.cross(self.r_cp_starboard_canard, F_cp_starboard_canard)
        tau_cp_port_fin         = np.cross(        self.r_cp_port_fin,         F_cp_port_fin)
        tau_cp_starboard_fin    = np.cross(   self.r_cp_starboard_fin,    F_cp_starboard_fin)

        tau = tau_cp_port_canard + tau_cp_port_fin + tau_cp_starboard_canard + tau_cp_starboard_fin
        # print(tau)
        angles = self.toFinAngles.calc_def(tau, x)
        
        for i in range(0,len(angles)):
            if angles[i] < -np.pi:
                angles[i] = -np.pi
            elif angles[i] > -np.pi/2:
                angles[i] = -np.pi/2
        # print(np.rad2deg(angles))

        # print(x[2])
        return F_E, angles,self.x_r.copy()

        # if self.rigid_body:
        #     return F_E, F_cp, tau, self.x_r.copy()
        # else:
        #     angles = self.toFinAngles.calc_def(tau, x)
        #     for i in range(0,len(angles)):
        #         if np.abs(angles[i]) < 0:
        #             angles[i] = 0
        #         elif np.abs(angles[i]) > np.abs(np.deg2rad(-180)):
        #             angles[i] = np.deg2rad(-180)

        #     print(np.rad2deg(angles))
        #     return F_E, angles, self.x_r.copy()


class DescentCP(DescentStateSpaceCP, BaseFullStateFeedBack):
    def __init__(self, compute_gains=False, add_integrator=False) -> None:
        BaseFullStateFeedBack.__init__(self, "descentCP")
        DescentStateSpaceCP.__init__(self)

        # self.Q_diagonal = 1/(1**2)*np.ones(self.A.shape[0])
        self.Q_diagonal = np.array([1E-6,1E-6,1E-6,1E-6,1,1E5,1,1,1,1])
        #self.Q_diagonal[0]*=1
        self.R_diagonal = 1E-6*np.ones(self.B.shape[1])
        # self.R_diagonal = np.array([1E5,1E5,1E5,1E-5,1E-5,1E-5,1,1])

        self.generateGains(compute_gains=compute_gains, add_integrator=add_integrator)

    # def update(self, _x, _x_r):
    #     # x_vars = p_n,p_e,u,v,e_1,e_2,e_3,p,q,r # order matters
    #     error   = self.global_x_r_to_local_x_r @ (_x -     _x_r)
    #     x_tilde = self.global_x_to_local_x     @ (_x - self.x_e)
    #     # if abs(x_tilde.item(0)) > 1 and abs(x_tilde.item(1)) < 1:
    #     #     self.x_e[6] = Euler2Quaternion(0, 0, np.deg2rad(90)).item(0)
    #     #     self.x_e[7] = Euler2Quaternion(0, 0, np.deg2rad(90)).item(1)
    #     #     self.x_e[8] = Euler2Quaternion(0, 0, np.deg2rad(90)).item(2)
    #     #     self.x_e[9] = Euler2Quaternion(0, 0, np.deg2rad(90)).item(3)
    #     # else:
    #     #     self.x_e[6] = Euler2Quaternion(0, 0,  np.deg2rad(0)).item(0)
    #     #     self.x_e[7] = Euler2Quaternion(0, 0,  np.deg2rad(0)).item(1)
    #     #     self.x_e[8] = Euler2Quaternion(0, 0,  np.deg2rad(0)).item(2)
    #     #     self.x_e[9] = Euler2Quaternion(0, 0,  np.deg2rad(0)).item(3)

    #     if self.has_integrator:
    #         self.integrateError(error)
    #         x_tilde_with_integral = np.hstack((x_tilde, self.integrator))
    #         u = self.local_u_to_global_u @ (-1*self.K @ x_tilde_with_integral)
    #     else:
    #         u = self.local_u_to_global_u @ (-1*self.K @ x_tilde)
    #     return self.u_e + u


class FlipCP(FlipStateSpaceCP, BaseFullStateFeedBack):
    def __init__(self, compute_gains=False, add_integrator=False) -> None:
        BaseFullStateFeedBack.__init__(self, "flipCP")
        FlipStateSpaceCP.__init__(self)

        if add_integrator:
            self.Q_diagonal = np.ones(self.A.shape[0] + self.C.shape[0])
            self.R_diagonal = np.ones(self.B.shape[0] + self.C.shape[0])
        else:
            self.Q_diagonal = np.ones(self.A.shape[0])
            self.R_diagonal = np.ones(self.B.shape[0])

        self.generateGains(compute_gains=compute_gains, add_integrator=add_integrator)


class Landing(LandingStateSpace, BaseFullStateFeedBack):
    def __init__(self, compute_gains=False, add_integrator=False) -> None:
        BaseFullStateFeedBack.__init__(self, "landing")
        LandingStateSpace.__init__(self)

        # self.Q_diagonal = 1E12/(0.1**2)*np.ones(self.A.shape[0])
        self.Q_diagonal = np.array([1E3,1E-2,1E3,1E-8,1E-8,1E-8,1E11,1,1,1])
        self.R_diagonal = 5E1*np.ones(self.B.shape[1])
        self.computeQR()
        if add_integrator:
            self.Q_diagonal = np.ones(self.A.shape[0] + self.C.shape[0])
            self.R_diagonal = np.ones(self.B.shape[0] + self.C.shape[0])
        else:
            # print("     No integrator")
            self.Q_diagonal = 1/(0.001**2)*np.ones(self.A.shape[0])
            self.R_diagonal = np.ones(self.B.shape[0])


        self.generateGains(compute_gains=compute_gains, add_integrator=add_integrator)