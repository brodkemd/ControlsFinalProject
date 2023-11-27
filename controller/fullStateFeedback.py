from controller.stateSpace import LandingStateSpace, FlipStateSpaceCP, DescentStateSpaceCP, BaseStateSpace
from controller.baseFullStateFeedback import BaseFullStateFeedBack
from tools.rotations import Euler2Quaternion
from parameters.baseClass import Base
from parameters import simulation_parameters as SIM
from tools.loadMathModule import LoadModule

import numpy as np
import os, warnings
warnings.filterwarnings("error")


class FullStateFeedBack(Base):
    global_u_to_total_u = None

    def __init__(self, compute_gains=False) -> None:
        super().__init__()
        print("\nFullStateFeedBack:")
        # self.landing    = Landing(compute_gains=False, add_integrator=False)
        self.descentCP  = DescentCP(compute_gains=compute_gains, add_integrator=False)
        self.flipCP     = FlipCP(compute_gains=False, add_integrator=False)

        print("    Loading Transforms: ", end="")
        cwd = os.path.dirname(__file__)
        LoadModule(os.path.join(cwd, "total.h"), self)
        print("done")

        self.state = 2 # 0 = landing, 1 = flip, 2 = descent

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
        # if self.state == 2:
        #     if x.item(2) > -7000: # determining from position
        #         self.state = 1
        #         self.check = 0
        #         phi_e = 0
        #         theta_e = np.deg2rad(90)
        #         psi_e = 0
        #         self.x_r[6] = Euler2Quaternion(phi_e, theta_e, psi_e).item(0) # e_0
        #         self.x_r[7] = Euler2Quaternion(phi_e, theta_e, psi_e).item(1) # e_1
        #         self.x_r[8] = Euler2Quaternion(phi_e, theta_e, psi_e).item(2) # e_2
        #         self.x_r[9] = Euler2Quaternion(phi_e, theta_e, psi_e).item(3) # e_3

        #     else:
        #         phi_e = 0
        #         theta_e = 0
        #         psi_e = 0
        #         self.x_r[6] = Euler2Quaternion(phi_e, theta_e, psi_e).item(0) # e_0
        #         self.x_r[7] = Euler2Quaternion(phi_e, theta_e, psi_e).item(1) # e_1
        #         self.x_r[8] = Euler2Quaternion(phi_e, theta_e, psi_e).item(2) # e_2
        #         self.x_r[9] = Euler2Quaternion(phi_e, theta_e, psi_e).item(3) # e_3
        
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
        
        # elif self.state == 0:
        #     if x.item(2) > -250 and self.check == 0 and np.abs(x.item(3)) < 5:
        #             self.x_r[0] = -x.item(0)
        #             self.x_r[1] = x.item(1)
        #             self.check = 1
    
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

        return F_E, F_cp, tau, self.x_r.copy()


class DescentCP(DescentStateSpaceCP, BaseFullStateFeedBack):
    def __init__(self, compute_gains=False, add_integrator=False) -> None:
        BaseFullStateFeedBack.__init__(self, "descentCP")
        DescentStateSpaceCP.__init__(self)

        self.zetas    = 0.8*np.ones(len(self.A)//2)
        # self.omega_ns = np.array([2, 1, 0.8, 0.1, 0.1], dtype=float)
        self.omega_ns = np.array([0.2, 0.9, 0.8, 0.1, 0.1], dtype=float)
        self.computePoles()
        if add_integrator:
            #print(self.poles)
            self.poles = self.poles + [-0.0001, -0.0002, -0.0003, -0.0001]
        #print(self.poles)
        # exit()
        self.generateGains(compute_gains=compute_gains, add_integrator=add_integrator)


class FlipCP(FlipStateSpaceCP, BaseFullStateFeedBack):
    def __init__(self, compute_gains=False, add_integrator=False) -> None:
        BaseFullStateFeedBack.__init__(self, "flipCP")
        FlipStateSpaceCP.__init__(self)

        self.zetas    =  0.707*np.ones(len(self.A)//2)
        self.omega_ns = 0.05*np.array([3.0, 1.5, 2.5, 1.25, 1.75, 1.15])

        self.computePoles()
        self.generateGains(compute_gains=compute_gains, add_integrator=add_integrator)


class Landing(LandingStateSpace, BaseFullStateFeedBack):
    def __init__(self, compute_gains=False, add_integrator=False) -> None:
        BaseFullStateFeedBack.__init__(self, "landing")
        LandingStateSpace.__init__(self)

        self.zetas    =  0.97*np.ones(len(self.A)//2)
        self.omega_ns = 0.0251*np.array([3.0, 1.5, 2.5, 1.25, 1.75]) # [3.0, 1.5, 2.5, 1.25, 1.75]

        self.computePoles()
        self.generateGains(compute_gains=compute_gains, add_integrator=add_integrator)