#

#delta_sc: (-21755.5291261093*AOA*A_c*Mz - 909465.715980001*A_c**2*V**2*rho + 1155468.73751558*A_c*A_f*V**2*rho + 75640.0*A_c*Fy - 21755.5291261093*A_c*My + 3904.0*A_c*Mz - 96100.0*A_f*Fy - 4960.0*A_f*Mz)/(424667.928541654*A_c**2*V**2*rho - 539537.122327511*A_c*A_f*V**2*rho)
#delta_pc: (25.0*AOA*Mz - 488.0*A_c*V**2*rho + 620.0*A_f*V**2*rho + 25.0*My)/(488.0*A_c*V**2*rho - 620.0*A_f*V**2*rho)
#delta_sf: (28434476.5678249*AOA**2*A_c*Mz + 316815354.09096*AOA*A_c**2*V**2*rho - 1303051520.57752*AOA*A_c*A_f*V**2*rho - 9243940.0*AOA*A_c*Fy + 28434476.5678249*AOA*A_c*My - 5186464.0*AOA*A_c*Mz + 1144128955.39408*AOA*A_f**2*V**2*rho + 11744350.0*AOA*A_f*Fy + 6589360.0*AOA*A_f*Mz - 13517600.0*A_c*Mx + 17174000.0*A_f*Mx)/(573301703.531233*AOA*A_c*A_f*V**2*rho - 728375115.14214*AOA*A_f**2*V**2*rho)
#delta_pf: (-28434476.5678249*AOA**2*A_c*Mz - 316815354.09096*AOA*A_c**2*V**2*rho - 498028899.526717*AOA*A_c*A_f*V**2*rho + 89617540.0*AOA*A_c*Fy - 28434476.5678249*AOA*A_c*My - 83936.0*AOA*A_c*Mz + 1144128955.39408*AOA*A_f**2*V**2*rho - 113858350.0*AOA*A_f*Fy + 106640.0*AOA*A_f*Mz + 13517600.0*A_c*Mx - 17174000.0*A_f*Mx)/(573301703.531233*AOA*A_c*A_f*V**2*rho - 728375115.14214*AOA*A_f**2*V**2*rho)
import sys
sys.path.append('.')
import parameters.body_parameters as BODY
import numpy as np

class deflection_calc:
    def __init__(self):
        self.A_c = BODY.AsurfC
        self.A_f = BODY.AsurfF
        self.rho = BODY.rhoAvg

    def calc_def(self,M,state):
        AOA = np.arctan2(state[5],state[3])
        V = np.sqrt(state[3]**2 + state[4]**2 + state[5]**2)
        Fy = 0
        A_c = self.A_c
        A_f = self.A_f
        rho = self.rho
        AOA -= np.deg2rad(90)
        Mx = M[0]
        My = M[1]
        Mz = M[2]
        if AOA == 0 or (np.abs(AOA)-np.deg2rad(0.05)) < 0:
            delta_sc = (-21755.5291261093*AOA*self.A_c*Mz - 909465.715980001*self.A_c**2*V**2*self.rho + 1155468.73751558*self.A_c*self.A_f*V**2*self.rho + 75640.0*self.A_c*Fy - 21755.5291261093*self.A_c*My + 3904.0*self.A_c*Mz - 96100.0*self.A_f*Fy - 4960.0*self.A_f*Mz)/(424667.928541654*self.A_c**2*V**2*self.rho - 539537.122327511*self.A_c*self.A_f*V**2*self.rho)
            delta_pc = (25.0*AOA*Mz - 488.0*self.A_c*V**2*self.rho + 620.0*self.A_f*V**2*self.rho + 25.0*My)/(488.0*self.A_c*V**2*self.rho - 620.0*self.A_f*V**2*self.rho)
            delta_sf = 0
            delta_pf = 0
            # delta_sc = np.abs(delta_pc + np.deg2rad(90))
            # delta_pc += np.deg2rad(90)
        else:
            delta_sc = (-21755.5291261093*AOA*A_c*Mz - 909465.715980001*A_c**2*V**2*rho + 1155468.73751558*A_c*A_f*V**2*rho + 75640.0*A_c*Fy - 21755.5291261093*A_c*My + 3904.0*A_c*Mz - 96100.0*A_f*Fy - 4960.0*A_f*Mz)/(424667.928541654*A_c**2*V**2*rho - 539537.122327511*A_c*A_f*V**2*rho)
            delta_pc = (25.0*AOA*Mz - 488.0*A_c*V**2*rho + 620.0*A_f*V**2*rho + 25.0*My)/(488.0*A_c*V**2*rho - 620.0*A_f*V**2*rho)
            delta_sf = (28434476.5678249*AOA**2*A_c*Mz + 316815354.09096*AOA*A_c**2*V**2*rho - 1303051520.57752*AOA*A_c*A_f*V**2*rho - 9243940.0*AOA*A_c*Fy + 28434476.5678249*AOA*A_c*My - 5186464.0*AOA*A_c*Mz + 1144128955.39408*AOA*A_f**2*V**2*rho + 11744350.0*AOA*A_f*Fy + 6589360.0*AOA*A_f*Mz - 13517600.0*A_c*Mx + 17174000.0*A_f*Mx)/(573301703.531233*AOA*A_c*A_f*V**2*rho - 728375115.14214*AOA*A_f**2*V**2*rho)
            delta_pf = (-28434476.5678249*AOA**2*A_c*Mz - 316815354.09096*AOA*A_c**2*V**2*rho - 498028899.526717*AOA*A_c*A_f*V**2*rho + 89617540.0*AOA*A_c*Fy - 28434476.5678249*AOA*A_c*My - 83936.0*AOA*A_c*Mz + 1144128955.39408*AOA*A_f**2*V**2*rho - 113858350.0*AOA*A_f*Fy + 106640.0*AOA*A_f*Mz + 13517600.0*A_c*Mx - 17174000.0*A_f*Mx)/(573301703.531233*AOA*A_c*A_f*V**2*rho - 728375115.14214*AOA*A_f**2*V**2*rho)
            # delta_sc = np.abs(delta_pc + np.deg2rad(90))
            # delta_pc += np.deg2rad(90)
            # delta_sf = np.abs(delta_pf + np.deg2rad(90))
            # delta_pf = np.abs(delta_pf + np.deg2rad(90))
        
        return np.array([delta_pc,delta_sc,delta_pf,delta_sf])

