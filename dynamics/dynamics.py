import numpy as np
from numpy import sin, cos, tan
import parameters.simulation_parameters as SIM
import parameters.body_parameters as BODY

class Dynamics:
    def __init__(self, initial_condition:np.ndarray):
        """
        state = [
            pn, # position north
            pe, # position east
            pd, # position down
            u, # velocity measured along i body
            v, # velocity measured along j body
            w, # velocity measured along k body
            phi, # roll angle
            theta, # pitch angle
            psi, # yaw angle
            p, # roll rate measured along i body
            q, # pitch rate measured along j body
            r # yaw rate measured along k body
        ]
        
        """
        # Initial state conditions, see state format
        self.state = initial_condition
        self.Ts = np.float64(SIM.ts_simulation) # sample rate of system
        self.limit = 1.0 # input saturation limit

        self.m = BODY.mass
        self.J_xx = MAV.jx
        self.J_yy = MAV.jy
        self.J_zz = MAV.jz
        self.J_xz = MAV.jxz
        self.Gamma = self.J_xx*self.J_zz - self.J_xz**2
        
        self.Gammas = np.array([
            (self.J_xz*(self.J_xx - self.J_yy + self.J_zz))/self.Gamma,
            (self.J_zz*(self.J_zz - self.J_yy) + self.J_xz**2)/self.Gamma,
            self.J_zz/self.Gamma,
            self.J_xz/self.Gamma,
            (self.J_zz - self.J_xx)/self.J_yy,
            self.J_xz/self.J_yy,
            ((self.J_xx - self.J_yy)*self.J_xx + self.J_xz**2)/self.Gamma,
            self.J_xx/self.Gamma
        ])


    def f(self, _state, quaternion, _u):
        # for system xdot = f(x,u), return f(x,u)
        # u  = [[f_x], [f_y], [f_z], [l], [m], [n]]
        # force in body frame: <f_x, f_y, f_z>
        # moment in body frame: <l, m, n>
        
        _state = _state.T[0]
        _u = _u.T[0]

        pn      = _state.item(0) # position north
        pe      = _state.item(1) # position east
        pd      = _state.item(2) # position down
        u       = _state.item(3) # velocity measured along i body
        v       = _state.item(4) # velocity measured along j body
        w       = _state.item(5) # velocity measured along k body
        p       = _state.item(6) # roll rate measured along i body
        q       = _state.item(7) # pitch rate measured along j body
        r       = _state.item(8) # yaw rate measured along k body

        f_x     = _u.item(0)
        f_y     = _u.item(1)
        f_z     = _u.item(2)
        l       = _u.item(3)
        m       = _u.item(4)
        n       = _u.item(5)

        e_0 = quaternion.item(0)
        e_1 = quaternion.item(1)
        e_2 = quaternion.item(2)
        e_3 = quaternion.item(3)

        M1 = np.array([
            [e_1**2 + e_0**2 - e_2**2 - e_3**2, 2*(e_1*e_2 - e_3*e_0), 2*(e_1*e_3 + e_2*e_0)],
            [2*(e_1*e_2 + e_3*e_0), e_2**2 + e_0**2 - e_1**2 - e_3**2, 2*(e_2*e_3 - e_1*e_0)],
            [2*(e_1*e_3 - e_2*e_0), 2*(e_2*e_3 + e_1*e_0), e_3**2 + e_0**2 - e_1**2 - e_2**2]
        ])
        M2 = np.array([
            [r*v - q*w],
            [p*w - r*u],
            [q*u - p*v]
        ])
        M3 = np.array([
            [0, -p, -q, -r],
            [p, 0, r, -q],
            [q, -r, 0, p],
            [r, q, -p, 0]
        ])
        M4 = np.array([
            [self.Gammas[0]*p*q - self.Gammas[1]*q*r],
            [self.Gammas[4]*p*r - self.Gammas[5]*(p**2 - r**2)],
            [self.Gammas[6]*p*q - self.Gammas[0]*q*r]
        ])

        x1 = np.array([[u], [v], [w]])
        x2 = 1
        x3 = np.array([[p], [q], [r]])
        x4 = 1

        B1 = np.zeros((3, 1))
        B2 = 1/self.m*np.array([
            [f_x],
            [f_y],
            [f_z]
        ])
        B3 = np.zeros((3, 1))
        B4 = np.array([
            [self.Gammas[2]*l + self.Gammas[3]*n],
            [1/self.J_yy*m],
            [self.Gammas[3]*l + self.Gammas[7]*n]
        ])

        b1 = np.dot(M1, x1) + B1
        b2 = np.dot(M2, x2) + B2
        b3 = np.dot(M3, x3) + B3
        b4 = np.dot(M4, x4) + B4

        result = b1[:,0]
        result = np.append(result, b2[:,0])
        result = np.append(result, b3[:,0])
        result = np.append(result, b4[:,0])

        return np.atleast_2d(result.copy()).T


    def h(self, u):
        return self.state.copy()


    def update(self, u):
        # This is the external method that takes the input u(t)
        # and returns the output y(t).
        self.rk4_step(u) # propagate the state by one time step
        return self.h(u) # compute the output at the current state


    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)