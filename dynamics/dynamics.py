import numpy as np
from numpy import sin, cos, tan
import parameters.simulation_parameters as SIM
import parameters.body_parameters as BODY

class Dynamics:
    def __init__(self, initial_condition:np.ndarray):
        # Initial state conditions, see state format
        self.state = initial_condition
        self.Ts    = SIM.ts_simulation # sample rate of system
        self.limit = 1.0 # input saturation limit
        self.m     = BODY.mass
        self.J     = BODY.J
        self.J_inv = np.linalg.inv(self.J)


    def f(self, _state, _u):
        # for system xdot = f(x,u), return f(x,u)
        # u  = [[f_x], [f_y], [f_z], [l], [m], [n]]
        # force in body frame: <f_x, f_y, f_z>
        # moment in body frame: <l, m, n>
        # pn      = _state.item(0) # position north
        # pe      = _state.item(1) # position east
        # pd      = _state.item(2) # position down
        # u       = _state.item(3) # velocity measured along i body
        # v       = _state.item(4) # velocity measured along j body
        # w       = _state.item(5) # velocity measured along k body
        e_0     = _state.item(6)
        e_1     = _state.item(7)
        e_2     = _state.item(8)
        e_3     = _state.item(9)
        p       = _state.item(10) # roll rate measured along i body
        q       = _state.item(11) # pitch rate measured along j body
        r       = _state.item(12) # yaw rate measured along k body
        # f_x     = _u.item(0)
        # f_y     = _u.item(1)
        # f_z     = _u.item(2)
        # l       = _u.item(3)
        # m       = _u.item(4)
        # n       = _u.item(5)

        V     = _state[3:6]
        omega = _state[10:13]
        e     = _state[6:10]

        F     = _u[0:3]
        tau   = _u[3:6]

        M1 = np.array([
            [e_1**2 + e_0**2 - e_2**2 - e_3**2, 2*(e_1*e_2 - e_3*e_0), 2*(e_1*e_3 + e_2*e_0)],
            [2*(e_1*e_2 + e_3*e_0), e_2**2 + e_0**2 - e_1**2 - e_3**2, 2*(e_2*e_3 - e_1*e_0)],
            [2*(e_1*e_3 - e_2*e_0), 2*(e_2*e_3 + e_1*e_0), e_3**2 + e_0**2 - e_1**2 - e_2**2]
        ])

        M2 = 1/2*np.array([
            [0, -p, -q, -r],
            [p, 0, r, -q],
            [q, -r, 0, p],
            [r, q, -p, 0]
        ])

        P_g_dot   =  M1 @ V
        V_dot     = -1*np.cross(omega, V) + 1/self.m*F
        e_dot     =  M2 @ e
        omega_dot = -1*self.J_inv @ (np.cross(omega, self.J @ omega)) + self.J_inv @ tau

        return np.hstack((P_g_dot, V_dot, e_dot, omega_dot))


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