from sympy import *
import os
#from tools.rotations import Euler2Quaternion
# from sympy.printing.latex import LatexPrinter
# from sympy.core.function import UndefinedFunction

# class MyLatexPrinter(LatexPrinter):
#     """Print derivative of a function of symbols in a shorter form.
#     """
#     def _print_Derivative(self, expr):
#         function, *vars = expr.args
#         if not isinstance(type(function), UndefinedFunction) or \
#            not all(isinstance(i, Symbol) for i in vars):
#             return super()._print_Derivative(expr)

#         # If you want the printer to work correctly for nested
#         # expressions then use self._print() instead of str() or latex().
#         # See the example of nested modulo below in the custom printing
#         # method section.
#         return "{}".format(
#             self._print(Symbol(function.func.__name__)))

# #init_printing(derivative='dot')
# # Initialize the pretty printing system
# init_printing()

# # Define a function to customize derivative printing using dot notation
# def custom_derivative_printer(expr, **kwargs):
#     return latex(expr, symbol_names={Derivative: lambda x: f'{x.args[0]}\\dot{{{x.args[1]}}}'})

# # Change the way derivatives are represented globally
# init_printing(latex_printer=custom_derivative_printer)

cwd = os.path.join(os.path.dirname(__file__), "sympy_output")
def write(s, _file:str):
    if not isinstance(s, str):
        s = latex(s)

    _file = os.path.join(cwd, "2_"+_file)
    with open(_file, "w") as file:
        file.write(s)

def Euler2Quaternion(phi, theta, psi):
    """
    Converts an euler angle attitude to a quaternian attitude
    :param euler: Euler angle attitude in a np.matrix(phi, theta, psi)
    :return: Quaternian attitude in np.array(e0, e1, e2, e3)
    """

    e0 = cos(psi/2.0) * cos(theta/2.0) * cos(phi/2.0) + sin(psi/2.0) * sin(theta/2.0) * sin(phi/2.0)
    e1 = cos(psi/2.0) * cos(theta/2.0) * sin(phi/2.0) - sin(psi/2.0) * sin(theta/2.0) * cos(phi/2.0)
    e2 = cos(psi/2.0) * sin(theta/2.0) * cos(phi/2.0) + sin(psi/2.0) * cos(theta/2.0) * sin(phi/2.0)
    e3 = sin(psi/2.0) * cos(theta/2.0) * cos(phi/2.0) - cos(psi/2.0) * sin(theta/2.0) * sin(phi/2.0)

    return Matrix([e0, e1, e2, e3])

def Diff(f):
    return diff(f, t)

t = symbols("t")
f_E_x, f_E_y, f_E_z, f_g_x, f_g_y, f_g_z, f_cp_x, f_cp_y, f_cp_z, r_E_x, r_E_y, r_E_z, r_cp_x, r_cp_y, r_cp_z = symbols("f_E_x, f_E_y, f_E_z, f_g_x, f_g_y, f_g_z, f_cp_x, f_cp_y, f_cp_z, r_E_x, r_E_y, r_E_z, r_cp_x, r_cp_y, r_cp_z")
p_n,p_e,p_d,e_0,e_1,e_2,e_3,u,v,w,p,q,r = symbols("p_n,p_e,p_d,e_0,e_1,e_2,e_3,u,v,w,p,q,r")
p_n_f = Function("p_n", real=True)(t)
p_e_f = Function("p_e", real=True)(t)
p_d_f = Function("p_d", real=True)(t)
e_0_f = Function("e_0", real=True)(t)
e_1_f = Function("e_1", real=True)(t)
e_2_f = Function("e_2", real=True)(t)
e_3_f = Function("e_3", real=True)(t)
u_f   = Function("u",   real=True)(t)
v_f   = Function("v",   real=True)(t)
w_f   = Function("w",   real=True)(t)
p_f   = Function("p",   real=True)(t)
q_f   = Function("q",   real=True)(t)
r_f   = Function("r",   real=True)(t)

f_subs = []
for item in ("p_n,p_e,p_d,e_0,e_1,e_2,e_3,u,v,w,p,q,r".split(",")):
    f_subs.append((eval(item), eval(f"{item}_f")))

m, g, J_xx, J_yy, J_zz, J_xy, J_xz, J_yz = symbols("m, g, J_xx, J_yy, J_zz, J_xy, J_xz, J_yz")

J = Matrix([
    [ J_xx, -J_xy, -J_xz],
    [-J_xy,  J_yy, -J_yz],
    [-J_xz, -J_yz,  J_zz]
])

V     = Matrix([u, v, w])
omega = Matrix([p, q, r])
e     = Matrix([e_0, e_1, e_2, e_3])

F_E = Matrix([f_E_x, f_E_y, f_E_z])
F_g = Matrix([f_g_x, f_g_y, f_g_z])
F_cp = Matrix([f_cp_x, f_cp_y, f_cp_z])
r_E = Matrix([r_E_x, r_E_y, r_E_z])
r_cp = Matrix([r_cp_x, r_cp_y, r_cp_z])

F = F_E + F_g + F_cp
write(F, "F")

tau = r_E.cross(F_E) + r_cp.cross(F_cp)
write(tau, "tau")

F_g = m*g*Matrix([
    2*(e_1*e_3 - e_2*e_0),
    2*(e_2*e_3 + e_1*e_0),
    e_3**2 + e_0**2 - e_1**2 - e_2**2
])

F = F_E + F_g + F_cp

# vectors
M1 = Matrix([
    [e_1**2 + e_0**2 - e_2**2 - e_3**2, 2*(e_1*e_2 - e_3*e_0), 2*(e_1*e_3 + e_2*e_0)],
    [2*(e_1*e_2 + e_3*e_0), e_2**2 + e_0**2 - e_1**2 - e_3**2, 2*(e_2*e_3 - e_1*e_0)],
    [2*(e_1*e_3 - e_2*e_0), 2*(e_2*e_3 + e_1*e_0), e_3**2 + e_0**2 - e_1**2 - e_2**2]
])

M2 = Rational(1, 2)*Matrix([
    [0, -p, -q, -r],
    [p, 0, r, -q],
    [q, -r, 0, p],
    [r, q, -p, 0]
])


P_g_dot   =  M1*V
V_dot     = -1*omega.cross(V) + 1/m*F
e_dot     =  M2*e
omega_dot = -1*J.inv()*(omega.cross(J*omega)) + J.inv()*tau

dot_x_e = Matrix([
    0, # dot u
    0, # dot v
    0, # dot w
    0, # dot e_0
    0, # dot e_1
    0, # dot e_2
    0, # dot e_3
    0, # dot p
    0, # dot q
    0  # dot r
])

f = zeros(len(dot_x_e), 1)
shift = 3
# f[0]  = P_g_dot[0]
# f[1]  = P_g_dot[1]
# f[2]  = P_g_dot[2]
f[3-shift]  = V_dot[0]
f[4-shift]  = V_dot[1]
f[5-shift]  = V_dot[2]
f[6-shift]  = e_dot[0]
f[7-shift]  = e_dot[1]
f[8-shift]  = e_dot[2]
f[9-shift]  = e_dot[3]
f[10-shift] = omega_dot[0]
f[11-shift] = omega_dot[1]
f[12-shift] = omega_dot[2]

write(f, "f")

#### for landing 
phi_e   = rad(0) # dont care, pick 0
theta_e = rad(0) # has to be this
psi_e   = rad(0) # dont care, pick 0
x_e = Matrix([
    0, # p_n
    0, # p_e
    0, # p_d
    0, # u
    0, # v
    0, # w
    Euler2Quaternion(phi_e, theta_e, psi_e)[0], # e_0
    Euler2Quaternion(phi_e, theta_e, psi_e)[1], # e_1
    Euler2Quaternion(phi_e, theta_e, psi_e)[2], # e_2
    Euler2Quaternion(phi_e, theta_e, psi_e)[3], # e_3
    0, # p
    0, # q
    0  # r
])

u_e = Matrix([
    m*g, 0, 0, # f_E_x, f_E_y, f_E_z
    0, 0, 0,   # f_cp_x, f_cp_y, f_cp_z
    0, 0, 0    # r_cp_x, r_cp_y, r_cp_z
])

# u = [f_E_x, f_E_y, f_E_z, f_cp_x, f_cp_y, f_cp_z, r_cp_x, r_cp_y, r_cp_z]

write(phi_e, "phi_e")
write(theta_e, "theta_e")
write(psi_e, "psi_e")

x_e_subs = []
for i, item in enumerate(("e_0,e_1,e_2,e_3,u,v,w,p,q,r".split(","))):
    x_e_subs.append((eval(item), x_e[i]))

u_e_subs = []
for i, item in enumerate(("f_E_x,f_E_y,f_E_z,f_cp_x,f_cp_y,f_cp_z,r_cp_x,r_cp_y,r_cp_z".split(","))):
    u_e_subs.append((eval(item), u_e[i]))

df_dx = f.jacobian([e_0,e_1,e_2,e_3,u,v,w,p,q,r])
write(df_dx, "df_dx")

df_du = f.jacobian([f_E_x, f_E_y, f_E_z, f_cp_x, f_cp_y, f_cp_z, r_cp_x, r_cp_y, r_cp_z])

# f_e = f.subs(x_e_subs + u_e_subs)
# write(f_e, "f_e")

A = df_dx.subs(x_e_subs + u_e_subs)
write(A, "A")

B = df_du.subs(x_e_subs + u_e_subs)
write(B, "B")