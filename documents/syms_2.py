from sympy import *
import os
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

def Diff(f):
    return diff(f, t)

t = symbols("t")
f_E_x, f_E_y, f_E_z, f_g_x, f_g_y, f_g_z, f_cp_x, f_cp_y, f_cp_z, r_E_x, r_E_y, r_E_z, r_cp_x, r_cp_y, r_cp_z = symbols("f_E_x, f_E_y, f_E_z, f_g_x, f_g_y, f_g_z, f_cp_x, f_cp_y, f_cp_z, r_E_x, r_E_y, r_E_z, r_cp_x, r_cp_y, r_cp_z")
e_0,e_1,e_2,e_3,u,v,w,p,q,r = symbols("e_0,e_1,e_2,e_3,u,v,w,p,q,r")
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
for item in ("e_0,e_1,e_2,e_3,u,v,w,p,q,r".split(",")):
    f_subs.append((eval(item), eval(f"{item}_f")))

m, J_xx, J_yy, J_zz, J_xy, J_xz, J_yz = symbols("m, J_xx, J_yy, J_zz, J_xy, J_xz, J_yz")

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