import sympy as sp
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
        #s = MyLatexPrinter().doprint(s)
        s = sp.latex(s)
        print

    _file = os.path.join(cwd, _file)
    with open(_file, "w") as file:
        file.write(s)

def diff(f):
    return sp.diff(f, t)

t = sp.symbols("t")
f_E_x, f_E_y, f_E_z, f_g_x, f_g_y, f_g_z, f_cp_x, f_cp_y, f_cp_z, r_E_x, r_E_y, r_E_z, r_cp_x, r_cp_y, r_cp_z = sp.symbols("f_E_x, f_E_y, f_E_z, f_g_x, f_g_y, f_g_z, f_cp_x, f_cp_y, f_cp_z, r_E_x, r_E_y, r_E_z, r_cp_x, r_cp_y, r_cp_z")
e_0,e_1,e_2,e_3,u,v,w,p,q,r = sp.symbols("e_0,e_1,e_2,e_3,u,v,w,p,q,r")
e_0_f = sp.Function("e_0", real=True)(t)
e_1_f = sp.Function("e_1", real=True)(t)
e_2_f = sp.Function("e_2", real=True)(t)
e_3_f = sp.Function("e_3", real=True)(t)
u_f = sp.Function("u", real=True)(t)
v_f = sp.Function("v", real=True)(t)
w_f = sp.Function("w", real=True)(t)
p_f = sp.Function("p", real=True)(t)
q_f = sp.Function("q", real=True)(t)
r_f = sp.Function("r", real=True)(t)

f_subs = []
for item in ("e_0,e_1,e_2,e_3,u,v,w,p,q,r".split(",")):
    f_subs.append((eval(item), eval(f"{item}_f")))

Gamma_1, Gamma_2, Gamma_3, Gamma_4, Gamma_5, Gamma_6, Gamma_7, Gamma_8, m, J_xx, J_yy, J_zz, J_xz = sp.symbols("Gamma_1, Gamma_2, Gamma_3, Gamma_4, Gamma_5, Gamma_6, Gamma_7, Gamma_8, m, J_xx, J_yy, J_zz, J_xz")

J = Matrix([
    [ J_xx,    0, -J_xz],
    [    0, J_yy,     0],
    [-J_xz,    0,  J_zz]
])

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
x1 = Matrix([u, v, w])
x2 = Matrix([1])
x3 = Matrix([e_0, e_1, e_2, e_3])
x4 = Matrix([1])

M1 = Matrix([
    [e_1**2 + e_0**2 - e_2**2 - e_3**2, 2*(e_1*e_2 - e_3*e_0), 2*(e_1*e_3 + e_2*e_0)],
    [2*(e_1*e_2 + e_3*e_0), e_2**2 + e_0**2 - e_1**2 - e_3**2, 2*(e_2*e_3 - e_1*e_0)],
    [2*(e_1*e_3 - e_2*e_0), 2*(e_2*e_3 + e_1*e_0), e_3**2 + e_0**2 - e_1**2 - e_2**2]
])
M2 = Matrix([
    [r*v - q*w],
    [p*w - r*u],
    [q*u - p*v]
])
M3 = sp.Rational(1, 2)*Matrix([
    [0, -p, -q, -r],
    [p, 0, r, -q],
    [q, -r, 0, p],
    [r, q, -p, 0]
])
x_acc = Matrix([p, q, r])
M4 =  Matrix([
    [Gamma_1*p*q - Gamma_2*q*r],
    [Gamma_5*p*r - Gamma_6*(p**2 - r**2)],
    [Gamma_7*p*q - Gamma_1*q*r]
])

for item in ["M1", "M2", "M3", "M4"]:
    write(eval(item), item)



for item in ["x1", "x2", "x3", "x4"]:
    write(eval(item), item)

# tau = [l, m, n]

B1_F = sp.zeros(3, 3)
B1_tau = sp.zeros(3, 3)

B2_F = 1/m*sp.eye(3)
B2_tau = sp.zeros(3, 3)

B3_F = sp.zeros(4, 3)
B3_tau = sp.zeros(4, 3)

B4_F = sp.zeros(3, 3)
B4_tau = Matrix([
    [Gamma_3, 0, Gamma_4],
    [0, 1/J_yy, 0],
    [Gamma_4, 0, Gamma_8]
])

for item in ["B1", "B2", "B3", "B4"]:
    for suf in ["F", "tau"]:
        temp = f"{item}_{suf}"
        write(eval(temp), temp)

# write(J.inv(), "J_inv")
# x_acc = Matrix([p, q, r])

# M4_acc = -1*J.inv()*x_acc.cross(J*x_acc)
# write(M4_acc, "M4_acc")

# B4 = Matrix([
#     [Gamma_3*l + Gamma_4*n],
#     [1/J_yy*m],
#     [Gamma_4*l + Gamma_8*n]
# ])
# state = [p_n p_e p_d  u  v  w  e_0  e_1  e_2  e_3  p  q  r]
#b1 = M1*x1 + B1_F*F + B1_tau*tau
#b2 = M2*x2 + B2_F*F + B2_tau*tau
#b3 = M3*x3 + B3_F*F + B3_tau*tau
#b4 = M4*x4 + B4_F*F + B4_tau*tau