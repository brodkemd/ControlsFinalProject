from sympy import symbols, Matrix, Function

t, m, g, J_xx, J_yy, J_zz, J_xy, J_xz, J_yz = symbols("t, m, g, J_xx, J_yy, J_zz, J_xy, J_xz, J_yz")
f_E_x, f_E_y, f_E_z, f_g_x, f_g_y, f_g_z, f_cp_x, f_cp_y, f_cp_z, r_E_x, r_E_y, r_E_z, r_cp_x, r_cp_y, r_cp_z = symbols("f_E_x, f_E_y, f_E_z, f_g_x, f_g_y, f_g_z, f_cp_x, f_cp_y, f_cp_z, r_E_x, r_E_y, r_E_z, r_cp_x, r_cp_y, r_cp_z")
p_n,p_e,p_d,e_0,e_1,e_2,e_3,u,v,w,p,q,r = symbols("p_n,p_e,p_d,e_0,e_1,e_2,e_3,u,v,w,p,q,r")

Phi, Theta, T = symbols("Phi, Theta, T")

J = Matrix([
    [ J_xx, -J_xy, -J_xz],
    [-J_xy,  J_yy, -J_yz],
    [-J_xz, -J_yz,  J_zz]
])

J_xy_val = 0
J_xz_val = 0
J_yz_val = 0
r_E_y = 0
r_E_z = 0


vals = [(J_xy, J_xy_val), (J_xz, J_xz_val), (J_yz, J_yz_val)]
J_val = J.subs(vals)

V     = Matrix([u, v, w])
omega = Matrix([p, q, r])
e     = Matrix([e_0, e_1, e_2, e_3])

F_E  = Matrix([ f_E_x,  f_E_y,  f_E_z])
F_g  = Matrix([ f_g_x,  f_g_y,  f_g_z])
F_cp = Matrix([f_cp_x, f_cp_y, f_cp_z])
r_E  = Matrix([ r_E_x,  r_E_y,  r_E_z])
r_cp = Matrix([r_cp_x, r_cp_y, r_cp_z])


p_n_f = Function("p_n", real=True)(t)
p_e_f = Function("p_e", real=True)(t)
p_d_f = Function("p_d", real=True)(t)
u_f   = Function("u",   real=True)(t)
v_f   = Function("v",   real=True)(t)
w_f   = Function("w",   real=True)(t)
e_0_f = Function("e_0", real=True)(t)
e_1_f = Function("e_1", real=True)(t)
e_2_f = Function("e_2", real=True)(t)
e_3_f = Function("e_3", real=True)(t)
p_f   = Function("p",   real=True)(t)
q_f   = Function("q",   real=True)(t)
r_f   = Function("r",   real=True)(t)

f_subs = []
for item in ("p_n,p_e,p_d,e_0,e_1,e_2,e_3,u,v,w,p,q,r".split(",")):
    f_subs.append((eval(item), eval(f"{item}_f")))