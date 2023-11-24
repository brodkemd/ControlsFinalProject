from sympy import symbols, Matrix

t, m, g, J_xx, J_yy, J_zz, J_xy, J_xz, J_yz, R                          = symbols("t, m, g, J_xx, J_yy, J_zz, J_xy, J_xz, J_yz, R")
f_E_x, f_E_y, f_E_z, f_g_x, f_g_y, f_g_z, r_E_x, r_E_y, r_E_z           = symbols("f_E_x, f_E_y, f_E_z, f_g_x, f_g_y, f_g_z, r_E_x, r_E_y, r_E_z")
p_n,p_e,p_d,e_0,e_1,e_2,e_3,u,v,w,p,q,r                                 = symbols("p_n,p_e,p_d,e_0,e_1,e_2,e_3,u,v,w,p,q,r")

f_cp_port_canard_x,f_cp_port_canard_y,f_cp_port_canard_z                = symbols("f_cp_port_canard_x,f_cp_port_canard_y,f_cp_port_canard_z")
f_cp_starboard_canard_x,f_cp_starboard_canard_y,f_cp_starboard_canard_z = symbols("f_cp_starboard_canard_x,f_cp_starboard_canard_y,f_cp_starboard_canard_z")

f_cp_port_fin_x,f_cp_port_fin_y,f_cp_port_fin_z                         = symbols("f_cp_port_fin_x,f_cp_port_fin_y,f_cp_port_fin_z")
f_cp_starboard_fin_x,f_cp_starboard_fin_y,f_cp_starboard_fin_z          = symbols("f_cp_starboard_fin_x,f_cp_starboard_fin_y,f_cp_starboard_fin_z")

r_cp_port_canard_x,r_cp_port_canard_y,r_cp_port_canard_z                = symbols("r_cp_port_canard_x,r_cp_port_canard_y,r_cp_port_canard_z")
r_cp_starboard_canard_x,r_cp_starboard_canard_y,r_cp_starboard_canard_z = symbols("r_cp_starboard_canard_x,r_cp_starboard_canard_y,r_cp_starboard_canard_z")

r_cp_port_fin_x,r_cp_port_fin_y,r_cp_port_fin_z                         = symbols("r_cp_port_fin_x,r_cp_port_fin_y,r_cp_port_fin_z")
r_cp_starboard_fin_x,r_cp_starboard_fin_y,r_cp_starboard_fin_z          = symbols("r_cp_starboard_fin_x,r_cp_starboard_fin_y,r_cp_starboard_fin_z")

J = Matrix([
    [ J_xx, -J_xy, -J_xz],
    [-J_xy,  J_yy, -J_yz],
    [-J_xz, -J_yz,  J_zz]
])

J_xy_val = 0
J_xz_val = 0
J_yz_val = 0
r_E_y_val = 0
r_E_z_val = 0


vals = [(J_xy, J_xy_val), (J_xz, J_xz_val), (J_yz, J_yz_val), (r_E_y, r_E_y_val), (r_E_z, r_E_z_val)]

V     = Matrix([u, v, w])
omega = Matrix([p, q, r])
e     = Matrix([e_0, e_1, e_2, e_3])

F_E  = Matrix([ f_E_x,  f_E_y,  f_E_z])
F_g  = Matrix([ f_g_x,  f_g_y,  f_g_z])
r_E  = Matrix([ r_E_x,  r_E_y,  r_E_z])

default_values_x = {
    p_n : 0,
    p_e : 0,
    p_d : 0,
    u   : 0,
    v   : 0,
    w   : 0,
    e_0 : 1,
    e_1 : 0,
    e_2 : 0,
    e_3 : 0,
    p   : 0,
    q   : 0,
    r   : 0
}


# F_cp_body = Matrix([f_cp_body_x,f_cp_body_y,f_cp_body_z])

F_cp_port_canard      = Matrix([      f_cp_port_canard_x,      f_cp_port_canard_y,      f_cp_port_canard_z])
F_cp_starboard_canard = Matrix([ f_cp_starboard_canard_x, f_cp_starboard_canard_y, f_cp_starboard_canard_z])
F_cp_port_fin         = Matrix([         f_cp_port_fin_x,         f_cp_port_fin_y,         f_cp_port_fin_z])
F_cp_starboard_fin    = Matrix([    f_cp_starboard_fin_x,    f_cp_starboard_fin_y,    f_cp_starboard_fin_z])

F_cp_canard = F_cp_port_canard + F_cp_starboard_canard
F_cp_fin    = F_cp_port_fin + F_cp_starboard_fin
F_cp        = F_cp_fin + F_cp_canard

# r_cp_body = Matrix([r_cp_body_x,r_cp_body_y,r_cp_body_z])
r_cp_port_canard        = Matrix([      r_cp_port_canard_x,      r_cp_port_canard_y,      r_cp_port_canard_z])
r_cp_starboard_canard   = Matrix([ r_cp_starboard_canard_x, r_cp_starboard_canard_y, r_cp_starboard_canard_z])
r_cp_port_fin           = Matrix([         r_cp_port_fin_x,         r_cp_port_fin_y,         r_cp_port_fin_z])
r_cp_starboard_fin      = Matrix([    r_cp_starboard_fin_x,    r_cp_starboard_fin_y,    r_cp_starboard_fin_z])


tau_cp_port_canard      = r_cp_port_canard.cross(F_cp_port_canard)
tau_cp_starboard_canard = r_cp_starboard_canard.cross(F_cp_starboard_canard)
tau_cp_port_fin         = r_cp_port_fin.cross(F_cp_port_fin)
tau_cp_starboard_fin    = r_cp_starboard_fin.cross(F_cp_starboard_fin)


tau_control = tau_cp_port_canard + tau_cp_port_fin + tau_cp_starboard_canard + tau_cp_starboard_fin

# f_cp_port_canard_x_val      = 0
# f_cp_starboard_canard_x_val = 0
# f_cp_port_fin_x_val         = 0
# f_cp_starboard_fin_x_val    = 0

r_cp_port_canard_x_val      =  1
r_cp_port_canard_y_val      = -1
r_cp_port_canard_z_val      =  0
r_cp_starboard_canard_x_val =  1
r_cp_starboard_canard_y_val =  1
r_cp_starboard_canard_z_val =  0
r_cp_port_fin_x_val         = -1
r_cp_port_fin_y_val         = -1
r_cp_port_fin_z_val         =  0
r_cp_starboard_fin_x_val    = -1
r_cp_starboard_fin_y_val    =  1
r_cp_starboard_fin_z_val    =  0

names = "f_cp_port_canard_x,f_cp_starboard_canard_x,f_cp_port_fin_x,f_cp_starboard_fin_x,r_cp_body_x,r_cp_body_y,r_cp_body_z,r_cp_port_canard_x,r_cp_port_canard_y,r_cp_port_canard_z,r_cp_starboard_canard_x,r_cp_starboard_canard_y,r_cp_starboard_canard_z,r_cp_port_fin_x,r_cp_port_fin_y,r_cp_port_fin_z,r_cp_starboard_fin_x,r_cp_starboard_fin_y,r_cp_starboard_fin_z".split(",")

for name in names:
    if f"{name}_val" in vars():
        vals.append((eval(name), eval(f"{name}_val")))

default_values_u = {
    f_E_x                   : 0,
    f_E_y                   : 0,
    f_E_z                   : 0,
    f_cp_port_canard_x      : 0,
    f_cp_port_canard_y      : 0,
    f_cp_port_canard_z      : 0,
    f_cp_starboard_canard_x : 0,
    f_cp_starboard_canard_y : 0,
    f_cp_starboard_canard_z : 0,
    f_cp_port_fin_x         : 0,
    f_cp_port_fin_y         : 0,
    f_cp_port_fin_z         : 0,
    f_cp_starboard_fin_x    : 0,
    f_cp_starboard_fin_y    : 0,
    f_cp_starboard_fin_z    : 0
}

x_names = "p_n,p_e,p_d,u,v,w,e_0,e_1,e_2,e_3,p,q,r".split(",")
u_names = "f_E_x,f_E_y,f_E_z,f_cp_port_canard_x,f_cp_port_canard_y,f_cp_port_canard_z,f_cp_starboard_canard_x,f_cp_starboard_canard_y,f_cp_starboard_canard_z,f_cp_port_fin_x,f_cp_port_fin_y,f_cp_port_fin_z,f_cp_starboard_fin_x,f_cp_starboard_fin_y,f_cp_starboard_fin_z".split(",")