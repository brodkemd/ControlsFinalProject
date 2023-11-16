from syms_vars import *
from sympy import Rational, rad, zeros, solve
from syms_tools import write, Euler2Quaternion

F = F_E + F_g + F_cp
write(F, "F")

tau = r_E.cross(F_E) + r_cp.cross(F_cp)
write(tau, "tau")

F_g = m*g*Matrix([
    2*(e_1*e_3 - e_2*e_0),
    2*(e_2*e_3 + e_1*e_0),
    e_3**2 + e_0**2 - e_1**2 - e_2**2
])

# F_E = Matrix([
#     T*cos(Phi),
#     T*sin(Phi)*cos(Theta),
#     T*sin(Phi)*cos(Theta)
# ])

# update the force to include actual force of gravity
F = F_E + F_g + F_cp
write(F, "F_subbed")

# tau = r_E.cross(F_E) + r_cp.cross(F_cp)
# write(tau, "tau_subbed")

# V to V_g
M1 = Matrix([
    [e_1**2 + e_0**2 - e_2**2 - e_3**2, 2*(e_1*e_2 - e_3*e_0), 2*(e_1*e_3 + e_2*e_0)],
    [2*(e_1*e_2 + e_3*e_0), e_2**2 + e_0**2 - e_1**2 - e_3**2, 2*(e_2*e_3 - e_1*e_0)],
    [2*(e_1*e_3 - e_2*e_0), 2*(e_2*e_3 + e_1*e_0), e_3**2 + e_0**2 - e_1**2 - e_2**2]
])
write(M1, "M1")

# e to dot e
M2 = Rational(1, 2)*Matrix([
    [0, -p, -q, -r],
    [p, 0, r, -q],
    [q, -r, 0, p],
    [r, q, -p, 0]
])
write(2*M2, "M2")

# dynamics equation
P_g_dot   =  M1*V
V_dot     = -1*omega.cross(V) + 1/m*F
e_dot     =  M2*e
omega_dot = -1*J_val.inv()*(omega.cross(J*omega)) + J_val.inv()*tau

dot_p_n  = P_g_dot[0]
dot_p_e  = P_g_dot[1]
dot_p_d  = P_g_dot[2]
dot_u    = V_dot[0]
dot_v    = V_dot[1]
dot_w    = V_dot[2]
dot_e_0  = e_dot[0]
dot_e_1  = e_dot[1]
dot_e_2  = e_dot[2]
dot_e_3  = e_dot[3]
dot_p    = omega_dot[0]
dot_q    = omega_dot[1]
dot_r    = omega_dot[2]

x_names = "p_n,p_e,p_d,e_0,e_1,e_2,e_3,u,v,w,p,q,r".split(",")
u_names = "f_E_x,f_E_y,f_E_z,f_cp_x,f_cp_y,f_cp_z,r_cp_x,r_cp_y,r_cp_z".split(",")

f_general = zeros(len(x_names), 1)
for i in range(len(x_names)):
    f_general[i] = eval(f"dot_{x_names[i]}")

f_general = 

write(f_general, "f_general")

exit()
def landing():
    space = "  "
    ############ landing state space
    # f in dot x = f(x,u)
    x_vars = "p_n,p_e,p_d,e_0,e_2,u,v,w,p,r".split(",")
    u_vars = "f_E_x,f_E_y,f_E_z,f_cp_x,f_cp_y,f_cp_z,r_cp_x,r_cp_y,r_cp_z".split(",")

    x_vec = [eval(item) for item in x_vars]
    u_vec = [eval(item) for item in u_vars]
    n = len(x_vec)
    f = zeros(n, 1)

    vals = [(f_cp_x, 0), (f_cp_y, 0), (f_cp_z, 0), (r_cp_x, 0), (r_cp_y,0), (r_cp_z,0)]

    for i in range(n): f[i] = eval(f"dot_{x_vars[i]}")

    write(f, "f_x_u")
    #### for landing 

    phi_e   = rad(0) # dont care, pick 0
    theta_e = rad(90) # has to be this
    psi_e   = rad(0) # dont care, pick 0
    print(space+"Quaternion:",Euler2Quaternion(phi_e, theta_e, psi_e))
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

    # u_e = Matrix([
    #     m*g, 0, 0, # f_E_x, f_E_y, f_E_z
    #     0, 0, 0,   # f_cp_x, f_cp_y, f_cp_z
    #     0, 0, 0    # r_cp_x, r_cp_y, r_cp_z
    # ])

    x_e_subs = []
    for i, item in enumerate(x_names):
        x_e_subs.append((eval(item), x_e[i]))

    f_e = f.subs(x_e_subs)
    print(f_e)
    print(solve(f_e, u_vec, dict=True))
    exit()

    u_e_subs = []
    for i, item in enumerate(u_names):
        u_e_subs.append((eval(item), u_e[i]))

    df_dx = f.jacobian(x_vec)
    write(df_dx, "df_dx")

    df_du = f.jacobian(u_vec)
    write(df_du, "df_du")

    A = df_dx.subs(x_e_subs + u_e_subs)
    write(A, "A")

    B = df_du.subs(x_e_subs + u_e_subs)
    write(B, "B")

    C = []
    for i in range(0, len(x_e)):
        C.append((A**i)*B)

    CC = Matrix.hstack(*C)
    print(space+"Rows:", CC.shape[0])
    print(space+"Rank:", CC.rank())


def roll():
    space = "  "
    print("Roll:")
    ############ landing state space
    # f in dot x = f(x,u)
    x_vars = "e_1,e_3,q".split(",")
    u_vars = "f_cp_x,f_cp_y,f_cp_z,r_cp_x,r_cp_y,r_cp_z".split(",")

    x_vec = [eval(item) for item in x_vars]
    u_vec = [eval(item) for item in u_vars]
    n = len(x_vec)
    f = zeros(n, 1)

    vals = [(f_E_x, 0), (f_E_y, 0), (f_E_z, 0)]

    for i in range(n): f[i] = eval(f"dot_{x_vars[i]}")

    write(f, "f_x_u")
    #### for landing

    phi_e   = rad(0) # has to be, this is omitted
    theta_e = rad(0) # has to be
    psi_e   = rad(0) # has to be
    print(space+"Quaternion:",Euler2Quaternion(phi_e, theta_e, psi_e))
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

    # u_e = Matrix([
    #     0,    0, 0, # f_E_x, f_E_y, f_E_z
    #     0, -m*g, 0, # f_cp_x, f_cp_y, f_cp_z
    #     0,    0, 0  # r_cp_x, r_cp_y, r_cp_z
    # ])

    x_e_subs = []
    for i, item in enumerate(x_names):
        x_e_subs.append((eval(item), x_e[i]))

    f_e = f.subs(x_e_subs)
    sol = solve(f_e, u_vec, dict=True)[0]
    for param in sol:
        print(2*space+f"{param} = {sol[param]}")
        # exec(f"{param} = {sol[param]}")
    
    # # print(solve(f_e, u_vec, dict=True))

    #print(f_e)
    # exit()

    u_e_subs = []
    for i, item in enumerate(u_names):
        u_e_subs.append((eval(item), u_e[i]))

    df_dx = f.jacobian(x_vec)
    write(df_dx, "df_dx")

    df_du = f.jacobian(u_vec)
    write(df_du, "df_du")

    A = df_dx.subs(x_e_subs + u_e_subs)
    write(A, "A")

    B = df_du.subs(x_e_subs + u_e_subs)
    write(B, "B")

    C = []
    for i in range(0, len(x_e)):
        C.append((A**i)*B)

    CC = Matrix.hstack(*C)
    print(space+"Rows:", CC.shape[0])
    print(space+"Rank:", CC.rank())




def descent():
    space = "  "
    print("Descent:")
    ############ landing state space
    # f in dot x = f(x,u)
    x_vars = "p_n,p_e,p_d,e_0,e_2,u,v,w,p,r".split(",")
    u_vars = "f_cp_x,f_cp_y,f_cp_z,r_cp_x,r_cp_y,r_cp_z".split(",")

    x_vec = [eval(item) for item in x_vars]
    u_vec = [eval(item) for item in u_vars]
    n = len(x_vec)
    f = zeros(n, 1)

    vals = [(f_E_x, 0), (f_E_y, 0), (f_E_z, 0)]

    for i in range(n): f[i] = eval(f"dot_{x_vars[i]}")

    write(f, "f_x_u")
    #### for landing

    phi_e   = rad(0) # has to be, this is omitted
    theta_e = rad(0) # has to be
    psi_e   = rad(0) # has to be
    print(space+"Quaternion:",Euler2Quaternion(phi_e, theta_e, psi_e))
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

    # u_e = Matrix([
    #     0,    0, 0, # f_E_x, f_E_y, f_E_z
    #     0, -m*g, 0, # f_cp_x, f_cp_y, f_cp_z
    #     0,    0, 0  # r_cp_x, r_cp_y, r_cp_z
    # ])

    x_e_subs = []
    for i, item in enumerate(x_names):
        x_e_subs.append((eval(item), x_e[i]))

    f_e = f.subs(x_e_subs)
    sol = solve(f_e, u_vec, dict=True)[0]
    for param in sol:
        print(2*space+f"{param} = {sol[param]}")
        # exec(f"{param} = {sol[param]}")
    
    # # print(solve(f_e, u_vec, dict=True))

    #print(f_e)
    exit()

    u_e_subs = []
    for i, item in enumerate(u_names):
        u_e_subs.append((eval(item), u_e[i]))

    df_dx = f.jacobian(x_vec)
    write(df_dx, "df_dx")

    df_du = f.jacobian(u_vec)
    write(df_du, "df_du")

    A = df_dx.subs(x_e_subs + u_e_subs)
    write(A, "A")

    B = df_du.subs(x_e_subs + u_e_subs)
    write(B, "B")

    C = []
    for i in range(0, len(x_e)):
        C.append((A**i)*B)

    CC = Matrix.hstack(*C)
    print(space+"Rows:", CC.shape[0])
    print(space+"Rank:", CC.rank())

# roll()
# landing()
descent()
