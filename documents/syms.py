from syms_vars import *
from sympy import Rational, rad, zeros, solve, eye, Symbol, Poly, ones
from syms_tools import write, Euler2Quaternion, writeMathModule


F = F_E + F_g + F_cp
write(F, "F")

tau = r_E.cross(F_E) + r_cp.cross(F_cp)
write(tau, "tau")


F_g = m*g*Matrix([
    2*(e_1*e_3 - e_2*e_0),
    2*(e_2*e_3 + e_1*e_0),
    e_3**2 + e_0**2 - e_1**2 - e_2**2
])

# update the force to include actual force of gravity
F = F_E + F_g + F_cp
write(F, "F_subbed")

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
omega_dot = -1*J.inv()*(omega.cross(J*omega)) + J.inv()*tau

dot_p_n  = P_g_dot[0].subs(vals)
dot_p_e  = P_g_dot[1].subs(vals)
dot_p_d  = P_g_dot[2].subs(vals)
dot_u    = V_dot[0].subs(vals)
dot_v    = V_dot[1].subs(vals)
dot_w    = V_dot[2].subs(vals)
dot_e_0  = e_dot[0].subs(vals)
dot_e_1  = e_dot[1].subs(vals)
dot_e_2  = e_dot[2].subs(vals)
dot_e_3  = e_dot[3].subs(vals)
dot_p    = omega_dot[0].subs(vals)
dot_q    = omega_dot[1].subs(vals)
dot_r    = omega_dot[2].subs(vals)

x_names = "p_n,p_e,p_d,u,v,w,e_0,e_1,e_2,e_3,p,q,r".split(",")
u_names = "f_E_x,f_E_y,f_E_z,f_cp_x,f_cp_y,f_cp_z,r_cp_x,r_cp_y,r_cp_z".split(",")

f_general = zeros(len(x_names), 1)
for i in range(len(x_names)):
    f_general[i] = eval(f"dot_{x_names[i]}")

write(f_general, "f_general")

def computeStateSpace(name, x_vars, u_vars, x_e, x_dot_e = None):
    space = "  "
    if x_dot_e is None:
        x_dot_e_vec = [0 for i in range(len(x_names))]
    else:
        # setting equilibrium dot x values to 0 if they are not provided
        for item in x_names:
            item = eval(item)
            if item not in x_dot_e:
                print(item)
                # print(2*space + "Adding:", item)
                x_dot_e[item] = 0
        
        x_dot_e_vec  = []
        for item in x_names:
            item = eval(item)
            x_dot_e_vec.append(x_dot_e[item])
    x_dot_e_vec = Matrix(x_dot_e_vec)
    print(space+"x_dot_e:")
    for i in range(len(x_names)):
        print(2*space + f"{x_names[i]} = {x_dot_e_vec[i]}")

    # getting the x symbolic variables
    x_vec = [eval(item) for item in x_vars]

    # getting the u symbolic variables
    u_vec = [eval(item) for item in u_vars]

    # the number of control variables
    n = len(x_vec)
    f = zeros(n, 1) # setting up the f in dot x = f(x,u)

    # grabbing the equations
    for i in range(n): f[i] = eval(f"dot_{x_vars[i]}")
    write(f, f"{name}_f")

    # setting equilibrium x values to default if they are not provided
    for item in x_names:
        item = eval(item)
        if item not in x_e:
            # print(2*space + "Adding:", item)
            x_e[item] = default_values_x[item]

    # making a list of values to substitute for in the equation
    x_e_subs = []
    x_e_vec  = []
    for item in x_names:
        item = eval(item)
        x_e_subs.append((item, x_e[item]))
        x_e_vec.append(x_e[item])
    write(Matrix(x_e_vec), f"{name}_x_e")

    # printing out equilibrium state
    print(space+"x_e:")
    for item in x_e_subs:
        print(2*space + f"{item[0]} = ", item[1], sep="")
    
    # filling in inputs that are not the ones being controlled with the default values
    u_subs = []
    for item in u_names:
        item = eval(item)
        if item not in u_vec:
            # print(2*space + f"{item} = ", default_values_u[item], sep="")
            u_subs.append((item, default_values_u[item]))

    # subbing the values into the general f equation
    f_e = f_general.subs(x_e_subs + u_subs)
    print(f_e)
    write(f_e, f"{name}_f_e")
    print(space+"Computing u_e")
    # finding the rest of the inputs
    u_e = solve(f_e-x_dot_e_vec, u_vec, dict=True)

    if not len(u_e):
        print("Could not find u_e")
        exit()
    else:
        u_e = u_e[0]

    # filling out the rest of the input solution with the default values
    for item in u_names:
        item = eval(item)
        if item not in u_e:
            print(2*space + "Adding:", item)
            u_e[item] = default_values_u[item]

    # making a list of values to substitute for in the equation
    u_e_subs = []
    u_e_vec  = []
    for item in u_names:
        item = eval(item)
        u_e_subs.append((item, u_e[item]))
        u_e_vec.append(u_e[item])
    write(Matrix(u_e_vec), f"{name}_u_e")

    # printing out equilibrium input
    print(space+"u_e:")
    for item in u_e_subs:
        print(2*space + f"{item[0]} = ", item[1], sep="")

    # computing jacobian of the system with respect to the state, first step to A
    df_dx = f.jacobian(x_vec)
    write(df_dx, f"{name}_df_dx")

    # computing jacobian of the system with respect to the input, first step to B
    df_du = f.jacobian(u_vec)
    write(df_du, f"{name}_df_du")

    # evaluating at equilibrium to get A
    A = df_dx.subs(x_e_subs + u_e_subs)
    write(A, f"{name}_A")

    # evaluating at equilibrium to get B
    B = df_du.subs(x_e_subs + u_e_subs)
    write(B, f"{name}_B")

    # computing the controllability matrix
    C = [(A**i)*B for i in range(0, len(x_e))]
    CC = Matrix.hstack(*C)

    # printing out the number of rows and the rank (hopefully they are equal)
    rank = CC.rank()
    
    if rank == CC.shape[0]:
        print(space+"System is controlable")
        print(2*space+"Rows:", CC.shape[0])
        print(2*space+"Rank:", rank)
    else:
        print(space+"System is NOT controlable, exiting")
        print(2*space+"Rows:", CC.shape[0])
        print(2*space+"Rank:", rank)
        exit()
    
    return A.copy(), B.copy(), CC.copy(), x_e_vec.copy(), u_e_vec.copy()

def landing():
    print("Landing")
    # for printing's sake
    space = "  "

    # variables to control
    x_vars = "p_n,p_e,p_d,u,v,w,e_0,q".split(",")
    u_vars = "f_E_x,f_E_y,f_E_z".split(",")
    n = len(x_vars)

    # write(f, "f_x_u")
    #### for landing 

    phi_e   = rad(0)  # dont care, pick 0
    theta_e = rad(90) # has to be this
    psi_e   = rad(0)  # dont care, pick 0
    print(space+"Quaternion:",Euler2Quaternion(phi_e, theta_e, psi_e))

    x_e = {
        e_0 : Euler2Quaternion(phi_e, theta_e, psi_e)[0],
        e_1 : Euler2Quaternion(phi_e, theta_e, psi_e)[1],
        e_2 : Euler2Quaternion(phi_e, theta_e, psi_e)[2],
        e_3 : Euler2Quaternion(phi_e, theta_e, psi_e)[3]
    }

    A, B, CC, x_e, u_e = computeStateSpace("landing", x_vars, u_vars, x_e)
    C = eye(n)
    print("A =", A)
    print("B =", B)
    print("x_e =", x_e)
    print("u_e =", u_e)
    
    writeMathModule("landing", A=A, B=B, C=C, x_e=x_e, u_e=u_e)


landing()

    # K = (alpha - a_A)*A_A.inv()


def flip():
    print("Descent")
    # for printing's sake
    space = "  "

    # variables to control
    x_vars = "p_n,p_e,u,v,w,e_0,e_3,q,r".split(",")
    u_vars = "f_E_x,f_E_y,f_E_z,r_cp_x,r_cp_y,r_cp_z".split(",")
    n = len(x_vars)

    # write(f, "f_x_u")
    #### for landing 

    phi_e   = rad(0)  # dont care, pick 0
    theta_e = rad(90) # has to be this
    psi_e   = rad(0)  # dont care, pick 0
    print(space+"Quaternion:",Euler2Quaternion(phi_e, theta_e, psi_e))

    x_e = {
        e_0 : Euler2Quaternion(phi_e, theta_e, psi_e)[0],
        e_1 : Euler2Quaternion(phi_e, theta_e, psi_e)[1],
        e_2 : Euler2Quaternion(phi_e, theta_e, psi_e)[2],
        e_3 : Euler2Quaternion(phi_e, theta_e, psi_e)[3]
    }
    print(x_e)


    A, B, CC, x_e, u_e = computeStateSpace("flip", x_vars, u_vars, x_e)
    C = eye(n)
    print("A =", A)
    print("B =", B)
    print("x_e =", x_e)
    print("u_e =", u_e)

    writeMathModule("flip", A=A, B=B, C=C, x_e=x_e, u_e=u_e)

flip()


def descent():
    print("Descent")
    # for printing's sake
    space = "  "

    # variables to control
    x_vars = "p_n,p_e,u,v,e_1,e_2,p,q".split(",")
    u_vars = "f_cp_x,f_cp_y,f_cp_z,r_cp_x,r_cp_y,r_cp_z".split(",")
    n = len(x_vars)

    # write(f, "f_x_u")
    #### for landing 

    phi_e   = rad(0)  # dont care, pick 0
    theta_e = rad(0) # has to be this
    psi_e   = rad(0)  # dont care, pick 0
    print(space+"Quaternion:",Euler2Quaternion(phi_e, theta_e, psi_e))

    x_e = {
        w   : 50,
        e_0 : Euler2Quaternion(phi_e, theta_e, psi_e)[0],
        e_1 : Euler2Quaternion(phi_e, theta_e, psi_e)[1],
        e_2 : Euler2Quaternion(phi_e, theta_e, psi_e)[2],
        e_3 : Euler2Quaternion(phi_e, theta_e, psi_e)[3]
    }
    print(x_e)
    # exit()
    x_dot_e = {
        p_d : 50
    }


    A, B, CC, x_e, u_e = computeStateSpace("descent", x_vars, u_vars, x_e, x_dot_e)
    C = eye(n)
    print("A =", A)
    print("B =", B)
    print("x_e =", x_e)
    print("u_e =", u_e)

    writeMathModule("descent", A=A, B=B, C=C, x_e=x_e, u_e=u_e)

descent()



"""
s = Symbol("s")
    sI_minus_A = s*eye(n) - A
    write(sI_minus_A, "sI_minus_A")

    Delta_ol = Poly(sI_minus_A.det(), s)
    a_A = Delta_ol.all_coeffs()
    print(a_A)
    if a_A[0] != 1:
        print(space+"Normalizing Delta_ol to get leading coef. of 1, current coef. is", a_A[0])
        Delta_ol = Delta_ol/a_A[0]
        a_A = Delta_ol.all_coeffs()
    
    A_A = eye(n)
    for i in range(n):
        for j in range(i+1, n):
            A_A[i,j] = a_A[-(j-i)]
    # print(A_A)

    # must have length n
    p = -1*Rational(1, 2)*ones(n, 1)
    # print(p)

    if len(p) != n:
        print("p (list of poles) does not have the right length")
        exit()
    
    Delta_cl = 1
    for i in range(n):
        Delta_cl*=(s - p[i])
    Delta_cl = Poly(Delta_cl, s)

    alpha = Delta_cl.all_coeffs()
    print(alpha)
    if alpha[0] != 1:
        print(space+"Normalizing Delta_cl to get leading coef. of 1, current coef. is", alpha[0])
        Delta_cl = Delta_cl/alpha[0]
        alpha = Delta_cl.all_coeffs()
    
    a_A   = Matrix(a_A).T
    alpha = Matrix(alpha).T
    print(alpha, a_A)

"""