from syms_vars import *
from sympy import Rational, rad, zeros, solve, sqrt
from syms_tools import write, Euler2Quaternion, writeMathModule
import sys
sys.path.append(".")
import tools.Logging as log
import time


tau = (r_E.cross(F_E) + tau_control).subs(vals)
write(tau, "tau")

F_g = m*g*Matrix([
    2*(e_1*e_3 - e_2*e_0),
    2*(e_2*e_3 + e_1*e_0),
    e_3**2 + e_0**2 - e_1**2 - e_2**2
])

# update the force to include actual force of gravity
F = (F_E + F_g + F_cp).subs(vals)
write(F, "F")

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


f_general = zeros(len(x_names), 1)
for i in range(len(x_names)):
    f_general[i] = eval(f"dot_{x_names[i]}")

write(f_general, "f_general")


# control calculation here
compute_descent = 1
compute_flip    = 1
compute_landing = 1


def descentCP():
    log.info("DescentCP")

    # variables to control
    x_vars = "e_1,e_2,e_3,p,q,r".split(",")
    u_vars = "f_cp_port_canard_y,f_cp_port_canard_z,f_cp_starboard_canard_x,f_cp_starboard_canard_y,f_cp_starboard_canard_z,f_cp_port_fin_x,f_cp_port_fin_y,f_cp_port_fin_z,f_cp_starboard_fin_x,f_cp_starboard_fin_y,f_cp_starboard_fin_z".split(",")

    # reference_params = "p_n,p_e,u,v,w,e_0,e_1,e_2,e_3,p,q,r".split(",") # must be same length as u_vars
    reference_params = "u,v,w,e_0,e_1,e_2,e_3,p,q,r".split(",") # must be same length as u_vars

    phi_e   = rad(0) # has to be this
    theta_e = rad(0) # has to be this
    psi_e   = rad(0) # dont care, pick 0
    log.info(1,"Quaternion:", Euler2Quaternion(phi_e, theta_e, psi_e))

    x_e = {
        #w   : 50,
        e_0 : Euler2Quaternion(phi_e, theta_e, psi_e)[0],
        e_1 : Euler2Quaternion(phi_e, theta_e, psi_e)[1],
        e_2 : Euler2Quaternion(phi_e, theta_e, psi_e)[2],
        e_3 : Euler2Quaternion(phi_e, theta_e, psi_e)[3]
    }
    fraction_gravity_acceleration = 2
    u_e = {
        # f_cp_port_canard_z      : (-m*g/4)/fraction_gravity_acceleration,
        # f_cp_starboard_canard_z : (-m*g/4)/fraction_gravity_acceleration,
        # f_cp_port_fin_z         : (-m*g/4)/fraction_gravity_acceleration,
        # f_cp_starboard_fin_z    : (-m*g/4)/fraction_gravity_acceleration
    }

    if fraction_gravity_acceleration != 1:
        x_dot_e = {
            #p_d : 50,
            #w   :  g/fraction_gravity_acceleration
        }
    else:
        x_dot_e = {
            # p_d : 50
        }

    export("descentCP", x_vars, u_vars, reference_params, x_e, x_dot_e, u_e) # computes the state space and saves it to a file


def flipCP():
    print("flipCP")

    # variables to control
    x_vars = "p_n,p_e,p_d,u,v,w,e_1,e_2,e_3,p,q,r".split(",")
    u_vars = "f_E_x,f_E_y,f_E_z,f_cp_port_canard_y,f_cp_port_canard_z,f_cp_starboard_canard_y,f_cp_starboard_canard_z".split(",")
    n = len(x_vars)

    reference_params:list = "p_n,p_e,u,v,e_1,e_2,e_3".split(",") # must be same length as u_vars

    phi_e   = rad(0)  # dont care, pick 0
    theta_e = rad(90) # has to be this
    psi_e   = rad(0)  # dont care, pick 0
    log.info(1,"Quaternion:",Euler2Quaternion(phi_e, theta_e, psi_e))

    x_e = {
        e_0 : Euler2Quaternion(phi_e, theta_e, psi_e)[0],
        e_1 : Euler2Quaternion(phi_e, theta_e, psi_e)[1],
        e_2 : Euler2Quaternion(phi_e, theta_e, psi_e)[2],
        e_3 : Euler2Quaternion(phi_e, theta_e, psi_e)[3]
    }

    u_e = {
        f_E_x : m*g
    }

    export("flipCP", x_vars, u_vars, reference_params, x_e, u_e=u_e) # computes the state space and saves it to a file

def landing():
    print("Landing")

    # variables to control
    x_vars = "p_n,p_e,p_d,u,v,w,e_0,e_3,q,r".split(",") # order matters
    u_vars = "f_E_x,f_E_y,f_E_z".split(",")
    n      = len(x_vars)

    reference_params:list = "p_n,p_e,p_d".split(",")

    #### for landing 
    phi_e   = rad(0)  # dont care, pick 0
    theta_e = rad(90) # has to be this
    psi_e   = rad(0)  # dont care, pick 0
    log.info(1,"Quaternion:",Euler2Quaternion(phi_e, theta_e, psi_e))

    x_e = {
        e_0 : Euler2Quaternion(phi_e, theta_e, psi_e)[0],
        e_1 : Euler2Quaternion(phi_e, theta_e, psi_e)[1],
        e_2 : Euler2Quaternion(phi_e, theta_e, psi_e)[2],
        e_3 : Euler2Quaternion(phi_e, theta_e, psi_e)[3]
    }

    export("landing", x_vars, u_vars, reference_params, x_e) # computes the state space and saves it to a file



##### Begin computational functions

def computeStateSpace(name, x_vars, u_vars, x_e, x_dot_e = None, u_e = None):
    # sets up the dot x at equilibrium vector
    if x_dot_e is None:
        x_dot_e_vec = [0 for i in range(len(x_names))]
    else:
        # setting equilibrium dot x values to 0 if they are not provided
        for item in x_names:
            item = eval(item)
            if item not in x_dot_e:
                x_dot_e[item] = 0
        
        x_dot_e_vec  = []
        for item in x_names:
            item = eval(item)
            x_dot_e_vec.append(x_dot_e[item])
    
    x_dot_e_vec = Matrix(x_dot_e_vec)
    log.info(1, "x_dot_e:")
    for i in range(len(x_names)):
        log.info(2, f"{x_names[i]} = {x_dot_e_vec[i]}")

    # getting the x symbolic variables
    x_vec = [eval(item) for item in x_vars]

    # getting the u symbolic variables
    u_vec = [eval(item) for item in u_vars]

    # the number of control variables
    n = len(x_vec)
    f = zeros(n, 1) # setting up the f in dot x = f(x,u)

    # grabbing the equations that pertain to these parameters
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
    log.info(1,"x_e:")
    for item in x_e_subs:
        log.info(2, f"{item[0]} = " + str(item[1]))
    
    
    if u_e is None:
        # filling in inputs that are not the ones being controlled with the default values
        u_subs = []
        for item in u_names:
            item = eval(item)
            if item not in u_vec:
                # print(2*space + f"{item} = ", default_values_u[item], sep="")
                u_subs.append((item, default_values_u[item]))

        # subbing the values into the general f equation
        f_e = f_general.subs(x_e_subs + u_subs)
        log.info(2, "f_e equation:")
        for i in range(len(f_e)):
            log.info(2, x_dot_e_vec[i], "=", f_e[i])

        write(f_e, f"{name}_f_e")
        log.info(1, "Computing u_e")

        # finding the rest of the inputs
        u_e = solve(f_e-x_dot_e_vec, u_vec, dict=True)

        if not len(u_e): log.error("Could not find u_e")
        else:            u_e = u_e[0]

        # filling out the rest of the input solution with the default values
        for item in u_names:
            item = eval(item)
            if item not in u_e:
                u_e[item] = default_values_u[item]

    else:
        for item in u_names:
            item = eval(item)
            if item not in u_e:
                u_e[item] = default_values_u[item]

        u_check_subs = []
        for item in u_names:
            item = eval(item)
            u_check_subs.append((item, u_e[item]))
        f_e = f_general.subs(x_e_subs + u_check_subs)
        if f_e == x_dot_e_vec:
            log.info(1,"Provided u_e satisfies constraints")
        else:
            log.info(1, "Provided u_e DOES NOT satisfy constraints")
            log.info(2, "f_e equation:")
            for i in range(len(f_e)):
                log.info(2, x_dot_e_vec[i], "=", f_e[i])
            exit()
        
    # making a list of values to substitute for in the equation
    u_e_subs = []
    u_e_vec  = []
    for item in u_names:
        item = eval(item)
        u_e_subs.append((item, u_e[item]))
        u_e_vec.append(u_e[item])
    write(Matrix(u_e_vec), f"{name}_u_e")

    # printing out equilibrium input
    log.info(2, "u_e:")
    for item in u_e_subs:
        log.info(3, f"{item[0]} = " + str(item[1]))

    log.info(1,"Computing x jacobian:", end=" ")
    df_dx = f.jacobian(x_vec)
    write(df_dx, f"{name}_df_dx")
    log.info("done")

    # computing jacobian of the system with respect to the input, first step to B
    log.info(1, "Computing u jacobian:", end=" ")
    df_du = f.jacobian(u_vec)
    write(df_du, f"{name}_df_du")
    log.info("done")

    # evaluating at equilibrium to get A
    A = df_dx.subs(x_e_subs + u_e_subs)
    write(A, f"{name}_A")

    # evaluating at equilibrium to get B
    B = df_du.subs(x_e_subs + u_e_subs)
    write(B, f"{name}_B")

    log.info(1, "Computing Controllability Matrix:", end=" ")
    # computing the controllability matrix
    C = [(A**i)*B for i in range(n)]
    CC = Matrix.hstack(*C)
    log.info("done")

    # printing out the number of rows and the rank (hopefully they are equal)
    log.info("Computing Controllability Matrix Rank:")
    time.sleep(0.5)
    rank = CC.rank()
    log.info(2, "done")
    
    if rank == CC.shape[0]:
        log.info(1,"System is controlable")
        log.info(2, "Rows:", CC.shape[0])
        log.info(2, "Rank:", rank)
    else:
        log.info(1, "System is NOT controlable, exiting")
        log.info(2, "Rows:", CC.shape[0])
        log.info(2, "Rank:", rank)
        exit()
    
    return A.copy(), B.copy(), CC.copy(), Matrix(x_e_vec).copy(), Matrix(u_e_vec).copy()


def computeTransforms(x_vars, u_vars, reference_params, x_e):
    n = len(x_vars)
    if len(reference_params) != len(u_vars):
        log.info("Error: len(reference_params) != len(u_vars)")
        log.info(1, "len(reference_params):", len(reference_params))
        log.info(1, "len(u_vars):", len(u_vars))
        exit()

    global_x_to_local_x = zeros(len(x_vars), len(x_names))
    for i in range(len(x_vars)):
        for j in range(len(x_names)):
            if x_vars[i] == x_names[j]:
                global_x_to_local_x[i,j] = 1
    
    global_x_r_to_local_x_r = zeros(len(reference_params), len(x_names))
    for i in range(len(reference_params)):
        for j in range(len(x_names)):
            if reference_params[i] == x_names[j]:
                global_x_r_to_local_x_r[i,j] = 1

    local_u_to_global_u = zeros(len(u_names), len(u_vars))
    for i in range(len(u_names)):
        for j in range(len(u_vars)):
            if u_names[i] == u_vars[j]:
                local_u_to_global_u[i,j] = 1

    C_r    =  zeros(len(reference_params), n)
    for i in range(len(x_vars)):
        if x_vars[i] in reference_params:
            C_r[reference_params.index(x_vars[i]),i] = 1
    
    y_re = zeros(len(reference_params), 1)
    for i in range(len(reference_params)):
        item = eval(reference_params[i])
        if item in x_e:
            y_re[i] = x_e[item]

    return global_x_to_local_x, global_x_r_to_local_x_r, local_u_to_global_u, C_r, y_re

def export(name, x_vars, u_vars, reference_params, x_e, x_dot_e = None, u_e = None):
    A, B, CC, x_e, u_e = computeStateSpace(name, x_vars, u_vars, x_e, x_dot_e, u_e)
    global_x_to_local_x, global_x_r_to_local_x_r, local_u_to_global_u, C_r, y_re = computeTransforms(x_vars, u_vars, reference_params, x_e)

    writeMathModule(name, A=A, B=B, C_r=C_r, CC=CC, x_e=x_e, u_e=u_e, y_re=y_re, global_x_to_local_x=global_x_to_local_x, global_x_r_to_local_x_r=global_x_r_to_local_x_r, local_u_to_global_u=local_u_to_global_u)


if compute_descent: descentCP()
if compute_flip:    flipCP()
if compute_landing: landing()