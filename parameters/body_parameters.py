from tools.rotations import Euler2Quaternion
import numpy as np

mass = 100.0
# gravity = 9.81
Aref = 1413.72 #Starship reference area, m^2
AsurfC = 15.81 #Canard ref area, m^2
AsurfF = 45.39 #Fin ref area, m^2
maxThrust = 9578.342 #N

# Atmospheric Parameters
rhoAvg = 0.016164897 #kg/m^3
gamma = 1.29
gravity = 3.73 #m/s^2
R = 192.0092 # Martian atmospheric gas constant, J/kgK
Tavg = 210 #Martian atmospheric temp, K

max_ground_incident_velocity = 10.0 # m/s

J_xx = 12814.379
J_yy = 940258.811
J_zz = 940258.811 
J_xy = 0.0
J_xz = 0.0
J_yz = 0.0

J = np.array([
    [ J_xx, -J_xy, -J_xz],
    [-J_xy,  J_yy, -J_yz],
    [-J_xz, -J_yz,  J_zz]
], dtype=float)

# use these euler angles to initial quaternion (easiest method)
phi   = 0
theta = np.deg2rad(0) # np.deg2rad(90)
psi   = 0

# defining the initial state
initial_state = np.array([
    -10, # p_n
    0, # p_e
   -11000, # p_d
    0, # u
    0, # v
    182.45, # w, M=0.8 for entry
    Euler2Quaternion(phi, theta, psi).item(0), # e_0
    Euler2Quaternion(phi, theta, psi).item(1), # e_1
    Euler2Quaternion(phi, theta, psi).item(2), # e_2
    Euler2Quaternion(phi, theta, psi).item(3), # e_3
    0, # p
    0, # q
    0 # r
], dtype=float)

r_E_x = -10.0
r_E_y = 0.0
r_E_z = 0.0

r_cp_port_canard_x      =  1
r_cp_port_canard_y      = -1
r_cp_port_canard_z      =  0
r_cp_starboard_canard_x =  1
r_cp_starboard_canard_y =  1
r_cp_starboard_canard_z =  0
r_cp_port_fin_x         = -1
r_cp_port_fin_y         = -1
r_cp_port_fin_z         =  0
r_cp_starboard_fin_x    = -1
r_cp_starboard_fin_y    =  1
r_cp_starboard_fin_z    =  0

r_cp_port_canard = np.array([r_cp_port_canard_x,r_cp_port_canard_y,r_cp_port_canard_z])
r_cp_starboard_canard = np.array([r_cp_starboard_canard_x,r_cp_starboard_canard_y,r_cp_starboard_canard_z])

r_cp_port_fin = np.array([r_cp_port_fin_x,r_cp_port_fin_y,r_cp_port_fin_z])
r_cp_starboard_fin = np.array([r_cp_starboard_fin_x,r_cp_starboard_fin_y,r_cp_starboard_fin_z])

r_E = np.array([r_E_x, r_E_y, r_E_z], dtype=float)

