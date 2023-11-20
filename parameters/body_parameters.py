import numpy as np

mass = 100.0
gravity = 9.81
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

r_E_x = 10.0
r_E_y = 0.0
r_E_z = 0.0

r_E = np.array([r_E_x, r_E_y, r_E_z], dtype=float)

