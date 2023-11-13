import numpy as np

mass = 100
gravity = 9.81

J_xx = 10
J_yy = 10
J_zz = 10
J_xy = 10
J_xz = 10
J_yz = 10

J = np.array([
    [ J_xx, -J_xy, -J_xz],
    [-J_xy,  J_yy, -J_yz],
    [-J_xz, -J_yz,  J_zz]
], dtype=float)

r_E_x = 10
r_E_y = 0
r_E_z = 0

r_E = np.array([r_E_x, r_E_y, r_E_z], dtype=float)

