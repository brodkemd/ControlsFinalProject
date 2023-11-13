import numpy as np

mass = 100

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