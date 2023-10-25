from viewers.animation import Animation
from tools.rotations import Euler2Quaternion
import numpy as np


animation = Animation("meshes/Body.vtu")
animation.update(0, 0, 0, Euler2Quaternion(0, 0, np.deg2rad(0)))
animation.update(1, 0, 0, Euler2Quaternion(0, 0, np.deg2rad(15)))
animation.update(2, 0, 0, Euler2Quaternion(0, 0, np.deg2rad(30)))
animation.update(3, 0, 0, Euler2Quaternion(0, 0, np.deg2rad(45)))