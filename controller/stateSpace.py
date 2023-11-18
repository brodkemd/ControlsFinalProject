import numpy as np
import os
from numpy import sqrt
from parameters.baseClass import Base
from tools.loadMathModule import LoadModule


class DescentStateSpaceCP(Base):
    def __init__(self) -> None:
        super().__init__()
        print(dir(self))
        cwd = os.path.dirname(__file__)
        LoadModule(os.path.join(cwd, "descent.h"), self)


class FlipStateSpaceCP(Base):
    def __init__(self) -> None:
        super().__init__()
        print(dir(self))
        cwd = os.path.dirname(__file__)
        LoadModule(os.path.join(cwd, "flip.h"), self)


class LandingStateSpace(Base):
    def __init__(self) -> None:
        super().__init__()
        print(dir(self))
        cwd = os.path.dirname(__file__)
        LoadModule(os.path.join(cwd, "landing.h"), self)

        # self.A = np.array([
        #     [0, 0, 0, 0, 0, 1, 0, 0],
        #     [0, 0, 0, 0, 1, 0, 0, 0],
        #     [0, 0, 0, -1, 0, 0, 0, 0],
        #     [0, 0, 0, 0, 0, 0, -sqrt(2)*self.g, 0],
        #     [0, 0, 0, 0, 0, 0, 0, 0],
        #     [0, 0, 0, 0, 0, 0, sqrt(2)*self.g, 0],
        #     [0, 0, 0, 0, 0, 0, 0, -sqrt(2)/4],
        #     [0, 0, 0, 0, 0, 0, 0, 0]
        # ], dtype=float)
        
        # self.B = np.array([
        #     [0, 0, 0],
        #     [0, 0, 0],
        #     [0, 0, 0],
        #     [1/self.m, 0, 0],
        #     [0, 1/self.m, 0],
        #     [0, 0, 1/self.m],
        #     [0, 0, 0],
        #     [0, 0, -self.r_E_x/self.J_yy]
        # ], dtype=float)

        # self.C = np.eye(len(self.A), dtype=float)

        # self.x_e = np.array([0, 0, 0, 0, 0, 0, sqrt(2)/2, 0, sqrt(2)/2, 0, 0, 0, 0], dtype=float)
        # self.u_e = np.array([self.g*self.m, 0, 0, 0, 0, 0, 0, 0, self.R],            dtype=float)
