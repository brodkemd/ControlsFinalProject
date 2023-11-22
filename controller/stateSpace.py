import numpy as np
import os
from numpy import sqrt
from parameters.baseClass import Base
from tools.loadMathModule import LoadModule


class BaseStateSpace(Base):
    # where the file is located
    cwd = os.path.dirname(__file__)

    # initializing for convenience
    A, B, C_r, x_e, u_e, y_re = None, None, None, None, None, None

    def __init__(self, module_file:str) -> None:
        super().__init__()
        print(f"    Loading module: '{module_file}'")
        LoadModule(os.path.join(self.cwd, module_file), self)


class DescentStateSpaceCP(BaseStateSpace):
    def __init__(self) -> None:
        super().__init__("descentCP.h")


class FlipStateSpaceCP(BaseStateSpace):
    def __init__(self) -> None:
        super().__init__("flipCP.h")


class LandingStateSpace(BaseStateSpace):
    def __init__(self) -> None:
        super().__init__("landing.h")
