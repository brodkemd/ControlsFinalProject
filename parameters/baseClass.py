import parameters.body_parameters as BODY

class Base:
    def __init__(self) -> None:
        self.m     = BODY.mass
        self.g     = BODY.gravity
        self.Aref  = BODY.Aref
        self.J_xx  = BODY.J_xx
        self.J_yy  = BODY.J_yy
        self.J_zz  = BODY.J_zz
        self.J_xy  = BODY.J_xy
        self.J_xz  = BODY.J_xz
        self.J_yz  = BODY.J_yz
        self.J     = BODY.J
        self.r_E_x = BODY.r_E_x
        self.r_E_y = BODY.r_E_y
        self.r_E_z = BODY.r_E_z
        self.r_E   = BODY.r_E
        self.R     = BODY.R