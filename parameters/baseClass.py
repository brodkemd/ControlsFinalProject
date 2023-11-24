import parameters.body_parameters as BODY
import numpy as np

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

        self.r_cp_port_canard_x      = BODY.r_cp_port_canard_x
        self.r_cp_port_canard_y      = BODY.r_cp_port_canard_y
        self.r_cp_port_canard_z      = BODY.r_cp_port_canard_z
        self.r_cp_starboard_canard_x = BODY.r_cp_starboard_canard_x
        self.r_cp_starboard_canard_y = BODY.r_cp_starboard_canard_y
        self.r_cp_starboard_canard_z = BODY.r_cp_starboard_canard_z
        self.r_cp_port_fin_x         = BODY.r_cp_port_fin_x
        self.r_cp_port_fin_y         = BODY.r_cp_port_fin_y
        self.r_cp_port_fin_z         = BODY.r_cp_port_fin_z
        self.r_cp_starboard_fin_x    = BODY.r_cp_starboard_fin_x
        self.r_cp_starboard_fin_y    = BODY.r_cp_starboard_fin_y
        self.r_cp_starboard_fin_z    = BODY.r_cp_starboard_fin_z
        self.r_cp_port_canard        = BODY.r_cp_port_canard
        self.r_cp_starboard_canard   = BODY.r_cp_starboard_canard
        self.r_cp_port_fin           = BODY.r_cp_port_fin
        self.r_cp_starboard_fin      = BODY.r_cp_starboard_fin