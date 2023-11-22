from sympy import symbols, Matrix

f_body_x,f_body_y,f_body_z = symbols("f_body_x,f_body_y,f_body_z")

f_port_canard_x,f_port_canard_y,f_port_canard_z = symbols("f_port_canard_x,f_port_canard_y,f_port_canard_z")
f_starboard_canard_x,f_starboard_canard_y,f_starboard_canard_z = symbols("f_starboard_canard_x,f_starboard_canard_y,f_starboard_canard_z")

f_port_fin_x,f_port_fin_y,f_port_fin_z = symbols("f_port_fin_x,f_port_fin_y,f_port_fin_z")
f_starboard_fin_x,f_starboard_fin_y,f_starboard_fin_z = symbols("f_starboard_fin_x,f_starboard_fin_y,f_starboard_fin_z")

r_body_x,r_body_y,r_body_z = symbols("r_body_x,r_body_y,r_body_z")

r_port_canard_x,r_port_canard_y,r_port_canard_z = symbols("r_port_canard_x,r_port_canard_y,r_port_canard_z")
r_starboard_canard_x,r_starboard_canard_y,r_starboard_canard_z = symbols("r_starboard_canard_x,r_starboard_canard_y,r_starboard_canard_z")

r_port_fin_x,r_port_fin_y,r_port_fin_z = symbols("r_port_fin_x,r_port_fin_y,r_port_fin_z")
r_starboard_fin_x,r_starboard_fin_y,r_starboard_fin_z = symbols("r_starboard_fin_x,r_starboard_fin_y,r_starboard_fin_z")

F_body = Matrix([f_body_x,f_body_y,f_body_z])

F_port_canard = Matrix([f_port_canard_x,f_port_canard_y,f_port_canard_z])
F_starboard_canard = Matrix([f_starboard_canard_x,f_starboard_canard_y,f_starboard_canard_z])

F_port_fin = Matrix([f_port_fin_x,f_port_fin_y,f_port_fin_z])
F_starboard_fin = Matrix([f_starboard_fin_x,f_starboard_fin_y,f_starboard_fin_z])

F_canard = F_port_canard + F_starboard_canard
F_fin = F_port_fin + F_starboard_fin

F = F_body + F_fin + F_canard

r_body = Matrix([r_body_x,r_body_y,r_body_z])

r_port_canard = Matrix([r_port_canard_x,r_port_canard_y,r_port_canard_z])
r_starboard_canard = Matrix([r_starboard_canard_x,r_starboard_canard_y,r_starboard_canard_z])

r_port_fin = Matrix([r_port_fin_x,r_port_fin_y,r_port_fin_z])
r_starboard_fin = Matrix([r_starboard_fin_x,r_starboard_fin_y,r_starboard_fin_z])


tau_port_canard      = r_port_canard.cross(F_port_canard)
tau_starboard_canard = r_starboard_canard.cross(F_starboard_canard)
tau_port_fin         = r_port_fin.cross(F_port_fin)
tau_starboard_fin    = r_starboard_fin.cross(F_starboard_fin)
tau_body             = r_body.cross(F_body)

tau = tau_port_canard + tau_port_fin + tau_starboard_canard + tau_starboard_fin + tau_body

r_cp = (F.cross(tau))/(F.dot(F))
print(r_cp)


f_body_x_val = 0
f_body_y_val = 0
f_body_z_val = 0

f_port_canard_x_val = 0
f_port_canard_y_val = 0
f_port_canard_z_val = 0

f_starboard_canard_x_val = 0
f_starboard_canard_y_val = 0
f_starboard_canard_z_val = 0

f_port_fin_x_val = 0
f_port_fin_y_val = 0
f_port_fin_z_val = 1

f_starboard_fin_x_val = 0
f_starboard_fin_y_val = 0
f_starboard_fin_z_val = 1

r_body_x_val = 0
r_body_y_val = 0
r_body_z_val = 0

r_port_canard_x_val =  1
r_port_canard_y_val = -1
r_port_canard_z_val =  0

r_starboard_canard_x_val = 1
r_starboard_canard_y_val = 1
r_starboard_canard_z_val = 0

r_port_fin_x_val = -1
r_port_fin_y_val = -1
r_port_fin_z_val =  0

r_starboard_fin_x_val = -1
r_starboard_fin_y_val =  1
r_starboard_fin_z_val =  0


names = "f_body_x,f_body_y,f_body_z,f_port_canard_x,f_port_canard_y,f_port_canard_z,f_starboard_canard_x,f_starboard_canard_y,f_starboard_canard_z,f_port_fin_x,f_port_fin_y,f_port_fin_z,f_starboard_fin_x,f_starboard_fin_y,f_starboard_fin_z,r_body_x,r_body_y,r_body_z,r_port_canard_x,r_port_canard_y,r_port_canard_z,r_starboard_canard_x,r_starboard_canard_y,r_starboard_canard_z,r_port_fin_x,r_port_fin_y,r_port_fin_z,r_starboard_fin_x,r_starboard_fin_y,r_starboard_fin_z".split(",")

vals = []
for name in names:
    vals.append((eval(name), eval(f"{name}_val")))

r_cp_val = r_cp.subs(vals)
print(r_cp_val)
