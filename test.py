from scipy import integrate

# Define the density function for the shape (replace this with the actual density function)
def density_function(x, y, z):
    # Example density function (replace this with the actual function)
    return x**2 + y**2 + z**2

# Define the integration limits based on the shape's geometry
x_min, x_max = 0, 1
y_min, y_max = 0, 1
z_min, z_max = 0, 1

# Perform triple integration to calculate moments of inertia
Ix, err_x = integrate.tplquad(lambda x, y, z: density_function(x, y, z) * (y**2 + z**2),
                               x_min, x_max,
                               lambda x: y_min, lambda x: y_max,
                               lambda x, y: z_min, lambda x, y: z_max)

Iy, err_y = integrate.tplquad(lambda x, y, z: density_function(x, y, z) * (x**2 + z**2),
                               y_min, y_max,
                               lambda y: x_min, lambda y: x_max,
                               lambda x, y: z_min, lambda x, y: z_max)

Iz, err_z = integrate.tplquad(lambda x, y, z: density_function(x, y, z) * (x**2 + y**2),
                               z_min, z_max,
                               lambda z: x_min, lambda z: x_max,
                               lambda x, y: y_min, lambda x, y: y_max)

print("Moment of inertia about the x-axis:", Ix)
print("Moment of inertia about the y-axis:", Iy)
print("Moment of inertia about the z-axis:", Iz)