from mayavi import mlab
import numpy as np
import time
# Create some data (example data)
x, y, z = np.mgrid[-5:5:50j, -5:5:50j, -5:5:50j]
scalar_field = x**2 + y**2 + z**2

# Create a Mayavi figure
fig = mlab.figure(size=(800, 600))

# Create an initial visualization
contour = mlab.contour3d(x, y, z, scalar_field, contours=10, colormap='viridis')

# Show the initial visualization
mlab.draw()
mlab.process_ui_events()
time.sleep(5)

# Update the visualization (for example, change contour levels) without closing the window
contour.contour.number_of_contours = 5
mlab.draw()

# Continue the script
input("Press Enter to exit.")