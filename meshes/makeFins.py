from tools.makeSimpleMeshes import makeSimpleMeshes as Mesher
import os
import meshes.makeBody as P
from math import acos, atan2, degrees, pi, cos, sin, sqrt

# size of rectangle
w = 1 
l = 2
t = 0.1

delta = 0.5*sqrt(P.tip_height**2 + (P.radius/2)**2) # position along nose starting from bottom of nose

cwd       = os.path.dirname(__file__)
body_file = os.path.join(cwd, "Body.vtu")
file      = os.path.join(cwd, "canard.vtu")

theta = atan2(P.tip_height, P.radius/2)
print(degrees(theta))
phi = pi/2 - theta

x = P.radius - delta*cos(theta) + w/2*cos(phi)
y = P.body_height/2 + delta*sin(theta) + w/2*sin(phi)
print(x, y)
mesher = Mesher()

body = mesher.loadVTU(body_file)
canard = mesher.rectangularPrism(w, l, t, 10)
canard = mesher.transform(canard, -x, y, 0, 0, 0, -30)
whole = mesher.appendObjects(body, canard)
#mesher.render(whole)
mesher.exportToVtu(file, whole)