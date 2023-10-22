from tools.loadVTU import loadVTU

points, tris = loadVTU("meshes/Cylinder.vtu")
print(points)
print(tris)