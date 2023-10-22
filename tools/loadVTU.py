import meshio
import numpy as np

def loadVTU(file:str):
    # Load VTU file
    mesh = meshio.read(file)

    # Access points
    points = mesh.points
    # for i, point in enumerate(points):
    #     print(f"Point {i}: {point}")

    # Access triangles (cell type 'triangle')
    triangles = mesh.cells_dict.get("triangle")

    # # Iterate through triangles
    # if triangles is not None:
    #     for i, triangle in enumerate(triangles):
    #         print(f"Triangle {i}: {triangle}")
    # else:
    #     print("No triangles found in the VTU file.")
    
    return np.array(points, dtype=float), np.array(triangles, dtype=int)