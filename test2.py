import vtk

file = "Cylinder.vtu"

# Create a 3D mesh (example: a cylinder)
cylinder = vtk.vtkCylinderSource()
cylinder.SetResolution(40)
cylinder.SetHeight(3.0)
cylinder.SetRadius(1.0)
cylinder.Update()

# Get points and cells from the mesh
mesh_polydata = cylinder.GetOutput()
points = mesh_polydata.GetPoints()
#cells = mesh_polydata.GetCells()

# Create a vtkUnstructuredGrid
unstructured_grid = vtk.vtkUnstructuredGrid()
unstructured_grid.SetPoints(points)

# Populate vtkUnstructuredGrid with cells
for i in range(mesh_polydata.GetNumberOfCells()):
    cell = mesh_polydata.GetCell(i)
    unstructured_grid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())

# Write vtkUnstructuredGrid to a VTU file
vtu_writer = vtk.vtkXMLUnstructuredGridWriter()
vtu_writer.SetFileName(file)
vtu_writer.SetInputData(unstructured_grid)
vtu_writer.Write()

print(f"Mesh exported to \"{file}\" successfully.")
