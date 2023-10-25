import vtk

file = "Cone.vtu"

# Create a cone source
cone = vtk.vtkConeSource()
cone.SetRadius(1.0)  # Set the radius of the base of the cone
cone.SetHeight(2.0)  # Set the height of the cone
cone.SetResolution(50)  # Set the resolution of the cone (number of sides)

triangle_filter = vtk.vtkTriangleFilter()
triangle_filter.SetInputConnection(cone.GetOutputPort())
triangle_filter.Update()

# Get points and cells from the mesh
mesh_polydata = triangle_filter.GetOutput()
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
