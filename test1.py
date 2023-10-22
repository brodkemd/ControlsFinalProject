import vtk

# Create a cylinder
cylinder = vtk.vtkCylinderSource()
cylinder.SetResolution(20)  # Set the resolution of the cylinder
cylinder.SetHeight(3.0)      # Set the height of the cylinder
cylinder.SetRadius(1.0)      # Set the radius of the cylinder

# Create a mapper
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(cylinder.GetOutputPort())

# Create an actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)

# Create a renderer
renderer = vtk.vtkRenderer()

# Create a render window
render_window = vtk.vtkRenderWindow()
render_window.SetWindowName("Cylinder Mesh")
render_window.AddRenderer(renderer)

# Create a render window interactor
render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)

# Add the actor to the renderer
renderer.AddActor(actor)

# Set the background color
renderer.SetBackground(1.0, 1.0, 1.0)  # white

# Reset the camera
renderer.ResetCamera()

# Create a sample vtkPolyData (replace this with your own vtkPolyData)
input_polydata = mapper.GetInput()

# Extract surface geometry using vtkGeometryFilter
geometry_filter = vtk.vtkGeometryFilter()
geometry_filter.SetInputData(input_polydata)
geometry_filter.Update()

# Create vtkUnstructuredGrid from the extracted geometry
unstructured_grid = vtk.vtkUnstructuredGrid()
unstructured_grid.SetPoints(geometry_filter.GetOutput().GetPoints())
unstructured_grid.SetCells(geometry_filter.GetOutput().GetCellTypesArray(), geometry_filter.GetOutput().GetCellLocationsArray(), geometry_filter.GetOutput().GetPolys())

# Now you have unstructured_grid as vtkUnstructuredGrid based on input_polydata's geometry.

# Create a VTU writer
vtu_writer = vtk.vtkXMLUnstructuredGridWriter()
vtu_writer.SetFileName("cylinder.vtu")
print("here")
vtu_writer.SetInputData(unstructured_grid)  # Use the input of the mapper (cylinder)
print("here")
# Write the VTU file
vtu_writer.Write()
print("here")
# Render the scene
render_window.Render()

# Start the render loop
render_window_interactor.Start()

print("Cylinder VTU file written successfully.")