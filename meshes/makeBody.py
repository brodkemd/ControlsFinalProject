from tools.makeSimpleMeshes import makeSimpleMeshes as Mesher
import os

cwd = os.path.dirname(__file__)
file = os.path.join(cwd, "Body.vtu")

mesher = Mesher()

radius      = 4
tip_height  = 6
body_height = 20
resolution  = 20

cone = mesher.cone(radius, tip_height, resolution)
cone = mesher.transform(cone, 0, body_height/2+tip_height/2, 0, 0, 0, 90)
cylinder = mesher.cylinder(radius, body_height, resolution)
body = mesher.appendObjects(cone, cylinder)
mesher.render(body)
mesher.exportToVtu(file, body)


# import vtk

# # Create the cylinder
# cylinder = vtk.vtkCylinderSource()
# cylinder.SetRadius(1.0)
# cylinder.SetHeight(2.0)
# cylinder.SetResolution(50)
# cylinder.Update()

# # Create the cone
# cone = vtk.vtkConeSource()
# cone.SetRadius(1.0)
# cone.SetHeight(2.0)
# cone.SetResolution(50)
# cone.Update()

# # Create a transform to position the cone on top of the cylinder
# transform = vtk.vtkTransform()
# transform.Translate(0, 0, 0)  # Translate the cone along the z-axis to place it on top of the cylinder
# transform.RotateZ(90)

# # Apply the transform to the cone
# transformFilter = vtk.vtkTransformPolyDataFilter()
# transformFilter.SetInputConnection(cone.GetOutputPort())
# transformFilter.SetTransform(transform)
# transformFilter.Update()

# # Append the cylinder and transformed cone
# appendFilter = vtk.vtkAppendPolyData()
# appendFilter.AddInputData(cylinder.GetOutput())
# appendFilter.AddInputData(transformFilter.GetOutput())
# appendFilter.Update()

# # Mapper
# mapper = vtk.vtkPolyDataMapper()
# mapper.SetInputData(appendFilter.GetOutput())

# # Actor
# actor = vtk.vtkActor()
# actor.SetMapper(mapper)

# # Renderer
# renderer = vtk.vtkRenderer()

# # Render Window
# renderWindow = vtk.vtkRenderWindow()
# renderWindow.SetWindowName("Cylinder with Cone on Top")
# renderWindow.AddRenderer(renderer)

# # Render Window Interactor
# renderWindowInteractor = vtk.vtkRenderWindowInteractor()
# renderWindowInteractor.SetRenderWindow(renderWindow)

# # Add actor to the renderer
# renderer.AddActor(actor)
# renderer.SetBackground(0.1, 0.1, 0.1)  # Set the background color

# # Set up the camera and render
# renderer.ResetCamera()
# renderer.GetActiveCamera().Azimuth(30)
# renderer.GetActiveCamera().Elevation(30)

# # Start the rendering loop
# renderWindow.Render()
# renderWindowInteractor.Start()




# # Create first object (for example, a cylinder)
# cylinder = vtk.vtkCylinderSource()
# cylinder.SetRadius(1.0)
# cylinder.SetHeight(2.0)
# cylinder.SetResolution(50)

# # Create second object (for example, a cone)
# cone = vtk.vtkConeSource()
# cone.SetRadius(1.0)
# cone.SetHeight(2.0)
# cone.SetResolution(50)

# # Append the two objects
# append_filter = vtk.vtkAppendPolyData()
# append_filter.AddInputData(cylinder.GetOutput())
# append_filter.AddInputData(cone.GetOutput())
# append_filter.Update()

# # Mapper
# mapper = vtk.vtkPolyDataMapper()
# mapper.SetInputConnection(append_filter.GetOutputPort())

# # Actor
# actor = vtk.vtkActor()
# actor.SetMapper(mapper)

# # Renderer
# renderer = vtk.vtkRenderer()

# # Render Window
# render_window = vtk.vtkRenderWindow()
# render_window.SetWindowName("Joined Objects Example")
# render_window.AddRenderer(renderer)

# # Render Window Interactor
# render_window_interactor = vtk.vtkRenderWindowInteractor()
# render_window_interactor.SetRenderWindow(render_window)

# # Add actor to the renderer
# renderer.AddActor(actor)
# renderer.SetBackground(0.1, 0.1, 0.1)  # Set the background color

# # Set up the camera and render
# renderer.ResetCamera()
# renderer.GetActiveCamera().Azimuth(30)
# renderer.GetActiveCamera().Elevation(30)

# # Start the rendering loop
# render_window.Render()
# render_window_interactor.Start()