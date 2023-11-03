import vtk

class makeSimpleMeshes:
    def __init__(self) -> None:
        pass


    def appendObjects(self, obj_1, obj_2):
        append_filter = vtk.vtkAppendPolyData()
        append_filter.AddInputData(obj_1)
        append_filter.AddInputData(obj_2)
        append_filter.Update()
        return append_filter.GetOutput()


    def transform(self, obj, dx:float, dy:float, dz:float, rotate_X:float, rotate_Y:float, rotate_Z:float):
        # Create a transform to position the cone on top of the cylinder
        transform = vtk.vtkTransform()
        transform.Translate(dx, dy, dz)  # Translate the cone along the z-axis to place it on top of the cylinder
        transform.RotateX(rotate_X)
        transform.RotateY(rotate_Y)
        transform.RotateZ(rotate_Z)

        # Apply the transform to the cone
        transformFilter = vtk.vtkTransformPolyDataFilter()
        transformFilter.SetInputData(obj)
        transformFilter.SetTransform(transform)
        transformFilter.Update()
        return transformFilter.GetOutput()


    def exportToVtu(self, file:str, obj):
        triangle_filter = vtk.vtkTriangleFilter()
        triangle_filter.SetInputData(obj)
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


    def cone(self, radius:float, height:float, resolution:int, center=(0, 0, 0)):
        # Create a cone source
        cone = vtk.vtkConeSource()
        cone.SetRadius(radius)  # Set the radius of the base of the cone
        cone.SetHeight(height)  # Set the height of the cone
        cone.SetResolution(resolution)  # Set the resolution of the cone (number of sides)
        cone.SetCenter(*center)
        cone.Update()
        return cone.GetOutput()


    def cylinder(self, radius:float, height:float, resolution:int, center=(0, 0, 0)):
        cylinder = vtk.vtkCylinderSource()
        cylinder.SetHeight(height)
        cylinder.SetRadius(radius)
        cylinder.SetResolution(resolution)  # Set the resolution of the cone (number of sides)
        cylinder.SetCenter(*center)
        cylinder.Update()
        return cylinder.GetOutput()


    def rectangularPrism(self, x_width:float, y_width:float, z_width:float, resolution:int, center=(0, 0, 0)):
        # Create a cube source (rectangular prism)
        cube = vtk.vtkCubeSource()
        cube.SetXLength(x_width)  # Set the length of the rectangular prism along the x-axis
        cube.SetYLength(y_width)  # Set the length of the rectangular prism along the y-axis
        cube.SetZLength(z_width)  # Set the length of the rectangular prism along the z-axis
        #cube.SetResolution(resolution)  # Set the resolution of the cone (number of sides)
        cube.SetCenter(*center)
        cube.Update()
        return cube.GetOutput()


    def loadVTU(self, file:str):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(file)
        reader.Update()
    
        geometry_filter = vtk.vtkGeometryFilter()
        geometry_filter.SetInputData(reader.GetOutput())
        geometry_filter.Update()
        return geometry_filter.GetOutput()


    def render(self, obj):
        # # Mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(obj)

        # Actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # Renderer
        renderer = vtk.vtkRenderer()

        # Render Window
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.SetWindowName("Cylinder with Cone on Top")
        renderWindow.AddRenderer(renderer)

        # Render Window Interactor
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderWindowInteractor.SetRenderWindow(renderWindow)

        # Add actor to the renderer
        renderer.AddActor(actor)
        renderer.SetBackground(0.1, 0.1, 0.1)  # Set the background color

        # Set up the camera and render
        renderer.ResetCamera()
        renderer.GetActiveCamera().Azimuth(30)
        renderer.GetActiveCamera().Elevation(30)

        # Start the rendering loop
        renderWindow.Render()
        renderWindowInteractor.Start()
        print("here")
