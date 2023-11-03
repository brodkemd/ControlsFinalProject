import sys
sys.path.append('.')# one directory up
import numpy as np

from tools.rotations import Euler2Rotation, Quaternion2Rotation
from tools.loadVTU import loadVTU
import meshio
from mayavi import mlab


class Animation():

    centroid = []
    trail_points = []

    points    = []
    triangles = []

    show_trail = False
    interactive = False
    write_meshes = False
    fig = None

    def __init__(self, file:str, scale:np.ndarray=np.eye(3), shift:np.ndarray=np.zeros(3), name="animation", interactive=False, write_meshes=True):
        print("\nAnimation:")
        self.count = 0
        self.name = name
        self.createObject(file, scale, shift)


        # rotate for plotting north=x east=y h=-z
        self.R_plot=np.array([
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1]
        ]) # for plotting's sake

        self.interactive = False# interactive
        self.write_meshes = write_meshes
        self.init = False

        print("  Rendering Setup:")
        print("    interactive:", self.interactive)
        print("    writing meshes:", self.write_meshes)



    def vertices(self, 
            position_north, # position north
            position_east, # position east
            position_down, # position down
            quaternion
        ):
        plot_points = []
        
        shift = np.array([position_north, position_east, position_down])
        
        R = Quaternion2Rotation(quaternion)
        
        shift = self.R_plot @ shift
        self.centroid = shift
        
        T = self.R_plot @ R # making the total transformation matrix

        for i in range(len(self.points)):
            plot_points.append((np.matmul(T, self.points[i]) + shift).tolist())
        
        return np.array(plot_points)


    def update(self, 
            position_north, # position north
            position_east, # position east
            position_down, # position down
            quaternion
        ):
        # draw plot elements: cart, bob, rod
        self.drawObject(position_north, position_east, position_down, quaternion)



    def drawObject(self, 
            position_north, # position north
            position_east, # position east
            position_down, # position down
            quaternion
        ):
        
        verts=self.vertices(position_north, position_east, position_down, quaternion)
        
        if self.write_meshes:
            meshio.write_points_cells(
                f"{self.name}_t{self.count}.vtu",
                verts,
                {"triangle" : self.triangles}
            )
        # if self.interactive:
        #     self.fig.mlab_source.set(x=self.points[:,0])
        #     self.fig.mlab_source.set(y=self.points[:,1])
        #     self.fig.mlab_source.set(z=self.points[:,2])
        #     self.fig.draw()

        self.count += 1


    def createObject(self, file:str, scale, shift):
        self.points, self.triangles = loadVTU(file)

        for i in range(len(self.points)):
            self.points[i] = np.matmul(scale, self.points[i]) + shift

        self.analyze()
        # if self.interactive:
        #     print("Making Visual Animation")
        #     self.fig = mlab.triangular_mesh(self.points[:,0], self.points[:,1], self.points[:,2], self.triangles)


    def analyze(self):
        round_to = 4
        x_diff = max(self.points[:,0]) - min(self.points[:,0])
        y_diff = max(self.points[:,1]) - min(self.points[:,1])
        z_diff = max(self.points[:,2]) - min(self.points[:,2])

        space = "  "
        print(space+"Mesh Details:", end="")
        greatest_side = max([x_diff, y_diff, z_diff])
        if greatest_side == x_diff:
            print(" length oriented along x")
        elif greatest_side == y_diff:
            print(" length oriented along y")
        else:
            print(" length oriented along z")
        print(space+"  x details:")
        print(space+"    - max x:", round(max(self.points[:,0]), round_to))
        print(space+"    - min x:", round(min(self.points[:,0]), round_to))
        print(space+"    - x spread:", round(x_diff, round_to))
        print(space+"  y details:")
        print(space+"    - max y:", round(max(self.points[:,1]), round_to))
        print(space+"    - min y:", round(min(self.points[:,1]), round_to))
        print(space+"    - y spread:", round(y_diff, round_to))
        print(space+"  z details:")
        print(space+"    - max z:", round(max(self.points[:,2]), round_to))
        print(space+"    - min z:", round(min(self.points[:,2]), round_to))
        print(space+"    - z spread:", round(z_diff, round_to))
        



    def savefig(self, file, *args, **kwargs):
        self.fig.savefig(file, *args, **kwargs)