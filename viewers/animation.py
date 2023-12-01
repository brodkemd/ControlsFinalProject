import sys, os
sys.path.append('.')# one directory up
import numpy as np

from tools.rotations import Euler2Rotation, Quaternion2Rotation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from tools.loadVTU import loadVTU
import meshio
import matplotlib.pyplot as plt


class Animation():
    centroid = []
    trail_points = []

    points    = []
    triangles = []

    interactive = False
    write_meshes = False

    def __init__(self, file:str, scale:np.ndarray=np.eye(3), shift:np.ndarray=np.zeros(3), name="animation", interactive=False, write_meshes=True, width=10, height=9):
        print("\nAnimation:")
        self.count = 0
        self.name = name

        self.fig = plt.figure(figsize=(width, height))
        self.ax:plt.Axes = self.fig.add_subplot(1, 1, 1, projection='3d')

        self.ax.set_ylabel('North(m)')
        self.ax.set_xlabel('East(m)')
        self.ax.set_zlabel('Height(m)')
        self.ax.set_ylim([-5000,250])
        self.ax.set_xlim([-5250/2,5250/2])
        self.ax.set_zlim([-50,12000])

        self.createObject(file, scale, shift)


        self.fig.tight_layout()

        # rotate for plotting north=x east=y h=z
        self.R_plot=np.array([
                [0, 1, 0],
                [1, 0, 0],
                [0, 0, -1]
        ]) # for plotting's sake

        self.interactive = interactive
        self.write_meshes = write_meshes
        self.init = True

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
        # self.centroid[2] += 78
        
        T = self.R_plot @ R # making the total transformation matrix

        for i in range(len(self.points)):
            plot_points.append((np.matmul(T, self.points[i]) + shift).tolist())
        
        return np.array(plot_points)


    def update(self, 
            _state
        ):
        pn = _state.item(0) # position north
        pe = _state.item(1) # position east
        pd = _state.item(2) # position down
        e  = _state[6:10]

        

        # draw plot elements: cart, bob, rod
        self.drawObject(pn, pe, pd, e)

        self.ax.set_ylim([-1500,10])
        self.ax.set_xlim([-5250/10,5250/10])        
        if -pd-500 > 0:
            self.ax.set_zlim([-pd-500,-pd+500])
        else:
            self.ax.set_zlim([0,1000])

        # Set initialization flag to False after first call
        if self.init:
            self.init = False


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
            self.count += 1

        buffer = 50
        self.ax.set_xlim([-buffer + self.centroid[0], buffer + self.centroid[0]])
        self.ax.set_ylim([-buffer + self.centroid[1], buffer + self.centroid[1]])
        self.ax.set_zlim([-buffer + self.centroid[2], buffer + self.centroid[2]])
        
        plot_points = []
        for i in range(len(self.triangles)):
            plot_points.append([])
            for j in range(3):
                plot_points[-1].append((verts[self.triangles[i][j]]).tolist())

        if self.init:
            poly = Poly3DCollection(plot_points, alpha=1.0)
            self.cube =self.ax.add_collection3d(poly)#
        else:
            self.cube.set_verts(plot_points)

    def createObject(self, file:str, scale, shift):
        self.points, self.triangles = loadVTU(file)

        for i in range(len(self.points)):
            self.points[i] = np.matmul(scale, self.points[i]) + shift

        self.analyze()
        # if self.interactive:
        #     print("Making Visual Animation")
        #     self.fig = mlab.triangular_mesh(self.points[:,0], self.points[:,1], self.points[:,2], self.triangles)


    def analyze(self):
        round_to = 8
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

        x_c, y_c, z_c = self.computeCentroid()

        print(space+"Centroid:")
        print(space+"  - x:", round(x_c, round_to))
        print(space+"  - y:", round(y_c, round_to))
        print(space+"  - z:", round(z_c, round_to))
        print(space+"Normalizing about centroid:", end=" ")
        self.normalizeAboutCentriod()
        print("done")
        x_c, y_c, z_c = self.computeCentroid()
        print(space+"New Centroid:")
        print(space+"  - x:", round(x_c, round_to))
        print(space+"  - y:", round(y_c, round_to))
        print(space+"  - z:", round(z_c, round_to))


    def computeCentroid(self):
        # assumes constant density, so just average the points
        x_c = np.average(self.points[:,0])
        y_c = np.average(self.points[:,1])
        z_c = np.average(self.points[:,2])
        self.centroid = [x_c, y_c, z_c]
        return x_c, y_c, z_c

    # def computeMoments(self):
    #     theta = np.arange(0, np.pi/2, 0.1)
    #     phi = np.arange(0, 2*np.pi, 0.1)

    #     for i in range(len(theta)):
    #         for j in range(len(phi)):
    #             r = 1
    #             vec = r*np.array([np.sin(theta[i])*np.cos(phi[j]), np.sin(theta[i])*np.sin(phi[i]), np.cos(theta[i])])





    def normalizeAboutCentriod(self):
        for i in range(3):
            self.points[:,i] = self.points[:,i] - self.centroid[i]


    def savefig(self, file, *args, **kwargs):
        self.fig.savefig(file, *args, **kwargs)