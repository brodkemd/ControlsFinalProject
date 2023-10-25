"""
Class for plotting a little plane

Author: Marek Brodke # 
"""
import sys
sys.path.append('.')# one directory up
import numpy as np

from tools.rotations import Euler2Rotation, Quaternion2Rotation
from tools.loadVTU import loadVTU
import meshio

class Animation():

    centroid = []
    trail_points = []

    points    = []
    triangles = []

    show_trail = False

    def __init__(self, file:str, scale:np.ndarray=np.eye(3), shift:np.ndarray=np.zeros(3), name="animation"):
        self.count = 0
        self.name = name
        self.createObject(file, scale, shift)


        # rotate for plotting north=x east=y h=-z
        self.R_plot=np.array([
                [0, 1, 0],
                [1, 0, 0],
                [0, 0, -1]
        ]) # for plotting's sake

        
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
        
        return plot_points


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
        meshio.write_points_cells(
            f"{self.name}_t{self.count}.vtu",
            verts,
            {"triangle" : self.triangles}
        )
        self.count += 1


    def createObject(self, file:str, scale, shift):
        self.points, self.triangles = loadVTU(file)

        for i in range(len(self.points)):
            self.points[i] = np.matmul(scale, self.points[i]) + shift

    
    def savefig(self, file, *args, **kwargs):
        self.fig.savefig(file, *args, **kwargs)