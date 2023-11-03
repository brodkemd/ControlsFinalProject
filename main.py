from tools.rotations import Euler2Quaternion
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from mayavi import mlab

import parameters.simulation_parameters as SIM
from viewers.animation                  import Animation
from viewers.dataPlotter                import DataPlotter


make_output        = 1
show_figures       = 0
include_animation  = 0
include_plotter    = 0
write_data_to_file = 1


if include_animation:
    # read mesh and oriente along x-axis, that is what the scale matrix does
    animation = Animation("meshes/Body.vtu", scale=np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]))

if include_plotter or write_data_to_file:
    plotter = DataPlotter(width=SIM.fig_width, height=SIM.fig_height, plot=include_plotter)

# turn of screen rendering if show animation is false
if not show_figures:
    plt.close('all')
    mpl.use('Agg')

t = SIM.start_time
inputs = np.zeros((4, 1))
count = 0

print("\nProgress:")
while t < SIM.end_time:
    # setting next time to make a plot, exits inner loop at that time
    t_next_plot = t + SIM.ts_plot

    # loop that runs the dynamics
    while t < t_next_plot:
        y = [0, 0, 0, Euler2Quaternion(0, 0, np.deg2rad(t))]
        t += SIM.ts_simulation
        # Propagate dynamics

    # updating animation from result of dynamics
    if include_animation:
        animation.update(*y)

    # plotting the state variables and forces with respect to time
    if include_plotter or write_data_to_file:
        plotter.update(t, t, t**(1/2))

    # saving image of figure if on a valid timestep
    if include_plotter and make_output:
        plotter.savefig(SIM.output_plot_file_format.format(num=count))

    if show_figures:
        plt.pause(0.001) # allows animation to draw

    # -------increment time-------------
    count +=1
    
    # progresses bar with stats
    print(round(t, 6), "/", SIM.end_time, end="      \r")
print()

if write_data_to_file:
    plotter.saveData(SIM.data_file)