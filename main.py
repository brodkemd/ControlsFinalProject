from tools.rotations import Euler2Quaternion
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import platform, os

import parameters.simulation_parameters as SIM
from viewers.animation                  import Animation
from viewers.dataPlotter                import DataPlotter
from dynamics.dynamics                  import Dynamics
from dynamics.forcesMoments             import ForcesMomentsFromCP
from controller.fullStateFeedback       import FullStateFeedBack

# use these euler angles to initial quaternion (easiest method)
phi   = 0
theta = np.deg2rad(90)
psi   = 0

# defining the initial state
state = np.array([
    0, # p_n
    0, # p_e
   -1000, # p_d
    0, # u
    0, # v
    0, # w
    Euler2Quaternion(phi, theta, psi).item(0), # e_0
    Euler2Quaternion(phi, theta, psi).item(1), # e_1
    Euler2Quaternion(phi, theta, psi).item(2), # e_2
    Euler2Quaternion(phi, theta, psi).item(3), # e_3
    0, # p
    0, # q
    0 # r
], dtype=float)

state_reference    = state.copy()
state_reference[2] = 0

make_output        = 0
show_figures       = 0
include_animation  = 0
include_plotter    = 0
write_data_to_file = 1

if include_animation:
    # read mesh and oriente along x-axis, that is what the scale matrix does
    animation = Animation(
        file         = "meshes/Body.vtu",
        scale        = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]), 
        interactive  = show_figures,
        write_meshes = False,
        width        = SIM.fig_width,
        height       = SIM.fig_height
    )

if include_plotter or write_data_to_file:
    plotter = DataPlotter(
        width       = SIM.fig_width,
        height      = SIM.fig_height, 
        plot        = include_plotter,
        interactive = show_figures
    )

# turn of screen rendering if show animation is false
if not show_figures:
    plt.close('all')
    mpl.use('Agg')


controller = FullStateFeedBack(compute_gains=False)
dynamics   = Dynamics(state)
forces     = ForcesMomentsFromCP()


t = SIM.start_time
count = 0

print("\nProgress:")
while t < SIM.end_time:
    # setting next time to make a plot, exits inner loop at that time
    t_next_plot = t + SIM.ts_plot

    # loop that runs the dynamics
    while t < t_next_plot:
        F_E, F_cp, r_cp = controller.update(state, state_reference)
        # print(F_E)
        u        = forces.update(state, F_E, F_cp, r_cp)
        y, crash = dynamics.update(u) #[0, 0, 0, Euler2Quaternion(0, 0, np.deg2rad(t))]
        t += SIM.ts_simulation
        if crash: break

    # updating animation from result of dynamics
    if include_animation:
        animation.update(y)

    # plotting the state variables and forces with respect to time
    if include_plotter or write_data_to_file:
        response = [[y.item(i), state_reference.item(i)] for i in range(len(y))]
        plotter.update(t, response, F_E.tolist() + F_cp.tolist() + r_cp.tolist())

    # saving image of figure if on a valid timestep
    if include_plotter and make_output:
        plotter.savefig(SIM.output_plot_file_format.format(num=count), SIM.input_plot_file_format.format(num=count))
    if include_animation and make_output:
        animation.savefig(SIM.animation_file_format.format(num=count))
    if show_figures:
        plt.pause(0.001) # allows animation to draw

    # -------increment time-------------
    count +=1

    if crash:
        print("Detected crash, exiting simulation loop")
        break

    # progresses bar with stats
    print(round(t, 6), "/", SIM.end_time, end="      \r")
print()

if write_data_to_file: plotter.saveData(SIM.data_file)

if "linux" in platform.system().lower():
    if (include_animation or include_plotter) and make_output:
        # makes a video from the saved images
        #print((f"convert {animation_file_format} {plot_file_format} +append {combined_file_format}").format(num="*"))
        #exit()
        if include_plotter:
            if include_animation:
                for i in range(count):
                    # print(i)
                    # print((f"convert {animation_file_format} {plot_file_format} +append {combined_file_format}").format(num=i))
                    os.system((f"convert {SIM.animation_file_format} {SIM.input_plot_file_format} -append {SIM.intermediate_file_format}").format(num=i))
                    os.system(f"rm {SIM.animation_file_format} {SIM.input_plot_file_format}".format(num=i))
                
                for i in range(count):
                    # print(i)
                    # print((f"convert {animation_file_format} {plot_file_format} +append {combined_file_format}").format(num=i))
                    os.system((f"convert {SIM.intermediate_file_format} {SIM.output_plot_file_format} +append {SIM.combined_file_format}").format(num=i))
                    os.system(f"rm {SIM.intermediate_file_format} {SIM.output_plot_file_format}".format(num=i))
            else:
                for i in range(count):
                    os.system((f"convert {SIM.output_plot_file_format} {SIM.input_plot_file_format} -append {SIM.combined_file_format}").format(num=i))
                    os.system(f"rm {SIM.output_plot_file_format} {SIM.input_plot_file_format}".format(num=i))
                    

        else:
            if include_animation:
                for i in range(count):
                    os.system(f"mv {SIM.animation_file_format} {SIM.combined_file_format}".format(num=i))

        os.system((f"ffmpeg -framerate {SIM.video_framerate} -i {SIM.combined_file_format} -c:v libx264 -r 30 -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -y {SIM.video_file_name}").format(num="%d"))

        for i in range(count):
            os.system(f"rm {SIM.combined_file_format}".format(num=i))
        
        plt.close()
else:
    print("Can not render video on non-linux machines, use the standalone data plotter in tools")