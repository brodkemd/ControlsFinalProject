import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import platform, os, traceback, warnings
warnings.filterwarnings("error")

from tools.rotations import Euler2Rotation
import parameters.simulation_parameters as SIM
import parameters.body_parameters       as P
from viewers.animation                  import Animation
from viewers.dataPlotter                import DataPlotter
from dynamics.dynamics                  import Dynamics
from dynamics.forcesMoments             import ForcesMoments
from controller.LQR                     import LQR

make_output        = 0
show_figures       = 1
include_animation  = 1
include_plotter    = 0
write_meshes       = 0
write_data_to_file = 1

if include_animation:
    # read mesh and oriente along x-axis, that is what the scale matrix does
    animation = Animation(
        file         = "meshes/Body.vtu",
        scale        = 3*np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]), 
        interactive  = show_figures,
        write_meshes = write_meshes,
        width        = SIM.fig_width,
        height       = SIM.fig_height,
    )
    
    # animation = Animation(
    #     file         = "meshes/starship.vtu",
    #     scale        = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]), 
    #     interactive  = show_figures,
    #     write_meshes = write_meshes,
    #     width        = SIM.fig_width,
    #     height       = SIM.fig_height,
    # )

# animation.update(P.initial_state)
# plt.show(block=True)
# exit()

if include_plotter or write_data_to_file:
    plotter = DataPlotter(width = SIM.fig_width, height = SIM.fig_height, plot = include_plotter, interactive = show_figures)

# turn of screen rendering if show animation is false
if not show_figures:
    plt.close('all')
    mpl.use('Agg')

state      = P.initial_state.copy()
#controller = FullStateFeedBack(compute_gains=True)
controller = LQR(compute_gains=False)
dynamics   = Dynamics(state)
forces     = ForcesMoments()

try:
    print("\nProgress:")
    t = SIM.start_time
    count = 0
    while t < SIM.end_time:
        # setting next time to make a plot, exits inner loop at that time
        t_next_plot = t + SIM.ts_plot

        # loop that runs the dynamics
        while t < t_next_plot:
            F_E, angles, x_r    = controller.update(state)
            u                   = forces.update(state, F_E, angles)
            y, crash, landed    = dynamics.update(u)
            t += SIM.ts_simulation
            if crash or landed: break

        # updating animation from result of dynamics
        if include_animation: animation.update(y)

        # plotting the state variables and forces with respect to time
        if include_plotter or write_data_to_file:
            response = [[y.item(i), x_r.item(i)] for i in range(len(y))]
            plotter.update(t, response, F_E.tolist() + angles.tolist())

        # saving image of figure if on a valid timestep
        if include_plotter and make_output:
            plotter.savefig(SIM.output_plot_file_format.format(num=count), SIM.input_plot_file_format.format(num=count))
        if include_animation and make_output:
            animation.savefig(SIM.animation_file_format.format(num=count))
        if show_figures: plt.pause(0.001) # allows animation to draw

        # -------increment time-------------
        count+=1

        if crash:
            print("Detected crash :(\nexiting simulation loop")
            break
        if landed:
            print("Detected landing, yay!!!!!\nexiting simulation loop")
            break

        # progresses bar with stats
        #print(round(t, 6), "/", SIM.end_time, round(np.rad2deg(Quaternion2Euler(y[6:10])[1]),2), end="      \r")
        print(round(t, 6), "/", SIM.end_time, end="      \r")
    print()

except KeyboardInterrupt:
    # catches user interrupts, stops wordy error messages that occur due a usual user interruption
    print("\nCaught keyboard interrupt, Exiting ...")
    exit()

except Exception:
    # saves data if there is an error
    traceback.print_exc()
    print("\nDumping collected data")
    if write_data_to_file: plotter.saveData(SIM.data_file)
    exit()






################
#
# The rest is just post processing
#
################



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
    if make_output:
        print("Can not render video on non-linux machines, use the standalone data plotter in tools")