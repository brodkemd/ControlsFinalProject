import time
import numpy as np

# sample times, etc
ts_simulation = 0.01  # smallest time step for simulation
start_time    = 0.0  # start time for simulation
end_time      = 80.0  # end time for simulation
ts_plot       = 0.2  # refresh rate for plots
ts_control    = ts_simulation  # sample rate for the controller

# control figures
fig_width = 10
fig_height = 9

# for rendering purposes
prefix                   = str(round(time.time()*10))
output_plot_file_format  = f"{prefix}_"+"plot_{num}_output.jpg"
video_framerate          = 15
data_file                = f"{prefix}_"+"data.json"

# video_file_name = lambda : f"Va_{round(V_A)}_Y_{round(np.rad2deg(Y))}_R_{round(R) if not isinstance(R, str) else R}.mp4"
video_file_name = lambda : f"lon_mode_demo.mp4"