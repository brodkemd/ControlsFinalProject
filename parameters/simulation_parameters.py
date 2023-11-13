import time, os
import numpy as np

# sample times, etc
ts_simulation = 0.01  # smallest time step for simulation
start_time    = 0.0  # start time for simulation
end_time      = 2.0  # end time for simulation
ts_plot       = 0.2  # refresh rate for plots
ts_control    = ts_simulation  # sample rate for the controller

# control figures
fig_width = 10
fig_height = 9

data_dir                 = os.path.join(os.getcwd(), "data")

# for rendering purposes
prefix                   = str(round(time.time()*10))
output_plot_file_format  = os.path.join(data_dir, f"{prefix}_"+"plot_{num}_output.jpg")
input_plot_file_format   = os.path.join(data_dir, f"{prefix}_"+"plot_{num}_input.jpg")
animation_file_format    = os.path.join(data_dir, f"{prefix}_"+"animation_{num}.jpg")
intermediate_file_format = os.path.join(data_dir, f"{prefix}_"+"temp_{num}.jpg")
combined_file_format     = os.path.join(data_dir, f"{prefix}_"+"{num}.jpg")
video_framerate          = 15
data_file                = os.path.join(data_dir, f"{prefix}_"+"data.json")

video_file_name          = os.path.join(data_dir, "output.mp4")