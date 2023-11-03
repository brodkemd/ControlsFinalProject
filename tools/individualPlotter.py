import sys
sys.path.append(".")
import parameters.simulation_parameters as SIM
import matplotlib.pyplot as plt
import numpy as np
import json, sys

from viewers.dataPlotter import DataPlotter

plotter = DataPlotter(width=SIM.fig_width, height=SIM.fig_height, plot=True)

file = sys.argv[1]

with open(file, 'r') as f:
    data = json.loads(f.read())

plotter.t = data["t"]
del data["t"]
for name in data:
    for item in data[name]:
        if isinstance(data[name][item], list):
            exec(f"plotter.{item} = {data[name][item]}")
        else:
            for dim in data[name][item]:
                exec(f"plotter.{item}_{dim} = {data[name][item][dim]}")
plotter.createPlot()
plotter.orientePlot()
plt.show(block=True)