import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import numpy as np
import json

plt.ion()  # enable interactive drawing

# def wrapString(s):
#    return f"\'{s}\'"

class DataPlotter:
    special_words       = ["phi", "theta", "psi", "alpha", "beta", "delta"]

    vars                = ["x", "y"]
    vars_lengths        = np.zeros(len(vars), dtype=int).tolist()
    layout              = (2, 1)

    plot = True
    dimension_names = ["North", "East", "Down"]
    colors          = ["tab:blue", "tab:orange", "tab:green"]

    def __init__(self, width=10, height=8, plot=True) -> None:
        print("\nData Plotter:", end="")
        self.plot = plot

        if self.plot:
            print()
            print("  Creating plot grid:", f"{self.layout[0]}x{self.layout[1]}")
            self.fig, axes_temp           = plt.subplots(*self.layout, sharex=True)
            if isinstance(axes_temp, plt.Axes):
                axes_temp = np.array([axes_temp])
            column_length_response  = len(axes_temp)

            self.axes:list[plt.Axes]      = []

            if isinstance(axes_temp[0], list|np.ndarray):
                for j in range(len(axes_temp[0])):
                    for i in range(len(axes_temp)): self.axes.append(axes_temp[i][j])
            else: self.axes = axes_temp.copy().tolist()
            del axes_temp

            for i in range(column_length_response-1, len(self.axes), column_length_response):
                self.axes[i].set_xlabel("Time (s)")

            self.fig.subplots_adjust(hspace=0.0)

            for ax in self.axes:
                ax.minorticks_on()
                ax.grid(True)

            for i in range(len(self.vars)):
                self.axes[i].set_ylabel(f"$\\{self.vars[i]}$" if self.vars[i] in self.special_words else f"${self.vars[i]}$")

            self.init = False
            
            self.fig.suptitle("Dynamic Responses")
            self.fig.supylabel("")
            self.fig.set_figheight(height)
            self.fig.set_figwidth(width)
            self.fig.tight_layout()
            self.fig.subplots_adjust(hspace=0.06)
        else:
            print(" not rendering plots")

        self.t = []
        lengths = self.vars_lengths
        for i, var in enumerate(self.vars):
            if lengths[i]:
                for j in range(lengths[i]):
                    exec(f"self.{var}_{self.dimension_names[j]} = []")
            else: exec(f"self.{var} = []")


    def createPlot(self):
        for i, var in enumerate(self.vars):
            if self.vars_lengths[i]:
                for j in range(self.vars_lengths[i]):
                    exec(f"self.line_{var}_{self.dimension_names[j]} = Line2D(self.t, self.{var}_{self.dimension_names[j]}, color=\'{self.colors[j]}\')")
                    exec(f"self.axes[{i}].add_line(self.line_{var}_{self.dimension_names[j]})")

            else:
                exec(f"self.line_{var} = Line2D(self.t, self.{var}, color=\'{self.colors[0]}\')")
                exec(f"self.axes[{i}].add_line(self.line_{var})")
        self.init = True


    def updatePlot(self):
        lengths = self.vars_lengths
        for i, var in enumerate(self.vars):
            if lengths[i]:
                for j in range(lengths[i]):
                    exec(f"self.line_{var}_{self.dimension_names[j]}.set_data(self.t, self.{var}_{self.dimension_names[j]})")
            else: exec(f"self.line_{var}.set_data(self.t, self.{var})")


    def orientePlot(self):
        for ax in self.axes:
            ax.relim()
            ax.autoscale()


    def update(self, t, *dynamic_response):
        self.t.append(t)
                
        for i, var in enumerate(self.vars):
            if self.vars_lengths[i]:
                for j in range(self.vars_lengths[i]):
                    exec(f"self.{var}_{self.dimension_names[j]}.append(inputs[{i}][{j}])")
            else:
                exec(f"self.{var}.append(dynamic_response[{i}])")


        if self.plot:
            if not self.init: self.createPlot()
            else: self.updatePlot()
            self.orientePlot()


    def show(self): self.fig.show()


    def savefig(self, file, *args, **kwargs):
        if self.plot: self.fig.savefig(file, *args, **kwargs)


    def saveData(self, file_for_data):
        data = {}
        data["t"] = self.t

        temp = {}
        for i, var in enumerate(self.vars):
            if self.vars_lengths[i]:
                temp[var] = {}
                for j in range(self.vars_lengths[i]):
                    temp[var][self.dimension_names[j]] = eval(f"self.{var}_{self.dimension_names[j]}")
                    if isinstance(temp[var][self.dimension_names[j]][0], np.ndarray):
                        for i in range(len(temp[var][self.dimension_names[j]])):
                            temp[var][self.dimension_names[j]][i] = temp[var][self.dimension_names[j]][i][0]
            else: temp[var] = eval(f"self.{var}")

        data["vars"] = temp.copy()

        with open(file_for_data, "w") as write_file: json.dump(data, write_file)