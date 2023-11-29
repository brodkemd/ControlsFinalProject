import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from collections import OrderedDict
from math import degrees
import numpy as np
import json, re

plt.ion()  # enable interactive drawing

def wrapString(s):
    return f"\'{s}\'"

class DataPlotter:
    special_words       = ["phi", "theta", "psi", "alpha", "beta", "delta", "Omega"]

    inputs              = ["f_E_x","f_E_y","f_E_z", "f_cp_x", "f_cp_y", "f_cp_z", "tau_cp_x", "tau_cp_y", "tau_cp_z"]
    inputs_lengths      = np.zeros(len(inputs)).tolist()

    vars                = ["p_n", "p_e", "p_d", "u", "v", "w", "e_0", "e_1", "e_2", "e_3", "p", "q", "r"]
    vars_lengths        = (2*np.ones(len(vars), dtype=int)).tolist()

    layout           = (5, 3)
    dashboard_layout = (3, 3)

    plot = True
    dimension_names = ["Actual", "Reference"]#, "Down"]
    colors          = ["tab:blue", "tab:orange"]#, "tab:green"]


    def __init__(self, width=10, height=8, plot=True, interactive=True) -> None:
        self.plot = plot
        print("\nPlotter:")
        print("  Rendering Setup:")
        print("    interactive:", interactive)
        print("    plotting:",    plot)

        if self.plot:
            self.fig, axes_temp           = plt.subplots(*self.layout, sharex=True)
            self.fig_dashboard, dashboard = plt.subplots(*self.dashboard_layout, sharex=True)

            column_length_response  = len(axes_temp)
            column_length_dashboard = len(dashboard)

            self.axes:list[plt.Axes]      = []
            self.dashboard:list[plt.Axes] = []

            if isinstance(axes_temp[0], np.ndarray) or isinstance(axes_temp[0], list):
                for j in range(len(axes_temp[0])):
                    for i in range(len(axes_temp)): self.axes.append(axes_temp[i][j])
            else: self.axes = axes_temp.copy().tolist()
            
            if isinstance(dashboard[0], np.ndarray) or isinstance(dashboard[0], list):
                for j in range(len(dashboard)):
                    for i in range(len(dashboard[0])): self.dashboard.append(dashboard[j][i])
            else: self.dashboard = dashboard.copy().tolist()

            del axes_temp, dashboard
           

            for i in range(column_length_response-1, len(self.axes), column_length_response):
                self.axes[i].set_xlabel("Time (s)")

            for i in range(len(self.dashboard)//column_length_dashboard):
                self.dashboard[-1-i].set_xlabel("Time (s)")

            self.fig.subplots_adjust(hspace=0.0)

            for ax in self.axes+self.dashboard:
                ax.minorticks_on()
                ax.grid(True)

            for i in range(len(self.vars)):
                # if self.vars_lengths[i]:
                #     plt.legend(self.dimension_names[:self.vars_lengths[i]])
                self.axes[i].set_ylabel(self.formatStringForLatex(self.vars[i]))

            for i in range(len(self.inputs)):
                # if self.inputs_lengths[i]:
                #     plt.legend(self.dimension_names[:self.inputs_lengths[i]])
                self.dashboard[i].set_ylabel(self.formatStringForLatex(self.inputs[i]))

            # Create a figure and axis
            if sum(self.inputs_lengths + self.vars_lengths):
                patches = [mpatches.Patch(color=color, label=label, linestyle="-") for label, color in zip(self.dimension_names, self.colors)]
                self.fig_dashboard.legend(patches, self.dimension_names, ncols=len(self.dimension_names), frameon=False)

            # plt.subplots_adjust(hspace=0.0)
            # self.fig_dashboard.subplots_adjust(hspace=0)
            self.init = False
            
            self.fig.suptitle("Dynamic Responses")
            self.fig_dashboard.suptitle("Inputs from user and Flight Conditions")
            self.fig.supylabel("")
            self.fig_dashboard.supylabel("")

            self.fig_dashboard.set_figheight(height/3)
            self.fig_dashboard.set_figwidth(width)
            self.fig_dashboard.tight_layout()
            self.fig_dashboard.subplots_adjust(hspace=0.08)

            self.fig.set_figheight(height+self.fig_dashboard.get_figheight())
            self.fig.set_figwidth(width)
            self.fig.tight_layout()
            self.fig.subplots_adjust(hspace=0.06)

            self.fig.canvas.mpl_connect("close_event", self.onClose)
            self.fig_dashboard.canvas.mpl_connect("close_event", self.onClose)


        for i in range(len(self.vars)):
            self.vars[i] = re.sub(r"\s+", "_", self.vars[i])

        for i in range(len(self.inputs)):
            self.inputs[i] = re.sub(r"\s+", "_", self.inputs[i])

        self.t = []
        lengths = self.vars_lengths+self.inputs_lengths
        for i, var in enumerate(self.vars+self.inputs):
            if lengths[i]:
                for j in range(lengths[i]):
                    exec(f"self.{var}_{self.dimension_names[j]} = []")
            else:
                exec(f"self.{var} = []")


    def onClose(self, event):
        plt.close(self.fig)
        plt.close(self.fig_dashboard)


    def formatStringForLatex(self, s:str):
        return s
        # components = re.split(r"\s+", s)
        # for j, item in enumerate(components):
        #     underscores = re.split(r"_+", item)
        #     for i, name in enumerate(underscores):
        #         if name in self.special_words:
        #             underscores[i] = f"\\{name}"
        #         else:
        #             underscores[i] = f"\\text{{{name}}}"
        #     components[j] = "_{".join(underscores) + "}"*(len(underscores)-1)
        # return "$" + "\;".join(components) + "$"


    def createPlot(self):
        for i, var in enumerate(self.vars):
            if self.vars_lengths[i]:
                for j in range(self.vars_lengths[i]):
                    exec(f"self.line_{var}_{self.dimension_names[j]} = Line2D(self.t, self.{var}_{self.dimension_names[j]}, color=\'{self.colors[j]}\')")
                    exec(f"self.axes[{i}].add_line(self.line_{var}_{self.dimension_names[j]})")

            else:
                exec(f"self.line_{var} = Line2D(self.t, self.{var}, color=\'{self.colors[0]}\')")
                exec(f"self.axes[{i}].add_line(self.line_{var})")
        
        for i, var in enumerate(self.inputs):
            if self.inputs_lengths[i]:
                for j in range(self.inputs_lengths[i]):
                    exec(f"self.line_{var}_{self.dimension_names[j]} = Line2D(self.t, self.{var}_{self.dimension_names[j]}, color=\'{self.colors[j]}\')")
                    exec(f"self.dashboard[{i}].add_line(self.line_{var}_{self.dimension_names[j]})")

            else:
                exec(f"self.line_{var} = Line2D(self.t, self.{var}, color=\'{self.colors[0]}\')")
                exec(f"self.dashboard[{i}].add_line(self.line_{var})")

        self.init = True


    def updatePlot(self):
        lengths = self.vars_lengths+self.inputs_lengths
        for i, var in enumerate(self.vars+self.inputs):
            if lengths[i]:
                for j in range(lengths[i]):
                    exec(f"self.line_{var}_{self.dimension_names[j]}.set_data(self.t, self.{var}_{self.dimension_names[j]})")
            else:
                exec(f"self.line_{var}.set_data(self.t, self.{var})")


    def orientePlot(self):
        for ax in self.axes+self.dashboard:
            #self.axes[i].set_xlim(max([self.t[-1]-self.time_to_display, 0]), max([self.t[-1]+self.time_to_display, self.time_to_display]))
            ax.relim()
            ax.autoscale()


    def update(self, t, dynamic_response, inputs):
        self.t.append(t)
                
        for i, var in enumerate(self.vars):
            if self.vars_lengths[i]:
                for j in range(self.vars_lengths[i]):
                    exec(f"self.{var}_{self.dimension_names[j]}.append(dynamic_response[{i}][{j}])")
            else:
                exec(f"self.{var}.append(dynamic_response[{i}])")

        for i, var in enumerate(self.inputs):
            if self.inputs_lengths[i]:
                for j in range(self.inputs_lengths[i]):
                    exec(f"self.{var}_{self.dimension_names[j]}.append(inputs[{i}][{j}])")
            else:
                exec(f"self.{var}.append(inputs[{i}])")

        if self.plot:
            if not self.init:
                self.createPlot()
                
            else:
                self.updatePlot()

            self.orientePlot()


    def show(self): self.fig.show()


    def savefig(self, image_file_for_outputs_plot=None, image_file_for_inputs_plot=None, *args, **kwargs):
        if image_file_for_inputs_plot is not None and self.plot:
            self.fig_dashboard.savefig(image_file_for_inputs_plot, *args, **kwargs)
        
        if image_file_for_outputs_plot is not None and self.plot:
            self.fig.savefig(image_file_for_outputs_plot, *args, **kwargs)
    

    def saveData(self, file_for_data):
        data = {}
        data["t"] = self.t

        temp = {}
        for i, var in enumerate(self.inputs):
            if self.inputs_lengths[i]:
                temp[var] = {}
                for j in range(self.inputs_lengths[i]):
                    temp[var][self.dimension_names[j]] = eval(f"self.{var}_{self.dimension_names[j]}")
                    if isinstance(temp[var][self.dimension_names[j]][0], np.ndarray):
                        for i in range(len(temp[var][self.dimension_names[j]])):
                            temp[var][self.dimension_names[j]][i] = temp[var][self.dimension_names[j]][i][0]

            else:
                temp[var] = eval(f"self.{var}")

        data["inputs"] = temp.copy()
        temp.clear()
        
        for i, var in enumerate(self.vars):
            if self.vars_lengths[i]:
                temp[var] = {}
                for j in range(self.vars_lengths[i]):
                    temp[var][self.dimension_names[j]] = eval(f"self.{var}_{self.dimension_names[j]}")
                    if isinstance(temp[var][self.dimension_names[j]][0], np.ndarray):
                        for i in range(len(temp[var][self.dimension_names[j]])):
                            temp[var][self.dimension_names[j]][i] = temp[var][self.dimension_names[j]][i][0]
                    
            else:
                temp[var] = eval(f"self.{var}")

        data["vars"] = temp.copy()
        temp.clear()

        # for item in data:
        #     if isinstance(data[item], dict|OrderedDict):
        #         print(item+":")
        #         for val in data[item]:
        #             print("-->", val, type(data[item][val]))
        #     else:
        #         print("->",item, type(data[item]))

        with open(file_for_data, "w") as write_file:
            json.dump(data, write_file)

# import matplotlib.pyplot as plt
# from matplotlib.lines import Line2D
# import matplotlib.patches as mpatches
# import numpy as np
# import json

# plt.ion()  # enable interactive drawing

# # def wrapString(s):
# #    return f"\'{s}\'"

# class DataPlotter:
#     special_words       = ["phi", "theta", "psi", "alpha", "beta", "delta"]

#     vars                = ["x", "y"]
#     vars_lengths        = np.zeros(len(vars), dtype=int).tolist()
#     layout              = (2, 1)

#     plot = True
#     dimension_names = ["North", "East", "Down"]
#     colors          = ["tab:blue", "tab:orange", "tab:green"]

#     def __init__(self, width=10, height=8, plot=True) -> None:
#         print("\nData Plotter:", end="")
#         self.plot = plot

#         if self.plot:
#             print()
#             print("  Creating plot grid:", f"{self.layout[0]}x{self.layout[1]}")
#             self.fig, axes_temp           = plt.subplots(*self.layout, sharex=True)
#             if isinstance(axes_temp, plt.Axes):
#                 axes_temp = np.array([axes_temp])
#             column_length_response  = len(axes_temp)

#             self.axes:list[plt.Axes]      = []

#             if isinstance(axes_temp[0], list|np.ndarray):
#                 for j in range(len(axes_temp[0])):
#                     for i in range(len(axes_temp)): self.axes.append(axes_temp[i][j])
#             else: self.axes = axes_temp.copy().tolist()
#             del axes_temp

#             for i in range(column_length_response-1, len(self.axes), column_length_response):
#                 self.axes[i].set_xlabel("Time (s)")

#             self.fig.subplots_adjust(hspace=0.0)

#             for ax in self.axes:
#                 ax.minorticks_on()
#                 ax.grid(True)

#             for i in range(len(self.vars)):
#                 self.axes[i].set_ylabel(f"$\\{self.vars[i]}$" if self.vars[i] in self.special_words else f"${self.vars[i]}$")

#             self.init = False
            
#             self.fig.suptitle("Dynamic Responses")
#             self.fig.supylabel("")
#             self.fig.set_figheight(height)
#             self.fig.set_figwidth(width)
#             self.fig.tight_layout()
#             self.fig.subplots_adjust(hspace=0.06)
#         else:
#             print(" not rendering plots")

#         self.t = []
#         lengths = self.vars_lengths
#         for i, var in enumerate(self.vars):
#             if lengths[i]:
#                 for j in range(lengths[i]):
#                     exec(f"self.{var}_{self.dimension_names[j]} = []")
#             else: exec(f"self.{var} = []")


#     def createPlot(self):
#         for i, var in enumerate(self.vars):
#             if self.vars_lengths[i]:
#                 for j in range(self.vars_lengths[i]):
#                     exec(f"self.line_{var}_{self.dimension_names[j]} = Line2D(self.t, self.{var}_{self.dimension_names[j]}, color=\'{self.colors[j]}\')")
#                     exec(f"self.axes[{i}].add_line(self.line_{var}_{self.dimension_names[j]})")

#             else:
#                 exec(f"self.line_{var} = Line2D(self.t, self.{var}, color=\'{self.colors[0]}\')")
#                 exec(f"self.axes[{i}].add_line(self.line_{var})")
#         self.init = True


#     def updatePlot(self):
#         lengths = self.vars_lengths
#         for i, var in enumerate(self.vars):
#             if lengths[i]:
#                 for j in range(lengths[i]):
#                     exec(f"self.line_{var}_{self.dimension_names[j]}.set_data(self.t, self.{var}_{self.dimension_names[j]})")
#             else: exec(f"self.line_{var}.set_data(self.t, self.{var})")


#     def orientePlot(self):
#         for ax in self.axes:
#             ax.relim()
#             ax.autoscale()


#     def update(self, t, *dynamic_response):
#         self.t.append(t)
                
#         for i, var in enumerate(self.vars):
#             if self.vars_lengths[i]:
#                 for j in range(self.vars_lengths[i]):
#                     exec(f"self.{var}_{self.dimension_names[j]}.append(inputs[{i}][{j}])")
#             else:
#                 exec(f"self.{var}.append(dynamic_response[{i}])")


#         if self.plot:
#             if not self.init: self.createPlot()
#             else: self.updatePlot()
#             self.orientePlot()


#     def show(self): self.fig.show()


#     def savefig(self, file, *args, **kwargs):
#         if self.plot: self.fig.savefig(file, *args, **kwargs)


#     def saveData(self, file_for_data):
#         data = {}
#         data["t"] = self.t

#         temp = {}
#         for i, var in enumerate(self.vars):
#             if self.vars_lengths[i]:
#                 temp[var] = {}
#                 for j in range(self.vars_lengths[i]):
#                     temp[var][self.dimension_names[j]] = eval(f"self.{var}_{self.dimension_names[j]}")
#                     if isinstance(temp[var][self.dimension_names[j]][0], np.ndarray):
#                         for i in range(len(temp[var][self.dimension_names[j]])):
#                             temp[var][self.dimension_names[j]][i] = temp[var][self.dimension_names[j]][i][0]
#             else: temp[var] = eval(f"self.{var}")

#         data["vars"] = temp.copy()

#         with open(file_for_data, "w") as write_file: json.dump(data, write_file)