import mesh
import tracking
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
from os import path
from os.path import dirname
import numpy as np

def make_t1_swap_figure():

    correctly_tracked_cells_data = np.loadtxt("correctly_tracked_cells.csv")
    incorrectly_tracked_cells_data = np.loadtxt("incorrectly_tracked_cells.csv")
    success_ratio_data = np.loadtxt("success_ratio.csv")

    correctly_tracked_cells_percentage = correctly_tracked_cells_data[:,0]
    incorrectly_tracked_cells_percentage = incorrectly_tracked_cells_data[:,0]
    success_ratio_percentage = success_ratio_data[:,0]

    correctly_tracked_cells_values = correctly_tracked_cells_data[:,1:]
    incorrectly_tracked_cells_values = incorrectly_tracked_cells_data[:,1:]
    success_ratios = success_ratio_data[:,1]

    font = {'size'   : 10}
    plt.rc('font', **font)
    mpl.rcParams['mathtext.default'] = 'regular'
    
    figuresize = (4.35,2.75)
    this_figure = plt.figure( figsize = figuresize )
#     ax = this_figure.add_axes([0.12, 0.15, 0.52, 0.75])
    plt.errorbar(correctly_tracked_cells_percentage, 
             np.mean(correctly_tracked_cells_values, axis = 1)*100,
             np.std(correctly_tracked_cells_values, axis = 1)*100/np.sqrt(10.0),
             ls = '--', dashes = (7,2),
             color = 'blue', label = 'Correctly tracked')
    for counter, percentage in enumerate(correctly_tracked_cells_percentage):
        data_points = correctly_tracked_cells_values[counter]*100
        x_values = np.array([percentage + 0.1]*len(data_points))
        plt.scatter(x_values, data_points, marker = '.', s = 30, 
                    alpha = 0.3, lw =0, color = 'blue')
    plt.errorbar(incorrectly_tracked_cells_percentage, 
             np.mean(incorrectly_tracked_cells_values, axis = 1)*100,
             np.std(incorrectly_tracked_cells_values, axis = 1)*100/np.sqrt(10.0),
             color = 'red', label = 'Incorrectly tracked')
    for counter, percentage in enumerate(incorrectly_tracked_cells_percentage):
        data_points = incorrectly_tracked_cells_values[counter]*100
        x_values = np.array([percentage - 0.1]*len(data_points))
        plt.scatter(x_values, data_points, marker = '.', s = 20, 
                    alpha = 0.3, lw=0, color = 'red')

#     plt.plot(success_ratio_percentage, success_ratios)
    
    plt.xlim(0,15)
    plt.ylim(0,100)
    
    plt.xlabel("Percentage of edges undergoing T1 transitions")
    plt.ylabel("Percentage of cells")
    
#     plt.legend(loc = "center left", bbox_to_anchor=(1, 0.85))
    this_legend = plt.legend(loc = "upper right")
    this_legend.set_zorder(0)
    plt.setp(plt.gca().get_legend().get_texts(), fontsize = '10')
    plt.tight_layout()
    plt.savefig("t1_swap_analysis_full.pdf")

if __name__ == "__main__":
    make_t1_swap_figure()
    
