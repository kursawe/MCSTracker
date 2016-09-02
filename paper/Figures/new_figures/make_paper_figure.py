import mesh
import tracking
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
from os import path
from os.path import dirname
import numpy as np

def make_paper_figure():
    
    
    figuresize = (6.5,5)
    this_figure = plt.figure( figsize = figuresize )

    this_figure.add_subplot(221)
    tracking_data = np.loadtxt("tracking_statistics.csv")
    
    x_values = np.arange(0, len(tracking_data[:,0]), 1)
    x_values *= 5
    plt.plot(x_values, tracking_data[:,0], color = 'lightgray', label = 'total' )
    plt.plot(x_values, tracking_data[:,1], color = 'darkslategray', label = 'tracked' )
    plt.legend(loc = "upper left" )
    plt.setp(plt.gca().get_legend().get_texts(), fontsize = '10')
    plt.xlabel("Time (min)")
    plt.ylabel("Number of cells")

    this_figure.add_subplot(222)
    rearrangement_data = np.loadtxt("rearrangement_statistics.csv")
    plt.plot( x_values, rearrangement_data, color = 'darkslategray' )
    plt.xlabel("Time (min)")
    plt.ylabel("Cell rearrangements (/5 min)")
    
    this_figure.add_subplot(223)
    area_data = np.loadtxt("area_statistics.csv")
    area_x_values = np.arange(0, len(area_data), 1)
    area_x_values *= 5
    area_data = np.loadtxt("area_statistics.csv")
    plt.plot( area_x_values, area_data, color = 'darkslategray' )
    plt.xlabel("Time (min)")
    plt.ylabel("Average cell area")
    
    this_figure.add_subplot(224)

    tracking_data = np.loadtxt("tracking_statistics.csv")

    new_data = np.divide( tracking_data[:,1].astype('float'), tracking_data[:, 0])
#     ax = this_figure.add_axes([0.12, 0.15, 0.52, 0.75])
    plt.plot( x_values, new_data, color = 'darkslategray' )
    plt.xlabel("Time (min) ")
    plt.ylabel("Percentage of tracked cells")
 
#     plt.legend(loc = "center left", bbox_to_anchor=(1, 0.85))

    this_figure.tight_layout()
    plt.savefig("tracking_figure.pdf")

if __name__ == "__main__":
    # font = {'size'   : 10,
    #         'sans-serif' : 'Helvetica'}
    font = {'size'   : 10}
    plt.rc('font', **font)
    mpl.rcParams['mathtext.default'] = 'regular'

    make_paper_figure()
