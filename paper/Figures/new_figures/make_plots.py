import mesh
import tracking
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
from os import path
from os.path import dirname
import numpy as np

def make_rearrangement_figure():

    rearrangement_data = np.loadtxt("rearrangement_statistics.csv")

    figuresize = (4.35,2.75)
    this_figure = plt.figure( figsize = figuresize )
#     ax = this_figure.add_axes([0.12, 0.15, 0.52, 0.75])
    plt.plot( rearrangement_data )
    plt.xlabel("tracking step")
    plt.ylabel("number of rearrangements")
    
#     plt.legend(loc = "center left", bbox_to_anchor=(1, 0.85))
    plt.tight_layout()
    plt.savefig("rearrangement_statistics.pdf")

def make_tracking_data_figure():

    tracking_data = np.loadtxt("tracking_statistics.csv")

    figuresize = (4.35,2.75)
    this_figure = plt.figure( figsize = figuresize )
#     ax = this_figure.add_axes([0.12, 0.15, 0.52, 0.75])
    plt.plot( tracking_data[:,0] )
    plt.plot( tracking_data[:,1] )
    plt.xlabel("tracking step")
    plt.ylabel("number of cells")
    
#     plt.legend(loc = "center left", bbox_to_anchor=(1, 0.85))
    plt.tight_layout()
    plt.savefig("tracking_statistics.pdf")

def make_area_data_figure():

    area_data = np.loadtxt("area_statistics.csv")

    figuresize = (4.35,2.75)
    this_figure = plt.figure( figsize = figuresize )
#     ax = this_figure.add_axes([0.12, 0.15, 0.52, 0.75])
    plt.plot( area_data )
    plt.xlabel("tracking step")
    plt.ylabel("cell area")
    
#     plt.legend(loc = "center left", bbox_to_anchor=(1, 0.85))
    plt.tight_layout()
    plt.savefig("area_statistics.pdf")
    
def make_tracking_ratio_statistics():
    tracking_data = np.loadtxt("tracking_statistics.csv")

    new_data = np.divide( tracking_data[:,1].astype('float'), tracking_data[:, 0])
    figuresize = (4.35,2.75)
    this_figure = plt.figure( figsize = figuresize )
#     ax = this_figure.add_axes([0.12, 0.15, 0.52, 0.75])
    plt.plot( new_data )
    plt.xlabel("tracking step")
    plt.ylabel("number of cells")
    
#     plt.legend(loc = "center left", bbox_to_anchor=(1, 0.85))
    plt.tight_layout()
    plt.savefig("tracking_ratio_statistics.pdf")

def make_cell_growth_plot():
    area_data = np.loadtxt("cell_area_statistics.csv")
    figuresize = (4.35,2.75)
    this_figure = plt.figure( figsize = figuresize )
    area_ratio = []
    for trajectory in area_data:
            length = np.count_nonzero(~np.isnan( trajectory ))
            area_ratio.append( np.nanmax(trajectory)/np.nanmin(trajectory) )
            if length < 5:
#                 plt.plot(trajectory, color = 'black', alpha = 0.2, zorder = 4)
                pass
            elif length < 10:
#                 plt.plot(trajectory, color = 'red', alpha = 0.1, zorder = 3)
                pass
            elif length < 15:
#                 plt.plot(trajectory, color = 'blue', alpha = 0.1, zorder = 2)
                pass
            else:
                plt.plot(trajectory, color = 'green', alpha = 0.1, zorder = 1)
    plt.xlabel("tracking step")
    plt.ylabel("cell area")
    plt.ylim(0,3500)
    plt.tight_layout()
    plt.savefig("cell_area_data.pdf")
    
    print 'max area ratio'
    print np.nanmax(area_ratio)
    print 'min area_ratio'
    print np.nanmin(area_ratio)
    print 'average area ratio'
    print np.nanmean(area_ratio)

def make_dying_cells_area_change_plot():
    area_data = np.loadtxt("cell_area_statistics.csv")
    death_ids = np.loadtxt("dying_cells.csv")
    figuresize = (4.35,2.75)
    this_figure = plt.figure( figsize = figuresize )
    area_ratio = []
    for death_id in death_ids: 
        trajectory = area_data[death_id]
        plt.plot(trajectory, color = 'black', alpha = 0.2, zorder = 4)
    plt.xlabel("tracking step")
    plt.ylabel("cell area")
#     plt.ylim(0,3500)
    plt.tight_layout()
    plt.savefig("cell_death_data.pdf")
 
if __name__ == "__main__":
    make_rearrangement_figure()
    make_tracking_data_figure()
    make_area_data_figure()
    make_tracking_ratio_statistics()
    make_cell_growth_plot()
    make_dying_cells_area_change_plot()
    
