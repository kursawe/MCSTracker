# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

import mesh
import tracking
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys

def make_cell_network_figure():
    """We load a mesh where two maximum common subgraphs will be identified, and track the division event"""

    mesh_one = mesh.creation.generate_random_tesselation( 9, 9 )
    
    mesh_one.assign_frame_ids_in_order()
    
    first_polygon_collection = mesh_one.get_polygon_collection( color_by_global_id = 'white' )
    first_polygon_collection.set_edgecolor('grey')
#     first_polygon_collection.set_dashes([(0,[2,1])])
    first_polygon_collection.set_linewidth(0.2)
    all_edges_by_elements = mesh_one.collect_elements_of_inner_edges()

    figuresize = (4,2.75)
    figure_1 = plt.figure(figsize = figuresize)
    for edge_elements in all_edges_by_elements:
        first_centroid = mesh_one.get_element_with_frame_id(edge_elements.pop()).calculate_centroid()
        second_centroid = mesh_one.get_element_with_frame_id(edge_elements.pop()).calculate_centroid()
        plt.plot( [ first_centroid[0], second_centroid[0]],
                  [ first_centroid[1], second_centroid[1]], color = 'black', solid_capstyle = 'round' )
    figure_1.gca().add_collection(first_polygon_collection)
    figure_1.gca().set_aspect('equal')
    figure_1.gca().autoscale_view()
    plt.axis('off')
    this_filename = 'cell_network.pdf'
    figure_1.savefig(this_filename, bbox_inches = 'tight')
    plt.close(figure_1)

if __name__ == "__main__":

    make_cell_network_figure()