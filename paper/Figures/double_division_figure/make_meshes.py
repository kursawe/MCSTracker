# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

import mesh
import tracking
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys

def make_double_division_plot():
    """We load a mesh where two maximum common subgraphs will be identified, and track the division event"""

    mesh_one = mesh.creation.generate_random_tesselation(9,9)
    mesh_two = copy.deepcopy(mesh_one)

    mesh_one.assign_frame_ids_in_order()
    mesh_two.assign_frame_ids_randomly()

    # build ground truth for testing the mapping
    ground_truth = {}
    for element_index, element in enumerate(mesh_one.elements):
        ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame
 
    # pick the element closest to the centre
    most_central_element = mesh_two.find_most_central_element()
    
    # get adjacent element ids
    adjacent_elements = most_central_element.get_ids_of_adjacent_elements()
    
    # pick one
    
    second_dividing_cell = adjacent_elements[0]
            
    mesh_two.divide_element_with_frame_id(most_central_element.id_in_frame)
    mesh_two.divide_element_with_frame_id(second_dividing_cell)
    
    tracked_ids = tracking.track(mesh_one, mesh_two)
     
    plotname = 'double_division'

    figuresize = (4,2.75)

    figure_1 = plt.figure(figsize = figuresize)
    first_polygon_collection = mesh_one.get_polygon_collection( color_by_global_id = True, 
                               total_number_of_global_ids = len( tracked_ids ) )
    figure_1.gca().add_collection(first_polygon_collection)
    figure_1.gca().set_aspect('equal')
    figure_1.gca().autoscale_view()
    plt.axis('off')
    this_filename = plotname + '_before_all_colors.pdf'
    figure_1.savefig(this_filename, bbox_inches = 'tight')
    plt.close(figure_1)

    figure_2 = plt.figure(figsize = figuresize)
    second_polygon_collection = mesh_two.get_polygon_collection( color_by_global_id = True, 
                               total_number_of_global_ids = len( tracked_ids ) )
    figure_2.gca().add_collection(second_polygon_collection)
    figure_2.gca().set_aspect('equal')
    figure_2.gca().autoscale_view()
    plt.axis('off')
    this_filename = plotname + '_after_all_colors.pdf'
    figure_2.savefig(this_filename, bbox_inches = 'tight')
    plt.close(figure_2)

    figure_3 = plt.figure(figsize = figuresize)
    third_polygon_collection = mesh_one.get_polygon_collection( color_by_global_id = True, 
                               total_number_of_global_ids = len( tracked_ids ) )
    polygon_list = []
    for element in mesh_one.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
        if element.global_id == None:
            this_polygon.set_facecolor([1.0, 1.0, 1.0])
        else:
            if element.id_in_frame in ground_truth:
                second_element = mesh_two.get_element_with_global_id(element.global_id).id_in_frame
                if ground_truth[element.id_in_frame] == second_element:
                    this_polygon.set_facecolor('#2dec34')
                else: 
                    this_polygon.set_facecolor('red')
            else:
                this_polygon.set_facecolor('red')
        polygon_list.append(this_polygon)
    third_polygon_collection = mpl.collections.PatchCollection(polygon_list, match_original = True)
    figure_3.gca().add_collection(third_polygon_collection)
    figure_3.gca().set_aspect('equal')
    figure_3.gca().autoscale_view()
    plt.axis('off')
    this_filename = plotname + '_before.pdf'
    figure_3.savefig(this_filename, bbox_inches = 'tight')
    plt.close(figure_3)

    figure_4 = plt.figure(figsize = figuresize)
    fourth_polygon_collection = mesh_two.get_polygon_collection( color_by_global_id = True, 
                               total_number_of_global_ids = len( tracked_ids ) )
    polygon_list = []
    for element in mesh_two.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
        if element.global_id == None:
            this_polygon.set_facecolor([1.0, 1.0, 1.0])
        else:
            if element.id_in_frame in ground_truth.values():
                first_element = mesh_one.get_element_with_global_id(element.global_id).id_in_frame
                if first_element in ground_truth:
                    if ground_truth[first_element] == element.id_in_frame:
                        this_polygon.set_facecolor('#2dec34')
                    else: 
                        this_polygon.set_facecolor('red')
                else: 
                    this_polygon.set_facecolor('red')
            else:
                this_polygon.set_facecolor('red')
        polygon_list.append(this_polygon)
    fourth_polygon_collection = mpl.collections.PatchCollection(polygon_list, match_original = True)
    figure_4.gca().add_collection( fourth_polygon_collection )
    figure_4.gca().set_aspect('equal')
    figure_4.gca().autoscale_view()
    plt.axis('off')
    this_filename = plotname + '_after.pdf'
    figure_4.savefig(this_filename, bbox_inches = 'tight')
    plt.close(figure_4)

if __name__ == "__main__":

    make_double_division_plot()