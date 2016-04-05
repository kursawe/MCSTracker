# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/LICENSE.

import mesh
import tracking
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys

def make_death_plots(figuresize):
    mesh_one = mesh.creation.generate_random_tesselation(9,9)
    mesh_two = copy.deepcopy(mesh_one)

    # First pick a cell in the centre
    mesh_centre = mesh_two.calculate_centre()
        
    mesh_one.assign_frame_ids_in_order()
    mesh_two.assign_frame_ids_randomly()

    most_central_element = mesh_two.find_most_central_element()
               
    mesh_two.kill_element_with_frame_id(most_central_element.id_in_frame)

    tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )

    for element in mesh_one.elements:
        element.global_id = None
  
    for element in mesh_two.elements:
        element.global_id = None
  
    tracked_ids_all = tracking.track( mesh_one, mesh_two )
    
    plotname = 'cell_death'
    
    make_difference_plot(mesh_one, mesh_two, tracked_ids, tracked_ids_all, plotname, figuresize)

def make_difference_plot(mesh_one, mesh_two, tracked_ids, tracked_ids_all, plotname, figuresize):

    # plot before cell death with white filled in
    # figure properties
    mesh_figure = plt.figure(figsize = figuresize)
    polygon_list = []
    for element in mesh_one.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
        if element.global_id == None:
           this_polygon.set_facecolor([1.0, 1.0, 1.0])
        elif (element.global_id in tracked_ids) and (element.global_id in tracked_ids_all):
           this_polygon.set_facecolor('g')
        elif (element.global_id in tracked_ids_all):
           this_polygon.set_facecolor('r')
        polygon_list.append(this_polygon)
    polygon_collection = mpl.collections.PatchCollection(polygon_list, match_original = True)
    mesh_figure.gca().add_collection(polygon_collection)
    mesh_figure.gca().set_aspect('equal')
    mesh_figure.gca().autoscale_view()
    plt.axis('off')
    this_filename = plotname + '_before.pdf'
    mesh_figure.savefig(this_filename, bbox_inches = 'tight')

    # plot before cell death with white filled in
    # figure properties
    mesh_figure = plt.figure(figsize = figuresize)
    polygon_list = []
    for element in mesh_two.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
        if element.global_id == None:
           this_polygon.set_facecolor([1.0, 1.0, 1.0])
        elif (element.global_id in tracked_ids) and (element.global_id in tracked_ids_all):
           this_polygon.set_facecolor('g')
        elif (element.global_id in tracked_ids_all):
           this_polygon.set_facecolor('r')
        polygon_list.append(this_polygon)
    polygon_collection = mpl.collections.PatchCollection(polygon_list, match_original = True)
    mesh_figure.gca().add_collection(polygon_collection)
    mesh_figure.gca().set_aspect('equal')
    mesh_figure.gca().autoscale_view()
    plt.axis('off')
    this_filename = plotname + '_after.pdf'
    mesh_figure.savefig(this_filename, bbox_inches = 'tight')

def make_t1_plots(figuresize):

    mesh_one = mesh.creation.generate_random_tesselation(9,9)
    mesh_two = copy.deepcopy(mesh_one)

    mesh_one.assign_frame_ids_in_order()
    mesh_two.assign_frame_ids_randomly()

    # build ground truth for testing the mapping
    ground_truth = {}
    for element_index, element in enumerate(mesh_one.elements):
        ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

    # Perform T1 swap on mesh two 
    # First pick a node in the centre
    mesh_centre = mesh_two.calculate_centre()

    # pick the node closest to the centre
    min_distance = 3*mesh_two.calculate_height()
    for node in mesh_two.nodes:
        distance = np.linalg.norm(node.position - mesh_centre)
        if distance < min_distance:
           min_distance = distance
           most_central_node = node 
               
    # pick a node that shares an edge with this central node
    for local_index, element_node in enumerate(most_central_node.adjacent_elements[0].nodes):
        if element_node.id == most_central_node.id:
            num_nodes_this_element = most_central_node.adjacent_elements[0].get_num_nodes()
            one_edge_node = most_central_node.adjacent_elements[0].nodes[(local_index+1)%num_nodes_this_element]
            break
        
    mesh_two.perform_t1_swap( most_central_node.id, one_edge_node.id )
         
    tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )

    for element in mesh_one.elements:
        element.global_id = None
  
    for element in mesh_two.elements:
        element.global_id = None
  
    tracked_ids_all = tracking.track( mesh_one, mesh_two )
    
    plotname = 't1_swap'
    
    make_difference_plot(mesh_one, mesh_two, tracked_ids, tracked_ids_all, plotname, figuresize)

def make_translation_plots(figuresize):
    sys.setrecursionlimit(40000)

    large_mesh = mesh.creation.generate_random_tesselation(15,8)
    
    large_mesh.assign_frame_ids_in_order()

    mesh_one = large_mesh.copy_region_of_interest((3,10),(-1,9))
    
    mesh_two = large_mesh.copy_region_of_interest((5,12),(-1,9))
    
    ground_truth_indices = {}
    for element_index, element in enumerate(mesh_one.elements):
        if element.id_in_frame in mesh_two.frame_id_dictionary:
            ground_truth_indices[element_index] = mesh_two.frame_id_dictionary[element.id_in_frame]
            
    mesh_one.assign_frame_ids_in_order()
    mesh_two.assign_frame_ids_randomly()

    ground_truth = {}
    for element_index, element in enumerate(mesh_one.elements):
        if element_index in ground_truth_indices:
            ground_truth[element.id_in_frame] = mesh_two.elements[ground_truth_indices[element_index]]
    
    tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
     
    for element in mesh_one.elements:
        element.global_id = None
   
    for element in mesh_two.elements:
        element.global_id = None
   
    mesh_one.index_global_ids()
    mesh_two.index_global_ids()

    tracked_ids_all = []
    tracked_ids_all = tracking.track( mesh_one, mesh_two )
    
    plotname = 'translation'
    make_difference_plot(mesh_one, mesh_two, tracked_ids, tracked_ids_all, plotname, figuresize)

def make_division_resolved_plot(figuresize):
    """We load a mesh where two maximum common subgraphs will be identified, and track the division event"""

    mesh_one = mesh.load('standard_ambiguous_division_one.mesh')
    mesh_two = copy.deepcopy(mesh_one)

    mesh_one.assign_frame_ids_in_order()
    mesh_two.assign_frame_ids_randomly()

    # build ground truth for testing the mapping
    ground_truth = {}
    for element_index, element in enumerate(mesh_one.elements):
        ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

    # pick the element closest to the centre
    most_central_element = mesh_two.find_most_central_element()
            
    mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, [1.0, 1.0])
     
    tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
     
    for element in mesh_one.elements:
        element.global_id = None
  
    for element in mesh_two.elements:
        element.global_id = None
 
    tracked_ids_all = tracking.track(mesh_one, mesh_two)        

    plotname = 'division_all'

    make_difference_plot(mesh_one, mesh_two, tracked_ids, tracked_ids_all, plotname, figuresize)

def make_division_I_plot(figuresize):

    """We load a mesh where two maximum common subgraphs will be identified, and track the division event"""

    mesh_one = mesh.load('standard_ambiguous_division_one.mesh')
    mesh_two = copy.deepcopy(mesh_one)

    mesh_one.assign_frame_ids_in_order()
    mesh_two.assign_frame_ids_randomly()

    # build ground truth for testing the mapping
    ground_truth = {}
    for element_index, element in enumerate(mesh_one.elements):
        ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

    # pick the element closest to the centre
    most_central_element = mesh_two.find_most_central_element()
            
    mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, [1.0, 1.0])

    subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
    subgraph_finder.find_maximum_common_subgraph()

    post_processor = tracking.PostProcessor(mesh_one, mesh_two, [subgraph_finder.largest_mappings[0]])
    post_processor.index_global_ids_from_largest_mappings()

    tracked_ids = post_processor.mapped_ids

    tracked_ids_all = tracked_ids
        
    plotname = 'division_I'

    make_difference_plot(mesh_one, mesh_two, tracked_ids, tracked_ids_all, plotname, figuresize)

def make_division_II_plot(figuresize):

    """We load a mesh where two maximum common subgraphs will be identified, and track the division event"""

    mesh_one = mesh.load('standard_ambiguous_division_one.mesh')
    mesh_two = copy.deepcopy(mesh_one)

    mesh_one.assign_frame_ids_in_order()
    mesh_two.assign_frame_ids_randomly()

    # build ground truth for testing the mapping
    ground_truth = {}
    for element_index, element in enumerate(mesh_one.elements):
        ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

    # pick the element closest to the centre
    most_central_element = mesh_two.find_most_central_element()
            
    mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, [1.0, 1.0])

    subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
    subgraph_finder.find_maximum_common_subgraph()

    post_processor = tracking.PostProcessor(mesh_one, mesh_two, [subgraph_finder.largest_mappings[1]])
    post_processor.index_global_ids_from_largest_mappings()

    tracked_ids = post_processor.mapped_ids

    tracked_ids_all = tracked_ids
        
    plotname = 'division_II'

    make_difference_plot(mesh_one, mesh_two, tracked_ids, tracked_ids_all, plotname, figuresize)

def make_general_division_plot(figuresize):
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
            
    mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, [1.0, 1.0])
     
    tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
    
    for element in mesh_one.elements:
        element.global_id = None
  
    for element in mesh_two.elements:
        element.global_id = None
 
    tracked_ids_all = tracking.track(mesh_one, mesh_two)        

    plotname = 'division_general'

    make_difference_plot(mesh_one, mesh_two, tracked_ids, tracked_ids_all, plotname, figuresize)

def make_permutation_plot(figuresize):
    mesh_one = mesh.creation.generate_random_tesselation(9,9)
    mesh_two = copy.deepcopy(mesh_one)

    mesh_one.assign_frame_ids_in_order()
    mesh_two.assign_frame_ids_randomly()

    # build ground truth for testing the mapping
    ground_truth = {}
    for element_index, element in enumerate(mesh_one.elements):
        ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

    tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
    
    for element in mesh_one.elements:
        element.global_id = None
  
    for element in mesh_two.elements:
        element.global_id = None
 
    tracked_ids_all = tracking.track(mesh_one, mesh_two)        

    plotname = 'permutation'

    make_difference_plot(mesh_one, mesh_two, tracked_ids, tracked_ids_all, plotname, figuresize)

if __name__ == "__main__":

    figuresize = (2,1.7)
    
    make_death_plots(figuresize)
    make_t1_plots(figuresize)
    make_translation_plots(figuresize)
    make_general_division_plot(figuresize)
#     make_division_resolved_plot(figuresize)
#     make_division_I_plot(figuresize)
#     make_division_II_plot(figuresize)
    make_permutation_plot(figuresize)