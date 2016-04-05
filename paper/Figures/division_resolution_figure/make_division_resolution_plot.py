# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This tests our first tracking example
"""
import mesh
import tracking
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from os import path
from os.path import dirname
import sys
import pickle
import copy
 
def make_simple_division_plot():
    """We load a mesh where only one maximum common subgraph should be identified, and track the division event"""
 
    mesh_one = mesh.load(path.join(dirname(__file__),'..','..','..',
                                   'test','test_tracking','data',
                                   'standard_simple_division_one.mesh'))
    mesh_two = mesh.load(path.join(dirname(__file__),'..','..','..',
                                   'test','test_tracking','data',
                                   'standard_simple_division_two.mesh'))
 
    # First pick a cell in the centre
    mesh_centre = mesh_two.calculate_centre()
     
    mesh_one.assign_frame_ids_in_order()
    mesh_two.assign_frame_ids_in_order()
 
    most_central_element = mesh_two.find_most_central_element()
    
    most_central_element_id = most_central_element.id_in_frame
            
    mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame,
                                                       [1.0,1.2])
     
    daughter_cell_id_one = mesh_one.get_maximal_frame_id() + 1
    daughter_cell_id_two = mesh_one.get_maximal_frame_id() + 2

    cells_adjacent_to_daughter_cell_one = mesh_two.get_element_with_frame_id(daughter_cell_id_one).get_ids_of_adjacent_elements()
    cells_adjacent_to_daughter_cell_two = mesh_two.get_element_with_frame_id(daughter_cell_id_two).get_ids_of_adjacent_elements()
    
    bordering_cells = set.intersection(set(cells_adjacent_to_daughter_cell_one),
                                       set(cells_adjacent_to_daughter_cell_two))
    
    first_mesh_polygon_list = []
    for element in mesh_one.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
         
        if element.id_in_frame in bordering_cells:
            this_polygon.set_facecolor('slategrey')
        elif element.id_in_frame == most_central_element_id:
                this_polygon.set_facecolor('lightskyblue')
        else: 
            this_polygon.set_facecolor([1.0, 1.0, 1.0])
        first_mesh_polygon_list.append(this_polygon)

    first_mesh_polygon_collection = mpl.collections.PatchCollection(first_mesh_polygon_list, match_original = True)
    first_mesh_figure = plt.figure()
    first_mesh_figure.gca().add_collection(first_mesh_polygon_collection)
    first_mesh_figure.gca().set_aspect('equal')
    first_mesh_figure.gca().autoscale_view()
    plt.axis('off')
    first_mesh_figure.savefig('simple_division_first.pdf', bbox_inches = 'tight')
    plt.close(first_mesh_figure)
 
    second_mesh_polygon_list = []
    for element in mesh_two.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
         
        if element.id_in_frame in bordering_cells:
            this_polygon.set_facecolor('slategrey')
        elif ( element.id_in_frame == daughter_cell_id_one or
               element.id_in_frame == daughter_cell_id_two): 
            this_polygon.set_facecolor('lightskyblue')
        else:
            this_polygon.set_facecolor('white')
        second_mesh_polygon_list.append(this_polygon)
 
    second_mesh_polygon_collection = mpl.collections.PatchCollection(second_mesh_polygon_list, match_original = True)
    second_mesh_figure = plt.figure()
    second_mesh_figure.gca().add_collection(second_mesh_polygon_collection)
    second_mesh_figure.gca().set_aspect('equal')
    second_mesh_figure.gca().autoscale_view()
    plt.axis('off')
    second_mesh_figure.savefig('simple_division_second.pdf', bbox_inches = 'tight')
    plt.close(second_mesh_figure)

def make_four_sided_division_plot():
    """We load a mesh where only one maximum common subgraph should be identified, and track the division event"""
 
    mesh_one = mesh.load(path.join(dirname(__file__),'..','..','..',
                                   'test','test_tracking','data',
                                   'standard_ambiguous_division_one.mesh'))
    mesh_two = mesh.load(path.join(dirname(__file__),'..','..','..',
                                   'test','test_tracking','data',
                                   'standard_ambiguous_division_two.mesh'))
 
    # First pick a cell in the centre
    mesh_centre = mesh_two.calculate_centre()
     
    mesh_one.assign_frame_ids_in_order()
    mesh_two.assign_frame_ids_in_order()
 
    most_central_element = mesh_two.find_most_central_element()
    
    most_central_element_id = most_central_element.id_in_frame
            
    mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame,
                                                       [1.0,1.0])
     
    daughter_cell_id_one = mesh_one.get_maximal_frame_id() + 1
    daughter_cell_id_two = mesh_one.get_maximal_frame_id() + 2

    cells_adjacent_to_daughter_cell_one = mesh_two.get_element_with_frame_id(daughter_cell_id_one).get_ids_of_adjacent_elements()
    cells_adjacent_to_daughter_cell_two = mesh_two.get_element_with_frame_id(daughter_cell_id_two).get_ids_of_adjacent_elements()
    
    bordering_cells = set.intersection(set(cells_adjacent_to_daughter_cell_one),
                                       set(cells_adjacent_to_daughter_cell_two))
    
    first_mesh_polygon_list = []
    for element in mesh_one.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
         
        if element.id_in_frame in bordering_cells:
            this_polygon.set_facecolor('slategrey')
        elif element.id_in_frame == most_central_element_id:
                this_polygon.set_facecolor('lightskyblue')
        else: 
            this_polygon.set_facecolor([1.0, 1.0, 1.0])
        first_mesh_polygon_list.append(this_polygon)

    first_mesh_polygon_collection = mpl.collections.PatchCollection(first_mesh_polygon_list, match_original = True)
    first_mesh_figure = plt.figure()
    first_mesh_figure.gca().add_collection(first_mesh_polygon_collection)
    first_mesh_figure.gca().set_aspect('equal')
    first_mesh_figure.gca().autoscale_view()
    plt.axis('off')
    first_mesh_figure.savefig('four_sided_division_first.pdf', bbox_inches = 'tight')
    plt.close(first_mesh_figure)
 
    second_mesh_polygon_list = []
    for element in mesh_two.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
         
        if element.id_in_frame in bordering_cells:
            this_polygon.set_facecolor('slategrey')
        elif ( element.id_in_frame == daughter_cell_id_one or
               element.id_in_frame == daughter_cell_id_two): 
            this_polygon.set_facecolor('lightskyblue')
        else:
            this_polygon.set_facecolor('white')
        second_mesh_polygon_list.append(this_polygon)
 
    second_mesh_polygon_collection = mpl.collections.PatchCollection(second_mesh_polygon_list, match_original = True)
    second_mesh_figure = plt.figure()
    second_mesh_figure.gca().add_collection(second_mesh_polygon_collection)
    second_mesh_figure.gca().set_aspect('equal')
    second_mesh_figure.gca().autoscale_view()
    plt.axis('off')
    second_mesh_figure.savefig('four_sided_division_second.pdf', bbox_inches = 'tight')
    plt.close(second_mesh_figure)

def make_three_sided_division_plot():
    """We load a mesh where only one maximum common subgraph should be identified, and track the division event"""
 
    mesh_one = mesh.load(path.join(dirname(__file__),'..','..','..',
                                   'test','test_tracking','data',
                                   'division_special_mesh_six.mesh'))
    mesh_two = mesh.load(path.join(dirname(__file__),'..','..','..',
                                   'test','test_tracking','data',
                                   'division_special_mesh_six.mesh'))
 
    # First pick a cell in the centre
    mesh_centre = mesh_two.calculate_centre()
     
    mesh_one.assign_frame_ids_in_order()
    mesh_two.assign_frame_ids_in_order()
 
    most_central_element = mesh_two.find_most_central_element()
    
    most_central_element_id = most_central_element.id_in_frame

    file_to_read = open(path.join(dirname(__file__),'..','..','..',
                                  'test','test_tracking','data',
                                  'division_special_direction_six.pickle'), 'r')
    division_direction = pickle.load(file_to_read)
    file_to_read.close()

    mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame,
                                                       division_direction)
     
    daughter_cell_id_one = mesh_one.get_maximal_frame_id() + 1
    daughter_cell_id_two = mesh_one.get_maximal_frame_id() + 2

    cells_adjacent_to_daughter_cell_one = mesh_two.get_element_with_frame_id(daughter_cell_id_one).get_ids_of_adjacent_elements()
    cells_adjacent_to_daughter_cell_two = mesh_two.get_element_with_frame_id(daughter_cell_id_two).get_ids_of_adjacent_elements()
    
    bordering_cells = set.intersection(set(cells_adjacent_to_daughter_cell_one),
                                       set(cells_adjacent_to_daughter_cell_two))
    
    first_mesh_polygon_list = []
    for element in mesh_one.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
         
        if element.id_in_frame in bordering_cells:
            this_polygon.set_facecolor('slategrey')
        elif element.id_in_frame == most_central_element_id:
                this_polygon.set_facecolor('lightskyblue')
        else: 
            this_polygon.set_facecolor([1.0, 1.0, 1.0])
        first_mesh_polygon_list.append(this_polygon)

    first_mesh_polygon_collection = mpl.collections.PatchCollection(first_mesh_polygon_list, match_original = True)
    first_mesh_figure = plt.figure()
    first_mesh_figure.gca().add_collection(first_mesh_polygon_collection)
    first_mesh_figure.gca().set_aspect('equal')
    first_mesh_figure.gca().autoscale_view()
    plt.axis('off')
    first_mesh_figure.savefig('three_sided_division_first.pdf', bbox_inches = 'tight')
    plt.close(first_mesh_figure)
 
    second_mesh_polygon_list = []
    for element in mesh_two.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
         
        if element.id_in_frame in bordering_cells:
            this_polygon.set_facecolor('slategrey')
        elif ( element.id_in_frame == daughter_cell_id_one or
               element.id_in_frame == daughter_cell_id_two): 
            this_polygon.set_facecolor('lightskyblue')
        else:
            this_polygon.set_facecolor('white')
        second_mesh_polygon_list.append(this_polygon)
 
    second_mesh_polygon_collection = mpl.collections.PatchCollection(second_mesh_polygon_list, match_original = True)
    second_mesh_figure = plt.figure()
    second_mesh_figure.gca().add_collection(second_mesh_polygon_collection)
    second_mesh_figure.gca().set_aspect('equal')
    second_mesh_figure.gca().autoscale_view()
    plt.axis('off')
    second_mesh_figure.savefig('three_sided_division_second.pdf', bbox_inches = 'tight')
    plt.close(second_mesh_figure)

if __name__=="__main__":
    make_simple_division_plot()
    make_four_sided_division_plot()
    make_three_sided_division_plot()
