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
import copy
 
def make_initial_mapping_plot():

    first_mesh = mesh.load('first_mesh.mesh')
    second_mesh = mesh.load('second_mesh.mesh')
     
    subgraph_finder = tracking.LocalisedSubgraphFinder(first_mesh, second_mesh)
    initial_tracking_state = subgraph_finder.create_starting_tracking_state()

    cell_id_in_first_mesh = initial_tracking_state.id_map.keys()[0]
    cell_id_in_second_mesh = initial_tracking_state.id_map.values()[0]

    vertex = subgraph_finder.pick_next_vertex(initial_tracking_state, is_global = True )
    possible_images = subgraph_finder.get_mappable_vertices(initial_tracking_state, vertex)
    assert(len(possible_images) == 1)
    possible_image = possible_images[0]
    
    next_cell = subgraph_finder.network_one_index_lookup[vertex]
    next_image = subgraph_finder.network_two_index_lookup[possible_image]

    neighbours_in_first_mesh = nx.single_source_shortest_path_length(subgraph_finder.network_one, 
                                                                     next_cell, 
                                                                     cutoff=2 ).keys()

    first_mesh_polygon_list = []
    for element in first_mesh.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
         
        if element.id_in_frame in neighbours_in_first_mesh:
            this_polygon.set_facecolor('lightgrey')
            if element.id_in_frame == next_cell:
                this_polygon.set_facecolor('black')
        else: 
            this_polygon.set_facecolor([1.0, 1.0, 1.0])
        if element.id_in_frame == cell_id_in_first_mesh:
            this_polygon.set_facecolor('blue')
        first_mesh_polygon_list.append(this_polygon)

    first_mesh_polygon_collection = mpl.collections.PatchCollection(first_mesh_polygon_list, match_original = True)
    first_mesh_figure = plt.figure()
    first_mesh_figure.gca().add_collection(first_mesh_polygon_collection)
    first_mesh_figure.gca().set_aspect('equal')
    first_mesh_figure.gca().autoscale_view()
    plt.axis('off')
    first_mesh_figure.savefig('first_mesh_second_match.pdf', bbox_inches = 'tight')
    plt.close(first_mesh_figure)
 
    second_mesh_polygon_list = []
    for element in second_mesh.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
         
        if element.id_in_frame == next_image:
            this_polygon.set_facecolor('black')
        else: 
            this_polygon.set_facecolor([1.0, 1.0, 1.0])
        if element.id_in_frame == cell_id_in_second_mesh:
            this_polygon.set_facecolor('blue')
        second_mesh_polygon_list.append(this_polygon)
 
    second_mesh_polygon_collection = mpl.collections.PatchCollection(second_mesh_polygon_list, match_original = True)
    second_mesh_figure = plt.figure()
    second_mesh_figure.gca().add_collection(second_mesh_polygon_collection)
    second_mesh_figure.gca().set_aspect('equal')
    second_mesh_figure.gca().autoscale_view()
    plt.axis('off')
    second_mesh_figure.savefig('second_mesh_second_match.pdf', bbox_inches = 'tight')
    plt.close(second_mesh_figure)

    first_mesh_removed_polygon_list = []
    for element in first_mesh.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
         
        if element.id_in_frame in neighbours_in_first_mesh:
            this_polygon.set_facecolor('lightgrey')
            if element.id_in_frame == next_cell:
                this_polygon.set_facecolor('white')
        else: 
            this_polygon.set_facecolor([1.0, 1.0, 1.0])
        if element.id_in_frame == cell_id_in_first_mesh:
            this_polygon.set_facecolor('blue')
        first_mesh_removed_polygon_list.append(this_polygon)

    first_mesh_removed_polygon_collection = mpl.collections.PatchCollection(first_mesh_removed_polygon_list, match_original = True)
    first_mesh_removed_figure = plt.figure()
    first_mesh_removed_figure.gca().add_collection(first_mesh_removed_polygon_collection)
    first_mesh_removed_figure.gca().set_aspect('equal')
    first_mesh_removed_figure.gca().autoscale_view()
    plt.axis('off')
    first_mesh_removed_figure.savefig('first_mesh_second_match_removed.pdf', bbox_inches = 'tight')
    plt.close(first_mesh_removed_figure)
 
    second_mesh_removed_polygon_list = []
    for element in second_mesh.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
         
        if element.id_in_frame == next_image:
            this_polygon.set_facecolor('white')
        else: 
            this_polygon.set_facecolor([1.0, 1.0, 1.0])
        if element.id_in_frame == cell_id_in_second_mesh:
            this_polygon.set_facecolor('blue')
        second_mesh_removed_polygon_list.append(this_polygon)
 
    second_mesh_removed_polygon_collection = mpl.collections.PatchCollection(second_mesh_removed_polygon_list, match_original = True)
    second_mesh_removed_figure = plt.figure()
    second_mesh_removed_figure.gca().add_collection(second_mesh_removed_polygon_collection)
    second_mesh_removed_figure.gca().set_aspect('equal')
    second_mesh_removed_figure.gca().autoscale_view()
    plt.axis('off')
    second_mesh_removed_figure.savefig('second_mesh_removed_second_match.pdf', bbox_inches = 'tight')
    plt.close(second_mesh_removed_figure)

if __name__=="__main__":
    print 'I am printing, too'
    make_initial_mapping_plot()