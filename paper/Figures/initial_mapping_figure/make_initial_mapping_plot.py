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

#     print 'I am printing'
#     mesh_sequence = mesh.read_sequence_from_data(path.join(dirname(__file__),'..','..','..',
#                                                            'test','test_on_data',
#                                                            'data','first_few_frames')) 
#     
#     first_mesh = mesh_sequence[0]
#     
#     second_mesh = mesh_sequence[1]

    sys.setrecursionlimit(40000)
    first_mesh = mesh.creation.generate_random_tesselation(9,9)
    second_mesh = copy.deepcopy(first_mesh)
     
    first_mesh.assign_frame_ids_in_order()
    second_mesh.assign_frame_ids_randomly()
 
    first_mesh.save(path.join(dirname(__file__),'first_mesh.mesh'))
    second_mesh.save(path.join(dirname(__file__),'second_mesh.mesh'))
    
    subgraph_finder = tracking.LocalisedSubgraphFinder(first_mesh, second_mesh)
    
    initial_tracking_state = subgraph_finder.create_starting_tracking_state()
    cell_id_in_first_mesh = initial_tracking_state.id_map.keys()[0]
    cell_id_in_second_mesh = initial_tracking_state.id_map.values()[0]

    neighbours_in_first_mesh = nx.single_source_shortest_path_length(subgraph_finder.network_one, 
                                                                     cell_id_in_first_mesh, 
                                                                     cutoff=2 ).keys()

    neighbours_in_second_mesh = nx.single_source_shortest_path_length(subgraph_finder.network_two, 
                                                                      cell_id_in_second_mesh, 
                                                                      cutoff=2 ).keys()

    first_mesh_polygon_list = []
    for element in first_mesh.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
        
        if element.id_in_frame in neighbours_in_first_mesh:
            this_polygon.set_facecolor('lightgrey')
            if element.id_in_frame == cell_id_in_first_mesh:
                this_polygon.set_facecolor('blue')
        else: 
            this_polygon.set_facecolor([1.0, 1.0, 1.0])
        first_mesh_polygon_list.append(this_polygon)

    first_mesh_polygon_collection = mpl.collections.PatchCollection(first_mesh_polygon_list, match_original = True)
    first_mesh_figure = plt.figure()
    first_mesh_figure.gca().add_collection(first_mesh_polygon_collection)
    first_mesh_figure.gca().set_aspect('equal')
    first_mesh_figure.gca().autoscale_view()
    plt.axis('off')
    first_mesh_figure.savefig('first_mesh.pdf', bbox_inches = 'tight')
    plt.close(first_mesh_figure)

    second_mesh_polygon_list = []
    for element in second_mesh.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
        
        if element.id_in_frame in neighbours_in_second_mesh:
            this_polygon.set_facecolor('lightgrey')
            if element.id_in_frame == cell_id_in_second_mesh:
                this_polygon.set_facecolor('blue')
        else: 
            this_polygon.set_facecolor([1.0, 1.0, 1.0])
        second_mesh_polygon_list.append(this_polygon)

    second_mesh_polygon_collection = mpl.collections.PatchCollection(second_mesh_polygon_list, match_original = True)
    second_mesh_figure = plt.figure()
    second_mesh_figure.gca().add_collection(second_mesh_polygon_collection)
    second_mesh_figure.gca().set_aspect('equal')
    second_mesh_figure.gca().autoscale_view()
    plt.axis('off')
    second_mesh_figure.savefig('second_mesh.pdf', bbox_inches = 'tight')
    plt.close(second_mesh_figure)

if __name__=="__main__":
    print 'I am printing, too'
    make_initial_mapping_plot()