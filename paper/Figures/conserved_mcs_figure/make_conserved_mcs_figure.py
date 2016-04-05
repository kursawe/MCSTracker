import mesh
import tracking
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
from os import path
from os.path import dirname

def make_all_plots():
    first_mesh = mesh.load(path.join(dirname(__file__),
                                   '..','..','..',
                                   'test', 'test_tracking',
                                   'output','multiple_mcs_before.mesh'))
    second_mesh = mesh.load(path.join(dirname(__file__),
                                   '..','..','..',
                                   'test', 'test_tracking',
                                   'output','multiple_mcs_before.mesh'))

    mesh_centre = second_mesh.calculate_centre()
    # pick the node closest to the centre
    min_distance = 3*second_mesh.calculate_height()
    for node in second_mesh.nodes:
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
    
    second_mesh.perform_t1_swap( most_central_node.id, one_edge_node.id )
    
    first_mesh.plot('first_mesh.pdf')
    second_mesh.plot('second_mesh.pdf')

    subgraph_finder = tracking.KrissinelMaximumCommonSubgraphFinder( first_mesh, 
                                                                     second_mesh,
                                                                     0.4 )
    subgraph_finder.find_maximum_common_subgraph()
    
    for counter, largest_mapping in enumerate(subgraph_finder.largest_mappings):
        tracked_ids = []
        for element_one in first_mesh.elements:
            element_one.global_id = None

        for element_two in second_mesh.elements:
            element_two.global_id = None

        for global_id, frame_one_id in enumerate(largest_mapping):
            first_mesh.get_element_with_frame_id(frame_one_id).global_id = global_id
            second_mesh.get_element_with_frame_id(largest_mapping[frame_one_id]).global_id = global_id
            tracked_ids.append(global_id)

        first_mesh.index_global_ids()
        second_mesh.index_global_ids()

        first_mesh.plot('first_mesh_tracked_' + str(counter) + '.pdf', color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        second_mesh.plot('second_mesh_tracked_' + str(counter) + '.pdf', color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        
        first_mesh_polygon_list = []
        for element in first_mesh.elements:
            this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                                fill = True)
            if element.id_in_frame in largest_mapping:
                this_polygon.set_facecolor('slategrey')
            else: 
                this_polygon.set_facecolor('white')
            first_mesh_polygon_list.append(this_polygon)
        first_mesh_polygon_collection = mpl.collections.PatchCollection(first_mesh_polygon_list,
                                                                        match_original = True)

        first_figure = plt.figure()
        first_figure.gca().add_collection(first_mesh_polygon_collection)
        first_figure.gca().set_aspect('equal')
        first_figure.gca().autoscale_view()
        plt.axis('off')
        first_figure.savefig('maping_nr_' + str(counter) +
                                  '_first_mesh.pdf', bbox_inches = 'tight')
        plt.close(first_figure)
 
        second_mesh_polygon_list = []
        for element in second_mesh.elements:
            this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                                fill = True)
            if element.id_in_frame in largest_mapping.values():
                this_polygon.set_facecolor('slategrey')
            else: 
                this_polygon.set_facecolor('white')
            second_mesh_polygon_list.append(this_polygon)
        second_mesh_polygon_collection = mpl.collections.PatchCollection( second_mesh_polygon_list,
                                                                          match_original = True)

        second_figure = plt.figure()
        second_figure.gca().add_collection(second_mesh_polygon_collection)
        second_figure.gca().set_aspect('equal')
        second_figure.gca().autoscale_view()
        plt.axis('off')
        second_figure.savefig('maping_nr_' + str(counter) +
                                  '_second_mesh.pdf', bbox_inches = 'tight')
        plt.close(second_figure)

if __name__ == "__main__":

    make_all_plots()