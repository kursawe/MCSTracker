# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__),'..','..','src'))

import mesh
import tracking
import copy
import numpy as np
import random
import numpy as np
import sys
from os import path
from os.path import dirname
sys.setrecursionlimit(40000)

def make_multiple_mcs():
    
    counter = 0
    example_not_found = True
    while example_not_found:
        print('running example number')
        print(counter)
        example_not_found = run_single_repeat()
        counter += 1
       
def run_single_repeat():
    first_mesh = mesh.creation.generate_random_tesselation(6,6)
    second_mesh = copy.deepcopy(first_mesh)
    
    first_mesh.assign_frame_ids_in_order()
    second_mesh.assign_frame_ids_randomly()

    first_mesh.save(path.join(dirname(__file__),'output','multiple_mcs_before.mesh'))
    second_mesh.save(path.join(dirname(__file__),'output','multiple_mcs_after.mesh'))
    
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

    subgraph_finder = tracking.KrissinelMaximumCommonSubgraphFinder( first_mesh, 
                                                                     second_mesh,
                                                                     0.4 )

    subgraph_finder.find_maximum_common_subgraph()

    if len(subgraph_finder.largest_mappings)<=1: 
        return True
    else:
        return False
    
if __name__ == "__main__":
    make_multiple_mcs()
