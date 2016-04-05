# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

import mesh
import tracking
import copy
import sys
import numpy as np
import random
import numpy as np
from os import path
from os.path import dirname
import multiprocessing as mp
sys.setrecursionlimit(40000)

"""This is a simple modified copy of make_t1_analysis.py that helps us
to count how many cells are involved for a given percentage of T1 swaps."""

def make_t1_analysis():
    """What we want is: three tables
    
    table 1
    percentage of t1 swaps | percentage_of_correctly_tracked_cells | 

    table 2
    percentage of t1 swaps | percentage_of_incorrectly_tracked_cells | 
    
    table 3
    percentage of t1 swaps | percentage_of_correctly_tracked_tissues | 
    """
    
    repetitions_number = 10
    success_ratio_matrix = np.zeros((15,2), dtype = 'int')
    correctly_tracked_cells_matrix = np.zeros((15,repetitions_number + 1), dtype = 'int')
    incorrectly_tracked_cells_matrix = np.zeros((15,repetitions_number + 1), dtype = 'int')
    for percentage_counter in range(15):
        success_ratio, correctly_tracked_cells, incorrectly_tracked_cells = get_success_ratio( percentage_counter, 
                                                                             repetitions_number)
        success_ratio_matrix[percentage_counter,0] = percentage_counter
        success_ratio_matrix[percentage_counter,1] = success_ratio
        correctly_tracked_cells_matrix[percentage_counter,0] = percentage_counter
        correctly_tracked_cells_matrix[percentage_counter,1:] = correctly_tracked_cells
        incorrectly_tracked_cells_matrix[percentage_counter,0] = percentage_counter
        incorrectly_tracked_cells_matrix[percentage_counter,1:] = incorrectly_tracked_cells

#     np.savetxt(path.join(dirname(__file__),
#                          'output','success_ratio.csv'), 
#                success_ratio_matrix)
# 
#     np.savetxt(path.join(dirname(__file__),
#                          'output','correctly_tracked_cells.csv'), 
#                correctly_tracked_cells_matrix)
# 
#     np.savetxt(path.join(dirname(__file__),
#                          'output','incorrectly_tracked_cells.csv'), 
#                incorrectly_tracked_cells_matrix)

def get_success_ratio(t1_percentage, repetitions_number):
    """Get success ratio of performing t1_number of t1 swaps repeatedly on
    repetitions_number different meshes.
    
    Parameters
    ----------
    
    t1_number : int
        how many t1 swaps should be performed
    
    repetitions_number : int
        how many repetitions should be performed
        
    Returns
    -------
    
    success_ratio : float between 0 and 1
        ratio of success on the repetitions_number meshes
        
    average_number_tracked_cells : float
        average number of correctly tracked cells
    """
    
    pool = mp.Pool(processes=4)
    
    arguments = [t1_percentage]*repetitions_number

    results = [pool.apply_async(try_test_success, args=(x,)) 
               for x in arguments]

    print 'started pool for t1 swap percentage'
    print t1_percentage

    results_list = []
    for result in results:
        try:
            this_result = result.get(10*60)
            results_list.append(this_result)
        except mp.TimeoutError:
            results_list.append([0, 0, 0])

    output = np.array(results_list)
    
    success_ratio = sum(output[:,0])/len(output)
    
    correctly_tracked_cells = output[:,1]
    incorrectly_tracked_cells = output[:,2]
    
    return success_ratio, correctly_tracked_cells, incorrectly_tracked_cells
   
def try_test_success(t1_percentage, output = False):
    try:
        success, correctly_tracked_cells, incorrectly_tracked_cells = test_success(t1_percentage, output)
        return success, correctly_tracked_cells, incorrectly_tracked_cells
    except:
        return 0, 0, 0

def test_success(t1_percentage, output = False):
    """Perform tracking test on a mesh where t1_number of t1 swaps are randomly performed
    
    Parameters
    ----------
    
    t1_number : int
        How many T1 swaps should be performed
        
    output : bool
        whether the meshes should be written to a file. Defaults to False.
        
    Returns
    -------
    
    success_boolean : bool
       True if tracking successful
       
    average_number_tracked_cells : float
       average number of cells that were tracked correctly
    """
    
    print 'making meshes'
    first_mesh = mesh.creation.generate_random_tesselation(20,20)
    second_mesh = copy.deepcopy(first_mesh)
    
    first_mesh.assign_frame_ids_in_order()
    second_mesh.assign_frame_ids_randomly()

    # build ground truth for testing the mapping
    print 'making ground truth'
    ground_truth = {}
    for element_index, element in enumerate(first_mesh.elements):
        ground_truth[element.id_in_frame] = second_mesh.elements[element_index].id_in_frame

    total_number_edges = len(first_mesh.collect_elements_of_inner_edges())
    edges_to_swap = int( 0.01 * t1_percentage * total_number_edges)
    print 'applying t1 swaps'
    apply_t1_swaps_to_mesh( second_mesh, edges_to_swap )
    
    counter = 0
    for cell in ground_truth:
        this_cell = first_mesh.get_element_with_frame_id(cell)
        image_cell = second_mesh.get_element_with_frame_id(ground_truth[cell])
        if np.linalg.norm( this_cell.calculate_centroid() -
                           image_cell.calculate_centroid()) > 0.0001:
            counter += 1
            
    print 'tissue with t1 percentage' 
    print t1_percentage
    print 'has this rearranged_cell percentage'
    print float(counter)/first_mesh.get_num_elements()
        
#     if output:
#         first_mesh.save(path.join(dirname(__file__),'output','multiple_t1_before.mesh'))
#         second_mesh.save(path.join(dirname(__file__),'output','multiple_t1_after.mesh'))
#     
#     try:
#         tracking.track(first_mesh, second_mesh)
#     except:
#         return 0, 0, 0

    if output:
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_before.pdf'), 
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_after.pdf'),
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
    
    print 'evaluate tracking'
    tracking_success, number_tracked_cells, number_incorrectly_tracked_cells = tracking.evaluate_tracking(first_mesh, 
                                                                                                          second_mesh, 
                                                                                                          ground_truth)
    
    percentage_correctly_tracked_cells = float(number_tracked_cells)/first_mesh.get_num_elements()
    percentage_incorrectly_tracked_cells = float(number_incorrectly_tracked_cells)/first_mesh.get_num_elements()

    return int(tracking_success), percentage_correctly_tracked_cells, percentage_incorrectly_tracked_cells
    
def apply_t1_swaps_to_mesh(sample_mesh, t1_number):
    """Apply t1_number T1 swaps on sample mesh.
    
    Parameters
    ----------
    
    sample_mesh : mesh instance
    
    t1_number : int
        number of t1 swaps that should be performed.
    """
    
    unswapped_edges = t1_number 
    
    while unswapped_edges > 0:
        all_inner_edges_by_elements = sample_mesh.collect_elements_of_inner_edges()
        if len(all_inner_edges_by_elements) < unswapped_edges:
            unswapped_edges = len(all_inner_edges_by_elements)
        edges_to_swap = random.sample(all_inner_edges_by_elements, unswapped_edges)
        unswapped_edges = 0
        for counter, edge in enumerate( edges_to_swap ):
            polygon_numbers_at_edge = [sample_mesh.get_element_with_frame_id(frame_id).get_num_nodes() 
                                       for frame_id in edge]
            if polygon_numbers_at_edge[0] > 3 and polygon_numbers_at_edge[1] > 3:
                nodes_at_edge = list(sample_mesh.get_nodes_shared_by_elements( list(edge) ))
                sample_mesh.perform_t1_swap(nodes_at_edge[0], nodes_at_edge[1])
            else:
                unswapped_edges += 1
        
if __name__ == "__main__":
    make_t1_analysis()
