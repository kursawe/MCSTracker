# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This file contains preliminary work on making an analysis similar to the T1 analysis,
except that here we randomly remove cells rather than randomly swap edges.
This work is not finished yet"""

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__),'..','..','src'))

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

def make_death_analysis():
    """What we want is: three tables
    
    table 1
    percentage of t1 swaps | percentage_of_correctly_tracked_cells | 

    table 2
    percentage of t1 swaps | percentage_of_incorrectly_tracked_cells | 
    
    table 3
    percentage of t1 swaps | percentage_of_correctly_tracked_tissues | 
    """
    
    repetitions_number = 2
    success_ratio_matrix = np.zeros((15,2), dtype = 'int')
    correctly_tracked_cells_matrix = np.zeros((15,repetitions_number + 1), dtype = 'int')
    incorrectly_tracked_cells_matrix = np.zeros((15,repetitions_number + 1), dtype = 'int')
    for percentage_counter in range(2):
        success_ratio, correctly_tracked_cells, incorrectly_tracked_cells = get_success_ratio( percentage_counter, 
                                                                             repetitions_number)
        success_ratio_matrix[percentage_counter,0] = percentage_counter
        success_ratio_matrix[percentage_counter,1] = success_ratio
        correctly_tracked_cells_matrix[percentage_counter,0] = percentage_counter
        correctly_tracked_cells_matrix[percentage_counter,1:] = correctly_tracked_cells
        incorrectly_tracked_cells_matrix[percentage_counter,0] = percentage_counter
        incorrectly_tracked_cells_matrix[percentage_counter,1:] = incorrectly_tracked_cells

    np.savetxt(path.join(dirname(__file__),
                         'output','death_success_ratio.csv'), 
               success_ratio_matrix)

    np.savetxt(path.join(dirname(__file__),
                         'output','death_correctly_tracked_cells.csv'), 
               correctly_tracked_cells_matrix)

    np.savetxt(path.join(dirname(__file__),
                         'output','death_incorrectly_tracked_cells.csv'), 
               incorrectly_tracked_cells_matrix)

def get_success_ratio(death_percentage, repetitions_number):
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
    
    arguments = [death_percentage]*repetitions_number

    results = [pool.apply_async(try_test_success, args=(x,)) 
               for x in arguments]

    print('started pool for death swap percentage')
    print(death_percentage)

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
   
def try_test_success(death_percentage, output = False):
    try:
        success, correctly_tracked_cells, incorrectly_tracked_cells = test_success(death_percentage, output)
        return success, correctly_tracked_cells, incorrectly_tracked_cells
    except:
        return 0, 0, 0

def test_success(death_percentage, output = False):
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
    
    print('making meshes')
    first_mesh = mesh.creation.generate_random_tesselation(20,20)
    second_mesh = copy.deepcopy(first_mesh)
    
    first_mesh.assign_frame_ids_in_order()
    second_mesh.assign_frame_ids_randomly()

    # build ground truth for testing the mapping
    print('making ground truth')
    ground_truth = {}
    for element_index, element in enumerate(first_mesh.elements):
        ground_truth[element.id_in_frame] = second_mesh.elements[element_index].id_in_frame

    total_number_elements = first_mesh.get_num_elements()
    cells_to_kill = int( 0.01 * death_percentage * total_number_elements)
    print('killing cells')
    kill_cells( second_mesh, cells_to_kill, ground_truth )
    
    if output:
        first_mesh.save(path.join(dirname(__file__),'output','multiple_death_before.mesh'))
        second_mesh.save(path.join(dirname(__file__),'output','multiple_death_after.mesh'))
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_death_before.pdf') )
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_death_after.pdf') )
    
#     try:
#         tracking.track(first_mesh, second_mesh)
#     except:
#         return 0, 0, 0

    tracking.track(first_mesh, second_mesh)


    if output:
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_death_before.pdf'), 
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_death_after.pdf'),
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
    
    print('evaluate tracking')
    tracking_success, number_tracked_cells, number_incorrectly_tracked_cells = tracking.evaluate_tracking(first_mesh, 
                                                                                                          second_mesh, 
                                                                                                          ground_truth)
    
    percentage_correctly_tracked_cells = float(number_tracked_cells)/second_mesh.get_num_elements()
    percentage_incorrectly_tracked_cells = float(number_incorrectly_tracked_cells)/second_mesh.get_num_elements()

    return int(tracking_success), percentage_correctly_tracked_cells, percentage_incorrectly_tracked_cells
    
def kill_cells(sample_mesh, number_cells_to_kill, ground_truth):

    inverse_ground_truth = { value : key for key,value in ground_truth.items() }

    while number_cells_to_kill > 0:
        cells_to_kill = random.sample(sample_mesh.elements, number_cells_to_kill)
        number_cells_to_kill = 0
        for counter, cell in enumerate( cells_to_kill ):
            adjacent_element_ids = cell.get_ids_of_adjacent_elements()
            cell_can_be_killed = True
            for frame_id in adjacent_element_ids:
                neighbour_number = sample_mesh.get_element_with_frame_id(frame_id).get_num_nodes()
                if neighbour_number < 4:
                    cell_can_be_killed = False
                    break
            
            if cell_can_be_killed:
                sample_mesh.kill_element_with_frame_id(cell.id_in_frame)
                del ground_truth[inverse_ground_truth[cell.id_in_frame]]
            else:
                number_cells_to_kill +=1 

if __name__ == "__main__":
    make_death_analysis()
