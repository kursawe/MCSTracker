# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

import unittest
import mesh
import tracking
import copy
import sys
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import pickle
from os import path
from os.path import dirname
import multiprocessing as mp

def make_timing_array():
    """Obtain averages for the times it takes to run mcs finding
    This test needs to be run BEFORE plot_timings"""
    
    # make a dictionary with mesh sizes as keys and computation times as Values
    sizes_and_times = {}
    
    number_of_runs_per_size = 10
    # sizes = [5,10,20,30,40,70,100]
    sizes = [5,10,20,30,40]

    argument_tuples= []
    for size in sizes:
        for index in range(number_of_runs_per_size):
            argument_tuples.append(size)

    pool = mp.Pool(processes=3)
    results = [pool.apply_async(get_single_mcs_time, args=(x,)) 
               for x in argument_tuples]
    output = [result.get() for result in results]
    
    size_matrix = []
    counter = 0
    for size in sizes:
        this_line = [size]
        for index in range(number_of_runs_per_size):
            this_line.append(output[counter])
            counter += 1
        size_matrix.append(this_line)
        
    size_matrix_np = np.array(size_matrix)

    np.savetxt(path.join(dirname(__file__),
                         'output','parallel_timings_array.csv'), 
               size_matrix_np)
    
#     plot_figure = plt.figure(figsize = figuresize)
#     plt.xlabel('Mesh size')
#     plt.ylabel('Calculation time')
#     plt.plot(sizes_and_times[:,0], sizes_and_times[:,1])
#     plot_figure.tight_layout()
#     file_name = "Calculation_times.pdf"
#     plot_figure.savefig(path.join(dirname(__file__),'data',file_name))

def get_mcs_times(size):
    
    calculation_times = []
    
    for instance in range(10):
        calculation_times.append( self.get_single_mcs_time(size) )
        
    calculation_times_np = np.array(calculation_times)
    
    return calculation_times

def get_single_mcs_time(size):

    sys.setrecursionlimit(40000)

    mesh_one = mesh.creation.generate_random_tesselation(size,size)
    mesh_two = copy.deepcopy(mesh_one)
     
    mesh_one.assign_frame_ids_in_order()
    mesh_two.assign_frame_ids_randomly()

    subgraph_finder = tracking.LocalisedSubgraphFinder(mesh_one, mesh_two)
    subgraph_finder.turn_timing_on()
    subgraph_finder.find_maximum_common_subgraph()
    print 'finished subraph of size'
    print size
    
    return subgraph_finder.total_execution_time

if __name__ == "__main__":
    make_timing_array()
