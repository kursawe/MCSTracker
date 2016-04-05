# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This tests our first tracking example
"""
import unittest
import mesh
import tracking
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import os
import copy
import sys
from os import path
from os.path import dirname
from nose.plugins.attrib import attr

class TestTracking(unittest.TestCase):
                                 
    @attr(level = 'standard')
    def test_visualise_search(self):
        """generate a small random mesh and track it.
        """
        sys.setrecursionlimit(40000)
        mesh_one = mesh.creation.generate_random_tesselation(5,5)
        
        mesh_two = copy.deepcopy(mesh_one)
         
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
         
        tracked_ids = tracking.track(mesh_one, mesh_two)
 
        mesh_two.plot(path.join(dirname(__file__),'output','visualisation_mesh_after_tracking.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = mesh_one.get_num_elements())
         
    @attr(level = 'standard')
    def test_plot_collection_into_separate_figure(self):
        """plot a polygon collection into a figure designed by us 
        """
        sys.setrecursionlimit(40000)
        mesh_one = mesh.creation.generate_random_tesselation(5,5)
       
        print 'start to copy'
        mesh_two = copy.deepcopy(mesh_one)
        print 'stop to copy'
         
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        tracked_ids = tracking.track(mesh_one, mesh_two)

        # figure properties
        figuresize = (4,2.75)
        font = {'size'   : 10}
        plt.rc('font', **font)

        mesh_figure = plt.figure(figsize = figuresize)

        ax = plt.axes([0,0,1,1])
        polygon_collection = mesh_two.get_polygon_collection(color_by_global_id = True,
                                           total_number_of_global_ids = mesh_one.get_num_elements())
        mesh_figure.gca().add_collection(polygon_collection)
        mesh_figure.gca().set_aspect('equal')
        mesh_figure.gca().autoscale_view()
        plt.axis('off')
        filename = path.join(dirname(__file__),'output','own_visualisation_figure.pdf')
        mesh_figure.savefig(filename, bbox_inches = 'tight')
        
    def test_plot_all_global_ids_same_color(self):
        """plot all global ids same color
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.creation.generate_random_tesselation(8,8)
        
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # Perform T1 swap on mesh two 
        # First pick a node in the centre
        mesh_centre = mesh_two.calculate_centre()

        # pick the node closest to the centre
        most_central_node = mesh_two.find_most_central_node()
               
        # pick a node that shares an edge with this central node
        for local_index, element_node in enumerate(most_central_node.adjacent_elements[0].nodes):
            if element_node.id == most_central_node.id:
                num_nodes_this_element = most_central_node.adjacent_elements[0].get_num_nodes()
                one_edge_node = most_central_node.adjacent_elements[0].nodes[(local_index+1)%num_nodes_this_element]
                break
        
        mesh_two.perform_t1_swap( most_central_node.id, one_edge_node.id )
         
        tracked_ids = tracking.find_maximum_common_subgraph(mesh_one, mesh_two)

        mesh_two.plot(path.join(dirname(__file__),'output','visualisation_in_green.pdf'), color_by_global_id = 'g', 
                      total_number_of_global_ids = mesh_one.get_num_elements())