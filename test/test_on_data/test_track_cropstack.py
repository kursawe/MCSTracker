# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This tests our first tracking example
"""
import unittest
import mesh
import tracking
import time
from os import path
from os.path import dirname
from nose.plugins.attrib import attr
import sys
import threading
import resource
 
class TestTrackCropstack(unittest.TestCase):
                                 
#    @attr(level = 'standard')
#    def xest_track_simple_case(self):
#        """testing something"""
#
#        threading.stack_size(67108864) # 64MB stack
#        sys.setrecursionlimit(2 ** 20)
#        # only new threads get the redefined stack size
#        thread = threading.Thread(target=self.maximum_common_subgraph_for_data_pair)
#        thread.start()
#        thread.join()

        
    def test_track_and_write_full_sequence(self):

#         mesh_sequence = mesh.read_sequence_from_data(path.join(dirname(__file__),
#                                                                '..','test_mesh','output','converted')) 
        
        tracking.track_and_write_sequence( path.join(dirname(__file__), 
                                                     '..','test_mesh','output','converted'),
                                           path.join(dirname(__file__),
                                                     'output','crostack_track'),
                                           start_number = 3, number_meshes = 9 ) 
        
    def test_visualize_sequence(self):
        tracking.plot_tracked_sequence( path.join(dirname(__file__),
                                                  'output','crostack_track'),
                                        path.join(dirname(__file__),
                                                  '..', 'test_mesh', 'data', 'image_data'),
                                        path.join(dirname(__file__),
                                                  '..', 'test_mesh', 'output', 'converted' ),
                                        path.join(dirname(__file__),
                                                  'output','cropstack_track_visualized') )
 
    def xest_post_process_sequence(self):

        data_collector = tracking.analyse_tracked_sequence(path.join(dirname(__file__),
                                                                     'output','first_few_frames'))
                                                           
        for counter, this_mesh in enumerate(data_collector.mesh_sequence):
            print 'mesh number'
            print counter
            print 'no of cells'
            print this_mesh.get_num_elements()
            print 'average cell area'
            print this_mesh.calculate_average_element_area()
            
        for mesh_data in data_collector.steps:
            print 'step no'
            print mesh_data.step_number
            print 'no of cells in first mesh'
            print mesh_data.mesh_one.get_num_elements()
            print 'no of tracked cells'
            print mesh_data.number_of_tracked_cells
            print 'global ids of tracked cells'
            print mesh_data.global_ids_of_tracked_cells
            print 'no of dying cells'
            print mesh_data.number_dying_cells
            print 'global ids of dying cells'
            print mesh_data.global_ids_of_dying_cells
            print 'average centroid displacement'
            print mesh_data.average_centroid_displacement
            print 'maximal_centroid_displacement'
            print mesh_data.maximal_centroid_displacement
            print 'minimal_centroid_displacement'
            print mesh_data.minimal_centroid_displacement
            print 'print number of cells gaining edges'
            print mesh_data.number_of_cells_gaining_edges
            print 'print number of cells losing edges'
            print mesh_data.number_of_cells_loosing_edges
        print 'total_no_of_cell deaths'
        print data_collector.number_dying_cells
        print 'global ids of dying cells'
        print data_collector.global_ids_of_dying_cells
        print 'total average centroid displacement'
        print data_collector.average_centroid_displacement
        print 'total maximal_centroid_displacement'
        print data_collector.maximal_centroid_displacement
        print 'total minimal_centroid_displacement'
        print data_collector.minimal_centroid_displacement
        print 'total print number of cells gaining edges'
        print data_collector.number_of_cells_gaining_edges
        print 'total print number of cells losing edges'
        print data_collector.number_of_cells_loosing_edges
        print 'no of shared tracked cells'
        print data_collector.number_of_tracked_cells
        print 'global ids of shared tracked cells'
        print data_collector.global_ids_of_tracked_cells
        
