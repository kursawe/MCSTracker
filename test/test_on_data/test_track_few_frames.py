# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This tests our first tracking example
"""
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__),'..','..','src'))

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
 
class TestTrackFewFrames(unittest.TestCase):
                                 
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

    def test_maximum_common_subgraph_for_data_pair(self):
        """Load two segmented data frames and track them"""

        mesh_sequence = mesh.read_sequence_from_data(path.join(dirname(__file__),'data','first_few_frames')) 
        
        first_mesh = mesh_sequence[0]
        
        second_mesh = mesh_sequence[1]

        print('I am testing...')
        subgraph_finder = tracking.LocalisedSubgraphFinder(first_mesh, second_mesh)
        print('Can I find a subgraph?')
        subgraph_finder.find_maximum_common_subgraph()
        print(len(subgraph_finder.largest_mappings))
        post_processor = tracking.PostProcessor(first_mesh, second_mesh, subgraph_finder.largest_mappings)
        post_processor.index_global_ids_from_largest_mappings()

#         tracked_ids = tracking.find_maximum_common_subgraph( first_mesh, second_mesh )

        first_mesh.plot(path.join(dirname(__file__),'output','first_mesh_tracked.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( first_mesh.elements ))

        second_mesh.plot(path.join(dirname(__file__),'output','second_mesh_tracked.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(first_mesh.elements ))

        first_mesh.save(path.join(dirname(__file__),'output','first_mesh_tracked.mesh'))
        second_mesh.save(path.join(dirname(__file__),'output','second_mesh_tracked.mesh'))

        post_processor.tidy_current_mapping()

        first_mesh.plot(path.join(dirname(__file__),'output','first_mesh_tracked_cleaned.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())
        
        second_mesh.plot(path.join(dirname(__file__),'output','second_mesh_tracked_cleaned.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())

        first_mesh.save(path.join(dirname(__file__),'output','first_mesh_tracked_cleaned.mesh'))
        second_mesh.save(path.join(dirname(__file__),'output','second_mesh_tracked_cleaned.mesh'))

    def test_maximum_common_subgraph_for_second_data_pair(self):
        """Load two segmented data frames and track them"""
        
        mesh_sequence = mesh.read_sequence_from_data(path.join(dirname(__file__),'data','first_few_frames')) 
        
        first_mesh = mesh_sequence[1]
        
        second_mesh = mesh_sequence[2]
        
        first_mesh.plot(path.join(dirname(__file__),'output','second_mesh_not_tracked.pdf') )
        second_mesh.plot(path.join(dirname(__file__),'output','third_mesh_not_tracked.pdf') )

        subgraph_finder = tracking.LocalisedSubgraphFinder(first_mesh, second_mesh)
        subgraph_finder.find_maximum_common_subgraph()
        post_processor = tracking.PostProcessor(first_mesh, second_mesh, subgraph_finder.largest_mappings)
        post_processor.index_global_ids_from_largest_mappings()

#         tracked_ids = tracking.find_maximum_common_subgraph( first_mesh, second_mesh )

        first_mesh.plot(path.join(dirname(__file__),'output','second_mesh_tracked_with_third_mesh.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())
        
        second_mesh.plot(path.join(dirname(__file__),'output','third_mesh_tracked.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())

        first_mesh.save(path.join(dirname(__file__),'output','second_mesh_tracked_with_third_mesh.mesh'))
        second_mesh.save(path.join(dirname(__file__),'output','third_mesh_tracked.mesh'))
        post_processor.tidy_current_mapping()

        first_mesh.plot(path.join(dirname(__file__),'output','second_mesh_tracked_with_third_mesh_cleaned.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())
        
        second_mesh.plot(path.join(dirname(__file__),'output','third_mesh_tracked_cleaned.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())
        
        first_mesh.save(path.join(dirname(__file__),'output','second_mesh_tracked_with_third_mesh_cleaned.mesh'))
        second_mesh.save(path.join(dirname(__file__),'output','third_mesh_tracked_cleaned.mesh'))

    def test_post_processing_for_first_data_pair(self):
        """Load two segmented data frames and track them"""
        
        first_mesh = mesh.load(path.join(dirname(__file__),'output','first_mesh_tracked_cleaned.mesh'))
        second_mesh = mesh.load(path.join(dirname(__file__),'output','second_mesh_tracked_cleaned.mesh'))

        largest_mapping = {}
        
        for element in first_mesh.elements:
            if element.global_id != None:
                second_frame_id = second_mesh.get_element_with_global_id(element.global_id).id_in_frame
                largest_mapping[element.id_in_frame] = second_frame_id
        
        largest_mappings = [largest_mapping]

        post_processor = tracking.PostProcessor(first_mesh, second_mesh, largest_mappings)
        post_processor.index_global_ids_from_largest_mappings()

        post_processor.post_process_with_data()

        first_mesh.plot(path.join(dirname(__file__),'output','first_mesh_tracked_processed.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())
        
        second_mesh.plot(path.join(dirname(__file__),'output','second_mesh_tracked_processed.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())
        
        first_mesh.save(path.join(dirname(__file__),'output','first_mesh_tracked_processed.mesh'))
        second_mesh.save(path.join(dirname(__file__),'output','second_mesh_tracked_processed.mesh'))
        global_ids = []
        for element in second_mesh.elements:
            if element.global_id != None:
                if element.global_id in global_ids:
                    print('found double global id')
                else:
                    global_ids.append(element.global_id)

        global_ids_one = []
        for element in first_mesh.elements:
            if element.global_id != None:
                if element.global_id in global_ids_one:
                    print('found double global id')
                else:
                    global_ids_one.append(element.global_id)
                    
        self.assertEqual(len(global_ids_one), len(global_ids))

    def test_post_processing_for_second_data_pair(self):
        """Load two segmented data frames and track them"""
        
        first_mesh = mesh.load(path.join(dirname(__file__),'output','second_mesh_tracked_with_third_mesh_cleaned.mesh'))
        second_mesh = mesh.load(path.join(dirname(__file__),'output','third_mesh_tracked_cleaned.mesh'))

        largest_mapping = {}
        
        for element in first_mesh.elements:
            if element.global_id != None:
                second_frame_id = second_mesh.get_element_with_global_id(element.global_id).id_in_frame
                largest_mapping[element.id_in_frame] = second_frame_id
        
        largest_mappings = [largest_mapping]

        post_processor = tracking.PostProcessor(first_mesh, second_mesh, largest_mappings)
        post_processor.index_global_ids_from_largest_mappings()

        post_processor.post_process_with_data()

        first_mesh.plot(path.join(dirname(__file__),'output','second_mesh_tracked_with_third_mesh_processed.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())
        
        second_mesh.plot(path.join(dirname(__file__),'output','third_mesh_tracked_processed.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())
        
        first_mesh.save(path.join(dirname(__file__),'output','second_mesh_tracked_with_third_mesh_tracked_processed.mesh'))
        second_mesh.save(path.join(dirname(__file__),'output','third_mesh_tracked_processed.mesh'))

    def test_track_function_on_first_data(self):
        """See whether we can fully track the mesh using the track function of the tracking module"""

        mesh_sequence = mesh.read_sequence_from_data(path.join(dirname(__file__),'data','first_few_frames')) 
        
        first_mesh = mesh_sequence[0]
        second_mesh = mesh_sequence[1]

        start_time = time.clock()
        tracked_ids = tracking.track( first_mesh, second_mesh )
        end_time = time.clock()

        print('total time for first data pair is:')
        print(end_time - start_time)

        first_mesh.plot( path.join(dirname(__file__),'output','first_mesh_tracked_total.pdf'), color_by_global_id = True, 
                         total_number_of_global_ids = first_mesh.get_num_elements())
        
        second_mesh.plot(path.join(dirname(__file__),'output','second_mesh_tracked_total.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())
        
    def test_track_function_on_second_data(self):
        """Test tracking function on second data pair"""

        mesh_sequence = mesh.read_sequence_from_data(path.join(dirname(__file__),'data','first_few_frames')) 
        
        first_mesh = mesh_sequence[1]
        second_mesh = mesh_sequence[2]

        start_time = time.clock()
        tracked_ids = tracking.track( first_mesh, second_mesh )
        end_time = time.clock()

        print('total time for second data pair is:')
        print(end_time - start_time)

        first_mesh.plot(path.join(dirname(__file__),'output','second_mesh_tracked_with_third_total.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())
        
        second_mesh.plot(path.join(dirname(__file__),'output','third_mesh_tracked_total.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = first_mesh.get_num_elements())
 
        
    def test_track_and_write_sequence(self):

        mesh_sequence = mesh.read_sequence_from_data(path.join(dirname(__file__),'data','first_few_frames')) 
        
        tracking.track_and_write_sequence(path.join(dirname(__file__),'data','first_few_frames'),
                                          path.join(dirname(__file__),'output','first_few_frames'))
 
    def test_post_process_sequence(self):

        data_collector = tracking.analyse_tracked_sequence(path.join(dirname(__file__),
                                                                     'output','first_few_frames'))
                                                           
        for counter, this_mesh in enumerate(data_collector.mesh_sequence):
            print('mesh number')
            print(counter)
            print('no of cells')
            print(this_mesh.get_num_elements())
            print('average cell area')
            print(this_mesh.calculate_average_element_area())
            
        for mesh_data in data_collector.steps:
            print('step no')
            print(mesh_data.step_number)
            print('no of cells in first mesh')
            print(mesh_data.mesh_one.get_num_elements())
            print('no of tracked cells')
            print(mesh_data.number_of_tracked_cells)
            print('global ids of tracked cells')
            print(mesh_data.global_ids_of_tracked_cells)
            print('no of dying cells')
            print(mesh_data.number_dying_cells)
            print('global ids of dying cells')
            print(mesh_data.global_ids_of_dying_cells)
            print('average centroid displacement')
            print(mesh_data.average_centroid_displacement)
            print('maximal_centroid_displacement')
            print(mesh_data.maximal_centroid_displacement)
            print('minimal_centroid_displacement')
            print(mesh_data.minimal_centroid_displacement)
            print('print number of cells gaining edges')
            print(mesh_data.number_of_cells_gaining_edges)
            print('print number of cells losing edges')
            print(mesh_data.number_of_cells_loosing_edges)
        print('total_no_of_cell deaths')
        print(data_collector.number_dying_cells)
        print('global ids of dying cells')
        print(data_collector.global_ids_of_dying_cells)
        print('total average centroid displacement')
        print(data_collector.average_centroid_displacement)
        print('total maximal_centroid_displacement')
        print(data_collector.maximal_centroid_displacement)
        print('total minimal_centroid_displacement')
        print(data_collector.minimal_centroid_displacement)
        print('total print number of cells gaining edges')
        print(data_collector.number_of_cells_gaining_edges)
        print('total print number of cells losing edges')
        print(data_collector.number_of_cells_loosing_edges)
        print('no of shared tracked cells')
        print(data_collector.number_of_tracked_cells)
        print('global ids of shared tracked cells')
        print(data_collector.global_ids_of_tracked_cells)
        
