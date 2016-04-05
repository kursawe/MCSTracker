# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/LICENSE.

"""This tests our first tracking example
"""
import unittest
import mesh
import tracking
import copy
import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle
from os import path
from os.path import dirname
from nose.plugins.attrib import attr

class TestTrackDivision(unittest.TestCase):
                                 
    @attr(level = 'standard')
    def test_maximum_common_subgraph_for_division(self):
        """generate a random mesh, copy it, perform a division event, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.creation.generate_random_tesselation(7,7)
         
        mesh_two = copy.deepcopy(mesh_one)
#         mesh_one = mesh.load('standard_ambiguous_division_one.mesh')
#         mesh_two = mesh.load('standard_ambiguous_division_two.mesh')

        # Perform T1 swap on mesh two 
        # First pick a cell in the centre
        mesh_centre = mesh_two.calculate_centre()
        
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # pick the_central_element = element 
        most_central_element = mesh_two.find_most_central_element()
               
        mesh_two.divide_element_with_frame_id(most_central_element.id_in_frame)
         
        tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
 
        mesh_one.plot(path.join(dirname(__file__),'output','mesh_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','mesh_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
         
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()
 
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
     
        plt.close('all')

    @attr(level = 'standard')
    def test_track_division_standard_case(self):
        """We load a mesh where only one maximum common subgraph should be identified, and track the division event"""
 
        mesh_one = mesh.load(path.join(dirname(__file__),'data','standard_simple_division_one.mesh'))
        mesh_two = mesh.load(path.join(dirname(__file__),'data','standard_simple_division_two.mesh'))
 
        # First pick a cell in the centre
        mesh_centre = mesh_two.calculate_centre()
         
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
 
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame
 
        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
                
        mesh_two.divide_element_with_frame_id(most_central_element.id_in_frame)
         
        tracked_ids = tracking.track( mesh_one, mesh_two )
  
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_standard_mesh_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_standard_mesh_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
          
        # make sure that the entire mesh was tracked (except for the dead cell)
        self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() - 2 )
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )
             
        plt.close('all')

    @attr(level = 'standard')
    def test_track_division_standard_ambiguous_case(self):
        """We load a mesh where two maximum common subgraphs will be identified, and track the division event"""
 
        mesh_one = mesh.load(path.join(dirname(__file__),'data','standard_ambiguous_division_one.mesh'))
        mesh_two = mesh.load(path.join(dirname(__file__),'data','standard_ambiguous_division_two.mesh'))
 
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame
 
        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
                
        mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, [1.0, 1.0])
         
        tracked_ids = tracking.track( mesh_one, mesh_two )
  
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_ambiguous_mesh_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_ambiguous_mesh_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
          
        # make sure that the entire mesh was tracked (except for the dead cell)
        self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() - 2 )
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )
             
        plt.close('all')       

    @attr(level = 'standard')
    def test_track_special_division(self):
        """read a special mesh, perform a division event, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','division_special_mesh_one.mesh'))
        mesh_two = mesh.load(path.join(dirname(__file__),'data','division_special_mesh_two.mesh'))

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        file_to_read = open(path.join(dirname(__file__),'data','last_division_direction.pickle'), 'r')
        division_direction = pickle.load(file_to_read)
        file_to_read.close()

        mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, division_direction)
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
#         largest_mappings = subgraph_finder.largest_mappings
# 
#         for mapping_index, large_mapping in enumerate(largest_mappings):
# 
#             tracked_ids = []
#             
#             for element_one in mesh_one.elements:
#                 element_one.global_id = None
# 
#             for element_two in mesh_two.elements:
#                 element_two.global_id = None
# 
#             for global_id, frame_one_id in enumerate(large_mapping):
#                 mesh_one.get_element_with_frame_id(frame_one_id).global_id = global_id
#                 mesh_two.get_element_with_frame_id(large_mapping[frame_one_id]).global_id = global_id
#                 tracked_ids.append(global_id)
# 
#             mesh_one.index_global_ids()
#             mesh_two.index_global_ids()
# 
#             mesh_one.plot('mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
         
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_special_mesh_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_special_mesh_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
         
        # make sure that the entire mesh was tracked (except for the dead cell)
        self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() - 1 )

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            if (ground_truth[element_one.id_in_frame] != element_two.id_in_frame ):
                print element_one.calculate_centroid()
                print element_two.calculate_centroid()
            
        plt.close('all')
        
    @attr(level = 'standard')
    def test_track_special_division_one(self):
        """read a special mesh, perform a division event, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','division_special_one_mesh_one.mesh'))
        mesh_two= mesh.load(path.join(dirname(__file__),'data','division_special_one_mesh_two.mesh'))

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        file_to_read = open(path.join(dirname(__file__),'data','division_special_direction_one.pickle'), 'r')
        division_direction = pickle.load(file_to_read)
        file_to_read.close()

        mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, division_direction)
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
#         largest_mappings = subgraph_finder.largest_mappings
# 
#         for mapping_index, large_mapping in enumerate(largest_mappings):
# 
#             tracked_ids = []
#             
#             for element_one in mesh_one.elements:
#                 element_one.global_id = None
# 
#             for element_two in mesh_two.elements:
#                 element_two.global_id = None
# 
#             for global_id, frame_one_id in enumerate(large_mapping):
#                 mesh_one.get_element_with_frame_id(frame_one_id).global_id = global_id
#                 mesh_two.get_element_with_frame_id(large_mapping[frame_one_id]).global_id = global_id
#                 tracked_ids.append(global_id)
# 
#             mesh_one.index_global_ids()
#             mesh_two.index_global_ids()
# 
#             mesh_one.plot('mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
         
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_special_mesh_one_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_special_mesh_one_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
         
        # make sure that the entire mesh was tracked (except for the dead cell)
        self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() - 2 )

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(ground_truth[element_one.id_in_frame] , element_two.id_in_frame )
            if (ground_truth[element_one.id_in_frame] != element_two.id_in_frame ):
                print element_one.calculate_centroid()
                print element_two.calculate_centroid()
            
        plt.close('all')

    @attr(level = 'standard')
    def test_track_special_division_two(self):
        """read a special mesh, perform a division event, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','division_special_mesh_three.mesh'))
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        file_to_read = open(path.join(dirname(__file__),'data','division_special_direction_three.pickle'), 'r')
        division_direction = pickle.load(file_to_read)
        file_to_read.close()

        mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, division_direction)
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
#         largest_mappings = subgraph_finder.largest_mappings
# # 
#         for mapping_index, large_mapping in enumerate(largest_mappings):
#  
#             tracked_ids = []
#              
#             for element_one in mesh_one.elements:
#                 element_one.global_id = None
#  
#             for element_two in mesh_two.elements:
#                 element_two.global_id = None
#  
#             for global_id, frame_one_id in enumerate(large_mapping):
#                 mesh_one.get_element_with_frame_id(frame_one_id).global_id = global_id
#                 mesh_two.get_element_with_frame_id(large_mapping[frame_one_id]).global_id = global_id
#                 tracked_ids.append(global_id)
#  
#             mesh_one.index_global_ids()
#             mesh_two.index_global_ids()
#  
#             mesh_one.plot('debug_mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('debug_mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#          
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_special_mesh_two_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_special_mesh_two_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
          
        # make sure that the entire mesh was tracked (except for the dead cell)
        self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() - 2 )

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual( ground_truth[element_one.id_in_frame] , element_two.id_in_frame )
            if ( ground_truth[element_one.id_in_frame] != element_two.id_in_frame ):
                print 'hello?'
                print element_one.calculate_centroid()
                print element_two.calculate_centroid()
            
        plt.close('all')

    @attr(level = 'standard')
    def test_track_special_division_three(self):
        """read a special mesh, perform a division event, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','division_special_mesh_four.mesh'))
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_in_reverse_order()
        
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        file_to_read = open(path.join(dirname(__file__),'data','division_special_direction_four.pickle'), 'r')
        division_direction = pickle.load(file_to_read)
        file_to_read.close()

        mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, division_direction)
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
# #         tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
#         largest_mappings = subgraph_finder.largest_mappings
# 
#         for mapping_index, large_mapping in enumerate(largest_mappings):
#   
#             tracked_ids = []
#               
#             for element_one in mesh_one.elements:
#                 element_one.global_id = None
#   
#             for element_two in mesh_two.elements:
#                 element_two.global_id = None
#   
#             for global_id, frame_one_id in enumerate(large_mapping):
#                 mesh_one.get_element_with_frame_id(frame_one_id).global_id = global_id
#                 mesh_two.get_element_with_frame_id(large_mapping[frame_one_id]).global_id = global_id
#                 tracked_ids.append(global_id)
#   
#             mesh_one.index_global_ids()
#             mesh_two.index_global_ids()
#   
#             mesh_one.plot('debug_mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('debug_mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             if mapping_index == 2:
#                 import pdb; pdb.set_trace()
#           
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_special_mesh_three_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_special_mesh_three_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
           
        # make sure that the entire mesh was tracked (except for the dead cell)
        self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() - 1 )

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual( ground_truth[element_one.id_in_frame] , element_two.id_in_frame )
            
        plt.close('all')

    @attr(level = 'standard')
    def test_track_special_division_four(self):
        """read a special mesh, perform a division event, and track it.
           This is the mesh where a four-sided cell divides into two four-sided cells
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','division_special_mesh_five.mesh'))
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_in_reverse_order()
        
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        file_to_read = open(path.join(dirname(__file__),'data','division_special_direction_five.pickle'), 'r')
        division_direction = pickle.load(file_to_read)
        file_to_read.close()

        mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, division_direction)
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
#         subgraph_finder = tracking.LocalisedSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
#         largest_mappings = subgraph_finder.largest_mappings
#    
#         for mapping_index, large_mapping in enumerate(largest_mappings):
#      
#             tracked_ids = []
#                  
#             for element_one in mesh_one.elements:
#                 element_one.global_id = None
#      
#             for element_two in mesh_two.elements:
#                 element_two.global_id = None
#      
#             for global_id, frame_one_id in enumerate(large_mapping):
#                 mesh_one.get_element_with_frame_id(frame_one_id).global_id = global_id
#                 mesh_two.get_element_with_frame_id(large_mapping[frame_one_id]).global_id = global_id
#                 tracked_ids.append(global_id)
#      
#             mesh_one.index_global_ids()
#             mesh_two.index_global_ids()
#      
#             mesh_one.plot('debug_mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('debug_mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#              
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_special_mesh_four_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_special_mesh_four_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
#             
#         # make sure that the entire mesh was tracked (except for the dead cell)
#         self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() - 1 )
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual( ground_truth[element_one.id_in_frame] , element_two.id_in_frame )
             
        plt.close('all')

    @attr(level = 'standard')
    def test_track_special_division_five(self):
        """read a special mesh, perform a division event, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','division_special_mesh_six.mesh'))
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_in_reverse_order()
        
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        file_to_read = open(path.join(dirname(__file__),'data','division_special_direction_six.pickle'), 'r')
        division_direction = pickle.load(file_to_read)
        file_to_read.close()

        mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, division_direction)
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
#         largest_mappings = subgraph_finder.largest_mappings
#   
#         for mapping_index, large_mapping in enumerate(largest_mappings):
#     
#             tracked_ids = []
#                 
#             for element_one in mesh_one.elements:
#                 element_one.global_id = None
#     
#             for element_two in mesh_two.elements:
#                 element_two.global_id = None
#     
#             for global_id, frame_one_id in enumerate(large_mapping):
#                 mesh_one.get_element_with_frame_id(frame_one_id).global_id = global_id
#                 mesh_two.get_element_with_frame_id(large_mapping[frame_one_id]).global_id = global_id
#                 tracked_ids.append(global_id)
#     
#             mesh_one.index_global_ids()
#             mesh_two.index_global_ids()
#     
#             mesh_one.plot('debug_mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('debug_mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_special_mesh_five_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_special_mesh_five_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
             
        # make sure that the entire mesh was tracked (except for the dead cell)
        self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() - 1 )
  
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual( ground_truth[element_one.id_in_frame] , element_two.id_in_frame )
#              
        plt.close('all')

    @attr(level = 'standard')
    def test_track_special_division_six(self):
        """read a special mesh, perform a division event, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','division_special_mesh_seven.mesh'))
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_in_reverse_order()
        
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        file_to_read = open(path.join(dirname(__file__),'data','division_special_direction_seven.pickle'), 'r')
        division_direction = pickle.load(file_to_read)
        file_to_read.close()

        mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, division_direction)
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
#         largest_mappings = subgraph_finder.largest_mappings
#    
#         for mapping_index, large_mapping in enumerate(largest_mappings):
#      
#             tracked_ids = []
#                  
#             for element_one in mesh_one.elements:
#                 element_one.global_id = None
#      
#             for element_two in mesh_two.elements:
#                 element_two.global_id = None
#      
#             for global_id, frame_one_id in enumerate(large_mapping):
#                 mesh_one.get_element_with_frame_id(frame_one_id).global_id = global_id
#                 mesh_two.get_element_with_frame_id(large_mapping[frame_one_id]).global_id = global_id
#                 tracked_ids.append(global_id)
#      
#             mesh_one.index_global_ids()
#             mesh_two.index_global_ids()
#      
#             mesh_one.plot('debug_mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('debug_mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#              
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_special_mesh_six_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_special_mesh_six_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
             
        # make sure that the entire mesh was tracked (except for the dead cell)
        self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() - 1 )
  
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual( ground_truth[element_one.id_in_frame] , element_two.id_in_frame )
#              
        plt.close('all')

    @attr(level = 'known_to_fail')
    def test_track_special_division_seven(self):
        """read a special mesh, perform a division event, and track it.
           This mesh is too small to find an initial mapping.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','division_special_mesh_eight.mesh'))
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_in_reverse_order()
        
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        file_to_read = open(path.join(dirname(__file__),'data','division_special_direction_eight.pickle'), 'r')
        division_direction = pickle.load(file_to_read)
        file_to_read.close()

        mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, division_direction)
        
#         mesh_one.plot(path.join(dirname(__file__),'output','tracked_special_mesh_seven_before_division.pdf') )
#         mesh_two.plot(path.join(dirname(__file__),'output','tracked_special_mesh_seven_after_division.pdf'))
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
#         largest_mappings = subgraph_finder.largest_mappings
#    
#         for mapping_index, large_mapping in enumerate(largest_mappings):
#      
#             tracked_ids = []
#                  
#             for element_one in mesh_one.elements:
#                 element_one.global_id = None
#      
#             for element_two in mesh_two.elements:
#                 element_two.global_id = None
#      
#             for global_id, frame_one_id in enumerate(large_mapping):
#                 mesh_one.get_element_with_frame_id(frame_one_id).global_id = global_id
#                 mesh_two.get_element_with_frame_id(large_mapping[frame_one_id]).global_id = global_id
#                 tracked_ids.append(global_id)
#      
#             mesh_one.index_global_ids()
#             mesh_two.index_global_ids()
#      
#             mesh_one.plot('debug_mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('debug_mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#              
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_special_mesh_seven_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_special_mesh_seven_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
             
        # make sure that the entire mesh was tracked (except for the dead cell)
        self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() - 1 )
  
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual( ground_truth[element_one.id_in_frame] , element_two.id_in_frame )
#              
        plt.close('all')

    @attr(level = 'known_to_fail')
    def test_track_special_division_eight(self):
        """read a special mesh, perform a division event, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','division_special_mesh_nine.mesh'))
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_in_reverse_order()
        
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        file_to_read = open(path.join(dirname(__file__),'data','division_special_direction_nine.pickle'), 'r')
        division_direction = pickle.load(file_to_read)
        file_to_read.close()

        mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, division_direction)
        
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_special_mesh_eight_before_division.pdf') )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_special_mesh_eight_after_division.pdf'))
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
#         largest_mappings = subgraph_finder.largest_mappings
#    
#         for mapping_index, large_mapping in enumerate(largest_mappings):
#      
#             tracked_ids = []
#                  
#             for element_one in mesh_one.elements:
#                 element_one.global_id = None
#      
#             for element_two in mesh_two.elements:
#                 element_two.global_id = None
#      
#             for global_id, frame_one_id in enumerate(large_mapping):
#                 mesh_one.get_element_with_frame_id(frame_one_id).global_id = global_id
#                 mesh_two.get_element_with_frame_id(large_mapping[frame_one_id]).global_id = global_id
#                 tracked_ids.append(global_id)
#      
#             mesh_one.index_global_ids()
#             mesh_two.index_global_ids()
#      
#             mesh_one.plot('debug_mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('debug_mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#              
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_special_mesh_eight_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_special_mesh_eight_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
             
        # make sure that the entire mesh was tracked (except for the dead cell)
        self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() - 1 )
  
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual( ground_truth[element_one.id_in_frame] , element_two.id_in_frame )
#              
        plt.close('all')

    @attr(level = 'standard')
    def xest_track_division(self):
        """generate a random mesh, copy it, perform a division event, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.creation.generate_random_tesselation(9,9)
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.save(path.join(dirname(__file__),'output','division_mesh_one.mesh'))
        mesh_two.save(path.join(dirname(__file__),'output','division_mesh_two.mesh'))

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        division_angle = np.random.uniform(0.0, 2.0*np.pi)
        random_direction = np.array( [np.cos(division_angle), np.sin(division_angle)] )

        file_to_write = open(path.join(dirname(__file__),'output','last_division_direction.pickle'), 'wb')
        pickle.dump(random_direction, file_to_write)
        file_to_write.close()

        mesh_two.divide_element_with_frame_id_in_direction(most_central_element.id_in_frame, random_direction)
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
 
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_mesh_before_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_mesh_after_division.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
         
        # make sure that the entire mesh was tracked (except for the dead cell)

        dangling_elements_one = mesh_one.get_dangling_element_ids(ground_truth.keys()) 
        self.assertGreaterEqual(len(tracked_ids), 
                         max(mesh_one.get_num_elements() - 4, 
                             mesh_one.get_num_elements() - len(dangling_elements_one) - 1))

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )
            
        plt.close('all')
        
    @attr(level = 'weekly')
    def test_track_division_multiple_times(self):
        for n in range(100):
            print 'division run ' + str(n)
            self.xest_track_division()
