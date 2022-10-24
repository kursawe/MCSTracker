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
import copy
import sys
import numpy as np
import matplotlib.pyplot as plt
from os import path
from os.path import dirname
from nose.plugins.attrib import attr

class TestTrackDeath(unittest.TestCase):
                                 
    @attr(level = 'standard')
    def test_maximum_common_subgraph_for_cell_death(self):
        """generate a random mesh, copy it, perform a death event, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.creation.generate_random_tesselation(9,9)
        mesh_one.save(path.join(dirname(__file__),'output','subgraph_death_mesh_one.mesh'))
        mesh_two = copy.deepcopy(mesh_one)
        mesh_one.save(path.join(dirname(__file__),'output','subgraph_death_mesh_two.mesh'))
#         mesh_one = mesh.load('special_death_mesh_one.mesh')
#         mesh_two = mesh.load('special_death_mesh_two.mesh')

        # First pick a cell in the centre
        mesh_centre = mesh_two.calculate_centre()
        
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # pick the element closest to the centre
        
        most_central_element = mesh_two.find_most_central_element()
               
        mesh_two.kill_element_with_frame_id(most_central_element.id_in_frame)

        tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
 
        mesh_one.plot(path.join(dirname(__file__),'output','mesh_before_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','mesh_after_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
         
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
     
    @attr(level = 'standard')
    def test_track_special_case(self):
        """track a cell death event and identify the cells that are involved in rearrangement"""

        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_one.mesh'))
        mesh_two = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_two.mesh'))

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        mesh_two.kill_element_with_frame_id(most_central_element.id_in_frame)

        mesh_one.plot(path.join(dirname(__file__),'output','not_tracked_mesh_before_death.pdf' ))
        mesh_two.plot(path.join(dirname(__file__),'output','not_tracked_mesh_after_death.pdf' ))

        tracked_ids = tracking.track( mesh_one, mesh_two )
#   
#         mesh_one.plot('special_tracked_mesh_before_death.pdf', color_by_global_id = True, 
#                       total_number_of_global_ids = len( tracked_ids ) )
#         mesh_two.plot('special_tracked_mesh_after_death.pdf', color_by_global_id = True, 
#                       total_number_of_global_ids = len( tracked_ids ) )

#         tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
#         subgraph_finder = tracking.LocalisedSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
#         tracked_ids = tracking.find_maximum_common_subgraph(mesh_one, mesh_two)

        mesh_one.plot(path.join(dirname(__file__),'output','special_tracked_mesh_before_death.pdf'), color_by_global_id = True, 
                    total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','special_tracked_mesh_after_death.pdf'), color_by_global_id = True, 
                    total_number_of_global_ids = len( tracked_ids ) )

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
#             mesh_one.plot('mesh_special_before_death_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('mesh_special_after_death_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#   


        # make sure that the entire mesh was tracked (except for the dead cell)
        self.assertEqual( len(tracked_ids), mesh_two.get_num_elements() - 2 )

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )
            
        plt.close('all')

    @attr(level = 'standard')
    def test_track_special_case_one(self):
        """track a cell death event and identify the cells that are involved in rearrangement"""

        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_1_one.mesh'))
        mesh_two = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_1_two.mesh'))

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        mesh_two.kill_element_with_frame_id(most_central_element.id_in_frame)

        mesh_one.plot(path.join(dirname(__file__),'output'),'not_tracked_mesh_before_death.pdf' )
        mesh_two.plot(path.join(dirname(__file__),'output'),'not_tracked_mesh_after_death.pdf' )

        tracked_ids = tracking.track( mesh_one, mesh_two )
#   
#         mesh_one.plot('special_tracked_mesh_before_death.pdf', color_by_global_id = True, 
#                       total_number_of_global_ids = len( tracked_ids ) )
#         mesh_two.plot('special_tracked_mesh_after_death.pdf', color_by_global_id = True, 
#                       total_number_of_global_ids = len( tracked_ids ) )

#         tracked_ids = tracking.find_maximum_common_subgraph( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()

        mesh_one.plot(path.join(dirname(__file__),'output','special_tracked_mesh_1_before_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','special_tracked_mesh_1_after_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )

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
#             mesh_one.plot('mesh_special_before_death_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('mesh_special_after_death_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#  


        # make sure that the entire mesh was tracked (except for the dead cell)
        self.assertEqual( len(tracked_ids), mesh_two.get_num_elements() - 1 )

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )
            
        plt.close('all')

    @attr(level = 'standard')
    def test_track_special_case_two(self):
        """generate a random mesh, copy it, perform a death event, and find maximum common subgraph
           not sure whether this one even works!
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','multiple_subgraphs_one.mesh'))
        mesh_two = mesh.load(path.join(dirname(__file__),'data','multiple_subgraphs_two.mesh'))

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_in_reverse_order()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        mesh_two.kill_element_with_frame_id(most_central_element.id_in_frame)
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
# 
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
#             mesh_one.plot('mesh_special_before_death_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                         total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('mesh_special_after_death_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                         total_number_of_global_ids = len( tracked_ids ) )

        mesh_one.plot(path.join(dirname(__file__),'output','special_tracked_mesh_2_before_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','special_tracked_mesh_2_after_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )

        self.assertEqual( len(tracked_ids), mesh_two.get_num_elements() )

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )
            
        plt.close('all')

    @attr(level = 'weekly')
    def test_track_special_case_three(self):
        """generate a random mesh, copy it, perform a death event, and find maximum common subgraph
           not sure whether this one even works!
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_three.mesh'))
        mesh_two = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_three.mesh'))

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        mesh_two.kill_element_with_frame_id(most_central_element.id_in_frame)

        tracked_ids = tracking.track( mesh_one, mesh_two )
# 
        mesh_one.plot(path.join(dirname(__file__),'output','special_tracked_mesh_three_before_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','special_tracked_mesh_three_after_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
 
        self.assertEqual( len(tracked_ids), mesh_two.get_num_elements() )
  
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


        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )
            
        plt.close('all')

    @attr(level = 'known_to_fail')
    def test_track_special_case_four(self):
        """generate a random mesh, copy it, perform a death event, and find maximum common subgraph
           not sure whether this one even works!
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_four.mesh'))
        mesh_two = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_four.mesh'))

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        mesh_two.kill_element_with_frame_id(most_central_element.id_in_frame)

        mesh_one.plot(path.join(dirname(__file__),'output','special_mesh_four_before_death_before_tracking.pdf'))
        mesh_two.plot(path.join(dirname(__file__),'output','special_mesh_four_after_death_before_tracking.pdf'))
        tracked_ids = tracking.track( mesh_one, mesh_two )
# 
        mesh_one.plot(path.join(dirname(__file__),'output','special_tracked_mesh_four_before_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','special_tracked_mesh_four_after_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
 
        self.assertEqual( len(tracked_ids), mesh_two.get_num_elements() )
  
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


        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )
            
        plt.close('all')

    @attr(level = 'standard')
    def test_track_special_case_five(self):
        """generate a random mesh, copy it, perform a death event, and find maximum common subgraph
           not sure whether this one even works!
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_five.mesh'))
        mesh_two = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_five.mesh'))

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_in_reverse_order()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        mesh_two.kill_element_with_frame_id(most_central_element.id_in_frame)
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
# 
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
#             mesh_one.plot('mesh_special_before_death_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                         total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('mesh_special_after_death_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                         total_number_of_global_ids = len( tracked_ids ) )

        mesh_one.plot(path.join(dirname(__file__),'output','special_tracked_mesh_5_before_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','special_tracked_mesh_5_after_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )

        dangling_elements_one = mesh_one.get_dangling_element_ids(ground_truth.keys()) 
        self.assertGreaterEqual(len(tracked_ids), 
                         max(mesh_one.get_num_elements() - 4, 
                             mesh_one.get_num_elements() - len(dangling_elements_one) -1 ))

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )
            
        plt.close('all')

    @attr(level = 'standard')
    def test_track_special_case_six(self):
        """generate a random mesh, copy it, perform a death event, and find maximum common subgraph
           not sure whether this one even works!
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_six.mesh'))
        mesh_two = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_six.mesh'))

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_in_reverse_order()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        mesh_two.kill_element_with_frame_id(most_central_element.id_in_frame)
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
# 
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
#             mesh_one.plot('mesh_special_before_death_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                         total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('mesh_special_after_death_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                         total_number_of_global_ids = len( tracked_ids ) )

        mesh_one.plot(path.join(dirname(__file__),'output','special_tracked_mesh_6_before_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','special_tracked_mesh_6_after_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )

        dangling_elements_one = mesh_one.get_dangling_element_ids(ground_truth.keys()) 
        self.assertGreaterEqual(len(tracked_ids), 
                         max(mesh_one.get_num_elements() - 4, 
                             mesh_one.get_num_elements() - len(dangling_elements_one) -1))

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )
            
        plt.close('all')

    @attr(level = 'known_to_fail')
    def test_track_special_case_seven(self):
        """generate a random mesh, copy it, perform a death event, and find maximum common subgraph
           not sure whether this one even works!
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_seven.mesh'))
        mesh_two = mesh.load(path.join(dirname(__file__),'data','special_death_mesh_seven.mesh'))

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_in_reverse_order()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        mesh_two.kill_element_with_frame_id(most_central_element.id_in_frame)
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         subgraph_finder = tracking.ReducedBacktrackingSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
# 
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
#             mesh_one.plot('mesh_special_before_death_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                         total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('mesh_special_after_death_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                         total_number_of_global_ids = len( tracked_ids ) )

        mesh_one.plot(path.join(dirname(__file__),'output','special_tracked_mesh_7_before_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','special_tracked_mesh_7_after_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )

        dangling_elements_one = mesh_one.get_dangling_element_ids(ground_truth.keys()) 
        self.assertGreaterEqual(len(tracked_ids), 
                         max(mesh_one.get_num_elements() - 4, 
                             mesh_one.get_num_elements() - len(dangling_elements_one) -1))

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )
            
        plt.close('all')


    @attr(level = 'standard')
    def xest_track_cell_death(self):
        """track a cell death event and identify the cells that are involved in rearrangement"""

        sys.setrecursionlimit(40000)

        mesh_one = mesh.creation.generate_random_tesselation(9,9)
        mesh_one.save(path.join(dirname(__file__),'output','death_mesh_one.mesh'))
        mesh_two = copy.deepcopy(mesh_one)
#        mesh_one = mesh.load('death_mesh_one.mesh')
#        mesh_two = mesh.load('death_mesh_two.mesh')

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # pick the element closest to the centre
        most_central_element = mesh_two.find_most_central_element()
               
        mesh_two.kill_element_with_frame_id(most_central_element.id_in_frame)

        mesh_one.plot(path.join(dirname(__file__),'output','not_tracked_mesh_before_death.pdf' ))
        mesh_two.plot(path.join(dirname(__file__),'output','not_tracked_mesh_after_death.pdf' ))

        tracked_ids = tracking.track( mesh_one, mesh_two )
 
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_mesh_before_death.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_mesh_after_death.pdf'), color_by_global_id = True, 
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
    def test_track_death_multiple_times(self):
        for n in range(100):
            print('death run ' + str(n))
            self.xest_track_cell_death()
