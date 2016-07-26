# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This tests our first tracking example
"""
import unittest
import mesh
import tracking
import copy
import sys
import numpy as np
import matplotlib.pyplot as plt
import make_t1_analysis
from os import path
from os.path import dirname
from nose.plugins.attrib import attr

class TestTrackT1Swap(unittest.TestCase):
                                 
    @attr(level = 'standard')
    def test_maximum_common_subgraph_for_t1_swap(self):
        """generate a random mesh, copy it, perform a t1 swap on it, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.creation.generate_random_tesselation(8,8)
         
        mesh_two = copy.deepcopy(mesh_one)
 
        mesh_one.save(path.join(dirname(__file__),'output','subgraph_t1_mesh_one.mesh'))
        mesh_one.save(path.join(dirname(__file__),'output','subgraph_t1_mesh_two.mesh'))

#         mesh_one = mesh.load('mesh_t1_special_one.mesh')
#         mesh_two = mesh.load('mesh_t1_special_two.mesh')

        # Perform T1 swap on mesh two 
        # First pick a node in the centre
        mesh_centre = mesh_two.calculate_centre()
        
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # pick the node closest to the centre
        min_distance = 3*mesh_two.calculate_height()

        for node in mesh_two.nodes:
            distance = np.linalg.norm(node.position - mesh_centre)
            if distance < min_distance:
               min_distance = distance
               most_central_node = node 
               
        for local_index, element_node in enumerate(most_central_node.adjacent_elements[0].nodes):
            if element_node.id == most_central_node.id:
                num_nodes_this_element = most_central_node.adjacent_elements[0].get_num_nodes()
                one_edge_node = most_central_node.adjacent_elements[0].nodes[(local_index+1)%num_nodes_this_element]
                break
        
        mesh_two.perform_t1_swap( most_central_node.id, one_edge_node.id )
         
        tracked_ids = tracking.find_maximum_common_subgraph(mesh_one, mesh_two)
 
        mesh_one.plot(path.join(dirname(__file__),'output','mesh_before_t1_swap.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len( tracked_ids ) )
        mesh_two.plot(path.join(dirname(__file__),'output','mesh_after_t1_swap.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
         
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())

        plt.close('all')

#     def test_track_special_t1_swap(self):
#         """generate a random mesh, copy it, perform a t1 swap on it, and track it.
#         """
#         sys.setrecursionlimit(40000)
# 
# #         mesh_one = mesh.creation.generate_random_tesselation(6,6)
# #         
# #         mesh_two = copy.deepcopy(mesh_one)
# # 
# #         mesh_one.assign_frame_ids_in_order()
# #         mesh_two.assign_frame_ids_randomly()
# # 
# #         mesh_one.save('mesh_t1_one.mesh')
# #         mesh_one.save('mesh_t1_two.mesh')
# 
#         mesh_one = mesh.load('mesh_t1_special_one.mesh')
#         mesh_two = mesh.load('mesh_t1_special_two.mesh')
# 
#         # build ground truth for testing the mapping
#         ground_truth = {}
#         for element_index, element in enumerate(mesh_one.elements):
#             ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame
# 
#         # Perform T1 swap on mesh two 
#         # First pick a node in the centre
#         mesh_centre = mesh_two.calculate_centre()
# 
#         # pick the node closest to the centre
#         min_distance = 3*mesh_two.calculate_height()
#         for node in mesh_two.nodes:
#             distance = np.linalg.norm(node.position - mesh_centre)
#             if distance < min_distance:
#                min_distance = distance
#                most_central_node = node 
#                
#         # pick a node that shares an edge with this central node
#         for local_index, element_node in enumerate(most_central_node.adjacent_elements[0].nodes):
#             if element_node.id == most_central_node.id:
#                 num_nodes_this_element = most_central_node.adjacent_elements[0].get_num_nodes()
#                 one_edge_node = most_central_node.adjacent_elements[0].nodes[(local_index+1)%num_nodes_this_element]
#                 break
#         
#         mesh_two.perform_t1_swap( most_central_node.id, one_edge_node.id )
#          
#         mesh_one.plot('tracked_mesh_before_t1_swap.pdf')
#         mesh_two.plot('tracked_mesh_after_t1_swap.pdf')
#         
#         tracked_ids = tracking.track( mesh_one, mesh_two )
# 
#         mesh_one.plot('tracked_mesh_before_t1_swap.pdf', color_by_global_id = True, 
#                       total_number_of_global_ids = len(tracked_ids) )
#         mesh_two.plot('tracked_mesh_after_t1_swap.pdf', color_by_global_id = True, 
#                       total_number_of_global_ids = len(tracked_ids) )
#          
#         # make sure that the entire mesh was tracked
#         self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() )
# 
#         for global_id in tracked_ids:
#             # and that the mapping coincides with the ground truth for all tracked ids
#             element_one = mesh_one.get_element_with_global_id(global_id)
#             element_two = mesh_two.get_element_with_global_id(global_id)
#             self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )
# 
#         plt.close('all')

    @attr(level = 'standard')
    def test_track_special_t1_swap(self):
        """generate a random mesh, copy it, perform a t1 swap on it, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__),'data','special_t1_mesh.mesh'))
        
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_in_reverse_order()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # Perform T1 swap on mesh two 
        # First pick a node in the centre
        mesh_centre = mesh_two.calculate_centre()

        # pick the node closest to the centre
        min_distance = 3*mesh_two.calculate_height()
        for node in mesh_two.nodes:
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
        
        mesh_two.perform_t1_swap( most_central_node.id, one_edge_node.id )
         
        tracked_ids = tracking.track( mesh_one, mesh_two )
  
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
#             mesh_one.plot('mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#   
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_special_mesh_before_t1_swap.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_special_mesh_after_t1_swap.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
          
        # make sure that the entire mesh was tracked
        dangling_elements_one = mesh_one.get_dangling_element_ids(ground_truth.keys()) 
        self.assertGreaterEqual(len(tracked_ids), 
                         max(mesh_one.get_num_elements() - 3, 
                             mesh_one.get_num_elements() - len(dangling_elements_one)))

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

        plt.close('all')

    @attr(level = 'standard')
    def xest_track_t1_swap(self):
        """generate a random mesh, copy it, perform a t1 swap on it, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.creation.generate_random_tesselation(9,9)
        
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        mesh_one.save(path.join(dirname(__file__),'output','mesh_t1_one.mesh'))
        mesh_one.save(path.join(dirname(__file__),'output','mesh_t1_two.mesh'))

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # Perform T1 swap on mesh two 
        # First pick a node in the centre
        mesh_centre = mesh_two.calculate_centre()

        # pick the node closest to the centre
        min_distance = 3*mesh_two.calculate_height()
        for node in mesh_two.nodes:
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
        
        mesh_two.perform_t1_swap( most_central_node.id, one_edge_node.id )
         
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_mesh_before_t1_swap.pdf'))
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_mesh_after_t1_swap.pdf'))
        
        tracked_ids = tracking.track( mesh_one, mesh_two )

        mesh_one.plot(path.join(dirname(__file__),'output','tracked_mesh_before_t1_swap.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_mesh_after_t1_swap.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
         
        # make sure that the entire mesh was tracked
        dangling_elements_one = mesh_one.get_dangling_element_ids(ground_truth.keys()) 
        self.assertGreaterEqual(len(tracked_ids), 
                         max(mesh_one.get_num_elements() - 3, 
                             mesh_one.get_num_elements() - len(dangling_elements_one)))

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

        plt.close('all')

    @attr(level = 'standard')
    def test_track_special_t1_swap_2(self):
        """generate a random mesh, copy it, perform a t1 swap on it, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__), 'data', 'special_t1_mesh_two.mesh'))
        
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # Perform T1 swap on mesh two 
        # First pick a node in the centre
        mesh_centre = mesh_two.calculate_centre()

        # pick the node closest to the centre
        min_distance = 3*mesh_two.calculate_height()
        for node in mesh_two.nodes:
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
        
        mesh_two.perform_t1_swap( most_central_node.id, one_edge_node.id )
         
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_mesh_before_special_t1_swap_two.pdf'))
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_mesh_after_special_t1_swap_two.pdf'))
        
        tracked_ids = tracking.track( mesh_one, mesh_two )
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
#             mesh_one.plot('mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
# 
#         mesh_one.plot(path.join(dirname(__file__),'output','tracked_mesh_before_special_t1_swap_two.pdf'), color_by_global_id = True, 
#                       total_number_of_global_ids = len(tracked_ids) )
#         mesh_two.plot(path.join(dirname(__file__),'output','tracked_mesh_after_special_t1_swap_two.pdf'), color_by_global_id = True, 
#                       total_number_of_global_ids = len(tracked_ids) )
#          
#         # make sure that the entire mesh was tracked
#         self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() )

        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

        plt.close('all')

    @attr(level = 'standard')
    def test_track_special_t1_swap_3(self):
        """generate a random mesh, copy it, perform a t1 swap on it, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__), 'data', 'special_t1_mesh_three.mesh'))
        
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # Perform T1 swap on mesh two 
        # First pick a node in the centre
        mesh_centre = mesh_two.calculate_centre()

        # pick the node closest to the centre
        min_distance = 3*mesh_two.calculate_height()
        for node in mesh_two.nodes:
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
        
        mesh_two.perform_t1_swap( most_central_node.id, one_edge_node.id )
         
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         subgraph_finder = tracking.LocalisedSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
#         post_processor = tracking.PostProcessor(mesh_one, mesh_two, subgraph_finder.largest_mappings)
#         post_processor.tidy_current_mapping()
#         post_processor.index_global_ids_from_largest_mappings()
#         network_one = mesh_one.generate_network_of_unidentified_elements()
#         print network_one.nodes()
#         post_processor.altered_fill_in_by_adjacency(network_one)
#         post_processor.index_global_ids()
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
#             mesh_one.plot('mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#  
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_mesh_before_special_t1_swap_three.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = mesh_one.get_num_elements() )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_mesh_after_special_t1_swap_three.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = mesh_one.get_num_elements() )
           
#         # make sure that the entire mesh was tracked
#         self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() )

#         mesh_one.plot(path.join(dirname(__file__),'output','tracked_mesh_before_special_t1_swap_three.pdf'))
#         mesh_two.plot(path.join(dirname(__file__),'output','tracked_mesh_after_special_t1_swap_three.pdf'))
#         
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

        plt.close('all')

    @attr(level = 'standard')
    def test_track_special_t1_swap_4(self):
        """generate a random mesh, copy it, perform a t1 swap on it, and track it.
        """
        sys.setrecursionlimit(40000)

        mesh_one = mesh.load(path.join(dirname(__file__), 'data', 'special_t1_mesh_four.mesh'))
        
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()

        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame

        # Perform T1 swap on mesh two 
        # First pick a node in the centre
        mesh_centre = mesh_two.calculate_centre()

        # pick the node closest to the centre
        min_distance = 3*mesh_two.calculate_height()
        for node in mesh_two.nodes:
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
        
        mesh_two.perform_t1_swap( most_central_node.id, one_edge_node.id )
         
        tracked_ids = tracking.track( mesh_one, mesh_two )
#         subgraph_finder = tracking.LocalisedSubgraphFinder(mesh_one, mesh_two)
#         subgraph_finder.find_maximum_common_subgraph()
#         post_processor = tracking.PostProcessor(mesh_one, mesh_two, subgraph_finder.largest_mappings)
#         post_processor.tidy_current_mapping()
#         post_processor.index_global_ids_from_largest_mappings()
#         network_one = mesh_one.generate_network_of_unidentified_elements()
#         print network_one.nodes()
#         post_processor.altered_fill_in_by_adjacency(network_one)
#         post_processor.index_global_ids()
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
#             mesh_one.plot('mesh_special_before_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('mesh_special_after_division_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#  
        mesh_one.plot(path.join(dirname(__file__),'output','tracked_mesh_before_special_t1_swap_four.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = mesh_one.get_num_elements() )
        mesh_two.plot(path.join(dirname(__file__),'output','tracked_mesh_after_special_t1_swap_four.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = mesh_one.get_num_elements() )
           
#         # make sure that the entire mesh was tracked
#         self.assertEqual( len(tracked_ids), mesh_one.get_num_elements() )

#         mesh_one.plot(path.join(dirname(__file__),'output','tracked_mesh_before_special_t1_swap_three.pdf'))
#         mesh_two.plot(path.join(dirname(__file__),'output','tracked_mesh_after_special_t1_swap_three.pdf'))
#         
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

        plt.close('all')

    @attr(level = 'known_to_fail')
    def test_track_multiple_t1_swaps_special_example(self):

        sys.setrecursionlimit(40000)
 
        first_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_before.mesh'))
        second_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_after.mesh'))
        
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(first_mesh.elements):
            ground_truth[element.id_in_frame] = second_mesh.elements[element_index].id_in_frame
        
        tracking.track(first_mesh, second_mesh)
#         tracking.find_maximum_common_subgraph(first_mesh, second_mesh)
    
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_before.pdf'), 
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_after.pdf'),
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
        
        tracking_success, number_tracked_cells = tracking.evaluate_tracking(first_mesh, 
                                                                   second_mesh, 
                                                                   ground_truth)
        
    @attr(level = 'known_to_fail')
    def test_track_multiple_t1_swaps_special_example_2(self):

        sys.setrecursionlimit(40000)
 
        first_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_example_2_before.mesh'))
        second_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_example_2_after.mesh'))
        
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(first_mesh.elements):
            ground_truth[element.id_in_frame] = second_mesh.elements[element_index].id_in_frame
        
        tracking.track(first_mesh, second_mesh)
#         tracking.find_maximum_common_subgraph(first_mesh, second_mesh)
    
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_2_before.pdf'), 
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_2_after.pdf'),
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
        
        tracking_success, number_tracked_cells = tracking.evaluate_tracking(first_mesh, 
                                                                   second_mesh, 
                                                                   ground_truth)
 
    @attr(level = 'known_to_fail')
    def test_track_multiple_t1_swaps_special_example_3(self):

        sys.setrecursionlimit(40000)
 
        first_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_example_3_before.mesh'))
        second_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_example_3_after.mesh'))
        
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(first_mesh.elements):
            ground_truth[element.id_in_frame] = second_mesh.elements[element_index].id_in_frame
        
        tracking.track(first_mesh, second_mesh)
#         tracking.find_maximum_common_subgraph(first_mesh, second_mesh)
    
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_3_before.pdf'), 
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_3_after.pdf'),
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
        
        tracking_success, number_tracked_cells = tracking.evaluate_tracking(first_mesh, 
                                                                   second_mesh, 
                                                                   ground_truth)

    @attr(level = 'known_to_fail')
    def test_track_multiple_t1_swaps_special_example_4(self):

        sys.setrecursionlimit(40000)
 
        first_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_example_4_before.mesh'))
        second_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_example_4_after.mesh'))
         
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(first_mesh.elements):
            ground_truth[element.id_in_frame] = second_mesh.elements[element_index].id_in_frame
         
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_4_before.pdf'))
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_4_after.pdf'))
 
        tracking.track(first_mesh, second_mesh)
#         tracking.find_maximum_common_subgraph(first_mesh, second_mesh)
     
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_4_before.pdf'), 
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_4_after.pdf'),
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
         
        tracking_success, number_tracked_cells = tracking.evaluate_tracking(first_mesh, 
                                                                   second_mesh, 
                                                                   ground_truth)
  
    @attr(level = 'known_to_fail')
    def test_track_multiple_t1_swaps_special_example_5(self):

        sys.setrecursionlimit(40000)
 
        first_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_example_5_before.mesh'))
        second_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_example_5_after.mesh'))
         
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(first_mesh.elements):
            ground_truth[element.id_in_frame] = second_mesh.elements[element_index].id_in_frame
         
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_5_before.pdf'))
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_5_after.pdf'))
 
        tracking.track(first_mesh, second_mesh)
#         tracking.find_maximum_common_subgraph(first_mesh, second_mesh)
     
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_5_before.pdf'), 
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_5_after.pdf'),
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
         
        print 'evaluate tracking'
        tracking_success, number_tracked_cells = tracking.evaluate_tracking(first_mesh, 
                                                                   second_mesh, 
                                                                   ground_truth)
 
    @attr(level = 'known_to_fail')
    def test_track_multiple_t1_swaps_special_example_6(self):

        sys.setrecursionlimit(40000)
 
        first_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_example_6_before.mesh'))
        second_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_example_6_after.mesh'))
         
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(first_mesh.elements):
            ground_truth[element.id_in_frame] = second_mesh.elements[element_index].id_in_frame
         
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_6_before.pdf'))
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_6_after.pdf'))
 
        tracking.track(first_mesh, second_mesh)
#         tracking.find_maximum_common_subgraph(first_mesh, second_mesh)
     
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_6_before.pdf'), 
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_6_after.pdf'),
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
         
        print 'evaluate tracking'
        tracking_success, number_tracked_cells = tracking.evaluate_tracking(first_mesh, 
                                                                   second_mesh, 
                                                                   ground_truth)
    @attr(level = 'known_to_fail')
    def test_track_multiple_t1_swaps_special_example_7(self):

        sys.setrecursionlimit(40000)
 
        first_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_example_7_before.mesh'))
        second_mesh = mesh.load(path.join(dirname(__file__), 'data', 'multiple_t1_example_7_after.mesh'))
         
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(first_mesh.elements):
            ground_truth[element.id_in_frame] = second_mesh.elements[element_index].id_in_frame
         
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_7_before.pdf'))
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_7_after.pdf'))
 
        tracking.track(first_mesh, second_mesh)
#         tracking.find_maximum_common_subgraph(first_mesh, second_mesh)
     
        first_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_7_before.pdf'), 
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
        second_mesh.plot(path.join(dirname(__file__),'output','multiple_t1_example_7_after.pdf'),
                        color_by_global_id = True, total_number_of_global_ids = first_mesh.get_num_elements())
         
        print 'evaluate tracking'
        tracking_success, number_tracked_cells = tracking.evaluate_tracking(first_mesh, 
                                                                   second_mesh, 
                                                                   ground_truth)
 
    @attr(level = 'weekly')
    def test_track_multiple_t1_swaps(self):

        sys.setrecursionlimit(40000)
 
        """Here, we run perform the tracking after multiple t1 swaps have been executed.
           We only test that the tracking algorithm actually finishes."""
           
        for n in range(150):
            print 'multiple t1 swap run ' + str(n)
            make_t1_analysis.test_success(20, output = True )
        
    @attr(level = 'weekly')
    def test_track_t1_swap_multiple_times(self):
        for n in range(100):
            print 't1 swap run ' + str(n)
            self.xest_track_t1_swap()  
