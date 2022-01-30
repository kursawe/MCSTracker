# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This tests our first tracking example
"""
import unittest
import mesh
import tracking
import sys
from os import path
from os.path import dirname
from nose.plugins.attrib import attr
import matplotlib.pyplot as plt

class TestTrackGlobalMovements(unittest.TestCase):
                                 
    @attr(level = 'standard')
    def xest_translation(self):
        """generate a random mesh, move it, and track it.
        """
        sys.setrecursionlimit(40000)
        large_mesh = mesh.creation.generate_random_tesselation(15,8)
        large_mesh.save(path.join(dirname(__file__),'output','large_transation_mesh.mesh'))

        large_mesh.assign_frame_ids_in_order()
        
        mesh_one = large_mesh.copy_region_of_interest((3,10),(-1,9))
        
        mesh_two = large_mesh.copy_region_of_interest((5,12),(-1,9))
        
        # Build ground truth
        ground_truth_indices = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element.id_in_frame in mesh_two.frame_id_dictionary:
                ground_truth_indices[element_index] = mesh_two.frame_id_dictionary[element.id_in_frame]
                
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
 
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element_index in ground_truth_indices:
                ground_truth[element.id_in_frame] = mesh_two.elements[ground_truth_indices[element_index]].id_in_frame
        
        mesh_one.plot(path.join(dirname(__file__),'output','mesh_before_translation.pdf'))
        mesh_two.plot(path.join(dirname(__file__),'output','mesh_after_translation.pdf'))
 
        print('call tracker')
        tracked_ids = tracking.track(mesh_one, mesh_two)
        print('tracker returned')
  
        mesh_one.plot(path.join(dirname(__file__),'output','larger_mesh_before_translation'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
        mesh_two.plot(path.join(dirname(__file__),'output','larger_mesh_after_translation'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
         
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()

        # Find any dangling elements that may not be identifiable through adjacency
        dangling_elements_counter = 0
        dangling_elements = []
        for element_one_id in ground_truth:
            element_one = mesh_one.get_element_with_frame_id(element_one_id)
            element_two = mesh_two.get_element_with_frame_id(ground_truth[element_one_id])
            ids_adjacent_to_element_one = element_one.get_ids_of_adjacent_elements()
            ids_adjacent_to_element_two = element_two.get_ids_of_adjacent_elements()
            mappable_elements_adjacent_to_elements_one = set(ids_adjacent_to_element_one).intersection(ground_truth.keys())
            mappable_elements_adjacent_to_elements_two = set(ids_adjacent_to_element_two).intersection(ground_truth.values())
            if (len(mappable_elements_adjacent_to_elements_one) == 1 or 
                len(mappable_elements_adjacent_to_elements_two) == 1):
                dangling_elements_counter += 1
                dangling_elements.append(element_one_id)
            elif (len(mappable_elements_adjacent_to_elements_one) == 2):
                neighbouring_polygon_number_is_the_same = False
                for element in mappable_elements_adjacent_to_elements_one:
                    element_polygon_number = mesh_one.get_element_with_frame_id(element).get_num_nodes()
                    if element_polygon_number == element_one.get_num_nodes():
                        dangling_elements.append(element_one_id)
                        dangling_elements_counter += 1
                        break
            elif (len(mappable_elements_adjacent_to_elements_two) == 2):
                neighbouring_polygon_number_is_the_same = False
                for element in mappable_elements_adjacent_to_elements_two:
                    element_polygon_number = mesh_two.get_element_with_frame_id(element).get_num_nodes()
                    if element_polygon_number == element_two.get_num_nodes():
                        dangling_elements.append(element_one_id)
                        dangling_elements_counter += 1
                        break

        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

        self.assertGreaterEqual( len(tracked_ids), max( len(ground_truth) - dangling_elements_counter, len(ground_truth) - 3 ) )

        plt.close('all')
 
    @attr(level = 'standard')
    def test_translation_special_case(self):
        """load first random mesh, move it, and track it.
        """
        sys.setrecursionlimit(40000)
        large_mesh = mesh.load(path.join(dirname(__file__),'data','special_translation_mesh_one.mesh'))

        large_mesh.assign_frame_ids_in_order()
         
        mesh_one = large_mesh.copy_region_of_interest((3,10),(-1,9))
         
        mesh_two = large_mesh.copy_region_of_interest((5,12),(-1,9))
         
        # Build ground truth
        ground_truth_indices = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element.id_in_frame in mesh_two.frame_id_dictionary:
                ground_truth_indices[element_index] = mesh_two.frame_id_dictionary[element.id_in_frame]
                 
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
 
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element_index in ground_truth_indices:
                ground_truth[element.id_in_frame] = mesh_two.elements[ground_truth_indices[element_index]].id_in_frame
         
        print('call tracker')
        tracked_ids = tracking.track(mesh_one, mesh_two)
#         tracked_ids = tracking.find_maximum_common_subgraph(mesh_one, mesh_two)
        print('tracker returned')
  
        mesh_one.plot(path.join(dirname(__file__),'output', 'special_translation_mesh_one.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
        mesh_two.plot(path.join(dirname(__file__),'output', 'special_translation_mesh_two.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
          
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()
  
        self.assertEqual( len(tracked_ids), len(ground_truth) )

#         for frame_id in ground_truth:
#             if mesh_one.get_element_with_frame_id(frame_id).global_id == None:
#                 print 'element with this frame id did not get tracked:'
#                 print frame_id
            
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

    @attr(level = 'standard')
    def test_translation_special_case_one(self):
        """load second random mesh, move it, and track it.
        """
        sys.setrecursionlimit(40000)
        large_mesh = mesh.load(path.join(dirname(__file__),'data','translation_mesh_one.mesh'))

        large_mesh.assign_frame_ids_in_order()
         
        mesh_one = large_mesh.copy_region_of_interest((3,10),(-1,9))
         
        mesh_two = large_mesh.copy_region_of_interest((5,12),(-1,9))
         
        # Build ground truth
        ground_truth_indices = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element.id_in_frame in mesh_two.frame_id_dictionary:
                ground_truth_indices[element_index] = mesh_two.frame_id_dictionary[element.id_in_frame]
                 
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
 
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element_index in ground_truth_indices:
                ground_truth[element.id_in_frame] = mesh_two.elements[ground_truth_indices[element_index]].id_in_frame
         
        print('call tracker')
        tracked_ids = tracking.track(mesh_one, mesh_two)
#         tracked_ids = tracking.find_maximum_common_subgraph(mesh_one, mesh_two)
        print('tracker returned')
  
        mesh_one.plot(path.join(dirname(__file__),'output', 'translation_mesh_two_before.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
        mesh_two.plot(path.join(dirname(__file__),'output', 'translation_mesh_two_after.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
          
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()
  
        self.assertEqual( len(tracked_ids), len(ground_truth) )

        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

    @attr(level = 'standard')
    def test_translation_special_case_two(self):
        """load third random mesh, move it, and track it.
        """
        sys.setrecursionlimit(40000)
        large_mesh = mesh.load(path.join(dirname(__file__),'data','special_translation_mesh_two.mesh'))

        large_mesh.assign_frame_ids_in_order()
         
        mesh_one = large_mesh.copy_region_of_interest((3,10),(-1,9))
         
        mesh_two = large_mesh.copy_region_of_interest((5,12),(-1,9))
         
        # Build ground truth
        ground_truth_indices = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element.id_in_frame in mesh_two.frame_id_dictionary:
                ground_truth_indices[element_index] = mesh_two.frame_id_dictionary[element.id_in_frame]
                 
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
 
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element_index in ground_truth_indices:
                ground_truth[element.id_in_frame] = mesh_two.elements[ground_truth_indices[element_index]].id_in_frame
         
#         tracked_ids = tracking.find_maximum_common_subgraph(mesh_one, mesh_two)
        tracked_ids = tracking.track(mesh_one, mesh_two)

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
#             mesh_one.plot('debug_mesh_special_before_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('debug_mesh_special_after_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#  
        mesh_one.plot(path.join(dirname(__file__),'output', 'special_translation_mesh_two_before.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
        mesh_two.plot(path.join(dirname(__file__),'output', 'special_translation_mesh_two_after.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
          
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()
  
        self.assertEqual( len(tracked_ids), len(ground_truth) - 1 )

        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

    @attr(level = 'standard')
    def test_translation_special_case_three(self):
        """load third random mesh, move it, and track it.
        """
        sys.setrecursionlimit(40000)
        large_mesh = mesh.load(path.join(dirname(__file__),'data','special_translation_mesh_three.mesh'))

        large_mesh.assign_frame_ids_in_order()
         
        mesh_one = large_mesh.copy_region_of_interest((3,10),(-1,9))
         
        mesh_two = large_mesh.copy_region_of_interest((5,12),(-1,9))
         
        # Build ground truth
        ground_truth_indices = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element.id_in_frame in mesh_two.frame_id_dictionary:
                ground_truth_indices[element_index] = mesh_two.frame_id_dictionary[element.id_in_frame]
                 
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_in_reverse_order()
 
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element_index in ground_truth_indices:
                ground_truth[element.id_in_frame] = mesh_two.elements[ground_truth_indices[element_index]].id_in_frame
         
#         tracked_ids = tracking.find_maximum_common_subgraph(mesh_one, mesh_two)
        tracked_ids = tracking.track(mesh_one, mesh_two)

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
#             mesh_one.plot('debug_mesh_special_before_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#             mesh_two.plot('debug_mesh_special_after_' + str(mapping_index) + '.pdf', color_by_global_id = True, 
#                           total_number_of_global_ids = len( tracked_ids ) )
#  
        mesh_one.plot(path.join(dirname(__file__),'output', 'special_translation_mesh_three_before.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
        mesh_two.plot(path.join(dirname(__file__),'output', 'special_translation_mesh_three_after.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
          
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()
  
        self.assertEqual( len(tracked_ids), len(ground_truth) )

        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

    @attr(level = 'standard')
    def test_translation_special_case_six(self):
        """generate a random mesh, move it, and track it.
        """
        sys.setrecursionlimit(40000)
        large_mesh = mesh.load(path.join(dirname(__file__),'data','special_translation_mesh_eight.mesh'))

        large_mesh.assign_frame_ids_in_order()
        
        mesh_one = large_mesh.copy_region_of_interest((3,10),(-1,9))
        
        mesh_two = large_mesh.copy_region_of_interest((5,12),(-1,9))
        
        # Build ground truth
        ground_truth_indices = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element.id_in_frame in mesh_two.frame_id_dictionary:
                ground_truth_indices[element_index] = mesh_two.frame_id_dictionary[element.id_in_frame]
                
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
 
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element_index in ground_truth_indices:
                ground_truth[element.id_in_frame] = mesh_two.elements[ground_truth_indices[element_index]].id_in_frame
        
 
        print('call tracker')
#         tracked_ids = tracking.find_maximum_common_subgraph(mesh_one, mesh_two)
        tracked_ids = tracking.track(mesh_one, mesh_two)
        print('tracker returned')
  
        mesh_one.plot(path.join(dirname(__file__),'output','special_translation_mesh_six_before.pdf'),
                      color_by_global_id = True, total_number_of_global_ids = len(tracked_ids))
        mesh_two.plot(path.join(dirname(__file__),'output','special_translation_mesh_six_after.pdf'),
                      color_by_global_id = True, total_number_of_global_ids = len(tracked_ids))
         
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()

        # Find any dangling elements that may not be identifiable through adjacency
        dangling_elements_counter = 0
        dangling_elements = []
        for element_one_id in ground_truth:
            element_one = mesh_one.get_element_with_frame_id(element_one_id)
            element_two = mesh_two.get_element_with_frame_id(ground_truth[element_one_id])
            ids_adjacent_to_element_one = element_one.get_ids_of_adjacent_elements()
            ids_adjacent_to_element_two = element_two.get_ids_of_adjacent_elements()
            mappable_elements_adjacent_to_elements_one = set(ids_adjacent_to_element_one).intersection(ground_truth.keys())
            mappable_elements_adjacent_to_elements_two = set(ids_adjacent_to_element_two).intersection(ground_truth.values())
            if (len(mappable_elements_adjacent_to_elements_one) == 1 or 
                len(mappable_elements_adjacent_to_elements_two) == 1):
                dangling_elements_counter += 1
                dangling_elements.append(element_one_id)
               
            
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

        self.assertGreaterEqual( len(tracked_ids), max( len(ground_truth) - dangling_elements_counter, len(ground_truth) - 3 ) )

    @attr(level = 'standard')
    def test_translation_special_case_seven(self):
        """generate a random mesh, move it, and track it.
        """
        sys.setrecursionlimit(40000)
        large_mesh = mesh.load(path.join(dirname(__file__),'data','special_translation_mesh_nine.mesh'))

        large_mesh.assign_frame_ids_in_order()
        
        mesh_one = large_mesh.copy_region_of_interest((3,10),(-1,9))
        
        mesh_two = large_mesh.copy_region_of_interest((5,12),(-1,9))
        
        # Build ground truth
        ground_truth_indices = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element.id_in_frame in mesh_two.frame_id_dictionary:
                ground_truth_indices[element_index] = mesh_two.frame_id_dictionary[element.id_in_frame]
                
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
 
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element_index in ground_truth_indices:
                ground_truth[element.id_in_frame] = mesh_two.elements[ground_truth_indices[element_index]].id_in_frame
        
 
        print('call tracker')
#         tracked_ids = tracking.find_maximum_common_subgraph(mesh_one, mesh_two)
        tracked_ids = tracking.track(mesh_one, mesh_two)
        print('tracker returned')
  
        mesh_one.plot(path.join(dirname(__file__),'output','special_translation_mesh_seven_before.pdf'),
                      color_by_global_id = True, total_number_of_global_ids = len(tracked_ids))
        mesh_two.plot(path.join(dirname(__file__),'output','special_translation_mesh_seven_after.pdf'),
                      color_by_global_id = True, total_number_of_global_ids = len(tracked_ids))
         
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()

        # Find any dangling elements that may not be identifiable through adjacency
        dangling_elements_counter = 0
        dangling_elements = []
        for element_one_id in ground_truth:
            element_one = mesh_one.get_element_with_frame_id(element_one_id)
            element_two = mesh_two.get_element_with_frame_id(ground_truth[element_one_id])
            ids_adjacent_to_element_one = element_one.get_ids_of_adjacent_elements()
            ids_adjacent_to_element_two = element_two.get_ids_of_adjacent_elements()
            mappable_elements_adjacent_to_elements_one = set(ids_adjacent_to_element_one).intersection(ground_truth.keys())
            mappable_elements_adjacent_to_elements_two = set(ids_adjacent_to_element_two).intersection(ground_truth.values())
            if (len(mappable_elements_adjacent_to_elements_one) == 1 or 
                len(mappable_elements_adjacent_to_elements_two) == 1):
                dangling_elements_counter += 1
                dangling_elements.append(element_one_id)
            elif (len(mappable_elements_adjacent_to_elements_one) == 2):
                neighbouring_polygon_number_is_the_same = False
                for element in mappable_elements_adjacent_to_elements_one:
                    element_polygon_number = mesh_one.get_element_with_frame_id(element).get_num_nodes()
                    if element_polygon_number == element_one.get_num_nodes():
                        dangling_elements.append(element_one_id)
                        dangling_elements_counter += 1
                        break
            elif (len(mappable_elements_adjacent_to_elements_two) == 2):
                neighbouring_polygon_number_is_the_same = False
                for element in mappable_elements_adjacent_to_elements_two:
                    element_polygon_number = mesh_two.get_element_with_frame_id(element).get_num_nodes()
                    if element_polygon_number == element_two.get_num_nodes():
                        dangling_elements.append(element_one_id)
                        dangling_elements_counter += 1
                        break
               
            
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

        self.assertGreaterEqual( len(tracked_ids), max( len(ground_truth) - dangling_elements_counter, len(ground_truth) - 3 ) )

    @attr(level = 'known_to_fail')
    def test_translation_special_case_eight(self):
        """generate a random mesh, move it, and track it.
           In this example, the local neighbourhood of the first vertex has some symmetry
           that allows it to be mapped incorrectly.
        """
        sys.setrecursionlimit(40000)
        large_mesh = mesh.load(path.join(dirname(__file__),'data','special_translation_mesh_ten.mesh'))

        large_mesh.assign_frame_ids_in_order()
        
        mesh_one = large_mesh.copy_region_of_interest((3,10),(-1,9))
        
        mesh_two = large_mesh.copy_region_of_interest((5,12),(-1,9))
        
        # Build ground truth
        ground_truth_indices = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element.id_in_frame in mesh_two.frame_id_dictionary:
                ground_truth_indices[element_index] = mesh_two.frame_id_dictionary[element.id_in_frame]
                
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
 
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element_index in ground_truth_indices:
                ground_truth[element.id_in_frame] = mesh_two.elements[ground_truth_indices[element_index]].id_in_frame
        
 
        print('call tracker')
#         tracked_ids = tracking.find_maximum_common_subgraph(mesh_one, mesh_two)
        tracked_ids = tracking.track(mesh_one, mesh_two)
        print('tracker returned')
  
        mesh_one.plot(path.join(dirname(__file__),'output','special_translation_mesh_ten_before.pdf'),
                      color_by_global_id = True, total_number_of_global_ids = len(tracked_ids))
        mesh_two.plot(path.join(dirname(__file__),'output','special_translation_mesh_ten_after.pdf'),
                      color_by_global_id = True, total_number_of_global_ids = len(tracked_ids))
         
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()

        # Find any dangling elements that may not be identifiable through adjacency
        dangling_elements_counter = 0
        dangling_elements = []
        for element_one_id in ground_truth:
            element_one = mesh_one.get_element_with_frame_id(element_one_id)
            element_two = mesh_two.get_element_with_frame_id(ground_truth[element_one_id])
            ids_adjacent_to_element_one = element_one.get_ids_of_adjacent_elements()
            ids_adjacent_to_element_two = element_two.get_ids_of_adjacent_elements()
            mappable_elements_adjacent_to_elements_one = set(ids_adjacent_to_element_one).intersection(ground_truth.keys())
            mappable_elements_adjacent_to_elements_two = set(ids_adjacent_to_element_two).intersection(ground_truth.values())
            if (len(mappable_elements_adjacent_to_elements_one) == 1 or 
                len(mappable_elements_adjacent_to_elements_two) == 1):
                dangling_elements_counter += 1
                dangling_elements.append(element_one_id)
            elif (len(mappable_elements_adjacent_to_elements_one) == 2):
                neighbouring_polygon_number_is_the_same = False
                for element in mappable_elements_adjacent_to_elements_one:
                    element_polygon_number = mesh_one.get_element_with_frame_id(element).get_num_nodes()
                    if element_polygon_number == element_one.get_num_nodes():
                        dangling_elements.append(element_one_id)
                        dangling_elements_counter += 1
                        break
            elif (len(mappable_elements_adjacent_to_elements_two) == 2):
                neighbouring_polygon_number_is_the_same = False
                for element in mappable_elements_adjacent_to_elements_two:
                    element_polygon_number = mesh_two.get_element_with_frame_id(element).get_num_nodes()
                    if element_polygon_number == element_two.get_num_nodes():
                        dangling_elements.append(element_one_id)
                        dangling_elements_counter += 1
                        break
               
            
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

        self.assertGreaterEqual( len(tracked_ids), max( len(ground_truth) - dangling_elements_counter, len(ground_truth) - 3 ) )

    @attr(level = 'standard')
    def test_translation_special_case_nine(self):
        """generate a random mesh, move it, and track it.
        """
        sys.setrecursionlimit(40000)
        large_mesh = mesh.load(path.join(dirname(__file__),'data','special_translation_mesh_eleven.mesh'))

        large_mesh.assign_frame_ids_in_order()
        
        mesh_one = large_mesh.copy_region_of_interest((3,10),(-1,9))
        
        mesh_two = large_mesh.copy_region_of_interest((5,12),(-1,9))
        
        # Build ground truth
        ground_truth_indices = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element.id_in_frame in mesh_two.frame_id_dictionary:
                ground_truth_indices[element_index] = mesh_two.frame_id_dictionary[element.id_in_frame]
                
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
 
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element_index in ground_truth_indices:
                ground_truth[element.id_in_frame] = mesh_two.elements[ground_truth_indices[element_index]].id_in_frame
        
 
        print('call tracker')
        tracked_ids = tracking.track(mesh_one, mesh_two)
        print('tracker returned')
  
        mesh_one.plot(path.join(dirname(__file__),'output','special_translation_mesh_eleven_before.pdf'),
                      color_by_global_id = True, total_number_of_global_ids = len(tracked_ids))
        mesh_two.plot(path.join(dirname(__file__),'output','special_translation_mesh_eleven_after.pdf'),
                      color_by_global_id = True, total_number_of_global_ids = len(tracked_ids))
         
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()

        # Find any dangling elements that may not be identifiable through adjacency
        dangling_elements_counter = 0
        dangling_elements = []
        for element_one_id in ground_truth:
            element_one = mesh_one.get_element_with_frame_id(element_one_id)
            element_two = mesh_two.get_element_with_frame_id(ground_truth[element_one_id])
            ids_adjacent_to_element_one = element_one.get_ids_of_adjacent_elements()
            ids_adjacent_to_element_two = element_two.get_ids_of_adjacent_elements()
            mappable_elements_adjacent_to_elements_one = set(ids_adjacent_to_element_one).intersection(ground_truth.keys())
            mappable_elements_adjacent_to_elements_two = set(ids_adjacent_to_element_two).intersection(ground_truth.values())
            if (len(mappable_elements_adjacent_to_elements_one) == 1 or 
                len(mappable_elements_adjacent_to_elements_two) == 1):
                dangling_elements_counter += 1
                dangling_elements.append(element_one_id)
            elif (len(mappable_elements_adjacent_to_elements_one) == 2):
                neighbouring_polygon_number_is_the_same = False
                for element in mappable_elements_adjacent_to_elements_one:
                    element_polygon_number = mesh_one.get_element_with_frame_id(element).get_num_nodes()
                    if element_polygon_number == element_one.get_num_nodes():
                        dangling_elements.append(element_one_id)
                        dangling_elements_counter += 1
                        break
            elif (len(mappable_elements_adjacent_to_elements_two) == 2):
                neighbouring_polygon_number_is_the_same = False
                for element in mappable_elements_adjacent_to_elements_two:
                    element_polygon_number = mesh_two.get_element_with_frame_id(element).get_num_nodes()
                    if element_polygon_number == element_two.get_num_nodes():
                        dangling_elements.append(element_one_id)
                        dangling_elements_counter += 1
                        break
               
            
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

        self.assertGreaterEqual( len(tracked_ids), max( len(ground_truth) - dangling_elements_counter, len(ground_truth) - 3 ) )

    @attr(level = 'standard')
    def test_translation_special_case_ten(self):
        """generate a random mesh, move it, and track it.
        """
        sys.setrecursionlimit(40000)
        large_mesh = mesh.load(path.join(dirname(__file__),'data','special_translation_mesh_twelve.mesh'))

        large_mesh.assign_frame_ids_in_order()
        
        mesh_one = large_mesh.copy_region_of_interest((3,10),(-1,9))
        
        mesh_two = large_mesh.copy_region_of_interest((5,12),(-1,9))
        
        # Build ground truth
        ground_truth_indices = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element.id_in_frame in mesh_two.frame_id_dictionary:
                ground_truth_indices[element_index] = mesh_two.frame_id_dictionary[element.id_in_frame]
                
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
 
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element_index in ground_truth_indices:
                ground_truth[element.id_in_frame] = mesh_two.elements[ground_truth_indices[element_index]].id_in_frame
        
 
        print('call tracker')
        tracked_ids = tracking.track(mesh_one, mesh_two)
        print('tracker returned')
  
        mesh_one.plot(path.join(dirname(__file__),'output','special_translation_mesh_twelve_before.pdf'),
                      color_by_global_id = True, total_number_of_global_ids = len(tracked_ids))
        mesh_two.plot(path.join(dirname(__file__),'output','special_translation_mesh_twelve_after.pdf'),
                      color_by_global_id = True, total_number_of_global_ids = len(tracked_ids))
         
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()

        # Find any dangling elements that may not be identifiable through adjacency
        dangling_elements_counter = 0
        dangling_elements = []
        for element_one_id in ground_truth:
            element_one = mesh_one.get_element_with_frame_id(element_one_id)
            element_two = mesh_two.get_element_with_frame_id(ground_truth[element_one_id])
            ids_adjacent_to_element_one = element_one.get_ids_of_adjacent_elements()
            ids_adjacent_to_element_two = element_two.get_ids_of_adjacent_elements()
            mappable_elements_adjacent_to_elements_one = set(ids_adjacent_to_element_one).intersection(ground_truth.keys())
            mappable_elements_adjacent_to_elements_two = set(ids_adjacent_to_element_two).intersection(ground_truth.values())
            if (len(mappable_elements_adjacent_to_elements_one) == 1 or 
                len(mappable_elements_adjacent_to_elements_two) == 1):
                dangling_elements_counter += 1
                dangling_elements.append(element_one_id)
            elif (len(mappable_elements_adjacent_to_elements_one) == 2):
                neighbouring_polygon_number_is_the_same = False
                for element in mappable_elements_adjacent_to_elements_one:
                    element_polygon_number = mesh_one.get_element_with_frame_id(element).get_num_nodes()
                    if element_polygon_number == element_one.get_num_nodes():
                        dangling_elements.append(element_one_id)
                        dangling_elements_counter += 1
                        break
            elif (len(mappable_elements_adjacent_to_elements_two) == 2):
                neighbouring_polygon_number_is_the_same = False
                for element in mappable_elements_adjacent_to_elements_two:
                    element_polygon_number = mesh_two.get_element_with_frame_id(element).get_num_nodes()
                    if element_polygon_number == element_two.get_num_nodes():
                        dangling_elements.append(element_one_id)
                        dangling_elements_counter += 1
                        break
               
            
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

        self.assertGreaterEqual( len(tracked_ids), max( len(ground_truth) - dangling_elements_counter, len(ground_truth) - 3 ) )

    @attr(level = 'standard')
    def test_translation_special_case_eleven(self):
        """generate a random mesh, move it, and track it.
        """
        sys.setrecursionlimit(40000)
        large_mesh = mesh.load(path.join(dirname(__file__),'data','special_translation_mesh_thirteen.mesh'))

        large_mesh.assign_frame_ids_in_order()
        
        mesh_one = large_mesh.copy_region_of_interest((3,10),(-1,9))
        
        mesh_two = large_mesh.copy_region_of_interest((5,12),(-1,9))
        
        # Build ground truth
        ground_truth_indices = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element.id_in_frame in mesh_two.frame_id_dictionary:
                ground_truth_indices[element_index] = mesh_two.frame_id_dictionary[element.id_in_frame]
                
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
 
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            if element_index in ground_truth_indices:
                ground_truth[element.id_in_frame] = mesh_two.elements[ground_truth_indices[element_index]].id_in_frame
        
 
        print('call tracker')
#         tracked_ids = tracking.find_maximum_common_subgraph(mesh_one, mesh_two)
        tracked_ids = tracking.track(mesh_one, mesh_two)
        print('tracker returned')
  
        mesh_one.plot(path.join(dirname(__file__),'output','special_translation_mesh_thirteen_before.pdf'),
                      color_by_global_id = True, total_number_of_global_ids = len(tracked_ids))
        mesh_two.plot(path.join(dirname(__file__),'output','special_translation_mesh_thirteen_after.pdf'),
                      color_by_global_id = True, total_number_of_global_ids = len(tracked_ids))
         
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()

        # Find any dangling elements that may not be identifiable through adjacency
        dangling_elements_counter = 0
        dangling_elements = []
        for element_one_id in ground_truth:
            element_one = mesh_one.get_element_with_frame_id(element_one_id)
            element_two = mesh_two.get_element_with_frame_id(ground_truth[element_one_id])
            ids_adjacent_to_element_one = element_one.get_ids_of_adjacent_elements()
            ids_adjacent_to_element_two = element_two.get_ids_of_adjacent_elements()
            mappable_elements_adjacent_to_elements_one = set(ids_adjacent_to_element_one).intersection(ground_truth.keys())
            mappable_elements_adjacent_to_elements_two = set(ids_adjacent_to_element_two).intersection(ground_truth.values())
            if (len(mappable_elements_adjacent_to_elements_one) == 1 or 
                len(mappable_elements_adjacent_to_elements_two) == 1):
                dangling_elements_counter += 1
                dangling_elements.append(element_one_id)
            elif (len(mappable_elements_adjacent_to_elements_one) == 2):
                neighbouring_polygon_number_is_the_same = False
                for element in mappable_elements_adjacent_to_elements_one:
                    element_polygon_number = mesh_one.get_element_with_frame_id(element).get_num_nodes()
                    if element_polygon_number == element_one.get_num_nodes():
                        dangling_elements.append(element_one_id)
                        dangling_elements_counter += 1
                        break
            elif (len(mappable_elements_adjacent_to_elements_two) == 2):
                neighbouring_polygon_number_is_the_same = False
                for element in mappable_elements_adjacent_to_elements_two:
                    element_polygon_number = mesh_two.get_element_with_frame_id(element).get_num_nodes()
                    if element_polygon_number == element_two.get_num_nodes():
                        dangling_elements.append(element_one_id)
                        dangling_elements_counter += 1
                        break
               
            
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
 
        for global_id in tracked_ids:
            # and that the mapping coincides with the ground truth for all tracked ids
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEquals( ground_truth[element_one.id_in_frame], element_two.id_in_frame )

        self.assertGreaterEqual( len(tracked_ids), max( len(ground_truth) - dangling_elements_counter, len(ground_truth) - 3 ) )


    @attr(level = 'weekly')
    def test_track_translation_multiple_times(self):
        for n in range(100):
            print('translation run ' + str(n))
            self.xest_translation()
