# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/LICENSE.

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
    def test_maximum_common_subgraph_finder_on_simple_mesh(self):
        mesh_one = mesh.creation.generate_hexagonal(2,2)

        mesh_two = mesh.creation.generate_hexagonal(2,2)

        mesh_one.assign_frame_ids_in_order()
        
        # assign ids to second mesh in reverse order
        for index, element in enumerate(mesh_two.elements):
            element.id_in_frame = 3 - index
        
        subgraph_finder = tracking.SlowMaximumCommonSubgraphFinder(mesh_one, mesh_two)

        assert(subgraph_finder.extendable({0:3, 1: 2, 2: 1})) 
        assert(not subgraph_finder.extendable({0:3, 1: 2, 2: 1, 3 :0})) 
        
        assert(subgraph_finder.pick_next_vertex({0:3, 1: 2}) in set([2,3]))
        
        assert(subgraph_finder.get_image_of_edges_from_vertex_in_network_one(0, {}) == set())
        assert(subgraph_finder.get_edges_from_vertex_in_network_two(3, {}) == set())
        assert(subgraph_finder.get_mappable_vertices(0, {}) == [3])

    @attr(level = 'known_failing')
    def xest_track_shuffled_2_by_2_mesh(self):
        """ This tests whether we can identify cells after they got shuffled (2 by 2)
            This test is known to fail: due to the symmetry of the problem cells tend to
            be mapped in the wrong order.
        """
        mesh_one = mesh.creation.generate_hexagonal(2,2)
        mesh_two = mesh.creation.generate_hexagonal(2,2)
          
        mesh_one.assign_frame_ids_in_order()
          
        mesh_two.assign_frame_ids_in_reverse_order()
  
        id_map = tracking.track(mesh_one, mesh_two)
          
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()
  
        subgraph_network = tracking.generate_subgraph_network(mesh_two, id_map)
  
        self.assertEqual(len(id_map), 4)
  
        for index_one in id_map:
            element_one = mesh_one.get_element_with_frame_id(index_one)
            element_two = mesh_two.get_element_with_global_id(id_map[index_one])
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertEqual(network_one.degree(index_one), network_two.degree(id_map[index_one]))
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
            np.testing.assert_almost_equal(element_one.calculate_centroid(), element_two.calculate_centroid())
  
    @attr(level = 'standard')
    def xest_track_shuffled_3_by_3_mesh(self):
        """ This tests whether we can identify cells after they got shuffled (3 by 3)
        """
        mesh_one = mesh.creation.generate_hexagonal(3,3)
        mesh_two = mesh.creation.generate_hexagonal(3,3)
   
        mesh_one.assign_frame_ids_in_order()
        mesh_one.plot(path.join(dirname(__file__),'output','mesh_one.pdf'))
          
        mesh_two.assign_frame_ids_randomly()
  
        mesh_two.plot(path.join(dirname(__file__),'output','mesh_two_before_tracking.pdf'))
  
        tracked_ids = tracking.track(mesh_one, mesh_two)
      
        mesh_two.plot(path.join(dirname(__file__),'output','mesh_two_after_tracking.pdf'), color_by_global_id = True, total_number_of_global_ids = 9)
  
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()
  
        subgraph_network = tracking.generate_subgraph_network(mesh_two, tracked_ids)
  
        self.assertEqual(len(tracked_ids), 9)
  
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_frame_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertEqual(network_one.degree(element_one.id_in_frame), network_two.degree(element_two.id_in_frame))
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
            np.testing.assert_almost_equal(element_one.calculate_centroid(), element_two.calculate_centroid())
  
    @attr(level = 'standard')
    def xest_track_shuffled_6_by_6_mesh(self):
        """ This tests whether we can identify cells after they got shuffled (6 by 6)
        """
        mesh_one = mesh.creation.generate_hexagonal(6,6)
        mesh_two = mesh.creation.generate_hexagonal(6,6)
   
        mesh_one.assign_frame_ids_in_order()
        mesh_one.plot(path.join(dirname(__file__),'output','mesh_one_hex_six.pdf'))
          
        mesh_two.assign_frame_ids_randomly()
  
        mesh_two.plot(path.join(dirname(__file__),'output','mesh_two_hex_six_before_tracking.pdf'))
  
        tracked_ids = tracking.track(mesh_one, mesh_two)
      
        mesh_two.plot(path.join(dirname(__file__),'output','mesh_two_hex_six_after_tracking.pdf'), color_by_global_id = True, total_number_of_global_ids = 36)
  
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()
  
        subgraph_network = tracking.generate_subgraph_network(mesh_two, tracked_ids)
  
        self.assertEqual(len(tracked_ids), 36)
  
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_frame_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertEqual(network_one.degree(element_one.id_in_frame), network_two.degree(element_two.id_in_frame))
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
            np.testing.assert_almost_equal(element_one.calculate_centroid(), element_two.calculate_centroid())
 
    @attr(level = 'standard')
    def test_track_3_by_3_random_mesh(self):
        """generate a small random mesh and track it.
        """
          
        mesh_one = mesh.creation.generate_random_tesselation(3,3,7)
        mesh_two = copy.deepcopy(mesh_one)
          
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
          
        mesh_one.plot(path.join(dirname(__file__),'output','random_mesh_one.pdf'))
        mesh_two.plot(path.join(dirname(__file__),'output','random_mesh_two_before_tracking.pdf'))
          
        tracked_ids = tracking.track(mesh_one, mesh_two)
  
        mesh_two.plot(path.join(dirname(__file__),'output','random_mesh_after_tracking.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = mesh_one.get_num_elements())
          
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()
  
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertEqual(network_one.degree(element_one.id_in_frame), network_two.degree(element_two.id_in_frame))
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
            np.testing.assert_almost_equal(element_one.calculate_centroid(), element_two.calculate_centroid())       
      
    @attr(level = 'standard')
    def xest_track_9_by_9_random_mesh(self):
        """generate a small random mesh and track it.
        """
        sys.setrecursionlimit(40000)
        mesh_one = mesh.creation.generate_random_tesselation(9,9)
        mesh_one.save(path.join(dirname(__file__),'output','nine_by_nine_permutation_mesh_one.mesh'))
        mesh_two = copy.deepcopy(mesh_one)
         
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
         
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame
 
        mesh_one.plot(path.join(dirname(__file__),'output','larger_random_mesh_one.pdf'))
        mesh_two.plot(path.join(dirname(__file__),'output','larger_random_mesh_two_before_tracking.pdf'))
         
        tracked_ids = tracking.track(mesh_one, mesh_two)

        dangling_elements_one = mesh_one.get_dangling_element_ids(ground_truth.keys()) 
        self.assertGreaterEqual(len(tracked_ids), 
                         max(mesh_one.get_num_elements() - 3, 
                             mesh_one.get_num_elements() - len(dangling_elements_one)))
 
        mesh_one.plot(path.join(dirname(__file__),'output','nine_by_nine_permuation_mesh_before.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = mesh_one.get_num_elements())
        mesh_two.plot(path.join(dirname(__file__),'output','nine_by_nine_permuation_mesh_after.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = mesh_one.get_num_elements())
         
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()
 
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertEqual(network_one.degree(element_one.id_in_frame), network_two.degree(element_two.id_in_frame))
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
            np.testing.assert_almost_equal(element_one.calculate_centroid(), element_two.calculate_centroid())       

        plt.close('all')
     
    @attr(level = 'standard')
    def test_track_special_permutation_mesh(self):
        """generate a small random mesh and track it.
        """
         
        mesh_one = mesh.load(path.join(dirname(__file__),'data','special_permutation_mesh.mesh'))
        mesh_two = copy.deepcopy(mesh_one)

        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
 
        # build ground truth for testing the mapping
        ground_truth = {}
        for element_index, element in enumerate(mesh_one.elements):
            ground_truth[element.id_in_frame] = mesh_two.elements[element_index].id_in_frame
 
        mesh_one.plot(path.join(dirname(__file__),'output','special_permutation_mesh_one.pdf'))
        mesh_two.plot(path.join(dirname(__file__),'output','special_permutation_mesh_two_before_tracking.pdf'))
         
        tracked_ids = tracking.track(mesh_one, mesh_two)

        mesh_one.plot(path.join(dirname(__file__),'output','special_permutation_mesh_before.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = mesh_one.get_num_elements())
        mesh_two.plot(path.join(dirname(__file__),'output','special_permutation_mesh_after.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = mesh_one.get_num_elements())

        dangling_elements_one = mesh_one.get_dangling_element_ids(ground_truth.keys()) 
        self.assertGreaterEqual(len(tracked_ids), 
                         max(mesh_one.get_num_elements() - 3, 
                             mesh_one.get_num_elements() - len(dangling_elements_one)))
         
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()
 
        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_global_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertEqual(network_one.degree(element_one.id_in_frame), network_two.degree(element_two.id_in_frame))
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
            np.testing.assert_almost_equal(element_one.calculate_centroid(), element_two.calculate_centroid())       
 
    @attr(level = 'weekly')
    def test_track_very_large_random_mesh(self):
        """generate a small random mesh and track it.
        """
        sys.setrecursionlimit(40000)
        mesh_one = mesh.creation.generate_random_tesselation(30,30)
        
        print 'start to copy'
        mesh_two = copy.deepcopy(mesh_one)
        print 'stop to copy'
         
        mesh_one.assign_frame_ids_in_order()
        mesh_two.assign_frame_ids_randomly()
         
        mesh_one.plot(path.join(dirname(__file__),'output','very_large_random_mesh_one.pdf'))
        mesh_two.plot(path.join(dirname(__file__),'output','very_large_random_mesh_two_before_tracking.pdf'))
         
        tracked_ids = tracking.track(mesh_one, mesh_two)
 
        mesh_two.plot(path.join(dirname(__file__),'output','very_large_random_mesh_after_tracking.pdf'), color_by_global_id = True, 
                      total_number_of_global_ids = mesh_one.get_num_elements())
         
        network_one = mesh_one.generate_network()
        network_two = mesh_two.generate_network()
 
        self.assertEqual( len(tracked_ids), mesh_one.get_num_elements())

        for global_id in tracked_ids:
            element_one = mesh_one.get_element_with_frame_id(global_id)
            element_two = mesh_two.get_element_with_global_id(global_id)
            self.assertEqual(element_one.get_num_nodes(), element_two.get_num_nodes())
            self.assertEqual(network_one.degree(element_one.id_in_frame), network_two.degree(element_two.id_in_frame))
            self.assertAlmostEqual(element_one.calculate_area(), element_two.calculate_area())
            np.testing.assert_almost_equal(element_one.calculate_centroid(), element_two.calculate_centroid())       
     
    @attr(level = 'weekly')
    def test_track_permutation_multiple_times(self):
        for n in range(100):
            print 'permutation run ' + str(n)
            self.xest_track_9_by_9_random_mesh()
