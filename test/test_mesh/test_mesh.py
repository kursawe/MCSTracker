# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This test tests our initial mesh creation framework.
"""
import unittest
import mesh
import numpy as np
import networkx as nx
import os
import copy
import sys
from os import path
from os.path import dirname
from nose.plugins.attrib import attr

@attr(level = 'standard')
class TestMesh(unittest.TestCase):
 
    def test_element(self):
        """
        Test the creation and methods of a single element
        """
        nodes = [mesh.Node([0.0,0.0]), mesh.Node([1.0, 0.0]), mesh.Node([1.0,1.0]), mesh.Node([0.0,1.0])]
        for node in nodes:
            self.assertEquals(node.adjacent_elements,[])
        element  = mesh.Element(nodes, 42)
        for node in element.nodes:
            self.assertEqual(len(node.adjacent_elements),1)
            self.assertEqual(node.adjacent_elements[0].id_in_frame, 42)
            self.assertEqual(node.get_adjacent_element_ids(), [42])
        self.assertEqual(element.get_num_nodes(),4)
        self.assertEqual(element.id_in_frame, 42)
        self.assertAlmostEqual(element.calculate_area(), 1.0)
        np.testing.assert_almost_equal(element.calculate_centroid(), [0.5,0.5])

    def test_mesh_methods(self):
        """create a 2 by 2 rectangular mesh and test the number
        of elements, their number of nodes, etc"""
        # make some vertices
        vertices = []
        vertices.append(mesh.Node([0.0, 0.0]))
        vertices.append(mesh.Node([1.0, 0.0]))
        vertices.append(mesh.Node([2.0, 0.0]))
        vertices.append(mesh.Node([0.0, 1.0]))
        vertices.append(mesh.Node([1.0, 1.0]))
        vertices.append(mesh.Node([2.0, 1.0]))
        vertices.append(mesh.Node([0.0, 2.0]))
        vertices.append(mesh.Node([1.0, 2.0]))
        vertices.append(mesh.Node([2.0, 2.0]))
        
        # make some elements
        list_of_elements = []
        list_of_elements.append(mesh.Element([vertices[0],
                                              vertices[1],
                                              vertices[4],
                                              vertices[3]]))
        list_of_elements.append(mesh.Element([vertices[1],
                                              vertices[2],
                                              vertices[5],
                                              vertices[4]]))
        list_of_elements.append(mesh.Element([vertices[3],
                                              vertices[4],
                                              vertices[7],
                                              vertices[6]]))
        list_of_elements.append(mesh.Element([vertices[4],
                                              vertices[5],
                                              vertices[8],
                                              vertices[7]]))
        
        simple_mesh = mesh.Mesh(vertices, list_of_elements)
        simple_mesh.assign_frame_ids_in_order()

        self.assertEqual(simple_mesh.get_num_nodes(), 9)
        self.assertEqual(simple_mesh.get_num_elements(), 4)
        for counter, element in enumerate( simple_mesh.elements ):
            self.assertEqual(element.get_num_nodes(), 4)
            self.assertEqual(element.id_in_frame, counter)
            
        for counter, node in enumerate(simple_mesh.nodes):
            self.assertEqual(node.id, counter)
            
        expected_edges = np.array([[0,1],
                                   [1,4],
                                   [4,3],
                                   [3,0],
                                   [1,2],
                                   [2,5],
                                   [5,4],
                                   [4,7],
                                   [7,6],
                                   [6,3],
                                   [5,8],
                                   [8,7]])

        edges = simple_mesh.collect_edges()
        np.testing.assert_equal(edges, expected_edges)
        
        inner_edges = simple_mesh.collect_elements_of_inner_edges()
        self.assertEqual(len(inner_edges), 4)
        assert( set([0,1]) in inner_edges )
        assert( set([0,2]) in inner_edges )
        assert( set([1,3]) in inner_edges )
        assert( set([2,3]) in inner_edges )

        nodes_between_element_0_and_1 = simple_mesh.get_nodes_shared_by_elements([0,1])
        nodes_between_element_0_and_2 = simple_mesh.get_nodes_shared_by_elements([0,2])
        nodes_between_element_0_and_3 = simple_mesh.get_nodes_shared_by_elements([0,3])
        self.assertEqual( set([1,4]), nodes_between_element_0_and_1 )
        self.assertEqual( set([3,4]), nodes_between_element_0_and_2 )
        self.assertEqual( set([4]), nodes_between_element_0_and_3 )
        
        simple_mesh.assign_frame_ids_in_order()
        network_representation = simple_mesh.generate_network()
        
        self.assertEqual(network_representation.number_of_nodes(), 4)
        self.assertEqual(network_representation.number_of_edges(), 4)
        
        assert(network_representation.has_edge(0,1))
        assert(network_representation.has_edge(0,2))
        assert(network_representation.has_edge(1,3))
        assert(network_representation.has_edge(2,3))
        
        # test width and height methods
        
        self.assertAlmostEquals(simple_mesh.calculate_width(), 2.0)
        self.assertAlmostEquals(simple_mesh.calculate_height(), 2.0)

        filename = 'simple_mesh.pdf'
        if os.path.exists(filename):
            os.remove(filename)
        simple_mesh.plot(path.join(dirname(__file__),'output',filename))
        
        assert(os.path.exists(path.join(dirname(__file__),'output',filename)))

        # test the calculate total area and calculate average area methods
        # as well as calculate centroid       

        self.assertEqual(simple_mesh.calculate_total_area(), 4.0)
        self.assertEqual(simple_mesh.calculate_average_element_area(), 1.0 )
        np.testing.assert_almost_equal(simple_mesh.calculate_centre(), [1.0, 1.0] )

    def test_assign_ids_to_elements_and_generate_network(self):
        """Test id_in_frame assignment functions"""
        my_mesh = mesh.creation.generate_hexagonal(5,5)
        my_mesh.assign_frame_ids_in_order()
        for counter, element in enumerate(my_mesh.elements):
            self.assertEquals(element.id_in_frame, counter)
        
        for id in range(25):
            self.assertEquals(my_mesh.get_element_with_frame_id(id).id_in_frame, id)
            
        my_mesh.assign_frame_ids_in_reverse_order()
        for counter, element in enumerate(my_mesh.elements):
            self.assertEquals(element.id_in_frame, 24 - counter)
        
        for id in range(25):
            self.assertEquals(my_mesh.get_element_with_frame_id(id).id_in_frame, id)
            
        # for completeness: test area calculation
        
        self.assertEqual(my_mesh.calculate_total_area(), 25.0)
        self.assertEqual(my_mesh.calculate_average_element_area(), 1.0)

        my_mesh.assign_frame_ids_randomly()
        all_ids = np.zeros(my_mesh.get_num_elements())
        for counter, element in enumerate(my_mesh.elements):
            all_ids[counter] = element.id_in_frame
        self.assertEquals(all_ids.max(), my_mesh.get_num_elements()-1)
        self.assertEquals(np.unique(all_ids).size, my_mesh.get_num_elements())
        
        for id in range(25):
            self.assertEquals(my_mesh.get_element_with_frame_id(id).id_in_frame, id)

        network_representation = my_mesh.generate_network()
        
        self.assertEquals(network_representation.number_of_nodes(), 25)
        # the number of connections between cells in this mesh will be 56
        self.assertEquals(network_representation.number_of_edges(), 56)
        
        for element in my_mesh.elements:
            np.testing.assert_almost_equal(network_representation.node[element.id_in_frame]['position'],
                                           element.calculate_centroid())
            self.assertEquals(network_representation.node[element.id_in_frame]['num_neighbours'], 6)
    
    def test_assign_global_ids(self):
        """test whether we can assign global ids correctly"""
        
        my_mesh = mesh.creation.generate_hexagonal(5,5)
        
        for counter, element in enumerate(my_mesh.elements):
            element.global_id = counter
            
        my_mesh.index_global_ids()
        
        for counter in range(my_mesh.get_num_elements()):
            self.assertEquals(my_mesh.get_element_with_global_id(counter).global_id, counter)
            
    def test_copy_mesh(self):
        """test whether we can use deepcopy() to copy a mesh"""
        
        random_mesh = mesh.creation.generate_random_tesselation(4, 5, 0)
        random_mesh.assign_frame_ids_in_order()

        mesh_copy = copy.deepcopy(random_mesh)
        
        node_5_position = random_mesh.nodes[5].position
        self.assertNotAlmostEqual(node_5_position[0], 1.01, msg = 'randomly ran into this, just run the test again!')
        self.assertNotAlmostEqual(node_5_position[1], 1.01, msg = 'randomly ran into this, just run the test again!')
        
        mesh_copy.nodes[5].position = np.array([1.01, 1.01])

        self.assertNotAlmostEqual(random_mesh.nodes[5].position[0], 1.01)
        self.assertNotAlmostEqual(random_mesh.nodes[5].position[1], 1.01)
        
        frame_id_of_7th_element = random_mesh.elements[7].id_in_frame
        
        mesh_copy.elements[7].id_in_frame = 485
        
        self.assertEqual(random_mesh.elements[7].id_in_frame, frame_id_of_7th_element)
        self.assertNotEqual(random_mesh.elements[7].id_in_frame, 485)
    
    def test_copy_region_of_interest(self):
        """test whether we can cut a smaller mesh from a bigger mesh"""
        
        large_mesh = mesh.creation.generate_random_tesselation(20,20)
        large_mesh.assign_frame_ids_randomly()
        sys.setrecursionlimit(40000)
        large_mesh_copy = copy.deepcopy(large_mesh)
        
        small_mesh = large_mesh.copy_region_of_interest((3,12), (5,15))
        
        small_mesh.plot(path.join(dirname(__file__),'output','cut_mesh.pdf'))
        
        small_mesh_width = small_mesh.calculate_width()
        small_mesh_height = small_mesh.calculate_height()
        
        self.assertGreater(small_mesh_width, 8)
        self.assertLess(small_mesh_width, 11)
        
        self.assertGreater(small_mesh_height, 9)
        self.assertLess(small_mesh_height, 12)
        
        ## Test that the copied nodes occupy a window of the correct size
        for node in small_mesh.nodes:
            self.assertGreaterEqual(node.position[0], 0.0)
            self.assertLess(node.position[0], 11)
            self.assertGreaterEqual(node.position[1], 0.0)
            self.assertLess(node.position[1], 11.5)

        # Test that the original mesh didn't get altered
        for element_index, copied_element in enumerate(large_mesh_copy.elements):
            original_element = large_mesh.elements[element_index]
            np.testing.assert_almost_equal(original_element.calculate_centroid(), copied_element.calculate_centroid())
            self.assertEqual(original_element.id_in_frame, copied_element.id_in_frame)
            
        # Test that the new elements have the same ids as their counterparts in the old mesh,
        # and that they are properly indexed
        for element_index, element in enumerate(small_mesh.elements):
            self.assert_(element.id_in_frame != None)
            # get element with same frame id from large mesh
            other_element = large_mesh.get_element_with_frame_id(element.id_in_frame)
            self.assertAlmostEqual(element.calculate_area(), other_element.calculate_area())
            self.assertAlmostEqual(element.get_num_nodes(), other_element.get_num_nodes())
            self.assert_( small_mesh.frame_id_dictionary[element.id_in_frame] == element_index )

    def test_check_mesh_equality(self):
        """tests whether we can find whether meshes are identical or not"""
        
        mesh_one = mesh.creation.generate_random_tesselation(5,5)
        mesh_two = copy.deepcopy(mesh_one)
        
        assert( mesh_one == mesh_two )
        assert( not (mesh_one != mesh_two) )
        
        mesh_one.elements[0].id_in_frame = 1000
        self.assertEqual( mesh_one.elements[0].id_in_frame, 1000)
        self.assertEqual( mesh_two.elements[0].id_in_frame, None)
        
        assert( mesh_one != mesh_two )
        assert( not (mesh_one == mesh_two) )
        
        mesh_one.elements[0].id_in_frame = None
        assert( mesh_one == mesh_two )
        
        mesh_one.nodes[0].id = 1000
        assert( mesh_one != mesh_two )
        assert( not (mesh_one == mesh_two) )
        
        mesh_one.nodes[0].id = mesh_two.nodes[0].id
        assert( mesh_one == mesh_two )

        mesh_one.nodes[0].position = [1000.0, 1000.0]
        assert( mesh_one != mesh_two )
        assert( not (mesh_one == mesh_two) )
        
        mesh_one.nodes[0].position = np.array( mesh_two.nodes[0].position )
        assert( mesh_one == mesh_two )

        mesh_one.nodes[0].adjacent_elements.append( mesh_one.elements[-1] )
        assert( mesh_one != mesh_two )
        assert( not (mesh_one == mesh_two) )
        
    def test_get_most_central_element_and_node(self):
        """tests whether we can get the element whose centroid is closest to the centre of the mesh, 
        as well as the most central node."""
        
        sample_mesh = mesh.creation.generate_random_tesselation(10,10)
        
        mesh_centre = sample_mesh.calculate_centre()
        
        central_element = sample_mesh.find_most_central_element()
        
        central_element_distance_to_centre = np.linalg.norm(central_element.calculate_centroid() - mesh_centre)
        
        for element in sample_mesh.elements:
            this_distance = np.linalg.norm(element.calculate_centroid() - mesh_centre)
            self.assertGreaterEqual(this_distance, central_element_distance_to_centre)
        
        central_node = sample_mesh.find_most_central_node()
        
        central_node_distance_to_centre = np.linalg.norm(central_node.position - mesh_centre)
        
        for node in sample_mesh.nodes:
            this_distance = np.linalg.norm(node.position - mesh_centre)
            self.assertGreaterEqual(this_distance, central_node_distance_to_centre)

    def test_remove_boundary_elements(self):
        """tests whether we can correctly remove all boundary elements in a mesh
        """
        
        vertices = []
        vertices.append(mesh.Node([-2.0, 0.0]))
        vertices.append(mesh.Node([ 0.0, 0.0]))
        vertices.append(mesh.Node([ 2.0, 0.0]))
        vertices.append(mesh.Node([-2.0, 2.0]))
        vertices.append(mesh.Node([-1.0, 2.0]))
        vertices.append(mesh.Node([ 0.0, 1.0]))
        vertices.append(mesh.Node([ 1.0, 2.0]))
        vertices.append(mesh.Node([ 2.0, 2.0]))
        vertices.append(mesh.Node([-2.0, 3.0]))
        vertices.append(mesh.Node([-1.0, 3.0]))
        vertices.append(mesh.Node([ 1.0, 3.0]))
        vertices.append(mesh.Node([ 2.0, 3.0]))
        vertices.append(mesh.Node([-2.0, 5.0]))
        vertices.append(mesh.Node([ 0.0, 4.0]))
        vertices.append(mesh.Node([ 2.0, 5.0]))
        vertices.append(mesh.Node([ 0.0, 5.0]))
        
        # make some elements
        list_of_elements = []
        list_of_elements.append(mesh.Element([vertices[0],
                                              vertices[1],
                                              vertices[5],
                                              vertices[4],
                                              vertices[3]]))
        list_of_elements.append(mesh.Element([vertices[1],
                                              vertices[2],
                                              vertices[7],
                                              vertices[6],
                                              vertices[5]]))
        list_of_elements.append(mesh.Element([vertices[3],
                                              vertices[4],
                                              vertices[9],
                                              vertices[8]]))
        list_of_elements.append(mesh.Element([vertices[4],
                                              vertices[5],
                                              vertices[6],
                                              vertices[10],
                                              vertices[13],
                                              vertices[9]]))
        list_of_elements.append(mesh.Element([vertices[6],
                                              vertices[7],
                                              vertices[11],
                                              vertices[10]]))
        list_of_elements.append(mesh.Element([vertices[8],
                                              vertices[9],
                                              vertices[13],
                                              vertices[15],
                                              vertices[12]]))
        list_of_elements.append(mesh.Element([vertices[10],
                                              vertices[11],
                                              vertices[14],
                                              vertices[15],
                                              vertices[13]]))
        
        simple_mesh = mesh.Mesh(vertices, list_of_elements)

        simple_mesh.assign_frame_ids_in_order()
        
        self.assert_( not simple_mesh.elements[3].check_if_on_boundary() )
        
        for element in simple_mesh.elements:
            if element.id_in_frame != 3:
                self.assert_( element.check_if_on_boundary() )
                
        simple_mesh.remove_boundary_elements()
        self.assertEqual(simple_mesh.get_num_elements(), 1)
        self.assertEqual(simple_mesh.get_num_nodes(), 6)
        self.assertEqual(simple_mesh.elements[0].id_in_frame, 3)
        
        for node in simple_mesh.nodes:
            self.assertEqual( len(node.get_adjacent_element_ids()) , 1)

#     def test_copy_large_mesh(self):
#         """test whether we can use deepcopy() to copy a mesh"""
#         
#         random_mesh = mesh.creation.generate_random_tesselation(30, 30, 0)
#         random_mesh.assign_frame_ids_randomly()
# 
#         mesh_copy = random_mesh.copy()
#         
#         node_5_position = random_mesh.nodes[5].position
#         self.assertNotAlmostEqual(node_5_position[0], 1.01, msg = 'randomly ran into this, just run the test again!')
#         self.assertNotAlmostEqual(node_5_position[1], 1.01, msg = 'randomly ran into this, just run the test again!')
# 
#         for node_index, node in enumerate(random_mesh.nodes):
#             np.testing.assert_almost_equal(node.position, mesh_copy.nodes[node_index].position)       
#             self.assertEqual(node.id, mesh_copy.nodes[node_index].id)
#             
#         for element_index, element in enumerate(random_mesh.elements):
#             self.assertEqual(element.id_in_frame, mesh_copy.elements[element_index].id_in_frame)
#             copied_element = mesh_copy.elements[element_index]
#             for node_index, node in enumerate(element.nodes):
#                 np.testing.assert_almost_equal(node.position, mesh_copy.nodes[node_index].position)       
#                 self.assertEqual(node.id, copied_element.nodes[node_index].id)
#                 
# 
#         mesh_copy.nodes[5].position = np.array([1.01, 1.01])
# 
#         self.assertNotAlmostEqual(random_mesh.nodes[5].position[0], 1.01)
#         self.assertNotAlmostEqual(random_mesh.nodes[5].position[1], 1.01)
#         
#         frame_id_of_7th_element = random_mesh.elements[7].id_in_frame
#         
#         mesh_copy.elements[7].id_in_frame = 5000
#         
#         self.assertEqual(random_mesh.elements[7].id_in_frame, frame_id_of_7th_element)
#         self.assertNotEqual(random_mesh.elements[7].id_in_frame, 5000)
#         
#         for node in enumerate(random_mesh.nodes):
#             self.assertEquals(random_mesh.node.position, mesh_copy.node.position)