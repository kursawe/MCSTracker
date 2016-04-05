# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/LICENSE.

"""This test tests our initial mesh creation framework.
"""
import unittest
import mesh
import numpy as np
from os import path
from os.path import dirname
from nose.plugins.attrib import attr

@attr(level = 'standard')
class TestMeshFurther(unittest.TestCase):
 
    def test_t1_swap(self):
        """Test whether we can perform a t1 swap"""

        vertices = []
        vertices.append(mesh.Node([ 0.0, 0.0]))
        vertices.append(mesh.Node([ 1.0, 1.0]))
        vertices.append(mesh.Node([ 0.0, 2.0]))
        vertices.append(mesh.Node([-1.0, 1.0]))
        vertices.append(mesh.Node([ 3.0, 2.5]))
        vertices.append(mesh.Node([-3.0, 2.5]))
        vertices.append(mesh.Node([ 0.0, 3.0]))
        vertices.append(mesh.Node([ 1.0, 4.0]))
        vertices.append(mesh.Node([-1.0, 4.0]))
        vertices.append(mesh.Node([ 0.0, 5.0]))
        
        # make some elements
        list_of_elements = []
        list_of_elements.append(mesh.Element([vertices[0],
                                              vertices[1],
                                              vertices[2],
                                              vertices[3]]))
        list_of_elements.append(mesh.Element([vertices[1],
                                              vertices[4],
                                              vertices[7],
                                              vertices[6],
                                              vertices[2]]))
        list_of_elements.append(mesh.Element([vertices[3],
                                              vertices[2],
                                              vertices[6],
                                              vertices[8],
                                              vertices[5]]))
        list_of_elements.append(mesh.Element([vertices[6],
                                              vertices[7],
                                              vertices[9],
                                              vertices[8]]))
        
        simple_mesh = mesh.Mesh(vertices, list_of_elements)

        self.assertEqual(simple_mesh.get_num_nodes(), 10)
        self.assertEqual(simple_mesh.get_num_elements(), 4)
        
        simple_mesh.assign_frame_ids_in_order()

        simple_mesh.plot(path.join(dirname(__file__), 'output', 't1_swap_before.pdf'))
        simple_mesh.perform_t1_swap(2,6)
        
        for element in simple_mesh.elements:
            if element.id_in_frame in [1,2]:
                self.assertEqual(element.get_num_nodes(), 4)
            if element.id_in_frame in [0,3]:
                self.assertEqual(element.get_num_nodes(), 5)
        
        distance_between_the_new_nodes = np.linalg.norm( simple_mesh.nodes[2].position - 
                                                         simple_mesh.nodes[6].position )
        
        centre_of_the_two_nodes = np.mean( np.vstack( (simple_mesh.nodes[2].position, 
                                                       simple_mesh.nodes[6].position) ), axis = 0)
        
#         self.assert_( simple_mesh.nodes[2].get_adjacent_element_ids()
        simple_mesh.plot(path.join(dirname(__file__), 'output', 't1_swap_after.pdf'))
        self.assertAlmostEqual(distance_between_the_new_nodes, 0.2)
        np.testing.assert_almost_equal(centre_of_the_two_nodes, [0.0, 2.5])
        
        for local_index, node in enumerate(simple_mesh.elements[0].nodes):
            if np.allclose( node.position, [0.1, 2.5]):
                np.testing.assert_almost_equal( simple_mesh.elements[0].nodes[(local_index+1)%5].position,
                                                [-0.1, 2.5])
                self.assert_( 0 in node.get_adjacent_element_ids())
                self.assert_( 1 in node.get_adjacent_element_ids())
                self.assert_( 3 in node.get_adjacent_element_ids())
                self.assert_( 2 not in node.get_adjacent_element_ids())
                break
            else:
                self.assertNotEqual(local_index, 4)

        for local_index, node in enumerate(simple_mesh.elements[3].nodes):
            if np.allclose( node.position, [-0.1, 2.5] ):
                np.testing.assert_almost_equal( simple_mesh.elements[3].nodes[(local_index+1)%5].position,
                                                [0.1, 2.5])
                self.assert_( 2 in node.get_adjacent_element_ids())
                self.assert_( 0 in node.get_adjacent_element_ids())
                self.assert_( 3 in node.get_adjacent_element_ids())
                self.assert_( 1 not in node.get_adjacent_element_ids())
                break
            else:
                self.assertNotEqual( local_index, 4 )
