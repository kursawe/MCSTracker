# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This test tests our initial mesh creation framework.
"""

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__),'..','..','src'))

import unittest
import mesh
import numpy as np
from os import path
from os.path import dirname
from nose.plugins.attrib import attr

@attr(level = 'standard')
class TestDeathEvent(unittest.TestCase):
 
    def test_kill_element(self):
        """Test whether we can perform a death event"""

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

        self.assertEqual(simple_mesh.get_num_nodes(), 16)
        self.assertEqual(simple_mesh.get_num_elements(), 7)
        
        simple_mesh.assign_frame_ids_in_order()
        
        simple_mesh.plot(path.join(dirname(__file__),'output','death_test_before.pdf'))
        
        simple_mesh.kill_element_with_frame_id(3)
        
        simple_mesh.plot(path.join(dirname(__file__),'output','death_test_after.pdf'))

        self.assertEqual(simple_mesh.get_num_nodes(), 11)
        self.assertEqual(simple_mesh.get_num_elements(), 6)

        self.assertEquals(simple_mesh.get_element_with_frame_id(0).get_num_nodes(), 4)
        self.assertEquals(simple_mesh.get_element_with_frame_id(1).get_num_nodes(), 4)
        self.assertEquals(simple_mesh.get_element_with_frame_id(2).get_num_nodes(), 3)
        self.assertEquals(simple_mesh.get_element_with_frame_id(4).get_num_nodes(), 3)
        self.assertEquals(simple_mesh.get_element_with_frame_id(5).get_num_nodes(), 4)
        self.assertEquals(simple_mesh.get_element_with_frame_id(6).get_num_nodes(), 4)
        
        new_node_exists = False
        for node in simple_mesh.nodes:
            if np.allclose( node.position, [0.0,2.5] ):
                new_node_exists = True
                new_node = node
                break
                
        assert(new_node_exists)
        
        self.assertEquals(len(new_node.adjacent_elements), 6)
        
        # We want to test that the new node got added in the right order to all elements. If it isn't added at the right position to the
        # elements, the elements will have a different area, and hence their total area will not add up to the total area 
        # that the mesh should have.
        
        self.assertAlmostEqual( simple_mesh.calculate_total_area(), 20.0) 