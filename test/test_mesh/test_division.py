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
class TestDivisionEvent(unittest.TestCase):
 
    def test_divide_element_in_direction(self):
        """Test whether we can perform a division event"""

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
        
        simple_mesh.plot(path.join(dirname(__file__),'output','division_test_before.pdf'))
        
        area_of_element_to_divide = simple_mesh.get_element_with_frame_id(3).calculate_area()

        simple_mesh.divide_element_with_frame_id_in_direction(3, [0.5, 1.0])
        
        simple_mesh.plot(path.join(dirname(__file__),'output','division_test_after.pdf'))

        self.assertEquals(simple_mesh.get_element_with_frame_id(0).get_num_nodes(), 6)
        self.assertEquals(simple_mesh.get_element_with_frame_id(1).get_num_nodes(), 5)
        self.assertEquals(simple_mesh.get_element_with_frame_id(2).get_num_nodes(), 4)
        self.assertEquals(simple_mesh.get_element_with_frame_id(4).get_num_nodes(), 4)
        self.assertEquals(simple_mesh.get_element_with_frame_id(5).get_num_nodes(), 5)
        self.assertEquals(simple_mesh.get_element_with_frame_id(6).get_num_nodes(), 6)
        
        first_new_node_exists = False
        second_new_node_exists = False
        for node in simple_mesh.nodes:
            if np.allclose( node.position, [0.5,3.5] ):
                first_new_node_exists = True
                first_new_node = node
            if np.allclose( node.position, [-0.5,1.5] ):
                second_new_node_exists = True
                second_new_node = node
                
        assert(first_new_node_exists)
        assert(second_new_node_exists)
        
        self.assertEquals(len(first_new_node.adjacent_elements), 3)
        self.assertEquals(len(second_new_node.adjacent_elements), 3)
        
        daughter_element_ids = set.intersection( set( first_new_node.get_adjacent_element_ids() ), 
                                                 set( second_new_node.get_adjacent_element_ids() ))
        
        self.assertEquals(len(daughter_element_ids), 2)
        
        first_daughter_element = simple_mesh.get_element_with_frame_id( daughter_element_ids.pop() )
        self.assertAlmostEqual( first_daughter_element.get_num_nodes(), 5)

        second_daughter_element = simple_mesh.get_element_with_frame_id( daughter_element_ids.pop() )
        self.assertAlmostEqual( second_daughter_element.get_num_nodes(), 5)
        
        self.assertAlmostEqual( first_daughter_element.calculate_area() + second_daughter_element.calculate_area(),
                                area_of_element_to_divide)

    def test_make_triangular_element_in_division(self):
        """Test whether we can perform a division event"""

        vertices = []
        vertices.append(mesh.Node([-2.0, 0.0]))
        vertices.append(mesh.Node([ 0.0, 0.0]))
        vertices.append(mesh.Node([ 2.0, 0.0]))
        vertices.append(mesh.Node([-2.0, 2.0]))
        vertices.append(mesh.Node([-1.0, 2.0]))
        vertices.append(mesh.Node([ 0.0, 1.0]))
        vertices.append(mesh.Node([ 1.0, 2.0]))
        vertices.append(mesh.Node([ 2.0, 2.0]))
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
        list_of_elements.append(mesh.Element([vertices[4],
                                              vertices[5],
                                              vertices[6],
                                              vertices[9]]))
        list_of_elements.append(mesh.Element([vertices[9],
                                              vertices[11],
                                              vertices[8],
                                              vertices[3],
                                              vertices[4]]))
        list_of_elements.append(mesh.Element([vertices[10],
                                              vertices[11],
                                              vertices[9],
                                              vertices[6],
                                              vertices[7]]))
        
        simple_mesh = mesh.Mesh(vertices, list_of_elements)

        self.assertEqual(simple_mesh.get_num_nodes(), 12)
        self.assertEqual(simple_mesh.get_num_elements(), 5)
        
        simple_mesh.assign_frame_ids_in_order()
        
        simple_mesh.plot(path.join(dirname(__file__),'output','triangular_division_test_before.pdf'))
        
        area_of_element_to_divide = simple_mesh.get_element_with_frame_id(3).calculate_area()

        simple_mesh.divide_element_with_frame_id_in_direction(2, [1.0, 0.0])
        
        simple_mesh.plot(path.join(dirname(__file__),'output','triangular_division_test_after.pdf'))

        self.assertEquals(simple_mesh.get_element_with_frame_id(0).get_num_nodes(), 5)
        self.assertEquals(simple_mesh.get_element_with_frame_id(1).get_num_nodes(), 5)
        self.assertEquals(simple_mesh.get_element_with_frame_id(3).get_num_nodes(), 6)
        self.assertEquals(simple_mesh.get_element_with_frame_id(4).get_num_nodes(), 6)
        self.assertEquals(simple_mesh.get_element_with_frame_id(5).get_num_nodes(), 3)
        self.assertEquals(simple_mesh.get_element_with_frame_id(6).get_num_nodes(), 5)
        
    def test_divide_element_in_random_direction(self):
        """test whether we can divide an element in a random direction"""
        
        my_mesh = mesh.creation.generate_random_tesselation(7,7)
        
        # First pick an element in the centre
        mesh_centre = my_mesh.calculate_centre()
        
        my_mesh.assign_frame_ids_in_order()

        # pick the element closest to the centre
        min_distance = 3*my_mesh.calculate_height()

        for element in my_mesh.elements:
            distance = np.linalg.norm( element.calculate_centroid() - mesh_centre )
            if distance < min_distance:
               min_distance = distance
               most_central_element = element 
               
        old_num_elements = my_mesh.get_num_elements()

        old_num_nodes = my_mesh.get_num_nodes()
        
        old_element_area = most_central_element.calculate_area()

        my_mesh.plot(path.join(dirname(__file__),'output','before_random_division.pdf'))
        my_mesh.divide_element_with_frame_id(most_central_element.id_in_frame)
        my_mesh.plot(path.join(dirname(__file__),'output','after_random_division.pdf'))
        
        self.assertEqual(old_num_elements + 1, my_mesh.get_num_elements())
        self.assertEqual(old_num_nodes + 2, my_mesh.get_num_nodes())
        
        daughter_element_one = most_central_element
        daughter_element_two = my_mesh.elements[old_num_elements]
        
        self.assertAlmostEqual( daughter_element_one.calculate_area() + daughter_element_two.calculate_area() ,
                                old_element_area )
        