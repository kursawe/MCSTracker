# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/LICENSE.

"""This test tests our initial mesh creation framework.
"""
import unittest
import mesh
from os import path
from os.path import dirname
from nose.plugins.attrib import attr

@attr(level = 'standard')
class TestMesh(unittest.TestCase):
 
    def setUp(self):
        pass
    
    def test_make_simple_hexagonal_mesh(self):
        """create a 2 by 2 hexagonal mesh and test the number
        of elements, their number of nodes, etc"""
        simple_mesh = mesh.creation.generate_hexagonal(2,2)

        # plotting doesn't work without ids assigned to it

        simple_mesh.plot(path.join(dirname(__file__),'output','hexagonal_mesh.pdf'))
        self.assertEqual(simple_mesh.get_num_nodes(), 16)
        self.assertEqual(simple_mesh.get_num_elements(), 4)
        for element in simple_mesh.elements:
            self.assertEqual(element.get_num_nodes(), 6)
        
    def test_nodes_are_ordered_counter_clockwise_in_each_element(self):
        """Create a simple mesh and make sure that the nodes are ordered clockwise
        in each element.
        """
        simple_mesh = mesh.creation.generate_hexagonal(7,10)
        for element in simple_mesh.elements:
            self.assertEqual(element.calculate_area(), 1.0)

        second_simple_mesh = mesh.creation.generate_hexagonal(3,3)
        for element in second_simple_mesh.elements:
            self.assertEqual(element.calculate_area(), 1.0)
            
    def test_generate_vornoi_padding_centroids(self):
        """test whether the random mesh generator draws the padding centroids correctly"""
        
        random_mesh_generator = mesh.creation.RandomTesselationGenerator(5,3)
        padding_centroids = random_mesh_generator.generate_padding_voronoi_centroids()
        
        for centroid in padding_centroids:
            assert( ( centroid[0] < -3.5 and centroid[0] > -5.5 ) or
                    ( centroid[0] > 9.5 and centroid[0] < 11.5 ) or
                    ( centroid[1] < -3.5 and centroid[1] > -5.5 ) or
                    ( centroid[1] > 7.5 and centroid[1] < 9.5 ) )
            
            if ( centroid[0] < -3.5 and centroid[0] > -5.5 ) :
                assert( centroid[1] > -5.5 and centroid[1] < 9.5 )
            elif ( ( centroid[0] > 9.5 and centroid[0] < 11.5 ) ):
                assert( centroid[1] > -5.5 and centroid[1] < 9.5 )
            elif ( ( centroid[1] < -3.5 and centroid[1] > -5.5 ) ):
                assert( centroid[0] > -5.505 and centroid[1] < 11.5 )
            elif ( ( centroid[1] > 7.5 and centroid[1] < 9.5 ) ):
                assert( centroid[0] > -5.505 and centroid[1] < 11.5 )
    
    def test_create_random_tesselation(self):
        """Create a mesh following Luisma Escuderos algorithm for
        mesh creation using voronoi tesselations"""
        
        random_mesh = mesh.creation.generate_random_tesselation(15,10, 5)

        random_mesh.plot(path.join(dirname(__file__),'output','random_mesh.pdf'))
        for element in random_mesh.elements:
            self.assertGreater(element.calculate_area(), 0.0)
            for node in element.nodes:
                self.assertNotAlmostEqual(node.position[0], -10.101)
                self.assertNotAlmostEqual(node.position[1], -10.101)
            this_centroid = element.calculate_centroid()
            self.assertGreater(this_centroid[0] , 0.0)
            self.assertGreater(this_centroid[1] , 0.0)
            self.assertLess(this_centroid[0] , 15.5)
            self.assertLess(this_centroid[1] , 10.5)

        average_cell_area = random_mesh.calculate_average_element_area()
        self.assertGreater(average_cell_area, 0.9)
        self.assertLess(average_cell_area, 1.1)
        
        total_area = random_mesh.calculate_total_area()
        self.assertGreater(total_area, 135)
        self.assertLess(total_area, 165.0)
        
        self.assertGreater(random_mesh.get_num_elements(), 135)
        self.assertLess(random_mesh.get_num_elements(), 165)
        

