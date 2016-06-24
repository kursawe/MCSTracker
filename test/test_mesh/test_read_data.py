# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This test tests our initial mesh creation framework.
"""
import unittest
import mesh
from os import path
from os.path import dirname
from nose.plugins.attrib import attr

@attr(level = 'standard')
class TestReadData(unittest.TestCase):
 
    def test_read_single_frame(self):
        filename = path.join(dirname(__file__),'data','first_few_frames', 'Segment_0_000.tif')
        this_mesh = mesh.read_frame_from_data(filename)
        this_network = this_mesh.generate_network()
        self.assertGreater( this_mesh.get_num_elements(), 10)

        for this_element in this_mesh.elements:
            self.assertGreater(this_element.get_num_nodes(), 2)
            self.assertGreater(this_element.id_in_frame, 1)
            if not this_element.check_if_on_boundary():
                self.assertEqual( this_element.get_num_nodes(), this_network.degree(this_element.id_in_frame))
            
        for node in this_mesh.nodes:
            self.assertGreaterEqual( node.position[0], 0)
            self.assertGreaterEqual( node.position[1], 0)
            self.assertGreater( len(node.adjacent_elements) , 0 )
        this_mesh.plot(path.join(dirname(__file__), 'output', 'first_data_mesh.pdf'))
 
    def test_read_sequence(self):
        """Test whether we can read the seedwater output"""
        
        mesh_sequence = mesh.read_sequence_from_data(path.join(dirname(__file__),'data','first_few_frames')) 
        
        self.assertEqual(len(mesh_sequence), 3)

        for frame_number, this_mesh in enumerate( mesh_sequence ):
            print 'reading frame number ' + str(frame_number)
            this_network = this_mesh.generate_network()
            this_mesh.plot(path.join(dirname(__file__),'output','first_frames_' + str(frame_number) + '.pdf'))
            this_network = this_mesh.generate_network()
            self.assertGreater( this_mesh.get_num_elements(), 10)
            for this_element in this_mesh.elements:
                self.assertGreater(this_element.get_num_nodes(), 2)
                self.assertGreater(this_element.id_in_frame, 1)
                self.assert_(this_network.node[this_element.id_in_frame]['position'] != None )
                if not this_element.check_if_on_boundary():
                    self.assertEqual( this_element.get_num_nodes(), this_network.degree(this_element.id_in_frame))

