# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__),'..','..','src'))

import unittest
import mesh
import tracking
from os import path
from os.path import dirname
from nose.plugins.attrib import attr
 
class TestPlotWithData(unittest.TestCase):
                                 
    def test_plot_with_data(self):
        """Test to plot an overlay of data"""

        first_filename = path.join(dirname(__file__),'data','first_few_frames', 'Segment_0_000.tif')
        second_filename = path.join(dirname(__file__),'data','first_few_frames', 'Segment_0_001.tif')
        third_filename = path.join(dirname(__file__),'data','first_few_frames', 'Segment_0_002.tif')

        first_mesh = mesh.read_frame_from_data(first_filename)
        second_mesh = mesh.read_frame_from_data(second_filename)
        third_mesh = mesh.read_frame_from_data(third_filename)
        
        first_real_data_filename = path.join(dirname(__file__), 'data', 'first_data_frame.tif')
        second_real_data_filename = path.join(dirname(__file__), 'data', 'second_data_frame.tif')
        third_real_data_filename = path.join(dirname(__file__), 'data', 'third_data_frame.tif')

        first_filename_to_save = path.join(dirname(__file__), 'output','first_overlay.pdf')
        second_filename_to_save = path.join(dirname(__file__), 'output','second_overlay.pdf')
        third_filename_to_save = path.join(dirname(__file__), 'output','third_overlay.pdf')
        

        first_mesh.plot_with_data( first_filename_to_save , first_real_data_filename )
        second_mesh.plot_with_data( second_filename_to_save , second_real_data_filename )
        third_mesh.plot_with_data( third_filename_to_save , third_real_data_filename )
