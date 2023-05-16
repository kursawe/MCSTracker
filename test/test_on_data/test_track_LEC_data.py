# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This tests our first tracking example
"""
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__),'..','..','src'))

import unittest
import mesh
import tracking
import time
from os import path
from os.path import dirname
from nose.plugins.attrib import attr
import sys
import threading
#import resource
 
class TestTrackLECData(unittest.TestCase):
                                 
      
    def xest_track_first_frames(self):

        my_path_1 = os.path.join(dirname(__file__),'..','test_tracking','data','AL151222_Job6_0.mesh')
        my_path_2 = os.path.join(dirname(__file__),'..','test_tracking','data','AL151222_Job6_1.mesh')

        mesh_one = mesh.load(my_path_1)
        mesh_two = mesh.load(my_path_2)

        tracked_ids = tracking.track( mesh_one, mesh_two )

        mesh_one.plot('tracked_mesh_before_t1_swap.pdf', color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
        mesh_two.plot('tracked_mesh_after_t1_swap.pdf', color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
         
    def test_track_raw_segmentation_1(self): 
        
        if not os.path.isdir(path.join(dirname(__file__),'output','tracked_LEC_segmentation_1/')):
            os.mkdir(path.join(dirname(__file__),'output','tracked_LEC_segmentation_1/'))
                     
        tracking.track_and_write_sequence(path.join(dirname(__file__),'data','LEC_segmentation_1'), path.join(dirname(__file__),'output','tracked_LEC_segmentation_1/'),
                                        use_geometry=True)
                                        # use_geometry=True, number_meshes= 2)

        tracking.plot_tracked_sequence(path.join(dirname(__file__),'output','tracked_LEC_segmentation_1/'), 
                                path.join(dirname(__file__),'data', 'LEC_segmentation_1_raw_images/'),
                                path.join(dirname(__file__),'data','LEC_segmentation_1/'),
                                path.join(dirname(__file__),'output','tracked_LEC_segmentation_1/'))

    def xest_track_other_frames(self):

        my_path_1 = os.path.join(dirname(__file__),'..','test_tracking','data','AL151222_Job6_0.mesh')
        my_path_2 = os.path.join(dirname(__file__),'..','test_tracking','data','AL151222_Job6_1.mesh')

        mesh_one = mesh.load(my_path_1)
        mesh_two = mesh.load(my_path_2)

        tracked_ids = tracking.track( mesh_one, mesh_two )

        mesh_one.plot('tracked_mesh_before_t1_swap.pdf', color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )
        mesh_two.plot('tracked_mesh_after_t1_swap.pdf', color_by_global_id = True, 
                      total_number_of_global_ids = len(tracked_ids) )        
