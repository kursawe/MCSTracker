# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This test tests our initial mesh creation framework.
"""
import unittest
import mesh
from os import path
from os.path import dirname
import os
import glob
from nose.plugins.attrib import attr
import matplotlib.pyplot as plt
import numpy as np
import cv2
from mesh.in_out import _natural_keys

@attr(level = 'standard')
class TestConversion(unittest.TestCase):
 
    def test_convert_from_ilastik_runs(self):
        "test that we can convert ilastik-type segmentation to seedwater-type segmentation"
        
        mesh.convert_from_ilastik( path.join(dirname(__file__), 'data', 'image_data' ),
                                   path.join(dirname(__file__), 'data', 'ilastik_data' ),
                                   out_path = path.join(dirname(__file__), 'output', 'converted' ) )
        
        files_in_output_dir = os.listdir(path.join(dirname(__file__), 'output', 'converted') )
        self.assertEqual(len(files_in_output_dir), 21)
        for filepath in files_in_output_dir:
            filename = os.path.split(filepath)[1]
            assert(filename.startswith('segmented_CropStack200'))
            assert(filename.endswith('.tif'))
            this_image = plt.imread( path.join(dirname(__file__), 'output', 'converted', filename) )
            self.assertEqual(this_image.dtype, np.dtype('uint16'))
 
    def test_create_seeds_from_image(self):
        "test whether our create seeds function works"
        first_image = plt.imread(path.join(dirname(__file__),
                                           'data/ilastik_data/CropStack20001_Simple Segmentation.tif'))
        
        seed_image = mesh.create_seeds_from_image( first_image )
        
        plt.imsave( path.join(dirname(__file__), 
                    'output/testseeds.tif'), seed_image )

    def test_create_segmentation_from_seeds(self):
        "test whether we can create a segmentation from given seeds."

        first_image = plt.imread(path.join(dirname(__file__),
                                           'data/ilastik_data/CropStack20001_Simple Segmentation.tif'))

        seed_image = mesh.create_seeds_from_image( first_image )
        
        actual_image = plt.imread(path.join(dirname(__file__),
                                           'data/image_data/CropStack20001.tif'))
        
        segmented_image = mesh.create_segmentation_from_seeds( actual_image, seed_image )
        
        self.assertEqual( segmented_image.dtype, np.dtype('uint16') )
        
        plt.imsave( path.join(dirname(__file__), 
                    'output/testsegmentation.tif'), segmented_image )
        
        cv2.imwrite( path.join(dirname(__file__), 
                     'output/testsegmentationwithcv2.tif'), segmented_image )
        
        reloaded_image = cv2.imread( path.join(dirname(__file__), 
                         'output/testsegmentationwithcv2.tif'), flags = -1 )
        
        np.testing.assert_equal( segmented_image, reloaded_image )

    def test_plot_first_contour_of_first_frame(selfs):

        folder_name = path.join(dirname(__file__), 'output','converted')

        first_filename = path.join(folder_name, 'segmented_CropStack20000.tif')
        
        first_image = cv2.imread( first_filename, flags = -1 ) 

        contour_list, cell_ids = mesh.get_contour_list(first_image)
        
        helper_image = np.zeros_like(first_image, dtype = 'uint8')
        
        cv2.drawContours( helper_image, np.array( contour_list[0] ), -1, color = 1, thickness = -1 )
    
        plt.imsave(path.join(dirname(__file__), 'output', 'testcontour.tif'), helper_image )
        plt.imsave(path.join(dirname(__file__), 'output', 'testfirstimage.tif'), first_image ) 

    def test_overlay(self):
        "overlay ilastik segmentation with real data for visual inspection"
        
        mesh_sequence = mesh.read_sequence_from_data(path.join(dirname(__file__), 'output','converted'))
        
        real_sequence = glob.glob(path.join(dirname(__file__), 'data', 'image_data', '*.tif'))
        
        real_sequence.sort(key=_natural_keys)

        output_dir = path.join(dirname(__file__),'output','overlays')
        if not path.isdir(output_dir):
            os.mkdir(output_dir)
        
        for mesh_counter, mesh_instance in enumerate( mesh_sequence ):
            print mesh_counter
            name_of_real_image = real_sequence[mesh_counter]
            mesh_instance.plot_with_data( path.join(output_dir, str(mesh_counter) + '.pdf' ), 
                                          name_of_real_image )
            