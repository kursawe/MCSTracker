# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__),'..','..','src'))

import unittest
import mesh
 
class TestColorMap(unittest.TestCase):
 
    def test_simple_color_map(self):
        my_colormap = mesh.core._get_distinct_colors(5)
        self.assertEquals(len(my_colormap),5)
        for color in my_colormap:
            self.assertEqual(len(color), 3)