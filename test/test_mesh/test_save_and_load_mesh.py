# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/LICENSE.

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
class TestSaveAndLoadMesh(unittest.TestCase):
 
    def test_save_and_load(self):
        """make a mesh, save, and load it"""
        
        mesh_to_save = mesh.creation.generate_random_tesselation(3,3)
        
        mesh_to_save.save(path.join(dirname(__file__),'output','test_mesh.mesh'))

        loaded_mesh = mesh.load(path.join(dirname(__file__),'output','test_mesh.mesh'))
        
        assert( loaded_mesh == mesh_to_save )
        
        self.assertAlmostEqual( loaded_mesh.calculate_total_area(), mesh_to_save.calculate_total_area() )