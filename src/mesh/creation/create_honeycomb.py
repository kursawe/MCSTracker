# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This module contains functions to generate meshes with defined shapes
"""

import numpy as np
from ..core import *
from pyhull.voronoi import VoronoiTess # in maths type easy_install --user pyhull to make this work

def generate_hexagonal(nx,ny):
    """Create a hexagonal mesh of dimension nx x ny.
        
    Parameters
    ----------
    nx : int
        Number of cells in x direction
    ny : int
        Number of cells in y direction
        
    Returns
    -------
    mesh : Mesh type
        hexagonal mesh of dimension nx, ny
    """
    # We specifiy the dimension of the grid in x and y direction
    number_of_cells = np.array((nx,ny))
    number_of_cells_and_dummies = number_of_cells + 2
    
    # We make a mesh of cell centroids
    centroid_positions = np.zeros((number_of_cells_and_dummies.prod(),2), 'double')
    origin = -number_of_cells_and_dummies/2.0+0.5;
    counter = 0;
    for x_index in range(number_of_cells_and_dummies[0]):
        for y_index in range(number_of_cells_and_dummies[1]):
            if y_index%2==0:
                centroid_positions[counter,0] = origin[0] + x_index*1.0;
            else:
                centroid_positions[counter,0] = origin[0] + 0.5 + x_index*1.0;
            centroid_positions[counter,1] = origin[1] + y_index*1.0;
            counter += 1
    
    # We generate a Voronoi tesselation of that
    voronoi_diagram = VoronoiTess(centroid_positions)
    
    # Now, let's get all vertices
    all_vertices = np.array(voronoi_diagram.vertices)
    
    # and all regions
    
    all_regions = voronoi_diagram.regions
    
    # now, let's make a new array of vertices and a new array of regions
    
    new_vertices = []
    new_elements = []
    index_map = [-1]*all_vertices.shape[0]
    current_number_of_new_vertices = 0
    # Loop over all regions
    for region in all_regions:
        # If this region is not a boundary region
        if (0 not in region) and (len(region) == 6):
            # then we want to keep it
            new_region = []
            for index in region:
                # if the index is not in the new vertices, add it
                if (index_map[index] == -1):
                    this_node = Node(all_vertices[index])
                    new_vertices.append(this_node)
                    index_map[index] = current_number_of_new_vertices
                    # and add it to the new region
                    new_region.append(this_node)
                    current_number_of_new_vertices += 1
                else:
                    # if it is in the new vertices though then add
                    # the right vertex to the new region
                    new_region.append(new_vertices[index_map[index]])
            # lets make an element from that region and add it to the list of elements
            new_elements.append(Element(new_region))
                    
    # ensure that all the vertices in each element are ordered counterclockwise
    for element in new_elements:
        area = element.calculate_area()
        if area < 0:
            element.nodes.reverse()
        area = element.calculate_area()
        assert( abs(area - 1.0) < 1e-10 )
            
    # now we can create a mesh and return it
    return Mesh(new_vertices, new_elements)
