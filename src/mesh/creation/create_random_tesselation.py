# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.


import numpy as np
from ..core import *
from pyhull.voronoi import VoronoiTess # in maths type easy_install --user pyhull to make this work

def generate_random_tesselation(nx, ny, number_of_relaxation_steps = 4):
    """Generate a random tesselation of roughly nx times ny cells
        
    Parameters
    ----------
    nx : int
        Approximate number of cells in x direction

    ny : int
        Approximate number of cells in y direction
        
    number_of_relaxation_steps : int
        In each relaxation step the centroids of the voronoi cells
        are used as seeds for a new voronoi tesselation
        
    Returns
    -------
    mesh : Mesh type
        random mesh of dimension nx and ny
        
    THIS DESCRIPTION NEEDS UPDATING, SEE CLASS DOCUMENTATION BELOW

    This mesh generator will distribute random points in a plane of size
    nx x ny and and create a voronoi tesselation from it. Then, number_of_relaxation_steps
    relaxation steps will be applied. In each relaxation step, the centroids of all voronoi cells
    will be collected and used as seeds for a new tesselation.
    """
    
    mesh_generator = RandomTesselationGenerator(nx, ny, number_of_relaxation_steps)
    
    return mesh_generator.generate_tesselation()
 
class RandomTesselationGenerator:
    def __init__(self, nx, ny, number_of_relaxation_steps = 4):
        """A generator class for the generation of random tesslations.
    
        For more detailed documentation of the input arguments see the
        generate_random_tesselation function in this module
        
        The constructor will not perform the voronoi tesselation and Lloyd relaxation.
        To get the mesh call the function generate_tesselation() of this object.
        """
        self.nx = nx
        """width of the sheet in x direction"""

        self.ny = ny
        """width of the sheet in y direction"""
        
        self.nx_with_dummies = nx + 8
        """ We make the voronoi mesh bigger by 4 cells on each side in order
        to avoid boundary effects in the mesh that we draw. This is then the
        overall mesh dimension in x directions"""

        self.ny_with_dummies = ny + 8
        """The overall mesh dimension in y direction"""

        self.number_of_relaxation_steps = number_of_relaxation_steps
        """The number of voronoi relaxation steps that the generator should perform on the tesselation"""

        self.current_mesh = None
        """after each iteration the generator will reset this internal variable"""
        
        self.current_centroids = None
        """will be reset after each iteration"""

        self.padding_centroid_positions = self.generate_padding_voronoi_centroids()
        """At each iteration, we put these additional centroids ad the outside of the mesh"""
        
        self.voronoi_diagram = None
        """the full voronoi diagram after the current iteration"""
        
    def generate_tesselation(self):
        """Generate and return the tesselation as prescribed by the arguments of the __init__() function
        
        Returns
        -------
        
        random_mesh : Mesh type
            a mesh object with based on random voronoi seeds of the requested size, 
            and where the requested number of Lloyd relaxation steps has been applied
        
        """

        # First, we generate an initial tesselation. This is where the random seeds are made and the internal
        # mesh variable are first written
        self.generate_initial_tesselation()
        
        # if the user does not request relaxation steps, just return this mesh
        if ( self.number_of_relaxation_steps == 0 ):
            # we need to crop the mesh before we return it, i.e. cut of the regions and cells that we added 
            # to hide boundary effects
            return self.crop_and_return_mesh()
        else:
            # if relaxation steps are required we perform them
            for iteration in range(self.number_of_relaxation_steps):

                # get the centroids of the existing mesh
                centroid_positions = []
                for element in self.mesh.elements:
                    centroid_positions.append(element.calculate_centroid())
                    
                # add our padding centroids, i.e. two rows of regularly spaced centroids
                # on the outside of the entire region of interest (which includes extra space and cells)
                all_centroid_positions = np.vstack((np.array(centroid_positions),
                                                    self.padding_centroid_positions))
    
                # We generate a Voronoi tesselation of that
                self.voronoi_diagram = VoronoiTess(all_centroid_positions)
                
                # now we clean up that voronoi tesselation and set the mesh in the generator
                self.remove_padding_cells_from_voronoi_diagram_and_create_mesh()

            # after we went through the iterations we can crop and return the mesh
            return self.crop_and_return_mesh()
                
    def crop_and_return_mesh(self):
        """We generated a (Lloyd relaxed) voronoi mesh that is larger than the requested mesh.
        Here, we cut off all the extra cells
        
        The method works by drawing a completely new mesh."""
        
        # Let's make a list of new elements for the new mesh
        new_elements = []
        
        # a dictionary where keys are ids of old nodes and values are new nodes
        # at the same position. New nodes are required because we currently don't have functionality to
        # delete elements from nodes (shouldn't be difficult to implement, but also didn't seem necessary
        # here)
        node_dictionary = {}
        
        for old_element in self.mesh.elements:

            this_centroid = old_element.calculate_centroid()

            # if this centroid is inside the box of interest, create an identical new element
            if ( not( this_centroid[0] < 0.5 or this_centroid[0] > ( self.nx + 0.5 )
                      or this_centroid[1] < 0.5 or this_centroid[1] > (self.ny + 0.5) )):
                     
                new_element_nodes = []
                
                for node in old_element.nodes:
                    if node.id not in node_dictionary:
                        # The new node has the same position as the old node
                        node_dictionary[node.id] = Node(node.position)

                    new_element_nodes.append(node_dictionary[node.id])
                
                new_elements.append(Element(new_element_nodes))
                
        return Mesh(node_dictionary.values(), new_elements)
   
    def generate_initial_tesselation(self):
        """Generate the initial tesselation.
        
        In this method we distribute nx_with_dummies*ny_with_dummies points randomly in a box of 
        dimension [nx_with_dummies, ny_with_dummies]. We then add the padding centroids on the outside
        of the box to `cut of' infinity. Finally, we take the voronoi tesselation of all these centroids,
        and remove the voronoi cells originating from the padding centroids from this tesselation.
        This is designed such that the average area of the remaining cells should be one.
        """

        total_number_of_random_points_to_create = np.prod(self.nx_with_dummies * self.ny_with_dummies)
    
        # We make a mesh of cell centroids.
        # First we get an array with random entries between 0 and 1 in x and y directions of correct length
        centroid_positions = np.random.rand(total_number_of_random_points_to_create, 2)
    
        # Then we blow this array up to the correct dimensions in x and y
    
        centroid_positions[:,0] *= self.nx_with_dummies
        centroid_positions[:,1] *= self.ny_with_dummies
        
        # and we shift the mesh such that the bottom and left added space is negative
        centroid_positions[:,0] -= 3.5
        centroid_positions[:,1] -= 3.5
        # the region with the centroids now spreads from -3.5 to (self.[nx/ny] + 4.5) so that we can add padding
        # cells at integer coordinates, i.e -4 and self.[nx/ny] + 5

        # Now we add some padding cells on the outside.
    
        centroid_positions = np.vstack((centroid_positions, self.padding_centroid_positions))
    
        # We generate a Voronoi tesselation of that
        self.voronoi_diagram = VoronoiTess(centroid_positions)
        
        # now we clean up this voronoi tesselation and set the mesh of the generator
        # cleaning up means we remove the voronoi cells originating from the padding seeds
        self.remove_padding_cells_from_voronoi_diagram_and_create_mesh()
    
    def remove_padding_cells_from_voronoi_diagram_and_create_mesh(self):
        """This function interrogates the internal voronoi_diagram and creates a Mesh class object instance from it.
        This mesh will not include the cells that were added for padding to cut off infinity.
        """

        # Let's get all vertices
        all_vertices = np.array(self.voronoi_diagram.vertices)
    
        # and all regions
        all_regions = self.voronoi_diagram.regions
    
        # now, let's make a new array of vertices and a new array of regions
    
        new_vertices = []
        new_elements = []
        index_map = [-1]*all_vertices.shape[0]
        current_number_of_new_vertices = 0

        # Loop over all regions and filter all regions that come from the dummy indices
        # this relies on the voronoi diagram indexing the regions and voronoi points in the same order
        
        for region_index, voronoi_centre in enumerate(self.voronoi_diagram.points):
            # if this voronoi centre is inside the box (including the dummies)
            if ( not( voronoi_centre[0] < -3.5 or voronoi_centre[0] > self.nx + 4 + 0.5 
                      or voronoi_centre[1] < -3.5 or voronoi_centre[1] > self.ny + 4 + 0.5 )):
                new_region = []
                for vertex_index in self.voronoi_diagram.regions[region_index]:
                    # if the index is not in the new vertices, add it
                    if (index_map[vertex_index] == -1):
                        this_node = Node(all_vertices[vertex_index])
                        new_vertices.append(this_node)
                        index_map[vertex_index] = current_number_of_new_vertices
                        # and add it to the new region
                        new_region.append(this_node)
                        current_number_of_new_vertices += 1
                    else:
                        # if it is in the new vertices though then add
                        # the right vertex to the new region
                        new_region.append(new_vertices[index_map[vertex_index]])
                # lets make an element from that region and add it to the list of elements
                this_element = Element(new_region)
                new_elements.append(this_element)

        # ensure that all the vertices in each element are ordered counterclockwise
        for element in new_elements:
            area = element.calculate_area()
            if area < 0:
                element.nodes.reverse()
            
        # now we have made all the nodes and elements for the mesh, and associated them.
        # Time to generate the mesh and set the member variable for further processing
        self.mesh = Mesh(new_vertices, new_elements)
 
    def generate_padding_voronoi_centroids(self):
        """Generate 2 rows or colomns on all sides of the nx_with_dummies x ny_with_dummies grid to seal
        the voronoi tesselation from infinity.
    
        The returned centroids will form two rows around the grid, the inner row
        will have integer coordinates, the outer row will be shifted by 0.5 in x direction
        on the top and the bottom, and shifted by 0.5 in y direction on the left and the right.

        Parameters
        ----------
    
        nx : size of the grid in x direction
        ny : size of the grid in y direction
    
        Returns
        -------
    
        padding_centroids : numpy array
            list of centroids of the points making up the padding. 
        """

        # 2 extra rows of cells at the top
        x_positions = np.arange( -5, self.nx + 7 )
    
        top_row_one = np.zeros( ( self.nx_with_dummies + 4, 2 ) )
        top_row_one[:, 0] = x_positions
        top_row_one[:, 1] = self.ny + 5

        top_row_two = np.zeros( ( self.nx_with_dummies + 4 + 1, 2 ) )
        top_row_two[:,1] = self.ny + 6 
        top_row_two[1:, 0] = (x_positions + 0.5)
        top_row_two[0, 0] = -5.5

        # 2 extra rows of cells at the bottom
        bottom_row_one = np.zeros( ( self.nx_with_dummies + 4, 2 ) )
        bottom_row_one[:, 0] = x_positions
        bottom_row_one[:, 1] = - 4

        bottom_row_two = np.zeros( ( self.nx_with_dummies + 4 + 1, 2 ) )
        bottom_row_two[:, 1] = - 5
        bottom_row_two[1:, 0] = (x_positions + 0.5)
        bottom_row_two[0,0] = -5.5
    
        # 2 extra rows left
    
        y_positions = np.arange( -3 , self.ny + 5 )

        left_row_one = np.zeros( ( self.ny_with_dummies , 2 ) )
        left_row_one[:,1] = y_positions
        left_row_one[:,0] = -4
    
        left_row_two = np.zeros( ( self.ny_with_dummies + 1, 2) )
        left_row_two[1:, 1] = y_positions + 0.5
        left_row_two[0, 1] = -4.5
        left_row_two[:, 0] = -5
    
        # and 2 extra rows on the right

        right_row_one = np.zeros( ( self.ny_with_dummies , 2) )
        right_row_one[:,1] = y_positions
        right_row_one[:,0] = self.nx + 5
        
        right_row_two = np.zeros( ( self.ny_with_dummies + 1, 2) )
        right_row_two[1:, 1] = y_positions + 0.5
        right_row_two[0, 1] = - 4.5
        right_row_two[:, 0] = self.nx + 6
    
        all_padding_centroids = np.vstack((right_row_one, right_row_two,
                                           left_row_one, left_row_two,
                                           bottom_row_one, bottom_row_two,
                                           top_row_one, top_row_two))
        
        return all_padding_centroids
    