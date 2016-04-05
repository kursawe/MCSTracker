# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/LICENSE.

"""This tests our first tracking example
"""
import mesh
from mesh import Node
from mesh import Element
from mesh import Mesh
import tracking
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
 
def make_step_one_plot():
    print 'hello'
    
    
    mesh_generator = mesh.creation.RandomTesselationGenerator(5, 5, 5)

    total_number_of_random_points_to_create = np.prod(mesh_generator.nx_with_dummies * 
                                                      mesh_generator.ny_with_dummies)
    
    # We make a mesh of cell centroids.
    # First we get an array with random entries between 0 and 1 in x and y directions of correct length
    centroid_positions = np.random.rand(total_number_of_random_points_to_create, 2)
    
    # Then we blow this array up to the correct dimensions in x and y
    
    centroid_positions[:,0] *= mesh_generator.nx_with_dummies
    centroid_positions[:,1] *= mesh_generator.ny_with_dummies

    # and we shift the mesh such that the bottom and left added space is negative
    centroid_positions[:,0] -= 3.5
    centroid_positions[:,1] -= 3.5
    # the region with the centroids now spreads from -3.5 to (self.[nx/ny] + 4.5) so that we can add padding
    # cells at integer coordinates, i.e -4 and self.[nx/ny] + 5

    # Now we add some padding cells on the outside.
    
    centroid_positions = np.vstack((centroid_positions, mesh_generator.padding_centroid_positions))
    
    # We generate a Voronoi tesselation of that
    mesh_generator.voronoi_diagram = mesh.creation.VoronoiTess(centroid_positions)
    
    # Let's get all vertices
    all_vertices = np.array(mesh_generator.voronoi_diagram.vertices)
    
    # and all regions
    all_regions = mesh_generator.voronoi_diagram.regions
    
    # now, let's make a new array of vertices and a new array of regions
    
    new_vertices = []
    new_elements = []
    index_map = [-1]*all_vertices.shape[0]
    current_number_of_new_vertices = 0

    # Loop over all regions and filter all regions that come from the dummy indices
    # this relies on the voronoi diagram indexing the regions and voronoi points in the same order
    
    region_centroids = []
    for region_index, voronoi_centre in enumerate(mesh_generator.voronoi_diagram.points):
        # if this voronoi centre is inside the box (including the dummies)
        new_region = []
        if ( not( voronoi_centre[0] < -4.5 or voronoi_centre[0] > mesh_generator.nx + 5 + 0.5 
                      or voronoi_centre[1] < -4.5 or voronoi_centre[1] > mesh_generator.ny + 5 + 0.5 )):
            new_region = []
            for vertex_index in mesh_generator.voronoi_diagram.regions[region_index]:
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
    mesh_generator.mesh = Mesh(new_vertices, new_elements)

    filtered_region_centroids = []
    polygon_list = []
    for element in mesh_generator.mesh.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
         
        this_polygon.set_facecolor('white')
        polygon_list.append(this_polygon)
        centroid = element.calculate_centroid()
        if ( not( centroid[0] < -3.5 or centroid[0] > mesh_generator.nx + 4 + 0.5 
             or centroid[1] < -3.5 or centroid[1] > mesh_generator.ny + 4 + 0.5 )):
            filtered_region_centroids.append(centroid)

    filtered_region_centroids_np = np.array(filtered_region_centroids)

    polygon_collection = mpl.collections.PatchCollection(polygon_list, match_original = True)
    polygon_collection.set_edgecolor('darkgrey')
    mesh_figure = plt.figure()
    mesh_figure.gca().add_collection(polygon_collection)
    plt.scatter(centroid_positions[:,0], centroid_positions[:,1], marker = '.', s = 5)
    plt.scatter(filtered_region_centroids_np[:,0], filtered_region_centroids_np[:,1], marker = 'x', s = 15, color = 'darkgrey')

    # outer line
    plt.plot((-3.5, mesh_generator.nx_with_dummies -3.5), (-3.5, -3.5), 'k-', lw = 2)
    plt.plot((-3.5, mesh_generator.nx_with_dummies -3.5), (mesh_generator.ny_with_dummies -3.5, mesh_generator.ny_with_dummies -3.5), 'k-', lw = 2)
    plt.plot((-3.5, -3.5), (-3.5, mesh_generator.ny_with_dummies -3.5), 'k-', lw = 2)
    plt.plot((mesh_generator.nx_with_dummies -3.5, mesh_generator.nx_with_dummies -3.5), (-3.5, mesh_generator.ny_with_dummies -3.5), 'k-', lw = 2)

    #inner line
#     plt.plot((0.5, mesh_generator.nx +0.5), (0.5, 0.5), 'k--', lw = 2)
#     plt.plot((0.5, mesh_generator.nx +0.5), (mesh_generator.ny +0.5, mesh_generator.ny +0.5), 'k--', lw = 2)
#     plt.plot((0.5, 0.5), (0.5, mesh_generator.ny +0.5), 'k--', lw = 2)
#     plt.plot((mesh_generator.nx +0.5, mesh_generator.nx +0.5), (+0.5, mesh_generator.ny +0.5), 'k--', lw = 2)
    mesh_figure.gca().set_aspect('equal')
    mesh_figure.gca().autoscale_view()
    plt.axis('off')
    mesh_figure.savefig('initial_tesselation.pdf', bbox_inches = 'tight')
    plt.close(mesh_figure)
    
#### SECOND LLOYDS RELAXATION STEP

    centroid_positions = np.vstack((filtered_region_centroids_np, mesh_generator.padding_centroid_positions))
    
    # We generate a Voronoi tesselation of that
    mesh_generator.voronoi_diagram = mesh.creation.VoronoiTess(centroid_positions)
    
    # Let's get all vertices
    all_vertices = np.array(mesh_generator.voronoi_diagram.vertices)
    
    # and all regions
    all_regions = mesh_generator.voronoi_diagram.regions
    
    # now, let's make a new array of vertices and a new array of regions
    
    new_vertices = []
    new_elements = []
    index_map = [-1]*all_vertices.shape[0]
    current_number_of_new_vertices = 0

    # Loop over all regions and filter all regions that come from the dummy indices
    # this relies on the voronoi diagram indexing the regions and voronoi points in the same order
    
    region_centroids = []
    for region_index, voronoi_centre in enumerate(mesh_generator.voronoi_diagram.points):
        # if this voronoi centre is inside the box (including the dummies)
        new_region = []
        if ( not( voronoi_centre[0] < -4.5 or voronoi_centre[0] > mesh_generator.nx + 5 + 0.5 
                      or voronoi_centre[1] < -4.5 or voronoi_centre[1] > mesh_generator.ny + 5 + 0.5 )):
            new_region = []
            for vertex_index in mesh_generator.voronoi_diagram.regions[region_index]:
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
    mesh_generator.mesh = Mesh(new_vertices, new_elements)

    filtered_region_centroids = []
    polygon_list = []
    for element in mesh_generator.mesh.elements:
        this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                           fill = True)
         
        this_polygon.set_facecolor('white')
        polygon_list.append(this_polygon)
        centroid = element.calculate_centroid()
        if ( not( centroid[0] < 0.5 or centroid[0] > ( mesh_generator.nx + 0.5 )
                  or centroid[1] < 0.5 or centroid[1] > (mesh_generator.ny + 0.5) )):
#             this_polygon.set_facecolor('#2dec34')
            this_polygon.set_facecolor('royalblue')
        if ( not( centroid[0] < -3.5 or centroid[0] > mesh_generator.nx + 4 + 0.5 
             or centroid[1] < -3.5 or centroid[1] > mesh_generator.ny + 4 + 0.5 )):
            filtered_region_centroids.append(centroid)

    filtered_region_centroids_np = np.array(filtered_region_centroids)

    polygon_collection = mpl.collections.PatchCollection(polygon_list, match_original = True)
    polygon_collection.set_edgecolor('darkgrey')
    mesh_figure = plt.figure()
    mesh_figure.gca().add_collection(polygon_collection)
    plt.scatter(centroid_positions[:,0], centroid_positions[:,1], marker = '.', s = 5)
#     plt.scatter(filtered_region_centroids_np[:,0], filtered_region_centroids_np[:,1], marker = 'x', s = 4.5, color = 'darkgrey')

    # outer line
    plt.plot((-3.5, mesh_generator.nx_with_dummies -3.5), (-3.5, -3.5), 'k-', lw = 2)
    plt.plot((-3.5, mesh_generator.nx_with_dummies -3.5), (mesh_generator.ny_with_dummies -3.5, mesh_generator.ny_with_dummies -3.5), 'k-', lw = 2)
    plt.plot((-3.5, -3.5), (-3.5, mesh_generator.ny_with_dummies -3.5), 'k-', lw = 2)
    plt.plot((mesh_generator.nx_with_dummies -3.5, mesh_generator.nx_with_dummies -3.5), (-3.5, mesh_generator.ny_with_dummies -3.5), 'k-', lw = 2)

    #inner line
    plt.plot((0.5, mesh_generator.nx +0.5), (0.5, 0.5), 'k--', lw = 2)
    plt.plot((0.5, mesh_generator.nx +0.5), (mesh_generator.ny +0.5, mesh_generator.ny +0.5), 'k--', lw = 2)
    plt.plot((0.5, 0.5), (0.5, mesh_generator.ny +0.5), 'k--', lw = 2)
    plt.plot((mesh_generator.nx +0.5, mesh_generator.nx +0.5), (+0.5, mesh_generator.ny +0.5), 'k--', lw = 2)
    mesh_figure.gca().set_aspect('equal')
    mesh_figure.gca().autoscale_view()
    plt.axis('off')
    mesh_figure.savefig('second_tesselation.pdf', bbox_inches = 'tight')
    plt.close(mesh_figure)

if __name__=="__main__":
    make_step_one_plot()
