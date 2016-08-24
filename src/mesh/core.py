# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""In this file the main mesh classes 'Mesh', 'Element', and 'Node' are defined.
"""

import numpy as np
import networkx as nx
import itertools
import colorsys
import warnings
import pickle
import sys
import cv2

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

class Mesh():
    """A mesh similar to a mutable vertex mesh in Chaste or an unstructured grid in VTK.
    """
    def __init__(self, nodes, elements):
        """The standard generator.
        
        The generator will number all nodes in the order they are provided,
        and it will index frame ids if all elements are provided with frame ids.
        (frame ids are identifies for elements in a given frame).
        
        Parameters
        ----------
        nodes : a list of nodes
            Each entry is a node.    
        elements : a list of element type objects
            Elements referencing the vertices vector.
            
        Returns
        -------
        mesh : the mesh
        """

        self.nodes = nodes
        """A list with all the nodes"""
        
        self.elements = elements
        """A list with all the elements"""
        
        self.frame_id_dictionary = {}
        """A dictionary with element frame ids as keys and the location
        in the element vector as values. Warning: currently, no checking is
        done as to whether the mapping between elements and frame ids is one-to-one"""

        self.global_id_dictionary = {}
        """A dictionary with global element ids as keys and the location
        in the element vector as values. The global ids have to be assigned
        to cells by a tracking algorithm, and the method 'index_global_ids' has to be called
        in order to get elements from the mesh by their global ids. Warning: currently, no checking is
        done as to whether the mapping between elements and frame ids is one-to-one.
        """
         
        self.assign_node_ids_in_order()
        self.build_frame_id_dictionary()
        
    def get_num_nodes(self):
        """Get the number of nodes in the mesh.
        
        Returns
        -------
        number_of_vertices : int
            Number of vertices in the mesh.
        """
        return len(self.nodes)

    def get_num_elements(self):
        """Get the number of elements in the mesh.
        
        Returns
        -------
        number_of_vertices : int
            Number of elements in the mesh.
        """
        return len(self.elements)
    
    def plot(self, filename, color_by_global_id = False, total_number_of_global_ids = None, reduced_mcs_only = False):
        """Plots the mesh to a file using matplotlib
        
        Parameters
        ----------
        
        filename : string
            name of the file to which we should plot
        
        color_by_global_id : bool or matplotlib color
            If this is true, the color will not be selected by the frame id, but by the global id.
            This will only work if the input parameter total_number_of_global_ids is also specified.
            Otherwise the function will not know how many colors the colormap will need to have.
            If this is false, the color will be random as specified by the frame ids.
            If this is a matplotlib color, cells that have a global id will show up in this color.
            
        total_number_of_global_ids : int
            The total number of global ids present under the current mapping. If color_by_global_id is true,
            he function will divide the available color space into total_number_of_global_ids and then for each 
            cell use the global_idth color in this discretization for plotting. The colors between frames 
            will not match if this parameter is incorrectly specified.
            
        Returns
        -------
        
        polygon_collection : matplotlib.collections.PatchCollection
            You can plot this particular polygon collection using the line
            my_figure.gca().add_collection(polygon_collection)
            where my_figure is an active matplotlib figure.

        """
        # figure properties
        figuresize = (4,2.75)
        font = {'size'   : 10}
        plt.rc('font', **font)

        mesh_figure = plt.figure(figsize = figuresize)
        polygon_collection = self.get_polygon_collection( color_by_global_id, total_number_of_global_ids, 
                                                          reduced_mcs_only )
        mesh_figure.gca().add_collection(polygon_collection)
        mesh_figure.gca().set_aspect('equal')
        mesh_figure.gca().autoscale_view()
        if filename.endswith('.pdf'):
            mesh_figure.savefig(filename, bbox_inches = 'tight')
        else:
            mesh_figure.savefig(filename, bbox_inches = 'tight', dpi = 300)
        
    def get_polygon_collection( self, color_by_global_id = False, total_number_of_global_ids = None,
                                reduced_mcs_only = False ):
        """Helper method for plot command, allows to include mesh representations
           in your own matplotlib plots
        
        Parameters
        ----------

        color_by_global_id : bool or matplotlib color
            If this is true, the color will not be selected by the frame id, but by the global id.
            This will only work if the input parameter total_number_of_global_ids is also specified.
            Otherwise the function will not know how many colors the colormap will need to have.
            If this is false, the color will be random as specified by the frame ids.
            If this is a matplotlib color, cells that have a global id will show up in this color.
        
        total_number_of_global_ids : int
            The total number of global ids present under the current mapping. If color_by_global_id is true,
            he function will divide the available color space into total_number_of_global_ids and then for each 
            cell use the global_idth color in this discretization for plotting. The colors between frames 
            will not match if this parameter is incorrectly specified.
            
        Returns
        -------
        
        polygon_collection : matplotlib.collections.PatchCollection
            You can plot this particular polygon collection using the line
            my_figure.gca().add_collection(polygon_collection)
            where my_figure is an active matplotlib figure.

        """
        polygon_list = []
        if color_by_global_id == True:
            if total_number_of_global_ids == None:
                raise Exception("You need to provide the total_number_of_global_ids argument!")
            else:
                color_collection = _get_distinct_colors( total_number_of_global_ids + 1 ) 
        elif color_by_global_id == False:
            if self.get_maximal_frame_id() is not None:
                color_collection = _get_distinct_colors(self.get_maximal_frame_id() + 1)

        for element in self.elements:
            this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                               fill = True)
            if color_by_global_id == True:

                if element.global_id == None:
                    this_polygon.set_facecolor([1.0, 1.0, 1.0])
                else:
                    if reduced_mcs_only: 
                        if element.is_in_reduced_mcs_next and not element.is_in_reduced_mcs_previous:
                            this_polygon.set_facecolor('yellow')
                        elif element.is_in_reduced_mcs_previous and not element.is_in_reduced_mcs_next:
                            this_polygon.set_facecolor('blue')
                        elif element.is_in_reduced_mcs_next and element.is_in_reduced_mcs_previous:
                            this_polygon.set_facecolor('green')
                        else:
                            this_polygon.set_facecolor('white')
                    else:
                        this_polygon.set_facecolor(color_collection[element.global_id])

            elif color_by_global_id == False:

                if element.id_in_frame == None:
                    this_polygon.set_facecolor([1.0, 1.0, 1.0])
                else:
                    this_polygon.set_facecolor( color_collection[element.id_in_frame] )

            else:
                if element.global_id == None:
                    this_polygon.set_facecolor([1.0, 1.0, 1.0])
                else:
                    this_polygon.set_facecolor(color_by_global_id)

            polygon_list.append(this_polygon)

        polygon_collection = mpl.collections.PatchCollection(polygon_list, match_original = True)
        
        return polygon_collection   

    def assign_node_ids_in_order(self):
        """Assign ids to nodes in the order they appear in the node vector
        """
        for counter, node in enumerate(self.nodes):
            node.id = counter

    def assign_frame_ids_in_order(self):
        """Assign frame_ids to elements in the order they appear in the elements vector
        """
        for counter, element in enumerate(self.elements):
            element.id_in_frame = counter
            self.frame_id_dictionary[counter] = counter

    def assign_frame_ids_in_reverse_order(self):
        """Assign ids to elements in the reverse order of how they appear in the elements vector
        """
        number_of_elements = len(self.elements)
        for counter, element in enumerate(self.elements):
            element.id_in_frame = number_of_elements - 1 - counter
            self.frame_id_dictionary[number_of_elements - 1 - counter] = counter

    def assign_frame_ids_randomly(self):
        """Assign ids to elements in random order
        """
        all_ids = np.linspace(0,self.get_num_elements()-1, self.get_num_elements() )
        all_ids = all_ids.astype('int')
        np.random.shuffle(all_ids)
        for counter, element in enumerate(self.elements):
            element.id_in_frame = all_ids[counter]
            self.frame_id_dictionary[all_ids[counter]] = counter
            
    def count_rosettes(self):
        """Count the number of rosettes in the mesh. A rosette is a node
        that belongs to more than three elements.
        
        Returns
        -------
        
        number_of_rosettes : int
            number of rosettes in the tissue
        """
        number_of_rosettes = 0
        for this_node in self.nodes:
            number_elements = len(this_node.adjacent_elements)
            if number_elements > 4:
                number_of_rosettes += 1
        
        return number_of_rosettes
                
    def build_frame_id_dictionary(self):
        """Index the frame ids in a dictionary.
           Writes empty dictionary if there is at least one frame id that is set to `None'
           in the mesh.
        """
        self.frame_id_dictionary = {}
        for element_index, element in enumerate(self.elements):
            if element.id_in_frame==None:
                self.frame_id_dictionary = {}
                break
            else:
                self.frame_id_dictionary[element.id_in_frame] = element_index 

    def get_element_with_frame_id(self, id_in_frame):
        """Returns the element that has the given id
        
        Parameters
        ----------
        
        id_in_frame : int
            frame id of the element that should be returned
            
        Returns
        -------
        
        element : the element that has the given frame id
        
        Currently, no checking is done as to whether the mapping elements <-> ids is one-to-one
        """
        assert ( id_in_frame in self.frame_id_dictionary )
        return self.elements[self.frame_id_dictionary[id_in_frame]]
    
    def collect_edges(self):
        """Crawl through the elements and collect all edges
        
        Returns
        -------
        
        edges : numpy integer array
            An numpy array where each row contains an edge.
            Each edge is characterised by two integer values, 
            the ID's of the nodes that it contains.
        """
        edge_list = []
        for element in self.elements:
            for index, node in enumerate(element.nodes):
                next_index = (index + 1)%element.get_num_nodes()
                next_node = element.nodes[next_index]
                this_edge = [node.id, next_node.id]
                #the list() command will make an identical copy of the list
                this_edge_reversed = list(this_edge)
                this_edge_reversed.reverse()
                if (this_edge not in edge_list) and\
                    (this_edge_reversed not in edge_list):
                    edge_list.append(this_edge)
        return np.array(edge_list)

    def collect_elements_of_inner_edges(self):
        """Collect all element pairs that share an inner edge of the mesh.
        
        Returns
        -------
        
        inner_edge_sets : list of sets
            list of sets. Each set contains two integers which are frame ids
            of elements shared by one edge.
        """

        all_edges = self.collect_edges()
        
        inner_edge_sets = []
        for edge in all_edges:
            first_node = self.get_node_with_id(edge[0])
            second_node = self.get_node_with_id(edge[1])
            element_ids_at_first_node = set( first_node.get_adjacent_element_ids() )
            element_ids_at_second_node = set( second_node.get_adjacent_element_ids() )
            shared_element_ids = set.intersection( element_ids_at_first_node, 
                                                   element_ids_at_second_node )
            if len(shared_element_ids) == 2:
                inner_edge_sets.append(shared_element_ids)
        
        return inner_edge_sets

    def get_nodes_shared_by_elements(self, element_list ):
        """Get nodes that are shared by the elements in element_set
        
        Parameters
        ----------
        
        element_set : list of ints
            list containing exactly two frame ids of elements in the mesh
            
        Returns
        -------
        
        shared_nodes : set of ints
            set containing all ids of nodes that are shared by the elements
        """
    
        first_element = self.get_element_with_frame_id( element_list[0] )
        second_element = self.get_element_with_frame_id( element_list[1] )
        
        first_element_node_ids = [ node.id for node in first_element.nodes ]
        second_element_node_ids = [ node.id for node in second_element.nodes ]
        
        shared_nodes = set.intersection( set(first_element_node_ids),
                                         set(second_element_node_ids))
        
        return shared_nodes

    def generate_network(self):
        """Return a network of cell ID's
        
        Returns
        -------
        
        network : networkx graph object
            Network nodes are element frame ids
            node attributes are:
                'position' : the centroid of the cell

                'num_neighbours' : the number of mesh nodes that the element has. Note that this is different to the
                    degree of the element node in the generated network
        """
        network_edges = []

#         for node in self.nodes:
#             new_connections_between_cells = itertools.combinations(node.get_adjacent_element_ids(), 2)
#             for this_tuple in new_connections_between_cells:
#                 connection = list(this_tuple)
#                 connection_reversed = list(connection)
#                 connection_reversed.reverse()
#                 if ( connection not in network_edges ) and\
#                     (connection_reversed not in network_edges):
#                     network_edges.append(connection)
#          
        edges = self.collect_edges()
        for edge in edges:
            cells_along_this_edge = self.find_elements_at_edge(edge)
            if len(cells_along_this_edge) > 1:
                assert( len(cells_along_this_edge) == 2)
                network_edges.append(cells_along_this_edge)
             
        network = nx.Graph()
        network.add_edges_from(network_edges)
        
        for element in self.elements:
            if network.has_node(element.id_in_frame):
                network.node[element.id_in_frame]['position'] = element.calculate_centroid()
                network.node[element.id_in_frame]['num_neighbours'] = element.get_num_nodes()
        
        return network
    
    def find_elements_at_edge(self, edge):
        """Find all cells that are at this edge.
        
        Parameters
        ----------
        
        edge : [int, int]
            list of node_ids along this edge.
            
        Returns
        -------
        
        elements_at_edge : list of ints
            list of element frame ids that share this edge
        """
        
        node_one = self.get_node_with_id(edge[0])
        node_two = self.get_node_with_id(edge[1])

        adjacent_element_set = set.intersection( set( node_one.get_adjacent_element_ids() ), 
                                                      set( node_two.get_adjacent_element_ids() ))
        
        elements_at_edge = list( adjacent_element_set )

        return elements_at_edge
        
    def generate_network_of_unidentified_elements(self, identified_elements = [] ):
        """Returns a network of all nodes that have not yet been identified
        
        Parameters
        ----------
        
        identified_elements : list of ints
            frame_ids in the network that have been identified but may not have
            a global id

        Returns
        -------
        
        unidentified_network : networkx Graph instance
            network of all elements that have global_id == None and are not contained in 
            identified_elements
            
        See also
        --------
        
        generate_network()
        """
        
        complete_network = self.generate_network()
        
        nodes_to_remove = []
        for node in complete_network.nodes():
            if ( self.get_element_with_frame_id(node).global_id != None or 
                 node in identified_elements):
                nodes_to_remove.append(node)
                
        complete_network.remove_nodes_from( nodes_to_remove )
        
        return complete_network

    def generate_network_of_identified_elements(self, identified_elements = [] ):
        """Returns a network of all nodes that have not yet been identified
        
        Parameters
        ----------
        
        identified_elements : list of ints
            frame_ids in the network that have been identified but may not have
            a global id

        Returns
        -------
        
        unidentified_network : networkx Graph instance
            network of all elements that have global_id == None and are not contained in 
            identified_elements
            
        See also
        --------
        
        generate_network()
        """
        
        complete_network = self.generate_network()
        
        nodes_to_remove = []
        for node in complete_network.nodes():
            if ( self.get_element_with_frame_id(node).global_id is None and 
                 node not in identified_elements):
                nodes_to_remove.append(node)
                
        complete_network.remove_nodes_from( nodes_to_remove )
        
        return complete_network
           
    def calculate_width(self):
        """Calculate the width of the network
        
        Returns
        -------
        
        width : double
            width of the network defined as the difference between the maximal and minimal node position
            in x direction
        """
        
        x_positions_vector = np.zeros(self.get_num_nodes())
        for index, node in enumerate(self.nodes):
            x_positions_vector[index] = node.position[0]
        
        width = x_positions_vector.max() - x_positions_vector.min()
        
        return width

    def calculate_height(self):
        """Calculate the height of the network
        
        Returns
        -------
        
        height : double
            height of the network defined as the difference between the maximal and minimal node position
            in y direction
        """
        
        y_positions_vector = np.zeros(self.get_num_nodes())
        for index, node in enumerate(self.nodes):
            y_positions_vector[index] = node.position[1]
        
        height = y_positions_vector.max() - y_positions_vector.min()
        
        return height
    
    def calculate_centre(self):
        """Calculate the centroid of the mesh as defined by the average of all node positions
        
        Returns
        -------
        
        centroid : 1x2 double np.array
            the average position of all nodes
        """
        
        node_positions = [node.position for node in self.nodes]
            
        node_positions = np.array(node_positions)
        
        centroid = np.mean(node_positions, axis = 0)
        
        return centroid
    
    def index_global_ids(self):
        """build a lookup table for which element belongs to given global id.
        
        After elements have been assigned global ids we need to run this method in order
        to be able to interrogate the mesh for elements belonging to a global id
        """
        self.global_id_dictionary = {}
        for index, element in enumerate(self.elements):
            if not element.global_id == None:
                self.global_id_dictionary[element.global_id] = index

    def index_frame_ids(self):
        """build a lookup table for which element belongs to a given frame id.
        
        After elements ids or the position of elements in the element vector have 
        changed we need to run this method in order to be able to interrogate the mesh 
        for elements belonging to a frame id
        """
        self.frame_id_dictionary = {}
        for index, element in enumerate(self.elements):
            if not element.id_in_frame == None:
                self.frame_id_dictionary[element.id_in_frame] = index
            
    def get_element_with_global_id(self, global_id):
        """Get the element corresponding to the given global id
        
        Parameters
        ----------
        
        global_id : int
            global id of the element that should be returned
            
        Returns
        -------
        
        element : Element
            The element object in this mesh that has the given global id
        """
        
        return self.elements[self.global_id_dictionary[global_id]]
    
    def get_max_global_id(self):
        """Finds the maximal global id in this mesh.
        
        Returns
        -------
        
        max_global_id : int
            maximal global id in the mesh
            
        Warnings
        --------
        
        redundant, not tested
        """
        max_global_id = np.max( np.array( self.global_id_dictionary.keys() ) )
        
        return max_global_id
    
    def calculate_total_area(self):
        """Returns the total area of this mesh
        
        Returns
        -------
        
        total_area : double
            The sum of the areas of all elements in the mesh
        """
        
        total_area = 0.0
        for element in self.elements:
            total_area += element.calculate_area()
        
        return total_area
            
    def calculate_average_element_area(self):
        """Returns the average area across all elements in the mesh
        
        Returns
        -------
        
        average_area : double
            The average area of all elements in the mesh, i.e. the total area
            divided by the total number of elements
        """
        
        average_area = self.calculate_total_area()/self.get_num_elements()
        
        return average_area
        
    def copy_region_of_interest(self, x_limits, y_limits):
        """Returns all cells within the region defined by x_limits, y_limits as separate mesh.
        
        This method creates a new mesh that looks like the parent mesh within the specified limits. Cells are considered
        to be inside the region of interest if their centroid is.
        The new mesh is a separate full copy and the parent mesh will not be altered by the cutting process.
        global_ids of the parent elements will not be copied, but frame_ids will.
        
        Parameters
        ----------
        
        x_limits : 2-tuple ( double, double)
            the limits of the region of interest in x_direction, has to be ordered.

        y_limits : 2-tuple ( double, double)
            the limits of the region of interest in y_direction, has to be ordered.
        
        Returns
        -------
        
        window_mesh : Mesh type
            a new mesh that has identical geometry to the parent mesh within the specified region of interest.
        """
        
        new_nodes = []

        new_elements = []
        
        # keys are old node ids and values are indices for the new_nodes vector
        node_map = {}

        for element in self.elements:
            element_centroid = element.calculate_centroid()
            if ( element_centroid[0] > x_limits[0] and element_centroid[0] < x_limits[1] and
                 element_centroid[1] > y_limits[0] and element_centroid[1] < y_limits[1] ):
                new_nodes_for_this_element = []
                for node in element.nodes:
                    if node.id in node_map:
                        new_nodes_for_this_element.append(node_map[node.id])
                    else:
                        new_node = Node(node.position)
                        new_nodes.append(new_node)
                        node_map[node.id] = new_node
                        new_nodes_for_this_element.append(new_node)
                new_elements.append(Element(new_nodes_for_this_element, element.id_in_frame))
                
        new_mesh = Mesh(new_nodes, new_elements)
        
        node_positions = []
        for node in new_mesh.nodes:
            node_positions.append(node.position)
            
        position_array = np.array(node_positions)
        
        min_x_position = np.min(position_array[:,0])
        min_y_position = np.min(position_array[:,1])
        
        for node in new_mesh.nodes:
            node.position[0] -= min_x_position
            node.position[1] -= min_y_position
            
        return new_mesh
    
    def divide_element_with_frame_id(self, frame_id):
        """Divides the element with the given frame id in a random direction.
        
        Parameters
        ----------
        
        frame_id : int
            id_in_frame of the element that we wish to divide
        """
        
        angle = np.random.uniform(0.0, 2.0*np.pi)
        
        random_direction = np.array( [np.cos(angle), np.sin(angle)] )
        
        self.divide_element_with_frame_id_in_direction(frame_id, random_direction)
    
    def divide_element_with_frame_id_in_direction(self, frame_id, direction):
        """Divide the element with given frame id in the given direction
        
        The element is divided by a new edge passing the centroid of the cell,
        and pointing in the given direction.
        
        Parameters
        ----------
        
        frame_id : int
            id_in_frame of the element to be divided
            
        direction : [double, double]
            a vector pointing in the direction of the new edge
        """
        
        element_to_divide = self.get_element_with_frame_id(frame_id)
        
        centroid = element_to_divide.calculate_centroid()
        
        num_nodes = element_to_divide.get_num_nodes()

        # entries of this list will be [node, node]
        intersecting_edges = []

        division_direction = np.array(direction)
        division_direction /= np.linalg.norm(division_direction)

        perpendicular_direction = np.array([-division_direction[1], division_direction[0]])

        # identify the edges that are being cut by the axis. These are the ones where the nodes are on different sides
        # of the division axis
        is_current_node_on_left = np.inner( element_to_divide.nodes[0].position - centroid, perpendicular_direction ) > 0
        for local_index, node in enumerate(element_to_divide.nodes):
            next_node = element_to_divide.nodes[(local_index+1)%num_nodes]
            is_next_node_on_left = np.inner(next_node.position - centroid, perpendicular_direction) >= 0
            if (is_current_node_on_left != is_next_node_on_left):
                intersecting_edges.append([node, next_node])
            is_current_node_on_left = is_next_node_on_left

        # If the axis of division does not cross two edges then we cannot proceed

        assert( len(intersecting_edges) == 2 )

        # calculate the intersection on edge one
        edge_one = intersecting_edges[0]
        
        edge_one_direction = edge_one[1].position - edge_one[0].position
        
        centroid_to_edge_node = edge_one[0].position - centroid

        determinant = division_direction[1]*edge_one_direction[0] - division_direction[0]*edge_one_direction[1] 

        intersection_one = edge_one[0].position + 1/determinant*(division_direction[0]*centroid_to_edge_node[1] -
                                                                 division_direction[1]*centroid_to_edge_node[0])*edge_one_direction 
        
        # and edge two
        edge_two = intersecting_edges[1]
        
        edge_two_direction = edge_two[1].position - edge_two[0].position
        
        centroid_to_edge_node = edge_two[0].position - centroid

        determinant = division_direction[1]*edge_two_direction[0] - division_direction[0]*edge_two_direction[1] 

        intersection_two = edge_two[0].position + 1/determinant*(division_direction[0]*centroid_to_edge_node[1] -
                                                                 division_direction[1]*centroid_to_edge_node[0])*edge_two_direction 
                                                                 
        # We calculated the intersections, these now need to get some nodes:
        
        # They are going to be new nodes at the end of the nodes vector
        id_of_new_node_one = self.get_maximal_node_id() + 1
        id_of_new_node_two = self.get_maximal_node_id() + 2

        new_node_one = Node(intersection_one, id_of_new_node_one)
        new_node_two = Node(intersection_two, id_of_new_node_two)
        
        # The new nodes will need to belong to the node vector
        self.nodes.append( new_node_one )
        self.nodes.append( new_node_two )
        
        # We will also need to draw the new elements, the old element is going to be turned into a daughter cell
        
        all_old_nodes_in_order = element_to_divide.nodes
        
        # the node lists of the new elements are going to start with the nodes on the new edge.
        # how do we find the order of these nodes? By finding out whether their difference vector points
        # in the same direction or opposite direction as the division axis!
        
        new_edge_direction = new_node_two.position - new_node_one.position

        if np.inner( new_edge_direction, division_direction ) > 0:
            nodes_for_element_one = [ new_node_one, new_node_two ]
            nodes_for_element_two = [ new_node_two, new_node_one ]
        else:
            nodes_for_element_one = [ new_node_two, new_node_one ]
            nodes_for_element_two = [ new_node_one, new_node_two ]
        
        # we find the local index of the the second node on the first edge
        for local_index, node in enumerate(all_old_nodes_in_order):
            if node.id == edge_one[1].id:
                first_old_node_local_index = local_index
                break
                
        # we loop over the old nodes starting at that node
        for local_index in range( first_old_node_local_index, first_old_node_local_index + len(all_old_nodes_in_order) ):
            current_node = all_old_nodes_in_order[ local_index%len(all_old_nodes_in_order) ]
            current_node_is_left_of_axis = np.inner( current_node.position - centroid, perpendicular_direction ) >= 0
            if current_node_is_left_of_axis:
                nodes_for_element_one.append(current_node)
            else:
                nodes_for_element_two.append(current_node)
                # the element_to_divide needs to be removed from the elements of this node
                # Not sure if the remove() function always works. If this ever breaks it's because we broke
                # identifiability of our elements.
                new_adjacent_elements_to_node = [element for element in current_node.adjacent_elements 
                                                 if element.id_in_frame != element_to_divide.id_in_frame]
                current_node.adjacent_elements = new_adjacent_elements_to_node
                
        # the old element will get the nodes for element one
        element_to_divide.nodes = nodes_for_element_one
        new_element = Element(nodes_for_element_two)
        
        # find the right frame ids for the new daughter cells
        new_frame_id_of_element_one = self.get_maximal_frame_id() + 1
        new_frame_id_of_element_two = self.get_maximal_frame_id() + 2
        
        # remove the old frame id (function argument) from the dictionary
        location_of_element_to_divide_in_list = self.frame_id_dictionary[frame_id]
        del self.frame_id_dictionary[frame_id]
        element_to_divide.id_in_frame = new_frame_id_of_element_one
        new_element.id_in_frame = new_frame_id_of_element_two
        new_element_location_in_list = self.get_num_elements()
        
        # add the new element to the list
        self.elements.append(new_element)
        
        # update the frame_id_dictionary
        self.frame_id_dictionary[new_frame_id_of_element_one] = location_of_element_to_divide_in_list
        self.frame_id_dictionary[new_frame_id_of_element_two] = new_element_location_in_list

        # The new nodes will automatically have the new_element indexed, since the element constructor
        # will have taken care of that. We still need to tell them about the element_to_divide and their
        # other neighboring element.
        
        new_node_one.adjacent_elements.append(element_to_divide)
        new_node_two.adjacent_elements.append(element_to_divide)
        
        # find the other elements at node_one and node_two, respectively
        # Attention: these are going to be sets of cardinality one
        frame_id_of_element_on_edge_one = set.intersection( set( edge_one[0].get_adjacent_element_ids() ), 
                                                            set( edge_one[1].get_adjacent_element_ids() ))
        assert(len(frame_id_of_element_on_edge_one) <= 1)
        frame_id_of_element_on_edge_two = set.intersection( set( edge_two[0].get_adjacent_element_ids() ), 
                                                            set( edge_two[1].get_adjacent_element_ids() ))
        assert(len(frame_id_of_element_on_edge_two) <= 1)
        
        # add the found elements to their corresponding new nodes
        other_element_at_node_one = self.get_element_with_frame_id(frame_id_of_element_on_edge_one.pop())
        other_element_at_node_two = self.get_element_with_frame_id(frame_id_of_element_on_edge_two.pop())
        new_node_one.adjacent_elements.append( other_element_at_node_one )
        new_node_two.adjacent_elements.append( other_element_at_node_two )

        # and also add the new nodes to their corresponding elements
        new_nodes_for_other_element_one = []
        for node in other_element_at_node_one.nodes:
            new_nodes_for_other_element_one.append(node)
            if node.id == edge_one[1].id:
                new_nodes_for_other_element_one.append(new_node_one)
        other_element_at_node_one.nodes = new_nodes_for_other_element_one
        
        new_nodes_for_other_element_two = []
        for node in other_element_at_node_two.nodes:
            new_nodes_for_other_element_two.append(node)
            if node.id == edge_two[1].id:
                new_nodes_for_other_element_two.append(new_node_two)
        other_element_at_node_two.nodes = new_nodes_for_other_element_two
 
    def kill_element_with_frame_id(self, frame_id):
        """Remove the element with the given frame id and replace it with a node.
        
        The new node will be a member of all elements that were adjacent to the killed element,
        and its location will be the centroid of the killed element.
        
        Parameters
        ----------
        
        frame_id : int
            id_in_frame of the element that is to be killed
        """

        element_to_be_killed = self.get_element_with_frame_id(frame_id)
        
        location_of_new_node = element_to_be_killed.calculate_centroid()
        
        new_node = Node(location_of_new_node)
        self.nodes.append(new_node)
        
        ids_of_adjacent_elements = element_to_be_killed.get_ids_of_adjacent_elements()
        
        for adjacent_frame_id in ids_of_adjacent_elements:
            this_element = self.get_element_with_frame_id(adjacent_frame_id)
            this_element.replace_adjacent_element_with_node(frame_id, new_node)
            
        self.remove_list_of_nodes(element_to_be_killed.nodes)
        self.remove_element_with_frame_id(frame_id)
        
    def remove_element_with_frame_id(self, frame_id):
        """Removes the element with the given frame id from self.elements.
        Does not alter any nodes or elements, but keeps the internal mapping
        frame_ids <--> location in element vector up to date.
        
        Parameters
        ----------
        
        frame_id : int
            the id of the element that should be removed from the internal list of elements
        """
        
        remaining_elements = []
        
        for element in self.elements:
            if element.id_in_frame != frame_id:
                remaining_elements.append(element)
                
        self.elements = remaining_elements
        self.index_frame_ids()

    def remove_list_of_nodes(self, node_list):
        """Remove the given nodes from the internal list of nodes. Will not alter
        any elements or nodes.
        
        This is a helper function for kill_element_with_frame_id()
        
        Parameters
        ----------
        
        node_list: list of Node instances
            all the nodes that should be removed from self.nodes()
            
        See also
        --------
        
        kill_element_with_frame_id()
        """
        list_of_node_ids = []
        for node in node_list:
            list_of_node_ids.append(node.id)
        
        remaining_nodes = []
        for node in self.nodes:
            if node.id not in list_of_node_ids:
                remaining_nodes.append(node)
                
        self.nodes = remaining_nodes
        
    def get_maximal_node_id(self):
        """Get the largest id of all nodes in the mesh.
        
        Returns
        -------
        
        maximal_node_id : int
            the largest id among all nodes in the mesh
        """
        
        all_node_ids = []

        for node in self.nodes:
            all_node_ids.append(node.id)
        
        max_node_id = np.max(np.array( all_node_ids ))
        
        return max_node_id

    def get_maximal_frame_id(self):
        """Get the largest id_in_frame among all elements in the mesh.
        
        Returns
        -------
        
        maximal_frame_id : int
            the largest id_in_frame among all elements in the mesh
        """
        
        all_frame_ids = []

        for element in self.elements:
            all_frame_ids.append(element.id_in_frame)
        
        max_element_id = np.max(np.array( all_frame_ids ))

        return max_element_id
 
    def perform_t1_swap(self, node_one_id, node_two_id):
        """Performs a t1 swap. On the nodes with the given ids.
        
        The distance between the new nodes will be 0.2

            \ /       
             |    -->   \__/
            / \         /  \

        Parameters
        ----------
        
        node_one_id : int
            id of a node on one end of the edge

        node_two_id : int
            id of a node on the other end of the edge
            
        Warning
        -------
        
        Frame ids have to have been assigned to all elements in the swap.

        """
        
        # get the nodes
        node_one = self.nodes[node_one_id]
        node_two = self.nodes[node_two_id]
                
        # and the ids of their adjacent element
        node_one_elements = set(node_one.get_adjacent_element_ids())
        node_two_elements = set(node_two.get_adjacent_element_ids())
        
        # the elements that they share
        shared_element_ids = set.intersection( node_one_elements, node_two_elements )
        assert( len( shared_element_ids ) == 2 ) 
        
        #and the elements that they only have individually
        elements_in_node_one_only = set.difference( node_one_elements, shared_element_ids )
#         assert( len(elements_in_node_one_only) == 1 )

        elements_in_node_two_only = set.difference( node_two_elements, shared_element_ids )
#         assert( len(elements_in_node_two_only) == 1 )
        
        # the centre and the direction vector between the two
        centre_of_the_two_nodes = np.mean( np.vstack( (node_one.position, node_two.position) ), axis = 0 )
        direction_between_old_nodes = (node_two.position - node_one.position)/np.linalg.norm(node_one.position - 
                                                                                             node_two.position )
        
        # we 'turn' the direction vector by 90 degrees to the right
        direction_between_new_nodes = np.array([-direction_between_old_nodes[1], direction_between_old_nodes[0]])
        
        # and based on that update the two new node locations
        new_node_one_location = centre_of_the_two_nodes - 0.1 * direction_between_new_nodes
        new_node_two_location = centre_of_the_two_nodes + 0.1 * direction_between_new_nodes
       
        node_one.position = new_node_one_location
        node_two.position = new_node_two_location
        
        # the nodes will have different adjacent elements now, so reset the element vector for each node
        node_one.adjacent_elements = []
        node_two.adjacent_elements = []
        
        # deal with the elements that previously shared the edge
        for element_id in shared_element_ids:
            this_element = self.get_element_with_frame_id( element_id )
            for local_index, node in enumerate(this_element.nodes):
                if node.id == node_one.id:
                    next_node = this_element.nodes[(local_index+1)%this_element.get_num_nodes()]
                    if next_node.id == node_two.id:
                        # this is the element that should now get node one deleted
                        old_nodes = list( this_element.nodes )
                        this_element.nodes = []
                        for old_node in old_nodes:
                            if old_node.id != node_one.id:
                                this_element.nodes.append(old_node)
                        node_two.adjacent_elements.append(this_element)
                    else:
                        # this is the element that should now get node two deleted
                        old_nodes = list( this_element.nodes )
                        this_element.nodes = []
                        for old_node in old_nodes:
                            if old_node.id != node_two.id:
                                this_element.nodes.append(old_node)
                        node_one.adjacent_elements.append(this_element)
                    break
        
        if elements_in_node_one_only != set():
            element_in_node_one_only = self.get_element_with_frame_id(elements_in_node_one_only.pop())
            # this element will now have both nodes
        
            old_nodes = list(element_in_node_one_only.nodes)
            element_in_node_one_only.nodes = []
            for node in old_nodes:
                element_in_node_one_only.nodes.append(node)
                if node.id == node_one.id:
                    element_in_node_one_only.nodes.append(node_two)
                
            node_one.adjacent_elements.append(element_in_node_one_only)
            node_two.adjacent_elements.append(element_in_node_one_only)
        
        if elements_in_node_two_only != set():
            element_in_node_two_only = self.get_element_with_frame_id(elements_in_node_two_only.pop())
            # this element will now also have both nodes
            
            old_nodes = list(element_in_node_two_only.nodes)
            element_in_node_two_only.nodes = []
            for node in old_nodes:
                element_in_node_two_only.nodes.append(node)
                if node.id == node_two.id:
                    element_in_node_two_only.nodes.append( node_one )
                    
            node_one.adjacent_elements.append( element_in_node_two_only )
            node_two.adjacent_elements.append( element_in_node_two_only )
        
    def get_already_mapped_adjacent_element_ids(self, frame_id, already_mapped_ids = []):
        """Finds all elements that have already been mapped and are adjacent to the element with frame_id
        
        Elements that have been mapped are identified by having global_id != None.
        
        Parameters
        ----------
        
        frame_id : int
            id_in_frame of the element whose neighbours we are interested in
            
        already_mapped_ids : list of ints
            frame ids in mesh for which mapping exists although there global ids have not been set.
            
        Returns
        -------
        
        mapped_ids : list of ints
            list of frame ids of all elements that are adjacent to the element with frame_id and don not have None as global_id
        """
        
        all_adjacent_elements = self.get_element_with_frame_id(frame_id).get_ids_of_adjacent_elements()
        
        mapped_ids = []
        for adjacent_frame_id in all_adjacent_elements:
            this_element = self.get_element_with_frame_id(adjacent_frame_id)
            if this_element.global_id != None or adjacent_frame_id in already_mapped_ids:
                mapped_ids.append(adjacent_frame_id)
                
        return mapped_ids
 
    def get_not_yet_mapped_shared_neighbour_ids(self, list_of_frame_ids, already_mapped_ids = []):
        """Finds all elements that have not yet been mapped and are adjacent to all elements in the list.
        
        Elements that have not yet been mapped are identified by having global_id == None.
        
        Parameters
        ----------
        
        list_of_frame_ids : list of ints
            list of id_in_frames of the elements whose shared neighbours we are interested in
            
        already_mapped_ids : list of ints
            elements for which an entry exists but for which the global ids have not yet been set
            
        Returns
        -------
        
        not_yet_mapped_ids : list of ints
            list of frame ids of all elements that are adjacent to all elements in list_of_frame_ids and have None as global_id
        """
        
        if list_of_frame_ids != []:
            set_of_shared_adjacent_elements = set( self.get_element_with_frame_id(list_of_frame_ids[0]).get_ids_of_adjacent_elements() )
        else:
            set_of_shared_adjacent_elements = set()
        
        for frame_id in list_of_frame_ids:
            this_element = self.get_element_with_frame_id(frame_id)
            set_of_shared_adjacent_elements.intersection_update( self.get_element_with_frame_id(frame_id).get_ids_of_adjacent_elements() )
            
        not_yet_mapped_ids = []
        for frame_id in set_of_shared_adjacent_elements:
            this_element = self.get_element_with_frame_id(frame_id)
            if this_element.global_id == None and frame_id not in already_mapped_ids:
                not_yet_mapped_ids.append( frame_id )
                
        return not_yet_mapped_ids

    def get_inclusive_not_yet_mapped_shared_neighbour_ids(self, list_of_frame_ids):
        """Finds all elements that have not yet been mapped that are adjacent to all elements in the list,
        and are also part of the list.
        
        Finds elements in a set of cells that are adjacent to all other elements in the set.
        Elements that have not yet been mapped are identified by having global_id == None.
        
        Parameters
        ----------
        
        list_of_frame_ids : list of ints
            list of id_in_frames of the elements whose shared neighbours we are interested in
            
        Returns
        -------
        
        not_yet_mapped_ids : list of ints
            list of frame ids of all elements that are adjacent to all elements in list_of_frame_ids and have None as global_id
        """
        
        set_of_shared_adjacent_elements = set( self.get_element_with_frame_id(list_of_frame_ids[0]).get_ids_of_adjacent_elements() )
        set_of_shared_adjacent_elements.add(list_of_frame_ids[0])
        
        for frame_id in list_of_frame_ids:
            this_element = self.get_element_with_frame_id(frame_id)
            this_set_of_adjacent_neighbours = set(self.get_element_with_frame_id(frame_id).get_ids_of_adjacent_elements())
            this_set_of_adjacent_neighbours.add(frame_id)
            set_of_shared_adjacent_elements.intersection_update( this_set_of_adjacent_neighbours )
        
        set_of_shared_adjacent_elements.intersection_update(list_of_frame_ids)
            
        not_yet_mapped_ids = []
        for frame_id in set_of_shared_adjacent_elements:
            this_element = self.get_element_with_frame_id(frame_id)
            if this_element.global_id == None:
                not_yet_mapped_ids.append( frame_id )
                
        return not_yet_mapped_ids
     
    def __eq__(self, other):
        """Overridden equality method. Returns true if the meshes are identical.

        This method is called by the built in python operator "=="
        
        Parameters
        ----------
        
        other : any type
            object that we are comparing to

        Returns
        -------
        
        meshes_are_identical : bool
            True if meshes are identical
        """
        return ( type(other) is type(self) ) and self.__dict__ == other.__dict__

    def __ne__(self, other):
        """Overridden inequality method. Returns true if the meshes are not identical.

        This method is called by the built in python operator "!="
        
        Parameters
        ----------
        
        other : any type
            object that we are comparing to

        Returns
        -------
        
        meshes_are_not_identical : bool
            True if meshes are not identical
        """
        return not self.__eq__(other)
    
    def save(self, filename):
        """Save the mesh to a file.
        
        Parameters
        ----------
        
        filename : string
            the name under which you would like to save the mesh
        """
        current_recursion_limit = sys.getrecursionlimit()
        sys.setrecursionlimit(40000)
        file_to_write = open(filename, 'wb')
        pickle.dump(self, file_to_write)
        file_to_write.close()
        sys.setrecursionlimit(current_recursion_limit)
        
    def find_most_central_element(self):
        """Finds the most central element in the mesh.
        
        The most central element in the mesh is the element whose
        centroid is closest to the centre of the mesh
        
        Returns
        -------
        
        most_central_element : Element instance
            The element of the mesh whose centre is closest to the centre
            of the mesh
        """

        mesh_centre = self.calculate_centre()

        min_distance = 3*(self.calculate_height() + self.calculate_width())

        for element in self.elements:
            distance = np.linalg.norm( element.calculate_centroid() - mesh_centre )
            if distance < min_distance:
                min_distance = distance
                most_central_element = element 
               
        return most_central_element

    def find_most_central_node(self):
        """Finds the most central node in the mesh.
        
        The most central node in the mesh is the node whose
        position is closest to the centre of the mesh
        
        Returns
        -------
        
        most_central_node : Node instance
            The node of the mesh whose position is closest to the centre
            of the mesh
        """
        mesh_centre = self.calculate_centre()

        min_distance = 3*(self.calculate_height() + self.calculate_width())

        for node in self.nodes:
            distance = np.linalg.norm(node.position - mesh_centre)
            if distance < min_distance:
                min_distance = distance
                most_central_node = node 
               
        return most_central_node
    
    def remove_boundary_elements(self):
        """Removes all elements at the outside of the mesh
        
        All elements that have an edge that is not shared with other elements
        are removed.
        """
        
        frame_ids_to_delete = []
        for element in self.elements:
            if element.check_if_on_boundary():
                frame_ids_to_delete.append( element.id_in_frame )
                
        for frame_id in frame_ids_to_delete:
            self.delete_element_with_frame_id(frame_id)
    
    def delete_element_with_frame_id(self, frame_id):
        """Removes the element from the mesh entirely.
        
        Deletes the element from all shared nodes, remove all nodes belong to the element only.
        Deletes the element from element list, keeps internal frame id indexing updated.
        
        Parameters
        ----------
        
        frame_id : int
            id_in_frame of the element that is to be deleted
            
        See Also
        --------
        
        remove_element_with_frame_id : remove an alement without effecting the nodes
        """
        
        element_to_delete = self.get_element_with_frame_id( frame_id )
        
        nodes_to_remove = []

        for node in element_to_delete.nodes:
            if len(node.adjacent_elements) == 1:
                nodes_to_remove.append(node)
            else:
                node.remove_element_with_frame_id( frame_id )
        
        self.remove_list_of_nodes(nodes_to_remove)
        
        self.remove_element_with_frame_id( frame_id )
        
    def merge_short_edges(self, threshold_distance):
        """Merge short edges
        
        Merge short edges into single nodes. If multiple short edges are connected,
        merge them all into the same node.
        
        Parameters
        ----------
        
        threshold_distance : double
            the distance below which edges should be merged
        """

        # collect all edges that are shorter than the threshold distance
        edges = self.collect_edges()
        edges_to_merge = []
        for edge_index, edge in enumerate(edges):
            edge_length = np.linalg.norm( self.get_node_with_id( edge[0] ).position -
                                          self.get_node_with_id( edge[1] ).position )
            if edge_length < threshold_distance:
                edges_to_merge.append(edge)
        edges_np = np.array(edges_to_merge)
        
        # arrange short edges into clusters (pointslists)
        pointslists_to_merge = []
        edge_already_clustered = np.zeros(len(edges_to_merge),dtype = 'bool')

        for edge_index, edge in enumerate(edges_to_merge):
            if not edge_already_clustered[edge_index]:
                node_indices_to_merge = list(edge)
                still_extendable = True
                indices_to_keep = []
                while still_extendable:
                    still_extendable = False
                    for node_index in node_indices_to_merge:
                        occurences = np.where( edges_np == node_index )
                        edges_in_mersion = occurences[0]
                        new_node_indices_to_merge = np.unique(edges_np[edges_in_mersion])
                        for new_node_index in new_node_indices_to_merge:
                            if new_node_index not in node_indices_to_merge:
                                indices_to_keep.append(new_node_index)
                                still_extendable = True
                        edge_already_clustered[edges_in_mersion] = True
                    node_indices_to_merge += indices_to_keep
                pointslists_to_merge.append(node_indices_to_merge)

        # for each cluster, merge all nodes to their average position
        for pointslist in pointslists_to_merge:

            list_of_positions = []
            for point in pointslist:
                list_of_positions.append( self.get_node_with_id(point).position )
            new_node_position = np.mean(list_of_positions, axis = 0)

            # make new node
            maximal_node_id = self.get_maximal_node_id()
            new_node = Node(new_node_position, maximal_node_id + 1)
            self.nodes.append(new_node)

            # for each adjacent element of these nodes remove these nodes
            # from the element and add the new node instead
            adjacent_element_ids = []
            for node_id in pointslist:
                this_node = self.get_node_with_id(node_id)
                for element in this_node.adjacent_elements:
                    adjacent_element_ids.append(element.id_in_frame)
                    new_element_nodes = []
                    node_already_appended = False
                    for node in element.nodes:
                        if node.id != ( maximal_node_id + 1 ):
                            if node.id not in pointslist:
                                new_element_nodes.append(node)
                            else:
                                if not node_already_appended: 
                                    new_element_nodes.append(new_node)
                                    node_already_appended = True
                        else:
                            if not node_already_appended: 
                                new_element_nodes.append(new_node)
                                node_already_appended = True
                    element.nodes = new_element_nodes
            
            unique_adjacent_element_ids = np.unique(adjacent_element_ids)
            
            for element_id in unique_adjacent_element_ids:
                this_element = self.get_element_with_frame_id(element_id)
                new_node.adjacent_elements.append(this_element)
                
            self.remove_list_of_nodes( [self.get_node_with_id(node_id) for node_id in pointslist] )
            
    def plot_with_data(self, filename, real_frame_path):
        """Overlays the mesh with the actual data from the image
        
        Parameters
        ----------
        
        filename : string
            path to the file where the overlay should be saved
            
        real_frame_matrix : string
            path to the image that is associated with this frame.
        """
        # figure properties
        figuresize = (4,2.75)
        font = {'size'   : 10}
        plt.rc('font', **font)

        real_image = plt.imread(real_frame_path)

        rosette_list = []
        node_positions_list = []
        for node in self.nodes:
            node_positions_list.append(node.position)
            if len(node.adjacent_elements) > 4:
                rosette_list.append(node.position)

      
        node_positions = np.array(node_positions_list)
        rosette_positions = np.array(rosette_list)
        
        node_positions[:,1] *= -1
        node_positions[:,1] += len(real_image[:,0])
        try:
            rosette_positions[:,1] *= -1
            rosette_positions[:,1] += len(real_image[:,0])
        except:
            pass
        #transform vertex coordinates to image coordinates
        polygon_collection = self.get_polygon_collection()
        for path in polygon_collection.properties()['paths']:
            new_vertices = np.zeros_like(path.vertices)
            path.vertices[:,1] *= -1
            path.vertices[:,1] += len(real_image[:,0])

        polygon_collection.set_facecolor('None')
        polygon_collection.set_edgecolor('red')
        polygon_collection.set_alpha(0.25)
        polygon_collection.set_linewidth(3)
        
        overlay_figure = plt.figure()
        plt.imshow(real_image, cmap = plt.cm.Greys_r)
        overlay_figure.gca().add_collection(polygon_collection)
        plt.scatter(node_positions[:,0], node_positions[:,1], s = 1, color = 'black', alpha = 0.3, marker = 'o', linewidths = 0.0)
        try:
            plt.scatter(rosette_positions[:,0], rosette_positions[:,1], s = 100, color = 'black', alpha = 0.3, marker = 'o', linewidths = 0.0)
        except:
            pass
        overlay_figure.savefig(filename, bbox_inches = 'tight')
        
    def plot_tracked_data(self, filename, image_path, segmented_path, max_global_id):
        """Overlays the mesh with the actual data from the image
        
        Parameters
        ----------
        
        filename : string
            path to the file where the overlay should be saved
            
        real_frame_matrix : string
            path to the image that is associated with this frame.
        """
        # figure properties
        figuresize = (4,2.75)
        font = {'size'   : 10}
        plt.rc('font', **font)

        real_image = plt.imread(image_path)
        segmented_image = cv2.imread( segmented_path, flags = -1 )
        
        color_selection = np.array( _get_distinct_colors( max_global_id + 1 ) )*255
        np.random.seed(12)
        np.random.shuffle(color_selection)
        
        real_image_converted = np.zeros( ( real_image.shape[0], real_image.shape[1], 3 ),
                                           dtype = 'uint8' )
        
        real_image_converted[:,:,0] = real_image*255.0/65535.0
        real_image_converted[:,:,1] = real_image*255.0/65535.0
        real_image_converted[:,:,2] = real_image*255.0/65535.0
       
        overlay_image = np.zeros( ( real_image.shape[0], real_image.shape[1], 3 ),
                                    dtype = 'uint8' )

        overlay_image[:] = real_image_converted[:]
        helper_image = np.zeros( ( real_image.shape[0], real_image.shape[1], 3 ),
                                  dtype = 'uint8' )

        cell_counter = 0
        for element in self.elements:
            if element.global_id != None:
                cell_counter+=1
                this_color = color_selection[element.global_id]
                this_frame_id = element.id_in_frame
                this_mask = segmented_image == this_frame_id
                helper_image[this_mask] = this_color
                overlay_image[this_mask] = real_image_converted[this_mask]*0.5 + helper_image[this_mask]*0.5
                if getattr(element, 'is_new', False):
                    try:
                        centroid = element.calculate_centroid()
                        centroid_image_coordinates = centroid.copy()
#                         centroid_image_coordinates[0] = real_image.shape[0] - centroid[1]
                        centroid_image_coordinates[1] = real_image.shape[0] - centroid[1]
                        centroid_image_coordinates[0] = centroid[0]
#                         import pdb; pdb.set_trace()
#                         import pdb; pdb.set_trace()
                        cv2.circle(overlay_image, tuple( centroid_image_coordinates.astype('int') ),
                                   radius = 2, color = (0,0,255), thickness = -1)
                    except:
                        pass
                
                try:
                    centroid = element.calculate_centroid()
                    centroid_image_coordinates = centroid.copy()
                    centroid_image_coordinates[1] = real_image.shape[0] - centroid[1]
                    centroid_image_coordinates[0] = centroid[0]
                    centroid_image_coordinates += [-5,5]
                    cv2.putText(overlay_image, str(element.global_id), tuple(centroid_image_coordinates.astype('int')),
                                cv2.FONT_HERSHEY_SCRIPT_SIMPLEX, fontScale = 0.2, 
                                color = (0,0,0) )
#                     cv2.putText( overlay_image, str(element.id_in_frame), tuple(centroid_image_coordinates.astype('int')),
#                                  cv2.FONT_HERSHEY_SCRIPT_SIMPLEX, fontScale = 0.2, 
#                                  color = (0,0,0) )
                except:
                    pass
#                 overlay_image[this_mask] = real_image_converted[this_mask]
#                 import pdb; pdb.set_trace();
                helper_image[:] = 0
                
        cv2.imwrite( filename, overlay_image )

    def get_node_with_id(self, node_id):
        """Get the node instance that corresponds to this node id.
        
        Parameters
        ----------
        
        node_id : int
        
        Returns
        -------
        
        node_instance : node with id node_id
        """
        
        for node in self.nodes:
            if node.id == node_id:
                node_instance = node
                break
            
        return node_instance
    
    def get_dangling_element_ids(self, list_of_mapped_elements):
        """Returns the global ids of all `dangling' elements.
        
        An element is considered `dangling' if it is in list_of_mapped_elements and if
        it has only one neighbour within the list, or, if it has two neighbours within the list
        but one of them has the same neighbour number as the element.
        
        Parameters
        ----------
        
        list_of_mapped_elements : list of ints
            frame ids of all elements that are considered to be mapped

        Returns
        -------
        
        dangling_elements : list of ints
            entries in the list are global ids of elements that are `dangling'
        """

        dangling_elements_counter = 0
        dangling_elements = []
        for element in self.elements:
            # find how many 'mapped' neighbours the element has
            if element.id_in_frame in list_of_mapped_elements:
                ids_adjacent_to_element = element.get_ids_of_adjacent_elements()
                reduced_ids_adjacent_to_element = []
                for frame_id in ids_adjacent_to_element:
                    if frame_id in list_of_mapped_elements:
                        reduced_ids_adjacent_to_element.append(frame_id)
            # if it's only one such neighbour, then it's obviosly dangling
                if (len(reduced_ids_adjacent_to_element) == 1 ):
                    dangling_elements_counter += 1
                    dangling_elements.append(element.id_in_frame)
                # if it's two, then check the neighbours for their
                # polygon number, and their number of neighbours
                elif (len(reduced_ids_adjacent_to_element) == 2):
                    neighbouring_polygon_number_is_the_same = False
                    for frame_id in reduced_ids_adjacent_to_element:
                        neighbouring_element = self.get_element_with_frame_id(frame_id)
                        element_polygon_number = neighbouring_element.get_num_nodes()
                        if element_polygon_number == element.get_num_nodes():
                            ids_adjacent_to_neighbouring_element = element.get_ids_of_adjacent_elements()
                            reduced_ids_adjacent_to_neighbouring_element = []
                            for adjacent_frame_id in ids_adjacent_to_neighbouring_element:
                                if adjacent_frame_id in list_of_mapped_elements:
                                    reduced_ids_adjacent_to_neighbouring_element.append(adjacent_frame_id)
                            if len(reduced_ids_adjacent_to_neighbouring_element) <= 2:
                                dangling_elements.append(element.id_in_frame)
                                dangling_elements_counter += 1
                                break
        
        return dangling_elements
    
class Element():
    """Elements are members of a mesh.
    """
    def __init__(self,nodes, id_in_frame = None):
        """The element generator.
            
        Parameters
        ----------
        nodes : list
            Each entry is a node    
        id_in_frame : int
            id of this element, defaults to None
            
        Returns
        -------
        element : the element
        """
        self.nodes = nodes
        """ A list containing all the nodes of this element"""
        
        self.id_in_frame = id_in_frame
        """ The id of this element within this particular frame"""
        
        self.global_id = None
        """ The id of this element globally"""
        
        # We add this element to all its nodes
        self.add_element_to_nodes()
        
    def add_element_to_nodes(self):
        """Adds the element to all its nodes, make sure to not do this twice!"""

        for node in self.nodes:
            node.adjacent_elements.append(self)

    def get_num_nodes(self):
        """Returns the number of nodes that this element shares.
        
        Returns
        -------
        
        number_of_shared_nodes : int
            number of nodes that are members of this element
        """

        return len(self.nodes)
    
    def calculate_area(self):
        """A simple method to return the area of the element
        
        Returns
        -------
        
        area : double
            area of the element
        """
        area = 0.0
        for index, node in enumerate(self.nodes):
            next_index = (index + 1)%self.get_num_nodes()
            next_node = self.nodes[next_index]
            area += (node.position[0] * next_node.position[1]) - (next_node.position[0] * node.position[1])
        area = area/2.
        return area
    
    def calculate_centroid(self):
        """Calculate the centroid of the vertex element
        
        Returns
        -------
        
        centroid :  1x2 np_array
            position of the centroid of the element
        """
        centroid = np.zeros(2)
        for index, node in enumerate(self.nodes):
            next_index = (index + 1)%self.get_num_nodes()
            next_node = self.nodes[next_index]
            centroid[0] += (node.position[0] + next_node.position[0]) * ((node.position[0] * next_node.position[1]) - 
                                                                         (next_node.position[0] * node.position[1]))
            centroid[1] += (node.position[1] + next_node.position[1]) * ((node.position[0] * next_node.position[1]) - 
                                                                         (next_node.position[0] * node.position[1]))
        area = self.calculate_area()
        assert(area > 0)
        centroid = centroid/(6.0*area)
        assert(len(centroid) == 2)
        return centroid
    
    def replace_adjacent_element_with_node(self, frame_id, node):
        """Finds the edge to the element with given id, and replaces the edge with the given node.
        
        Throws an error if the two elements share a rosette. Does not alter the adjacent element.
        Adds this element (self) to the node. Designed as helper function for mesh.kill_element_with_frame_id().
        
        Parameters
        ----------
        
        frame_id : int
            id_in_frame of the element that shares the edge which we would like to replace.
            
        node : Node instance
            the node that replaces the shared edge.
            
        See also
        --------
        
        Mesh.kill_element_with_frame_id()
        """
        
        # Identify the local indices of nodes on the edge
        local_indices_on_edge = []
        for local_index, element_node in enumerate(self.nodes):
            if frame_id in element_node.get_adjacent_element_ids():
                local_indices_on_edge.append(local_index)
                
        if not (len(local_indices_on_edge) == 2):
            # replace this warning with an exception if you don't want this to happen
            warnings.warn("Merging rosette or double edge. This is not fully tested.")

        new_nodes = []
        for local_index, element_node in enumerate(self.nodes):
            if local_index not in local_indices_on_edge:
                new_nodes.append(element_node)
            elif local_index == local_indices_on_edge[0]:
                new_nodes.append(node)
                
        self.nodes = new_nodes
        node.adjacent_elements.append(self)
        
#     def get_ids_of_adjacent_elements(self):
#         """Returns a list of frame_ids of all adjacent elements.
#           
#         Assembles the list by collecting adjacent elements of contained nodes.
#           
#         Returns
#         -------
#           
#         adjacent_element_ids : np.array with integer values
#             all ids that are adjacent to this element
#           
#         See also
#         --------
#           
#           
#         Mesh.kill_element_with_frame_id() : uses this function
#         """
#           
#         # the following list is going to include double entries
#         all_adjacent_element_ids = []
#         for node in self.nodes:
#             all_adjacent_element_ids += node.get_adjacent_element_ids()
#               
#         unique_adjacent_element_ids = np.unique(np.array(all_adjacent_element_ids))
#           
#         # remove self from the element ids
#         unique_adjacent_element_ids = unique_adjacent_element_ids[unique_adjacent_element_ids != self.id_in_frame]
#           
#         return unique_adjacent_element_ids
    
    def get_ids_of_adjacent_elements(self):
        """Returns a list of frame_ids of all adjacent elements.
           
        Elements are adjacent with this element if they share an edge.
           
        Returns
        -------
           
        adjacent_element_ids : np.array with integer values
            all ids that are adjacent to this element
           
        See also
        --------
           
        Mesh.kill_element_with_frame_id() : uses this function
        """

        all_adjacent_element_ids = []
   
        for node_index, node in enumerate( self.nodes ):
            next_node_index = (node_index + 1)%self.get_num_nodes()
            next_node = self.nodes[next_node_index]
            shared_elements = set.intersection( set(node.get_adjacent_element_ids() ), set(next_node.get_adjacent_element_ids()) )
            shared_elements.remove( self.id_in_frame )
            if len( shared_elements ) > 0:
                assert(len(shared_elements) == 1)
                all_adjacent_element_ids.append(shared_elements.pop())
               
        return all_adjacent_element_ids

    def __eq__(self, other):
        """Overridden equality method. Returns true if the elements are identical.

        This method is called by the built in python operator "=="
        
        Parameters
        ----------
        
        other : any type
            object that we are comparing to

        Returns
        -------
        
        meshes_are_identical : bool
            True if elements are identical
        """
        return ( type(other) is type(self) ) and self.__dict__ == other.__dict__

    def __ne__(self, other):
        """Overridden inequality method. Returns true if the elements are not identical.

        This method is called by the built in python operator "!="
        
        Parameters
        ----------
        
        other : any type
            object that we are comparing to

        Returns
        -------
        
        elements_are_not_identical : bool
            True if elements are not identical
        """
        return not self.__eq__(other)

    def check_if_on_boundary(self):
        """Return true if this element is on the boundary of the mesh.
        
        The element is on the boundary of the mesh if it has an edge that it
        does not share with another element. The check is done `in place' and
        this variable is not indexed.
        
        Returns
        -------
        
        element_is_on_boundary : bool
            True if element is on boundary, False if otherwise.
        """
        
        element_is_on_boundary = False
        
        if self.get_num_nodes() <= 2:
            element_is_on_boundary = True

        for local_index, this_node in enumerate(self.nodes):
            next_local_index = (local_index + 1)%self.get_num_nodes()
            next_node = self.nodes[next_local_index]
            frame_ids_this_node = this_node.get_adjacent_element_ids()
            frame_ids_next_node = next_node.get_adjacent_element_ids()
            shared_frame_ids = set.intersection( set(frame_ids_this_node), set(frame_ids_next_node) )
            if len(shared_frame_ids) == 1:
                element_is_on_boundary = True
                break
        
        return element_is_on_boundary
                
class Node():
    """Nodes are points in a mesh.
    """
    def __init__(self, position, id = None):
        """The node generator.
        
        Parameters
        ----------
        
        position : double array like
            position of the Node
        
        id : int
            id of the node, defaults to None
        
        Returns
        -------
        
        node : the node
        
        Warnings
        --------
        
        The equality method is hardcoded. That means, if you change or add node members you will need to 
        also implement that their equality is checked in the __eq__() method.
        """

        self.position = np.array(position, dtype = 'double')
        """The position of the node"""

        self.id = id
        """The id of the node"""

        self.adjacent_elements = []
        """A list with all the elements this node belongs to"""
        
    def get_adjacent_element_ids(self):
        """returns a list of ids of all containing elements
        
        Returns
        -------
        
        id_list : int list
            list containing ids of all adjacent elements
        """
        
        id_list = []

        for element in self.adjacent_elements:
            id_list.append(element.id_in_frame)
        
        return id_list

    def remove_element_with_frame_id(self, id_in_frame ):
        """Remove element with id_in_frame from list of adjacent elements.
        
        Parameters
        ----------
        
        id_in_frame : int
            id of element that is to be removed from the list of adjacent elements
            for this node
        """
        self.adjacent_elements = [element for element in self.adjacent_elements if 
                                  element.id_in_frame != id_in_frame]
        
    def __eq__(self, other):
        """Overridden equality method. Returns true if the Nodes are identical.
 
        This method is called by the built in python operator "==". Identity for all
        member variables is hardcoded. Adjacent elements are compared by id only.
         
        Parameters
        ----------
         
        other : any type
            object that we are comparing to
 
        Returns
        -------
         
        nodes_are_identical : bool
            True if meshes are identical
        """
        if type(other) is type(self):
            position_is_equal = np.array_equal(self.position, other.position) 
            id_is_equal = ( self.id == other.id )
            element_ids_are_equal = ( self.get_adjacent_element_ids() == 
                                      other.get_adjacent_element_ids() )
            
            nodes_are_identical = ( position_is_equal and id_is_equal and element_ids_are_equal )
        else:
            nodes_are_identical = False

        return nodes_are_identical
 
    def __ne__(self, other):
        """Overridden inequality method. Returns true if the nodes are not identical.
 
        This method is called by the built in python operator "!="
         
        Parameters
        ----------
         
        other : any type
            object that we are comparing to
 
        Returns
        -------
         
        nodes_are_not_identical : bool
            True if nodes are not identical
            
        See also
        --------
        
        self.__eq__() equality method
        """
        return not self.__eq__(other)
    
def _get_distinct_colors(number_of_colors):
    """Returns number_of_colors different colors
    
    Parameters
    ----------
    
    number_of_colors : int
    
    Returns
    -------
    
    colormap : list
        a list where each entry is a different color
    """
    HSV_tuples = [(color_index*1.0/number_of_colors, 1.0, 0.7) for color_index in range(number_of_colors)]
    RGB_tuples = map(lambda color: colorsys.hsv_to_rgb(*color), HSV_tuples)

    return RGB_tuples
