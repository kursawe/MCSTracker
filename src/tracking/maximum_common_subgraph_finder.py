# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""This file contains the algorithm for finding maximum common subgraphs
"""
import numpy as np
import networkx as nx
import itertools
import time
import math
   
class KrissinelMaximumCommonSubgraphFinder:
    """A helper class for finding maximum common subgraphs
    
    The class implements the algorithm described in
    
    Evgeny B. Krissinel, Kim Henrick. Common subgraph isomorphism detection
    by backtracking search. Softw. Pract. Exper. 34, pp 591-607 (2007)
    DOI: 10.1002/spe.588
    """
    def __init__(self, mesh_one, mesh_two, minimum_size = 0.7):
        """Initialise algorithm with two Mesh-type objects
        
        Parameters
        ----------
        
        mesh_one: Mesh instance
            One of the two meshes to which we would like to apply the maximum
            common subgraph algorithm
            
        mesh_two: Mesh instance
            second of the two meshes
            
        minimum_size : double
            fraction of the number of cells in the first network that we expect to 
            be matched by the algorithm, has to be between zero and one. No smaller
            mappings will be considered
            
        The constructor will setup the class but not run the algorithm. Run the algorithm
        by calling find_maximum_common_subgraph
        """
        
        self.mesh_one = mesh_one
        """First mesh"""
        
        self.mesh_two = mesh_two
        """Second mesh"""
        
        self.network_one = mesh_one.generate_network()
        """Network of first mesh"""

        self.network_two = mesh_two.generate_network()
        """Network of first mesh"""

        self.width_one = mesh_one.calculate_width()
        """Width of the first mesh"""

        self.height_one = mesh_one.calculate_height()
        """height of the first mesh"""
        
        self.min_subgraph_size = int( minimum_size*len( self.network_one.nodes() ) )
        """Minimum number of matches that we expect for the maximum common subgraph"""

        self.max_found_subgraph_size = 0
        """maximum size of currently found subgraph"""
        
        self.largest_mappings = []
        """list of all tracking states that have the current maximal number of edges"""

        self.branches_counter = 0
        """keeps track of the number of branches we have analysed"""
        
        self.initial_tracking_state = self.create_initial_tracking_state()
        """The initial tracking state.
        
        See also
        --------
        
        create_initial_tracking_state
        TrackingState
        """

        self.timing_is_on = False
        """If this is true we turn on the collection of time statistics"""
        
        self.start_time = None
        """Start time of the search, for reference."""

        self.total_execution_time = None
        """Total time of the search"""
        
        self.network_one_node_iterator = self.network_one.nodes()
        self.network_two_node_iterator = self.network_two.nodes()
        
    def turn_timing_on(self):
        """turn on the collection of timing statistics"""

        self.timing_is_on = True
        
    def create_initial_tracking_state(self):
        """Creates an instance of TrackingState. In this instance all potential matchings
        between vertices in the first and second mesh are identified.
        
        Returns
        -------
        
        initial_tracking_state : TrackingState instance
            The id_map of this state will be empty, the vertex_matching_matrix will be filled
            with all possible matches, and the counts_of_mappable_vertices will be correspondingly initialised.
        
        See also
        --------
        
        TrackingState
        """
        
        # make a tracking state
        tracking_state = TrackingState()

        # initialize lists that will later form counts_of_mappable_vertices
        # and vertex_matching_matrix
        vertex_matches = []
        numbers_of_matches = []

        cutoff_distance = self.determine_cutoff_distance()

        # loop over all nodes in network one
        for node_one in self.network_one.nodes():

            matches_for_node_one = []

            # identify which nodes in network two could be corresponding
            for node_two_index, node_two in enumerate(self.network_two.nodes()):
                node_two_is_mappable_to_node_one = ( np.linalg.norm(self.network_one.node[node_one]['position'] -
                                                                    self.network_two.node[node_two]['position'] ) < cutoff_distance 
                                                     and self.network_one.node[node_one]['num_neighbours'] == 
                                                     self.network_two.node[node_two]['num_neighbours'])

                # if a node is mappable, remember that
                if node_two_is_mappable_to_node_one:
                    matches_for_node_one.append(node_two_index)
            
            # sort the matches by their distance to the original node
            matches_for_node_one.sort(key = lambda x: np.linalg.norm(self.network_two.node[self.network_two.nodes()[x]]['position']-
                                                                     self.network_one.node[node_one]['position']))

            # put what we learned in the lists
            vertex_matches.append(matches_for_node_one)
            numbers_of_matches.append(len(matches_for_node_one))
        
        # This allows us to initialise the counts_of_mappable_vertices
        # np.arrays are fast, hopefully
        tracking_state.counts_of_mappable_vertices = np.array(numbers_of_matches, dtype = 'int')

        # Find out how many columns we maximally need for the vertex matching matrix - it's size won't change later on
        maximal_number_of_matches = np.max(tracking_state.counts_of_mappable_vertices)

        # the vertex_matching_matrix has as many rows as the first network has nodes, maximal_number_of_matches columns 

        tracking_state.vertex_matching_matrix = np.zeros( ( self.network_one.number_of_nodes(), maximal_number_of_matches ), dtype = 'int')
        
        # The number of connections to the maximum common subgraph are initialised here but only used in child classes
        
        tracking_state.connections_to_current_subgraph = np.zeros( self.network_one.number_of_nodes(), dtype = 'int')
        
        # fill the entries of the vertex_matching_matrix
        for node_one_index, matching_vertices in enumerate(vertex_matches):
            for index_counter, node_two_index in enumerate(matching_vertices):
                tracking_state.vertex_matching_matrix[node_one_index, index_counter] = node_two_index
                
        return tracking_state
    
    def determine_cutoff_distance(self):
        """Estimate the maximal physical distance that we expect a cell and its image to be away from each other.
        
        This distance is currently 10 times the average cell distance in mesh_one                
        
        Returns
        -------
        
        cutoff_distance : float
            theshold for when cells are mappable to each other
        """
        
        cutoff_distance = 10.0*math.sqrt(self.mesh_one.calculate_average_element_area())
        
        return cutoff_distance
        
    def find_maximum_common_subgraph(self):
        """Runs the backtracking algorithm to find the maximum
        common subgraph
        
        The maximum common subgraph is found by recursive application
        of the backtrack() function to instances_of_tracking_state. 
        a tracking state keeps track of the current subgraph, and the remaining
        possible matchings between vertices.
        
        See also
        --------
        
        backtrack
        """

        if self.timing_is_on:
            self.start_time = time.clock()
        
        self.backtrack(self.initial_tracking_state)

        if self.timing_is_on:
            end_time = time.clock()
            self.total_execution_time = end_time - self.start_time

    def backtrack(self, tracking_state):
        """The recursive backtrack function
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance 
            The tracking state keeps track of the subgraph that is currently analyzed. And the
            possible matchings that remain to be analyzed.
        
        This function will analyze whether the current search branch can be expanded,
        and if not then check whether the branch is the largest currently found common subgraph.
        If the branch is expandable it will pick a new vertex in the first network, find all possible
        matches for this vertex, and set up a new branch for each of these matches - on which backtrack 
        will be called again. Since the picked vertex does not have to be part of the maximal common subgraph
        it will then set up a branch where this vertex is not in present in the first network, and attempt to call 
        backtrack on this branch.
        """

        if self.is_extendable(tracking_state):
            current_vertex = self.pick_next_vertex(tracking_state)
            possible_images = self.get_mappable_vertices(tracking_state, current_vertex)
        
            for image in possible_images:
                # this method will also update the vertex_matching_matrix in the tracking state
                new_tracking_state = self.create_extended_tracking_state(tracking_state, current_vertex, image)
                self.backtrack(new_tracking_state)
        
            # take the current vertex out of the possible extensions for the current tracking state
            if tracking_state.inverse_sparse_matrix_lookup == None:
                tracking_state.vertex_matching_matrix[current_vertex,:] = 0
                tracking_state.counts_of_mappable_vertices[current_vertex] = 0
            else:
                reduced_current_vertex = tracking_state.inverse_sparse_matrix_lookup[current_vertex]
                tracking_state.vertex_matching_matrix[reduced_current_vertex,:] = 0
                tracking_state.counts_of_mappable_vertices[reduced_current_vertex] = 0
            
            # and run backtrack on this
            self.backtrack(tracking_state)
            
        else:
            self.branches_counter += 1
            if len( tracking_state.id_map.keys() ) > self.max_found_subgraph_size:
                self.max_found_subgraph_size = len( tracking_state.id_map.keys() )
                self.largest_mappings = [tracking_state.id_map]
            elif len( tracking_state.id_map.keys() ) == self.max_found_subgraph_size:
                self.largest_mappings.append( tracking_state.id_map )
                
    def create_extended_tracking_state(self, tracking_state, vertex_one, vertex_two):
        """Create an new tracking state, where tracking state is extended by the mapping from current_vertex 
        to image_vertex, and the possible matches are refined
        
        This method corresponds to the Refine() method in Krissinel et al.
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            the tracking state that we are planning on copying and extending
            
        vertex_one : int
            index of a node in the network_one.nodes vector
            
        vertex_two : int
            index of a node in the network_two.nodes vector
            
        The method extends tracking_state by the mapping between network_one.nodes[current_vertex]->network_two.nodes[vertex_two],
        and refines the vertex_matching_matrix, such that only mappings are considered that have no edge conflicts with the current nodes.
        """
        
        # np.nonzero will return a tuple with one entry, we are interested in the entry, which is an array with all
        # the indices where the mappable vertices count is not zero. Hence, we take the zero'th element of the nonzero return value
        list_of_still_mappable_vertices_in_old_state = np.nonzero( tracking_state.counts_of_mappable_vertices )[0]
        
        new_tracking_state = TrackingState()
        
        new_tracking_state.id_map = tracking_state.id_map.copy()
        new_tracking_state.id_map[self.network_one_node_iterator[vertex_one]] = self.network_two_node_iterator[vertex_two]
        
        new_tracking_state.vertex_matching_matrix = np.zeros_like( tracking_state.vertex_matching_matrix, dtype = 'int')
        new_tracking_state.counts_of_mappable_vertices = np.zeros_like( tracking_state.counts_of_mappable_vertices, dtype = 'int')
        
        new_tracking_state.sparse_matrix_lookup = tracking_state.sparse_matrix_lookup
        new_tracking_state.inverse_sparse_matrix_lookup = tracking_state.inverse_sparse_matrix_lookup
        
        for reduced_vertex_index in list_of_still_mappable_vertices_in_old_state:
            if tracking_state.inverse_sparse_matrix_lookup == None:
                vertex_index = reduced_vertex_index
            else:
                vertex_index = tracking_state.sparse_matrix_lookup[reduced_vertex_index]
            if vertex_index != vertex_one:
                count_of_matching_vertices = 0
                for column_index_in_vertex_matching_matrix in range(tracking_state.counts_of_mappable_vertices[reduced_vertex_index]):
                    image_index = tracking_state.vertex_matching_matrix[reduced_vertex_index, column_index_in_vertex_matching_matrix]
                    if image_index != vertex_two:
                        # check whether this potential new pairing would have a conflicting edge with the pairing
                        # that we just added
                        if self.pairing_does_not_introduce_edge_conflicts((vertex_one, vertex_index), (vertex_two, image_index)):
                            new_tracking_state.vertex_matching_matrix[reduced_vertex_index][count_of_matching_vertices] = image_index
                            count_of_matching_vertices += 1
                new_tracking_state.counts_of_mappable_vertices[reduced_vertex_index] = count_of_matching_vertices
                
        return new_tracking_state
                        
    def pairing_does_not_introduce_edge_conflicts(self, network_one_tuple, network_two_tuple):
        """Returns False if the tuples have an edge in one of the networks but not in the other
        
        Parameters
        ----------
        
        network_one_tuple : ( int, int )
            a tuple of indices of network_one.nodes

        network_two_tuple : ( int, int )
            a tuple of indices of network_two.nodes
            
        Returns
        -------
        
        pairing_has_no_conflicting_edge : bool
            Is False if the tuple in one of the networks has an edge but the tuple in the other network doesn't.
            Is True otherwise
        """
        
        pair_one_has_edge = self.network_one.has_edge( self.network_one_node_iterator[network_one_tuple[0]], 
                                                       self.network_one_node_iterator[network_one_tuple[1]])

        pair_two_has_edge = self.network_two.has_edge( self.network_two_node_iterator[network_two_tuple[0]], 
                                                       self.network_two_node_iterator[network_two_tuple[1]])

        pairing_has_no_conflicting_edge = ( pair_one_has_edge == pair_two_has_edge )
        
        return pairing_has_no_conflicting_edge

    def is_extendable(self, tracking_state):
        """Check whether the current subgraph is extendable to a suitable maximum common subgraph
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            see backtrack() function
            
        Returns
        -------
        
        is_extendable : bool
            The current subgraph is extendable if the set of remaing possible matches is large enough, i.e.
            if it largest possible subgraph on the current branch is larger than our expectation, and than the current
            maximal common subgraph
        """
        number_vertices_still_mappable_in_network_one = np.sum(tracking_state.counts_of_mappable_vertices > 0)

        maximal_possible_size_of_this_subgraph = number_vertices_still_mappable_in_network_one + len(tracking_state.id_map)
        
        necessary_possible_size_of_this_subgraph_to_continue = max( self.min_subgraph_size, self.max_found_subgraph_size )
        
        is_extendable = ( ( number_vertices_still_mappable_in_network_one > 0 ) and 
                          ( maximal_possible_size_of_this_subgraph >= necessary_possible_size_of_this_subgraph_to_continue ) )

        return is_extendable
    
    def pick_next_vertex(self, tracking_state):
        """Pick the next vertex in network_one. Picks one of the vertices with the minimal number of possible matches
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            see backtrack() function, or TrackingState
            
        Returns
        -------
            
        vertex_index : int
            an index for an entry network_one.nodes
            
        Pick a node from network_one that hasn't been mapped yet. Pick one that has the minimum number of possible matches
        """
        minimum_number_of_mappable_vertices = np.min(tracking_state.counts_of_mappable_vertices[tracking_state.counts_of_mappable_vertices > 0])
        # np.nonzero will give us a tuple with the array that we are looking for as single entry. Let's get the array directly.
        possible_next_vertex_indices = np.nonzero( tracking_state.counts_of_mappable_vertices == minimum_number_of_mappable_vertices )[0]
        next_vertex_index = possible_next_vertex_indices[0]

        return next_vertex_index

    def get_mappable_vertices(self, tracking_state, vertex_index):
        """Get all vertices in network_two that could be images to the node in network one that is indexed by vertex_index.
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            see TrackingState or backtrack()

        vertex_index : int
            index for an entry in network_one.nodes
        
        Returns
        -------
        
        mappable_vertices : np.array, dtype = int
            a set of vertices in network_two that have the same degree as vertex,
            are within half the mesh width distance, and don't add an edge to the existing nodes
            that isn't there in network_one, while also not deleting or adding an edge that is or is not there
            in network_one
        """
        
        if tracking_state.inverse_sparse_matrix_lookup == None:
            mappable_vertices = tracking_state.vertex_matching_matrix[vertex_index,0:tracking_state.counts_of_mappable_vertices[vertex_index]]
        else:
            reduced_index = tracking_state.inverse_sparse_matrix_lookup[vertex_index]
            mappable_vertices = tracking_state.vertex_matching_matrix[reduced_index,0:tracking_state.counts_of_mappable_vertices[reduced_index]]

        return mappable_vertices

    def index_and_return_global_ids(self):
        """Write global id entries into the elements of mesh_one and mesh_two.
        
        Currently assigns global ids for matching cells in the first entry of the largest mappings.
        
        Not suitable if there are more than one largest mappings - use the PostProcessor class in this case.
        
        Returns
        -------
        
        mapped_ids : list of ints
            list of all global ids that have been assigned
        """
        
        assert( len(self.largest_mappings) == 1 )

        mapped_ids = []
        for global_id, frame_one_id in enumerate(self.largest_mappings[0]):
            self.mesh_one.get_element_with_frame_id(frame_one_id).global_id = global_id
            self.mesh_two.get_element_with_frame_id(self.largest_mappings[0][frame_one_id]).global_id = global_id
            mapped_ids.append(global_id)

        self.mesh_two.index_global_ids()
        self.mesh_one.index_global_ids() 
        
        return mapped_ids
    
class TrackingState:
    def __init__(self):
        """A helper class for maximum common subgraph finding. One instance of this class will describe
        one position in the search tree for maximum common subgraph finding. It keeps track of all the currently
        set mappings between the two networks, and all the possible remaining matches. This class is initialized
        in MaximumCommonSubgraphFinder.initialize_tracking_state(), modified in MaximumCommonSubgraphFinder.bactrack(), and
        copied and modified in MaximumCommonSubgraphFinder.create_extended_tracking_state()
        """

        self.id_map = {}
        """Keys of this dictionary are frame_ids in network_one, values are frame_ids in network two. This is the currently found
        subgraph_isomorphism"""

        self.vertex_matching_matrix = None
        """One row for each entry in network_one.nodes, has identical ordering. Entries are indices for network_two.nodes.
        Each row has entries for all possible matches in network_two.nodes that could correspond to the network_one.node that is identified
        by the row index of this matrix. This matrix is square and in most rows many of the entries will be unrelevant.
        Is zero in all redundant entries."""
         
        self.counts_of_mappable_vertices = None
        """Each entry in this vector contains the number of possible matches for a node in network one. Follows the same ordering
        as network_one.nodes. This vector is necessary to know which entries in vertex_matching_matrix we need to consider.
        """
        
        self.adjacent_vertices = None
        """List of elements that are adjacent to the current maximum common subgraph. Only used in ConnectedMaximumCommonSubgraphFinder,
        not KrissinalMaximumCommonSubgraphFinder"""
        
        self.connections_to_current_subgraph = None
        """One row for each entry in network one.nodes, has identical ordering. Each entry denotes how many connections this vertex
        has to the current subgraph"""

        self.not_blacklisted_vector = None
        """One row for each entry in network one.nodes, has identical ordering. Each entry denotes how many connections this vertex
        has to the current subgraph"""

        self.sparse_matrix_lookup = None
        """Used for reduced tracking states. Vector of indices self.network_one.nodes()"""
        
        self.inverse_sparse_matrix_lookup = None
        """Used for reduced tracking states. Dictionary of ints. Keys are indices of the self.network_one.nodes() vector, 
           values are indices for the reduced matrix lookups"""

class ConnectedMaximumCommonSubgraphFinder(KrissinelMaximumCommonSubgraphFinder):
    def __init__(self, mesh_one, mesh_two, minimum_size = 0.7):
        """A helper class to find maximum common subgraphs. Inherits from
        KrissinelMaximumCommonSubgraphFinder and adapts it to deal with connected subgraphs
        of planar topologies.
        
        Parameters for __init__ are the same as in the parent class.
        
        See also
        --------
        
        KrissinelMaximumCommonSubgraphFinder.__init__()
        """
        KrissinelMaximumCommonSubgraphFinder.__init__( self, mesh_one, mesh_two, minimum_size )

        self.initial_tracking_state.not_blacklisted_vector = np.ones(len(self.network_one.nodes()), dtype = 'bool')
        
        # create a lookup dict for indices in network one.
        network_one_node_dict = {}
        for counter, node in enumerate(self.network_one.nodes()):
            network_one_node_dict[counter] = node
        self.network_one_index_lookup = {node: index for index, node in network_one_node_dict.items()}
        """A dictionary to know which index belongs to a node in network_one.nodes()"""

        # create a lookup dict for indices in network one.
        network_two_node_dict = {}
        for counter, node in enumerate(self.network_two.nodes()):
            network_two_node_dict[counter] = node
        self.network_two_index_lookup = {node: index for index, node in network_two_node_dict.items()}
        """A dictionary to know which index belongs to a node in network_two.nodes()"""
        
    def pick_next_vertex(self, tracking_state, is_global = False):
        """Pick the next vertex in network_one. Picks one of the vertices with the minimal number of possible matches
        Overrides KrissinelMaximuCommonSubgraphFinder.pick_next_vertex().
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            see backtrack() function, or TrackingState
            
        is_global : bool
            if True, pick vertex for global scope (while loop) in LocalisedSubgraphFinder.
            This will avoid vertices which return false in tracking_state.not_blacklisted_vector.
            
        Returns
        -------
            
        vertex_index : int
            an index for an entry network_one.nodes
            
        Pick a node from network_one that hasn't been mapped yet. Pick one that has the minimum number of possible matches
        """
        if tracking_state.adjacent_vertices == None:
            next_vertex_index = KrissinelMaximumCommonSubgraphFinder.pick_next_vertex(self, tracking_state)
        else:
            if not is_global:
                mappable_mask = np.logical_and( tracking_state.counts_of_mappable_vertices > 0,
                                                tracking_state.connections_to_current_subgraph > 0)
            else: 
                mappable_mask = np.logical_and( np.logical_and( tracking_state.counts_of_mappable_vertices > 0,
                                                tracking_state.connections_to_current_subgraph > 0 ),
                                                tracking_state.not_blacklisted_vector )
                
            minimum_number_of_mappable_vertices = np.min(tracking_state.counts_of_mappable_vertices[mappable_mask])
            possible_vertices_with_minimal_matches = np.logical_and( mappable_mask,
                                                                     tracking_state.counts_of_mappable_vertices ==
                                                                     minimum_number_of_mappable_vertices )
            
            maximum_number_of_connections = np.max( tracking_state.connections_to_current_subgraph[possible_vertices_with_minimal_matches])

            possible_vertices_with_maximal_connection = np.logical_and ( possible_vertices_with_minimal_matches, 
                                                                         tracking_state.connections_to_current_subgraph ==
                                                                         maximum_number_of_connections )
            possible_next_vertex_indices = np.nonzero( possible_vertices_with_maximal_connection )
            
            if tracking_state.inverse_sparse_matrix_lookup == None:
                next_vertex_index = possible_next_vertex_indices[0][0]
            else:
                reduced_next_vertex_index = possible_next_vertex_indices[0][0]
                next_vertex_index = tracking_state.sparse_matrix_lookup[reduced_next_vertex_index]
            
        return next_vertex_index

    def is_extendable(self, tracking_state, is_global = False):
        """Check whether the current subgraph is extendable to a suitable maximum common subgraph
        Overrides and uses KrissinelMaximumCommonSubgraphFinder.is_extendable()
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            see backtrack() function
            
        Returns
        -------
        
        is_extendable : bool
            The current subgraph is extendable if the set of remaing possible matches is large enough, i.e.
            if it largest possible subgraph on the current branch is larger than our expectation, and than the current
            maximal common subgraph
        """
        if not is_global:
            if not np.sum( tracking_state.connections_to_current_subgraph ) == 0:
                mappable_mask = np.logical_and( tracking_state.counts_of_mappable_vertices > 0,
                                                tracking_state.connections_to_current_subgraph > 0)
            else:
                mappable_mask = ( tracking_state.counts_of_mappable_vertices > 0 )
        else: 
            mappable_mask = np.logical_and( np.logical_and( tracking_state.counts_of_mappable_vertices > 0,
                                            tracking_state.connections_to_current_subgraph > 0 ),
                                            tracking_state.not_blacklisted_vector )
 
        number_of_adjacent_vertices_still_mappable_in_network_one = np.sum( mappable_mask )

        is_extendable = ( KrissinelMaximumCommonSubgraphFinder.is_extendable(self, tracking_state) and
                          ( number_of_adjacent_vertices_still_mappable_in_network_one > 0 ) )

        return is_extendable

    def create_extended_tracking_state(self, tracking_state, vertex_one, vertex_two):
        """Create an new tracking state, where tracking state is extended by the mapping from current_vertex 
        to image_vertex, and the possible matches are refined
        
        Overrides and uses KrissinelMaximumCommonSubgraphFinder.create_extended_tracking_state()
        
        This method corresponds to the Refine() method in Krissinel et al.
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            the tracking state that we are planning on copying and extending
            
        vertex_one : int
            index of a node in the network_one.nodes vector
            
        vertex_two : int
            index of a node in the network_two.nodes vector
            
        The method extends tracking_state by the mapping between network_one.nodes[current_vertex]->network_two.nodes[vertex_two],
        and refines the vertex_matching_matrix, such that only mappings are considered that have no edge conflicts with the current nodes.
        """
        
        new_tracking_state = KrissinelMaximumCommonSubgraphFinder.create_extended_tracking_state(self, tracking_state, vertex_one, vertex_two)

        vertices_that_are_adjacent_to_vertex_one = self.network_one.neighbors(self.network_one.nodes()[vertex_one])
        indices_of_vertices_that_are_adjacent_to_vertex_one = np.array( [ self.network_one_index_lookup[vertex] 
                                                                          for vertex in vertices_that_are_adjacent_to_vertex_one])
        
        new_tracking_state.connections_to_current_subgraph = np.copy(tracking_state.connections_to_current_subgraph)
        if tracking_state.inverse_sparse_matrix_lookup == None:
            new_tracking_state.connections_to_current_subgraph[indices_of_vertices_that_are_adjacent_to_vertex_one] += 1
        else:
            for vertex_index in indices_of_vertices_that_are_adjacent_to_vertex_one:
                if vertex_index in tracking_state.inverse_sparse_matrix_lookup:
                    reduced_index = tracking_state.inverse_sparse_matrix_lookup[vertex_index] 
                    new_tracking_state.connections_to_current_subgraph[reduced_index] += 1

        new_tracking_state.not_blacklisted_vector = tracking_state.not_blacklisted_vector

        if (tracking_state.adjacent_vertices != None):
            new_tracking_state.adjacent_vertices = np.unique( np.hstack( ( tracking_state.adjacent_vertices,
                                                                           indices_of_vertices_that_are_adjacent_to_vertex_one ) ) )
        else:
            new_tracking_state.adjacent_vertices = indices_of_vertices_that_are_adjacent_to_vertex_one 
            
        return new_tracking_state

    
class ReducedBacktrackingSubgraphFinder(ConnectedMaximumCommonSubgraphFinder):
    def __init__(self, mesh_one, mesh_two, minimum_size = 0.65):
        """A helper class to find maximum common subgraphs. Inherits from
        ConnectedMaximumCommonSubgraphFinder and adapts it to deal with connected subgraphs
        of planar topologies.
        
        Parameters for __init__ are the same as in the parent class.
        
        See also
        --------
        
        ConnectedMaximumCommonSubgraphFinder.__init__()
        KrissinelMaximumCommonSubgraphFinder.__init__()
        """
        ConnectedMaximumCommonSubgraphFinder.__init__( self, mesh_one, mesh_two, minimum_size )

        self.safe_bets = {}
        """keys are indices in network_one.nodes(), and values are indices in network_two.nodes()"""

        self.preliminary_safe_bets = {}
        """keys are indices in network_one.nodes(), and values are indices in network_two.nodes()"""
        
        self.time_till_first_safe_bet = None
        """That's the time that passed from the start of the search till the first safe bet was found"""
        

    def backtrack(self, tracking_state):
        """The recursive backtrack function. Overrides KrissinelMaximumCommonSubgraphFinder.backtrack().
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance 
            The tracking state keeps track of the subgraph that is currently analyzed. And the
            possible matchings that remain to be analyzed.
        
        This function will analyze whether the current search branch can be expanded,
        and if not then check whether the branch is the largest currently found common subgraph.
        If the branch is expandable it will pick a new vertex in the first network, find all possible
        matches for this vertex, and set up a new branch for each of these matches - on which backtrack 
        will be called again. Since the picked vertex does not have to be part of the maximal common subgraph
        it will then set up a branch where this vertex is not in present in the first network, and attempt to call 
        backtrack on this branch.
        
        See also
        --------
        
        KrissinelMaximumCommonSubgraphFinder.backtrack()
        """
        if self.is_extendable(tracking_state):
            current_vertex = self.pick_next_vertex(tracking_state)
            possible_images = self.get_mappable_vertices(tracking_state, current_vertex)
        
            for image in possible_images:
                # only consider this image if it is not a safe bet
                if self.network_two.nodes()[image] not in self.safe_bets.values():
                    # this method will also update the vertex_matching_matrix in the tracking state
                    new_tracking_state = self.create_extended_tracking_state(tracking_state, current_vertex, image)
                    self.backtrack(new_tracking_state)
                # or if it is, than the current vertex better be as well
                elif self.network_one.nodes()[current_vertex] in self.safe_bets.keys():
                    new_tracking_state = self.create_extended_tracking_state(tracking_state, current_vertex, image)
                    self.backtrack(new_tracking_state)

                if self.check_is_safe_bet( current_vertex ):                    
                    break
        
            if self.network_one.nodes()[current_vertex] not in self.safe_bets:
                # take the current vertex out of the possible extensions for the current tracking state
                tracking_state.vertex_matching_matrix[current_vertex,:] = 0
                tracking_state.counts_of_mappable_vertices[current_vertex] = 0
                
                # and run backtrack on this
                self.backtrack(tracking_state)
            
        else:
            self.branches_counter += 1
            if len( tracking_state.id_map.keys() ) > self.max_found_subgraph_size:
                self.max_found_subgraph_size = len( tracking_state.id_map.keys() )
                self.largest_mappings = [tracking_state.id_map]
            elif len( tracking_state.id_map.keys() ) == self.max_found_subgraph_size:
                self.largest_mappings.append( tracking_state.id_map )
                
    def check_is_safe_bet(self, node_index):
        """This function checks whether the current vertex is a safe bet, and updates its status if necessary.
        
           If the node is in the safe_bets dictionary, then it just returns true. If it is in the
           preliminary safe bets dictionary, then it checks whether any of the neighbours are safe bets.
           In that case, it will get the correct image from the current maximum common subgraph and make this mapping
           into a safe bet. In any other case, returns false.
           
        Parameters
        ----------
        
            node_index : int
                index of to the self.network_one.nodes() vector.
                
        Returns
        -------
        
            is_safe_bet : bool
                Whether this node is a safe bet.
        """
        this_frame_id = self.network_one.nodes()[node_index]
        if this_frame_id in self.safe_bets:
            return True
        else:
            if self.largest_mappings != []:
                if this_frame_id in self.largest_mappings[0] and this_frame_id in self.preliminary_safe_bets:
                    adjacent_vertices = self.network_one.neighbors(self.network_one.nodes()[node_index])
                    two_adjacent_vertices_are_safe_bet = False
                    no_of_adjacent_safe_bets = 0
                    for vertex in adjacent_vertices:
                        if ( vertex in self.safe_bets or vertex in self.preliminary_safe_bets ):
                            no_of_adjacent_safe_bets += 1
                            if no_of_adjacent_safe_bets > 1:
                                self.safe_bets[this_frame_id] = self.largest_mappings[0][this_frame_id]
                                two_adjacent_vertices_are_safe_bet = True
                                break
                    return two_adjacent_vertices_are_safe_bet 
            # This line is only reached if we haven't returned True yet
            return False
                
    def get_mappable_vertices(self, tracking_state, vertex_index):
        """Get all vertices in network_two that could be images to the node in network one that is indexed by vertex_index.
        Overrides KrissinelMaximumCommonSubgraphFinder.get_mappable_vertices().
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            see TrackingState or backtrack()

        vertex_index : int
            index for an entry in network_one.nodes
        
        Returns
        -------
        
        mappable_vertices : np.array, dtype = int
            a set of vertices in network_two that have the same degree as vertex,
            are within half the mesh width distance, and don't add an edge to the existing nodes
            that isn't there in network_one, while also not deleting or adding an edge that is or is not there
            in network_one
            
        See also
        --------
        
        KrissinelMaximumCommonSubgraphFinder.get_mappable_vertices()
        """
        
        if self.network_one.nodes()[vertex_index] in self.safe_bets:
            mappable_vertices = [ self.network_two_index_lookup[ self.safe_bets[vertex_index] ] ]
        else:
            mappable_vertices = KrissinelMaximumCommonSubgraphFinder.get_mappable_vertices(self, tracking_state, vertex_index)

        return mappable_vertices
    
    def pick_next_vertex(self, tracking_state):
        """Pick the next vertex in network_one. Picks one of the vertices with the minimal number of possible matches
        Overrides ConnectedMaximuCommonSubgraphFinder.pick_next_vertex().
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            see backtrack() function, or TrackingState
            
        Returns
        -------
            
        vertex_index : int
            an index for an entry network_one.nodes
            
        Pick a node from network_one that hasn't been mapped yet. Pick one that has the minimum number of possible matches.
        If not all safe bets are mapped yet, pick the safe bets first.
        """
        
        unmapped_safe_bets = self.find_unmapped_safe_bets( tracking_state )

        if unmapped_safe_bets == set():
            next_vertex_index = ConnectedMaximumCommonSubgraphFinder.pick_next_vertex(self, tracking_state )
        else:
            next_vertex_index = unmapped_safe_bets.pop()

        return next_vertex_index
    
    def find_unmapped_safe_bets(self, tracking_state):
        """Identify all safe bets that haven't been mapped yet
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            see backtrack() function, or TrackingState
            
        Returns
        -------
        
            unmapped safe bets : set of ints
                all safe bets that have not yet been mapped
        """
        
        mapped_indices = set( [ self.network_one_index_lookup[key] for key in tracking_state.id_map.keys() ] )
        unmapped_safe_bets = set(self.safe_bets).difference( mapped_indices )
        
        return unmapped_safe_bets
    
    def create_extended_tracking_state(self, tracking_state, vertex_one, vertex_two):
        """Create an new tracking state, where tracking state is extended by the mapping from current_vertex 
        to image_vertex, and the possible matches are refined
        
        Overrides and uses ConnectedMaximumCommonSubgraphFinder.create_extended_tracking_state()
        
        This method corresponds to the Refine() method in Krissinel et al.
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            the tracking state that we are planning on copying and extending
            
        vertex_one : int
            index of a node in the network_one.nodes vector
            
        vertex_two : int
            index of a node in the network_two.nodes vector
            
        The method extends tracking_state by the mapping between network_one.nodes[current_vertex]->network_two.nodes[vertex_two],
        and refines the vertex_matching_matrix, such that only mappings are considered that have no edge conflicts with the current nodes.
        
        See also
        --------
        
        ConnectedMaximumCommonSubgraphFinder.create_extended_tracking_state()
        KrissinelMaximumCommonSubgraphFinder.create_extended_tracking_state()
        """
        
        new_tracking_state = ConnectedMaximumCommonSubgraphFinder.create_extended_tracking_state(self, tracking_state, vertex_one, vertex_two)

        # The mapping of vertex one might have just made one of it's neighbours a safe bet
        neighbour_ids_of_vertex_one = self.network_one.neighbors(self.network_one.nodes()[vertex_one])
        
        for neighbour_id in neighbour_ids_of_vertex_one:
            # Proceed if this index is not yet a safe bet but has already been mapped.
            if ( neighbour_id not in self.safe_bets and neighbour_id in new_tracking_state.id_map ):
                neighbours_neighbour_ids = self.network_one.neighbors(neighbour_id)
                mapped_neighbour_counter = 0
                for mapped_neighbour_counter, neighbours_neighbour_id in enumerate( neighbours_neighbour_ids ):
                    if neighbours_neighbour_id not in new_tracking_state.id_map:
                        break
                    else:
                        mapped_neighbour_counter += 1
                if ( mapped_neighbour_counter == len( neighbours_neighbour_ids ) and 
                     self.network_one.degree(neighbour_id) == self.network_two.degree(new_tracking_state.id_map[neighbour_id]) ):
                    if ( self.check_neighbours_are_mapped_in_order(new_tracking_state, neighbour_id)
                         and self.network_one.degree( neighbour_id ) > 3 ):
                        if ( self.check_not_too_many_hexagons(new_tracking_state, neighbour_id) ):
                            if self.timing_is_on:
                                if self.safe_bets == {}:
                                    first_safe_bet = time.clock()
                                    self.time_till_first_safe_bet = first_safe_bet - self.start_time
                            self.safe_bets[ neighbour_id ] = tracking_state.id_map[ neighbour_id ]
                        else:
                            self.preliminary_safe_bets[ neighbour_id ] = tracking_state.id_map[ neighbour_id ]
                    
        return new_tracking_state

    def check_not_too_many_hexagons(self, tracking_state, frame_id):
        """True if less than half of the neighbours (including the centre cell)
           are hexagons.
        
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            the tracking state that we are planning on copying and extending
            
        frame_id : int
            frame_id in network one for which we would like to check whether the neighbours got mapped in the right 
            order.
        """
        
        this_element = self.mesh_one.get_element_with_frame_id(frame_id)
        
        total_no_of_neighbours_including_centre = 1.0
        no_of_hexagons = 0.0
        if this_element.get_num_nodes() == 6:
            no_of_hexagons += 1            
        else:
            return True

        neighbour_element_ids = []

        for local_index, node in enumerate(this_element.nodes):
            next_node = this_element.nodes[ (local_index+1)%this_element.get_num_nodes() ]
            shared_element_ids =  ( set.intersection( set( node.get_adjacent_element_ids() ),
                                                          set( next_node.get_adjacent_element_ids() ) ) )
            shared_element_ids.remove(frame_id)
            if len(shared_element_ids) > 0:
                neighbour_element_ids.append( shared_element_ids.pop() )
                total_no_of_neighbours_including_centre += 1

        for this_id in neighbour_element_ids:
            this_polygon_no = self.mesh_one.get_element_with_frame_id(this_id).get_num_nodes()
            if this_polygon_no == 6:
                no_of_hexagons += 1
                
        not_too_many_hexagons = ( no_of_hexagons/total_no_of_neighbours_including_centre < 0.42 )

        return not_too_many_hexagons
   
    def check_neighbours_are_mapped_in_order(self, tracking_state, frame_id):
        """True if all neighbours of frame_id are mapped in the correct order
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            the tracking state that we are planning on copying and extending
            
        frame_id : int
            frame_id in network one for which we would like to check whether the neighbours got mapped in the right 
            order.
            
        Returns
        -------
        
        neighbours_are_mapped_in_order: bool
            true if all neighbours are mapped in the correct circular order.
        """
        
        this_element = self.mesh_one.get_element_with_frame_id(frame_id)
        
        neighbour_element_ids = []

        for local_index, node in enumerate(this_element.nodes):
            next_node = this_element.nodes[ (local_index+1)%this_element.get_num_nodes() ]
            shared_element_ids =  ( set.intersection( set( node.get_adjacent_element_ids() ),
                                                          set( next_node.get_adjacent_element_ids() ) ) )
            shared_element_ids.remove(frame_id)
            if len(shared_element_ids) > 0:
                neighbour_element_ids.append( shared_element_ids.pop() )

        mapped_neighbour_element_ids = []
        for this_id in neighbour_element_ids:
            mapped_neighbour_element_ids.append( tracking_state.id_map[this_id] )

        second_element = self.mesh_two.get_element_with_frame_id( tracking_state.id_map[frame_id] )
        
        second_neighbour_element_ids = []

        for local_index, node in enumerate(second_element.nodes):
            next_node = second_element.nodes[ (local_index+1)%second_element.get_num_nodes() ]
            shared_element_ids =  ( set.intersection( set( node.get_adjacent_element_ids() ),
                                                          set( next_node.get_adjacent_element_ids() ) ) )
            shared_element_ids.remove( second_element.id_in_frame )
            if len(shared_element_ids) > 0:
                second_neighbour_element_ids.append( shared_element_ids.pop() )

        # The next bit is from stackoverflow: a quick comparision whether two integer lists are in same circular order 
        # http://stackoverflow.com/questions/26924836/how-to-check-whether-two-lists-are-circularly-identical-in-python

        str1 = ' '.join(map(str, second_neighbour_element_ids))
        str2 = ' '.join(map(str, mapped_neighbour_element_ids))

        return str1 in str2 + ' ' + str2
    
class LocalisedSubgraphFinder(ConnectedMaximumCommonSubgraphFinder):
    def __init__(self, mesh_one, mesh_two):
        """A helper class to find maximum common subgraphs. Inherits from
        ConnectedMaximumCommonSubgraphFinder and adapts it to deal with connected subgraphs
        of planar topologies.
        
        Parameters for __init__ are the same as in the parent class.
        
        See also
        --------
        
        ConnectedMaximumCommonSubgraphFinder.__init__()
        KrissinelMaximumCommonSubgraphFinder.__init__()
        """
        ConnectedMaximumCommonSubgraphFinder.__init__( self, mesh_one, mesh_two, 0.0 )


    def find_maximum_common_subgraph(self):
        """Runs the backtracking algorithm to find the maximum
        common subgraph
        
        The maximum common subgraph is found by recursive application
        of the backtrack() function to instances_of_tracking_state. 
        a tracking state keeps track of the current subgraph, and the remaining
        possible matchings between vertices.
        
        See also
        --------
        
        backtrack
        """
        if self.timing_is_on:
            self.start_time = time.clock()
        
        self.linear_backtrack()
        
        if self.timing_is_on:
            end_time = time.clock()
            self.total_execution_time = end_time - self.start_time
            
    def create_starting_tracking_state(self):
        """Create the starting tracking state for the linear subgraph search
        
        Finds an initial mapping that we can be certain about.
        
        Returns:
            first_tracking_state
        """
        
#         print 'start searching for initialisation'

        first_vertices = self.get_vertices_ordered_by_number_of_matches()
        
        first_match_not_yet_found = True
        entry_counter = 0

        while first_match_not_yet_found:
            vertex = first_vertices[ entry_counter ]
#             print 'consider vertex with frame id, position, and degree'
#             vertex_id = self.network_one.nodes()[vertex]
#             print vertex_id
#             print self.network_one.node[vertex_id]['position']
#             print self.network_one.node[vertex_id]['num_neighbours']
            possible_images = self.get_mappable_vertices(self.initial_tracking_state, vertex)
            self.largest_mappings = []
            self.max_found_subgraph_size = 0
            enlarge_local_region = False
            max_reduced_size = 0
            for image in possible_images:
                new_tracking_state, reduced_size = self.create_localised_tracking_state(self.initial_tracking_state, vertex, image, 1)
                self.backtrack( new_tracking_state )
                if self.max_found_subgraph_size >= reduced_size:
                    enlarge_local_region = True
#                     first_mapping = [vertex, image]
                    break
                
#             vertex_polygon_number = self.network_one.node[self.network_one.nodes()[vertex]]['num_neighbours']
#             if enlarge_local_region and vertex_polygon_number > 4:
            if enlarge_local_region:
                for image in possible_images:
                    new_tracking_state, reduced_size = self.create_localised_tracking_state(self.initial_tracking_state, vertex, 
                                                                                            image)
                    self.backtrack( new_tracking_state )
                    if reduced_size > max_reduced_size:
                        max_reduced_size = reduced_size

                if self.max_found_subgraph_size >= max_reduced_size:
                    vertex_is_mapped_equally_in_all_largest_mappings, image_id = self.vertex_is_mapped_equally_in_all_largest_mappings(vertex)
                    if vertex_is_mapped_equally_in_all_largest_mappings:
                        first_match_not_yet_found = False
                        first_mapping = [vertex, self.network_two_index_lookup[image_id]]
                        break
            entry_counter += 1
                
        first_tracking_state = self.create_extended_tracking_state(self.initial_tracking_state, 
                                                                   first_mapping[0], first_mapping[1])
        
#         print 'needed so many tries for first mapping'
#         print entry_counter
#         print 'initial mapping is'
#         print self.network_one.node[self.network_one.nodes()[first_mapping[0]]]['position']
#         print self.network_two.node[self.network_two.nodes()[first_mapping[1]]]['position']

#         print 'finished searching for initialisation'
        return first_tracking_state
    
    def get_vertices_ordered_by_number_of_matches(self):
        """Get the order in which we consider vertices for the search of the first match.
        
        Returns
        -------
        
        ordered_vertices : list of ints
            entries are indices of network_one.nodes(), The vertices in network one are ordered by
            their number of possible matches.
        """
        sorted_indices_of_vertices_count = np.argsort( self.initial_tracking_state.counts_of_mappable_vertices )

        ordered_vertices = []
        for index in sorted_indices_of_vertices_count:
            if ( ( self.initial_tracking_state.counts_of_mappable_vertices[index] > 0 )
                 and ( self.network_one.degree( self.network_one.nodes()[index] ) > 3 ) ):
                ordered_vertices.append( index )

        return ordered_vertices
    
    def create_localised_tracking_state(self, tracking_state, vertex_one, vertex_two = None, cutoff_distance = 2 ):
        """Create a tracking state in which all vertices that are too far away from the current vertex are not
           mappable in network_one.
           
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            the tracking state that we are planning to edit
            
        vertex_one : int
            index in network_one.nodes(), the vertex around which we wish to localise.

        vertex_two : int
            index in network_two.nodes(), the image we wish to consider.
            
        cutoff_distance : int
            size of the localised subgraph in network step length.
            
        Returns
        -------
        
        localised_tracking_state : TrackingState instance
            the edited tracking state. Entries in counts_of_mappable_vertices and vertex_matching_matrix of all
            vertices that are beyond a threshold distance to the `vertex' argument 
            
        reduced_size : int
            The number of vertices in network_one that are not removed in localised_tracking_state (including
            vertices that are not removed and are not mappable).      
        """
        vertex_id = self.network_one.nodes()[vertex_one]
        second_class_neighbour_ids = nx.single_source_shortest_path_length(self.network_one, vertex_id, 
                                                                           cutoff=cutoff_distance ).keys()
        second_class_neighbour_indices = [ self.network_one_index_lookup[ frame_id ] for
                                           frame_id in second_class_neighbour_ids]

        vertices_that_are_adjacent_to_vertex_one = self.network_one.neighbors(self.network_one.nodes()[vertex_one])
        indices_of_vertices_that_are_adjacent_to_vertex_one = np.array( [ self.network_one_index_lookup[vertex] 
                                                                          for vertex in vertices_that_are_adjacent_to_vertex_one])
       
        if vertex_two != None:
            # here is the copy-pasted code for making extended trackingstates
            # here is the copy-pasted code for making extended Krissinel tracking states

            # np.nonzero will return a tuple with one entry, we are interested in the entry, which is an array with all
            # the indices where the mappable vertices count is not zero. Hence, we take the zero'th element of the nonzero return value
            list_of_still_mappable_vertices_in_old_state = np.nonzero( tracking_state.counts_of_mappable_vertices )[0]
            localised_tracking_state = TrackingState()
            localised_tracking_state.id_map = tracking_state.id_map.copy()
            localised_tracking_state.id_map[self.network_one.nodes()[vertex_one]] = self.network_two.nodes()[vertex_two]
            
            # Ok, here we start introducing a reduced structure
            localised_tracking_state.vertex_matching_matrix = np.zeros( (len(second_class_neighbour_indices), 
                                                                        tracking_state.vertex_matching_matrix.shape[1]), dtype = 'int')
            localised_tracking_state.counts_of_mappable_vertices = np.zeros( len(second_class_neighbour_indices), dtype = 'int')
            localised_tracking_state.connections_to_current_subgraph = np.zeros( len(second_class_neighbour_indices), dtype = 'int')
            
            localised_tracking_state.sparse_matrix_lookup = np.array(second_class_neighbour_indices)
            for reduced_index, vertex_index in enumerate(second_class_neighbour_indices):
                if vertex_index != vertex_one:
                    count_of_matching_vertices = 0
                    if vertex_index in indices_of_vertices_that_are_adjacent_to_vertex_one:
                        localised_tracking_state.connections_to_current_subgraph[reduced_index] = tracking_state.connections_to_current_subgraph[vertex_index] + 1
                    else:
                        localised_tracking_state.connections_to_current_subgraph[reduced_index] = tracking_state.connections_to_current_subgraph[vertex_index]
                    for column_index_in_vertex_matching_matrix in range(tracking_state.counts_of_mappable_vertices[vertex_index]):
                        image_index = tracking_state.vertex_matching_matrix[vertex_index, column_index_in_vertex_matching_matrix]
                        if image_index != vertex_two:
                            if self.pairing_does_not_introduce_edge_conflicts((vertex_one, vertex_index), (vertex_two, image_index)):
                                localised_tracking_state.vertex_matching_matrix[reduced_index][count_of_matching_vertices] = image_index
                                count_of_matching_vertices += 1
                    localised_tracking_state.counts_of_mappable_vertices[reduced_index] = count_of_matching_vertices

            localised_tracking_state.not_blacklisted_vector = tracking_state.not_blacklisted_vector

            if (tracking_state.adjacent_vertices != None):
                localised_tracking_state.adjacent_vertices = np.unique( np.hstack( ( tracking_state.adjacent_vertices,
                                                                               indices_of_vertices_that_are_adjacent_to_vertex_one ) ) )
            else:
                localised_tracking_state.adjacent_vertices = indices_of_vertices_that_are_adjacent_to_vertex_one 
        else:
            localised_tracking_state = TrackingState()
            localised_tracking_state.id_map = tracking_state.id_map.copy()

            localised_tracking_state.vertex_matching_matrix = np.zeros( (len(second_class_neighbour_indices), 
                                                                        tracking_state.vertex_matching_matrix.shape[1]), dtype = 'int')
            localised_tracking_state.counts_of_mappable_vertices = np.zeros( len(second_class_neighbour_indices), dtype = 'int')
            localised_tracking_state.connections_to_current_subgraph = np.zeros( len(second_class_neighbour_indices), dtype = 'int')
            
            localised_tracking_state.sparse_matrix_lookup = np.array( second_class_neighbour_indices )
            for reduced_index, vertex_index in enumerate(second_class_neighbour_indices):
                if vertex_index != vertex_one:
                    localised_tracking_state.connections_to_current_subgraph[reduced_index] = tracking_state.connections_to_current_subgraph[vertex_index]
                    count_of_matching_vertices = 0
                    localised_tracking_state.vertex_matching_matrix[reduced_index] = tracking_state.vertex_matching_matrix[vertex_index]
                    localised_tracking_state.counts_of_mappable_vertices[reduced_index] = tracking_state.counts_of_mappable_vertices[vertex_index]
            localised_tracking_state.adjacent_vertices = np.copy( tracking_state.adjacent_vertices )
            localised_tracking_state.not_blacklisted_vector = tracking_state.not_blacklisted_vector
                
        localised_tracking_state.inverse_sparse_matrix_lookup = { vertex_index : counter for counter, vertex_index 
                                                           in enumerate(localised_tracking_state.sparse_matrix_lookup)}
        reduced_size = len(second_class_neighbour_indices)

        return localised_tracking_state, reduced_size

    def linear_backtrack(self):
        """The recursive backtrack function
        
        This function replaces the recursive backtrack function with a linear search through the mesh
        and iteratively extends the current maximum common subgraph optimally.
        """
        
        tracking_state = self.create_starting_tracking_state()

        subgraph_is_still_changing = True
        subgraph_is_extendible = True
        while subgraph_is_still_changing:
#             print 'blacklist loop'
            old_id_map = tracking_state.id_map.copy()
            while subgraph_is_extendible:
                vertex = self.pick_next_vertex(tracking_state, is_global = True )
#                 if self.network_one.nodes()[vertex] == 0:
#                     print 'trying to map frame id 0'
                    
    #             print 'mapping vertex at position'
    #             print self.network_one.node[self.network_one.nodes()[vertex]]['position']
                possible_images = self.get_mappable_vertices(tracking_state, vertex)
#                 if self.network_one.nodes()[vertex] == 0:
#                     print 'possible images are'
#                     print possible_images

                self.max_found_subgraph_size = 0
                self.largest_mappings = []
                
                for image in possible_images:
                    new_tracking_state,_ = self.create_localised_tracking_state( tracking_state, vertex, image )
                    self.backtrack( new_tracking_state )
                    
                reduced_localised_tracking_state,_ = self.create_localised_tracking_state(tracking_state, vertex)
                self.backtrack( reduced_localised_tracking_state )
                
                tracking_state = self.extend_tracking_state_and_largest_mapping_if_feasible( tracking_state, vertex )
                
                subgraph_is_extendible = self.is_extendable(tracking_state, is_global = True)
                
            subgraph_is_still_changing = not old_id_map == tracking_state.id_map
            tracking_state.not_blacklisted_vector[:] = True

            subgraph_is_extendible = self.is_extendable(tracking_state, is_global = True)

#         while subgraph_is_extendible:
#             vertex = self.pick_next_vertex(tracking_state, is_global = True )
# #             print 'mapping vertex at position'
# #             print self.network_one.node[self.network_one.nodes()[vertex]]['position']
#             possible_images = self.get_mappable_vertices(tracking_state, vertex)
# 
#             self.max_found_subgraph_size = 0
#             self.largest_mappings = []
#             
#             for image in possible_images:
#                 new_tracking_state,_ = self.create_localised_tracking_state( tracking_state, vertex, image )
#                 self.backtrack( new_tracking_state )
#                 
#             reduced_localised_tracking_state,_ = self.create_localised_tracking_state(tracking_state, vertex)
#             self.backtrack( reduced_localised_tracking_state )
#             
#             tracking_state = self.extend_tracking_state_and_largest_mapping_if_feasible( tracking_state, vertex )
#             
#             subgraph_is_extendible = self.is_extendable(tracking_state, is_global = True)
#        
#        for vertex in blacklisted_vertices:
#
#            possible_images = self.get_mappable_vertices(tracking_state, vertex)
#
#            self.max_found_subgraph_size = 0
#            self.largest_mappings = []
#
#            for image in possible_images:
#                new_tracking_state,_ = self.create_localised_tracking_state( tracking_state, vertex, image )
#                self.backtrack( new_tracking_state )
#                
#            reduced_localised_tracking_state,_ = self.create_localised_tracking_state(tracking_state, vertex)
#            self.backtrack( reduced_localised_tracking_state )
#            
#            tracking_state = self.extend_tracking_state_and_largest_mapping_if_feasible( tracking_state, vertex )

        self.largest_mappings = [tracking_state.id_map]
            
    def extend_tracking_state_and_largest_mapping_if_feasible(self, tracking_state, vertex):
        """Extend the current tracking state with the current vertex if it is part of all largest mappings
        and had the same image in all largest mappings.
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            The tracking state representing the current search for a maximum common_subgraph
        
        vertex : int
            index in self.network_one.nodes(), the vertex by which the current tracking state is
            considered to be extended.
            
        Returns
        -------
        
        extended_tracking_state : TrackingState instance
            the tracking state is either extended or reduced by the current vertex. In the latter case
            the tracking state is not copied, but edited and returned.
        """
        vertex_id = self.network_one_node_iterator[vertex]

        vertex_is_mapped_equally_in_all_largest_mappings, image = self.vertex_is_mapped_equally_in_all_largest_mappings(vertex)
                
#         if self.network_one.nodes()[vertex] == 0:
#             print 'vertex 0 has position'
#             print self.network_one.node[vertex_id]['position']
#             print 'the current id map has this size'
#             print len(tracking_state.id_map)
#             print 'vertices that in the largest mappings'
#             for largest_mapping in self.largest_mappings:
#                 print 'inspecting a largest mapping'
#                 for frame_id in largest_mapping:
#                     if frame_id not in tracking_state.id_map:
#                         print 'vertex at this position got mapped additionally'
#                         print self.network_one.node[frame_id]['position']
#                         print 'onto:'
#                         print self.network_two.node[largest_mapping[frame_id]]['position']
                        
        if vertex_is_mapped_equally_in_all_largest_mappings:
#             if self.network_one.nodes()[vertex] == 0:
#                 print 'vertex 0 is mapped equally in all largest mappings'
            new_tracking_state = self.extend_tracking_state(tracking_state, vertex, self.network_two_index_lookup[image] )
            new_tracking_state.not_blacklisted_vector[vertex] = True
#             print 'extend tracking state with vertices'
#             print 'origin:'
#             print self.network_one.node[vertex_id]['position']
#             print 'image'
#             print self.network_two.node[image]['position']
        else:
#             if self.network_one.nodes()[vertex] == 0:
#                 print 'vertex 0 is not mapped equally in all largest mappings'
            new_tracking_state = tracking_state
            tracking_state.not_blacklisted_vector[vertex] = False
        
        return new_tracking_state

    def extend_tracking_state(self, tracking_state, vertex_one, vertex_two):
        """Create an new tracking state, where tracking state is extended by the mapping from current_vertex 
        to image_vertex, and the possible matches are refined
        
        Overrides and uses KrissinelMaximumCommonSubgraphFinder.create_extended_tracking_state()
        
        This method corresponds to the Refine() method in Krissinel et al.
        
        Parameters
        ----------
        
        tracking_state : TrackingState instance
            the tracking state that we are planning on copying and extending
            
        vertex_one : int
            index of a node in the network_one.nodes vector
            
        vertex_two : int
            index of a node in the network_two.nodes vector
            
        The method extends tracking_state by the mapping between network_one.nodes[current_vertex]->network_two.nodes[vertex_two],
        and refines the vertex_matching_matrix, such that only mappings are considered that have no edge conflicts with the current nodes.
        """
        
        # np.nonzero will return a tuple with one entry, we are interested in the entry, which is an array with all
        # the indices where the mappable vertices count is not zero. Hence, we take the zero'th element of the nonzero return value
        # list_of_still_mappable_vertices_in_old_state = np.nonzero( tracking_state.counts_of_mappable_vertices )[0]

        vertices_that_are_adjacent_to_vertex_one = self.network_one.neighbors(self.network_one_node_iterator[vertex_one])
        indices_of_vertices_that_are_adjacent_to_vertex_one = np.array( [ self.network_one_index_lookup[vertex] 
                                                                          for vertex in vertices_that_are_adjacent_to_vertex_one])
 
        vertices_that_are_adjacent_to_vertex_two = self.network_two.neighbors(self.network_two_node_iterator[vertex_two])
        indices_of_vertices_that_are_adjacent_to_vertex_two = np.array( [ self.network_two_index_lookup[vertex] 
                                                                          for vertex in vertices_that_are_adjacent_to_vertex_two])


        reverse_indices = []
        for index_adjacent_to_vertex_two in indices_of_vertices_that_are_adjacent_to_vertex_two:
            reverse_indices += np.where(tracking_state.vertex_matching_matrix == index_adjacent_to_vertex_two)[0].tolist()
 
        indices_to_be_inspected = set(reverse_indices + indices_of_vertices_that_are_adjacent_to_vertex_one.tolist())

        new_tracking_state = tracking_state
        
        new_tracking_state.id_map[self.network_one_node_iterator[vertex_one]] = self.network_two_node_iterator[vertex_two]
        
        for vertex_index in indices_to_be_inspected:
            if vertex_index != vertex_one:
                counts_of_mappable_vertices_this_vertex = tracking_state.counts_of_mappable_vertices[vertex_index]
                mappable_vertices_this_vertex = np.array(tracking_state.vertex_matching_matrix[vertex_index])
                tracking_state.vertex_matching_matrix[vertex_index,:] = 0
                count_of_matching_vertices = 0
                for column_index_in_vertex_matching_matrix in range(counts_of_mappable_vertices_this_vertex):
                    image_index = mappable_vertices_this_vertex[column_index_in_vertex_matching_matrix]
                    if image_index != vertex_two:
                        # check whether this potential new pairing would have a conflicting edge with the pairing
                        # that we just added
                        if self.pairing_does_not_introduce_edge_conflicts((vertex_one, vertex_index), (vertex_two, image_index)):
                            new_tracking_state.vertex_matching_matrix[vertex_index][count_of_matching_vertices] = image_index
                            count_of_matching_vertices += 1
                new_tracking_state.counts_of_mappable_vertices[vertex_index] = count_of_matching_vertices
                
        new_tracking_state.counts_of_mappable_vertices[vertex_one] = 0
        new_tracking_state.vertex_matching_matrix[vertex_one, :] = 0
        
        new_tracking_state.connections_to_current_subgraph[indices_of_vertices_that_are_adjacent_to_vertex_one] += 1

        if (tracking_state.adjacent_vertices != None):
            new_tracking_state.adjacent_vertices = np.unique( np.hstack( ( tracking_state.adjacent_vertices,
                                                                           indices_of_vertices_that_are_adjacent_to_vertex_one ) ) )
        else:
            new_tracking_state.adjacent_vertices = indices_of_vertices_that_are_adjacent_to_vertex_one 
            
        return new_tracking_state

    def vertex_is_mapped_equally_in_all_largest_mappings(self, vertex):
        """Check whether the element with network index vertex is mapped equally in all current largest mappings.
        
        Parameters
        ----------
        
        vertex : int
            index of self.network_one.nodes()
            
        Returns
        -------
        
        vertex_is_mapped_equally_in_all_largest_mappings : bool
            True if this vertex has a mapping in all largest mappings, and if these mappings are all the same.
            False otherwise
            
        image_id : int
            frame id of the image that is found. None if the vertex is not mapped equally in all largest mappings.
        """ 

        vertex_id = self.network_one.nodes()[vertex]
        vertex_is_mapped_equally_in_all_largest_mappings = True
#         if self.network_one.nodes()[vertex] == 0:
#             print 'vertex 0 is in first largest mapping'
#             print vertex_id in self.largest_mappings[0]
#             print 'the largest mapping has this size'
#             print len(self.largest_mappings[0])
#             print 'and there are so many largest mappings:'
#             print len(self.largest_mappings)
        image = None
        if vertex_id in self.largest_mappings[0]:
            image = self.largest_mappings[0][vertex_id]
            vertex_is_in_first_mapping = True
        else:
            vertex_is_mapped_equally_in_all_largest_mappings = False
            vertex_is_in_first_mapping = False
#             print 'found vertex not in first mapping'
#             print self.network_one.node[vertex_id]['position']
#             print 'vertex has id'
#             print vertex_id
#             print 'number of largest mappings'
#             print len(self.largest_mappings)
#             inverse_dictionary_one = { v: k for k, v in self.largest_mappings[0].items() }
#             inverse_image = inverse_dictionary_one[46]
#             print 'the vertex that got mapped onto this one in this image is'
#             print inverse_image
#             print self.network_two.node[inverse_image]['position']

        if vertex_is_in_first_mapping:
            for large_mapping in self.largest_mappings:
                if vertex_id in large_mapping:
                    if large_mapping[vertex_id] != image:
                        vertex_is_mapped_equally_in_all_largest_mappings = False
#                         print 'found vertex with different images in the mappings'
#                         print self.network_one.node[vertex_id]['position']
                        break
                else:
                    vertex_is_mapped_equally_in_all_largest_mappings = False
#                     print 'found vertex that is not in all mappings'
#                     print self.network_one.node[vertex_id]['position']
                    break
        
        return vertex_is_mapped_equally_in_all_largest_mappings, image
 
class SlowMaximumCommonSubgraphFinder:
    """A helper class for finding maximum common subgraphs
    
    The class implements the algorithm described in
    
    Evgeny B. Krissinel, Kim Henrick. Common subgraph isomorphism detection
    by backtracking search. Softw. Pract. Exper. 34, pp 591-607 (2007)
    DOI: 10.1002/spe.588
    """
    def __init__( self, mesh_one, mesh_two ):
        """Initialise algorithm with two Mesh-type objects
        
        Parameters
        ----------
        
        mesh_one: Mesh type
            One of the two meshes to which we would like to apply the maximum
            common subgraph algorithm
            
        mesh_two: Mesh type
            second of the two meshes
            
        This method will setup the class but not run the algorithm. Run the algorithm
        by calling find_maximum_common_subgraph
        """
        
        self.network_one = mesh_one.generate_network()
        """Network of first mesh"""

        self.network_two = mesh_two.generate_network()
        """Network of first mesh"""

        self.nodes_one = set(self.network_one.nodes())
        """all nodes in first, will get repeatedly reduced and increased again during
        the search"""
        
        self.width_one = mesh_one.calculate_width()
        """Width of the first mesh"""

        self.height_one = mesh_one.calculate_height()
        """Width of the second mesh"""

        self.n_max = 0
        """maximum size of currently found subgraph"""
        
        self.largest_mappings = []
        """list of all graphs that have the current maximal number of edges"""

        self.branches_counter = 0
        """keeps track of the number of branches we have analysed"""

    def find_maximum_common_subgraph(self):
        """Runs the backtracking algorithm to find the maximum
        common subgraph
        
        The maximum common subgraph is found by recursive application
        of the backtrack() function to copies of id_map. An id_map
        is a dictionary whose keys are element ids in mesh_one, whereas the 
        values are element ids in mesh_two.
        """
        initial_id_map = {}
        self.backtrack(initial_id_map)

    def backtrack(self, id_map):
        """The recursive backtrack function
        
        Parameters
        ----------
        
        id_map: dict 
            dictionary with integer keys and values. keys represent nodes
            in network_one, whereas values are the corresponding nodes in
            network_two. id_map represents the node mapping at the current
            location in the search tree.
        
        This function will analyse whether the current search branch can be expanded,
        and if not then check whether the branch is the largest currently found common subgraph.
        If the branch is expandable it will pick a new vertex in the first network, find all possible
        matches for this vertex, and set up a new branch for each of these matches - on which backtrack 
        will be called again. Since the picked vertex does not have to be part of the maximal common subgraph
        it will then set up a branch where this vertex is not in present in the first network, and attempt to call 
        backtrack on this branch.
        """
        if self.extendable(id_map):
            current_vertex = self.pick_next_vertex(id_map)
            possible_images = self.get_mappable_vertices(current_vertex, id_map)
        
            for image in possible_images:
                new_id_map_branch = id_map.copy()
                new_id_map_branch[current_vertex] = image
                self.backtrack(new_id_map_branch)
        
            self.nodes_one.remove(current_vertex)
            self.backtrack(id_map)
            self.nodes_one.add(current_vertex)
        else:
            self.branches_counter += 1
            if len(id_map.keys()) > self.n_max:
                self.n_max = len(id_map.keys())
                self.largest_mappings = [id_map]
            elif len(id_map.keys()) == self.n_max:
                self.largest_mappings.append(id_map)
    
    def extendable(self, id_map):
        """Check whether the current subgraph is extendable
        
        Parameters
        ----------
        
        id_map : dict
            see backtrack() function
            
        Returns
        -------
        
        is_extendable : bool
            currently, this only checks whether there are any unmatched vertices left in network one
        """
        vertices_left_in_network_one = set.difference(self.nodes_one, id_map.keys())
        return len(vertices_left_in_network_one) > 0
    
    def pick_next_vertex(self, id_map):
        """Pick the next vertex in network_one
        
        Parameters
        ----------
        
        id_map : dict
            see backtrack() function
            
        Returns
        -------
            
        vertex : int
            The id of a network node that hasn't been mapped yet
            
        Pick a node from network_one that hasn't been mapped yet
        """
        set_of_already_mapped_vertices = set(id_map.keys())
        set_of_not_yet_mapped_vertices = set.difference(self.nodes_one, set_of_already_mapped_vertices)
        return set_of_not_yet_mapped_vertices.pop()

    def get_mappable_vertices(self, vertex, id_map):
        """Get all vertices in network_two that could be images the argument vertex in network_one
        
        Parameters
        ----------
        
        vertex : int
            id of a node in network_one
        
        id_map : dict
            see backtrack() function    

        Returns
        -------
        
        mappable_vertices : set of ints
            a set of vertices in network_two that have the same degree as vertex,
            are within half the mesh width distance, and don't add an edge to the existing nodes
            that isn't there in network_one, while also not deleting an edge that is there
            in network_one
        """

        image_of_vertex_edges = self.get_image_of_edges_from_vertex_in_network_one(vertex, id_map)
    
        set_of_already_mapped_vertices_in_two = set(id_map.values())
        set_of_not_yet_mapped_vertices_in_two = set.difference(set(self.network_two.nodes()), set_of_already_mapped_vertices_in_two)

        mappable_vertices = []
        for vertex_two in set_of_not_yet_mapped_vertices_in_two:
            vertex_two_is_mappable = \
                (self.network_two.node[vertex_two]['position'][0] - self.network_one.node[vertex]['position'][0]) < 0.25*self.width_one and \
                (self.network_two.node[vertex_two]['position'][1] - self.network_one.node[vertex]['position'][1]) < 0.25*self.height_one and \
                self.network_one.degree(vertex) == self.network_two.degree(vertex_two)

            if vertex_two_is_mappable:
       
                vertex_two_edges = self.get_edges_from_vertex_in_network_two(vertex_two, id_map)
                
                edge_difference  = set.symmetric_difference(image_of_vertex_edges, vertex_two_edges)
                if edge_difference == set():
                    mappable_vertices.append(vertex_two)
                
        return mappable_vertices
    
    def get_image_of_edges_from_vertex_in_network_one(self, vertex, id_map):
        """Find the image of all vertices in the already mapped nodes of network_one to which vertex connects
        
        Parameters
        ----------
        
        vertex : int
            index of vertex for which we would like to know the connectivity
            
        id_map : dict
           see backtrack() function
           
        Returns
        -------
        
        set_of_edge_images : set of ints
            indices of vertices that are images in network_two of all the nodes in network_one that vertex connects to
            
        """

        # find all endpoints of edges from vertex
        vertex_edges = set(np.unique(np.array(self.network_one.edges(vertex))[:]))
        vertex_edges.remove(vertex)

        # all edges connecting to already mapped nodes
        edges_to_remove = []
        for end_point in vertex_edges:
            if end_point not in id_map.keys():
                edges_to_remove.append(end_point) 
        vertex_edges.difference_update(edges_to_remove)
            
        # image of the edge end points in network_two
        image_of_vertex_edges = set()
        for end_point in vertex_edges:
            image_of_vertex_edges.add(id_map[end_point])
            
        return image_of_vertex_edges
 
    def get_edges_from_vertex_in_network_two(self, vertex_two, id_map):
        """For vertex_two in network_two, this will find the already mapped vertices in network_two that it
        connects to.
        
        Parameters
        ----------
        
        vertex_two : int
            index of vertex in network_two for which we find the connected vertices
        
        Returns
        -------
        
        edge_connections : set of ints
            all the vertices in network_two that this node connects to
            
        """
        
        # all endpoints of edges from vertex
        vertex_two_edges = set(np.unique(np.array(self.network_two.edges(vertex_two))[:]))
        vertex_two_edges.remove(vertex_two)

        # all edges connecting to existing nodes
        edges_to_remove = []
        for end_point in vertex_two_edges:
            if end_point not in id_map.values():
                edges_to_remove.append(end_point) 
        vertex_two_edges.difference_update(edges_to_remove)

        return vertex_two_edges
    
   
def generate_subgraph_network(mesh_two, id_map):
    """generate a network in the second mesh with labels from
    the first mesh, as provided by id_map - experimental, not yet tested
    """
    
    network_edges = []

    for node in mesh_two.nodes:
        subgraph_adjacent_element_ids = []
        for id in node.get_adjacent_element_ids():
            subgraph_adjacent_element_ids.append(mesh_two.get_element_with_frame_id(id).global_id)
        new_connections_between_cells = itertools.combinations(subgraph_adjacent_element_ids, 2)
        for tuple in new_connections_between_cells:
            connection = list(tuple)
            connection_reversed = list(connection)
            connection_reversed.reverse()
            if ( connection not in network_edges ) and\
                (connection_reversed not in network_edges):
                network_edges.append(connection)
        
    network = nx.Graph()
    network.add_edges_from(network_edges)
        
    for element in mesh_two.elements:
        network.node[element.global_id]['position'] = element.calculate_centroid()
        
    return network
