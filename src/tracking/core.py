# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""In this the main tracking functions are defined
"""
import sys
import os
from .maximum_common_subgraph_finder import *
import mesh
from mesh.in_out import _natural_keys
import glob
import copy
import warnings
#from networkx.algorithms.components.connected import connected_component_subgraphs
# # # # # # # # # from networkx import connected_component_subgraphs

def track(mesh_one, mesh_two):
    """Find a mapping between the cell ids in both frames and assigns the global ids accordingly.
    
    Parameters
    ----------
    
    mesh_one : Mesh type
        First mesh
    
    mesh_two : Mesh type
        Second mesh
        
    Returns
    -------
    
    mapped_ids : the ids of elements that were identified in both meshes
    """
    
    subgraph_finder = LocalisedSubgraphFinder(mesh_one, mesh_two)
    subgraph_finder.find_maximum_common_subgraph()

    post_processor = PostProcessor(mesh_one, mesh_two, subgraph_finder.largest_mappings)
    post_processor.index_global_ids_from_largest_mappings()
    
    post_processor.tidy_current_mapping()
    mapped_ids = post_processor.post_process_with_data()
    
    return mapped_ids

def track_and_write_sequence(input_path, output_path, start_number = 1, number_meshes = None):
    """Reads a sequence and writes the tracked data into consecutive meshes
    
    Cells that are present in multiple frames will have the same global ids,
    and each other cell will have a distinct non-recurring global id.

    Parameters
    ----------
    
    input_path : string
        filename of seedwater-segmented data frames, without the file-endings
        and numberings
        
    output_path : string
        filename where the output should be saved, without file ending
        this name will be extended with a number and .mesh for each segmented
        frame
        
    start_number : int
        mesh number to be started with (indexing starts at one)
        
    number_meshes : int
        index of the last mesh we want to track (indexing starts at one)
    """

    mesh_sequence = mesh.read_sequence_from_data(input_path, start_number, number_meshes)
    previous_sequence = mesh.read_sequence_from_data(input_path, start_number, number_meshes)   
    next_sequence = mesh.read_sequence_from_data(input_path, start_number, number_meshes)
    
    # track all consecutive time frames individually
    step_sequence = []
    for counter, this_mesh in enumerate(mesh_sequence):
        if counter > 0:
            previous_mesh = previous_sequence[counter -1]
            corresponding_mesh = next_sequence[counter]
            try:
                track(previous_mesh, corresponding_mesh)
            except FirstIndexException:
                print("Could not find first index in tracking step " + str(counter))
            step_sequence.append([previous_mesh, corresponding_mesh])
    
    # give global ids to the first mesh
    global_ids = []
    for counter, element in enumerate(mesh_sequence[0].elements):  
        element.global_id = counter
        global_ids.append(counter)
        element.is_in_reduced_mcs_previous = False
    mesh_sequence[0].index_global_ids()

    # trace global ids through all the meshes, making new ones if necessary
    for counter, this_mesh in enumerate(mesh_sequence):
        if counter == 0:
            corresponding_mesh_next_step = step_sequence[counter][0]
            for element_counter, element in enumerate(this_mesh.elements):
                element.is_in_reduced_mcs_next = corresponding_mesh_next_step.elements[element_counter].is_in_reduced_mcs_next
        if counter > 0:
            previous_mesh = step_sequence[counter - 1][0]
            corresponding_mesh = step_sequence[counter - 1][1]

            if counter < len(step_sequence):
                corresponding_mesh_next_step = step_sequence[counter][0]

            for element_counter, element in enumerate(this_mesh.elements):
                corresponding_element = corresponding_mesh.get_element_with_frame_id(element.id_in_frame)
                this_global_id = corresponding_element.global_id
                if this_global_id is None:
                    new_global_id = max(global_ids) + 1
                    global_ids.append( max(global_ids) + 1 )
                    element.global_id = new_global_id
                    element.is_new = True
                else:
                    previous_frame_id = previous_mesh.get_element_with_global_id(this_global_id).id_in_frame
                    previous_global_id = mesh_sequence[counter - 1].get_element_with_frame_id(previous_frame_id).global_id
                    element.global_id = previous_global_id
                
                try:
                    element.is_in_reduced_mcs_previous = corresponding_element.is_in_reduced_mcs_previous
                except:
                    element.is_in_reduced_mcs_previous = False

                if counter < len(step_sequence):
                    try:
                        element.is_in_reduced_mcs_next = corresponding_mesh_next_step.elements[element_counter].is_in_reduced_mcs_next
                    except(AttributeError):
                        element.is_in_reduced_mcs_next = False
                else:
                    element.is_in_reduced_mcs_next = False
            this_mesh.index_global_ids()
        
    #now, save the mesh sequence
    for counter, this_mesh in enumerate(mesh_sequence):
        this_file_name = output_path + str(start_number + counter - 1) + '.mesh'
        this_mesh.save(this_file_name)

def analyse_tracked_sequence(input_path):
    """collect summary statistics on tracked data
    
    Parameters
    ----------
    
    input_path : string
        Path to the sequence that should be analysed. 
        Sequences are numbered already tracked meshes.
        
    Returns
    -------
    
    data_collector : DataCollector instance
        This object has member variables for various summary statistics
    """
    mesh_sequence = mesh.load_sequence(input_path)
    
    return DataCollector(mesh_sequence)

def plot_tracked_sequence( sequence_path, image_path, segmented_path, out_path ):
    """Plot a tracked sequence of meshes.
    
    This creates three types of plots for the entire sequence.  
    The first type of plot overlays the experimental data, the segmentation, and the 
    tracking outcome. Each tracked cell is given an individual colour and an id that is 
    included in the overlay.
    
    The second type of plots illustrates the maximum common subgraphs.
    
    The third type of plots shows the tracked tesselation of polygons
    
    Parameters
    ----------
    
    sequence_path : string
        path to the tracked mesh sequence (contains a series of .mesh files) 
        
    image_path : string
        path to the sequence of original images
        
    segmented_path : string
        path where the overlay should be saved. Will be created if required.
    """

    mesh_sequence = mesh.load_sequence( sequence_path )
    
    list_of_image_files = glob.glob( os.path.join( image_path , '*.tif') )
    list_of_image_files.sort(key=_natural_keys)

    list_of_segmented_files = glob.glob( os.path.join( segmented_path , '*.tif') )
    list_of_segmented_files.sort(key=_natural_keys)

    # get maximal global id
    max_global_id = 0
    for mesh_instance in mesh_sequence:
        this_max_global_id = mesh_instance.get_max_global_id()
        if this_max_global_id > max_global_id:
            max_global_id = this_max_global_id
    
    if not os.path.isdir(out_path):
        os.mkdir( out_path )
        
    overlay_path = os.path.join(out_path, 'overlay')
    if not os.path.isdir(overlay_path):
        os.mkdir( overlay_path )

    polygon_path = os.path.join(out_path, 'polygons')
    if not os.path.isdir(polygon_path):
        os.mkdir( polygon_path )
        
    mcs_path = os.path.join(out_path, 'mcs')
    if not os.path.isdir(mcs_path):
        os.mkdir( mcs_path )

    for mesh_counter, mesh_instance in enumerate( mesh_sequence ):
        this_image_path = list_of_image_files[mesh_counter]
        print(list_of_segmented_files)
        this_segmented_path = list_of_segmented_files[mesh_counter]
        out_file_name = os.path.split( this_image_path.replace('.tif', '_overlay.png') )[1]

        overlay_file_path = os.path.join(overlay_path, out_file_name)
        mesh_instance.plot_tracked_data(overlay_file_path, this_image_path, this_segmented_path, max_global_id)

        polygon_file_name = os.path.join( polygon_path, out_file_name )
        mesh_instance.plot( polygon_file_name, color_by_global_id = True,
                            total_number_of_global_ids = max_global_id)
        
        mcs_file_name = os.path.join( mcs_path, out_file_name )
        mesh_instance.plot( mcs_file_name, color_by_global_id = True,
                            total_number_of_global_ids = max_global_id, reduced_mcs_only = True )

class DataCollector():
    """A class for analysing tracked sequences."""
    
    def __init__(self, mesh_sequence):
        """The constructor of the DataCollector
        
        Parameters
        ----------
        
        mesh_sequence : list of Mesh instances
            The entries should have global ids in them.
        """
        self.mesh_sequence = mesh_sequence
        self.collect_all_steps()
        self.calculate_average_cell_area()
        self.generate_death_statistics()
        self.generate_centroid_statistics()
        self.generate_edge_difference_statistics()
        self.generate_tracking_statistics()
        self.generate_rosette_statistics()
        self.output_directory = None
    
    def set_output_directory(self, output_dir):
        """Sets the output dir.
        
        Parameters
        ----------
        
        output_dir : string
        """
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        
        self.output_directory = output_dir

    def write_area_statistics(self):
        """Write the area statistics"""
        area_statistics = []
        for this_mesh in self.mesh_sequence:
            this_area = this_mesh.calculate_total_area()
            this_number_cells = this_mesh.get_num_elements()
            this_average = this_area/this_number_cells
            area_statistics.append( this_average )
        
        area_statistics_np = np.array(area_statistics)
        np.savetxt(os.path.join(self.output_directory, 'area_statistics.csv' ), area_statistics_np)

    def write_rearrangement_statistics(self):
        """Write the area statistics"""
        number_of_rearrangements = []
        for step in self.steps:
            number_of_rearrangements.append( step.number_of_cells_gaining_edges + step.number_of_cells_loosing_edges )
        rearrangement_statistics_np = np.array(number_of_rearrangements)
        np.savetxt(os.path.join(self.output_directory, 'rearrangement_statistics.csv' ), rearrangement_statistics_np)
    
    def write_tracked_cell_statistics(self):
        """Write tracked_cells_statistics"""
        number_of_tracked_cells = []
        number_of_total_cells = []
        these_data = np.zeros( (len(self.steps), 2 ), dtype = 'int')
        for step_counter, step in enumerate( self.steps ):
            these_data[step_counter, 0] = step.mesh_one.get_num_elements()
            these_data[step_counter, 1] = step.number_of_tracked_cells
        np.savetxt(os.path.join(self.output_directory, 'tracking_statistics.csv' ), these_data)
 
    def write_dying_cells(self):
        """make a list of all global ids that are removed"""
        np.savetxt( os.path.join(self.output_directory, 'dying_cells.csv'), 
                    self.global_ids_of_dying_cells )
        
    def write_cell_area_statistics(self):
        """write the area evolution for each global id"""
        maximal_global_id = 0
        for this_mesh in self.mesh_sequence:
            this_max_global_id = this_mesh.get_max_global_id()
            if this_max_global_id > maximal_global_id:
                maximal_global_id = this_max_global_id
            
        cell_area_data = np.zeros( (maximal_global_id + 1, len(self.mesh_sequence)) )
        for mesh_counter, this_mesh in enumerate(self.mesh_sequence):
            for global_id in range(maximal_global_id + 1):
                try:
                    this_element = this_mesh.get_element_with_global_id( global_id )
                    this_area = this_element.calculate_area()
                except KeyError:
                    this_area = np.nan
                cell_area_data[global_id, mesh_counter] = this_area

#Added comma delimiter for reading in rows to make plots.  -Aimee
        np.savetxt(os.path.join(self.output_directory, 'cell_area_statistics.csv' ), cell_area_data, delimiter = ',')
                
    def collect_all_steps(self):
        """Generate StepDataCollectors for each time step"""
        self.steps = []
        for counter, this_mesh in enumerate(self.mesh_sequence):
            if counter > 0:
                previous_mesh = self.mesh_sequence[counter - 1]
                self.steps.append(StepDataCollector(previous_mesh,
                                                    this_mesh,
                                                    counter))
                
    def generate_rosette_statistics(self):
        "Get the total number of rosettes in all meshes"
        self.number_of_rosettes = 0
        for this_mesh in self.mesh_sequence:
            self.number_of_rosettes += this_mesh.count_rosettes()

    def generate_death_statistics(self):
        """Get the total number of dying cells in the sequence"""
        self.number_dying_cells = 0
        self.global_ids_of_dying_cells = []
        for step in self.steps:
            self.number_dying_cells += step.number_dying_cells
            self.global_ids_of_dying_cells += step.global_ids_of_dying_cells
            
    def generate_centroid_statistics(self):
        """Get statistics on centroid displacement"""
        self.centroid_displacements = self.steps[0].centroid_displacements
        for step in self.steps[1:]:
            step.centroid_displacements = np.hstack((self.centroid_displacements,
                                                     step.centroid_displacements))
        
        self.centroid_displacements /= np.sqrt(self.average_cell_area)
        self.maximal_centroid_displacement = np.max(self.centroid_displacements)
        self.minimal_centroid_displacement = np.min(self.centroid_displacements)
        self.average_centroid_displacement = np.mean(self.centroid_displacements)
        
    def calculate_average_cell_area(self):
        "Calculate the average area of all cells of all meshes in the sequence"
        total_area = 0
        total_number_of_cells = 0
        for this_mesh in self.mesh_sequence:
            total_area += this_mesh.calculate_total_area()
            total_number_of_cells += this_mesh.get_num_elements()
        
        self.average_cell_area = total_area/total_number_of_cells
            
    def generate_edge_difference_statistics(self):
        """Collect statistics on how many cells gain vs loose edges in this step"""
        self.number_of_cells_gaining_edges = 0
        self.number_of_cells_loosing_edges = 0
        
        for step in self.steps:
            self.number_of_cells_gaining_edges += step.number_of_cells_gaining_edges
            self.number_of_cells_loosing_edges += step.number_of_cells_loosing_edges

    def generate_tracking_statistics(self):
        """Generate statistics about number of tracked cells"""
        
        shared_global_ids = set(self.mesh_sequence[0].global_id_dictionary.keys())
        for this_mesh in self.mesh_sequence[1:]:
            shared_global_ids.intersection_update(set(this_mesh.global_id_dictionary.keys()))
        
        self.number_of_tracked_cells = len(shared_global_ids)
        self.global_ids_of_tracked_cells = list(shared_global_ids)


class StepDataCollector():
    """A class to analyse two consecutive tracked meshes"""
    def __init__(self, mesh_one, mesh_two, step_number = 0):
        """The constructor of the StepDataCollector
        
        Parameters
        ----------
        
        mesh_one : Mesh instance
            first mesh
            
        mesh_two : Mesh instance
            second_mesh
            
        step_number : int
           number of this step in the sequence
        """
        self.mesh_one = mesh_one
        self.mesh_two = mesh_two
        self.step_number = step_number
        self.generate_tracking_statistics()
        self.generate_death_statistics()
        self.generate_centroid_statistics()
        self.generate_edge_difference_statistics()
        
    def generate_tracking_statistics(self):
        """Generate statistics about number of tracked cells"""
        
        mesh_one_global_ids = self.mesh_one.global_id_dictionary.keys()
        mesh_two_global_ids = self.mesh_two.global_id_dictionary.keys()
        
        shared_global_ids = set.intersection(set(mesh_one_global_ids),
                                             set(mesh_two_global_ids))
        
        self.number_of_tracked_cells = len(shared_global_ids)
        self.global_ids_of_tracked_cells = list(shared_global_ids)

    def generate_death_statistics(self):
        """Collect the number of dying cells in this step
        """
        self.number_dying_cells = 0
        self.global_ids_of_dying_cells = []
        for element in self.mesh_one.elements:
            if element.global_id not in self.mesh_two.global_id_dictionary.keys():
                element_dyed = True
                if element.check_if_on_boundary():
                    element_dyed = False
                else:
                    adjacent_element_ids = element.get_ids_of_adjacent_elements()
                    for frame_id in adjacent_element_ids:
                        adjacent_global_id = self.mesh_one.get_element_with_frame_id(frame_id).global_id
                        if adjacent_global_id not in self.mesh_two.global_id_dictionary.keys():
                            element_dyed = False
                            break
                if element_dyed:
                    self.number_dying_cells +=1
                    self.global_ids_of_dying_cells.append(element.global_id)

    def generate_centroid_statistics(self):
        """Collect statistics on how much centroids move"""
        centroid_displacements = []
        for element in self.mesh_one.elements:
            if element.global_id in self.mesh_two.global_id_dictionary.keys():
                second_element_centroid = self.mesh_two.get_element_with_global_id(element.global_id).calculate_centroid()
                centroid_displacements.append(np.linalg.norm(second_element_centroid - 
                                                             element.calculate_centroid()))
        
        centroid_displacements_np = np.array(centroid_displacements)
        self.centroid_displacements = centroid_displacements_np
        centroid_displacements_rescaled = centroid_displacements_np/np.sqrt(self.mesh_one.calculate_average_element_area())
        self.maximal_centroid_displacement = np.max(centroid_displacements_rescaled)
        self.minimal_centroid_displacement = np.min(centroid_displacements_rescaled)
        self.average_centroid_displacement = np.mean(centroid_displacements_rescaled)
        
    def generate_edge_difference_statistics(self):
        """Collect statistics on how many cells gain vs loose edges in this step"""
        self.number_of_cells_gaining_edges = 0
        self.number_of_cells_loosing_edges = 0
        for element in self.mesh_one.elements:
            if element.global_id in self.mesh_two.global_id_dictionary.keys():
                second_element = self.mesh_two.get_element_with_global_id(element.global_id)
                if element.get_num_nodes() > second_element.get_num_nodes():
                    self.number_of_cells_gaining_edges += 1
                elif element.get_num_nodes() < second_element.get_num_nodes():
                    self.number_of_cells_loosing_edges += 1

class PostProcessor():
    """An object to postprocess a maximum common subgraph and identify rearrangements"""
    def __init__(self, mesh_one, mesh_two, largest_mappings ):
        """The constructor of the post processor
        
        Parameters
        ----------
        
        mesh_one : Mesh instance
            the first frame represented as mesh
            
        mesh_two : Mesh instance
            the second frame represented as mesh
            
        largest_mappings : list of dictionaries
            the list of equivalent largest mappings that the subgraph finder returned
        """
        self.largest_mappings = largest_mappings
        self.mapped_ids = []
        """All currently present global ids"""

        self.mesh_one = mesh_one
        self.network_one = mesh_one.generate_network()

        self.mesh_two = mesh_two
        self.network_two = mesh_two.generate_network()

        self.preliminary_mappings = {}
        """A dictionary of the same style as TrackingState.id_map. Keys are mesh_one frame ids
        and values are mesh_two frame_ids"""
        
    def get_multiple_images( self, list_of_arguments, preliminary_mapping = {} ):
        """Get a list of all images of the given arguments.
        
        Parameters
        ----------
        
        list_of_arguments : list of ints
            list containing frame_ids in mesh_one
            
        preliminary_mapping : dict
            mapping of cells between the two frames for which the global ids have not yet been set
            
        Returns
        -------
        
        list_of_images : list of ints
            list containing all frame_ids in mesh_two of elements that are images of frame_ids
            in list_of_arguments
        """
        
        list_of_images = []
        for frame_id in list_of_arguments:
            global_id = self.mesh_one.get_element_with_frame_id(frame_id).global_id
            if global_id is not None:
                list_of_images.append(self.mesh_two.get_element_with_global_id(global_id).id_in_frame )
            else:
                list_of_images.append(preliminary_mapping[frame_id])
            
        return list_of_images
    
    def post_process_with_data(self):
        """Post process the maximum common subgraph, 'fill in the gaps',
        and return the full list of global ids
        
        Identifies T1 Swaps and maps the involved cells
        
        Returns
        -------
        global_ids : list if ints
            list of all global ids present after post-processing
        """
#         self.index_global_ids_from_largest_mappings()

        network_one = self.mesh_one.generate_network_of_unidentified_elements()
    
        self.stable_fill_in_by_adjacency()

        self.resolve_division_events()
        self.index_global_ids()
        
        return self.mapped_ids
   
    def stable_fill_in_by_adjacency(self):
        """Fill in untracked elements. 
        
        This method sets up a registry of untracked cells and how many tracked neighbours they have.
        This registry is saved under self.connectivity_vector, which safes for each element
        in the first mesh the number of tracked neighbours in the first mesh.
        Based on this registry it will attempt to map cells in a way that maximises the number of preserved
        neighbours upon tracking.
        This is achieved by combining self.connectivity_vector with a boolean vector
        self.actual_connectivy_tested that saves whether the number of preserved neighbours under the best
        possible mapping has been found. The method also uses a current_best_match that has the
        connectivity self.maximal actual connectivity.
        """
        
        self.make_connectivity_vector()

        extension_found_with_relaxed_condition = True
        while extension_found_with_relaxed_condition:
            mapping_has_changed = True
            while mapping_has_changed:
                old_mapping = self.preliminary_mappings.copy()
                self.already_inspected_cells = np.zeros_like(self.connectivity_vector, dtype = 'bool')
                while self.check_mapping_is_extendible():
                    self.maximal_actual_connectivity = 0
                    self.current_best_match = None
                    self.actual_connectivity_tested = np.zeros_like( self.connectivity_vector, dtype = 'bool' )
                    while ( self.get_maximal_connectivity() > self.maximal_actual_connectivity and
                            self.get_maximal_connectivity() > 1 ):
                        next_frame_id = self.pick_next_cell()
                        mapping_candidate, actual_connectivity = self.alternative_find_safe_mapping_candidate_for_single_cell( next_frame_id )
                        element_index = self.mesh_one.frame_id_dictionary[next_frame_id]
                        self.actual_connectivity_tested[element_index] = True
                        if mapping_candidate is not None:
                            if actual_connectivity > self.maximal_actual_connectivity:
                                self.maximal_actual_connectivity = actual_connectivity
                                self.current_best_match = ( next_frame_id, mapping_candidate )
                        else:
                            self.already_inspected_cells[element_index] = True
                    if self.current_best_match is not None:
                        self.extend_preliminary_mapping( self.current_best_match[0], self.current_best_match[1] )
  
                if self.preliminary_mappings == old_mapping:
                    mapping_has_changed = False
                else:
                    mapping_has_changed = True
                  
            self.already_inspected_cells = np.zeros_like(self.connectivity_vector, dtype = 'bool')
            self.maximal_actual_connectivity = 0
            self.current_best_match = None
            self.actual_connectivity_tested = np.zeros_like( self.connectivity_vector, dtype = 'bool' )
            while ( self.get_maximal_connectivity() > self.maximal_actual_connectivity and
                    self.get_maximal_connectivity() > 1 ):
                next_frame_id = self.pick_next_cell()
                mapping_candidate, actual_connectivity = self.alternative_find_safe_mapping_candidate_for_single_cell( next_frame_id, relaxed_condition = True )
                element_index = self.mesh_one.frame_id_dictionary[next_frame_id]
                self.actual_connectivity_tested[element_index] = True
                if mapping_candidate is not None:
                        if actual_connectivity >= 2:
                            self.maximal_actual_connectivity = actual_connectivity
                            self.current_best_match = ( next_frame_id, mapping_candidate )
                else:
                    self.already_inspected_cells[element_index] = True
            if self.current_best_match is not None:
                self.extend_preliminary_mapping( self.current_best_match[0], self.current_best_match[1] )
                extension_not_yet_found = False
                extension_found_with_relaxed_condition = True
            else:
                extension_found_with_relaxed_condition = False
 
    def get_maximal_connectivity(self):
        """Helper method for stable_fill_in_by_adjacency. It it returns
        the maximal connectivity to the mcs among cells that have not yet been inspected
        for actual connectivity, i.e. the possible number of preserved neighbours under the
        best-possible mapping.
        
        Returns
        -------
        
        maximal_connectivity : int
            maximal connectivity among not yet inspected cells.
        """
        not_yet_visited_cells = np.logical_and( self.already_inspected_cells == False, self.actual_connectivity_tested == False )
        maximal_connectivity = np.max( self.connectivity_vector[not_yet_visited_cells])
        return maximal_connectivity

    def pick_next_cell(self):
        """Pick a next cell for inspection for actual connectivity
        
        Returns a cell that has not yet been inspected and for which the actual connectivity has not yet
        been tested.
        
        Returns
        -------
        
        next_cell : int 
            frame id of the cell that is to be inspected next.
        """
        maximal_connectivity = self.get_maximal_connectivity()
        assert(maximal_connectivity > 1)
        not_yet_visited_cells = np.logical_and( self.already_inspected_cells == False, self.actual_connectivity_tested == False )
        possible_indices = np.where( np.logical_and(self.connectivity_vector == maximal_connectivity,
                                                    not_yet_visited_cells ) )
        next_frame_id = self.mesh_one.elements[possible_indices[0][0]].id_in_frame

        return next_frame_id
    
    def check_mapping_is_extendible(self):
        """Returns whether the current mapping is extendible.
        
        Returns True if there are any cells that have not yet been inspected and for which
        the connectivity is larger than one
        
        Returns
        -------
        
        mapping_is_extendible : bool
            True if the mapping is extendible.
        """
        mapping_is_extendible = np.sum(np.logical_and( self.already_inspected_cells == False, 
                                                       self.connectivity_vector > 1 )) > 0

        return mapping_is_extendible
    
    def make_connectivity_vector(self):
        """Make a connectivity vector. The connectivity vector is used
        throughout the method stable_fill_in_by_adjacency. For each cell in the first
        mesh it saves an integer number denoting how many tracked neighbours that cell has. 
        The connectivity vector is stored as a member variable of the post processor.
        """
        connectivity_vector = np.zeros(self.mesh_one.get_num_elements(), dtype = 'int')
        for counter, element in enumerate(self.mesh_one.elements):
            if element.global_id is None:
                full_set_of_currently_mapped_neighbours = self.mesh_one.get_already_mapped_adjacent_element_ids( element.id_in_frame ) 
                connectivity_vector[counter] = len(full_set_of_currently_mapped_neighbours)
            else:
                connectivity_vector[counter] = 0
        
        self.connectivity_vector = connectivity_vector

    def extend_preliminary_mapping(self, next_frame_id, mapping_candidate):
        """Extend the preliminary mapping.
        
        Once stable_fill_in_by_adjacency has found a new mapping this method is called to
        add the mapping to the preliminary mapping.
        
        It will update the connectivity vector for any cells around the cell corresponding
        to next_frame_id and reset their already_inspected vector.
        
        Parameters
        ----------
        
        next_frame_id : int 
            frame id of cell in first mesh that is to be mapped 
        
        mapping_candidate : int
            frame id of cell in second mesh that is to be mapped
        """

        centroid_position = self.mesh_two.get_element_with_frame_id(mapping_candidate).calculate_centroid()
        new_centroid_position = np.array(centroid_position)
        new_centroid_position[1] = 326 - centroid_position[1]

        assert(next_frame_id not in self.preliminary_mappings)
        self.preliminary_mappings[next_frame_id] = mapping_candidate
        new_neighbour_ids = self.mesh_one.get_not_yet_mapped_shared_neighbour_ids( [next_frame_id],
                                                                                   self.preliminary_mappings.keys() )

        element_index = self.mesh_one.frame_id_dictionary[next_frame_id]
        self.connectivity_vector[element_index] = 0
        for neighbour_id in new_neighbour_ids:
            element_index = self.mesh_one.frame_id_dictionary[neighbour_id]
            self.connectivity_vector[element_index] += 1
            self.already_inspected_cells[element_index] = False

    def alternative_find_safe_mapping_candidate_for_single_cell(self, frame_id, relaxed_condition = False ):
        """This method finds a possible mapping candidate for a single cell.
        It is a helper method of stable_fill_in_by_adjacency.
        
        It returns a mapping candidate if the number of gained tracked neighbours is less than
        the number of preserved neighbours - 1. If relaxed_condition is True,
        it returns a mapping candidate if the number of gained tracked neighbours is
        less than the number of preserved tracked neighbours.
        
        Parameters
        ----------
        
        frame_id : int
            integer of the cell for which we try to find a mapping candidate
            
        relaxed_condition : bool
            If True, the number of gained tracked neighbours must be less than
            the number of preserved tracked neighbours. If False, the number
            of gained tracked neighbours must be less than the number
            of preserved tracked neighbours - 1.
            
        Returns
        -------
        
        mapping_candidate : int
            frame id in second mesh that indicates the mapping candidate
            
        current_neighbour_number : int
            number of preserved neighbours
        """

        mapping_candidate = None
        element_one = self.mesh_one.get_element_with_frame_id(frame_id)
        if ( frame_id not in self.preliminary_mappings ):
            full_set_of_currently_mapped_neighbours = self.mesh_one.get_already_mapped_adjacent_element_ids( frame_id, 
                                                                                                             self.preliminary_mappings.keys() )
            # get mapping candidates by all shared neighbours of currently mapped neighbours
            images_of_already_mapped_neighbours = self.get_multiple_images( full_set_of_currently_mapped_neighbours, 
                                                                            self.preliminary_mappings )
            mapping_candidates = self.mesh_two.get_not_yet_mapped_shared_neighbour_ids( images_of_already_mapped_neighbours,
                                                                                        self.preliminary_mappings.values() )
            full_neighbour_number = len( full_set_of_currently_mapped_neighbours )
            current_neighbour_number = len( full_set_of_currently_mapped_neighbours )
            if len(mapping_candidates) == 0:
                mapping_candidates = set()
                old_reduced_image_sets = [images_of_already_mapped_neighbours]
                while ( ( len(mapping_candidates) == 0 ) and 
                        ( current_neighbour_number > 2 ) ):
                    # They don't have a shared neighbour, see whether we can get better mapping candidates if we take one of the
                    # mapped neighbours out to allow for rearrangement 
                    new_reduced_image_sets = []
                    for image_set in old_reduced_image_sets:
                        for image in image_set:
                            reduced_images_of_already_mapped_neighbours = [item for item in image_set
                                                                           if item != image ]
                            mapping_candidates.update( self.mesh_two.get_not_yet_mapped_shared_neighbour_ids( reduced_images_of_already_mapped_neighbours,
                                                                                                              self.preliminary_mappings.values() ))
                            new_reduced_image_sets.append(list(reduced_images_of_already_mapped_neighbours))
                    current_neighbour_number = current_neighbour_number - 1
                    old_reduced_image_sets = list(new_reduced_image_sets)
        
            filtered_mapping_candidates = []
            for candidate in mapping_candidates:
                additional_neighbour_count = self.get_additional_neighbour_count( candidate, images_of_already_mapped_neighbours,
                                                                                  self.preliminary_mappings.values() )
                element_two = self.mesh_two.get_element_with_frame_id(candidate)
                if relaxed_condition:
                    if additional_neighbour_count < full_neighbour_number:
                        filtered_mapping_candidates.append( candidate )
                else:
                    if additional_neighbour_count < full_neighbour_number - 1:
                        filtered_mapping_candidates.append( candidate )

            if len(filtered_mapping_candidates) == 1:
                mapping_candidate = filtered_mapping_candidates[0]
                    
        return mapping_candidate, current_neighbour_number
    
    def find_safe_mapping_candidate_for_single_cell(self, frame_id, preliminary_mapping, min_neighbour_number = 3 ):
        """Finds a mapping candidate for the cell with frame_id
        
        Helper to altered_fill_in_by_adjacency which only gets calles upon division resolution.
        
        Parameters
        ----------
        
        frame_id : int
            frame_id of cell in network one for which a mapping candidate is needed
            
        preliminary_mapping : dict
            existing mappings from network one to network 2
            
        min_neighbour_number : int
            minimal number or connections to already mapped neighbours that the new mapping needs to preserve
        
        Returns
        -------
        
        mapping_candidate : int
            frame_id in network two that has minimal_number_of_connections to already mapped neighbours
            of the element in mesh_one with frame_id. Returns None if no mapping candidate could be found.
        """

        mapping_candidate = None
        # loop over the nodes in the connected component_one
        element_one = self.mesh_one.get_element_with_frame_id(frame_id)
        if ( frame_id not in preliminary_mapping ):
            full_set_of_currently_mapped_neighbours = self.mesh_one.get_already_mapped_adjacent_element_ids( frame_id, 
                                                                                                             preliminary_mapping.keys() )
            if len( full_set_of_currently_mapped_neighbours ) >= min_neighbour_number:
                # get mapping candidates by all shared neighbours of currently mapped neighbours
                images_of_already_mapped_neighbours = self.get_multiple_images( full_set_of_currently_mapped_neighbours, 
                                                                                preliminary_mapping )
                mapping_candidates = self.mesh_two.get_not_yet_mapped_shared_neighbour_ids( images_of_already_mapped_neighbours,
                                                                                            preliminary_mapping.values() )

                if len(mapping_candidates) == 0:
                    mapping_candidates = set()
                    current_neighbour_number = len( full_set_of_currently_mapped_neighbours )
                    old_reduced_image_sets = [images_of_already_mapped_neighbours]
                    while ( len(mapping_candidates) == 0 and current_neighbour_number > min_neighbour_number ):
                        # They don't have a shared neighbour, see whether we can get better mapping candidates if we take one of the
                        # mapped neighbours out to allow for rearrangement 
                        new_reduced_image_sets = []
                        for image_set in old_reduced_image_sets:
                            for image in image_set:
                                reduced_images_of_already_mapped_neighbours = [item for item in image_set
                                                                               if item != image ]
                                assert( len( reduced_images_of_already_mapped_neighbours ) >= min_neighbour_number )
                                mapping_candidates.update( self.mesh_two.get_not_yet_mapped_shared_neighbour_ids( reduced_images_of_already_mapped_neighbours,
                                                                                                                  preliminary_mapping.values() ))
                                new_reduced_image_sets.append(list(reduced_images_of_already_mapped_neighbours))
                        current_neighbour_number = current_neighbour_number - 1
                        old_reduced_image_sets = list(new_reduced_image_sets)
            
                filtered_mapping_candidates = []
                for candidate in mapping_candidates:
                    additional_neighbour_count = self.get_additional_neighbour_count( candidate, images_of_already_mapped_neighbours,
                                                                                      preliminary_mapping.values() )
                    element_two = self.mesh_two.get_element_with_frame_id(candidate)
                    polygon_numbers_add_up = element_two.get_num_nodes() < ( element_one.get_num_nodes() + additional_neighbour_count + 2 )
                    if additional_neighbour_count < 3 and additional_neighbour_count < min_neighbour_number and polygon_numbers_add_up:
                        filtered_mapping_candidates.append( candidate )

                if len(filtered_mapping_candidates) == 1:
                    mapping_candidate = filtered_mapping_candidates[0]
                    
        return mapping_candidate
    
    def get_additional_neighbour_count(self, candidate_id, expected_neighbours, mapped_cells):
        """See how many additional neighbours the cell with candidate_id in mesh_two has (within all already mapped cells).
        
        Parameters
        ----------
        
        candidate_id : int
            id_in_frame of cell in mesh_two
           
        expected_neighbours : list of ints
            cells in mesh two that we expect to be neighbours of candidate
        
        mapped_cells : list of ints
            frame ids in mesh two that have been mapped but whose global ids have not been set
        
        Returns
        -------
        
        additional_neighbour_count : int
            number of mapped neighbours of element with candidate_id that are not in expected_neighbours
        """
        additional_neighbour_count = 0
        candidates_mapped_neighbours = self.mesh_two.get_already_mapped_adjacent_element_ids( candidate_id, 
                                                                                              mapped_cells )
        for neighbour in candidates_mapped_neighbours:
            if neighbour not in expected_neighbours:
                additional_neighbour_count += 1
                
        return additional_neighbour_count
   
    def altered_fill_in_by_adjacency(self, network_one):
        """Fill in unmapped cells by adjacency to existing mapping.
               
        Takes a network of unmapped cells in the first mesh, 
        and fills in the cell-to-cell mapping between them based on adjacency
        with already mapped cells. This method has been replaced by stable_fill_in_by_adjacency
        and is now only used in the division resolution step.
        
        Parameters
        ----------
        
        network_one : networkx Graph instance
            subgraph of the network corresponding to mesh_one
        """

        preliminary_mappings = self.altered_get_mappings_by_adjacency(network_one)
        
        for node in preliminary_mappings:
            self.preliminary_mappings[node] = preliminary_mappings[node]
        
    def altered_get_mappings_by_adjacency(self, connected_component_one):
        """Gets a preliminary mapping based on the adjacency to already mapped nodes.
        
        Helper method for fill_in_by_adjacency and identify_division_event.
        
        Same as altered_fill_in_by_adjacency this method is now only used in the division resolution step

        Parameters
        ----------
        
        connected_component_one : networkx Graph instance
            subgraph of the network corresponding to mesh_one. network of ummapped cells
            
        Returns
        -------
        
        preliminary_mapping : dict
            keys are frame ids in mesh_one, values are frame_ids in mesh_two
        """

        preliminary_mapping = {}
        
        self.extend_current_preliminary_mapping(connected_component_one, preliminary_mapping, minimal_number_of_neighbours=4)
        self.extend_current_preliminary_mapping(connected_component_one, preliminary_mapping, minimal_number_of_neighbours=3)
        self.extend_current_preliminary_mapping(connected_component_one, preliminary_mapping, minimal_number_of_neighbours=2)
#         self.extend_current_preliminary_mapping(connected_component_one, preliminary_mapping, minimal_number_of_neighbours=1)
                
        return preliminary_mapping

    def extend_current_preliminary_mapping(self, network_one, preliminary_mapping, minimal_number_of_neighbours=3):
        """This fills in any unmapped nodes in network one into preliminary mapping, ensuring
           that any new mapping has at least minimal_number_of_neighbours tracked neighbours.
           
        As submethod to altered_fill_in_by_adjacency this method only gets called upon division resolution.
           
        Parameters
        ----------
        
        network_one : networkx.Graph instance
            network of unmapped frame ids in mesh one
            
        preliminary_mapping : dict int->int
            already known mappings from network one
        
        minimal_number_of_neighbours : int
            the minimum number of connections to already mapped cells that the mapping needs to preserve.
        """ 

        attempted_fill_in_counter = {}
        for node in network_one.nodes():
            attempted_fill_in_counter[node] = 0

        not_all_neighbours_mapped = True
        while not_all_neighbours_mapped:
            not_all_neighbours_mapped = False
            for node in network_one.nodes():
                if node not in preliminary_mapping and node not in self.preliminary_mappings:
                    mapping_candidate = self.find_safe_mapping_candidate_for_single_cell( node, preliminary_mapping,
                                                                                          minimal_number_of_neighbours )       
                    if mapping_candidate is not None and mapping_candidate not in preliminary_mapping.values():
                        preliminary_mapping[node] = mapping_candidate
                    else:
                        # this element is still not uniquely identifiable. If all its neighbours have been mapped, then
                        # this means that it actually does not exist in mesh 2, so we stop looking for a match.
                        # otherwise, try again.
                        if len(self.mesh_one.get_element_with_frame_id(node).get_ids_of_adjacent_elements() ) > 2:
                            not_yet_mapped_neighbours = self.mesh_one.get_not_yet_mapped_shared_neighbour_ids([ node ])
                            no_not_yet_mapped_neighbours = 0
                            for neighbour_id in not_yet_mapped_neighbours:
                                if neighbour_id not in preliminary_mapping:
                                    no_not_yet_mapped_neighbours += 1
                            if no_not_yet_mapped_neighbours > 0 and attempted_fill_in_counter[node] < 5:
                                not_all_neighbours_mapped = True
                                attempted_fill_in_counter[node] += 1
            
    def tidy_current_mapping(self):
        """This function resets all global id's that only have one connection to the current maximum common subgraph, or
           two isolated connections, or or members of a small extension to the mcs that contains maximally three cells and 
           has only one connection to the mcs, or connected components of less than ten members.
        """
        isolated_vector = np.zeros( len(self.mesh_one.elements), dtype = 'bool' )
        for element_counter, element in enumerate( self.mesh_one.elements ):
            if element.global_id is not None:
#                 if element.global_id == 166:
#                     import pdb; pdb.set_trace()
                if self.is_isolated( element ):
                    isolated_vector[ element_counter ] = True
                mapped_neighbours = self.mesh_one.get_already_mapped_adjacent_element_ids( element.id_in_frame )
                if len(mapped_neighbours) == 2:
                    first_neighbour_element = self.mesh_one.get_element_with_frame_id( mapped_neighbours[0] )
                    second_neighbour_element = self.mesh_one.get_element_with_frame_id( mapped_neighbours[1] )
                    if self.is_isolated(first_neighbour_element) or self.is_isolated(second_neighbour_element):
                        isolated_vector[element_counter] = True 
                     
        self.remove_global_ids_by_boolean_mask(isolated_vector)
             
        isolated_vector[:] = False
        # Now, let's deal with connected components

        network_one = self.mesh_one.generate_network_of_identified_elements()

        connected_components_in_network_one = nx.connected_components(network_one) 
        
        for connected_component in connected_components_in_network_one:
            if len(connected_component) < 10:
                for frame_id in connected_component:
                    index = self.mesh_one.frame_id_dictionary[frame_id]
                    isolated_vector[index] = True
  
        self.remove_global_ids_by_boolean_mask(isolated_vector)

        self.reindex_global_ids()
#         
        # apply reduced_mcs flags:
        for element in self.mesh_one.elements:
            if element.global_id in self.mapped_ids:
                element.is_in_reduced_mcs_next = True
            else:
                element.is_in_reduced_mcs_next = False

        for element in self.mesh_two.elements:
            if element.global_id in self.mapped_ids:
                element.is_in_reduced_mcs_previous = True
            else:
                element.is_in_reduced_mcs_previous = False
                
    def reindex_global_ids(self):
        """Reindexes the global ids such that the maximal global id corresponds
        to the total number of tracked cells. This method ensures a contiuous count
        of global ids.
        """
        # currently, the mapped ids are not a continuous count, let's change that
        new_mapped_ids = []
        for counter, mapped_id in enumerate(self.mapped_ids):
            self.mesh_one.get_element_with_global_id(mapped_id).global_id = counter
            self.mesh_two.get_element_with_global_id(mapped_id).global_id = counter
            new_mapped_ids.append(counter)
     
        # index the change
        self.mesh_one.index_global_ids()
        self.mesh_two.index_global_ids()
         
        self.mapped_ids = new_mapped_ids

    def remove_global_ids_by_boolean_mask(self, boolean_mask):
        """Remove global ids from all elements for which boolean_map is True
        
        Parameters
        ----------
        
        boolean_map : nd_array, dtype = 'bool'
            mask for elements in the mesh_one elements vector for which we plan to remove the global ids
        """
        for element_counter, element in enumerate( self.mesh_one.elements ):
            if boolean_mask[ element_counter ]:
                this_global_id = element.global_id
                self.mesh_two.get_element_with_global_id(this_global_id).global_id = None
                element.global_id = None   
                del self.largest_mappings[0][element.id_in_frame]
                self.mapped_ids.remove(this_global_id)
                     
        # index the change
        self.mesh_one.index_global_ids()
        self.mesh_two.index_global_ids()

    def is_isolated(self, element):
        """This function determines whether the element is isolated in mesh_one or not.
        
        Parameters
        ----------
        
        element : mesh.Element instance
            a element in a mesh, has to be an element in mesh_one
            
        Returns
        -------
        
        is_isolated : bool
            True if the element is isolated
        """
        
        adjacent_elements = element.get_ids_of_adjacent_elements()
        
        already_mapped_adjacent_elements = []
        for element_id in adjacent_elements:
            if self.mesh_one.get_element_with_frame_id(element_id).global_id is not None:
                already_mapped_adjacent_elements.append(element_id)
        
        if len( already_mapped_adjacent_elements ) == 1 or len(already_mapped_adjacent_elements) == 0:
            is_isolated = True
        elif len( already_mapped_adjacent_elements ) == 2:
            if not self.network_one.has_edge( already_mapped_adjacent_elements[0], already_mapped_adjacent_elements[1]):
                is_isolated = True
            else:
                is_isolated = False
        elif len( already_mapped_adjacent_elements ) == 3:
            number_edges = 0
            if self.network_one.has_edge( already_mapped_adjacent_elements[0], already_mapped_adjacent_elements[1]):
                number_edges+=1
            if self.network_one.has_edge( already_mapped_adjacent_elements[1], already_mapped_adjacent_elements[2]):
                number_edges+=1
            if self.network_one.has_edge( already_mapped_adjacent_elements[0], already_mapped_adjacent_elements[2]):
                number_edges+=1
            if number_edges < 2:
                is_isolated = True
            else:
                is_isolated = False
        else: 
            is_isolated = False
                
        return is_isolated
        
    def index_global_ids(self):
        """add the preliminary mapping to the meshes, i.e. fill in the global ids
        for all mapped cells"""

#         import pdb; pdb.set_trace()
        for element_one_id in self.preliminary_mappings:
            current_maximal_global_id = max( self.mapped_ids )
            new_global_id = current_maximal_global_id + 1
 
            element_one = self.mesh_one.get_element_with_frame_id(element_one_id)
            element_one.global_id = new_global_id
            element_two = self.mesh_two.get_element_with_frame_id(self.preliminary_mappings[element_one_id])
            element_two.global_id = new_global_id

            self.mapped_ids.append(new_global_id)
        
        self.mesh_one.index_global_ids()
        self.mesh_two.index_global_ids()
        
        self.reindex_global_ids()
        
    def index_global_ids_from_largest_mappings(self):
        """Index global ids using all mappings that are contained in all largest mappings"""
        
        preserved_mappings = {}
        
        for key in self.largest_mappings[0]:
            pair_is_in_other_mappings = True
            value = self.largest_mappings[0][key]
            for mapping in self.largest_mappings:
                if key not in mapping:
                    pair_is_in_other_mappings = False
                    break
                elif mapping[key] != value:
                    pair_is_in_other_mappings = False
                    break
            
            if pair_is_in_other_mappings:
                preserved_mappings[key] = value
                
        for global_id, frame_one_id in enumerate(preserved_mappings):
            self.mesh_one.get_element_with_frame_id(frame_one_id).global_id = global_id
            self.mesh_two.get_element_with_frame_id(self.largest_mappings[0][frame_one_id]).global_id = global_id
#             if global_id == 166:
#                 import pdb; pdb.set_trace();
            self.mapped_ids.append(global_id)

        self.mesh_two.index_global_ids()
        self.mesh_one.index_global_ids() 
        
    def identify_division(self, connected_component_one, connected_component_two):
        """Identifies the mother and daughter cells of a division event, and adds the 
        remaining cells to the preliminary mapping.
        
        Parameters
        ----------
         
        connected_component_one : networkx Graph instance
            subgraph of the network corresponding to mesh_one
            
        connected_component_two : networkx Graph instance
            subgraph of the network corresponding to mesh_two
        """
        
        mappings_based_on_adjacency = self.altered_get_mappings_by_adjacency(connected_component_one)
 
#        mappings_based_on_adjacency = self.get_mappings_by_adjacency(connected_component_one, connected_component_two)
         
        bordering_cells_mapping = self.find_bordering_cells_of_division( mappings_based_on_adjacency )
         
        potential_mother_cells = self.mesh_one.get_not_yet_mapped_shared_neighbour_ids( bordering_cells_mapping.keys() )
 
        mother_cell = None
        daughter_cells = None
 
        if len(potential_mother_cells) == 0:
            # In this case one of the daughter cells is triangular.
            # In this case it is not possible to say by adjacency only which cell is the mother cell,
            # Need to make geometric argument
            new_potential_mother_cells = bordering_cells_mapping.keys()
            potential_daughter_cells = bordering_cells_mapping.values() 
            # add the triangular cell
            # this `+' is a list concatenation
            potential_daughter_cells += self.mesh_two.get_not_yet_mapped_shared_neighbour_ids( bordering_cells_mapping.values() ) 
            mother_cell, daughter_cells = self.identify_division_event(new_potential_mother_cells, potential_daughter_cells,
                                                                  mappings_based_on_adjacency)
   
            connected_component_one.remove_node( mother_cell )
            connected_component_two.remove_nodes_from( daughter_cells )
   
            self.altered_fill_in_by_adjacency( connected_component_one )
  
        elif len(potential_mother_cells) == 1:
            potential_mother_cell = potential_mother_cells[0]
            if potential_mother_cell in mappings_based_on_adjacency:
                del mappings_based_on_adjacency[potential_mother_cell]
            for frame_id in mappings_based_on_adjacency:
                self.preliminary_mappings[frame_id] = mappings_based_on_adjacency[frame_id]
        else:
            potential_daughter_cells = self.mesh_two.get_not_yet_mapped_shared_neighbour_ids( bordering_cells_mapping.values() )
#             assert ( len(potential_daughter_cells) > 1)
            if len( potential_daughter_cells ) <= 1 :
                raise Exception("could not resolve division event")
            elif len(potential_daughter_cells) == 3:
                mother_cell, daughter_cells = self.identify_division_event(potential_mother_cells, potential_daughter_cells,
                                                                           mappings_based_on_adjacency)
#   
                connected_component_one.remove_node( mother_cell )
                connected_component_two.remove_nodes_from( daughter_cells )
#   
                self.altered_fill_in_by_adjacency( connected_component_one )
            elif len(potential_daughter_cells) == 4 :
                self.altered_fill_in_by_adjacency( connected_component_one )
            else:
                raise Exception("could not resolve division event")
                 
#             if mother_cell is not None and daughter_cells is not None and daughter_cells != 12:
#                 division_resolved = True
#             else:
#                 division_resolved = False

    def find_bordering_cells_of_division(self, preliminary_mapping):
        """Find the bordering cells of a division in a preliminary mapping. Looks for cells that gain an edge
        in the mapping.
        
        Parameters
        ----------
        
        preliminary_mapping : dict
            keys are frame ids in mesh_one, values are frame_ids in mesh_two. This preliminary mapping must contain the cells
            adjacent to the dividing cell.
        
        Returns
        -------
        
        bordering_cells : dict
            mapping of the cells adjacent to the division
        """
        bordering_cells = {}
        for cell_one in preliminary_mapping:
            num_edges_one = self.mesh_one.get_element_with_frame_id(cell_one).get_num_nodes()
            num_edges_two = self.mesh_two.get_element_with_frame_id(preliminary_mapping[cell_one]).get_num_nodes()
            if num_edges_two == num_edges_one + 1:
                bordering_cells[cell_one] = preliminary_mapping[cell_one]
                
        return bordering_cells
    
    def identify_division_event(self, potential_mother_cells, potential_daughter_cells, preliminary_mapping ):
        """Identify which of the potential mother cells and potential daughter cells are 
        the actual mother and daughter cells of the division
        
        Parameters
        ----------
        
        potential_mother_cells : list
            list of frame ids in mesh_one of potential mother cells
        
        potential_daughter cells : list
            list of frame ids in mesh_two of potential daughter cells
            
        preliminary_mapping : dict
            preliminary mapping that contains at least the two mother cells
        
        Returns
        -------
        
        mother_cell : int
            frame_id of the mother cell in mesh_one
            
        daughter_cells : list
            list containing the frame ids of the two daughter cells of the division
        """
        definite_daughter_cell_set = self.mesh_two.get_inclusive_not_yet_mapped_shared_neighbour_ids(potential_daughter_cells)
        
        # following if statement is to cover case of triangular cells
        if len( definite_daughter_cell_set ) == 1:
            definite_daughter_cell = definite_daughter_cell_set.pop()
        elif len( definite_daughter_cell_set ) == 4:
            # Only one of the provided cells will be triangular
            # if you reached this position in the code
            for frame_id in definite_daughter_cell_set:
                if self.mesh_two.get_element_with_frame_id(frame_id).get_num_nodes() == 3:
                    definite_daughter_cell = frame_id
                    break
        else:
            raise Exception("could not resolve division event")
        
        if definite_daughter_cell is None or definite_daughter_cell == 0 :
            raise Exception("could not resolve division event")

        if len(potential_daughter_cells) <= 1 :
            raise Exception("could not resolve division event")

        potential_daughter_cells.remove( definite_daughter_cell )
        
        inverse_preliminary_mapping = { value: key for key, value in preliminary_mapping.items() }
        
        closest_centroid_distance = sys.float_info.max

        for frame_id in potential_daughter_cells:
            merged_element = self.merge_elements( self.mesh_two.get_element_with_frame_id(definite_daughter_cell),
                                                  self.mesh_two.get_element_with_frame_id(frame_id) )
            merged_centroid = merged_element.calculate_centroid()
            this_mother_cell = self.mesh_one.get_element_with_frame_id(inverse_preliminary_mapping[frame_id])
            this_distance = np.linalg.norm(merged_centroid - this_mother_cell.calculate_centroid())
            if this_distance < closest_centroid_distance:
                definite_mother_cell = this_mother_cell.id_in_frame
                second_definite_daughter_cell = frame_id
                closest_centroid_distance = this_distance

        return definite_mother_cell, [definite_daughter_cell, second_definite_daughter_cell]

    def resolve_division_events(self):
        """Resolve division events.
        
        This method will find all connected components of untracked cells in the second mesh.
        If a connected component is not at the boundary the mothod
        resolve_division_event_for_connected_component is called to attempt to resolve
        the division.
        """
        
#         for frame_one_id in self.largest_mappings[0]:
#             self.preliminary_mappings[frame_one_id] = self.largest_mappings[0][frame_one_id]

#         self.preliminary_mappings = copy.copy(self.largest_mappings[0])
        # first, identify any cells that are in network two but are not mapped
        network_two = self.mesh_two.generate_network_of_unidentified_elements(self.preliminary_mappings.values())
        connected_components_in_network_two = nx.connected_components(network_two) 

        for connected_component in connected_components_in_network_two:
            #check whether component is at mesh boundary:
            component_is_on_boundary = False
            for node in connected_component:
                if self.mesh_two.get_element_with_frame_id(node).check_if_on_boundary():
                    component_is_on_boundary = True
                    break
            if not component_is_on_boundary:
                self.resolve_division_event_for_connected_component(connected_component)

#         self.reindex_global_ids()
        # then, get all their neighbouring cells, and all inverse images of neighbouring cells
        # make a connected component out of both
        # remove both from preliminary mappings
        # identify division event on both connected components
        
    def resolve_division_event_for_connected_component(self, connected_component):
        """This method will extend the connected component in network two by all it's 
        first-order adjacent elements. It will then find the corresponding tracked elements
        to these adjacent elements in the first mesh. It will then construct a connected
        component of the corresponding elements in the first mesh and subsequently
        add any of their shared neighbours.
        
        Finally, it will remove all tracked cells in the first connected component from the
        preliminary mapping and pass both connected components to the method identify_division.
        
        If identify_division fails a warning is given the preliminary mapping is returned to it's
        oroginal state. This means that the preliminar mapping remains unaltered if division resolution
        fails.
        
        Parameters
        ----------
        
        connected_component : list of ints
            list of frame ids of elements in network two that form a connected component.
        """
        # collect_cells_for_connected_component_two
        adjacent_elements = []
        for node in connected_component:
            this_element = self.mesh_two.get_element_with_frame_id(node)
            if not this_element.check_if_on_boundary():
                adjacent_elements += this_element.get_ids_of_adjacent_elements() 
            else:
                print('element to remove is on boundary')
        
        unique_adjacent_elements = np.unique(np.array(adjacent_elements))
        
        preliminary_adjacent_elements = list(set(unique_adjacent_elements).intersection( self.preliminary_mappings.values() ))
        mcs_adjacent_elements = list(set(unique_adjacent_elements).intersection( self.largest_mappings[0].values() ))
        
        # collect cells for connected_component_one
        inverse_preliminary_mapping = { value : key for key, value in self.preliminary_mappings.items() }
        inverse_largest_mapping = { value : key for key, value in self.largest_mappings[0].items() }

        inverse_images_of_preliminary_adjacent_elements = [ inverse_preliminary_mapping[frame_id] for
                                                            frame_id in preliminary_adjacent_elements]
        
        inverse_images_of_mcs_adjacent_elements = [ inverse_largest_mapping[frame_id] for
                                                    frame_id in mcs_adjacent_elements]

        unmapped_elements_belonging_to_connected_component_in_network_one = []

        for element_id in inverse_images_of_preliminary_adjacent_elements:
            unmapped_elements_belonging_to_connected_component_in_network_one += self.mesh_one.get_not_yet_mapped_shared_neighbour_ids([element_id])

        for element_id in inverse_images_of_mcs_adjacent_elements:
            unmapped_elements_belonging_to_connected_component_in_network_one += self.mesh_one.get_not_yet_mapped_shared_neighbour_ids([element_id])

        unmapped_elements_belonging_to_connected_component_in_network_one = list(np.unique(np.array(unmapped_elements_belonging_to_connected_component_in_network_one)))
        
        unmapped_elements_belonging_to_connected_component_in_network_one += inverse_images_of_preliminary_adjacent_elements
        unmapped_elements_belonging_to_connected_component_in_network_one += inverse_images_of_mcs_adjacent_elements
        unmapped_elements_belonging_to_connected_component_in_network_two = [node for node in connected_component] + preliminary_adjacent_elements + mcs_adjacent_elements
        
        # remove the collected cells from the mapping
        old_mappings = dict()
        for frame_id in unmapped_elements_belonging_to_connected_component_in_network_one:
            if frame_id in self.preliminary_mappings:
                old_mappings[frame_id] = self.preliminary_mappings[frame_id]
            elif frame_id in self.largest_mappings[0]:
                old_mappings[frame_id] = self.largest_mappings[0][frame_id]
            global_id = self.mesh_one.get_element_with_frame_id(frame_id).global_id
            try:
                self.mesh_two.get_element_with_global_id(global_id).global_id = None
                self.mapped_ids.remove(global_id)
            except KeyError:
                pass
            self.mesh_one.get_element_with_frame_id(frame_id).global_id = None
            try:
                del( self.preliminary_mappings[frame_id] )
            except KeyError:
                pass
            
        self.mesh_one.index_global_ids()
        self.mesh_two.index_global_ids()
        # make the connected components
        connected_component_one = self.network_one.subgraph( unmapped_elements_belonging_to_connected_component_in_network_one )
        connected_component_two = self.network_one.subgraph( unmapped_elements_belonging_to_connected_component_in_network_two )

        # pass to our connected component function
        try:
            self.identify_division(connected_component_one, connected_component_two)
        except:
            warnings.warn("could not resolve division event")
            for frame_id in old_mappings:
                self.preliminary_mappings[frame_id] = old_mappings[frame_id]

    def merge_elements(self, element_one, element_two):
        """Merge two elements into a bigger element, taking out the shared nodes.
        
        This function will leave the nodes untouched, i.e. their information about elements will not be updated.
        The original elements will also not be affected.
        
        Parameters
        ----------
        
        element_one : Element instance
            first element that we would like to merge

        element_two : Element instance
            second element that we would like to merge
            
        Returns
        -------
        
        merged_element : Element instance
            A new element over the existing nodes. Is not part of the element vectors in the nodes.
        """
        
        new_element_nodes = []
        for local_index, node in enumerate(element_one.nodes):
            if ( element_one.id_in_frame in node.get_adjacent_element_ids() and
                 element_two.id_in_frame in node.get_adjacent_element_ids() ):
                next_node = element_one.nodes[ (local_index + 1)%element_one.get_num_nodes() ]
                if ( element_one.id_in_frame in next_node.get_adjacent_element_ids() and
                     element_two.id_in_frame in next_node.get_adjacent_element_ids() ):
                    new_element_nodes.append(node)
                    one_edge_id = node.id
                    break
                else:
                    previous_node = element_one.nodes[ element_one.get_num_nodes() - 1 ]
                    new_element_nodes.append(previous_node)
                    one_edge_id = previous_node.id
            else:
                new_element_nodes.append(node)
            
        # we find the local index of the found node in the other cell
        for local_index, node in enumerate(element_two.nodes):
            if node.id == one_edge_id:
                second_element_local_index = local_index
                break
            
        # loop through the second element nodes 
        reached_other_side = False
        while reached_other_side == False:
            second_element_local_index = ( second_element_local_index + 1 )%element_two.get_num_nodes()
            next_node = element_two.nodes[second_element_local_index]
            if ( element_one.id_in_frame in next_node.get_adjacent_element_ids() and
                 element_two.id_in_frame in next_node.get_adjacent_element_ids() ):
                new_element_nodes.append(next_node)
                second_edge_id = next_node.id
                reached_other_side = True
            else:
                new_element_nodes.append(next_node)

        # we find the local index of the found node in the other cell
        for local_index, node in enumerate(element_one.nodes):
            if node.id == second_edge_id:
                first_element_local_index = local_index
                break
            
        for local_index in range( first_element_local_index + 1, element_one.get_num_nodes() ):
            new_element_nodes.append(element_one.nodes[local_index])
            
        # We add the nodes to the element after instantiation, so that the element is not added to the node
        merged_element = mesh.Element([])
        merged_element.nodes = new_element_nodes
        
        assert( merged_element.calculate_area() > 0 )
        
        return merged_element
        
def evaluate_tracking(first_mesh, second_mesh, ground_truth):
    """Evaluate the tracking.
    
    Parameters
    ----------
    
    first_mesh : Mesh instance
        this is a mesh that has global ids in it
        
    second_mesh : Mesh instance
        another mesh with global ids in it
        
    ground truth : dictionary, keys and values are integers
        Keys are frame ids in first_mesh, values are 
        frame ids in second_mesh
        
    Returns
    -------
    
    success_boolean : bool
        True if less than four cells in ground_truth are not tracked,
        and if all tracked cells correspond to pairings in ground_truth
        
    number_tracked_cells : int
        Number of correctly tracked cells between first_mesh and
        second_mesh

    Warning
    -------
    
    This function is not tested!
    """   

    correctly_tracked_cells = []
    incorrectly_tracked_cells = []
    missing_cells = []
    for first_element in first_mesh.elements:
        # and that the mapping coincides with the ground truth for all tracked ids
        first_frame_id = first_element.id_in_frame
        if first_frame_id in ground_truth:
            if first_element.global_id is None:
                missing_cells.append(first_frame_id)
            else:
                this_global_id = first_element.global_id
                second_element = second_mesh.get_element_with_global_id(this_global_id)
                second_frame_id = second_element.id_in_frame
                if second_frame_id == ground_truth[first_frame_id]:
                    correctly_tracked_cells.append(first_frame_id)
                else:
                    incorrectly_tracked_cells.append(first_frame_id)

    success_boolean = ( len(missing_cells) < 4 and len(incorrectly_tracked_cells) == 0 )
    number_tracked_cells = len(correctly_tracked_cells)
    number_incorrectly_tracked_cells = len(incorrectly_tracked_cells)
    
    return success_boolean, number_tracked_cells, number_incorrectly_tracked_cells

def find_maximum_common_subgraph(mesh_one, mesh_two):
    """Find a mapping between the cell ids in both frames and assigns the global ids according 
    to their maximum common subgraph.
       
    Writes global_id entries for all identified elements in both meshes.
    
    Parameters
    ----------
    
    mesh_one : Mesh type
        First mesh
    
    mesh_two : Mesh type
        Second mesh
        
    Returns
    -------
    
    mapped_ids : dict (int->int)
        the ids of elements that were identified in both meshes
    """
    
    subgraph_finder = LocalisedSubgraphFinder(mesh_one, mesh_two)
    subgraph_finder.find_maximum_common_subgraph()
    post_processor = PostProcessor(mesh_one, mesh_two, subgraph_finder.largest_mappings)
    post_processor.tidy_current_mapping()
    post_processor.index_global_ids_from_largest_mappings()

    mesh_two.index_global_ids()
    mesh_one.index_global_ids()

    return post_processor.mapped_ids
