# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

"""In this file the main mesh classes 'Mesh', 'Element', and 'Node' are defined.
"""
import glob
import re
from os import path
import os
import numpy as np
import matplotlib.pyplot as plt
import cv2
import pickle
from .core import Mesh, Node, Element
import mahotas


def load(filename):
    """Reads a saved mesh back from a file.
    
    Parameters
    ----------
    
    filename : string
        file that contains the pickled object. The file has to contain a
        pickled Mesh instance.
    
    Returns
    -------
    
    mesh : mesh contained in file
    """
    file_to_read = open(filename, 'rb')
    mesh_to_read = pickle.load( file_to_read )
    file_to_read.close()
    
    return mesh_to_read

def convert_from_ilastik( picture_path, segmentation_path, out_path ):
    """Reads a sequence of ilastik type segmentations and converts them
    into seedwater type segmentation.
    
    Parameters
    ----------
    
    picture_path : string
        path to the images that should be segmented. This folder is expected to contain
        a numbered sequence of .tif files.

    segmentation_path : string
        path to the ilastik sequence. This is expected to be a folder that
        contains a numbered sequence of .tif files that was generated
        using ilastik segmentation from the images in picture_path.
        In particular, the segmented images are expected to have the same order
        and number as the ones in picture_path.
        
    out_path : string
        will store a seedwater-type segmentation that we can then convert to
        a vertex mesh.
    """
    
    list_of_image_files = glob.glob( os.path.join(picture_path , '*.tif') )
    list_of_image_files.sort(key=_natural_keys)
    
    list_of_segmented_files = glob.glob( os.path.join( segmentation_path , '*.tif') )
    list_of_segmented_files.sort(key=_natural_keys)
    
    if not os.path.isdir(out_path):
        os.mkdir(out_path)

    for image_counter, image_path in enumerate( list_of_image_files ):
        image_filename = os.path.split(image_path)[1]
        segmented_image_path = list_of_segmented_files[ image_counter ]
        out_file_name = os.path.join(out_path, 'segmented_' + image_filename)
        convert_single_file_from_ilastik(image_path, segmented_image_path, out_file_name)
        
def convert_single_file_from_ilastik( image_path, segmentation_path, out_path ):
    """Reads a file that was segmented with ilastik and turns it into a
    file as seedwater would create it.
    
    Parameters
    ----------
    
    image_path : string
        path to the image file
        
    segmentation_path : string
        path to the ilastik segmented image file

    out_path : string
        where the converted file should be stored
    """
    full_image = plt.imread(image_path)
    
    segmented_image = plt.imread(segmentation_path)

    seed_image = create_seeds_from_image(segmented_image)

    segmentation = create_segmentation_from_seeds( full_image, seed_image )
    
    cv2.imwrite(out_path, segmentation)

def create_segmentation_from_seeds( input_image, seed_image ):
    """Create a Seedwater style segmentation from a bunch of seeds in an image.
    
    Parameters
    ----------
    
    input_image : nd_array
        image that should be segmented
        
    seed_image : nd_array
        image that contains the seeds
        
    Returns
    -------
    
    segmentation : nd_array
        segmented image represented as 16bit integer image
    """
    segmented_image = mahotas.cwatershed( input_image, seed_image )
    segmented_image_as_16_bit = np.array( segmented_image, dtype = 'uint16' )
    
    return segmented_image_as_16_bit
    
def create_seeds_from_image( input_image ):
    """Create seeds from image to grow full segmentation.
    
    Parameters
    ----------
    
    input_image : nd array
        array representation of the image
        
    Returns
    -------
    seed_image : nd array
        image filled with zeros, each 'seed' has a different integer value
    """
    
    # This will find all the cell-inner parts
    this_image_binary = np.array( (input_image == 2), dtype = 'uint8' )
    
    # distinguish cv2 versions
    opencvversion = float(cv2.__version__.split(".")[0])
    if ( opencvversion< 3.0 or opencvversion>= 4.0):
        these_contours, hirarchy = cv2.findContours(this_image_binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    else:
        _,these_contours, hirarchy = cv2.findContours(this_image_binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    seed_image = np.zeros_like(input_image, dtype = 'uint')
    helper_image = np.zeros_like(input_image, dtype = 'uint8')
    seed_counter = 1
    for contour_counter, contour in enumerate(these_contours):
        helper_image[:] = 0
        cv2.drawContours( helper_image, these_contours, contour_counter, color = 1,
                          thickness = -1 )
        if np.count_nonzero(helper_image) > 1:
            seed_image[ helper_image == 1 ] = seed_counter
            seed_counter +=1

    return seed_image
    
def load_sequence(path):
    """Reads a sequence of meshes
    
    Parameters
    ----------
    
    path : string
       path to the sequence, including the beginning of the filenames of the sequence.
       All files will be read whose file names start with the path string, and ends with `.mesh'
       
    Returns
    -------
    
    mesh_sequence : list of Mesh instances
       The meshes that are saved as indicated by the path parameter
    """
    list_of_files = glob.glob(path + '*.mesh')
    list_of_files.sort(key=_natural_keys)
    mesh_sequence = []
    for filename in list_of_files:
        mesh_sequence.append(load(filename))

    return mesh_sequence

def read_sequence_from_data(folder_name, start_number = 1, number_meshes = None):
    """Reads a seedwater output folder of .tif files and transforms it into a sequence
       of mesh objects
    
    Parameters
    ----------
    
    folder_name : string
        folder that contains the segmented .tif files.    
    
    start_number : int
        number of first file to be considered

    number_meshes : int
        number of last file to be considered

    Returns
    -------
   
    mesh_sequence : list of Mesh instances
        elements of the list are instances of Mesh
    """
   
    list_of_files = get_segmentation_list_in_folder(folder_name, start_number, number_meshes)

    mesh_sequence = []
    for filename in list_of_files:
        print("reading" + str(filename))
        mesh_sequence.append( read_frame_from_data(filename) )
    return mesh_sequence   

def get_segmentation_list_in_folder(folder_name, start_number = 1, number_meshes = None):
    """Get a sorted list of .tif file names in the provided folder
    
    Parameters:
    -----------
    
    folder_name : string
        path to the folder that contains

    start_number : int
        number of first file to be considered

    number_meshes : int
        number of last file to be considered
    
    Returns:
    --------
    
    list_of_files : list of strings
        ordered list of all .tif files in the provided folder
    """
    list_of_files = glob.glob(path.join(folder_name, '*.tif*'))
    list_of_files.sort(key=_natural_keys)

    if number_meshes != None:
        list_of_files = list_of_files[(start_number -1):number_meshes]
    else:
        list_of_files = list_of_files[(start_number -1):]
        
    return list_of_files
 
def read_frame_from_data(filename):
    """Reads a seedwater output segmented .tif files and transforms
       it into a mesh object, discarding boundary cells
    
    Parameters
    ----------
    
    filename : string
        full path of a seedwater segmented .tif file    

    Returns
    -------
    
    cleaned_mesh : Mesh instance
        a mesh representation of this .tif frame
    """
    
    read_mesh = read_data_and_create_raw_mesh( filename )
    read_mesh.merge_short_edges(2)
    read_mesh.remove_boundary_elements()
    
    return read_mesh

def read_data_and_create_raw_mesh( filename ):
    """Reads a seedwater output segmented .tif files and transforms
       it into a mesh object including all fully segmented cells.
    
    Parameters
    ----------
    
    filename : string
        full path of a seedwater segmented .tif file    

    Returns
    -------
    
    this_mesh : Mesh instance
        a mesh representation of this .tif frame
    """

#     this_image = plt.imread(filename)
    this_image = cv2.imread( filename, flags = -1 )

    contour_list, cell_ids = get_contour_list(this_image)
    
    vertex_array, vertices_of_cells, number_vertices_of_cells = extract_vertex_data(this_image, contour_list, cell_ids)
    
    nodes = []
    for vertex_position in vertex_array:
        nodes.append( Node( vertex_position ) )
        
    elements = []
    for cell_index, row in enumerate( vertices_of_cells ):
        element_nodes = []
        for local_index in range( 0, number_vertices_of_cells[cell_index] ):
            element_nodes.append( nodes[row[local_index]] )
        elements.append( Element( element_nodes, cell_ids[ cell_index ] ) )
        
    this_mesh = Mesh( nodes, elements )
    

    
    return this_mesh 

def get_contour_list(this_image):
    """Get a list of contours around each cell in the image.
    
    Uses the opencv function cv2.findContours
    
    Parameters
    ----------
    
    this_image : ndarray
        an image as integer numpy array

    Returns
    -------
    
    contour_list : list
        each entry is a contour as returned by the opencv function
        
    cell_ids : list
        entries are the integer values of segmented objects in the data frame.
        order of this list is the same as in contour_list
    """

    cell_ids = np.unique( this_image )[1:] # Skip background

    # collect contours, i.e. for each region we get the coordinates of the outer-most set of
    # pixels
    contour_list = [None]*len(cell_ids)

    for cell_number, cell_id in enumerate(cell_ids):
        # Get an image with ones for the current cell and zeros outside
        cell_mask = ( this_image == cell_id )
        # Transform to 8 bit for openCV
        cell_mask = np.uint8( cell_mask )
        opencvversion = float(cv2.__version__.split(".")[0])
        if (opencvversion < 3.0 or opencvversion >= 4.0):
            contour_list[cell_number],_ = cv2.findContours( cell_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE )    
        else:
            _,contour_list[cell_number],_ = cv2.findContours( cell_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE )    

    return contour_list, cell_ids
    
def find_triple_junctions_at_pixel(this_image, x_index, y_index):
    """Finds triple junctions in a segmented image around the pixel with coordinates (x_index, y_index)
    
    In this function, a triple junction is indexed by the pixel to the lower left pixel corner of where three
    pixels of different colour meet. Coordinates in the input and output of this function are like this:
    
            _____y______>
           |
          x|   IMAGECONTENT
           |
           \/
           
    Coordinates are different in all other functions in this package: in extract_vertex_data we 
    change to cartesian coordinates.
    
    See Also
    --------
        extract_vertex_data

    Parameters
    ----------
    
    this_image : ndarray
        an image as integer numpy array
        
    x_index : int
        x_coordinate of this pixel
        
    y_index : int 
        y_coordinate of this pixel

    Returns
    -------
    
    odered_triple_junctions: list of [x,y] coordinates
        ordered list of triple junctions, ordered anticlockwise starting
        at the first corner of the 'first' anticlockwise edge. This is the
        anticlockwise edge after an edge that is shared with a pixel of the 
        same colour.
    """

    values_in_neighbourhood = this_image[(x_index -1):(x_index + 2)][:,(y_index - 1):(y_index + 2)]
    central_value = values_in_neighbourhood[1,1]
    triple_junctions = []

    if ( central_value != values_in_neighbourhood[0,0] and central_value != values_in_neighbourhood[0,1] and
         values_in_neighbourhood[0,0] != values_in_neighbourhood [0,1] ): # top left

        triple_junctions.append( np.array([x_index, y_index -1]) )

    if ( central_value != values_in_neighbourhood[0,0] and central_value != values_in_neighbourhood[1,0] and
         values_in_neighbourhood[0,0] != values_in_neighbourhood [1,0] ): # top left

        triple_junctions.append( np.array([x_index, y_index -1]) )

    if values_in_neighbourhood.shape[0] == 3 and values_in_neighbourhood.shape[1] ==3:
        if ( central_value != values_in_neighbourhood[1,0] and central_value != values_in_neighbourhood[2,0] and
             values_in_neighbourhood[1,0] != values_in_neighbourhood [2,0] ): # bottom left

            triple_junctions.append( np.array([x_index +1, y_index -1]) )

    if values_in_neighbourhood.shape[0] == 3:
        if ( central_value != values_in_neighbourhood[2,0] and central_value != values_in_neighbourhood[2,1] and
             values_in_neighbourhood[2,0] != values_in_neighbourhood [2,1] ): # bottom left

            triple_junctions.append( np.array([x_index +1, y_index -1]) )

        if values_in_neighbourhood.shape[1] == 3:
            if ( central_value != values_in_neighbourhood[2,1] and central_value != values_in_neighbourhood[2,2] and
                 values_in_neighbourhood[2,1] != values_in_neighbourhood [2,2] ): # bottom right

                triple_junctions.append( np.array([x_index +1, y_index]) )

            if ( central_value != values_in_neighbourhood[2,2] and central_value != values_in_neighbourhood[1,2] and
                 values_in_neighbourhood[2,2] != values_in_neighbourhood [1,2] ): # bottom right

                triple_junctions.append( np.array([x_index +1, y_index]) )

    if values_in_neighbourhood.shape[1] == 3:
        if ( central_value != values_in_neighbourhood[0,1] and central_value != values_in_neighbourhood[0,2] and
             values_in_neighbourhood[0,1] != values_in_neighbourhood [0,2]): # top right

            triple_junctions.append( np.array([x_index, y_index]) )

        if ( central_value != values_in_neighbourhood[0,2] and central_value != values_in_neighbourhood[1,2] and
             values_in_neighbourhood[0,2] != values_in_neighbourhood [1,2]): # top right

            triple_junctions.append( np.array([x_index, y_index]) )

    ### ORDERING corners here
        
    ordered_triple_junctions = []
    triple_junctions = np.array(triple_junctions)

    if ( len( triple_junctions ) > 1 ):
        
        if ( central_value == values_in_neighbourhood[0,1] and central_value != values_in_neighbourhood[1,0] ):

            if np.any(np.all( triple_junctions == [x_index, y_index -1] , axis = 1)):
                ordered_triple_junctions.append( [x_index, y_index -1] )
            if np.any(np.all( triple_junctions == [x_index +1, y_index -1] , axis = 1)):
                ordered_triple_junctions.append( [x_index +1, y_index -1])
            if np.any(np.all( triple_junctions == [x_index +1, y_index] , axis = 1)):
                ordered_triple_junctions.append( [x_index +1, y_index])
            if np.any(np.all( triple_junctions == [x_index, y_index] , axis = 1)):
                ordered_triple_junctions.append( [x_index, y_index])

        if values_in_neighbourhood.shape[0] == 3:
            if central_value == values_in_neighbourhood[1,0] and central_value != values_in_neighbourhood[2,1]:

                if np.any(np.all( triple_junctions == [x_index +1, y_index -1] , axis = 1)):
                    ordered_triple_junctions.append( [x_index +1, y_index -1])
                if np.any(np.all( triple_junctions == [x_index +1, y_index] , axis = 1)):
                    ordered_triple_junctions.append( [x_index +1, y_index])
                if np.any(np.all( triple_junctions == [x_index, y_index] , axis = 1)):
                    ordered_triple_junctions.append( [x_index, y_index])
                if np.any(np.all( triple_junctions == [x_index, y_index -1] , axis = 1)):
                    ordered_triple_junctions.append( [x_index, y_index -1] )

        if values_in_neighbourhood.shape[0] == 3 and values_in_neighbourhood.shape[1] == 3:
            if central_value == values_in_neighbourhood[2,1] and central_value != values_in_neighbourhood[1,2]:
    
                if np.any(np.all( triple_junctions == [x_index +1, y_index] , axis = 1)):
                    ordered_triple_junctions.append( [x_index +1, y_index])
                if np.any(np.all( triple_junctions == [x_index, y_index] , axis = 1)):
                    ordered_triple_junctions.append( [x_index, y_index])
                if np.any(np.all( triple_junctions == [x_index, y_index -1] , axis = 1)):
                    ordered_triple_junctions.append( [x_index, y_index -1] )
                if np.any(np.all( triple_junctions == [x_index +1, y_index -1] , axis = 1)):
                    ordered_triple_junctions.append( [x_index +1, y_index -1])

        if values_in_neighbourhood.shape[1] == 3:
            if central_value == values_in_neighbourhood[1,2] and central_value != values_in_neighbourhood[0,1]:
    
                if np.any(np.all( triple_junctions == [x_index, y_index] , axis = 1)):
                    ordered_triple_junctions.append( [x_index, y_index])
                if np.any(np.all( triple_junctions == [x_index, y_index -1] , axis = 1)):
                    ordered_triple_junctions.append( [x_index, y_index -1] )
                if np.any(np.all( triple_junctions == [x_index +1, y_index -1] , axis = 1)):
                    ordered_triple_junctions.append( [x_index +1, y_index -1])
                if np.any(np.all( triple_junctions == [x_index +1, y_index] , axis = 1)):
                    ordered_triple_junctions.append( [x_index +1, y_index])
    else:
        ordered_triple_junctions = triple_junctions
        
    return ordered_triple_junctions

def extract_vertex_data(this_image, contour_list, cell_ids):
    """Extract location of vertices and associate vertices with cells.
    
    Vertex coordinates are returned in this coordinate system
    
           /\
           |
          x|   IMAGECONTENT
           |_____________>
                 y
                 
    Parameters
    ----------
    
    this_image : ndarray
        an image as integer numpy array

    contour_list : list
        each entry is a contour as returned by the opencv function,
        i.e. each entry an anticlockwise list of integer pixel coordinates
        
    cell_ids : list
        integer identifiers of cells from the segmentation
            
    Returns
    -------
    
    vertex_array : list of vertex coordinates
        each entry is a [x,y] pair of floating point values
        denoting the coordinate of where 3 pixels of different colours meet in the
        segmented image.
    
    vertices_of_cells : ndarray
        each row contains indices for vertex_array and is a anti-clockwise list of 
        vertices for each cell with integer coordinates, with trailing zeros. The 
        number of vertices that need to be considered is indexed in number_vertices_of_cells.
    
    number_vertices_of_cells : list of ints
        same order as vertices_of_cells, each entry indexes how many values are to be considered
        in vertices_of_cells.
        
    cell_ids : list of ints
        same order as vertices_of_cells and number_vertices_of_cells. Each entry is the integer
        identifier if this cell as provided by the segmentation.
    """
    
    no_of_cells = len(cell_ids)
    # We overestimate the size of this array by approximately a factor of 2
    number_of_vertices_estimate = int(no_of_cells*10)
    vertex_array = np.zeros( ( number_of_vertices_estimate, 2 ), dtype=np.int64 )
    no_of_vertices = 0
    no_of_cells_per_vertex = np.zeros(number_of_vertices_estimate, dtype=np.int64)

    # We also overestimate the maximal number of vertices that a cell can have
    vertices_of_cells = np.zeros( (no_of_cells,1000), dtype=np.int64)
    number_vertices_of_cells = np.zeros(no_of_cells, dtype=np.int64)

    # We loop over each cell
    for contour_index in range(no_of_cells):
        no_of_vertices_this_cell = 0
        # We loop over each border pixel of that cell
#         for pixel_index, pixel_coordinates in enumerate( contour_list[contour_index] ):
# 
#             y_index = pixel_coordinates[pixel_index][0][0]
#             x_index = pixel_coordinates[pixel_index][0][1]

        for contour in contour_list[contour_index]:
            for pixel_coordinates in contour:
 
                y_index = pixel_coordinates[0][0]
                x_index = pixel_coordinates[0][1]

                # Get all the surrounding pixels of the current pixel
                values_in_neighbourhood = this_image[(x_index -1):(x_index + 2)][:,(y_index - 1):(y_index + 2)]
                no_values_in_neighbourhood = len(np.unique(values_in_neighbourhood))

                if (no_values_in_neighbourhood > 2):
                    ordered_triple_junctions = find_triple_junctions_at_pixel(this_image, x_index, y_index)
                    for triple_junction in ordered_triple_junctions:
                        if np.any( np.all( vertex_array == triple_junction, axis =1 ) ):
                            #This_junction has been found before, here is the index of that vertex
                            vertex_index = np.where( np.all( vertex_array == triple_junction, axis = 1 ))[0][0]
                            # The Triple Junction is potentially in the Neighbourhood of multiple pixels of a cell
                            # Check whether it is member of the cell already
                            if no_of_vertices_this_cell > 0:
                                if ( vertex_index not in 
                                     vertices_of_cells[contour_index,0:number_vertices_of_cells[contour_index]]):
                                    # Ok, it was not in the neighbourhood of the last pixel, so we can
                                    # assign the Vertex to this cell
                                    no_of_cells_per_vertex[ vertex_index ]  += 1
                                    vertices_of_cells[contour_index, no_of_vertices_this_cell] = vertex_index                                
                                    number_vertices_of_cells[contour_index] += 1
                                    no_of_vertices_this_cell += 1
                            else:
                                # This is the first vertex for this cell, so let's assign the vertex to
                                # this cell
                                no_of_cells_per_vertex[ vertex_index ] += 1
                                vertices_of_cells[contour_index, no_of_vertices_this_cell] = vertex_index
                                number_vertices_of_cells[contour_index] += 1
                                no_of_vertices_this_cell += 1
                        else:
                            # This is a new vertex. Let's write that vertex and assign it to this cell
                            vertex_array[ no_of_vertices ] = triple_junction
                            no_of_cells_per_vertex[no_of_vertices] +=1

                            vertices_of_cells[ contour_index, no_of_vertices_this_cell] = no_of_vertices
                            no_of_vertices_this_cell += 1
                            number_vertices_of_cells[contour_index] += 1

                            no_of_vertices +=1
    
    # We now have some vectors that are too long. Let's shorten them!
    # We also need move the vertex positions from the pixel to the actual
    # junction
    
    no_of_cells_per_vertex = no_of_cells_per_vertex[0:no_of_vertices]
    vertices_of_cells = vertices_of_cells[0:no_of_vertices]
    
    # adjust coordinates of actual vertex positions
    vertex_array = np.float64(vertex_array[0:no_of_vertices])
    vertex_array = vertex_array + [-0.5,0.5]

    # move to Cartesian coordinates
    new_vertex_array = np.zeros_like(vertex_array)
    new_vertex_array[:,0] = vertex_array[:,1]
    new_vertex_array[:,1] = len(this_image[:,0]) - vertex_array[:,0]
#     new_vertex_array = np.hstack(( vertex_array[:,1,np.newaxis], vertex_array[:,0,np.newaxis] ))
#     new_vertex_array[:,1] *= -1
#     new_vertex_array[:,1] -= ( np.min(new_vertex_array[:,1]) + np.max(new_vertex_array[:,1]) )

    return new_vertex_array, vertices_of_cells, number_vertices_of_cells
 
def _atoi(text):
    """Helper method for sorting in read_sequence_from_data"""
    return int(text) if text.isdigit() else text

def _natural_keys(text):
    """Helper method for sorting in read_sequence_from_data

    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    """

    return [ _atoi(c) for c in re.split('(\d+)', text) ]
