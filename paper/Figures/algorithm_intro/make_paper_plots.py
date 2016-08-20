# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

import mesh
import tracking
import copy
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from os import path
from os.path import dirname

def is_isolated(sample_mesh, element):
    adjacent_elements = element.get_ids_of_adjacent_elements()
    already_mapped_adjacent_elements = []
    for element_id in adjacent_elements:
        if sample_mesh.get_element_with_frame_id(element_id).global_id != None:
            already_mapped_adjacent_elements.append(element_id)
    
    if len( already_mapped_adjacent_elements ) == 1 or len(already_mapped_adjacent_elements) == 0:
        is_isolated = True
    elif len( already_mapped_adjacent_elements ) == 2:
        if not already_mapped_adjacent_elements[1] in sample_mesh.get_element_with_frame_id(already_mapped_adjacent_elements[0]).get_ids_of_adjacent_elements():
            is_isolated = True
        else:
            is_isolated = False
    else: 
        is_isolated = False
            
    return is_isolated

def find_weakly_connected_cells(sample_mesh):
    
    weakly_connected_global_ids = []
    for element in sample_mesh.elements:
        if element.global_id != None:
            if is_isolated(sample_mesh, element):
                this_global_id = element.global_id
                element.global_id = None   
                weakly_connected_global_ids.append(this_global_id)
    
    return weakly_connected_global_ids

first_filename = path.join(dirname(__file__),'..','..','..','test','test_on_data','data','first_few_frames', 'Segment_0_000.tif')
second_filename = path.join(dirname(__file__),'..','..','..','test','test_on_data','data','first_few_frames', 'Segment_0_001.tif')
third_filename = path.join(dirname(__file__),'..','..','..','test','test_on_data','data','first_few_frames', 'Segment_0_002.tif')


# first_mcs = mesh.load(path.join(dirname(__file__),'..','..','..','test','test_on_data','output','first_mesh_tracked.mesh'))
# second_mcs = mesh.load(path.join(dirname(__file__),'..','..','..','test','test_on_data','output','second_mesh_tracked.mesh'))
# 
# 
# first_mcs_cleaned = mesh.load(path.join(dirname(__file__),'..','..','..','test','test_on_data','output','first_mesh_tracked_cleaned.mesh'))
# second_mcs_cleaned = mesh.load(path.join(dirname(__file__),'..','..','..','test','test_on_data','output','second_mesh_tracked_cleaned.mesh'))
# 
# 
# first_mcs_tracked = mesh.load(path.join(dirname(__file__),'..','..','..','test','test_on_data','output','first_mesh_tracked_processed.mesh'))
# second_mcs_tracked = mesh.load(path.join(dirname(__file__),'..','..','..','test','test_on_data','output','second_mesh_tracked_processed.mesh'))

first_mcs = mesh.load(path.join(dirname(__file__),'..','..','..','test','test_on_data','output','second_mesh_tracked_with_third_mesh.mesh'))
second_mcs = mesh.load(path.join(dirname(__file__),'..','..','..','test','test_on_data','output','third_mesh_tracked.mesh'))

sys.setrecursionlimit(40000)
first_mcs_copy = copy.deepcopy(first_mcs)

weakly_connected_global_ids = find_weakly_connected_cells(first_mcs_copy)

first_mcs_cleaned = mesh.load(path.join(dirname(__file__),'..','..','..','test','test_on_data','output','second_mesh_tracked_with_third_mesh_cleaned.mesh'))
second_mcs_cleaned = mesh.load(path.join(dirname(__file__),'..','..','..','test','test_on_data','output','third_mesh_tracked_cleaned.mesh'))


first_mcs_tracked = mesh.load(path.join(dirname(__file__),'..','..','..','test','test_on_data','output','second_mesh_tracked_with_third_mesh_tracked_processed.mesh'))
second_mcs_tracked = mesh.load(path.join(dirname(__file__),'..','..','..','test','test_on_data','output','third_mesh_tracked_processed.mesh'))


first_mesh = mesh.read_frame_from_data(first_filename)
second_mesh = mesh.read_frame_from_data(second_filename)
third_mesh = mesh.read_frame_from_data(third_filename)

plt.close('all')

# First wireframe

all_edges_by_elements = first_mcs.collect_elements_of_inner_edges()

first_wireframe = first_mcs.get_polygon_collection()
first_wireframe.set_facecolor('None')
first_wireframe.set_edgecolor('black')
first_wireframe.set_linewidth(0.2)
first_wireframe_figure = plt.figure()
for edge_elements in all_edges_by_elements:
    first_centroid = first_mcs.get_element_with_frame_id(edge_elements.pop()).calculate_centroid()
    second_centroid = first_mcs.get_element_with_frame_id(edge_elements.pop()).calculate_centroid()
    plt.plot( [ first_centroid[0], second_centroid[0]],
              [ first_centroid[1], second_centroid[1]], color = 'black', solid_capstyle = 'round',
              linewidth=2 )
first_wireframe_figure.gca().add_collection(first_wireframe)
first_wireframe_figure.gca().set_aspect('equal')
first_wireframe_figure.gca().autoscale_view()
plt.axis('off')
first_wireframe_figure.savefig('first_wireframe.pdf', bbox_inches = 'tight')
plt.close(first_wireframe_figure)
 
# Second wireframe
all_edges_by_elements = second_mcs.collect_elements_of_inner_edges()
second_wireframe = second_mcs.get_polygon_collection()
second_wireframe.set_facecolor('None')
second_wireframe.set_edgecolor('black')
second_wireframe.set_linewidth(0.2)
second_wireframe_figure = plt.figure()
for edge_elements in all_edges_by_elements:
    first_centroid = second_mcs.get_element_with_frame_id(edge_elements.pop()).calculate_centroid()
    second_centroid = second_mcs.get_element_with_frame_id(edge_elements.pop()).calculate_centroid()
    plt.plot( [ first_centroid[0], second_centroid[0]],
              [ first_centroid[1], second_centroid[1]], color = 'black', solid_capstyle = 'round',
              linewidth=2 )
second_wireframe_figure.gca().add_collection(second_wireframe)
second_wireframe_figure.gca().set_aspect('equal')
second_wireframe_figure.gca().autoscale_view()
plt.axis('off')
second_wireframe_figure.savefig('second_wireframe.pdf', bbox_inches = 'tight')
plt.close(second_wireframe_figure)
 
# First MCS
first_mcs_polygon_list = []
for element in first_mcs.elements:
    this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                       fill = True)
    if element.global_id == None:
        this_polygon.set_facecolor([1.0, 1.0, 1.0])
    elif element.global_id in weakly_connected_global_ids:
            this_polygon.set_facecolor('purple')
    else:
        this_polygon.set_facecolor('#2dec34')
    first_mcs_polygon_list.append(this_polygon)
first_mcs_polygon_collection = mpl.collections.PatchCollection(first_mcs_polygon_list, match_original = True)
first_mcs_polygon_collection.set_linewidth(2)
first_mcs_figure = plt.figure()
first_mcs_figure.gca().add_collection(first_mcs_polygon_collection)
first_mcs_figure.gca().set_aspect('equal')
first_mcs_figure.gca().autoscale_view()
plt.axis('off')
first_mcs_figure.savefig('first_mcs.pdf', bbox_inches = 'tight')
plt.close(first_mcs_figure)

# Second MCS
second_mcs_polygon_list = []
for element in second_mcs.elements:
    print 'hello, I am an element'
    this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                       fill = True)
    if element.global_id == None:
        this_polygon.set_facecolor([1.0, 1.0, 1.0])
    elif element.global_id in weakly_connected_global_ids:
        this_polygon.set_facecolor('purple')
    else:
        this_polygon.set_facecolor('#2dec34')
    second_mcs_polygon_list.append(this_polygon)
second_mcs_polygon_collection = mpl.collections.PatchCollection(second_mcs_polygon_list, match_original = True)
second_mcs_polygon_collection.set_linewidth(2)
second_mcs_figure = plt.figure()
second_mcs_figure.gca().add_collection(second_mcs_polygon_collection)
second_mcs_figure.gca().set_aspect('equal')
second_mcs_figure.gca().autoscale_view()
plt.axis('off')
second_mcs_figure.savefig('second_mcs.pdf', bbox_inches = 'tight')
plt.close(second_mcs_figure)

# First MCS cleaned
first_mcs_cleaned_polygon_list = []
for element in first_mcs_cleaned.elements:
    this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                       fill = True)
    if element.global_id == None:
        this_polygon.set_facecolor([1.0, 1.0, 1.0])
    else:
        this_polygon.set_facecolor('#2dec34')
    first_mcs_cleaned_polygon_list.append(this_polygon)
first_mcs_cleaned_polygon_collection = mpl.collections.PatchCollection(first_mcs_cleaned_polygon_list, match_original = True)
first_mcs_cleaned_polygon_collection.set_linewidth(2)
first_mcs_cleaned_figure = plt.figure()
first_mcs_cleaned_figure.gca().add_collection(first_mcs_cleaned_polygon_collection)
first_mcs_cleaned_figure.gca().set_aspect('equal')
first_mcs_cleaned_figure.gca().autoscale_view()
plt.axis('off')
first_mcs_cleaned_figure.savefig('first_mcs_cleaned.pdf', bbox_inches = 'tight')
plt.close(first_mcs_cleaned_figure)
 
# Second MCS cleaned
second_mcs_cleaned_polygon_list = []
for element in second_mcs_cleaned.elements:
    this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                       fill = True)
    if element.global_id == None:
        this_polygon.set_facecolor([1.0, 1.0, 1.0])
    else:
        this_polygon.set_facecolor('#2dec34')
    second_mcs_cleaned_polygon_list.append(this_polygon)
second_mcs_cleaned_polygon_collection = mpl.collections.PatchCollection(second_mcs_cleaned_polygon_list, match_original = True)
second_mcs_cleaned_polygon_collection.set_linewidth(2)
second_mcs_cleaned_figure = plt.figure()
second_mcs_cleaned_figure.gca().add_collection(second_mcs_cleaned_polygon_collection)
second_mcs_cleaned_figure.gca().set_aspect('equal')
second_mcs_cleaned_figure.gca().autoscale_view()
plt.axis('off')
second_mcs_cleaned_figure.savefig('second_mcs_cleaned.pdf', bbox_inches = 'tight')
plt.close(second_mcs_cleaned_figure)
 
# First MCS tracked
first_mcs_tracked_polygon_list = []
for element in first_mcs_tracked.elements:
    this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                       fill = True)
    if element.global_id == None:
        this_polygon.set_facecolor([1.0, 1.0, 1.0])
    else:
        this_polygon.set_facecolor('#2dec34')
    first_mcs_tracked_polygon_list.append(this_polygon)
first_mcs_tracked_polygon_collection = mpl.collections.PatchCollection(first_mcs_tracked_polygon_list, match_original = True)
first_mcs_tracked_polygon_collection.set_linewidth(2)
first_mcs_tracked_figure = plt.figure()
first_mcs_tracked_figure.gca().add_collection(first_mcs_tracked_polygon_collection)
first_mcs_tracked_figure.gca().set_aspect('equal')
first_mcs_tracked_figure.gca().autoscale_view()
plt.axis('off')
first_mcs_tracked_figure.savefig('first_mcs_tracked.pdf', bbox_inches = 'tight')
plt.close(first_mcs_tracked_figure)
 
# Second MCS tracked
second_mcs_tracked_polygon_list = []
for element in second_mcs_tracked.elements:
    this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                       fill = True)
    if element.global_id == None:
        this_polygon.set_facecolor([1.0, 1.0, 1.0])
    else:
        this_polygon.set_facecolor('#2dec34')
    second_mcs_tracked_polygon_list.append(this_polygon)
second_mcs_tracked_polygon_collection = mpl.collections.PatchCollection(second_mcs_tracked_polygon_list, match_original = True)
second_mcs_tracked_polygon_collection.set_linewidth(2)
second_mcs_tracked_figure = plt.figure()
second_mcs_tracked_figure.gca().add_collection(second_mcs_tracked_polygon_collection)
second_mcs_tracked_figure.gca().set_aspect('equal')
second_mcs_tracked_figure.gca().autoscale_view()
plt.axis('off')
second_mcs_tracked_figure.savefig('second_mcs_tracked.pdf', bbox_inches = 'tight')
plt.close(second_mcs_tracked_figure)

