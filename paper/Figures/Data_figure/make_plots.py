# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/LICENSE.

import mesh
from mesh.core import _get_distinct_colors
import tracking
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
import random
from os import path
from os.path import dirname
import numpy as np
import colorsys

first_filename = path.join(dirname(__file__),'..','..','..','test','test_on_data','data','first_few_frames', 'Segment_0_000.tif')
second_filename = path.join(dirname(__file__),'..','..','..','test','test_on_data','data','first_few_frames', 'Segment_0_001.tif')
third_filename = path.join(dirname(__file__),'..','..','..','test','test_on_data','data','first_few_frames', 'Segment_0_002.tif')

first_mesh = mesh.read_frame_from_data(first_filename)
second_mesh = mesh.read_frame_from_data(second_filename)
third_mesh = mesh.read_frame_from_data(third_filename)

first_image = plt.imread(first_filename)
second_image = plt.imread(second_filename)

colored_first_image = np.ones( (first_image.shape[0], first_image.shape[1], 3))
colored_second_image = np.ones( (second_image.shape[0], second_image.shape[1], 3))

num_distinct_colors = 10
HSV_tuples = [(color_index*(1.0/num_distinct_colors), 0.75,0.8) for color_index in range(num_distinct_colors)]
distinct_colors = map(lambda color: colorsys.hsv_to_rgb(*color), HSV_tuples)

set_of_all_colors = set(range(1, num_distinct_colors +1))

first_color_map = np.zeros(len(first_mesh.elements))
first_frame_id_dict = {}
for counter, element in enumerate(first_mesh.elements):
    first_frame_id_dict[element.id_in_frame] = counter

for counter, element in enumerate(first_mesh.elements):
    adjacent_frame_ids = element.get_ids_of_adjacent_elements()
    adjacent_colors = []
    for adjacent_frame_id in adjacent_frame_ids:
        adjacent_colors.append(first_color_map[first_frame_id_dict[adjacent_frame_id]])
    remaining_adjacent_colors = set_of_all_colors.difference(set(adjacent_colors))
    if len(remaining_adjacent_colors) > 0:
        this_color = random.sample(remaining_adjacent_colors, 1)[0]
    else:
        this_color = num_distinct_colors
    first_color_map[counter] = this_color
    colored_first_image[first_image == element.id_in_frame] = distinct_colors[this_color-1]

second_frame_id_dict = {}
second_color_map = np.zeros(len(second_mesh.elements))
for counter, element in enumerate(second_mesh.elements):
    second_frame_id_dict[element.id_in_frame] = counter

for counter, element in enumerate(second_mesh.elements):
    adjacent_frame_ids = element.get_ids_of_adjacent_elements()
    adjacent_colors = []
    for adjacent_frame_id in adjacent_frame_ids:
        adjacent_colors.append(second_color_map[second_frame_id_dict[adjacent_frame_id]])
    remaining_adjacent_colors = set_of_all_colors.difference(set(adjacent_colors))
    if len(remaining_adjacent_colors) > 0:
        this_color = random.sample(remaining_adjacent_colors, 1)[0]
    else:
        this_color = num_distinct_colors
    second_color_map[counter] = this_color
    colored_second_image[second_image == element.id_in_frame] = distinct_colors[this_color -1]

# First Figure
figuresize = (4,2.75)

node_positions_list = []
for node in first_mesh.nodes:
    node_positions_list.append(node.position)

node_positions = np.array(node_positions_list)

node_positions[:,1] *= -1
node_positions[:,1] += len(colored_first_image[:,0,0])
#transform vertex coordinates to image coordinates

first_polygon_collection = first_mesh.get_polygon_collection()

for path in first_polygon_collection.properties()['paths']:
    new_vertices = np.zeros_like(path.vertices)
    path.vertices[:,1] *= -1
    path.vertices[:,1] += len(colored_first_image[:,0,0])

first_polygon_collection.set_facecolor('None')
first_polygon_collection.set_edgecolor('black')
first_polygon_collection.set_alpha(1.0)
first_polygon_collection.set_linewidth(3)

first_figure = plt.figure()
plt.imshow(colored_first_image)
plt.axis('off')
first_figure.gca().add_collection(first_polygon_collection)
first_figure.savefig('first_data_segmented.pdf', bbox_inches = 'tight')

# Second Figure

node_positions_list = []
for node in first_mesh.nodes:
    node_positions_list.append(node.position)

node_positions = np.array(node_positions_list)

node_positions[:,1] *= -1
node_positions[:,1] += len(colored_second_image[:,0,0])
#transform vertex coordinates to image coordinates

second_polygon_collection = second_mesh.get_polygon_collection()

for path in second_polygon_collection.properties()['paths']:
    new_vertices = np.zeros_like(path.vertices)
    path.vertices[:,1] *= -1
    path.vertices[:,1] += len(colored_second_image[:,0,0])

second_polygon_collection.set_facecolor('None')
second_polygon_collection.set_edgecolor('black')
second_polygon_collection.set_alpha(1.0)
second_polygon_collection.set_linewidth(1)

second_figure = plt.figure()
plt.imshow(colored_second_image)
plt.axis('off')
second_figure.gca().add_collection(second_polygon_collection)
second_figure.savefig('second_data_segmented.pdf', bbox_inches = 'tight')

plt.close('all')
# First wireframe
first_wireframe = first_mesh.get_polygon_collection()
first_wireframe.set_facecolor('None')
first_wireframe.set_edgecolor('black')
first_wireframe.set_linewidth(1)
first_wireframe_figure = plt.figure()
first_wireframe_figure.gca().add_collection(first_wireframe)
first_wireframe_figure.gca().set_aspect('equal')
first_wireframe_figure.gca().autoscale_view()
plt.axis('off')
first_wireframe_figure.savefig('first_wireframe.pdf', bbox_inches = 'tight')

# Second wireframe
second_wireframe = second_mesh.get_polygon_collection()
second_wireframe.set_facecolor('None')
second_wireframe.set_edgecolor('black')
second_wireframe.set_linewidth(1)
second_wireframe_figure = plt.figure()
second_wireframe_figure.gca().add_collection(second_wireframe)
second_wireframe_figure.gca().set_aspect('equal')
second_wireframe_figure.gca().autoscale_view()
plt.axis('off')
second_wireframe_figure.savefig('second_wireframe.pdf', bbox_inches = 'tight')
