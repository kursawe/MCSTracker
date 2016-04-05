# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/LICENSE.

import mesh
import tracking
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
from os import path
from os.path import dirname
import numpy as np

data_collector = tracking.analyse_tracked_sequence(path.join(dirname(__file__),'..','..','..','test','test_on_data',
                                                             'output','first_few_frames'))
                                                   
# First MCS
first_mcs_tracked = data_collector.mesh_sequence[0]
first_mcs_polygon_list = []
first_mcs_centroids = []
for element in first_mcs_tracked.elements:
    this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                       fill = True)
    if element.global_id in data_collector.global_ids_of_tracked_cells:
        this_polygon.set_facecolor('#2dec34')
        first_mcs_centroids.append(element.calculate_centroid())
    else:
        if element.global_id in data_collector.global_ids_of_dying_cells:
            this_polygon.set_facecolor('black')
        else:
            this_polygon.set_facecolor([1.0, 1.0, 1.0])
    first_mcs_polygon_list.append(this_polygon)
first_mcs_polygon_collection = mpl.collections.PatchCollection(first_mcs_polygon_list, match_original = True)
first_mcs_figure = plt.figure()
first_mcs_figure.gca().add_collection(first_mcs_polygon_collection)
first_mcs_figure.gca().set_aspect('equal')
first_mcs_figure.gca().autoscale_view()
plt.axis('off')
first_mcs_figure.savefig('first_reduced_tracking.pdf', bbox_inches = 'tight')
first_mcs_centroids_np = np.array(first_mcs_centroids)
plt.scatter(first_mcs_centroids_np[:,0], first_mcs_centroids_np[:,1],
            color = 'yellow', edgecolors = 'black', s=45, lw=1.5)
first_mcs_figure.savefig('first_reduced_tracking_w_centroids.pdf', bbox_inches = 'tight')
plt.close(first_mcs_figure)

# Second MCS
second_mcs_tracked = data_collector.mesh_sequence[1]
second_mcs_polygon_list = []
second_mcs_centroids = []
for element in second_mcs_tracked.elements:
    this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                       fill = True)
    if element.global_id in data_collector.global_ids_of_tracked_cells:
        previous_centroid = first_mcs_tracked.get_element_with_global_id(element.global_id).calculate_centroid()
        number_contained_points = 0
        for sample_centroid in first_mcs_centroids_np:
            if this_polygon.contains_point ( sample_centroid ):
                number_contained_points += 1 
        if ( this_polygon.contains_point( previous_centroid ) and
             number_contained_points == 1):
            this_polygon.set_facecolor('#2dec34')
        else: 
            this_polygon.set_facecolor('darkmagenta')
        second_mcs_centroids.append(element.calculate_centroid())
    else:
        this_polygon.set_facecolor([1.0, 1.0, 1.0])
    second_mcs_polygon_list.append(this_polygon)
second_mcs_polygon_collection = mpl.collections.PatchCollection(second_mcs_polygon_list, match_original = True)
second_mcs_figure = plt.figure()
second_mcs_figure.gca().add_collection(second_mcs_polygon_collection)
second_mcs_figure.gca().set_aspect('equal')
second_mcs_figure.gca().autoscale_view()
plt.axis('off')
second_mcs_figure.savefig('second_reduced_tracking.pdf', bbox_inches = 'tight')
# plt.scatter(first_mcs_centroids_np[:,0], first_mcs_centroids_np[:,1], lw = 0)
plt.scatter(first_mcs_centroids_np[:,0], first_mcs_centroids_np[:,1], 
            color = 'yellow', edgecolors = 'black', s=45, lw=1.5)
second_mcs_figure.savefig('second_reduced_tracking_w_centroids.pdf', bbox_inches = 'tight')
plt.close(second_mcs_figure)
second_mcs_centroids_np = np.array(second_mcs_centroids)

# Third MCS
third_mcs_tracked = data_collector.mesh_sequence[2]
third_mcs_polygon_list = []
for element in third_mcs_tracked.elements:
    this_polygon = mpl.patches.Polygon([node.position for node in element.nodes],
                                       fill = True)
    if element.global_id in data_collector.global_ids_of_tracked_cells:
        previous_centroid = second_mcs_tracked.get_element_with_global_id(element.global_id).calculate_centroid()
        number_contained_points = 0
        for sample_centroid in first_mcs_centroids_np:
            if this_polygon.contains_point ( sample_centroid ):
                number_contained_points += 1 
        if ( this_polygon.contains_point( previous_centroid ) and
             number_contained_points == 1):
            this_polygon.set_facecolor('#2dec34')
        else: 
            this_polygon.set_facecolor('darkmagenta')
        second_mcs_centroids.append(element.calculate_centroid())
    else:
        this_polygon.set_facecolor([1.0, 1.0, 1.0])
    third_mcs_polygon_list.append(this_polygon)
third_mcs_polygon_collection = mpl.collections.PatchCollection(third_mcs_polygon_list, match_original = True)
third_mcs_figure = plt.figure()
third_mcs_figure.gca().add_collection(third_mcs_polygon_collection)
third_mcs_figure.gca().set_aspect('equal')
third_mcs_figure.gca().autoscale_view()
plt.axis('off')
third_mcs_figure.savefig('third_reduced_tracking.pdf', bbox_inches = 'tight')
plt.scatter(second_mcs_centroids_np[:,0], second_mcs_centroids_np[:,1],
            color = 'yellow', edgecolors = 'black', s=45, lw=1.5)
third_mcs_figure.savefig('third_reduced_tracking_w_centroids.pdf', bbox_inches = 'tight')
plt.close(third_mcs_figure)
