#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 23:37:06 2020

@author: miladrv
"""
import pandas as pd
import numpy as np
import csv
import numpy.matlib
from scipy.spatial import distance
# from __future__ import print_function

###############################################################################
############################### Inputs ########################################
###############################################################################

# Number of cell types presented in data
noClasses = 15

# Path to ground truth graph path
# groundtruth_graph_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/OCT16_calprotectin/Results_IL6/Results_0_2_5/IL6_IL6R_and_IL6_IL6ST/graph_data/Links.csv'
# resulting_graph_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/NOV4_cellphonedb/0225/Links.csv'
resulting_graph_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/NOV5_BH/IL6/02375/IL6_IL6R_and_IL6_IL6ST/graph_data/Links.csv'

# Path to resulting graph paths
groundtruth_graph_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/OCT16_calprotectin/Results_IL6/ground_truth/Links_subunit.csv'


###############################################################################
############################### Outputs #######################################
###############################################################################

# Path to the intersection graph between the input graphs
output_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/NOV4_cellphonedb/0175/'

###############################################################################
###############################################################################
###############################################################################

graph_a = np.array(pd.read_csv(groundtruth_graph_path, header = 0, index_col = None))
graph_b = np.array(pd.read_csv(resulting_graph_path, header = 0, index_col = None))

# calculate the kl divergence
def kl_divergence(p, q):
	return sum(p[i] * np.log2(p[i]/q[i]) for i in range(len(p)))

distance_vec_a = []
distance_vec_b = []
for i in range(1,noClasses+1):
    for j in range(noClasses+1, 2*noClasses+1):
        index_location_a_x = np.argwhere(graph_a[:,0] == i)
        index_location_a_y = np.argwhere(graph_a[:,1] == j)
        index_intersect_a = np.intersect1d(index_location_a_x ,index_location_a_y)
        common_graph_a_gene_names = list(graph_a[index_intersect_a,3])
        common_graph_a_weights = list(graph_a[index_intersect_a,2])
        
        index_location_b_x = np.argwhere(graph_b[:,0] == i)
        index_location_b_y = np.argwhere(graph_b[:,1] == j)
        index_intersect_b = np.intersect1d(index_location_b_x ,index_location_b_y)
        common_graph_b_gene_names = list(graph_b[index_intersect_b,3])
        common_graph_b_weights = list(graph_b[index_intersect_b,2])
        
        for k in range(len(index_intersect_b)):
            if common_graph_b_gene_names[k] in common_graph_a_gene_names:
                index_element = common_graph_a_gene_names.index(common_graph_b_gene_names[k])
                # distance = f(common_graph_a_weights[index_element],common_graph_b_weights[k])
                distance_vec_a.append(common_graph_a_weights[index_element])
                distance_vec_b.append(common_graph_b_weights[k])
         

## Save results as csv file
scores_graph = pd.concat([pd.DataFrame(distance_vec_a), 
                     pd.DataFrame(distance_vec_b)],
                     axis = 1)
with open(output_path + 'Common_graph.csv','w',newline = '') as f:
     thewriter = csv.writer(f, delimiter = ',')
     thewriter.writerow( ('GT', 'Outcome') )
     # thewriter = csv.writer(f, dialect='exce')
     for row in range(scores_graph.shape[0]):
         thewriter.writerow([scores_graph.iloc[row,0], scores_graph.iloc[row,1]])
print()
###############################################################################
########################## Cell type-Cell type interaction score ##############
###############################################################################

celltype_vec_a = np.zeros(noClasses ** 2)
celltype_vec_b = np.zeros(noClasses ** 2)
points_size = np.zeros(noClasses ** 2)
counter = 0
for i in range(1,noClasses+1):
    for j in range(noClasses+1, 2*noClasses+1):
        index_location_a_x = np.argwhere(graph_a[:,0] == i)
        index_location_a_y = np.argwhere(graph_a[:,1] == j)
        index_intersect_a = np.intersect1d(index_location_a_x ,index_location_a_y)
        celltype_vec_a[counter] = len(index_intersect_a)

        
        index_location_b_x = np.argwhere(graph_b[:,0] == i)
        index_location_b_y = np.argwhere(graph_b[:,1] == j)
        index_intersect_b = np.intersect1d(index_location_b_x ,index_location_b_y)
        common_graph_b_gene_names = list(graph_b[index_intersect_b,3])
        common_graph_b_weights = list(graph_b[index_intersect_b,2])
        celltype_vec_b[counter] = len(index_intersect_b)
        counter = counter + 1


## Calculate_size_of points

## Save results as csv file
scores_cellType = pd.concat([pd.DataFrame(celltype_vec_a), 
                     pd.DataFrame(celltype_vec_b)],
                     axis = 1)
with open(output_path + 'Common_heatmap.csv','w',newline = '') as f:
     thewriter = csv.writer(f, delimiter = ',')
     thewriter.writerow( ('GT', 'Outcome') )
     # thewriter = csv.writer(f, dialect='exce')
     for row in range(scores_cellType.shape[0]):
         thewriter.writerow([scores_cellType.iloc[row,0], scores_cellType.iloc[row,1]])
print()

intersect = np.min(np.transpose(scores_cellType))
sensitivity = np.sum(intersect)/np.sum(scores_cellType.iloc[:,0])
precision = np.sum(intersect)/np.sum(scores_cellType.iloc[:,1])
accuracy = 0.5*(precision + sensitivity)

s = f"""
{'-'*40}
# Prediction scores
# Accuracy: {accuracy}
# Precision: {precision}
# Sensitivity: {sensitivity}

{'-'*40}
"""
print(s)
