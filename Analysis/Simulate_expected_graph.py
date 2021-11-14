#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 12:15:25 2020

@author: miladrv
"""

import pandas as pd
import numpy as np
import scipy.stats
import random 
import time
import csv
import numpy.matlib
import matplotlib.pyplot as plt
# from __future__ import print_function
from scipy.io import mmread
import h5py

###############################################################################
############################### Inputs #######################################
###############################################################################
    
# Path to input gene expression matrix
# scRNA_path = '/cloud-data/cloud-pipeline-milad-storage/IL6_data/matrix.mtx'
# scRNA_path = '/cloud-data/cloud-pipeline-milad-storage/TNF_data/matrix.mtx'
scRNA_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/For_Naouel/RA_Data/New_data/matrix.mtx'


# Path to names of all genes in gene expression matrix
gene_names_all_path = '/cloud-data/cloud-pipeline-milad-storage/IL6_data/genes.txt'
# gene_names_all_path = '/cloud-data/cloud-pipeline-milad-storage/TNF_data/genes.txt'
# gene_names_all_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/For_Naouel/RA_Data/New_data/genes.txt'

# Path to desired genes (Ligand and receptors of the desired pathway)
gene_names_path = "/cloud-data/cloud-pipeline-milad-storage/IL6_data/Selected_gene_names.csv"
# gene_names_path = "/cloud-data/cloud-pipeline-milad-storage/TNF_data/Selected_gene_names.csv"
# gene_names_path = "/cloud-data/cloud-pipeline-milad-storage/CATS/For_Naouel/RA_Data/New_data/Selected_gene_names.csv"

# Path to selected interations (which correspond to the desired signaling pathway)
selected_interactions_path = '/cloud-data/cloud-pipeline-milad-storage/IL6_data/selected_interactions_richa.csv'
# selected_interactions_path = '/cloud-data/cloud-pipeline-milad-storage/TNF_data/selected_interactions_richa.csv'
# selected_interactions_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/For_Naouel/RA_Data/New_data/selected_interactions_richa.csv'

# Desired pair of interaction
desired_interactions = [0,1]

# Path to metadata
# metadata_path = '/cloud-data/cloud-pipeline-milad-storage/IL6_data/categorical_coloring_data.json'
# metadata_path = '/cloud-data/cloud-pipeline-milad-storage/TNF_data/categorical_coloring_data.json'
metadata_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/For_Naouel/RA_Data/New_data/categorical_coloring_data.json'

# Path to selected ligand-(first) subunit interations 
subUnit_interactions_path1 = '/cloud-data/cloud-pipeline-milad-storage/IL6_data/Sub_unit_interactions_richa1.csv'
# subUnit_interactions_path1 = '/cloud-data/cloud-pipeline-milad-storage/TNF_data/Sub_unit_interactions_richa1.csv'
# subUnit_interactions_path1 = '/cloud-data/cloud-pipeline-milad-storage/CATS/For_Naouel/RA_Data/New_data/Sub_unit_interactions_richa1.csv'

# Desired pair of subunit interaction 1
desired_interactions_subUnit_1 = [0]

# Path to selected ligand-(second) subunit interations 
subUnit_interactions_path2 = '/cloud-data/cloud-pipeline-milad-storage/IL6_data/Sub_unit_interactions_richa2.csv'
# subUnit_interactions_path2 = '/cloud-data/cloud-pipeline-milad-storage/TNF_data/Sub_unit_interactions_richa2.csv'
# subUnit_interactions_path2 = '/cloud-data/cloud-pipeline-milad-storage/CATS/For_Naouel/RA_Data/New_data/Sub_unit_interactions_richa2.csv'

# Desired pair of subunit interaction 2
desired_interactions_subUnit_2 = [0]

threshold_numbers = 0.1

threshold_expressions = 0.1
# Threshold on p-value for filtering non-significant LR interactions
# threshold = 0.05

output_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/NOV4_cellphonedb/Groundtruth/01/'
###############################################################################
############################### Outputs #######################################
###############################################################################

# heatmap_map_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/Simulated_Groundtruth/IL6/Heatmap_subunit.csv'
# heatmap_map_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/Simulated_Groundtruth/TNF/Heatmap_subunit.csv'
heatmap_map_path = output_path + 'Heatmap_subunit.csv'

# graph_links_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/Simulated_Groundtruth/IL6/Links_subunit.csv'
# graph_links_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/Simulated_Groundtruth/TNF/Links_subunit.csv'
graph_links_path = output_path + 'Links_subunit.csv'

# graph_nodes_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/Simulated_Groundtruth/IL6/Nodes_subunit.csv'
# graph_nodes_path = '/cloud-data/cloud-pipeline-milad-storage/CATS/Simulated_Groundtruth/TNF/Nodes_subunit.csv'
graph_nodes_path = output_path + 'Nodes_subunit.csv'

###############################################################################
########################## Read scRNA-seq data ################################
###############################################################################
def filter_out_unannotated_cells(scRNA_path,
                                 metadata_path,
                                 select_track,
                                 exclude_tracks,
                                 subset):
    
    metadata = pd.read_json(metadata_path)
    cell_type_labels = metadata.loc['label_list', select_track]
    cell_type_keep = [i[0] for i in enumerate(cell_type_labels) if i[1] not in exclude_tracks]
    # if subset: 
    #     subset_track = list(subset.keys())[0]
    #     cell_type_track = metadata.loc['label_list',  subset_track]
    #     subset_cats  = subset[subset_track]
        
    #     not_in_system = set(subset_cats) - set(cell_type_track)
    #     if len(not_in_system) > 0:
    #         print(f'WARNING! {not_in_system} ARE NOT IN {subset_track}')
            
    #     cell_type_keep = [i[0] for i, j in zip(enumerate(cell_type_labels), cell_type_track)
    #                       if i[1] not in exclude_tracks and j in subset_cats]
        
    cell_type_labels = [i[1] for i in enumerate(cell_type_labels) if i[0] in cell_type_keep]
    scRNA_data = mmread(scRNA_path)
    scRNA_array = scRNA_data.toarray()
    scRNA_array = scRNA_array[cell_type_keep, ]
    scRNA_array = np.transpose(scRNA_array)
    
    return scRNA_array, cell_type_labels

# Select groups in the metadata track for exclusion, default="Unclassified, Other":
exclude_tracks = "Unclassified, Other"

# Select metadata track for analysis, default="CellStates":
select_track = "CellStates"

subset = True
    
scRNA_array, cell_type_labels = filter_out_unannotated_cells(scRNA_path,
                             metadata_path, select_track, 
                             exclude_tracks, subset)

# scRNA_data = mmread(scRNA_path)
# scRNA_array = scRNA_data.toarray()
# scRNA_array = np.transpose(scRNA_array)
noCells = scRNA_array.shape[1]
###############################################################################
############################## selected genes #################################
###############################################################################

gene_names_all = pd.read_csv(gene_names_all_path, header = None, index_col = None)
gene_names = pd.read_csv(gene_names_path, header = None, index_col = None)
gene_number = len(gene_names)
expressions_scRNA = np.zeros([gene_number,noCells])
for k in range(gene_number):
    gene_index = gene_names_all.isin([gene_names.iloc[k,:][0]])
    seriesObj = gene_index.any(axis = 1) 
    rowNames = list(seriesObj[seriesObj == True].index) 
    expressions_scRNA[k,:] = scRNA_array[rowNames,:]

noCells = expressions_scRNA.shape[1]
# expressions_scRNA = expressions_scRNA.values.astype(float)
expressions_scRNA = np.log2(expressions_scRNA + 1)
###############################################################################
############################## Meta data analysis #############################
###############################################################################
# metadata = pd.read_json(metadata_path)
# cell_type_labels = metadata.loc['label_list','CellStates']
unique_cell_type_labels = np.unique(cell_type_labels)
noClasses = len(unique_cell_type_labels)
cell_type_labels_df = pd.DataFrame(cell_type_labels)

###############################################################################
####################### Average expression matrix #############################
###############################################################################
average_expression = np.zeros([gene_number,noClasses])
cell_type_numbers = [0]*noClasses
for k in range(noClasses):
    cell_type_index = cell_type_labels_df.isin([unique_cell_type_labels[k]])
    seriesObj = cell_type_index.any(axis = 1) 
    columnNames = list(seriesObj[seriesObj == True].index) 
    cell_type_numbers[k] = len(columnNames)
    for i in range(gene_number):
        average_expression[i,k] = np.mean(expressions_scRNA[i,columnNames])

# average_expression = np.transpose(average_expression)
average_expression_df = pd.DataFrame(average_expression, index = gene_names, columns = unique_cell_type_labels)
 
p_values =  np.zeros([gene_number,noClasses]) 
for i in range(gene_number):
    for j in range(noClasses):
        boolean_index_selected =  average_expression_df.columns != average_expression_df.columns[j]
        t, p = scipy.stats.ttest_1samp(np.array(average_expression_df.iloc[i, boolean_index_selected]), average_expression_df.iloc[i,j]) 
        if t < 0:
           p_values[i,j] = p


###############################################################################
################## Prepare data for plotting ##################################
###############################################################################           

genes_vec = ['a'] * gene_number * noClasses
counter = 0
for i in range(gene_number):
    for j in range(noClasses):
        genes_vec[counter] = gene_names.iloc[i][0]
        counter = counter + 1
        
clusters_vec = ['a'] * gene_number * noClasses
cluster_names = pd.DataFrame(unique_cell_type_labels)
counter = 0
for i in range(gene_number):
    for j in range(noClasses):
        clusters_vec[counter] = cluster_names.iloc[j][0]
        counter = counter + 1

p_values_vec = np.zeros(gene_number * noClasses)
counter = 0
for i in range(gene_number):
    for j in range(noClasses):
        p_values_vec[counter] = p_values[i,j]
        counter = counter + 1

average_expression_vec = np.zeros(gene_number * noClasses)
counter = 0
for i in range(gene_number):
    for j in range(noClasses):
        average_expression_vec[counter] = average_expression[i,j]
        counter = counter + 1
        
results = pd.concat([pd.DataFrame(genes_vec), 
                     pd.DataFrame(clusters_vec),
                     pd.DataFrame(p_values_vec), 
                     pd.DataFrame(average_expression_vec)],
                     axis = 1)

# with open('/Users/miladrv/Desktop/CellPhoneDB/NOV5/GroundTruth/IL6/deg_analysis_results.csv','w',newline = '') as f:
#      thewriter = csv.writer(f, delimiter = ',')
#      thewriter.writerow( ('genes', 'clusters', 'pvalue', 'mean') )
#      # thewriter = csv.writer(f, dialect='exce')
#      for row in range(results.shape[0]):
#          thewriter.writerow([results.iloc[row,0], results.iloc[row,1], results.iloc[row,2], results.iloc[row,3]])
# print()

###############################################################################
########################## Graph plot #########################################
###############################################################################
node_a_list = []
node_b_list = []
link_list = []
link_weight_list_expressions = []
link_weight_list_pvalues = []
for i_c in range(noClasses):
    for j_c in range(noClasses):
        for i_g in range(gene_number):
            for j_g in range(gene_number):
                # if p_values[i_g,i_c] > 0 and p_values[j_g,j_c] > 0:
                if cell_type_specific_numbers[i_g,i_c] > threshold_numbers and cell_type_specific_expressions[i_g,i_c] >  threshold_expressions and cell_type_specific_numbers[j_g,j_c] > threshold_numbers and cell_type_specific_expressions[j_g,j_c] > threshold_expressions: 
                   node_a_list.append(i_c + 1)
                   node_b_list.append(noClasses + j_c + 1)
                   link_weight_list_expressions.append((average_expression[i_g,i_c] * average_expression[j_g,j_c]))
                   link_weight_list_pvalues.append(0.5 * (-np.log10(p_values[i_g,i_c]) + -np.log10(p_values[j_g,j_c])))
                   link_list.append(gene_names.iloc[i_g][0] + ' ' + '|' + ' ' + gene_names.iloc[j_g][0])


graph_data_links_expressions = pd.concat([pd.DataFrame(node_a_list), 
                     pd.DataFrame(node_b_list),
                     pd.DataFrame(link_weight_list_expressions), 
                     pd.DataFrame(link_list)],
                     axis = 1)

graph_data_links_pvalues = pd.concat([pd.DataFrame(node_a_list), 
                     pd.DataFrame(node_b_list),
                     pd.DataFrame(link_weight_list_pvalues), 
                     pd.DataFrame(link_list)],
                     axis = 1)              

node_data_id = list(range(1,2*noClasses+1))
node_data_type = [2]*len(unique_cell_type_labels) + [1]*len(unique_cell_type_labels) 
node_data_name =  list(unique_cell_type_labels) + list(unique_cell_type_labels) 
x = list(range(1,noClasses*5+1,5)) + list(range(1,noClasses*5+1,5))
y = node_data_type    
graph_data_nodes = pd.concat([pd.DataFrame(node_data_id),
                     pd.DataFrame(node_data_name), 
                     pd.DataFrame(node_data_type),
                     pd.DataFrame(y),
                     pd.DataFrame(x)],
                     axis = 1)                
with open(graph_nodes_path,'w',newline = '') as f:
     thewriter = csv.writer(f, delimiter = ',')
     thewriter.writerow( ('id', 'name', 'type', 'y', 'x') )
     # thewriter = csv.writer(f, dialect='exce')
     for row in range(graph_data_nodes.shape[0]):
         thewriter.writerow([graph_data_nodes.iloc[row,0], graph_data_nodes.iloc[row,1], graph_data_nodes.iloc[row,2], graph_data_nodes.iloc[row,3], graph_data_nodes.iloc[row,4]])
print() 

###############################################################################
######################### Filtering based on (Richa LR) #######################
###############################################################################
selected_interactions = pd.read_csv(selected_interactions_path, header = None, index_col = None)
selected_interactions = selected_interactions.iloc[desired_interactions]
filter_interactions = graph_data_links_expressions.iloc[:,3].isin(selected_interactions.iloc[:,0])
index_selected = list(filter_interactions[filter_interactions == True].index)
graph_data_links_pvalues_richa = graph_data_links_pvalues.iloc[index_selected, :]
graph_data_links_expressions_richa = graph_data_links_expressions.iloc[index_selected, :]

# with open('/Users/miladrv/Desktop/CellPhoneDB/NOV5/GroundTruth/IL6/Links_pvalues.csv','w',newline = '') as f:
#      thewriter = csv.writer(f, delimiter = ',')
#      thewriter.writerow( ('from', 'to', 'weight', 'name') )
#      # thewriter = csv.writer(f, dialect='exce')
#      for row in range(graph_data_links_pvalues_richa.shape[0]):
#          thewriter.writerow([graph_data_links_pvalues_richa.iloc[row,0], graph_data_links_pvalues_richa.iloc[row,1], graph_data_links_pvalues_richa.iloc[row,2], graph_data_links_pvalues_richa.iloc[row,3]])
# print()

with open(graph_links_path,'w',newline = '') as f:
     thewriter = csv.writer(f, delimiter = ',')
     thewriter.writerow( ('from', 'to', 'weight', 'name') )
     # thewriter = csv.writer(f, dialect='exce')
     for row in range(graph_data_links_expressions_richa.shape[0]):
         thewriter.writerow([graph_data_links_expressions_richa.iloc[row,0], graph_data_links_expressions_richa.iloc[row,1], graph_data_links_expressions_richa.iloc[row,2], graph_data_links_expressions_richa.iloc[row,3]])
print()

##############################################################################
######################## Cell-Cell interaction heatmap data ##################
##############################################################################
heatmap = np.zeros([noClasses, noClasses])
count = 0
for i in range(1,noClasses+1):
    for j in range(noClasses+1, 2*noClasses+1):
        index_location_x = np.argwhere(np.asarray(graph_data_links_expressions_richa.iloc[:,0] == i))
        index_location_y = np.argwhere(np.asarray(graph_data_links_expressions_richa.iloc[:,1] == j))
        index_intersect = np.intersect1d(index_location_x ,index_location_y)
        heatmap[i - 1,j - (noClasses+1)] = len(index_intersect)
        count = count + 1
heatmap_df = pd.DataFrame(heatmap, index = unique_cell_type_labels.tolist(), columns = unique_cell_type_labels.tolist())
# heatmap_df.to_csv('/Users/miladrv/Desktop/CellPhoneDB/NOV5/GroundTruth/IL6/heatmap.csv')

##############################################################################
######################## Sub unit analysis ###################################
##############################################################################
subUnit_interactions_data1 = pd.read_csv(subUnit_interactions_path1, header = None, index_col = None)
subUnit_interactions1 = subUnit_interactions_data1.values
subUnit_interactions1 = subUnit_interactions1[desired_interactions_subUnit_1]

subUnit_interactions_data2 = pd.read_csv(subUnit_interactions_path2, header = None, index_col = None)
subUnit_interactions2 = subUnit_interactions_data2.values
subUnit_interactions2 = subUnit_interactions2[desired_interactions_subUnit_2]

node_a_list = []
node_b_list = []
link_list = []
link_weight_list = []
count = 0
count2 = 0
for i in range(1,noClasses+1):
    for j in range(noClasses+1, 2*noClasses+1):
        index_location_x = np.argwhere(np.asarray(graph_data_links_expressions_richa.iloc[:,0] == i))
        index_location_y = np.argwhere(np.asarray(graph_data_links_expressions_richa.iloc[:,1] == j))
        index_intersect = np.intersect1d(index_location_x ,index_location_y)
        interaction_gene_names_list = list(graph_data_links_expressions_richa.iloc[index_intersect,3])
        interaction_list = list(graph_data_links_expressions_richa.iloc[index_intersect,2])
        # Check if both sub units are presented:
        both_subunit_presented = []
        both_subunit_presented_pvalues = []
        for s in range(subUnit_interactions1.shape[0]):
            if subUnit_interactions1[s,:][0] in interaction_gene_names_list and subUnit_interactions2[s,:][0] in interaction_gene_names_list:
                both_subunit_presented.append(subUnit_interactions1[s,:][0] + '+' + subUnit_interactions2[s,:][0])
                k1 = interaction_gene_names_list.index(subUnit_interactions1[s,:][0])
                k2 = interaction_gene_names_list.index(subUnit_interactions2[s,:][0])
                both_subunit_presented_pvalues.append(0.5*(interaction_list[k1] + interaction_list[k2]))
            
        number_of_interactions = len(both_subunit_presented)
        for k in range(number_of_interactions):
            node_a_list.append(i)
            node_b_list.append(j)
            link_weight_list.append(both_subunit_presented_pvalues[k])
            link_list.append(both_subunit_presented[k])
            count2 = count2 + 1
        count = count + 1

graph_data_links = pd.concat([pd.DataFrame(node_a_list), 
                     pd.DataFrame(node_b_list),
                     pd.DataFrame(link_weight_list), 
                     pd.DataFrame(link_list)],
                     axis = 1)                
with open(graph_links_path,'w',newline = '') as f:
     thewriter = csv.writer(f, delimiter = ',')
     thewriter.writerow( ('from', 'to', 'weight', 'name') )
     # thewriter = csv.writer(f, dialect='exce')
     for row in range(graph_data_links.shape[0]):
         thewriter.writerow([graph_data_links.iloc[row,0], graph_data_links.iloc[row,1], graph_data_links.iloc[row,2], graph_data_links.iloc[row,3]])
print()            

# Subunit heatmap
heatmap_subunit = np.zeros([noClasses, noClasses])
count = 0
for i in range(1,noClasses+1):
    for j in range(noClasses+1, 2*noClasses+1):
        index_location_x = np.argwhere(np.asarray(graph_data_links.iloc[:,0] == i))
        index_location_y = np.argwhere(np.asarray(graph_data_links.iloc[:,1] == j))
        index_intersect = np.intersect1d(index_location_x ,index_location_y)
        heatmap_subunit[i - 1,j - (noClasses+1)] = len(index_intersect)
        count = count + 1
heatmap_subunit_df = pd.DataFrame(heatmap_subunit, index = unique_cell_type_labels.tolist(), columns = unique_cell_type_labels.tolist())
heatmap_subunit_df.to_csv(heatmap_map_path)


