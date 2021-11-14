#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 21:15:42 2020

@author: Milad R. Vahid
"""

import pandas as pd
import numpy as np
import scipy.stats
import numpy.matlib
from scipy.io import mmread

from scipy.stats.mstats import rankdata

from math import ceil
from math import log2

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as pch
from matplotlib import rcParams
from statsmodels.stats.multitest import multipletests

mpl.use('Agg')
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
#rcParams.update({'font.size': 15})

#import seaborn as sns

#################################################################################
############################## Functions ########################################
#################################################################################

def rm_outlier_mean(df, percentile=97.5):
    '''to exclude outliers. currently not in use'''
    df = df[df < np.percentile(df, percentile)]
    if len(df)==0:
        return 0
    else:
        return np.mean(df)

    
def transform_to_rank(df):
    '''transfer expression values to rank; decided against it'''
    mat = df.values
    mat_rank = rankdata(mat, axis = 1)
    mat_rank[np.where(mat  == 0)] = 0
    return pd.DataFrame(mat_rank)
 
   
def parse_interactions(subunit_interactions_path):
    LR_info = {}
    subunit_info = {}
    subunit_interactions = pd.read_csv(subunit_interactions_path, header = None, index_col = None)
    for row in range(len(subunit_interactions)):
        gene_names = []
        inter_names = list(subunit_interactions.iloc[row, : ].dropna().values)
        subunit_info[row] = inter_names
        for inter in inter_names:
            inter_sp = inter.split(' | ')
            gene_1 = inter_sp[0]
            gene_2 = inter_sp[1]
            if gene_1 not in gene_names:
                gene_names.append(gene_1)
            if gene_2 not in gene_names:
                gene_names.append(gene_2)
        LR_info[row] = gene_names
        
    return LR_info, subunit_info


def filter_out_unannotated_cells(scRNA_path,
                                 metadata_path,
                                 select_track,
                                 exclude_tracks,
                                 subset):
    
    metadata = pd.read_json(metadata_path)
    cell_type_labels = metadata.loc['label_list', select_track]
    cell_type_keep = [i[0] for i in enumerate(cell_type_labels) if i[1] not in exclude_tracks]
    if subset: 
        subset_track = list(subset.keys())[0]
        cell_type_track = metadata.loc['label_list',  subset_track]
        subset_cats  = subset[subset_track]
        
        not_in_system = set(subset_cats) - set(cell_type_track)
        if len(not_in_system) > 0:
            print(f'WARNING! {not_in_system} ARE NOT IN {subset_track}')
            
        cell_type_keep = [i[0] for i, j in zip(enumerate(cell_type_labels), cell_type_track)
                          if i[1] not in exclude_tracks and j in subset_cats]
        
    cell_type_labels = [i[1] for i in enumerate(cell_type_labels) if i[0] in cell_type_keep]
    scRNA_data = mmread(scRNA_path)
    scRNA_array = scRNA_data.toarray()
    scRNA_array = scRNA_array[cell_type_keep, ]
    scRNA_array = np.transpose(scRNA_array)
    
    return scRNA_array, cell_type_labels
    

def gene_expressions_with_selected_genes(scRNA_array,
                                         gene_names_all_path,
                                         subunit_interactions_path):

    gene_names_all = list(pd.read_csv(gene_names_all_path, header = None, index_col = None)[0])
    LR_info, subunit_info = parse_interactions(subunit_interactions_path)
    
    gene_names = list(LR_info.values())
    gene_names = sorted(set([item for subl in gene_names for item in subl]))

    gene_index = [gene_names_all.index(i) for i in gene_names]
    
    expressions_scRNA = scRNA_array[gene_index, : ]

    # Data normalization (here, log normalization)
    expressions_scRNA = np.log2(expressions_scRNA + 1)
        
    return expressions_scRNA, gene_names, LR_info, subunit_info


def calculate_celltype_average_expressions(expressions_scRNA,
                                           gene_names,
                                           cell_type_labels):    
    # Read meta data (cell type labels)
    unique_cell_type_labels = np.unique(cell_type_labels)
    noClasses = len(unique_cell_type_labels)
    cell_type_labels_df = pd.DataFrame(cell_type_labels)
    
    gene_number = len(gene_names)
    
    # Calculate cell type average expression matrix
    average_expression = np.zeros([gene_number,noClasses])
    number_expression = np.zeros([gene_number,noClasses])
    fraction_expression = np.zeros([gene_number,noClasses])
    
    cell_type_specific_numbers = np.zeros([gene_number,noClasses])
    cell_type_numbers = [0]*noClasses
    for k in range(noClasses):
        cell_type_index = cell_type_labels_df.isin([unique_cell_type_labels[k]])
        seriesObj = cell_type_index.any(axis = 1) 
        columnNames = list(seriesObj[seriesObj == True].index) 
        cell_type_numbers[k] = len(columnNames)
        for i in range(gene_number):
            average_expression[i,k] = np.mean(expressions_scRNA[i,columnNames])
            number_expression[i,k] = len(expressions_scRNA[i,columnNames][expressions_scRNA[i,columnNames] > 0])
            fraction_expression[i,k]  = number_expression[i,k]/len(columnNames)           
            cell_type_specific_numbers[i,k] = np.count_nonzero(expressions_scRNA[i,columnNames])/cell_type_numbers[k]
    
    average_expression_df = pd.DataFrame(average_expression, index = gene_names, columns = unique_cell_type_labels)
    number_expression_df = pd.DataFrame(number_expression, index = gene_names, columns = unique_cell_type_labels)
    fraction_expression_df = pd.DataFrame(fraction_expression, index = gene_names, columns = unique_cell_type_labels)   
    
    cell_type_specific_expressions = average_expression/np.max(average_expression)
    return average_expression_df, number_expression_df, fraction_expression_df, unique_cell_type_labels, cell_type_numbers, cell_type_specific_numbers, cell_type_specific_expressions


def calculate_LR_interactions_matrix_prefiltered(average_expression_df,
                            gene_names,
                            unique_cell_type_labels,
                            cell_type_specific_numbers,
                            cell_type_specific_expressions,
                            threshold_numbers,
                            threshold_expressions):
    
    gene_number = len(gene_names)
    noClasses = len(unique_cell_type_labels)
        
    # Calculate p-value of belonging Ligand (or receptor) to cell types
    p_values =  1000 * np.ones([gene_number,noClasses]) 
    for i in range(gene_number):
        for j in range(noClasses):
            boolean_index_selected =  average_expression_df.columns != average_expression_df.columns[j]
            t, p = scipy.stats.ttest_1samp(np.array(average_expression_df.iloc[i, boolean_index_selected]), average_expression_df.iloc[i,j]) 
            if t < 0:
               p_values[i,j] = p
    
    interactions = np.zeros(gene_number ** 2 * noClasses ** 2) 
    interactions_class_names = ['a'] * (gene_number ** 2 * noClasses ** 2)
    interactions_gene_names = ['a'] * (gene_number ** 2 * noClasses ** 2)
    counter = 0
    for i_g in range(gene_number):
        for j_g in range(gene_number):
            for i_c in range(noClasses):
                for j_c in range(noClasses):
                    interactions_class_names[counter] = unique_cell_type_labels[i_c] + ' '+ '|' + ' ' + unique_cell_type_labels[j_c]
                    interactions_gene_names[counter] = gene_names[i_g] + ' '+ '|' + ' '  + gene_names[j_g]
                    if cell_type_specific_numbers[i_g,i_c] > threshold_numbers and cell_type_specific_expressions[i_g,i_c] >  threshold_expressions and cell_type_specific_numbers[j_g,j_c] > threshold_numbers and cell_type_specific_expressions[j_g,j_c] > threshold_expressions: 
                      # interactions[counter] = np.mean([average_expression_df.iloc[i_g,i_c], average_expression_df.iloc[j_g,j_c]])
                      interactions[counter] = average_expression_df.iloc[i_g,i_c] * average_expression_df.iloc[j_g,j_c]
                    counter = counter + 1

    return interactions, interactions_gene_names, interactions_class_names


def calculate_LR_interactions_pvalues(interactions,
                                   interactions_gene_names,
                                   interactions_class_names,
                                   expressions_scRNA,
                                   cell_type_numbers,
                                   gene_names,
                                   unique_cell_type_labels,
                                   iteration):
                                                                    
    gene_number = len(gene_names)
    noClasses = len(unique_cell_type_labels)
    
    # Shuffle cell labels between all cells for "iteration" times 
    interactions_shuffle = np.zeros([gene_number ** 2 * noClasses ** 2, iteration])  
    average_expression_shuffle = np.zeros([gene_number,noClasses])   
    for no_iter in range(iteration):
        for k in range(noClasses):
            columnNames = np.random.choice(range(expressions_scRNA.shape[1]), cell_type_numbers[k]).tolist()
            for i in range(gene_number):
                average_expression_shuffle[i,k] = np.mean(expressions_scRNA[i,columnNames])
    
        counter = 0
        for i_g in range(gene_number):
            for j_g in range(gene_number):
                for i_c in range(noClasses):
                    for j_c in range(noClasses):
                        # interactions_shuffle[counter,no_iter] = np.mean([average_expression_shuffle[i_g,i_c], average_expression_shuffle[j_g,j_c]])
                        interactions_shuffle[counter,no_iter] = average_expression_shuffle[i_g,i_c] * average_expression_shuffle[j_g,j_c]
                        counter = counter + 1
          
    p_values_matrix = 1000* np.ones([gene_number ** 2, noClasses ** 2])
    expression_celltype_matrix = np.zeros([gene_number ** 2, noClasses ** 2])
    p_values =  1000*np.ones([gene_number ** 2 * noClasses ** 2])
    counter = 0
    counter_1 = 0
    for i_g in range(gene_number):
        for j_g in range(gene_number):
            counter_2 = 0 
            for i_c in range(noClasses):
                for j_c in range(noClasses):
                    t, p = scipy.stats.ttest_1samp(interactions_shuffle[counter,:], interactions[counter]) 
                    if t < 0:
                       p_values_matrix[counter_1,counter_2] = p
                       expression_celltype_matrix[counter_1,counter_2] = interactions[counter]
                       p_values[counter] = p
                    counter = counter + 1
                    counter_2 = counter_2 + 1
            counter_1 = counter_1 + 1
    
    output_table_results = pd.concat([pd.DataFrame(interactions_gene_names), 
                         pd.DataFrame(interactions_class_names),
                         pd.DataFrame(p_values), 
                         pd.DataFrame(interactions)],
                         axis = 1)    
    p_values_gene_names = ['a'] * (gene_number ** 2)
    counter = 0
    for i_g in range(gene_number):
        for j_g in range(gene_number):
            p_values_gene_names[counter] = gene_names[i_g] + ' '+ '|' + ' '  + gene_names[j_g]
            counter = counter + 1
    
    p_values_celltype_names = ['a'] * (noClasses ** 2)
    counter = 0
    for i_c in range(noClasses):
        for j_c in range(noClasses):
            p_values_celltype_names[counter] = unique_cell_type_labels[i_c] + ' '+ '|' + ' '  + unique_cell_type_labels[j_c]
            counter = counter + 1
    return pd.DataFrame(p_values_matrix), pd.DataFrame(p_values_gene_names), pd.DataFrame(p_values_celltype_names), pd.DataFrame(expression_celltype_matrix), output_table_results    
       

def build_interactions_graph_subunit(p_values_matrix,
                             p_values_gene_names,
                             unique_cell_type_labels,
                             expression_celltype_matrix,
                             selected_interactions,
                             threshold):    

    noClasses = len(unique_cell_type_labels)

    matched_gene_names = p_values_gene_names.iloc[:,0].isin(selected_interactions)
    
    index_matched_gene_names = list(matched_gene_names[matched_gene_names == True].index)
    selected_matched_gene_names = [i for i in list(p_values_gene_names[0]) if i in selected_interactions]

    p_values_matrix = p_values_matrix.iloc[index_matched_gene_names, :]
    p_values_matrix = p_values_matrix.values
    p_values_matrix = np.array([multipletests(i, method='fdr_bh')[1] for i in p_values_matrix])

    expression_celltype_matrix = expression_celltype_matrix.iloc[index_matched_gene_names, :]
    expression_celltype_matrix = expression_celltype_matrix.values
    node_a_list = []
    node_b_list = []
    link_list = []
    link_weight_list = []
    count = 0
    count2 = 0
    for i in range(noClasses):
        for j in range(noClasses):
            significant_pvalues = p_values_matrix[:,count] <= threshold
            interaction_list = expression_celltype_matrix[significant_pvalues == True, count]
            interaction_gene_names_list = [i for i,j in zip(selected_matched_gene_names, list(significant_pvalues)) if j == True]
            # Check if both sub units are presented:
            both_subunit_presented = []
            both_subunit_presented_pvalues = []
            #print(subUnit_interactions_data, 'subUnit_interactions_data')
            if len(set(selected_interactions) & set(interaction_gene_names_list)) == len(selected_interactions):
                both_subunit_presented.append('+'.join(selected_interactions))
                both_subunit_presented_pvalues.append(np.mean([i for i,j in zip(interaction_list, 
                                                                                interaction_gene_names_list) if j in selected_interactions]))
            number_of_interactions = len(both_subunit_presented)
            for k in range(number_of_interactions):
                node_a_list.append(i + 1)
                node_b_list.append(noClasses + j + 1)
                link_weight_list.append(both_subunit_presented_pvalues[k])
                link_list.append(both_subunit_presented[k])
                count2 = count2 + 1
            count = count + 1
    
    graph_data_links = pd.concat([pd.DataFrame(node_a_list), 
                         pd.DataFrame(node_b_list),
                         pd.DataFrame(link_weight_list), 
                         pd.DataFrame(link_list)],
                         axis = 1)
    
    node_data_id = list(range(1,2*noClasses+1))
    node_data_type = [2]*len(unique_cell_type_labels) + [1]*len(unique_cell_type_labels) 
    node_data_name =  unique_cell_type_labels.tolist() + unique_cell_type_labels.tolist() 
    x = list(range(1,noClasses * 20 + 1, 20)) + list(range(1,noClasses * 20 + 1, 20))
    y = node_data_type    
    graph_data_nodes = pd.concat([pd.DataFrame(node_data_id),
                         pd.DataFrame(node_data_name), 
                         pd.DataFrame(node_data_type),
                         pd.DataFrame(y),
                         pd.DataFrame(x)],
                         axis = 1)                         
   
    return graph_data_nodes, graph_data_links 


def calculate_celltype_interactions_heatmap(graph_data_links,
                                            unique_cell_type_labels,
                                            average_expression_df,
                                            selected_interactions):
    heatmap_dict = {}
        
    noClasses = len(unique_cell_type_labels)
    
    receptors = [i.split(' | ')[1] for i in selected_interactions]
    ligands = [i.split(' | ')[0] for i in selected_interactions]

    heatmap_dict = {}
    heatmap_dict['heatmap'] = np.zeros([noClasses, noClasses])
    
    if not graph_data_links.empty:
        for i in range(1,noClasses+1):
            for j in range(noClasses+1, 2*noClasses+1):
                index_location_x = np.where(graph_data_links.iloc[:,0] == i)
                index_location_y = np.where(graph_data_links.iloc[:,1] == j)
                index_intersect = np.intersect1d(index_location_x ,index_location_y)
                if len(index_intersect) != 0:
                    heatmap_dict['heatmap'][i - 1,j - (noClasses+1)] = graph_data_links.iloc[index_intersect, 2]
            
    heatmap_dict['heatmap_all_interactions'] = np.zeros([noClasses, noClasses])
    for i, label_i in enumerate(unique_cell_type_labels):
        for j, label_j in enumerate(unique_cell_type_labels):
            expression_sum = np.mean([average_expression_df.loc[k, label_i]*average_expression_df.loc[m, label_j] for k,m in zip(ligands,receptors)])
            heatmap_dict['heatmap_all_interactions'][i,j] = expression_sum
        
    heatmap_dict['heatmap'] = pd.DataFrame(heatmap_dict['heatmap'], index = unique_cell_type_labels, columns = unique_cell_type_labels)
    heatmap_dict['heatmap_all_interactions'] = pd.DataFrame(heatmap_dict['heatmap_all_interactions'], index = unique_cell_type_labels, columns = unique_cell_type_labels)

    return heatmap_dict


def plot_cell_dist(average_expression_df, 
                   number_expression_df,
                   fraction_expression_df,
                   outdir):
        
    Ny = len(average_expression_df)
    Nx = len(average_expression_df.columns)
    
    fig = plt.figure(figsize=(Nx*1.2, Ny), dpi=1500)
    ax = plt.axes()
    
    ylabels = list(average_expression_df.index)
    xlabels = list(average_expression_df.columns)
    
    max_val = ceil(max(average_expression_df.max()))
    
    norm = mpl.colors.Normalize( vmin=0, vmax=max_val )
    mapper = cm.ScalarMappable(norm=norm, cmap='Reds')
    
    circles = []
    
    for i in range(Nx):
        for j in range(Ny):
            rad = log2(fraction_expression_df.iloc[j, i]*100+1)/13.28
            col = mapper.to_rgba(average_expression_df.iloc[j, i])
            patch = pch.Circle((i+1, j+0.5), rad, fc=col, ec='black')
            text = f'{int(number_expression_df.iloc[j, i])} ({round(fraction_expression_df.iloc[j, i]*100, 1)}%)'
            ax.text(x=i+0.55, y=j+0.9, s=text, size=7.5)
            ax.add_patch(patch)
            circles.append(patch)
    
    ax.set(xticks=np.arange(1, Nx+2), yticks=np.arange(0.5, Ny+0.5), yticklabels=ylabels)
    ax.set_ylim(0, Ny+0.25)
    ax.set_xticklabels(xlabels+[''], rotation=80)
    
    plt.colorbar(mapper, ax=ax, label='log(expression)', pad=0.025)
    
    fig.savefig(f'{outdir}/expression_info.pdf', bbox_inches='tight')
    plt.close()

    
def plot_heatmaps(heatmaps, cellstates, outdir):
    
    for ht in heatmaps:
        
        heatmatrix = heatmaps[ht]
        maxval = max(heatmatrix.max())
    
        fig = plt.figure(figsize=(7.5, 5.5))
        ax = plt.axes()
        pcm = ax.imshow(heatmatrix, cmap = 'YlGnBu', vmin=0, vmax=maxval)
        
        ax.set_xticks(range(len(cellstates)))
        ax.set_yticks(range(len(cellstates)))
        ax.set_xticklabels(cellstates, rotation=90)
        ax.set_yticklabels(cellstates)
        ax.set_xlabel('Receptor cell state')
        ax.set_ylabel('Ligand cell state')
        fig.colorbar(pcm, ax=ax, label='log expression')
        
        for val in np.arange(0.5, len(cellstates)):
            ax.axhline(val, color='0.8', lw=1)
            ax.axvline(val, color='0.8', lw=1)
    
        fig.savefig(f'{outdir}/{ht}.pdf', bbox_inches='tight')
        plt.close()

#################################################################################
#################################################################################
#################################################################################
