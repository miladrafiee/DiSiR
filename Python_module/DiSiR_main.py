#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 16:38:34 2021

@author: E0463430

Pipeline for DiSir
Code and method development: Milad Vahid
Pipeline and minor edits: Andre Kurlovs

"""

import argparse
import os
import pandas as pd
import sys
import subprocess
import csv
import numpy as np

###############################################################################
############################### MASTER INPUTS #################################
###############################################################################


PARSER = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

# INPUT DIRECTORY SHOULD CONTAIN SPRING OUTPUT:
    # categorical_coloring_data.json
    # genes.txt
    # matrix.mtx

# THE INTERACTIONS FILE
# format: ligand | receptor
# column1: first subunit relationship
# column2: second subunit relationship
#EXAMPLE:
#IL1B | IL1RAP,IL1B | IL1R1  
#IL1RN | IL1RAP,IL1RN | IL1R1
        

PARSER.add_argument("-i", "--indir", required=True,
                    help="input directory with spring data:"
                          "matrix.mtx, categorical_coloring_data.json, genes.txt")

PARSER.add_argument("-o", "--outdir", required=True,
                    help="output directory")

PARSER.add_argument("-n", "--interactions", required=True,
                    help="comma-separated interaction file"
                         "one interaction per line"
                         "example: IL1A | IL1R1,IL1A | IL1RAP")

PARSER.add_argument("-s", "--select", required=False,
                    help="Select metadata track for analysis",
                    default="CellStates")

PARSER.add_argument("-x", "--exclude", required=False,
                    help="Select groups in the metadata track for exclusion",
                    default='Unclassified,Other')

PARSER.add_argument("-u1", "--subset_track", required=False,
                    help="Select track e.g Disease")

PARSER.add_argument("-u2", "--subset_categories", required=False,
                    help="Select categories e.g. SLE")

PARSER.add_argument("-t1", "--threshold_numbers", required=False,
                    help="Threshold on number of cells expressed each ligand or receptor per cell type",
                    default=0)
                    
PARSER.add_argument("-t2", "--threshold_expressions", required=False,
                    help="Threshold on scaled (max-normalized) average expression of each ligand or receptor within a cell type",
                    default=0)                    

PARSER.add_argument("-p", "--threshold", required=False,
                    help="Threshold on p-value for filtering non-significant LR interactions",
                    default=0.05)

PARSER.add_argument("-r", "--iterations", required=False,
                    help="Number of iterations for permutating data in statistical test",
                    default=100)


ARGIES = PARSER.parse_args()
indir = ARGIES.indir
outdir = ARGIES.outdir
subunit_interactions_path = ARGIES.interactions

select_track = ARGIES.select

if ARGIES.subset_track and ARGIES.subset_categories:
    subset = {}
    subset[ARGIES.subset_track] = ARGIES.subset_categories.split(',')
else:
    subset = False


exclude_tracks = ARGIES.exclude.split(',')
threshold_numbers = float(ARGIES.threshold_numbers)
threshold_expressions = float(ARGIES.threshold_expressions)
threshold = float(ARGIES.threshold)
iteration = int(ARGIES.iterations)

append = '/'.join(os.path.abspath(__file__).split('/')[:-1]) + '/'
sys.path.append(append)

###############################################################################
############################### MORE INPUTS ##################################
###############################################################################

# Path to input gene expression matrix
scRNA_path = indir + '/matrix.mtx'
# scRNA_path = '/cloud-data/cloud-pipeline-milad-storage/IL1RAP_data/matrix.mtx'

# Path to names of all genes in gene expression matrix
gene_names_all_path = indir + '/genes.txt'
# gene_names_all_path = '/cloud-data/cloud-pipeline-milad-storage/IL1RAP_data/genes.txt'

# Path to metadata
metadata_path = indir + '/categorical_coloring_data.json'
# metadata_path = '/cloud-data/cloud-pipeline-milad-storage/IL1RAP_data/categorical_coloring_data.json'

###################################################################################################
###################################### RUN MOD CELLPHONEDB ########################################
###################################################################################################

import DiSiR_functions as ds

np.random.seed(0)

scRNA_array, cell_type_labels = ds.filter_out_unannotated_cells(scRNA_path,
                             metadata_path, select_track, 
                             exclude_tracks, subset)


expressions_scRNA_allLR, gene_names_allLR, LR_info, subunit_info = ds.gene_expressions_with_selected_genes(scRNA_array,
                                            gene_names_all_path,
                                            subunit_interactions_path)

for LR_set in LR_info:
    
    gene_names = LR_info[LR_set]
    gene_index = [gene_names_allLR.index(i) for i in gene_names]
    expressions_scRNA = expressions_scRNA_allLR[gene_index, : ]
    
    selected_interactions = subunit_info[LR_set]
    
    average_expression_df, number_expression_df, fraction_expression_df, unique_cell_type_labels, cell_type_numbers, cell_type_specific_numbers, cell_type_specific_expressions = ds.calculate_celltype_average_expressions(expressions_scRNA,
                                               gene_names,
                                               cell_type_labels)
    
    ######################### Transfer to rank space #############################
    # average_expression_df = ds.transform_to_rank(average_expression_df)
    
    ##############################################################################
    interactions_prefiltered, interactions_gene_names, interactions_class_names = ds.calculate_LR_interactions_matrix_prefiltered(average_expression_df,
                                gene_names,
                                unique_cell_type_labels,
                                cell_type_specific_numbers,
                                cell_type_specific_expressions,
                                threshold_numbers,
                                threshold_expressions)
    
    
    p_values_matrix, p_values_gene_names, p_values_celltype_names, expression_celltype_matrix, output_table_results = ds.calculate_LR_interactions_pvalues(interactions_prefiltered,
                                       interactions_gene_names,
                                       interactions_class_names,
                                       expressions_scRNA,
                                       cell_type_numbers,
                                       gene_names,
                                       unique_cell_type_labels,
                                       iteration)
    

    graph_data_nodes, graph_data_links = ds.build_interactions_graph_subunit(p_values_matrix,
                                 p_values_gene_names,
                                 unique_cell_type_labels,
                                 expression_celltype_matrix,
                                 selected_interactions,
                                 threshold)
    
    heatmap_dict = ds.calculate_celltype_interactions_heatmap(graph_data_links,
                                                unique_cell_type_labels,
                                                average_expression_df,
                                                selected_interactions)
    
    ###############################################################################
    #################### Define and save outputs and plots #######################
    ###############################################################################
    
    relation = '_and_'.join(selected_interactions).replace(' | ', '_').replace(' ', '_')
    
    if not os.path.isdir(f'{outdir}'):
        subprocess.call(f'mkdir {outdir}', shell=True)
        
    if not os.path.isdir(f'{outdir}/{relation}'):
        subprocess.call(f'mkdir {outdir}/{relation}', shell=True)
        
    if not os.path.isdir(f'{outdir}/{relation}/heatmaps'):
        subprocess.call(f'mkdir {outdir}/{relation}/heatmaps', shell=True)
        
    heatmap_dict['heatmap'].to_csv(f'{outdir}/{relation}/heatmaps/Heatmap.csv')
    heatmap_dict['heatmap_all_interactions'].to_csv(f'{outdir}/{relation}/heatmaps/Heatmap_all_interactions.csv')
    
    ds.plot_heatmaps(heatmap_dict, unique_cell_type_labels, outdir=f'{outdir}/{relation}/heatmaps/')
    
    enzymes = [i.split('_') for i in relation.split('_and_')]
    enzymes = sorted(set([item for sublist in enzymes for item in sublist]))
    ds.plot_cell_dist(average_expression_df.loc[enzymes],
                      number_expression_df.loc[enzymes],
                      fraction_expression_df.loc[enzymes],
                      f'{outdir}/{relation}')
    
    
    if not os.path.isdir(f'{outdir}/{relation}/graph_data'):
        subprocess.call(f'mkdir {outdir}/{relation}/graph_data', shell=True)
        
    if not os.path.isdir(f'{outdir}/{relation}/network_plots'):
        subprocess.call(f'mkdir {outdir}/{relation}/network_plots', shell=True)
     
    with open(f'{outdir}/{relation}/graph_data/Links.csv','w',newline = '') as f:
         thewriter = csv.writer(f, delimiter = ',')
         thewriter.writerow( ('from', 'to', 'weight', 'name') )
         # thewriter = csv.writer(f, dialect='exce')
         for row in range(graph_data_links.shape[0]):
             thewriter.writerow([graph_data_links.iloc[row,0], graph_data_links.iloc[row,1], graph_data_links.iloc[row,2], graph_data_links.iloc[row,3]])
    print()
                 
    with open(f'{outdir}/{relation}/graph_data/Nodes.csv','w',newline = '') as f:
         thewriter = csv.writer(f, delimiter = ',')
         thewriter.writerow( ('id', 'name', 'type', 'y', 'x') )
         # thewriter = csv.writer(f, dialect='exce')
         for row in range(graph_data_nodes.shape[0]):
             thewriter.writerow([graph_data_nodes.iloc[row,0], graph_data_nodes.iloc[row,1], graph_data_nodes.iloc[row,2], graph_data_nodes.iloc[row,3], graph_data_nodes.iloc[row,4]])
    print()
    
    
    ###############################################################################
    ########################## Additional Plotting in R ###########################
    ###############################################################################
    
    
    subprocess.call(f'Rscript {append}DiSiR_subunit.R {outdir}/{relation}/graph_data/ {outdir}/{relation}/network_plots/', shell=True)
    

###############################################################################
###############################################################################
###############################################################################

