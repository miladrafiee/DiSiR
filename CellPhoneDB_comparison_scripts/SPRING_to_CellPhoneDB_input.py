#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 18:03:06 2021

this script preps our files for cellphonedb

@author: E0463430
"""

import json
import scipy.io
#import numpy as np
import os
import pandas as pd
import subprocess
import sys

############################################################
####################### INPUTS #############################
############################################################ 

indir = sys.argv[1] # spring files must contain genes, json annotation, matrix, and barcodes in folder
sel_genes = sys.argv[2].split(',') # which genes to select
outdir = sys.argv[3] # where to put the output
       
# synovium = IL6,IL6R,IL6ST,TNF,TNFRSF1A,TNFRSF1B
# ad = IL1B,IL1RL1,IL1RN,IL1R1,IL1RAP,IL33
# arizona = IL1A,IL1B,IL1RL1,IL1RL2,IL1RN,IL1R1,IL1RAP,IL33,IL36B

############################################################
################### PROCESS BARCODES #######################
############################################################ 

in_barcodes = indir + '/barcodes/'
in_matrix = indir + '/matrix.mtx'
in_genes = indir + '/genes.txt'
in_json = indir + '/categorical_coloring_data.json'
outfile_mat = f'{outdir}/counts.txt'
outfile_meta = f'{outdir}/meta.txt'

if not os.path.isdir(f'{outdir}'):
    subprocess.call(f'mkdir {outdir}', shell=True)

############################################################
################### PROCESS BARCODES #######################
############################################################

samples_order = []


with open(in_json) as in_json_handle:
    metadata = json.load(in_json_handle)
    for sample_name in metadata['Sample']['label_list']:
        if sample_name not in samples_order:
            samples_order.append(sample_name)


# if barcodes present, get them
# if not, label the cells sequentially
barcodes_all = []
if os.path.isdir(f'{outdir}/barcodes'):
    for sample in samples_order:
        with open(f'{in_barcodes}/{sample}') as bar_handle:
            barcodes = bar_handle.readlines()
            barcodes = [i[2:-2] for i in barcodes]
            barcodes_all += barcodes
else:
    for ncell in range(1, len(metadata['Sample']['label_list'])+1):
        barcodes_all.append(f'Cell_{ncell}')


cellstates_keep = []
barcodes_keep = []

if 'CellStates' in metadata:
    anno_track = 'CellStates'
elif 'CellStatesID' in metadata:
    anno_track = 'CellStatesID'
else:
    sys.exit('Cell state annotation tracks not found')

for anno, bar in zip(metadata[anno_track]['label_list'],
                     barcodes_all):
    if anno not in ['Other', 'Unclassified']:
        cellstates_keep.append(anno)
        barcodes_keep.append(bar)
        
############################################################
##################### PROCESS GENES ########################
############################################################

with open(in_genes) as genes_handle:
    genes_all = genes_handle.readlines()
    genes_all = [gene.rstrip() for gene in genes_all]
    
############################################################
##################### PROCESS MATRIX #######################
############################################################

mat = scipy.io.mmread(in_matrix)
dense = mat.todense().transpose()

densedf = pd.DataFrame(dense)
densedf.index = genes_all
densedf.columns = barcodes_all

select_mat = densedf.loc[sel_genes]
select_mat = select_mat[barcodes_keep]

# remove duplicate barcodes as they will cause issues in CellPhoneDB
no_duplicates = [i for i in barcodes_keep if barcodes_keep.count(i)==1]
select_mat = select_mat[no_duplicates]

select_mat.to_csv(outfile_mat, sep='\t', index_label='Gene')

############################################################
################### OUTPUT METADATA ########################
############################################################

meta_dict = {}

# remove duplicates here as well
for bar, anno in zip(barcodes_keep, cellstates_keep):
    if bar in no_duplicates:
        meta_dict[bar] = anno

meta_ser = pd.Series(meta_dict)
meta_ser.to_csv(outfile_meta, sep='\t', header=['cell_type'], index_label='Cell')

############################################################
############################################################
############################################################
