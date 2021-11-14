#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 12:27:25 2021

converts cellphonedb output
to nodes and links
kind of an awkward script
but it works

@author: E0463430
"""

import pandas as pd


workdirs = ['/Users/E0463430/Downloads/scripts/Milad_CATS_project/CellPhoneDB/ra_synovium_amp/cellphonedb_input/out',
            '/Users/E0463430/Downloads/scripts/Milad_CATS_project/CellPhoneDB/emma_gutman_milad/cellphonedb_input/out',
            '/Users/E0463430/Downloads/scripts/Milad_CATS_project/CellPhoneDB/fibrotic_lung_az/cellphonedb_input/out']


# FUNCTIONS TO APPLY TO COLUMNS TO GENERATE LINKS
def links_gen(df):
    reverse = 'False'
    if df.loc['receptor_a']:
        reverse='True'
        relationship = f'{reverse},{df.loc["gene_b"]} | {df.loc["partner_a"]}'
    else:
        relationship = f'{reverse},{df.loc["partner_a"]} | {df.loc["gene_b"]}'
    return relationship

# ANNOTATE INTERACTIONS
def anno_inter(interaction, anno, totes):
    interaction = interaction.split('|')
    interaction = [str(anno[interaction[0]]), str(anno[interaction[1]]+totes)]
    interaction = ','.join(interaction)
    return interaction
    

def master_script(workdir):
    # IMPORTATION OF DATA
    sig_matrix = pd.read_csv(f'{workdir}/significant_means.txt', sep='\t')
    sig_matrix = sig_matrix.T
    
    # GENERATION OF NODES -- AT LEAST ENOUGH TO PLOT
    inters = [i.split('|') for i in list(sig_matrix.index) if '|' in i]
    flat_inters = [item for sublist in inters for item in sublist]
    cell_states = sorted(set(flat_inters))
    totes = len(cell_states)
    w_cell_states = cell_states*2
    
    anno = {}
    with open(f'{workdir}/Nodes.csv', 'w') as node_handle:
        node_handle.write('id,name\n')
        for nstate, state in zip(range(1, len(w_cell_states)+1),
                                 w_cell_states):
            if state not in anno:
                anno[state] = nstate
            node_handle.write(f'{str(nstate)},{state}\n')
    
              
    # GENERATION OF LINKS
    with open(f'{workdir}/Links.csv', 'w') as links_handle:
        links_handle.write('from,to,weight,name\n')
        relats = sig_matrix.apply(links_gen)
        sig_matrix.columns = relats
        inter_matrix = sig_matrix.loc[[i for i in sig_matrix.index if '|' in i]]
        inter_matrix.fillna(0, inplace=True)
        for col in inter_matrix.columns:
            test = inter_matrix[col]
            if col.split(',')[0] == 'True':
                test.index = ['|'.join(i.split('|')[::-1]) for i in list(test.index)]
            test.index = [anno_inter(i, anno, totes) for i in list(test.index)]
            test = test.to_dict()
            for key in sorted(test):
                if test[key] > 0:
                    output = f'{key},{round(test[key], 3)},{col.split(",")[1]}'
                    links_handle.write(f'{output}\n')


###############################

if __name__ == '__main__':
    for workdir in workdirs:
        master_script(workdir)
        
###############################
