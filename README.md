# DiSiR: a Software framework to identify ligand-receptor interactions at subunit level from single-cell RNA-sequencing data

<img src="https://github.com/miladrafiee/DiSiR/blob/main/Data/ReadMe_data/Method.png" width="1000">

**DiSiR** (DiSiR: Dimer Signal Receptor Analysis) is an ultrafast and easy-to-use software framework to investigate how individual cells are interacting with each other by analyzing signaling pathways of multi-subunit ligand-activated receptors from scRNA-seq data, even for the genes that are not listed in available databases. DiSiR has been built upon the framework from existing expression permutation-based methods such as CellPhoneDB (Vento-Tormo, Efremova et al. 2018, Efremova, Vento-Tormo et al. 2020) and has been specifically optimized to analyze and also visualize complex pathways involving multi-subunit receptor complexes with multi-subunit ligands for specific individual user-defined ligand-receptor input.

This repository contains the codes used for the implementation and evaluation of our method and one case study in applying it.

The primary implementation is in Python 3. To see usage example of DiSiR keep reading, and also look in the Python-module and Analysis directories. Also, the scripts to plot graph representation of resulting ligand-receptor interactions at cell type level are located in  plotting folder. In order to identify ligand-receptor interactions, 'DiSiR_main.py' script in the Python-module directory needs to be run. Below is an usage example for analyzing IL6 signaling pathway in rheumatoid arthritis (RA) Synovium scRNA-seq data. 

## DiSiR test example use

As an example, we use [AMP consortium’s Phase 1 RA data from ](https://immunogenomics.io/ampra/) published in [Zhang, Wei et al. 2019](https://www.nature.com/articles/s41590-019-0378-1). Please find the processed expression matrix, genes list and meta data [here](https://drive.google.com/drive/folders/1j6TlTqEmrCeVED2Ld3GtNbvMdzmPjMMN?usp=sharing).

<img src="https://github.com/miladrafiee/DiSiR/blob/main/Data/ReadMe_data/Ease_of_use.png" width="1000">

DiSiR only needs three input files:  1) single-cell gene expression matrix 2) cell type annotations and 3) list of desired ligand-receptor interactions at subunit level. In the 'DiSiR_main.py' file, the paths to input and output files need to be set as well as few parameters as follows.

### Path to gene expression matrix file:

Path to input gene expression matrix which has been assumed in a format of MatrixMarket (or mtx) file. 

```python

# Path to input gene expression matrix
scRNA_path = indir + '/matrix.mtx'

```

### Path to genes list file:

Path to the text file that contains gene names assoctiated with the input gene expression matrix. 

```python

# Path to names of all genes in gene expression matrix
gene_names_all_path = indir + '/genes.txt'

```

### Path to meta data file:

Path to the text file that contains gene names assoctiated with the input gene expression matrix. 

```python

# Path to metadata file
metadata_path = indir + '/categorical_coloring_data.json'

```