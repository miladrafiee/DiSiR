# DiSiR: a Software framework to identify ligand-receptor interactions at subunit level from single-cell RNA-sequencing data

<img src="https://github.com/miladrafiee/DiSiR/blob/main/Data/ReadMe_data/Method.png" width="500">

**DiSiR** an ultrafast and easy-to-use software framework to investigate how individual cells are interacting with each other by analyzing signaling pathways of multi-subunit ligand-activated receptors from scRNA-seq data, even for the genes that are not listed in available databases.

This repository contains the codes used for the implementation and evaluation of our method and one case study in applying it.

The key innovations of our method are:

- In contrast to conventional methods, CytoSPACE dissects spatial organizations of cells in a given tissue at single cell level.

- Since our method maps single cells from scRNA-sequencing data, in which larger numbers of genes are sequenced per each cell compared to availavle spatial transcriptomics technology, our method imporves the gene coverage of a recontructed tissue significantly.

- Our method requires no prior knowledge about the cell types and cell states.

The primary implementation is in Python 3. To see usage example of CytoSPACE keep reading, and also look in the Python-module and Analysis directories. In order to identify cell-to-spot assignment, 'CytoSPACE.py' code in the Analysis directory needs to be run. Below is an usage example for Human Colorectal Cancer (CRC) 10X Visium spatial transcriptomics data. 

## CytoSPACE test example use

As an example, we use [Human Colorectal Cancer 10X Visium Spatial Transcriptomics data](https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_ColorectalCancer/Parent_Visium_Human_ColorectalCancer_web_summary.html). Also, the corresponding scRNA-seq data that we use here is from the Cancer Genome Atlas (TCGA) published in [Lee et al 2020](https://www.nature.com/articles/s41588-020-0636-z). Please find the full data [here](https://drive.google.com/drive/folders/1j6TlTqEmrCeVED2Ld3GtNbvMdzmPjMMN?usp=sharing).

In the 'CytoSPACE.py' file, the paths to input and output files need to be set as well as few parameters as follows.

### Paths to ST files:

```python

# Path to spatial transcriptomics data (gene expression matrix) 
ST_path = '.../TCGA_ST_crc_10x_tpm.fullIDs.remapped.csv'

# Path to spatial transcriptomics data (x-y coordinates) 
coordinates_path = '.../Coordinates.csv'

```
The formats of these two files are as below:

1. ST data:

