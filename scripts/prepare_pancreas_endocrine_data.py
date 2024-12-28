import numpy as np
import scanpy.api as sc
# import matplotlib.pyplot as pl
# import pandas as pd
# from scanpy.tools import rna_velocity
# from anndata import AnnData
# import seaborn as sns
# from scipy.sparse import csr_matrix
# import networkx as nx
# import xlsxwriter
# from matplotlib import rcParams
# import seaborn as sns
# import scipy as sci
# import gseapy as gp
# sc.settings.verbosity = 3
# sc.logging.print_versions()

# Read cellranger files for all four samples
filename = './E12_5_counts/mm10/matrix.mtx'
filename_genes = './E12_5_counts/mm10/genes.tsv'
filename_barcodes = './E12_5_counts/mm10/barcodes.tsv'

e125 = sc.read(filename).transpose()
e125.var_names = np.genfromtxt(filename_genes, dtype=str)[:, 1]
e125.obs_names = np.genfromtxt(filename_barcodes, dtype=str)

filename = './E13_5_counts/mm10/matrix.mtx'
filename_genes = './E13_5_counts/mm10/genes.tsv'
filename_barcodes = './E13_5_counts/mm10/barcodes.tsv'

e135 = sc.read(filename).transpose()
e135.var_names = np.genfromtxt(filename_genes, dtype=str)[:, 1]
e135.obs_names = np.genfromtxt(filename_barcodes, dtype=str)

filename = './E14_5_counts/mm10/matrix.mtx'
filename_genes = './E14_5_counts/mm10/genes.tsv'
filename_barcodes = './E14_5_counts/mm10/barcodes.tsv'

e145 = sc.read(filename).transpose()
e145.var_names = np.genfromtxt(filename_genes, dtype=str)[:, 1]
e145.obs_names = np.genfromtxt(filename_barcodes, dtype=str)

filename = './E15_5_counts/mm10/matrix.mtx'
filename_genes = './E15_5_counts/mm10/genes.tsv'
filename_barcodes = './E15_5_counts/mm10/barcodes.tsv'

e155 = sc.read(filename).transpose()
e155.var_names = np.genfromtxt(filename_genes, dtype=str)[:, 1]
e155.obs_names = np.genfromtxt(filename_barcodes, dtype=str)

# Add dev. timepoint label for each sample
e125.obs['day'] = '12.5'
e135.obs['day'] = '13.5'
e145.obs['day'] = '14.5'
e155.obs['day'] = '15.5'

# Create Concatenated anndata object for all timepoints
alldays = e125.concatenate(e135, e145, e155)

# Deleting individual day arrays
del e125
del e135
del e145
del e155
