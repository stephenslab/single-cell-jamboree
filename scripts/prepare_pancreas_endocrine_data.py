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

# Quality control - calculate QC covariates for all anndata objects
print(alldays.obs['day'].value_counts()) # number of cells per sample (day)

# #counts per cell
alldays.obs['n_counts'] = alldays.X.sum(1)
# #logcounts per cell
alldays.obs['log_counts'] = np.log(alldays.obs['n_counts'])
# #genes per cell
alldays.obs['n_genes'] = (alldays.X > 0).sum(1)
# mitochondrial gene fraction
mt_gene_mask = [gene.startswith('mt-') for gene in alldays.var_names]
mt_gene_index = np.where(mt_gene_mask)[0]
alldays.obs['mt_frac'] = alldays.X[:,mt_gene_index].sum(1) / alldays.X.sum(1)

#Sample quality plots
sc.pl.violin(alldays, ['n_counts', 'mt_frac'], groupby='day', size=1, log=False, cut=0)
#sc.pl.violin(alldays, 'mt_frac', groupby='day')

# Filter cells according to identified QC thresholds
print('Total number of cells: {:d}'.format(alldays.n_obs))
alldays = alldays[alldays.obs['mt_frac'] < 0.2]
print('Number of cells after MT filter: {:d}'.format(alldays.n_obs))

sc.pp.filter_cells(alldays, min_genes = 1200)
print('Number of cells after gene filter: {:d}'.format(alldays.n_obs))
#Filter genes:
print('Total number of genes: {:d}'.format(alldays.n_vars))

# Min 20 cells - filters out 0 count genes
sc.pp.filter_genes(alldays, min_cells=20)
print('Number of genes after cell filter: {:d}'.format(alldays.n_vars))

# Keep a copy of the raw, filtered data in a separate anndata object.
alldays_counts = alldays.copy()
