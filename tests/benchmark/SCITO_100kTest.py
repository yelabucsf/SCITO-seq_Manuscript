#!usr/bin/env python

from scito import preprocess
import seaborn as sns
import scanpy as sc
import pandas as pd

path = '/Users/antonogorodnikov/Documents/Work/DataSci/SCITO-seq/tests/100k_pbmc_filtered_feature_bc_matrix.h5'
verbose = True

adata = sc.read_10x_h5(path, gex_only=False)
corr_var = pd.read_csv("/Users/antonogorodnikov/Documents/Work/DataSci/SCITO-seq/tests/var.csv", index_col=0)

adata.var = corr_var
lol = preprocess.ScitoFrame(from_scanpy=True, adata=adata)
lol_2 = lol.detectMux(batchid_string="barcode",
                        positiveQuantile=0.99,
                        n_clust=None,
                        n_init=100,
                        kfunc="kmeans",
                        maxneighbor=100,
                        seed=42,
                        keep_input=True,
                        verbose=True)