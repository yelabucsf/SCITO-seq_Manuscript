'''
Data pre-processing for SCITO-seq
'''


import os
import scanpy as sc
import numpy as np
import multiprocessing

from sklearn.cluster import KMeans
from pyclustering.cluster.clarans import clarans
from utils import count_collapser

class ScitoFrame:
    '''
    Class storing the count matrix

    path (str): path to the h5 data
    '''

    def __init__(self, path):

        ext = os.path.splitext(path)[1].strip('.')

# TODO: add other format options (Kallisto?)
        if ext == "h5":
            self.adata = sc.read_10x_h5(path, gex_only=False)

        else:
            self.adata = None
            print("ERROR, unknown data format")



    def detectMux(self,
                  batchid_string="barcode",
                  positiveQuantile=0.99,
                  n_clust=None,
                  n_init=100,
                  kfunc="clarans",
                  maxneighbor=100,
                  seed=42,
                  keep_input=False,
                  verbose=False):
        '''
        Function to assign droplets to sample id and detect singlets vs multiplets
        :param batchid_string: string identifying batch barcode. Default: "barcode"
        :param positiveQuantile: The quantile of inferred 'negative' distribution for each batch -
                over which the cell is considered 'positive'. Default is 0.99
        :param n_clust: Initial number of clusters for batches.
                Default is the # of batch oligo names + 1 (to account for negatives)
        :param n_init: value for k-means clustering (for kfunc = "kmeans"). 100 by default
        :param kfunc: Clustering function for initial batch grouping. Default and only available now is "kmeans"
        :param maxneighbor: Max number of neighbors per CLARANS cluster, for kfunc = "clarans" (irrelevant for now)
        :param seed: Sets the random seed.
        :param keep_input: keep input previous step of data analysis
        :param verbose: Chatty
        :return: anndata object split to sample id's and marked as singlets or multiplets
        '''

        # TODO: random seed

        batches = self.adata.var.index.str.extract(r'(%s\d+)'%batchid_string).iloc[:,0].dropna().unique()
        nClust = n_clust if n_clust != None else len(batches)+1

        # extract only antibody counts
        ab_adata = self.adata[:,self.adata.var.index.str.contains(r'(%s\d+)'%batchid_string)]

        if keep_input:
            self.input = ab_adata
        ab_adata = None

        # collapse counts within batch
        batch_counts = np.transpose(np.array([count_collapser(ab_adata, bc) for bc in batches]))

        # clustering functions
        # TODO: re-implement CLARA or other k-medoid
        if kfunc == "kmeans":
            clust_model = KMeans(n_clusters=nClust, n_init=n_init, random_state=seed, n_jobs=int(multiprocessing.cpu_count()*0.6))
            clusters = clust_model.fit_predict(batch_counts)

        else:
            print("ERROR: unknown clustering function. Select from kmeans, ")

        self.adata.obs['batch_cluster'] = clusters







