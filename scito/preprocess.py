'''
Data pre-processing for SCITO-seq
'''


import os
import scanpy as sc
import numpy as np
import multiprocessing
import anndata
import pandas as pd

from sklearn.cluster import KMeans
from warnings import warn
from scipy import sparse
from scipy.stats import norm
from scito.utils import count_collapser, av_gene_expression, drop_assigner, drop_identifier

class ScitoFrame:
    '''
    Class storing the count matrix

    path (str): path to the h5 data
    '''

    def __init__(self, path=None, from_scanpy=None):

        if from_scanpy != None:
            self.adata = from_scanpy.copy()

        else:
            ext = os.path.splitext(path)[1].strip('.')
            self.ext = ext
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
                  kfunc="kmeans",
                  maxneighbor=100,
                  seed=33,
                  keep_input=False,
                  verbose=False):
        '''
        Function to assign droplets to sample id and detect singlets vs multiplets. Antibody counts are expected to be
        normalized and log scaled (e.g. using sc.pp.normalize_per_cell(), sc.pp.log1p())
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


        batches = self.adata.var_names.str.extract(r'(%s\d+)'%batchid_string).iloc[:,0].dropna().unique()
        nClust = n_clust if n_clust != None else len(batches)+1

        # extract only antibody counts
        ab_adata = self.adata[:,self.adata.var_names.str.contains(r'(%s\d+)'%batchid_string)]

        # collapse counts within batch
        batch_counts = np.transpose(np.array([count_collapser(ab_adata, bc) for bc in batches]))

        if keep_input:
            self.input = ab_adata
            if verbose:
                print("Keeping sparse matrix with antibody expression only. Target = self.input")
        ab_adata = None

        # create anndata object with collapsed counts per batch
        batch_adata = anndata.AnnData(X=sparse.csr_matrix(batch_counts),
                                      obs=self.adata.obs,
                                      var=pd.DataFrame(batches, columns=['batch']))
        batch_counts = None

        # Normalize and log scale data
        sc.pp.normalize_per_cell(batch_adata, counts_per_cell_after=1e4)
        if verbose:
            print("Keeping linear scale data for computing average")
        batch_adataNormLin = batch_adata.copy()
        sc.pp.log1p(batch_adata)

        # clustering functions
        # TODO: re-implement CLARA or other k-medoid
        if kfunc == "kmeans":
            if verbose:
                print("Performing k-means clustering")
            clust_model = KMeans(n_clusters=nClust, n_init=n_init, random_state=seed, n_jobs=int(multiprocessing.cpu_count()*0.6))
            clusters = clust_model.fit_predict(batch_adata.X.todense())

        else:
            print("ERROR: unknown clustering function. Select from kmeans, ")



        batch_adataNormLin.obs['batch_cluster'] = clusters
        batch_adata.obs['batch_cluster'] = clusters
        marker_dict = {batch_adataNormLin.var.keys()[0]: batch_adataNormLin.var['batch']} # for average expression of batch oligo

        av_batch_expr = av_gene_expression(batch_adataNormLin, marker_dict, gene_symbol_key='batch', partition_key='batch_cluster').iloc[:,:-1]

        # free up memory
        batch_adataNormLin = None

        if any(av_batch_expr.iloc[:,:-1]) == 0:
            warn("WARNING Cells with 0 counts exist as a cluster")

        # create a matrix to store classification result
        discrete = anndata.AnnData(X=sparse.csr_matrix(np.zeros(batch_adata.shape)),
                                   obs=self.adata.obs,
                                   var=pd.DataFrame(batches, columns=['batch']))


        # for each batch barcode, we will use the minimum cluster for fitting
        # NOTE: fitting normal distribution to normalized and log-transformed data. Thresholds will be more conservative
        # TODO implement nbinom fit
        for batch_name in av_batch_expr.index:
            values = batch_adata[:, batch_adata.var['batch'] == batch_name]
            values_use = values[values.obs['batch_cluster'] == np.argmin(av_batch_expr.loc[batch_name, :])]
            fitty = norm.fit(values_use.X.toarray())
            cutoff = np.quantile(norm.rvs(loc=fitty[0], scale=fitty[1], size=1000, random_state=seed),
                                 q=positiveQuantile)

            ox = [x[0] for x in np.argwhere(values.X > cutoff)]

            discrete.X[ox, int(values.var.index[0])] = 1
            if verbose:
                print("Cutoff for {}: {} reads".format(batch_name,
                                                       int(np.expm1(cutoff))))


        # assign whether drop is SNG, MTP or NEG
        n_positive = np.sum(discrete.X, axis=1)
        assignment = [drop_assigner(int(x)) for x in n_positive]
        batch_adata.obs['assignment'] = assignment

        if verbose:
            print("Assigning best guesses")

        # assign cells to HTO and get expression values of each
        best = [drop_identifier(a=batch_adata.X[x, :],
                                n_top=int(n_positive[x]),
                                bc_ids=batch_adata.var['batch'])
                for x in range(batch_adata.n_obs)]



        best_guess = [x['barcodes'] for x in best]
        best_exp = [x['expression'] for x in best]

        batch_adata.obs['best_guess'] = best_guess
        batch_adata.obs['expression'] = best_exp

        # Assemble some meta data
        n_cells_atLevel = [sum(n_positive == i).tolist()[0][0] for i in range(0,7)]
        n_cells_atLevel_df = pd.DataFrame({"cells_per_drop": ["{} cells per drop".format(x) for x in range(0,7)],
                                           "N_drops": n_cells_atLevel})


        self.meta = n_cells_atLevel_df

        self.n_positive = n_positive


        return batch_adata













