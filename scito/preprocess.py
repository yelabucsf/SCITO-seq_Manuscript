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
from scipy.stats import norm, nbinom
import statsmodels.api as sm
from scito.utils import count_collapser, av_gene_expression, drop_assigner, drop_identifier

class ScitoFrame:
    '''
    Class storing the count matrix
    :param path: path to count matrix
    :param from_scanpy: bool - is it going to be imported from already available scanpy AnnData? Use if path=None
    :param adata: give variable name used for storing scanpy AnnData. Use only if from_scanpy=True
    :return:
    '''

    def __init__(self, path=None, from_scanpy=False, adata=None):

        if from_scanpy:
            self.adata = adata.copy()

        else:
            ext = os.path.splitext(path)[1].strip('.')
            self.ext = ext
            # TODO: add other format options (Kallisto?)
            if ext == "h5":
                self.adata = sc.read_10x_h5(path, gex_only=False)
            else:
                self.adata = None
                print("ERROR, unknown data format")


    def resolveMux(self,
                  batchid_string="barcode",
                  positiveQuantile=0.99,
                  n_clust=None,
                  n_init=100,
                  kfunc="kmeans",
                  maxneighbor=100,
                  distr_fit="nbinom",
                  seed=33,
                  keep_input=False,
                  collapse=True,
                  verbose=False):
        '''
        Function to assign droplets to sample id and detect singlets vs multiplets. Antibody counts are expected to be
        normalized and log scaled (e.g. using sc.pp.normalize_per_cell(), sc.pp.log1p())
        :param batchid_string: string. Identifies batch barcode. Default: "barcode"
        :param positiveQuantile: float. The quantile of inferred 'negative' distribution for each batch -
                over which the cell is considered 'positive'. Default is 0.99
        :param n_clust: int. Initial number of clusters for batches.
                Default is the # of batch oligo names + 1 (to account for negatives)
        :param n_init: int. value for k-means clustering (for kfunc = "kmeans"). 100 by default
        :param kfunc: string. Clustering function for initial batch grouping. Default and only available now is "kmeans"
        :param maxneighbor: int. Max number of neighbors per CLARANS cluster, for kfunc = "clarans" (irrelevant for now)
        :param distr_fit: str. Which distribution to fit. Accepts ("nbinom", "norm"). Default: "nbinom"
        :param seed: int. Sets the random seed.
        :param keep_input: bool. keep input previous step of data analysis
        :param collapse: bool. Sum all antibody counts per batch
        :param verbose: Chatty
        :return: anndata object split to sample id's and marked as singlets or multiplets
        '''
        if distr_fit == "norm":
            warn("WARNING data fit to normal distribution will be deprecated in the next versions")

        # extract only antibody counts
        ab_adata = self.adata[:,self.adata.var_names.str.contains(r'(%s\d+)'%batchid_string)]

        if collapse:
            batches = self.adata.var_names.str.extract(r'(%s\d+)' % batchid_string).iloc[:, 0].dropna().unique()
            num_ab = set([ab_adata[:, ab_adata.var_names.str.contains(r'(%s$)' % x)].n_vars for x in batches])
            # collapse counts within batch
            batch_counts = np.transpose(np.array([count_collapser(ab_adata, bc) for bc in batches]))

            # count number of clusters
            nClust = n_clust if n_clust != None else len(batches) + 1
            # create anndata object with collapsed counts per batch
            batch_adata = anndata.AnnData(X=sparse.csr_matrix(batch_counts),
                                          obs=self.adata.obs,
                                          var=pd.DataFrame(batches, columns=['batch']))
            batch_counts = None

        else:
            batchid_string = None
            batches = ab_adata.var_names
            num_ab = {1} # each ab is a batch
            batch_adata = ab_adata.copy()
            # count number of clusters
            nClust = len(batch_adata.var_names)
            batch_adata.var = pd.DataFrame(batches, columns=['batch'])

        # test that every batch has same number of anitbodies

        # TODO: try to catch this exception and look for other identifiers
        assert (len(num_ab) == 1), "ERROR: different number of antibodies per batch. Program exit"

        # Normalize and log scale data
        if distr_fit == "nbinom":
            batch_adataRAW = batch_adata.copy()
        sc.pp.normalize_per_cell(batch_adata, counts_per_cell_after=1e4, min_counts=0)
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
        marker_dict = {batch_adataNormLin.var.keys()[0]: batch_adataNormLin.var['batch']} # for average expression of batch oligo

        av_batch_expr = av_gene_expression(batch_adataNormLin, marker_dict, gene_symbol_key='batch', partition_key='batch_cluster').iloc[:,:-1]

        if distr_fit == "norm":
        # free up memory
            batch_adata.obs['batch_cluster'] = clusters
            batch_adataNormLin = None

        elif distr_fit == "nbinom":
            batch_adataRAW.obs['batch_cluster'] = clusters
            batch_adata = batch_adataRAW.copy()
            batch_adataNormLin = None
            batch_adataRAW = None
        else:
            print("ERROR: unknown distribution to fit. Choose from 'norm', 'nbinom' ")

        if any(av_batch_expr.iloc[:,:-1]) == 0:
            warn("WARNING Cells with 0 counts exist as a cluster")

        # create a matrix to store classification result
        discrete = anndata.AnnData(X=sparse.csr_matrix(np.zeros(batch_adata.shape)),
                                   obs=self.adata.obs,
                                   var=pd.DataFrame(batches, columns=['batch']))

        # blank adata for resolved drops
        if collapse:
            result = anndata.AnnData(X=sparse.csr_matrix(np.zeros((1, list(num_ab)[0]))),
                                     var=pd.DataFrame(
                                         index=set([x[0] for x in ab_adata.var_names.str.split(pat=batchid_string)])),
                                     obs=None)
        else:
            result = anndata.AnnData(X=sparse.csr_matrix(np.zeros((1, len(batches)))),
                                     var=pd.DataFrame(
                                         index=batches),
                                     obs=None)
        # for each batch barcode, we will use the minimum cluster for fitting
        # NOTE: fitting normal distribution to normalized and log-transformed data. Thresholds will be more conservative
        for batch_name in av_batch_expr.index:
            values = batch_adata[:, batch_adata.var['batch'] == batch_name]
            values_use = values[values.obs['batch_cluster'] == np.argmin(av_batch_expr.loc[batch_name, :])]

            if distr_fit == "norm":
                fitty = norm.fit(values_use.X.toarray())
                cutoff = np.quantile(norm.rvs(loc=fitty[0], scale=fitty[1], size=1000, random_state=seed),
                                     q=positiveQuantile)
                if verbose:
                    print("Cutoff for {}: {} reads".format(batch_name,
                                                           int(np.expm1(cutoff))))

            elif distr_fit == "nbinom":
                endog = values_use.X.toarray()
                exog = [1] * len(values_use.X.toarray())
                fitty = sm.NegativeBinomial(endog, exog).fit(disp=False)
                mu = np.e ** fitty._results.params[0]
                n = 1 / fitty._results.params[1]
                p = n / (n + mu)
                cutoff = np.quantile(nbinom.rvs(n=n, p=p, size=1000, random_state=seed),
                                     q=positiveQuantile)
                if verbose:
                    print("Cutoff for {}: {} reads".format(batch_name,
                                                           int(cutoff)))



            ox = [x[0] for x in np.argwhere(values.X > cutoff)]
            res = ab_adata[ox, ab_adata.var_names.str.contains(r'(%s$)' % batch_name)]
            res.var_names = [x[0] for x in res.var_names.str.split(pat=batchid_string)]
            res.var.drop(res.var.columns, axis=1, inplace=True)
            res.obs['batch_name'] = batch_name

            # resolved adata
            result = result.concatenate(res, join='outer', index_unique=None)


            discrete.X[ox, int(values.var.index[0])] = 1


        result = result[1:,:]


        if keep_input:
            self.input = ab_adata
            if verbose:
                print("Keeping sparse matrix with antibody expression only. Target = self.input")
        ab_adata = None

        del(self.adata)

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

        # save drop assignment (for HTO demuxing only
        self.drop_assign = batch_adata

        # some meta data
        self.meta = n_cells_atLevel_df
        self.n_positive = n_positive
        self.discrete = discrete



        # demultiplexed and resolved bcs
        return result















