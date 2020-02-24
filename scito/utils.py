'''
Utilities for the package
'''

import numpy as np
import scanpy as sc

def count_collapser(data, bc):
    '''
    Collapses counts for all antibodies for a given batch barcode. Outputs dense data structure
    :param data: sparse or dense data (AnnData or DataFrame) with counts
    :param bc: batch barcode
    :return: .collapsed = Numpy array of cells by batch BC.
            .n_ab = number of antibodies per batch (quality check)
    '''
    # TODO add handle for dense data
    rel_data = data[:,data.var.index.str.contains(r'(%s$)'%bc)]

    # count number of antibodies per batch
    assertrel_data.n_vars

    # collapse counts
    collapsed = rel_data.X.toarray().sum(axis=1).astype("int")

    return collapsed


def av_gene_expression(anndata, marker_dict, gene_symbol_key=None, partition_key='batch_cluster'):
    """ Copied from https://github.com/theislab/scanpy/issues/181 - posted by one of scanpy developers
    A function go get mean expressions of feature per cluster (class).
    Data should be normalized to reads per cell, but kept in the linear scale
    #
    # Inputs:
    #    anndata         - An AnnData object containing the data set and a partition
    #    marker_dict     - A dictionary with cell-type markers. The markers should be stores as anndata.var_names or
    #                      an anndata.var field with the key given by the gene_symbol_key input
    #    gene_symbol_key - The key for the anndata.var field with gene IDs or names that correspond to the marker
    #                      genes
    #    partition_key   - The key for the anndata.obs field where the cluster IDs are stored. The default is
    #                      'batch_cluster' """

    # Test inputs
    if partition_key not in anndata.obs.columns.values:
        print('KeyError: The partition key was not found in the passed AnnData object.')
        print('   Have you done the clustering? If so, please tell pass the cluster IDs with the AnnData object!')
        raise

    if (gene_symbol_key != None) and (gene_symbol_key not in anndata.var.columns.values):
        print('KeyError: The provided gene symbol key was not found in the passed AnnData object.')
        print('   Check that your cell type markers are given in a format that your anndata object knows!')
        raise

    if gene_symbol_key:
        gene_ids = anndata.var[gene_symbol_key]
    else:
        gene_ids = anndata.var_names

    clusters = set(anndata.obs[partition_key])
    n_clust = len(clusters)
    marker_exp = pd.DataFrame(columns=clusters)
    marker_exp['cell_type'] = pd.Series({}, dtype='str')
    marker_names = []

    i = 0
    for group in marker_dict:
        # Find the corresponding columns and get their mean expression in the cluster
        for gene in marker_dict[group]:
            ens_idx = np.in1d(gene_ids, gene)  # Note there may be multiple mappings
            if np.sum(ens_idx) == 0:
                continue
            else:
                anndata.obs[ens_idx[0]] = anndata.X[:, ens_idx].mean(1)  # works for both single and multiple mapping
                ens_idx = ens_idx[0]

            clust_marker_exp = anndata.obs.groupby(partition_key)[ens_idx].apply(np.mean).tolist()
            clust_marker_exp.append(group)
            marker_exp.loc[i] = clust_marker_exp
            marker_names.append(gene)
            i += 1

    # Replace the rownames with informative gene symbols
    marker_exp.index = marker_names

    return (marker_exp)


def drop_assigner(x):
    '''
    Assigns droplet to be empty ("NEG"), singlet ("SNG") or multiplet ("MTP"). Accepts numeric values
    :param x: numeric value - number of batches detected above the threshold for this droplet
    :return: Assignment for a droplet
    '''
    if x == 0:
        res = "NEG"
    elif x == 1:
        res = "SNG"
    elif x > 1:
        res = "MTP"
    else:
        assert (x >= 0), "ERROR: wrong assignment"

    return res

def drop_identifier(a, n_top, bc_ids):
    '''
    Function to get the batch ID and barcode expression of the best guesses
    :param a: array or list of barcode expression values per cell
    :param n_top: number of top assignments to return - depends on how many cells are detected per multiplet
    :param bc_ids: adata.var with batch name
    :return: dictionary:
                "barcodes" - barcode IDs for best guesses in order of descending expression values
                "expression" - expression values in same order as barcode IDs
    '''
    a = a.toarray()[0].tolist()
    a = [round(x, 2) for x in a]
    bc_ids = bc_ids.tolist()
    a_ord, bc_ids_ord = zip(*sorted(zip(a, bc_ids), reverse=True))

    top_a = a_ord[:n_top]
    top_bc = bc_ids_ord[:n_top]

    results = {"barcodes": top_bc,
               "expression": top_a}

    return results
