'''
Utilities for the package
'''

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


def av_gene_expression(anndata, marker_dict, gene_symbol_key=None, partition_key='louvain_r1'):
    """ Copied from https://github.com/theislab/scanpy/issues/181 - posted by one of scanpy devs
    A function go get mean z-score expressions of feature per cluster (class)
    #
    # Inputs:
    #    anndata         - An AnnData object containing the data set and a partition
    #    marker_dict     - A dictionary with cell-type markers. The markers should be stores as anndata.var_names or
    #                      an anndata.var field with the key given by the gene_symbol_key input
    #    gene_symbol_key - The key for the anndata.var field with gene IDs or names that correspond to the marker
    #                      genes
    #    partition_key   - The key for the anndata.obs field where the cluster IDs are stored. The default is
    #                      'louvain_r1' """

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

    clusters = anndata.obs[partition_key].cat.categories
    n_clust = len(clusters)
    marker_exp = pd.DataFrame(columns=clusters)
    marker_exp['cell_type'] = pd.Series({}, dtype='str')
    marker_names = []

    z_scores = sc.pp.scale(anndata, copy=True)

    i = 0
    for group in marker_dict:
        # Find the corresponding columns and get their mean expression in the cluster
        for gene in marker_dict[group]:
            ens_idx = np.in1d(gene_ids, gene)  # Note there may be multiple mappings
            if np.sum(ens_idx) == 0:
                continue
            else:
                z_scores.obs[ens_idx[0]] = z_scores.X[:, ens_idx].mean(1)  # works for both single and multiple mapping
                ens_idx = ens_idx[0]

            clust_marker_exp = z_scores.obs.groupby(partition_key)[ens_idx].apply(np.mean).tolist()
            clust_marker_exp.append(group)
            marker_exp.loc[i] = clust_marker_exp
            marker_names.append(gene)
            i += 1

    # Replace the rownames with informative gene symbols
    marker_exp.index = marker_names

    return (marker_exp)
