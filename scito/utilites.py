'''
Utilities for the package
'''

def count_collapser(data, bc,):
    '''
    Collapses counts for all antibodies for a given batch barcode. Outputs dense data structure
    :param data: sparse or dense data (AnnData or DataFrame) with counts
    :param bc: batch barcode
    :return: Numpy array of cells by batch BC.
    '''

    rel_data = data