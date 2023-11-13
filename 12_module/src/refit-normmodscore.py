import sys
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

def main(group, nhvg, K, topn):
    prefix = './'

    ldsc = pd.read_csv(prefix+f'data/ldsc_{topn}/{group}_v{nhvg}_k{K}.log.bottomed_scaled_score_power_1_5.tsv', sep='\t')
    ldsc = ldsc.set_index('module_name')

    # The dataset for all cells in AnnData as counts.
    adata = sc.read(prefix+'../08_metacell/data/metacell_v2_cre.h5ad')

    refit_usage = pd.read_csv(prefix+f'data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_usage.csv', index_col='barcode')

    counts = refit_usage.sum(1)
    counts = np.ravel(counts)
    after = np.median(counts)
    counts += counts == 0
    counts = counts / after
    refit_usage = np.divide(refit_usage, counts[:, None])

    refit_usage.columns = ldsc.index
    m = refit_usage @ ldsc
    m.to_csv(prefix+f'data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_traits.csv', index_label='barcode')

if __name__ == '__main__':
    main(group = sys.argv[1], 
         nhvg  = sys.argv[2], 
         K     = sys.argv[3],
         topn  = sys.argv[4])
