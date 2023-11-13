import sys
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import zscore
import scipy

def norm_usage(refit_usage):
    counts = refit_usage.sum(1)
    counts = np.ravel(counts)
    after = np.median(counts)
    counts += counts == 0
    counts = counts / after
    refit_usage = np.divide(refit_usage, counts[:, None])
    return(refit_usage)

def main(group, nhvg, K, topn):
    
    prefix= './'
    n=100
    ldsc = pd.read_csv(prefix+f'data/ldsc_{topn}/{group}_v{nhvg}_k{K}.log.bottomed_scaled_score_power_1_5.tsv'
                       , sep='\t')
    ldsc = ldsc.set_index('module_name')
    adata = sc.read(prefix+'../08_metacell/data/metacell_v2_cre.h5ad')
    
    ncells = adata.shape[0]
    ntraits = ldsc.shape[1]
    res = np.zeros((n, ncells, ntraits))

    for i in range(n):
        refit_shuf = pd.read_csv(prefix+
                                    f'data/{topn}/refit_shuf/cre_{group}_v{nhvg}_K{K}_{i}_usage.csv',
                                    index_col='barcode')
        refit_shuf = norm_usage(refit_shuf)
        refit_shuf.columns = ldsc.index
        res[i] = refit_shuf @ ldsc
        
        
    res = np.dstack(res)
    refit = pd.read_csv(prefix+
                        f'data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_usage.csv'
                        , index_col='barcode')
    refit = norm_usage(refit)
    refit.columns = ldsc.index
    m = refit @ ldsc
    
    pval = np.zeros((ncells, ntraits))
    zs = np.zeros((ncells, ntraits))
    for j in range(ntraits):
        zs[:, j] = zscore(np.append(np.reshape(m.to_numpy()[:, j], (-1, 1)), res[:,j, :], axis=1))[:,0]
        pval[:, j] = scipy.stats.norm.sf(zs[:, j])
        
    pval = pd.DataFrame(pval, index=adata.obs_names, columns=ldsc.columns)    
    zs = pd.DataFrame(zs, index=adata.obs_names, columns=ldsc.columns)    
    pval.to_csv(prefix+f'data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_pval.csv', index_label='barcode')
    zs.to_csv(prefix+f'data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_zscore.csv', index_label='barcode')
    
if __name__ == '__main__':
    main(group = sys.argv[1], 
         nhvg  = sys.argv[2], 
         K     = sys.argv[3],
         topn  = sys.argv[4])
