import sys
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from sklearn.utils import shuffle
from cnmf import cNMF

def main(group, nhvg, K, topn, perm):
    print(perm)
    prefix = './'
    cnmf_obj = cNMF(output_dir=prefix+f'data/cNMF/', name=f'{group}_v{nhvg}')
    norm_hvg = sc.read(cnmf_obj.paths['normalized_counts']) 
    
    hvgs = list(norm_hvg.var_names)
    adata = sc.read(prefix+'../08_metacell/data/metacell_v2_cre.h5ad')
    adata = adata[:, hvgs]
    adata.var['sum'] = np.ravel(adata.X.sum(0))
    adata.var['expr_bin'] = pd.qcut(adata.var['sum'], 5, labels=list('ABCDE')) 
    levels = adata.var.expr_bin.unique()
    
    split_adata = [0] * len(levels)
    for i in range(len(levels)):
        split_adata[i] = adata[:, adata.var.expr_bin == levels[i]]
        split_adata[i].X = shuffle(split_adata[i].X, random_state=int(perm))
        
    x = ad.concat(split_adata, axis=1)
    x = x[:, adata.var_names]
    x.obs = adata.obs
    x.write(prefix+f'data/{topn}/shuffle/cre_{group}_v{nhvg}_K{K}/{perm}.h5ad')
    
if __name__ == '__main__':
    main(group = sys.argv[1], 
         nhvg  = sys.argv[2], 
         K     = sys.argv[3],
         topn  = sys.argv[4],
         perm  = sys.argv[5])
