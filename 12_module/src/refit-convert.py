import sys
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

def main(group, nhvg, K, topn):
    prefix = './'
    
    for score in ['usage', 'pval','zscore']: 
    
        adata = sc.read(prefix+'data/HCAJ_Gene_CRE.h5ad')
        m = pd.read_csv(prefix+f'data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_{score}.csv', index_col=0)
        m = m[sorted(m.columns)]
        if score == 'pval':
            m = -np.log10(m)
        
        cols = ['Broad_Celltype','Narrow_Celltype']
        x = ad.AnnData(adata.X, var = adata.var, obs = adata.obs[cols].join(m), obsm = adata.obsm)

        x.write(prefix+f'data/{topn}/refit-convert/cre_{group}_v{nhvg}_K{K}_{score}_Gene_CRE.h5ad')
        x = x[:, 0]
        x.write(prefix+f'data/{topn}/refit-convert/cre_{group}_v{nhvg}_K{K}_{score}.h5ad')

if __name__ == '__main__':
    main(group = sys.argv[1], 
         nhvg  = sys.argv[2], 
         K     = sys.argv[3],
         topn  = sys.argv[4])