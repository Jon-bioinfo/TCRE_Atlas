import sys
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from cnmf import cNMF

def main(group, nhvg, K, i, topn, ldthresh=0.2):
    prefix = './'
    cnmf_obj = cNMF(output_dir=prefix+f'data/cNMF/', name=f'{group}_v{nhvg}')
    norm_hvg = sc.read(cnmf_obj.paths['normalized_counts']) 
    usage_unnormed = pd.read_csv(cnmf_obj.paths['consensus_usages__txt'] % (int(K), str(ldthresh).replace('.','_')), sep='\t', index_col=0) 
    refit_spectra = cnmf_obj.refit_spectra(norm_hvg.X, usage_unnormed.astype(np.float32))
    
    hvgs = list(norm_hvg.var_names)
    adata = sc.read(prefix+f'data/{topn}/shuffle/cre_{group}_v{nhvg}_K{K}/{i}.h5ad')
    adata = adata[:, hvgs]
    sc.pp.scale(adata, zero_center=False)

    refit_usage = cnmf_obj.refit_usage(adata.X, refit_spectra)
    refit_usage = pd.DataFrame(refit_usage, index=adata.obs_names)
    refit_usage.columns = ['M'+"%03d" % x for x in  refit_usage.columns.values.astype(int) + 1]

    refit_usage.to_csv(prefix+f'data/{topn}/refit_shuf/cre_{group}_v{nhvg}_K{K}_{i}_usage.csv', index_label='barcode')

if __name__ == '__main__':
    main(group = sys.argv[1], 
         nhvg  = sys.argv[2], 
         K     = sys.argv[3],
         topn  = sys.argv[4],
         i     = sys.argv[5])
