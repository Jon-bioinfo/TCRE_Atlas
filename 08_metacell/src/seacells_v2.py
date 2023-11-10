import sys
import scanpy as sc
import pandas as pd
import SEACells
from SEACells.core import summarize_by_SEACell

def main(ct, n_SEACells):
    ad = sc.read('../04_subcluster_v2/data/03_subcluster.h5ad')
    ad = ad[ad.obs['Narrow_Celltype'] == ct]
    
    build_kernel_on = 'X_sub_umap_harmony'
    n_waypoint_eigs = 10
    model = SEACells.core.SEACells(ad, 
                  build_kernel_on=build_kernel_on, 
                  n_SEACells=n_SEACells, 
                  n_waypoint_eigs=n_waypoint_eigs,
                  convergence_epsilon = 1e-5)
    
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    model.fit(min_iter=10, max_iter=200)
    ad.obs[['SEACell']].to_csv(f'data/seacells_v2/{ct}.csv.gz', index_label='barcode', index=True)
    
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])


