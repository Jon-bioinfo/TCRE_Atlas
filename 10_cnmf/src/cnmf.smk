n_iter = 10
iters = [x for x in range(n_iter)]
nhvg = [2000]
k = [x for x in range(50,250+1, 10)]
k_str = ' '.join([str(x) for x in k])
jobs_tot = 36
jobs = [x for x in range(1,jobs_tot+1)]
dataset = 'mc_n3350'
out_dir = 'data/output'
cnmf_bin = 'cnmf'
#final_k = 150

rule all:
  input:
    expand('{out_dir}/{dataset}_v{nhvg}/{dataset}_v{nhvg}.usages.k_{k}.dt_0_2.consensus.txt', out_dir = out_dir, dataset = dataset, nhvg = nhvg, k=k)

rule prepare:
  input:
    h5ad = 'data/mc_3350_total.h5ad'
  output:
    h5ad = '{out_dir}/{dataset}_v{nhvg}/cnmf_tmp/{dataset}_v{nhvg}.norm_counts.h5ad'
  log:
    'logs/prepare/{out_dir}/{dataset}_v{nhvg}.log'
  params:
    cnmf = cnmf_bin,
    name = '{dataset}_v{nhvg}',
    n_iter = n_iter,
    k_str = k_str
  shell:
    '''
    {params.cnmf} prepare --output-dir {wildcards.out_dir} --name {params.name} -c {input.h5ad} -k {params.k_str} --numgenes {wildcards.nhvg} --n-iter {params.n_iter} --seed 1 > {log} 2>&1
    '''

rule factorize:
  input:
    rules.prepare.output.h5ad
  output:
    check = '{out_dir}/{dataset}_v{nhvg}/checks/factorize/{dataset}_v{nhvg}_{job}.check'
  log:
    'logs/factorize/{out_dir}/{dataset}_v{nhvg}_{job}.log'
  params:
    cnmf = cnmf_bin,
    name = '{dataset}_v{nhvg}',
    jobs = jobs_tot
  shell:
    '''
    {params.cnmf} factorize --output-dir {wildcards.out_dir} --name {params.name} --worker-index {wildcards.job} --total-workers {params.jobs} > {log} 2>&1
    echo 'done' > {output}
    '''

rule combine:
  input:
    expand('{out_dir}/{dataset}_v{nhvg}/checks/factorize/{dataset}_v{nhvg}_{job}.check', 
	out_dir = out_dir, dataset = dataset, nhvg = nhvg, k = k, job = jobs)
  output:
    check = '{out_dir}/{dataset}_v{nhvg}/checks/combine/{dataset}_v{nhvg}.check'
  log:
    'logs/combine/{out_dir}/{dataset}_v{nhvg}.log'
  params:
    cnmf = cnmf_bin,
    name = '{dataset}_v{nhvg}'
  shell:
    '''
    {params.cnmf} combine --output-dir {wildcards.out_dir} --name {params.name} > {log} 2>&1
    {params.cnmf} k_selection_plot --output-dir {wildcards.out_dir} --name {params.name} >> {log} 2>&1
    echo 'done' > {output} 
    '''

rule consensus:
  input:
    '{out_dir}/{dataset}_v{nhvg}/checks/combine/{dataset}_v{nhvg}.check'
  output:
    '{out_dir}/{dataset}_v{nhvg}/{dataset}_v{nhvg}.usages.k_{k}.dt_0_2.consensus.txt'
  log:
    'logs/consensus/{out_dir}/{dataset}_v{nhvg}_{k}.log'
  params:
    cnmf = cnmf_bin,
    name = '{dataset}_v{nhvg}'
  shell:
    '''
    {params.cnmf} consensus --output-dir {wildcards.out_dir} --name {params.name} --components {wildcards.k} --local-density-threshold 0.2 --show-clustering > {log} 2>&1
    '''
