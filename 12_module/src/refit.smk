import pandas as pd
dat = pd.read_csv('data/refit_params.tsv', sep='\t')
dat2 = pd.merge(dat, pd.DataFrame({'perm':range(100)}), how='cross')

rule all:
  input:
    expand('data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_usage.csv', zip, group = dat.group, nhvg=dat.nhvg, K=dat.K, topn = dat.topn),
    expand('data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_traits.csv', zip, group = dat.group, nhvg=dat.nhvg, K=dat.K, topn = dat.topn),
    expand('data/{topn}/shuffle/cre_{group}_v{nhvg}_K{K}/{i}.h5ad', zip, group = dat2.group, nhvg=dat2.nhvg, K=dat2.K, i=dat2.perm, topn = dat2.topn),
    expand('data/{topn}/refit_shuf/cre_{group}_v{nhvg}_K{K}_{i}_usage.csv', zip, group = dat2.group, nhvg=dat2.nhvg, K=dat2.K, i=dat2.perm, topn = dat2.topn),
    expand('data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_zscore.csv', zip, group = dat.group, nhvg=dat.nhvg, K=dat.K, topn = dat.topn),
    expand('data/{topn}/refit-convert/cre_{group}_v{nhvg}_K{K}_pval.h5ad', zip, group = dat.group, nhvg=dat.nhvg, K=dat.K, topn = dat.topn)

rule refit:
  input:
    h5ad = 'data/cNMF/{group}_v{nhvg}/cnmf_tmp/{group}_v{nhvg}.norm_counts.h5ad'
  output:
    usage = 'data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_usage.csv'
  log:
    'logs/{topn}/refit/cre_{group}_v{nhvg}_K{K}.log'
  shell:
    '''
    python src/refit.py {wildcards.group} {wildcards.nhvg} {wildcards.K} {wildcards.topn} > {log} 2>&1
    '''

rule refit_norm:
  input:
    ldsc = 'data/ldsc_{topn}/{group}_v{nhvg}_k{K}.log.bottomed_scaled_score_power_1_5.tsv',
    usage = 'data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_usage.csv'  
  output:
    traits = 'data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_traits.csv'
  log:
    'logs/{topn}/refit2/cre_{group}_v{nhvg}_K{K}.log'
  shell:
    '''
    python src/refit-normmodscore.py {wildcards.group} {wildcards.nhvg} {wildcards.K} {wildcards.topn} > {log} 2>&1
    '''
    
rule convert:
  input:
    'data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_zscore.csv'
  output:
    'data/{topn}/refit-convert/cre_{group}_v{nhvg}_K{K}_pval.h5ad'
  log:
    'logs/{topn}/refit-convert/cre_{group}_v{nhvg}_K{K}.log'
  shell:
    '''
    python src/refit-convert.py {wildcards.group} {wildcards.nhvg} {wildcards.K} {wildcards.topn} > {log} 2>&1
    '''

rule shuffle:
  input:
    'data/cNMF/{group}_v{nhvg}/cnmf_tmp/{group}_v{nhvg}.norm_counts.h5ad'
  output:
    'data/{topn}/shuffle/cre_{group}_v{nhvg}_K{K}/{i}.h5ad'
  log:
    'logs/{topn}/shuffle/cre_{group}_v{nhvg}_K{K}/{i}.log'
  threads: 2
  shell:
    '''
    python src/shuffle.py {wildcards.group} {wildcards.nhvg} {wildcards.K} {wildcards.topn} {wildcards.i} > {log} 2>&1
    '''
    
rule refit_shuf:
  input:
    'data/{topn}/shuffle/cre_{group}_v{nhvg}_K{K}/{i}.h5ad',
    'data/cNMF/{group}_v{nhvg}/cnmf_tmp/{group}_v{nhvg}.norm_counts.h5ad'
  output:
    'data/{topn}/refit_shuf/cre_{group}_v{nhvg}_K{K}_{i}_usage.csv'
  log:
    'logs/{topn}/refit_shuf/cre_{group}_v{nhvg}_K{K}_{i}.log'
  shell:
    '''
    python src/refit-shuf.py {wildcards.group} {wildcards.nhvg} {wildcards.K} {wildcards.topn} {wildcards.i} > {log} 2>&1
    '''

rule pval:
  input:
    ldsc = 'data/ldsc_{topn}/{group}_v{nhvg}_k{K}.log.bottomed_scaled_score_power_1_5.tsv',
  output:
    'data/{topn}/refit/cre_{group}_v{nhvg}_K{K}_zscore.csv'
  log:
    'logs/{topn}/refit_pval/cre_{group}_v{nhvg}_K{K}.log'
  shell:
    '''
    python src/refit-permute.py {wildcards.group} {wildcards.nhvg} {wildcards.K} {wildcards.topn} > {log} 2>&1
    '''
