import os
import pandas as pd
df = pd.read_csv('data/cell_counts.csv', index_col=1)

rule all:
  input:
    expand('data/cellbender/{sample}.h5', sample=df.index)

rule cellbender:
  input:
    h5 = 'data/raw_h5/{sample}.h5'
  output:
    h5 = 'data/cellbender/{sample}.h5'
  log:
    'logs/cellbender/{sample}.log'
  params:
    cells_good = lambda wildcards: round(df.loc[int(wildcards.sample)].cells * 0.6),
    cells_include = lambda wildcards: round(df.loc[int(wildcards.sample)].cells * 2.5)
  shell:
    """
    /home/jonathan/data/programs/miniconda_202002/envs/CellBender/bin/cellbender remove-background \
        --input {input.h5} \
        --output {output.h5} \
        --cuda \
        --expected-cells {params.cells_good} \
        --total-droplets-included {params.cells_include} \
        --fpr 0.01 \
        --epochs 150 > {log} 2>&1
    """
