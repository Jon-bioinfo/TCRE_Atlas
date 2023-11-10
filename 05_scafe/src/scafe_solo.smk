import pandas as pd
df = pd.read_csv('data/raw_paths.tsv', sep='\t')
df = df.set_index('id')
samples = df.index

rule all:
  input:
    expand('scafe/solo/count/{sample}/matrix/matrix.mtx', sample = samples)

rule scafe_solo:
  input:
    mtx_barcodes = lambda wildcards: df.loc[wildcards.sample].bc,
    bam = lambda wildcards: df.loc[wildcards.sample].bam
  output:
    mtx = 'scafe/solo/count/{sample}/matrix/matrix.mtx',
    unen = 'scafe/solo/bam_to_ctss/{sample}/bed/{sample}.unencoded_G.collapse.ctss.bed.gz',
    ctss = 'scafe/solo/remove_strand_invader/{sample}/bed/{sample}.pass.ctss.bed.gz'
  log:
    'logs/scafe_solo/{sample}.log'
  threads: 6
  shell:
    """
    scafe.check.dependencies >> {log} 2>&1

    ulimit -n 5000 >> {log} 2>&1

    scafe.workflow.sc.solo \
        --overwrite=no \
        --run_bam_path={input.bam} \
        --run_cellbarcode_path={input.mtx_barcodes} \
        --detect_TS_oligo=auto \
        --genome=hg38.gencode_v32 \
        --run_tag={wildcards.sample} \
        --run_outDir=scafe/solo \
        --max_thread={threads} \
        >> {log} 2>&1
    """

