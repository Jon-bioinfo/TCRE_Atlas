import pandas as pd

samples = pd.read_csv('data/count_data.tsv', sep='\t', header=None)
head = ['projectID', 'sampleID', 'CB_ctss_bed_path', 'cellBarcode_list_path']
samples.columns = head
samples['name'] = samples['projectID'] +'.'+ samples['sampleID']
samples = samples.set_index('name')
srr = samples.index
                  
rule all:
  input:
    expand('scafe/merge_all/count/{sample}/matrix/matrix.mtx', sample = srr)

rule count:
  input:
    countRegion = 'scafe/merge_all/annotate/HCAJ_all_optimal/bed/HCAJ_all_optimal.CRE.coord.bed.gz',
    ctssBed = 'scafe/solo/bam_to_ctss/{sample}/bed/{sample}.CB.ctss.bed.gz',
    barcodes = lambda wildcards: samples.loc[wildcards.sample]['cellBarcode_list_path'],
    #scope = 'scafe/merge_all/mask/stranded_count_mask.bed.gz',
    scope = 'scafe/merge_all/annotate/HCAJ_all_optimal/bed/HCAJ_all_optimal.cluster.coord.bed.gz'
  output:
    mtx = 'scafe/merge_all/count/{sample}/matrix/matrix.mtx'
  log:
    'logs/count/{sample}.log'
  shell:
    """
    scafe.tool.sc.count \
	--countRegion_bed_path={input.countRegion} \
	--ctss_bed_path={input.ctssBed} \
        --ctss_scope_bed_path={input.scope} \
	--cellBarcode_list_path={input.barcodes} \
	--outputPrefix={wildcards.sample} \
	--genome=hg38.gencode_v32 \
	--outDir=scafe/merge_all/count \
        > {log} 2>&1
    """
