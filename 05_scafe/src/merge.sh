#!/bin/bash

run_tag="HCAJ_all"
genome="hg38.gencode_v32"
ulimit -n 10000

scafe.check.dependencies

scafe.workflow.cm.aggregate \
 	--overwrite=no \
 	--lib_list_path=data/all_ctss.tsv \
 	--max_thread=20 \
 	--genome=$genome \
 	--run_tag=$run_tag \
 	--run_outDir=scafe/aggregate_all/

scafe.tool.cm.merge \
	--overwrite=no \
	--lib_list_path=data/00_lib_list.HCAJ_all.txt \
	--lib_region_filter_criteria_path=data/00_lib_cutoff.optimal.HCAJ_all.txt \
	--aggregated_collapse_ctss_bed_path=scafe/aggregate_all/aggregate/HCAJ_all/bed/HCAJ_all.aggregate.collapse.ctss.bed.gz \
	--aggregated_unencoded_G_collapse_ctss_bed_path=scafe/aggregate_all/aggregate/HCAJ_all/bed/HCAJ_all.aggregate.unencoded_G.collapse.ctss.bed.gz \
	--force_retain_CRE_bed=data/pooled.catlas.screen.bed.gz \
	--genome=hg38.gencode_v32 \
	--outputPrefix=HCAJ_all_optimal \
	--outDir=scafe/merge_all/merge/

scafe.tool.cm.annotate \
	--overwrite=no \
	--tssCluster_bed_path=scafe/merge_all/merge/HCAJ_all_optimal/bed/HCAJ_all_optimal.tssCluster.merged.filtered.bed.gz \
	--tssCluster_info_path=scafe/merge_all/merge/HCAJ_all_optimal/log/HCAJ_all_optimal.tssCluster.merged.filtered.log.tsv \
	--min_CRE_count=10 \
	--pseudocount_correction=0 \
	--genome=hg38.gencode_v32 \
	--outputPrefix=HCAJ_all_optimal \
	--outDir=scafe/merge_all/annotate/


