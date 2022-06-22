#!/usr/bin/bash

singularity exec \
-B in:/input_data \
-B out:/results \
R-modules.sif \
Rscript \
/opt/Rscripts/Reformat_pseudobulk.R \
--gene_counts=/input_data/salmon.merged.gene_counts_clip.tsv \
--samplesheet=/input_data/samplesheet_scrnaseq.csv \
--prefix="test" \
--results_dir=/results