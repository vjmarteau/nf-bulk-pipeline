#!/bin/bash
# Run custom bulk RNA-seq pipeline

nextflow run ./main.nf \
-w /path/to/workDir/work \
--outdir /path/to/resDir/results \
--samplesheet /path/to/input/data/samplesheet.csv \
--gene_expression_matrix /path/to/input/data/gene_expression_matrix.tsv \
-resume