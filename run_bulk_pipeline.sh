#!/bin/bash
# Run custom bulk RNA-seq pipeline

nextflow run ./main.nf -profile cluster \
-w /home/marteau/myScratch/nf-work-dir/nf-bulk-pipeline/work \
-resume