#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { Reformat_data } from "./modules/Reformat_data"
include { DESeq2_DGEA } from "./modules/DESeq2_DGEA"

workflow {
    // Retrieve and validate parameters
    assert params.gene_expression_matrix != null : "Please specify the `gene_expression_matrix` parameter"
    assert params.samplesheet != null : "Please specify the `samplesheet` parameter"
    gene_expression_matrix = file(params.gene_expression_matrix, checkIfExists: true)
    samplesheet = file(params.samplesheet, checkIfExists: true)
    prefix = params.prefix
    model = params.model
    treat_col = params.treat_col

    // start workflow
    Reformat_data(samplesheet, gene_expression_matrix, prefix)
    DESeq2_DGEA(Reformat_data.out.metadata, Reformat_data.out.count_mat, Reformat_data.out.gene_cnvan_key, model, treat_col, prefix)
}
