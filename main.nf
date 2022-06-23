#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { Reformat_data } from "./modules/Reformat_data"

workflow {
    // Retrieve and validate parameters
    assert params.gene_expression_matrix != null : "Please specify the `gene_expression_matrix` parameter"
    assert params.samplesheet != null : "Please specify the `samplesheet` parameter"
    gene_expression_matrix = file(params.gene_expression_matrix, checkIfExists: true)
    samplesheet = file(params.samplesheet, checkIfExists: true)
    prefix = params.prefix

    // start workflow
    Reformat_data(samplesheet, gene_expression_matrix, prefix)
}
