#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { Reformat_data } from "./modules/Reformat_data"
include { DESeq2_DGEA } from "./modules/DESeq2_DGEA"
include { Draw_volcano } from "./modules/Draw_volcano"
include { Plot_GOI_levels } from "./modules/Plot_GOI_levels"
include { variancePartition } from "./modules/variancePartition"

workflow {
    // Retrieve and validate parameters
    assert params.gene_expression_matrix != null : "Please specify the `gene_expression_matrix` parameter"
    assert params.samplesheet != null : "Please specify the `samplesheet` parameter"
    gene_expression_matrix = file(params.gene_expression_matrix, checkIfExists: true)
    samplesheet = file(params.samplesheet, checkIfExists: true)
    GOI = params.GOI ? file(params.GOI, checkIfExists: true) : []
    prefix = params.prefix
    model = params.model
    treat_col = params.treat_col
    pCutoff = params.pCutoff
    FCcutoff = params.FCcutoff

    // start workflow
    Reformat_data(samplesheet, gene_expression_matrix, prefix)
    DESeq2_DGEA(Reformat_data.out.metadata, Reformat_data.out.count_mat, Reformat_data.out.gene_cnvan_key, model, treat_col, prefix)
    Draw_volcano(DESeq2_DGEA.out.de_res.flatten(), GOI, pCutoff, FCcutoff, prefix)
    Plot_GOI_levels(Reformat_data.out.metadata, Reformat_data.out.count_mat, Reformat_data.out.gene_cnvan_key, GOI, model, prefix)
    variancePartition(Reformat_data.out.metadata, Reformat_data.out.count_mat, Reformat_data.out.gene_cnvan_key, model, prefix)
}

