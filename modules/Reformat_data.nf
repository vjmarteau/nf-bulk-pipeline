nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process Reformat_data {
    //Packages dependencies
    //conda "conda-forge::r-base=4.1.2 conda-forge::r-docopt=0.7.1 conda-forge::r-conflicted=1.1.0 conda-forge::dplyr=1.3.1"
    //container "./envs/nf-bulk-r.sif" // Can be specified in .config file
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(samplesheet)
	    path(gene_expression_matrix)
	    val(prefix)

    output:
        path("${prefix}_metadata.tsv"), emit: metadata
        path("${prefix}_count_mat.tsv"), emit: count_mat
        path("${prefix}_gene_cnvan_key.tsv"), emit: gene_cnvan_key


	script:
	"""
    Reformat_data.R \\
    --gene_counts=${gene_expression_matrix} \\
    --samplesheet=${samplesheet} \\
    --prefix=${prefix}
	"""
}
