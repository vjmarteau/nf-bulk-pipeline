nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process variancePartition {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(metadata)
	    path(count_mat)
        path(gene_cnvan_key)
        val(model)
	    val(prefix)

    output:
        path("${prefix}_variancePartition.pdf"), emit: variancePartition_plot
        path("${prefix}_variancePartition.tsv"), emit: variancePartition

	script:
	"""
    variancePartition.R \\
    --count_mat=${count_mat} \\
    --metadata=${metadata} \\
    --cnvan_key=${gene_cnvan_key} \\
    --model=${model} \\
    --prefix=${prefix}
	"""
}
