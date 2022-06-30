nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process Plot_GOI_levels {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(metadata)
	    path(count_mat)
        path(gene_cnvan_key)
        path(GOI)
        val(model)
	    val(prefix)

    output:
        path("${prefix}*.pdf"), emit: GOI_expression_level_plots, optional: true

	script:
	"""
    Plot_GOI_levels.R \\
    --count_mat=${count_mat} \\
    --metadata=${metadata} \\
    --cnvan_key=${gene_cnvan_key} \\
    --GOI=${GOI} \\
    --model="${model}" \\
    --prefix=${prefix}
	"""
}