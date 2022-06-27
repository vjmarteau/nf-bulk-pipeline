nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process DESeq2_DGEA {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(metadata)
	    path(count_mat)
        path(gene_cnvan_key)
        val(model)
        val(treat_col)
	    val(prefix)

    output:
        path("${prefix}_de_res_*.tsv"), emit: de_res


	script:
	"""
    DESeq2_DGEA.R \\
    --count_mat=${count_mat} \\
    --metadata=${metadata} \\
    --cnvan_key=${gene_cnvan_key} \\
    --model=${model} \\
    --treat_col=${treat_col} \\
    --prefix=${prefix}
	"""
}
