nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process Render_Report {
    publishDir "${out_dir}", mode: "$mode"
    stageInMode "copy"

    input:
        path(Report)
        path(metadata)
        path(de_res)
        path(volcano_plots)
        path(GOI_expression_level_plots)
        path(variancePartition)

    output:
        path("analysis/_build/html"), emit: Report

	script:
	"""
    jupyter-book build --all ${Report}
	"""
}