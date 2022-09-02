nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process Draw_volcano {
    //Packages dependencies
    //conda "conda-forge::r-base=4.1.2 conda-forge::r-docopt=0.7.1 conda-forge::r-conflicted=1.1.0 conda-forge::dplyr=1.3.1 bioconda::bioconductor-enhancedvolcano=1.12.0"
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(de_res)
	    path(GOI)
	    val(pCutoff)
	    val(FCcutoff)
	    val(prefix)

    output:
        path("${prefix}*.pdf"), emit: volcano_plots, optional: true

	script:
    contrast_name = de_res.getBaseName()
	contrast_name = contrast_name.replaceAll(prefix + "_de_res_", "")

	"""
    Draw_volcano.R \\
    --de_res=${de_res} \\
    --GOI=${GOI} \\
    --prefix=${prefix} \\
    --contrast_name=${contrast_name} \\
    --pCutoff=${pCutoff} \\
    --FCcutoff=${FCcutoff}
	"""
}