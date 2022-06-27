#!/usr/bin/env Rscript
'
Usage:
  Draw_volcano.R --de_res=<de_res> --GOI=<GOI> --prefix=<prefix> [options]

Mandatory arguments:
  --de_res=<de_res>           TopTable from DESeq2 in TSV format
  --GOI=<GOI>                 Genes of interest in txt format
  --prefix=<prefix>           Prefix for output filenames

Optional arguments:
  --contrast_name=<contrast_name>
  --pCutoff=<pCutoff>         Cut-off for statistical significance [default: 0.05]
  --FCcutoff=<FCcutoff>       Cut-off for absolute log2 fold-change [default: 2]
  --resDir=<resDir>           Output directory [default: ./]
' -> doc

library(conflicted)
library(docopt)
arguments <- docopt(doc, version = "0.1")
print(arguments)

library(tidyverse)
library(EnhancedVolcano)

# Load parameters
de_res <- read_tsv(arguments$de_res)
GOI <- read_lines(arguments$GOI)
pCutoff <- as.numeric(arguments$pCutoff)
FCcutoff <- as.numeric(arguments$FCcutoff)
prefix <- arguments$prefix
contrast_name <- as.character(arguments$contrast_name)
resDir <- arguments$resDir


# Make volcano plots using "EnhancedVolcano" package

message("Draw volcano plots ...")

p <- EnhancedVolcano(
  toptable = de_res,
  lab = NA,
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = pCutoff,
  FCcutoff = FCcutoff,
  title = paste0(prefix, "_", contrast_name, "_volcano_plot"),
  caption = paste0("fold change cutoff: ", FCcutoff, "; p-value cutoff: ", pCutoff)
)

ggsave(file.path(resDir, paste0(prefix, "_", contrast_name, "_volcano_plot.pdf")), plot = p, width = 297, height = 210, units = "mm")

p <- EnhancedVolcano(
  toptable = de_res,
  lab = NA,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = pCutoff,
  FCcutoff = FCcutoff,
  title = paste0(prefix, "_", contrast_name, "_volcano_plot"),
  caption = paste0("fold change cutoff: ", FCcutoff, "; adj.p-value cutoff: ", pCutoff)
)

ggsave(file.path(resDir, paste0(prefix, "_", contrast_name, "_volcano_padj.pdf")), plot = p, width = 297, height = 210, units = "mm")

p <- EnhancedVolcano(
  toptable = de_res,
  lab = de_res$symbol,
  selectLab = GOI,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = pCutoff,
  FCcutoff = FCcutoff,
  drawConnectors = TRUE,
  title = paste0(prefix, "_", contrast_name, "_volcano_plot_GOI"),
  caption = paste0("fold change cutoff: ", FCcutoff, "; adj.p-value cutoff: ", pCutoff)
)

ggsave(file.path(resDir, paste0(prefix, "_", contrast_name, "_volcano_padj_GoI.pdf")), plot = p, width = 297, height = 210, units = "mm")