#!/usr/bin/env Rscript
'
Usage:
  Plot_GOI_levels.R --count_mat=<count_mat> --metadata=<metadata> --cnvan_key=<cnvan_key> --GOI=<GOI> --model=<model> --prefix=<prefix> [options]

Mandatory arguments:
  --count_mat=<count_mat>   Count matrix
  --metadata=<metadata>     Experiment metadata
  --cnvan_key=<cnvan_key>   Gene conversion key
  --GOI=<GOI>               Genes of interest in txt format
  --model=<model>           Model design for DGEA
  --prefix=<prefix>         Prefix for output filenames

Optional arguments:
  --resDir=<resDir>           Output directory [default: ./]
' -> doc

library(conflicted)
library(docopt)
arguments <- docopt(doc, version = "0.1")
print(arguments)

library(tidyverse)
library(conflicted)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(cowplot)

conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Load parameters
count_mat <- read_tsv(arguments$count_mat)
metadata <- read_tsv(arguments$metadata)
cnvan_key <- read_tsv(arguments$cnvan_key)
GOI <- read_lines(arguments$GOI)
model <- as.formula(arguments$model)
prefix <- arguments$prefix
resDir <- arguments$resDir

# Data wrangling
count_mat <- count_mat |> column_to_rownames(var = "ensembl")
metadata <- metadata |>
  column_to_rownames(var = "sample") |>
  mutate_all(factor) |>
  mutate(treatment = fct_relevel(treatment, metadata |> pull(treatment) |> unique() )) |>
  mutate(label = fct_relevel(label, metadata |> pull(label) |> unique() ))

dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = metadata,
                              design = model)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Remove GOIs with no expression
GOI <- GOI[GOI %in% cnvan_key[cnvan_key$ensembl %in% rownames(dds), ]$symbol]
GOI <- GOI[order(GOI)]

# Set ggplot basic theme
theme_set(
  theme(panel.grid = element_blank(),
        text = element_text(size = 26),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.line = element_line(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_line(linetype = "dotted"))
        )

plot_GOI <- function (dds, GOI) {
  
gene <- cnvan_key[cnvan_key$symbol == GOI, ] |> pull("ensembl")
data <- plotCounts(dds, gene = gene,
                   intgroup = names(metadata),
                   returnData = TRUE)

ggplot(data, aes(x = patient, y = count)) +
  scale_y_log10() +
  geom_point(aes(fill = treatment),
             shape = 21, position = position_jitter(width = 0.1, height = 0)) +
  scale_fill_manual(values = brewer.pal(4, "RdYlBu")) + # Colour palette "RdYlBu" n = 11! is max number of levels in treatment factor!
  ggtitle(GOI) +
  theme_bw()
}

l <- lapply(seq_along(GOI), function(i) plot_GOI(dds = dds, GOI = GOI[i]) + theme(legend.position = "none"))

# Custom function to limit number of plots on single A4 page to 24 GOIs. If number > 24 plots are pushed to next page
plot_dim <- function(my_list) {
  
  n <- length(my_list)
  k <- 24 # max number of genes per page + legend = 5*5 = 25
  sub_l <- split(my_list, rep(1:ceiling(n/k), each = k)[1:n])
  return(sub_l)
}

l <- plot_dim(l)

lapply(seq_along(l), function(i){
    p <- plot_grid(plotlist = l[[i]],
    get_legend(plot_GOI(dds = dds, GOI = GOI[1])),
    ncol = 5, nrow = 5) # plot arrangement
    
    ggsave(file.path(resDir, paste0(prefix, "_GOI_expression_levels_", i, ".pdf")), plot = p, width = 297, height = 210, units = "mm")
  })