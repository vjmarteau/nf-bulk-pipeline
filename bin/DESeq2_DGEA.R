#!/usr/bin/env Rscript
'
Usage:
  DESeq2_DGEA.R --count_mat=<count_mat> --metadata=<metadata> --cnvan_key=<cnvan_key> --model=<model> --treat_col=<treat_col> --prefix=<prefix> [options]

Mandatory arguments:
  --count_mat=<count_mat>   Count matrix
  --metadata=<metadata>     Experiment metadata
  --cnvan_key=<cnvan_key>   Gene conversion key
  --model=<model>           Model design for DGEA
  --treat_col=<treat_col>   Column containing the treatment info
  --prefix=<prefix>         Prefix for output filenames

Optional arguments:
  --resDir=<resDir>         Output directory [default: ./]
' -> doc

# load required packages
library(docopt)
arguments <- docopt(doc, version = "0.1")
print(arguments)

library(conflicted)
library(tidyverse)
library(DESeq2)
library(IHW)

conflict_prefer("filter", "dplyr")

# Load parameters
count_mat <- read_tsv(arguments$count_mat)
metadata <- read_tsv(arguments$metadata)
cnvan_key <- read_tsv(arguments$cnvan_key)
model <- as.formula(arguments$model)
treat_col <- as.character(arguments$treat_col)
prefix <- arguments$prefix
resDir <- arguments$resDir

# Data wrangling
count_mat <- count_mat |> column_to_rownames(var = "ensembl")
metadata <- metadata |>
  column_to_rownames(var = "sample") |>
  mutate_all(factor) |>
  mutate(treatment = fct_relevel(treatment, metadata |> pull(treatment) |> unique() )) |>
  mutate(label = fct_relevel(label, metadata |> pull(label) |> unique() ))

# List all useful contrasts for DESeq2
contrast <- tidyr::crossing(treatment = treat_col,
                            var1 = metadata |> pull(treat_col),
                            var2 = metadata |> pull(treat_col)) |>
  filter(var1 != var2) |>
  filter(var2 != "Ctrl")
contrast <- asplit(contrast, 1) |> map(as.character)

# Function to run differential gene expression analysis and get all contrast
run_DESeq2 <- function(counts, meta, model, contrast, key) {
  DESeq2 <- lapply(contrast, function(cf) {
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = meta,
                                  design = model)
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    dds <- DESeq(dds)
    dds <- results(dds, filterFun = ihw, contrast = cf)
    dds <- as.data.frame(dds) |>
      rownames_to_column(var = "ensembl") |>
      left_join(key) |>
      relocate(symbol, .after = ensembl) |>
      as_tibble() |>
      arrange(padj)
  })
  
  # Generate comparison labels from contrasts  
  name <- unlist(contrast)
  name <- name[name!= "treatment"]
  name <- split(as.character(name), 1:length(name) %% 2 == 0)
  name <- paste(name[[1]], name[[2]], sep = "_vs_")
  names(DESeq2) <- name
  return(DESeq2)
}

res_list <- run_DESeq2(count_mat, metadata, model, contrast, cnvan_key)

lapply(seq_along(res_list), function(i) write_tsv(res_list[[i]], file.path(resDir, paste0(prefix, "_de_res_", names(res_list)[i], ".tsv"))))