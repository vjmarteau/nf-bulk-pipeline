#!/usr/bin/env Rscript
'
Usage:
  variancePartition.R --count_mat=<count_mat> --metadata=<metadata> --cnvan_key=<cnvan_key> --model=<model> --prefix=<prefix> [options]

Mandatory arguments:
  --count_mat=<count_mat>   Count matrix
  --metadata=<metadata>     Experiment metadata
  --cnvan_key=<cnvan_key>   Gene conversion key
  --model=<model>           Model design for DGEA
  --prefix=<prefix>         Prefix for output filenames

Optional arguments:
  --resDir=<resDir>         Output directory [default: ./]
' -> doc

library(conflicted)
library(docopt)
arguments <- docopt(doc, version = "0.1")
print(arguments)

library(tidyverse)
library(conflicted)
library(tidyverse)
library(DESeq2)
library(variancePartition)

conflict_prefer("expand", "tidyr")

# Load parameters
count_mat <- read_tsv(arguments$count_mat)
metadata <- read_tsv(arguments$metadata)
cnvan_key <- read_tsv(arguments$cnvan_key)
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
dds <- estimateSizeFactors(dds)

# variancePartition
isexpr <- rowSums(fpm(dds) > 1) >= 0.5 * ncol(dds) # identify genes that pass expression cutoff
quantLog <- assay(rlog(dds))[isexpr, ] + 1
varPart <- fitExtractVarPartModel(quantLog, model, metadata)
vp <- sortCols(varPart)
p <- plotVarPart(vp)

ggsave(file.path(resDir, paste0(prefix, "_variancePartition.pdf")), plot = p, width = 297, height = 210, units = "mm")

df <- as.data.frame(vp) |>
 rownames_to_column(var = "ensembl") |>
 left_join(cnvan_key) |>
 relocate(symbol, .after = ensembl) |>
 as_tibble() |>
 arrange(across(3))

write_tsv(df, file.path(resDir, paste0(prefix, "_variancePartition.tsv")))