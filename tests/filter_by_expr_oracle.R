#!/usr/bin/env Rscript
# edgeR filterByExpr() oracle.
# Usage: filter_by_expr_oracle.R <counts.tsv> [group_csv] [lib_size.tsv]
suppressMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
counts_path <- args[1]
group_csv <- if (length(args) >= 2 && nzchar(args[2])) args[2] else NA
libsize_path <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NA

raw <- read.delim(counts_path, check.names = FALSE)
gene_col <- colnames(raw)[1]
genes <- raw[[1]]
y <- as.matrix(raw[, -1, drop = FALSE])
rownames(y) <- genes

group <- if (!is.na(group_csv)) factor(strsplit(group_csv, ",")[[1]]) else NULL
lib.size <- if (!is.na(libsize_path)) {
  nf <- read.delim(libsize_path, header = FALSE)
  as.numeric(nf[[ncol(nf)]])
} else NULL

keep <- filterByExpr(y, group = group, lib.size = lib.size)

cat(gene_col, "keep", sep = "\t")
cat("\n")
write.table(
  data.frame(gene = genes, keep = ifelse(keep, "TRUE", "FALSE")),
  file = stdout(), sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = FALSE
)
