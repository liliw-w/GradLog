rm(list = ls())

library(WGCNA)
source("functions_for_additional_analysis.R")
options(stringsAsFactors = FALSE)

# expression input
datExpr = readRDS("./data/ex_var_regressed.rds")

# remove pseudogenes
pseudogenes = read.table("./result/pseudogenes.txt", sep ="\t", stringsAsFactors = F)
ind_pseudo = colnames(datExpr) %in% unique(pseudogenes$V1)

# remove low mappability genes
mappability = read.table("./result/mappability.txt", sep ="\t", stringsAsFactors = F, header = T)
ind_map = colnames(datExpr) %in% unique(mappability$gene_name)

# remove cross-mappable genes
cross_mappability = readRDS("./result/cross_mappable_genes.rds")
ind_cross_map = unique(names(cross_mappability[!is.na(cross_mappability)]))
ind_cross_map = colnames(datExpr) %in% ind_cross_map

# all removed genes
ind_remove = ind_pseudo | ind_map | ind_cross_map
