rm(list = ls())

options(stringsAsFactors = FALSE)

# dataDGN: dir_data =
# dataDGN: gene.name = colnames(datExpr)

dir_data = "/project2/xuanyao/llw/breastcancerTCGA/txt_file"
setwd(dir_data)

### Data preparation ###
datExpr = readRDS("./ex_var_regressed.rds")
gene.name = sapply(colnames(datExpr), function(x) strsplit(x, "\\|")[[1]][1])

# remove genes on chromosome X & Y
gene.meta = read.table("./tumor.gene_pos.txt",
                       header = TRUE,
                       stringsAsFactors = F)
rownames(gene.meta) = gene.meta$gene
ind_chrX = gene.meta[colnames(datExpr), "chr"] == "chrX" | gene.meta[colnames(datExpr), "chr"] == "chrY"

# remove pseudogenes
pseudogenes = read.table("./pseudogenes.txt", sep ="\t", stringsAsFactors = F)
ind_pseudo = gene.name %in% unique(pseudogenes$V1)

# remove low mappability genes
mappability = read.table("./mappability.txt", sep ="\t", stringsAsFactors = F, header = T)
ind_map = gene.name %in% unique(mappability$gene_name)

# remove cross-mappable genes
cross_mappability = readRDS("./cross_mappable_genes.rds")
cross_map_genes = unique(names(cross_mappability[!is.na(cross_mappability)]))
ind_cross_map = gene.name %in% cross_map_genes

ind_remove = ind_chrX | ind_pseudo | ind_map | ind_cross_map

sum(ind_chrX)
sum(ind_pseudo)
sum(ind_map)
sum(ind_cross_map)
sum(ind_remove)


write.table(as.matrix(data.frame("gene" = colnames(datExpr),
                                 "gene.name" = gene.name,
                                 "ind_remove" = ind_remove,
                                 "ind_chrX" = ind_chrX,
                                 "pseudo" = ind_pseudo,
                                 "lowmap" = ind_map,
                                 "crossmap" = ind_cross_map)),
            "genes_rm_info.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
