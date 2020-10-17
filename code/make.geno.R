rm(list = ls())

args = commandArgs(trailing = T)
chr = args[1]

# dataset DGN: dir_txt = "/project2/xuanyao/data/DGN/txt_file_of_DGN/"
# dataset DGN: colnames(geno) = unlist(lapply(strsplit(colnames(genotype)[7:ncol(genotype)],"_"),"[",1))

dir_txt = "/project2/xuanyao/llw/breastcancerTCGA/txt_file_of_DGN/"
genotype = read.table(paste0(dir_txt, "chr",chr,".raw"),
                 header = T, check.names = F)
snp.meta = read.table(paste0(dir_txt, "chr",chr,"_genotype_meta.txt"),
                      header = TRUE, stringsAsFactors = FALSE)

geno = as.matrix(genotype[,7:ncol(genotype)])
rownames(geno) = as.character(genotype$IID)
colnames(geno) = snp.meta$id


write.table(t(geno),
            paste0(dir_txt, "chr", chr,".genotype.matrix.eqtl.txt"),
            sep = "\t", quote = FALSE)

cat(paste0("chr", chr, " is done!"))

# check the samples in genotype data is same with that in tumor sample data.
# check if ".genotype.matrix.eqtl.txt" has sample name as header for chr1-22.
# extract sample names
