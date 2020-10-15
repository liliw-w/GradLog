rm(list = ls())

args = commandArgs(trailing = T)
chr = args[1]

dat = read.table(paste("/project2/xuanyao/data/DGN/txt_file_of_DGN/chr",chr,".raw",sep=""),header=T,check.names=F)
geno = as.matrix(dat[,7:ncol(dat)])
rownames(geno) = as.character(dat$IID)
colnames(geno) = unlist(lapply(strsplit(colnames(dat)[7:ncol(dat)],"_"),"[",1))
write.table(t(geno),
            paste("./project2/xuanyao/data/DGN/txt_file_of_DGN/chr",chr,".genotype.matrix.eqtl.txt",sep = ""),
            sep = "\t", quote = FALSE)