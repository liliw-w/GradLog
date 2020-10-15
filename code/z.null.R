rm(list = ls())

parameter = commandArgs(trailingOnly = T)
module_id = parameter[1]
chr = parameter[2]


library(foreach)
library(doParallel)

# data input
data_dir = "/project2/xuanyao/data/DGN/"
datExpr = readRDS(paste0(data_dir, "txt_file_of_DGN/ex_var_regressed.rds"))
gene_pos = read.table(paste0(data_dir, "txt_file_of_DGN/expression_meta_hg19.txt"),
                      sep = '\t',
                      header = T,
                      stringsAsFactors = F)
gene_result = readRDS(file = paste0(data_dir, "txt_file_of_DGN/DGN_coexpressed_gene_network.rds"))
header = ifelse(chr==1, TRUE, FALSE)
genotype = read.table(paste0(data_dir, "txt_file_of_DGN/chr", chr, ".genotype.matrix.eqtl.txt"),
                      header = header, row.names=NULL,
                      stringsAsFactors = F)
genotype = genotype[!duplicated(genotype[, 1]), ]
rownames(genotype) = genotype[, 1]
genotype = t(genotype[, -1])


# Use only the genes in the module and have postion info in the meta file
gene_in_cluster = data.frame("gene_name"=names(gene_result$moduleLabels)[gene_result$moduleLabels == module_id], stringsAsFactors = F)
gene_w_pos = merge(gene_in_cluster, gene_pos, by.x = 1, by.y = 1)
exp.data = datExpr[, gene_w_pos$gene_name]
print("Input part is done!")

# run parallel on snps
cores=(Sys.getenv("SLURM_NTASKS_PER_NODE"))
cl = makeCluster(as.numeric(cores), outfile="")
registerDoParallel(cl)

M = ncol(genotype); K = ncol(exp.data); N = nrow(exp.data); ind.perm = sample(1:N)
genotype = genotype[ind.perm, ]
chunk.size = M %/% as.numeric(cores)
system.time(
  z.null <- foreach(i=1:cores, .combine='rbind') %dopar%
    {
      res = matrix(nrow = chunk.size, ncol = K,
                   dimnames = list(colnames(genotype)[((i-1)*chunk.size+1):(i*chunk.size)], colnames(exp.data)))

      for(j in 1:chunk.size){
        snp = genotype[, j + (i-1)*chunk.size]
        res[j, ] = apply(exp.data, 2, function(x) summary(lm(x~snp))$coefficients[2, 3] )

        if(j %% 100 == 0){print(j+(i-1)*chunk.size)}
      }
      res
    }
)
stopCluster(cl)
print("Parallel part is done!")

# Run PCO on snps left in the parallel part
Nsnps_left = M %% as.numeric(cores)
Nsnps_done = (chunk.size)*as.numeric(cores)
if(Nsnps_left != 0){
  res = matrix(nrow = Nsnps_left, ncol = K,
               dimnames = list(colnames(genotype)[Nsnps_done+1:(Nsnps_left)], colnames(exp.data)) )

  for(j in 1:Nsnps_left){
    snp = genotype[, j + Nsnps_done]
    res[j, ] = apply(exp.data, 2, function(x) summary(lm(x~snp))$coefficients[2, 3] )

    print(Nsnps_done+j)
  }

  z.null = rbind(z.null, res)
}
print("All part is done!")


saveRDS(z.null, file = paste0("./NULL/z.null/z.null.module", module_id, ".chr", chr, ".rds"))
