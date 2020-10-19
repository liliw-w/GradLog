rm(list = ls())

parameter = commandArgs(trailingOnly = T)
module_id = parameter[1]
chr = parameter[2]
p.method = parameter[3] #"truncPCO" #"PCO"


setwd("/scratch/midway2/liliw1/TCGA/")

# dataset DGN: dir_data = "/project2/xuanyao/data/DGN/txt_file_of_DGN/"
# dataset DGN: gene_module = readRDS(file = paste0(dir_data, "DGN_coexpressed_gene_network.rds"))
# dataset DGN: gene_meta = read.table(paste0(dir_data, "expression_meta_hg19.txt"), ...)


source(paste0("./script/ModifiedPCOMerged.R"))
source(paste0("./script/liu.r"))
source(paste0("./script/liumod.R"))
source(paste0("./script/davies.R"))
source(paste0("./script/SigmaMetaEstimate.R"))
library(foreach)
library(doParallel)
library(mvtnorm)
library(MPAT)


# zscore input
z_mat = readRDS(file = paste0("./z/z.module", module_id, ".chr", chr, ".rds"))


# data input
dir_data = "/project2/xuanyao/llw/breastcancerTCGA/txt_file/"
datExpr = readRDS(paste0(dir_data, "ex_var_regressed.rds"))
gene_meta = read.table(paste0(dir_data, "tumor.gene_pos.txt"),
                      sep = '\t',
                      header = T,
                      stringsAsFactors = F)
gene_module = readRDS(file = paste0(dir_data, "coexp.module.rds"))


# Use only the genes in the module & have postion info in the meta file & not on chromosome chr
gene_in_cluster = data.frame("gene_name" = names(gene_module$moduleLabels[gene_module$moduleLabels == module_id]), stringsAsFactors = F)
gene_w_pos = merge(gene_in_cluster, gene_meta, by.x = 1, by.y = 1)
gene_trans = gene_w_pos[gene_w_pos[,2] != paste0("chr", chr), "gene_name"]


Sigma = cor(datExpr[, gene_trans])
SigmaO = ModifiedSigmaOEstimate(Sigma,simNum=1000)
z_mat_trans = z_mat[, gene_trans]; rm(z_mat)


# Do parallel
cores=(Sys.getenv("SLURM_NTASKS_PER_NODE"))
cl = makeCluster(as.numeric(cores), outfile="")
registerDoParallel(cl)

Nsnps = nrow(z_mat_trans); snps = rownames(z_mat_trans)
chunk.size = Nsnps %/% as.numeric(cores)
Nsnps_left = Nsnps %% as.numeric(cores)
Nsnps_done = (chunk.size)*as.numeric(cores)

system.time(
  res.all <- foreach(i=1:cores, .combine='c') %dopar%
    {
      library(MPAT)
      library(mvtnorm)

      res <- rep(NA, chunk.size); names(res) <- snps[((i-1)*chunk.size+1):(i*chunk.size)]
      for(x in 1:chunk.size) {
        res[x] <- ModifiedPCOMerged(Z.vec=z_mat_trans[x + (i-1)*chunk.size, ],
                                                       Sigma=Sigma, SigmaO=SigmaO,method = "davies")

        if(x %% 100 == 0){print(x + (i-1)*chunk.size)}
      }

      res
    }
)
stopCluster(cl)

print("Parallel part is done!")


if(Nsnps_left != 0){
  res_left = rep(NA, Nsnps_left); names(res_left) = snps[1:(Nsnps_left) + Nsnps_done]
  for(k in 1:(Nsnps_left)){
    res_left[k] = ModifiedPCOMerged(Z.vec=z_mat_trans[k + Nsnps_done, ],
                               Sigma=Sigma, SigmaO=SigmaO,method = "davies")

    print(k + Nsnps_done)
  }

  res.all = c(res.all, res_left)
}

print("All part is done!")


saveRDS(res.all, file = paste0("./p/p.module", module_id, ".chr", chr, ".", p.method, ".rds"))
