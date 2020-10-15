rm(list = ls())

parameter = commandArgs(trailingOnly = T)
module_id = parameter[1]
chr = parameter[2]
p.method = parameter[3] #"truncPCO" #"PCO"

script_dir = "./NULL/script/"
zscore_dir = "./NULL/z.null/"

# There files are modified PCO functions. Only use PC's with eigenvalue \lambda>1
source(paste0(script_dir, "ModifiedPCOMerged.R"))
source(paste0(script_dir, "liu.r"))
source(paste0(script_dir, "liumod.R"))
source(paste0(script_dir, "davies.R"))
source(paste0(script_dir, "SigmaMetaEstimate.R"))
library(foreach)
library(doParallel)
library(mvtnorm)
library(MPAT)


# zscore input
z.null = readRDS(file = paste0(zscore_dir, "z.null.module", module_id, ".chr", chr, ".rds"))


# input
data_dir = "/project2/xuanyao/data/DGN/txt_file_of_DGN/"
gene_all_pos = read.table(paste0(data_dir, "expression_meta_hg19.txt"),
                          sep = '\t',
                          header = T,
                          stringsAsFactors = F)
datExpr = readRDS(paste0(data_dir, "ex_var_regressed.rds"))
gene_module = readRDS(file = paste0(data_dir, "DGN_coexpressed_gene_network.rds"))


# Extract genes within the module that have gene meta info and are not on the chromosome for which we are calculating pvalues.
# This step makes sure the signals are due to trans genetic effects.
gene_cluster = data.frame("gene_name" = names(gene_module$moduleLabels[gene_module$moduleLabels == module_id]), stringsAsFactors = F)
gene_pos = merge(gene_cluster, gene_all_pos, by.x = 1, by.y = 1)
gene_trans = gene_pos[gene_pos$chrom != paste0("chr", chr), "gene_name"]

# inputs for PCO
Sigma = cor(datExpr[, gene_trans])
SigmaO = ModifiedSigmaOEstimate(Sigma,simNum=1000)
z.null_trans = z.null[, gene_trans]; rm(z.null)


# Do parallel across snps
cores=(Sys.getenv("SLURM_NTASKS_PER_NODE"))
cl = makeCluster(as.numeric(cores), outfile="")
registerDoParallel(cl)

Nsnps = nrow(z.null_trans); snps = rownames(z.null_trans)
chunk.size = Nsnps %/% as.numeric(cores)
Nsnps_left = Nsnps %% as.numeric(cores)
Nsnps_done = (chunk.size)*as.numeric(cores)

system.time(
  p.null <- foreach(i=1:cores, .combine='c') %dopar%
    {
      library(MPAT)
      library(mvtnorm)

      res <- rep(NA, chunk.size); names(res) <- snps[((i-1)*chunk.size+1):(i*chunk.size)]
      for(x in 1:chunk.size) {
        res[x] <- ModifiedPCOMerged(Z.vec=z.null_trans[x + (i-1)*chunk.size, ],
                                    Sigma=Sigma, SigmaO=SigmaO,method = "davies")

        if(x %% 100 == 0){print(x + (i-1)*chunk.size)}
      }

      res
    }
)
stopCluster(cl)

print("Parallel part is done!")


# Do PCO on snps that are left in the parallel part.
if(Nsnps_left != 0){
  res_left = rep(NA, Nsnps_left); names(res_left) = snps[1:(Nsnps_left) + Nsnps_done]
  for(k in 1:(Nsnps_left)){
    res_left[k] = ModifiedPCOMerged(Z.vec=z.null_trans[k + Nsnps_done, ],
                                    Sigma=Sigma, SigmaO=SigmaO,method = "davies")

    print(k + Nsnps_done)
  }

  p.null = c(p.null, res_left)
}

print("All part is done!")


saveRDS(pT.null, file = paste0("./NULL/p.null/p.null.module", module_id, ".chr", chr, ".", p.method, ".rds"))
