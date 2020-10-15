rm(list = ls())

# import expression and covariates data
load("./data/dgn.rdata")
covariates_tec = read.table("./data/Technical_factors.txt",
                        sep = '\t',
                        header = T, row.names = 1,
                        stringsAsFactors = F)
covariates_bio = read.table("./data/Biological_and_hidden_factors.txt",
                            sep = '\t',
                            header = T, row.names = 1,
                            stringsAsFactors = F)
sample_names = read.table("data/samples.txt",
                            header = T,
                            stringsAsFactors = F)


extract_residual <- function(y, x){
  return(lm(y ~ x)$residuals)
}


rownames(ex) = sample_names$sample; colnames(ex) = gnames
covariates_tec_used = as.matrix(covariates_tec[rownames(ex), ])
covariates_bio_used = as.matrix(covariates_bio[rownames(ex), ])


# regress out technical covariates
ex_tec_regressed = apply(ex, 2, function(y) extract_residual(y, covariates_tec_used))
ex_tec_bio_regressed = apply(ex_tec_regressed, 2, function(e) extract_residual(e, covariates_bio_used))

saveRDS(ex_tec_bio_regressed, "./data/ex_var_regressed.rds")

print("Done!")
