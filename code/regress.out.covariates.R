rm(list = ls())

dir_data = "/project2/xuanyao/llw/breastcancerTCGA/txt_file"
setwd(dir_data)

extract_residual <- function(y, x){
  return(lm(y ~ x)$residuals)
}


# import expression and covariates data
load("./tumor.rdata")
cov_exp = read.table("./tumor.exp_pc.txt",
                     header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
cov_snp = read.table("./tumor.snp_pc.txt",
                     header = FALSE, row.names = NULL, stringsAsFactors = FALSE)

colnames(ex) = gnames


# regress out covariates
ex_cov_regressed = apply(ex, 2, function(y) extract_residual(y, as.matrix(cbind(cov_exp, cov_snp))))
saveRDS(ex_cov_regressed, "./ex_var_regressed.rds")

print("Done!")
