rm(list = ls())

library(qvalue)

# estimate the p-value cut-off to control the false discovery rate at a certain level

parameter = commandArgs(trailingOnly = T)
fdr.level = parameter[1]
adjust.method = parameter[2]

# input pvalues for all 18 gene modules and 22 chromosomes
pvalue_pc = NULL
for(gene_cluster_id in 1:18){
  for(chr in 1:22){
    tmp_pc = readRDS(file = paste0("./TruncPCO/GeneCluster", gene_cluster_id, "/pvalue_pco_", chr, ".rds"))
    names(tmp_pc) = paste0("C", gene_cluster_id, ":", names(tmp_pc))
    pvalue_pc = c(pvalue_pc, tmp_pc)

    cat("Gene cluster:", gene_cluster_id, ". Chr:", chr, "\n")
    print(object.size(pvalue_pc))
  }
}

# histogram of all pvalues
png("./FDR/Hist of all pvalues")
hist(pvalue_pc)
dev.off()
print("Histogram of pvalues saved!")


# adjust pvalues
if(adjust.method == "BH"){
  p.adjusted = p.adjust(pvalue_pc, method = "BH")
  print("BH correction of pvalue done!")


  signal_q = p.adjusted[p.adjusted<fdr.level]
  signal_name = names(signal_q)
  signal_summary = cbind(pvalue_pc[signal_name], signal_q)
  colnames(signal_summary) = c("pvalue.PCO", "qvalue.PCO")
  write.table(signal_summary, file = "./FDR/transeQTL.Summary.BH.txt",
              quote = FALSE)
  print("Signals written in the file!")


  saveRDS(p.adjusted, file = "./FDR/qvalueofPC.BH.rds")
  print("qvalue results saved!")

}else if(adjust.method == "qvalue"){
  qobj_pc = qvalue(p = pvalue_pc, fdr.level=fdr.level)
  print("FDR analysis of pc has been done!")


  signal_q = qobj_pc$qvalues[qobj_pc$significant]
  signal_name = names(signal_q)
  signal_summary = cbind(qobj_pc$pvalues[signal_name], signal_q)
  colnames(signal_summary) = c("pvalue.PCO", "qvalue.PCO")
  write.table(signal_summary, file = "./FDR/transeQTL.Summary.q.txt",
              quote = FALSE)
  print("Signals written in the file!")


  saveRDS(qobj_pc, file = "./FDR/qvalueofPC.q.rds")
  print("qvalue results saved!")

}
