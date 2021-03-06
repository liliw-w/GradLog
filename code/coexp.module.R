rm(list = ls())

library(WGCNA)

dir_data = "/project2/xuanyao/llw/breastcancerTCGA/txt_file"
setwd(dir_data)

# Data preparation
datExpr = readRDS("./ex_var_regressed.rds")
rm_info = read.table("./genes_rm_info.txt",
                     header = TRUE, row.names = NULL,
                     stringsAsFactors = FALSE)
datExpr = datExpr[, !rm_info$ind_remove]


# Run WGCNA
### Parameter specification ###
minModuleSize = 30
MEDissThres = 0.15
if_plot_adjacency_mat_parameter_selection = F
if_plot_only_tree = F
if_plot_color_and_tree = F
if_plot_eigengene_heatmap_tree = F
if_heatmap_of_network = T


### Step1: network construction ###
# determine the paramter in adjacency function: pickSoftThreshold() OR pickHardThreshold()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower = sft$powerEstimate

# network construction
adjacency = adjacency(datExpr, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM


### Step2: module detection ###
# tree construction using hierarchical clustering based on TOM
geneTree = hclust(as.dist(dissTOM), method = "average")

# branch cutting using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)

# eigene genes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

moduleLabels = match(mergedColors, c("grey", standardColors(100)))-1
names(moduleLabels) = colnames(datExpr)

# Save results
result = list(moduleColors = mergedColors,
              moduleLabels = moduleLabels,
              MEs = mergedMEs,
              old_moduleColors = dynamicColors,
              old_moduleLabels = dynamicMods,
              old_MEs = MEs,
              geneTree = geneTree,
              ind_remove = rm_info$ind_remove)
saveRDS(result, file = "./coexp.module.rds")
