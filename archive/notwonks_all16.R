#### Network analysis of Barkai Data in all 4 experiments x 4 "organisms" (cer, par, hyc, hyp) ####
# Following WGCNA tutorials: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Following Jeremy Miller's tutorial on Meta-analysis with WGCNA: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/
sapply(c("tidyverse", "edgeR", "WGCNA", "matrixStats"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/paralogs/DDivergence/Barkai_data_analysis")
options(stringsAsFactors = FALSE)

#### Reading & Cleaning Data ####
# read in count data (un-normalized in order to properly interface with DESeq2, but we'll normalize in next step) and sample info 
load("Cleaned_Barkai_Data.RData")

# normalize counts to cpm (don't want rpkm b/c it's tag-seq data---1 read = 1 3' UTR)
counts<- cpm(counts)

# checking that there aren't any genes with lots of missing values or outliers
gsg <- goodSamplesGenes(counts, verbose = 3)
gsg$allOK # good as a whole set, but once we break into 16 sets, we'll have to check again

# separate out count matrices for each of the 4 organisms (cer, par, hyc, hyp) and 4 experiments (CC, LowN, LowPi, and HAP4) (also genes need to be columns aparently)
cat("percent of count data columns in the same order as sample_info rows (should be 1): ", sum(colnames(counts) == sample_info$sample_name)/ncol(counts))
# 16 organisms x experiments (oh my)
counts_CC_cer <- t(counts[, sample_info$experiment == "CC" & sample_info$organism == "cer"])
counts_CC_par <- t(counts[, sample_info$experiment == "CC" & sample_info$organism == "par"])
counts_CC_hyc <- t(counts[, sample_info$experiment == "CC" & sample_info$organism == "hyc"])
counts_CC_hyp <- t(counts[, sample_info$experiment == "CC" & sample_info$organism == "hyp"])
counts_LowN_cer <- t(counts[, sample_info$experiment == "LowN" & sample_info$organism == "cer"])
counts_LowN_par <- t(counts[, sample_info$experiment == "LowN" & sample_info$organism == "par"])
counts_LowN_hyc <- t(counts[, sample_info$experiment == "LowN" & sample_info$organism == "hyc"])
counts_LowN_hyp <- t(counts[, sample_info$experiment == "LowN" & sample_info$organism == "hyp"])
counts_LowPi_cer <- t(counts[, sample_info$experiment == "LowPi" & sample_info$organism == "cer"])
counts_LowPi_par <- t(counts[, sample_info$experiment == "LowPi" & sample_info$organism == "par"])
counts_LowPi_hyc <- t(counts[, sample_info$experiment == "LowPi" & sample_info$organism == "hyc"])
counts_LowPi_hyp <- t(counts[, sample_info$experiment == "LowPi" & sample_info$organism == "hyp"])
counts_HAP4_cer <- t(counts[, sample_info$experiment == "HAP4" & sample_info$organism == "cer"])
counts_HAP4_par <- t(counts[, sample_info$experiment == "HAP4" & sample_info$organism == "par"])
counts_HAP4_hyc <- t(counts[, sample_info$experiment == "HAP4" & sample_info$organism == "hyc"])
counts_HAP4_hyp <- t(counts[, sample_info$experiment == "HAP4" & sample_info$organism == "hyp"])

# Oh great, for downstream applications we need to make sure that none of the genes have zero variance
# but because we want the same set of genes in the end, we need to eliminate samples that have zero variance in ANY of the 16 networks
counts_all16 <- list(CC_cer = counts_CC_cer,
                     LowN_cer = counts_LowN_cer,
                     LowPi_cer = counts_LowPi_cer,
                     HAP4_cer = counts_HAP4_cer,
                     CC_par = counts_CC_par,
                     LowN_par = counts_LowN_par,
                     LowPi_par = counts_LowPi_par,
                     HAP4_par = counts_HAP4_par,
                     CC_hyc = counts_CC_hyc,
                     LowN_hyc = counts_LowN_hyc,
                     LowPi_hyc = counts_LowPi_hyc,
                     HAP4_hyc = counts_HAP4_hyc,
                     CC_hyp = counts_CC_hyp,
                     LowN_hyp = counts_LowN_hyp,
                     LowPi_hyp = counts_LowPi_hyp,
                     HAP4_hyp = counts_HAP4_hyp)
geneVar <- sapply(counts_all16, colVars)
geneVarIsZero <- geneVar == 0
dropGenes <- apply(geneVarIsZero, 1, any)
counts_all16 <- lapply(counts_all16, function(x) {
  x <- x[,!dropGenes]
  return(x)
})

# to avoid confusion, removing the 12 count matrices not in a list b/c we didn't prune their gene set for low var genes
rm(counts_LowN_cer, counts_CC_cer, counts_LowPi_cer, 
   counts_LowN_hyc, counts_LowPi_hyc, counts_CC_hyc,
   counts_LowN_hyp, counts_LowPi_hyp, counts_CC_hyp,
   counts_CC_par, counts_LowN_par, counts_LowPi_par,
   counts_HAP4_cer, counts_HAP4_hyc, counts_HAP4_hyp, counts_HAP4_par) # this order is chaos

# TODO: here's where I left off adapting the notwonks code to be specific to individual experiments' networks

#### Parameter Fitting ####
# pick power parameter ß for determining soft threshold of edge presence in our network
# I'm just doing this per organism for now, as it's just about picking a consistent parameter,
# but I should keep in mind that a more accurate version would probably be per experiment and organism (so 4^2 = 16 total categories)
powers <- c(1:20)
sfts <- lapply(counts_all4, pickSoftThreshold, powerVector = powers, verbose = 5)

# just to compare with a couple organism x experiment networks
sft_LowN_cer <- pickSoftThreshold(counts_LowN_cer, powerVector = powers, verbose = 5)
sft_LowPi_par <- pickSoftThreshold(counts_LowPi_par, powerVector = powers, verbose = 5)
sft_CC_hyp <- pickSoftThreshold(counts_CC_hyp, powerVector = powers, verbose = 5)
sft_LowN_cer_WT <- pickSoftThreshold(counts_LowN_cer_WT, powerVector = powers, verbose = 5)
sft_LowN_hyp <- pickSoftThreshold(counts_LowN_hyp, powerVector = powers, verbose = 5)
sft_LowN_hyp_WT <- pickSoftThreshold(counts_LowN_hyp_WT, powerVector = powers, verbose = 5)
sft_HAP4_cer <- pickSoftThreshold(counts_HAP4_cer, powerVector = powers, verbose = 5)

# Future Anna: It really just calculated R^2 for different log(k) vs log(p(k)) 
# degree distribution plots using different ß-based soft thresholds

# wrapping this WGCNA code in a function for better repeatability
plotSoftPower <- function(sft) {
  # plotting results to pick power
  par(mfrow = c(1,2))
  cex1 = 0.9
  # Scale-free topology fit index vs soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], # I think this is just to make sure that the slope is negative? Otherwise you could just use the 2nd column without the 3rd
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=0.9, col="red") # this line corresponds to using an R^2 cut-off of h apparently
  # Mean connectivity vs soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

par(mfrow = c(2,2))
lapply(sfts, plotSoftPower)
plotSoftPower(sft_cer)
plotSoftPower(sft_par)
plotSoftPower(sft_hyc)
plotSoftPower(sft_hyp)
plotSoftPower(sft_LowPi_par)
plotSoftPower(sft_CC_hyp)
plotSoftPower(sft_HAP4_cer) # HAP4 is definitely behaving weird, far from scale-free, might consider scrapping
plotSoftPower(sft_HAP4)
# hmm WT versus gene deletion doesn't actually affect much about parameter fit
plotSoftPower(sft_LowN_cer)
plotSoftPower(sft_LowN_cer_WT)
plotSoftPower(sft_LowN_hyp)
plotSoftPower(sft_LowN_hyp_WT)

# we'll choose ß = 10 for ours, as that's above an R^2 of 0.9 in the hybrids, and in the parents it's at least plateaued at that
softPower <- 10

# TODO: possibly repeat these with only the most highly connected genes (however you find those) to see if the data are more consistent between experiments
# but I did try this ^ a bit for building gene trees and it seemed to make things worse not better, so I'll not do this for now

##### Data exploration ####
# see how well average expression level of each gene and overall connectivity of the inferred networks correlate between data sets
# comparisons: each of the 4 experiments within one organism (ex: HAP4 vs LowPi in hyc), (4 choose 2) * 4 = 24 total comparisons
correlateExpressionAndConnectivity <- function(cts_list, ntwks_titles) {
  n_ntwks <- length(cts_list)
  n_comparisons <- choose(n_ntwks, 2)
  plotList <- vector(mode = "list", length = n_comparisons*2)
  plotCounter <- 1
  
  for (i in 1:(n_ntwks-1)) {
    for (j in (i+1):n_ntwks) {
      rankExpr1 <- rank(rowMeans(t(cts_list[[i]])))
      rankExpr2 <- rank(rowMeans(t(cts_list[[j]])))
      rankConn1 <- rank(softConnectivity(cts_list[[i]], type="signed", power=softPower))
      rankConn2 <- rank(softConnectivity(cts_list[[j]], type="signed", power=softPower))
      plotTibble <- tibble("rankExpr1" = rankExpr1,
                           "rankExpr2" = rankExpr2,
                           "rankConn1" = rankConn1,
                           "rankConn2" = rankConn2)
      exprCorr <- cor(rankExpr1, rankExpr2, method = "pearson")
      exprPvalue <- cor.test(rankExpr1, rankExpr2, method = "pearson")$p.value
      connCorr <- cor(rankConn1, rankConn2, method = "pearson")
      connPvalue <- cor.test(rankConn1, rankConn2, method = "pearson")$p.value
      p_expr <- ggplot(data = plotTibble, aes(x = rankExpr1, y = rankExpr2)) + 
        xlab(paste("Ranked Expression", ntwks_titles[[i]])) + 
        ylab(paste("Ranked Expression", ntwks_titles[[j]])) + 
        ggtitle(paste("Correlation = ", exprCorr, "\n", "p-value = ", exprPvalue)) +
        geom_point()
      p_conn <- ggplot(data = plotTibble, aes(x = rankConn1, y = rankConn2)) + 
        xlab(paste("Ranked Connectivity", ntwks_titles[[i]])) + 
        ylab(paste("Ranked Conenctivity", ntwks_titles[[j]])) + 
        ggtitle(paste("Correlation = ", connCorr, "\n", "p-value = ", connPvalue)) +
        geom_point()
      plotList[[plotCounter]] <- p_expr
      plotCounter <- plotCounter + 1
      plotList[[plotCounter]] <- p_conn
      plotCounter <- plotCounter + 1
    }
  }
  return(plotList)
}

# plotting the correlations for average expression and connectivity ("normalized" by rank, just seeing if the genes are in the same order in terms of their avg expression and connectivity (degree))
network_titles <- list("CC", "LowN", "LowPi", "HAP4") # this is our order now, so don't screw up

# cer (not just WT LowN)
counts_list_cer <- list(counts_CC_cer, counts_LowN_cer, counts_LowPi_cer, counts_HAP4_cer)
ExprConn_cer <- correlateExpressionAndConnectivity(counts_list_cer, network_titles) 
grid.arrange(grobs = ExprConn_cer, nrow = 6)

# par (not just WT LowN)
counts_list_par <- list(counts_CC_par, counts_LowN_par, counts_LowPi_par, counts_HAP4_par)
ExprConn_par <- correlateExpressionAndConnectivity(counts_list_par, network_titles) 
grid.arrange(grobs = ExprConn_par, nrow = 6)

# hybrid cer (not just WT LowN)
counts_list_hyc <- list(counts_CC_hyc, counts_LowN_hyc, counts_LowPi_hyc, counts_HAP4_hyc)
ExprConn_hyc <- correlateExpressionAndConnectivity(counts_list_hyc, network_titles) 
grid.arrange(grobs = ExprConn_hyc, nrow = 6)

# hybrid par (not just WT LowN)
counts_list_hyp <- list(counts_CC_hyp, counts_LowN_hyp, counts_LowPi_hyp, counts_HAP4_hyp)
ExprConn_hyp <- correlateExpressionAndConnectivity(counts_list_hyp, network_titles) 
grid.arrange(grobs = ExprConn_hyp, nrow = 6)

# Conclusion: Avg expression correlation looks good in all experiments,
# but connectiviy correlation is atrocious for HAP4. Further justification to drop HAP4 samples

#### Building networks ####
# might want to consider filtering for top 3000 or so highly connected genes if the data look too noisy
# something like this:
# keepGenesExpr_CC_cer <- rank(-softConnectivity(counts_CC_cer, type = "signed", power = softPower)) <= 3000 # rank gives smallest value the rankof #1, so to specify exactly 3000 genes we need to do negative connectivity
# keepGenesExpr <- (keepGenesExpr_1 & keepGenesExpr_2 & ... & keepGenesExpr_16)
# countsTop3000_CC_cer <- counts_CC_cer[,keepGenesExpr_CC_cer]
# I tried it but didn't seem to do much except remove a few outliers

# Module selection: Adjacency matrix, topical overlap (TOM), and heirachical clustering
# first on "control" dataset then assessing how well the modules from that align with modules
# from other, comparable datasets (so in my case, other experiments in the same species)

# set soft power (decided in parameter fitting section)
softPower <- 10

# getting adjacency matrix from count matrix (this takes a bit)
adjFromCounts <- function(cts) {
  adj <- abs(cor(cts, use = "p"))^softPower # equivalent to WGCNA::adjacency(x, type = "unsigned", power = softPower) as far as I can tell
  diag(adj) <- 0
  return(adj)
}

# all 12 at once! Takes a few minutes (note we build the gene trees below, after scaling the TOMs to be comparable)
adj_all12 <- lapply(counts_all12, adjFromCounts)
TOM_all12 <- lapply(adj_all12, TOMsimilarity, TOMType = "unsigned")

# all 4 species. Takes even longer than all 12 cause the networks are bigger
adj_all4 <- lapply(counts_all4, adjFromCounts)
TOM_all4 <- lapply(adj_all4, TOMsimilarity, TOMType = "unsigned")
dissTOM_all4 <- lapply(TOM_all4, function(x) {
  return(1-x)
})
geneTree_all4 <- lapply(dissTOM_all4, function(x) {
  return(hclust(as.dist(x), method = "average"))
})
par(mfrow = c(2,2))
for (i in 1:4) {
  plot(geneTree_all4[[i]], xlab="", sub="", main=names(counts_all4)[i],
       labels = FALSE, hang=0.04)
}

# just to look at hap4 samples for further justification that they're horrible
counts_HAP4 <- list(HAP4_cer = counts_HAP4_cer, 
                    HAP4_par = counts_HAP4_par, 
                    HAP4_hyc = counts_HAP4_hyc, 
                    HAP4_hyp = counts_HAP4_hyp)
dissTOM_HAP4 <- lapply(counts_HAP4, dissTOMFromCounts)
geneTree_HAP4 <- lapply(dissTOM_HAP4, function(x) {
  return(hclust(as.dist(x), method = "average"))
})
par(mfrow = c(2,2))
for (i in 1:4) {
  plot(geneTree_HAP4[[i]], xlab="", sub="", main=names(counts_HAP4)[i],
       labels = FALSE, hang=0.04) # HAP4 actually has better module separation than the other datasets, at least in the parents, but it's so different that I wouldn't want to include it in a consensus network with the others -- true biological coexpression should be robust to environment
}

# module selection
# testing to see if a consensus TOM or just conglomerate network leads to better module preservation in each of the 3 individual experiments datasets in each species (excluding HAP4)
minModuleSize <- 30

# consensus TOM
# scale the TOMs so that they're more comparable (95th percentile is equal in each)
scaleP <- 0.95
scaleQuant <- lapply(TOM_all12, quantile, probs = scaleP, type = 8)
scalePowers <- rep(1, 12)
unscaledTOM_all12 <- TOM_all12
for (set in 1:12) {
  if (set > 1) {
    scalePowers[set] <- log(scaleQuant[[1]])/log(scaleQuant[[set]])
  }
  TOM_all12[[set]] <- TOM_all12[[set]]^scalePowers[set]
}
unscaledDissTOM_all12 <- dissTOM_all12
dissTOM_all12 <- lapply(TOM_all12, function(x) {
  return(1-x)
})
geneTree_all12 <- lapply(dissTOM_all12, function(x) {
  return(hclust(as.dist(x), method = "average"))
})
# constructing consensus TOM
consensusTOM <- list(cer = pmin(TOM_all12$CC_cer, TOM_all12$LowN_cer, TOM_all12$LowPi_cer),
                     par = pmin(TOM_all12$CC_par, TOM_all12$LowN_par, TOM_all12$LowPi_par),
                     hyc = pmin(TOM_all12$CC_hyc, TOM_all12$LowN_hyc, TOM_all12$LowPi_hyc),
                     hyp = pmin(TOM_all12$CC_hyp, TOM_all12$LowN_hyp, TOM_all12$LowPi_hyp)) # pmin is "parallel" min, takes as many vectors as you give it of the same length and gives the min at each position

consensusDissTOM <- lapply(consensusTOM, function(x) {
  return(1-x)
})

# construct trees based on consensus TOM
consensusTrees <- lapply(consensusTOM, function(x) {
  return(hclust(as.dist(1-x), method = "average"))
})

# module selection
# getting individual module selections for each of the 12, to be used for module preservation evaluations
# trying out two different tree cutting algorithms (may also want to try out different values of deepSplit)
# all 12
dynamicColors_all12 <- map2(geneTree_all12, dissTOM_all12, function(x, y) {
  labels <- cutreeDynamic(dendro = x, 
                          distM = y,
                          deepSplit = 2, cutHeight = 0.995,
                          minClusterSize = minModuleSize,
                          pamRespectsDendro = FALSE)
  colors <- labels2colors(labels)
  return(colors)
})

hybridColors_all12 <- map2(geneTree_all12, dissTOM_all12, function(x, y) {
  chopped_tree <- cutreeHybrid(dendro = x, 
                               distM = y,
                               deepSplit = 2, 
                               cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamStage = FALSE)
  labels <- chopped_tree$labels
  colors <- labels2colors(labels)
  return(colors)
})

hybridColors_all12_ds1 <- map2(geneTree_all12, dissTOM_all12, function(x, y) {
  chopped_tree <- cutreeHybrid(dendro = x, 
                               distM = y,
                               deepSplit = 1, 
                               cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamStage = FALSE)
  labels <- chopped_tree$labels
  colors <- labels2colors(labels)
  return(colors)
})

hybridColors_all12_ds0 <- map2(geneTree_all12, dissTOM_all12, function(x, y) {
  chopped_tree <- cutreeHybrid(dendro = x, 
                               distM = y,
                               deepSplit = 0, 
                               cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamStage = FALSE)
  labels <- chopped_tree$labels
  colors <- labels2colors(labels)
  return(colors)
})

idx_12 <- sample(c(1:12), 1)
plotDendroAndColors(dendro = geneTree_all12[[idx_12]], colors = cbind(dynamicColors_all12[[idx_12]], hybridColors_all12[[idx_12]], hybridColors_all12_ds1[[idx_12]], hybridColors_all12_ds0[[idx_12]]),
                    groupLabels = c("dynamic", "hybrid deepSplit = 2", "hybrid deepSplit = 1", "hybrid deepSplit = 0"), main = names(geneTree_all12)[idx_12],
                    dendroLabels = FALSE)
# ooh hybrid seems like something you can be much more confident in! But maybe too confident? If the data looks fuzzy, the module definitions should probably also be fuzzy. DeepSplit also seems to be a more important parameter for the hybrid
dev.off()
par(mfcol = c(3,4))
for (i in 1:12) {
  plotDendroAndColors(dendro = geneTree_all12[[i]], colors = hybridColors_all12[[i]],
                      groupLabels = "hybrid deepSplit = 2", main = names(geneTree_all12)[i],
                      dendroLabels = FALSE)
}

# consensus
consensusDynamicColors <- map2(consensusTrees, consensusDissTOM, function(x, y) {
  labels <- cutreeDynamic(dendro = x, 
                          distM = y,
                          deepSplit = 2, cutHeight = 0.995,
                          minClusterSize = minModuleSize,
                          pamRespectsDendro = FALSE)
  colors <- labels2colors(labels)
  return(colors)
})

consensusHybridColors <- map2(consensusTrees, consensusDissTOM, function(x, y) {
  chopped_tree <- cutreeHybrid(dendro = x, 
                               distM = y,
                               deepSplit = 2, 
                               cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamStage = FALSE)
  labels <- chopped_tree$labels
  colors <- labels2colors(labels)
  return(colors)
})

consensusHybridColors_ds1 <- map2(consensusTrees, consensusDissTOM, function(x, y) {
  chopped_tree <- cutreeHybrid(dendro = x, 
                               distM = y,
                               deepSplit = 1, 
                               cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamStage = FALSE)
  labels <- chopped_tree$labels
  colors <- labels2colors(labels)
  return(colors)
})

consensusHybridColors_ds0 <- map2(consensusTrees, consensusDissTOM, function(x, y) {
  chopped_tree <- cutreeHybrid(dendro = x, 
                               distM = y,
                               deepSplit = 0, 
                               cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamStage = FALSE)
  labels <- chopped_tree$labels
  colors <- labels2colors(labels)
  return(colors)
})

for (i in 1:4) {
  plotDendroAndColors(dendro = consensusTrees[[i]], colors = cbind(consensusDynamicColors[[i]], consensusHybridColors[[i]], consensusHybridColors_ds1[[i]], consensusHybridColors_ds0[[i]]),
                      groupLabels = c("dynamic", "hybrid deepSplit = 2", "hybrid deepSplit = 1", "hybrid deepSplit = 0"), main = names(consensusTrees)[i],
                      dendroLabels = FALSE)
} # deepSplit = 2 is nice because it identifies all the modules that dynamic was most confident in without the noise. The downside is that it looks like some pairs of modules are probably supposed to be the same module, but we can merge them based on module eigengenes in the next section

# selecting dynamic deepSplit = 2 as our module definitions, fuzziness at this stage isn't a bad thing as we'll be merging modules
colors_all12 <- dynamicColors_all12
consensusColors <- consensusDynamicColors

# conglomerate
# we're just using the parameters (dynamic and deepSplit=2) we fit on the consensus network because the whole point of conglomerate is to do the same thing as consensus but without correcting for unequal sample numbers in the different experiments and seeing if it makes a difference
colors_all4 <- map2(geneTree_all4, dissTOM_all4, function(x, y) {
  labels <- cutreeDynamic(dendro = x, 
                          distM = y,
                          deepSplit = 2, 
                          cutHeight = 0.995,
                          minClusterSize = minModuleSize,
                          pamRespectsDendro = FALSE)
  colors <- labels2colors(labels)
  return(colors)
})

# here's the within-experiment module preservation---more of a data quality check than biologically meaningful
# note we haven't merged close modules or matched colors yet
for (i in 1:4) {
  idx_12 <- 3*(i-1)
  plotDendroAndColors(consensusTrees[[i]], cbind(colors_all12[[idx_12+1]], colors_all12[[idx_12+2]], colors_all12[[idx_12+3]]),
                      names(colors_all12)[c(idx_12+1, idx_12+2, idx_12+3)], main = names(consensusTrees)[i],
                      dendroLabels = FALSE) # this function unfortunately does not respect my desire to have a 2x2 grid of plots
} # I mean, who's to say anything from that

# merging close modules
# calculate eigengenes (eigenvector with largest eigenvalue aka 1st PC of expression matrix (woah way back to that thing) for each module)
# Note: colors is not an aesthetic argument, it's literally the only information on module membership
# individual experiment networks
MEs_all12 <- map2(counts_all12, colors_all12, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y)
  MEs <- MEoutput$eigengenes
  return(MEs)
})
MEDiss_all12 <- lapply(MEs_all12, function(x) {
  return(1-cor(x))
})
METree_all12 <- lapply(MEDiss_all12, function(x) {
  return(hclust(as.dist(x), method = "average"))
})
# consensus network
consensusMEs <- map2(counts_all4, consensusColors, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y)
  MEs <- MEoutput$eigengenes
  return(MEs)
})
consensusMEDiss <- lapply(consensusMEs, function(x) {
  return(1-cor(x)) # cor generates square matrix of every ME's correlation with ever other one (note that if there's a negative correlation, 1-cor(x) will be >1)
})
consensusMETree <- lapply(consensusMEDiss, function(x) {
  return(hclust(as.dist(x), method = "average"))
})
# conglomerate networks
MEs_all4 <- map2(counts_all4, colors_all4, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y)
  MEs <- MEoutput$eigengenes
  return(MEs)
})
MEDiss_all4 <- lapply(MEs_all4, function(x) {
  return(1-cor(x))
})
METree_all4 <- lapply(MEDiss_all4, function(x) {
  return(hclust(as.dist(x), method = "average"))
})

# plotting (relatively low effort)
# our somewhat arbitrary threshold of 0.75 correlation where we merge clusters
# (everything with a branch point at a lower height than the abline gets merged)
TestThresh <- 0.25
# all 12
par(mfcol = c(3,4))
for (i in 1:12) {
  plot(METree_all12[[i]], main = names(counts_all12)[i], xlab = "", sub = "")
  abline(h=TestThresh, col = "red")
}
# consensus
dev.off()
par(mfrow = c(2,2))
for (i in 1:4) {
  plot(consensusMETree[[i]], main = names(counts_all4)[i], xlab = "", sub = "")
  abline(h=TestThresh, col = "red")
}
# conglomerate
dev.off()
par(mfrow = c(2,2))
for (i in 1:4) {
  plot(METree_all4[[i]], main = names(counts_all4)[i], xlab = "", sub = "")
  abline(h=TestThresh, col = "red")
}

# A threshold of 0.25 looks a little too harsh, especially for consensus
MEDissThresh <- 0.1
dev.off()
par(mfrow = c(2,2))
for (i in 1:4) {
  plot(consensusMETree[[i]], main = names(counts_all4)[i], xlab = "", sub = "")
  abline(h=MEDissThresh, col = "red")
} # better! Now we're only merging modules if their eigengenes have 0.9 correlation

# actually merging
# individual experiments
unmergedColors_all12 <- colors_all12
unmergedMEs_all12 <- MEs_all12
merge_all12 <- map2(counts_all12, unmergedColors_all12, function(x,y) {
  return(mergeCloseModules(x, y, cutHeight = MEDissThresh, verbose = 3))
})
colors_all12 <- lapply(merge_all12, function(x) {
  return(x$colors)
})
MEs_all12 <- lapply(merge_all12, function(x) {
  return(x$newMEs)
})

# consensus
unmergedConsensusMEs <- consensusMEs
unmergedConsensusColors <- consensusColors
consensusMerge <- map2(counts_all4, unmergedConsensusColors, function(x,y) {
  return(mergeCloseModules(x, y, cutHeight = MEDissThresh, verbose = 3))
})
consensusColors <- list(cer = consensusMerge$cer$colors,
                        par = consensusMerge$par$colors,
                        hyc = consensusMerge$hyc$colors,
                        hyp = consensusMerge$hyp$colors)
consensusMEs <- list(cer = consensusMerge$cer$newMEs,
                     par = consensusMerge$par$newMEs,
                     hyc = consensusMerge$hyc$newMEs,
                     hyp = consensusMerge$hyp$newMEs)

# conglomerate
unmergedColors_all4 <- colors_all4
unmergedMEs_all4 <- MEs_all4
merge_all4 <- map2(counts_all4, unmergedColors_all4, function(x,y) {
  return(mergeCloseModules(x, y, cutHeight = MEDissThresh, verbose = 3))
})
colors_all4 <- lapply(merge_all4, function(x) {
  return(x$colors)
})
MEs_all4 <- lapply(merge_all4, function(x) {
  return(x$newMEs)
})

# matching module colors across networks
# matchLabels looks for modules with significantly overlapping gene content and makes them the same color (not as sophisticated as the permutation tests in modulePreservation, it's just a fisher exact test)
consensusColorsUnmatched <- consensusColors
colorsUnmatched_all12 <- colors_all12
colorsUnmatched_all4 <- colors_all4
consensusColors <- lapply(consensusColors, function(x) {
  return(matchLabels(source = x, reference = consensusColors$cer))
}) # note that using cer as the reference is arbitrary and I can try changing that and see if it affects downstream results -- it shouldn't but hey
colors_all12 <- lapply(colors_all12, function(x) {
  return(matchLabels(source = x, reference = consensusColors$cer))
}) # note that I'm still using the cer consensus as the reference so that the within-experiment module preservation below is also color-matched
colors_all4 <- lapply(colors_all4, function(x) {
  return(matchLabels(source = x, reference = consensusColors$cer))
})

# final network visualizations
# individual experiment networks
for (i in 1:12) {
  plotDendroAndColors(dendro = geneTree_all12[[i]], colors = colors_all12[[i]], 
                      dendroLabels = FALSE, main = names(geneTree_all12)[i])
}

# consensus networks
for (i in 1:4) {
  idx_12 <- 3*(i-1)
  plotDendroAndColors(consensusTrees[[i]], cbind(consensusColors[[i]], colors_all12[[idx_12+1]], colors_all12[[idx_12+2]], colors_all12[[idx_12+3]]),
                      c(names(consensusTrees)[i], names(colors_all12)[c(idx_12+1, idx_12+2, idx_12+3)]), main = names(consensusTrees)[i],
                      dendroLabels = FALSE)
}

# conglomerate networks
for (i in 1:4) {
  idx_12 <- 3*(i-1)
  plotDendroAndColors(geneTree_all4[[i]], cbind(colors_all4[[i]], colors_all12[[idx_12+1]], colors_all12[[idx_12+2]], colors_all12[[idx_12+3]]),
                      c(names(geneTree_all4)[i], names(colors_all12)[c(idx_12+1, idx_12+2, idx_12+3)]), main = names(geneTree_all4)[i],
                      dendroLabels = FALSE)
} # ok you know, the conglomerate networks don't actually look that bad
TOMplotSubset <- function(cts, dissTOM, moduleColors, sampleSize = 1000, plotTitle = "TOM plot") {
  nGenes <- ncol(cts)
  randomSample <- sample(nGenes, sampleSize)
  selectTOM <- dissTOM[randomSample, randomSample]
  plotTOM <- selectTOM^7 # I guess 7 is just a power they thought looked pretty
  diag(plotTOM) <- NA
  plotTree <-  hclust(as.dist(plotTOM), method = "average")
  plotColors <- moduleColors[randomSample]
  TOMplot(-plotTOM, plotTree, plotColors, main = plotTitle)
}

# TOMplot unfortunately doesn't work in a loop, so here are all the individual TOM plots
# conglomerate 
# cer
TOMplotSubset(cts = counts_all4[[1]], dissTOM = dissTOM_all4[[1]], 
              moduleColors = colors_all4[[1]], plotTitle = paste(names(counts_all4)[1], "consensus"))
# par
TOMplotSubset(cts = counts_all4[[2]], dissTOM = dissTOM_all4[[2]], 
              moduleColors = colors_all4[[2]], plotTitle = paste(names(counts_all4)[2], "consensus"))
# hyc
TOMplotSubset(cts = counts_all4[[3]], dissTOM = dissTOM_all4[[3]], 
              moduleColors = colors_all4[[3]], plotTitle = paste(names(counts_all4)[3], "consensus"))
# hyp
TOMplotSubset(cts = counts_all4[[4]], dissTOM = dissTOM_all4[[4]], 
              moduleColors = colors_all4[[4]], plotTitle = paste(names(counts_all4)[4], "consensus"))
# consensus 
# cer
TOMplotSubset(cts = counts_all4[[1]], dissTOM = consensusDissTOM[[1]], 
              moduleColors = consensusColors[[1]], plotTitle = paste(names(counts_all4)[1], "conglomerate"))
# par
TOMplotSubset(cts = counts_all4[[2]], dissTOM = consensusDissTOM[[2]], 
              moduleColors = consensusColors[[2]], plotTitle = paste(names(counts_all4)[2], "conglomerate"))
# hyc
TOMplotSubset(cts = counts_all4[[3]], dissTOM = consensusDissTOM[[3]], 
              moduleColors = consensusColors[[3]], plotTitle = paste(names(counts_all4)[3], "conglomerate"))
# hyp
TOMplotSubset(cts = counts_all4[[4]], dissTOM = consensusDissTOM[[4]], 
              moduleColors = consensusColors[[4]], plotTitle = paste(names(counts_all4)[4], "conglomerate"))
# conglomerate looks significantly better in terms of its final module definitions, we don't want to just end up with one single module afterall

# all 12
# Cer CC
TOMplotSubset(cts = counts_all12$CC_cer, dissTOM = dissTOM_all12$CC_cer, 
              moduleColors = colors_all12$CC_cer, plotTitle = paste(names(counts_all12)[1], "individual experiment"))
# Cer LowN
TOMplotSubset(cts = counts_all12$LowN_cer, dissTOM = dissTOM_all12$LowN_cer, 
              moduleColors = colors_all12$LowN_cer, plotTitle = paste(names(counts_all12)[2], "individual experiment"))
# Cer LowPi
TOMplotSubset(cts = counts_all12$LowPi_cer, dissTOM = dissTOM_all12$LowPi_cer, 
              moduleColors = colors_all12$LowPi_cer, plotTitle = paste(names(counts_all12)[3], "individual experiment"))
# Par LowN
TOMplotSubset(cts = counts_all12$LowN_par, dissTOM = dissTOM_all12$LowN_par, 
              moduleColors = colors_all12$LowN_par, plotTitle = paste(names(counts_all12)[5], "individual experiment"))
# Hyc LowN
TOMplotSubset(cts = counts_all12$LowN_hyc, dissTOM = dissTOM_all12$LowN_hyc, 
              moduleColors = colors_all12$LowN_hyc, plotTitle = paste(names(counts_all12)[8], "individual experiment"))
# Hyp LowN
TOMplotSubset(cts = counts_all12$LowN_hyp, dissTOM = dissTOM_all12$LowN_hyp, 
              moduleColors = colors_all12$LowN_hyp, plotTitle = paste(names(counts_all12)[11], "individual experiment"))

# The more I look, the more I think that meaningful network connections are being obscured in the
# consensus and conglomerate networks (also I do not for the life of me understand why
# the teal/brown module fully encompases the blue one even though it is clearly 
# doing something different on each side of the blue module)

# I'm going to proceed for now with only the individual experiment networks in order to 
# assess what info exactly is being obscurred in the larger networks

# We're done building our networks! Let's save all the relevant files
save(file = "Networks_all12.RData", counts_all12, adj_all12, TOM_all12, dissTOM_all12,
     geneTree_all12, colors_all12, MEs_all12)
save(file = "Networks_consensus.RData", counts_all4, adj_all4, consensusTOM, consensusDissTOM,
     consensusTrees, consensusColors, consensusMEs)
save(file = "Networks_conglomerate.RData", counts_all4, adj_all4, TOM_all4, dissTOM_all4,
     geneTree_all4, colors_all4, MEs_all4)

#### Evaluating network quality with module preservation ####
# goal: visualize module preservation in each experiment per species to see which network, consensus or conglomerate, 
# does a better job of preserving the modules present in each individual experiment

# visualizing
# consensus
par(mfrow = c(2,2))
for (i in 1:4) {
  idx_12 <- 3*(i-1)
  plotDendroAndColors(consensusTrees[[i]], cbind(consensusColors[[i]], colors_all12[[idx_12+1]], colors_all12[[idx_12+2]], colors_all12[[idx_12+3]]),
                      c(names(consensusTrees)[i], names(colors_all12)[c(idx_12+1, idx_12+2, idx_12+3)]), main = names(consensusTrees)[i],
                      dendroLabels = FALSE) # this function unfortunately does not respect my desire to have a 2x2 grid of plots
}
# conglomerate
for (i in 1:4) {
  idx_12 <- 3*(i-1)
  plotDendroAndColors(geneTree_all4[[i]], cbind(colors_all4[[i]], colors_all12[[idx_12+1]], colors_all12[[idx_12+2]], colors_all12[[idx_12+3]]),
                      c(names(geneTree_all4)[i], names(colors_all12)[c(idx_12+1, idx_12+2, idx_12+3)]), main = names(geneTree_all4)[i],
                      dendroLabels = FALSE) 
} # hard to say for sure which is better until we quantify

# streamlining some of the annoyance of the WGCNA modulePreservation function (whyyy does it need lists of lists named data?)
modulePreservationBySpecies <- function(cts_list, ref_colors) {
  nNetworks <- length(cts_list)
  cts <- list(ref = list(data = cts_list[[1]]))
  ref_colors <- list(ref = ref_colors)
  for (i in 2:nNetworks) {
    cts[[i]] <- list(data = cts_list[[i]])
    names(cts)[i] <- names(cts_list)[i]
  }
  mp <- modulePreservation(cts, ref_colors, 
                           referenceNetworks = 1, verbose = 3, 
                           networkType = "unsigned", nPermutations = 30, 
                           maxGoldModuleSize = 100, maxModuleSize = 400)
  return(mp)
}
# consensus: quantifying module preservation in the consensus network for each of the 3 experiments
# cer
mp_cer_consensus <- modulePreservationBySpecies(list(ref = counts_all4$cer,
                                                     CC_cer = counts_all12$CC_cer,
                                                     LowN_cer = counts_all12$LowN_cer,
                                                     LowPi_cer = counts_all12$LowPi_cer),
                                                ref = consensusColors$cer)
# par
mp_par_consensus <- modulePreservationBySpecies(list(ref = counts_all4$par,
                                                     CC_par = counts_all12$CC_par,
                                                     LowN_par = counts_all12$LowN_par,
                                                     LowPi_par = counts_all12$LowPi_par),
                                                ref = consensusColors$par)
# hyc
mp_hyc_consensus <- modulePreservationBySpecies(list(ref = counts_all4$hyc,
                                                     CC_hyc = counts_all12$CC_hyc,
                                                     LowN_hyc = counts_all12$LowN_hyc,
                                                     LowPi_hyc = counts_all12$LowPi_hyc),
                                                ref = consensusColors$hyc)
# hyp
mp_hyp_consensus <- modulePreservationBySpecies(list(ref = counts_all4$hyp,
                                                     CC_hyp = counts_all12$CC_hyp,
                                                     LowN_hyp = counts_all12$LowN_hyp,
                                                     LowPi_hyp = counts_all12$LowPi_hyp),
                                                ref = consensusColors$hyp)




# conglomerate: quantifying preservation of modules present in the 3 individual experiments in the conglomerate network of the corresponding species
# cer
mp_cer_conglomerate <- modulePreservationBySpecies(list(ref = counts_all4$cer,
                                                        CC_cer = counts_all12$CC_cer,
                                                        LowN_cer = counts_all12$LowN_cer,
                                                        LowPi_cer = counts_all12$LowPi_cer),
                                                   ref = colors_all4$cer)

# par
mp_par_conglomerate <- modulePreservationBySpecies(list(ref = counts_all4$par,
                                                        CC_par = counts_all12$CC_par,
                                                        LowN_par = counts_all12$LowN_par,
                                                        LowPi_par = counts_all12$LowPi_par),
                                                   ref = colors_all4$par)

# hyc
mp_hyc_conglomerate <- modulePreservationBySpecies(list(ref = counts_all4$hyc,
                                                        CC_hyc = counts_all12$CC_hyc,
                                                        LowN_hyc = counts_all12$LowN_hyc,
                                                        LowPi_hyc = counts_all12$LowPi_hyc),
                                                   ref = colors_all4$hyc)

# hyp
mp_hyp_conglomerate <- modulePreservationBySpecies(list(ref = counts_all4$hyp,
                                                        CC_hyp = counts_all12$CC_hyp,
                                                        LowN_hyp = counts_all12$LowN_hyp,
                                                        LowPi_hyp = counts_all12$LowPi_hyp),
                                                   ref = colors_all4$hyp)

# visualizing preservation by Zscore
# note this function assumes individual experiments were given to modulePreservationBySpecies in this order: CC, LowN, LowPi
plotPreservationZscores <- function(mp) {
  mp_CC <- mp$preservation$Z[[1]][[2]]
  mp_LowN <- mp$preservation$Z[[1]][[3]]
  mp_LowPi <- mp$preservation$Z[[1]][[4]]
  pres_df <- tibble(module = c(rownames(mp_CC), rownames(mp_LowN), rownames(mp_LowPi)),
                    Zscore = c(mp_CC$Zsummary.pres, mp_LowN$Zsummary.pres, mp_LowPi$Zsummary.pres),
                    modSize = c(mp_CC$Zsummary.pres, mp_LowN$Zsummary.pres, mp_LowPi$Zsummary.pres),
                    experiment = c(rep("CC", nrow(mp_CC)), rep("LowN", nrow(mp_LowN)), rep("LowPi", nrow(mp_LowPi))))
  p <- ggplot(data = pres_df, aes(x = module, y = Zscore)) + geom_point(aes(color = experiment, size = modSize))
  return(list(pres_df, p))
}

mps <- list(cer_conglomerate = mp_cer_conglomerate,
            cer_consensus = mp_cer_consensus,
            par_conglomerate = mp_par_conglomerate,
            par_consensus = mp_par_consensus,
            hyc_conglomerate = mp_hyc_conglomerate,
            hyc_consensus = mp_hyc_consensus,
            hyp_conglomerate = mp_hyp_conglomerate,
            hyp_consensus = mp_hyp_consensus)

pres <- lapply(mps, plotPreservationZscores)

comparePreservationZscores <- function(presdf1, presdf2, name1, name2) {
  modules1 <- paste(presdf1$module, name1, sep = "_")
  modules2 <- paste(presdf2$module, name2, sep = "_")
  presdf <- rbind(presdf1, presdf2)
  presdf$module <- c(modules1, modules2)
  p <- ggplot(data = presdf, aes(x = module, y = Zscore)) + 
    geom_point(aes(color = experiment, size = modSize)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
}

# note the module colors being the same between conglomerate and consensus doesn't necessarily mean they contain a lot of the same genes! (I think? But the doccumentation is making me start to think that they do somehow? Or maybe it's just that they tend to)
comparePreservationZscores(pres$cer_conglomerate[[1]],
                           pres$cer_consensus[[1]],
                           "conglomerate",
                           "consensus")
comparePreservationZscores(pres$par_conglomerate[[1]],
                           pres$par_consensus[[1]],
                           "conglomerate",
                           "consensus")
comparePreservationZscores(pres$hyc_conglomerate[[1]],
                           pres$hyc_consensus[[1]],
                           "conglomerate",
                           "consensus")
comparePreservationZscores(pres$hyp_conglomerate[[1]],
                           pres$hyp_consensus[[1]],
                           "conglomerate",
                           "consensus")
# not a resounding signal as to which network construction is better. But it does seem we'll want the one with the lowest variance in module preservation between the individual experiments
avg_var <- map_dfr(pres, function(x) {
  df <- x[[1]]
  df <- df[df$module != "grey",] # grey is genes not belonging to a module, so we don't care about variance in this category
  mods <- names(table(df$module))
  vs <- sapply(mods, function(mod) {
    return(var(df[df$module == mod, "Zscore"]))
  })
  return(list(avg = mean(vs), std = sd(vs)))
})
avg_var$name <- names(pres)
ggplot(data = avg_var) + geom_bar(aes(x = name, y = avg), stat = "identity") + 
  geom_errorbar(aes(x = name, ymin = avg - std, ymax = avg + std, width = 0.1)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Average variance in module preservation across experiments")
# consensus has lower average variance, but there is so much variance in the average variance (oh god) that it's hard to say for sure, might just be an artifact of the network construction method

#### Comparing Networks ####

# visualizing module preservation between species
# note I haven't run this yet so there's probably lots of errors---trying to 
# compare between-species module preservation within each of the 4 experiments (4 comparisons total)
for (i in 1:4) {
  idx_16 <- c(i+4, i+8, i+12)
  plotDendroAndColors(trees_all16[[i]], cbind(colors_all16[[i]], colors_all16[[idx_16]]),
                      names(trees_all16)[c(i, idx_16)], main = names(trees_all16)[i], 
                      dendroLabels=FALSE)
}

# TODO:
# comparing networks:
# quantify between-species module preservation in consensus vs conglomerate networks
# calculate kME (once you figure out what it is) and look for shared hub genes in different networks
# GO enrichment in modules
