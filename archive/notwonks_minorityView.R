#### Network analysis of Barkai Data ####
# building networks based on "minority view" (opposite of consensus) --- if there's a connection in at least one experiment, there's a conneciton in the combined network
# Following WGCNA tutorials: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Following Jeremy Miller's tutorial on Meta-analysis with WGCNA: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/
sapply(c("tidyverse", "edgeR", "WGCNA", "matrixStats"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/paralogs/DDivergence/Barkai_data_analysis")
options(stringsAsFactors = FALSE)

#### Reading & Cleaning Data ####
# read in count data (un-normalized) and sample info
load("Cleaned_Barkai_Data.RData")

# normalize counts to cpm (don't want rpkm b/c it's tag-seq data---1 read = 1 3' UTR)
counts <- cpm(counts)

# checking that there aren't any genes with lots of missing values or outliers
gsg <- goodSamplesGenes(counts, verbose = 3)
gsg$allOK

# separate out count matrices for each of the 4 organisms (cer, par, hyc, hyp) and 4 experiments (LowN, LowPi, CC, and HAP4) (also genes need to be columns aparently)
sum(colnames(counts) == sample_info$sample_name)
ncol(counts)
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
# but because we want the same set of genes in the end, we need to eliminate samples that have zero variance in ANY of the 12 networks
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
sum(dropGenes)
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

# removing outlier samples by Euclidean distance clustering
dev.off()
par(cex = 0.2);
par(mar = c(0,4,2,0))
for (i in 1:16) {
  sampleTree <- hclust(dist(counts_all16[[i]]), method = "average")
  plot(sampleTree, main = names(counts_all16)[i])
} # we're looking for samples that are all off on their own---good clustering has two initial branches both leading to a decent total number of samples
# seems both hybrid CC WT 37 samples are outliers as well as 2 samples in LowN_par and WT 28 in Cer_CC
# removing offending samples in the hybrid CC by hand
counts_all16$CC_hyp <- counts_all16$CC_hyp[rownames(counts_all16$CC_hyp) != "WT_37_hyp_CellCycle_rep2",]
counts_all16$CC_hyc <- counts_all16$CC_hyc[rownames(counts_all16$CC_hyc) != "WT_37_hyc_CellCycle_rep2",]
# removing a cluster of ~5 samples in paradoxus LowN using tree cutting
parTree <- hclust(dist(counts_all16$LowN_par), method = "average")
plot(parTree, main = "LowN par")
abline(h = 35000, col = "red");
# Determine cluster under the line
clust <- cutreeStatic(parTree, cutHeight = 35000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples <- (clust==1)
counts_all16$LowN_par <- counts_all16$LowN_par[keepSamples,]
dev.off()

# filtering on top 3000 most connected (highest degree) genes based on average connectivity across 
# we decide on softPower in the parameter fitting section, but it's needed for soft connectivity
softPower <- 4
avgConn <- map_dfr(counts_all16, softConnectivity, type = "unsigned", power = softPower) %>% 
  rowMeans() # I have no idea why it's rowMeans not colMeans... By my calculation map_dfr should get a 16 x nGenes matrix, so colMeans should return the mean connectivity of each gene but it gives a vector of length 16...
hist(avgConn)
keepGenes <- rank(-avgConn) <= 3000
counts_top3000 <- lapply(counts_all16, function(x) {
  return(x[,keepGenes])
})

#### Parameter Fitting ####
# pick power parameter ß for determining soft threshold of edge presence in our network
powers <- c(1:20)
sfts <- lapply(counts_top3000, pickSoftThreshold, powerVector = powers, verbose = 5) # It really just calculates R^2 for different log(k) vs log(p(k)) (aka the classic scale free network) plots
 
# degree distribution plots using different ß-based soft thresholds

# wrapping this WGCNA code in a function for better repeatability
plotSoftPower <- function(sft, sft_name) {
  # plotting results to pick power
  par(mfrow=c(1,2))
  cex1 = 0.9
  # Scale-free topology fit index vs soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], # I think this is just to make sure that the slope is negative? Otherwise you could just use the 2nd column without the 3rd
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence\n", sft_name))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=0.9, col="red") # this line corresponds to using an R^2 cut-off of h apparently
  # Mean connectivity vs soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}
map2(sfts, names(sfts), function(x, y) {
  return(plotSoftPower(x, y))
})

# 6 looks good for all but the HAP4 paradoxus, but nothing looks good for that batch so...
softPower <- 6

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

# can't do all 16 pairwise, that would be 16 choose 2 comparisons,
# but across-species comparisons would confound real evolution and experimental batch effects anyway
# cer
ExprConn_cer <- correlateExpressionAndConnectivity(cts_list = list(counts_top3000$CC_cer,
                                                               counts_top3000$LowN_cer,
                                                               counts_top3000$LowPi_cer,
                                                               counts_top3000$HAP4_cer),
                                               ntwks_titles = c("CC", "LowN", "LowPi", "HAP4"))
grid.arrange(grobs = ExprConn_cer, nrow = 6) # have to export multiplying height by 6 to get a good looking plot

# par
ExprConn_par <- correlateExpressionAndConnectivity(cts_list = list(counts_top3000$CC_par,
                                                                   counts_top3000$LowN_par,
                                                                   counts_top3000$LowPi_par,
                                                                   counts_top3000$HAP4_par),
                                                   ntwks_titles = c("CC", "LowN", "LowPi", "HAP4"))
grid.arrange(grobs = ExprConn_par, nrow = 6)

# hyc
ExprConn_hyc <- correlateExpressionAndConnectivity(cts_list = list(counts_top3000$CC_hyc,
                                                                   counts_top3000$LowN_hyc,
                                                                   counts_top3000$LowPi_hyc,
                                                                   counts_top3000$HAP4_hyc),
                                                   ntwks_titles = c("CC", "LowN", "LowPi", "HAP4"))
grid.arrange(grobs = ExprConn_hyc, nrow = 6)

# hyp
ExprConn_hyp <- correlateExpressionAndConnectivity(cts_list = list(counts_top3000$CC_hyp,
                                                                   counts_top3000$LowN_hyp,
                                                                   counts_top3000$LowPi_hyp,
                                                                   counts_top3000$HAP4_hyp),
                                                   ntwks_titles = c("CC", "LowN", "LowPi", "HAP4"))
grid.arrange(grobs = ExprConn_hyp, nrow = 6)

# Conclusion: Avg expression correlation looks good in all experiments,
# but connectiviy correlation is fairly bad.
# repeated with only top 3000 most connected genes to see if connectivity correlation improves, and it did not
# TODO: look into how big an issue this is. Why wasn't it an issue for my un-normalized counts? I guess becuase the noise from artifacts to do with library size in each sample were shared for all genes so correlations were driven by these artifacts?

#### Building networks ####
# Module selection: Adjacency matrix, topical overlap (TOM), and heirachical clustering

# set soft power (decided in parameter fitting section)
softPower <- 6

# getting adjacency matrix from count matrix (this takes a bit)
adjFromCounts <- function(cts) {
  adj <- abs(cor(cts, use = "p"))^softPower # equivalent to WGCNA::adjacency(x, type = "unsigned", power = softPower) as far as I can tell
  diag(adj) <- 0
  return(adj)
}

# creating adjacency matrix and TOM
adjs <- lapply(counts_top3000, adjFromCounts)
TOMs <- lapply(adjs, TOMsimilarity, TOMType = "unsigned")
dissTOMs <- lapply(TOMs, function(x) {
  return(1-x)
})
geneTrees <- lapply(dissTOMs, function(x) {
  return(hclust(as.dist(x), method = "average"))
})
for (i in 1:16) {
  plot(geneTrees[[i]], xlab="", sub="", main=names(counts_top3000)[i],
       labels = FALSE, hang=0.04)
}

# module selection
minModuleSize <- 20
# building one TOM for each "species" (cer, par, hyc, hyp) based on minority view
# first scaling the 16 individual TOMs so that the 95th percentile of each is equal
scaleP <- 0.95
scaleQuant <- lapply(TOMs, quantile, probs = scaleP, type = 8)
scalePowers <- rep(1, 16)
unscaledTOMs <- TOMs
for (set in 1:12) {
  if (set > 1) {
    scalePowers[set] <- log(scaleQuant[[1]])/log(scaleQuant[[set]])
  }
  TOMs[[set]] <- TOMs[[set]]^scalePowers[set]
}
unscaledDissTOMs <- dissTOMs
dissTOMs <- lapply(TOMs, function(x) {
  return(1-x)
})
geneTrees <- lapply(dissTOMs, function(x) {
  return(hclust(as.dist(x), method = "average"))
})
# constructing minority view TOM
mvTOM <- list(cer = pmax(TOMs$CC_cer, TOMs$LowN_cer, TOMs$LowPi_cer, TOMs$HAP4_cer),
              par = pmax(TOMs$CC_par, TOMs$LowN_par, TOMs$LowPi_par, TOMs$HAP4_par),
              hyc = pmax(TOMs$CC_hyc, TOMs$LowN_hyc, TOMs$LowPi_hyc, TOMs$HAP4_hyc),
              hyp = pmax(TOMs$CC_hyp, TOMs$LowN_hyp, TOMs$LowPi_hyp, TOMs$HAP4_hyp)) # pmax is "parallel" max, takes as many vectors as you give it of the same length and gives the max at each position

mvDissTOM <- lapply(mvTOM, function(x) {
  return(1-x)
})
# construct tree based on minority view TOM
mvTree <- lapply(mvTOM, function(x) {
  return(hclust(as.dist(1-x), method = "average"))
})

# module selection
# getting individual module selections for each of the 12, to be used for module preservation evaluations
# trying out two different tree cutting algorithms, dynamic and hybrid, and 3 values of deepSplit
colorList <- vector(mode = "list", length = 4)
for (i in 1:4) {
  for (ds in 0:2) {
    dynamic <- cutreeDynamic(dendro = mvTree[[i]], 
                              distM = mvDissTOM[[i]],
                              deepSplit = ds, cutHeight = 0.995,
                              minClusterSize = minModuleSize,
                              pamRespectsDendro = FALSE) %>% labels2colors()
      
    chopped_tree <- cutreeHybrid(dendro = mvTree[[i]], 
                           distM = mvDissTOM[[i]],
                           deepSplit = ds, cutHeight = 0.995,
                           minClusterSize = minModuleSize,
                           pamRespectsDendro = FALSE)
    hybrid <- chopped_tree$labels %>% labels2colors()
    colorList[[i]] <- cbind(colorList[[i]], dynamic, hybrid)
  }
}
# compare tree cutting (module selection) methods
for (i in 1:4) {
  plotDendroAndColors(dendro = mvTree[[i]], colors = cbind(unlist(colorList[[i]])),
                      groupLabels = c("dynamic deepSplit = 0", "hybrid deepSplit = 0",
                                      "dynamic deepSplit = 1", "hybrid deepSplit = 1",
                                      "dynamic deepSplit = 2", "hybrid deepSplit = 2"), 
                      main = names(mvTree)[i],
                      dendroLabels = FALSE)
}

# that turquoise module is a little out of control, but otherwise it's looking good!
# not much of a difference between dynamic and hybrid, or really the different deepSplit
# deepSplit=1 and hybrid seems like fair selections for now
mvColors <- map2(mvTree, mvDissTOM, function(x, y) {
  chopped_tree <- cutreeHybrid(dendro = x, 
                               distM = y,
                               deepSplit = 1, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE)
  colors <- chopped_tree$labels %>% labels2colors()
  return(colors)
})
colors_all16 <- map2(geneTrees, dissTOMs, function(x, y) {
  chopped_tree <- cutreeHybrid(dendro = x, 
                               distM = y,
                               deepSplit = 1, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE)
  colors <- chopped_tree$labels %>% labels2colors()
  return(colors)
})

# here's the within-experiment module preservation---more of a data quality check than biologically meaningful
# note we haven't merged close modules or matched colors yet
for (i in 1:4) {
  idx_16 <- 4*(i-1)
  plotDendroAndColors(mvTree[[i]], cbind(mvColors[[i]],colors_all16[[idx_16+1]], colors_all16[[idx_16+2]], colors_all16[[idx_16+3]], colors_all16[[idx_16+4]]),
                      c(names(mvTree)[i], names(colors_all16)[c(idx_16+1, idx_16+2, idx_16+3, idx_16+4)]), main = names(mvTree)[i],
                      dendroLabels = FALSE)
}

# merging close modules
# calculate eigengenes (eigenvector with largest eigenvalue aka 1st PC of expression matrix (woah way back to that thing) for each module)
# Note: colors is not an aesthetic argument, it's literally the only information on module membership

# creating count matrix of all 4 experiments per species (note that I don't believe 
# the sample order is going to matter, as long as it's the same for all of them)
counts_mv <- list("cer" = rbind(counts_top3000$CC_cer,
                                counts_top3000$LowN_cer,
                                counts_top3000$LowPi_cer,
                                counts_top3000$HAP4_cer),
                  "par" = rbind(counts_top3000$CC_par,
                                counts_top3000$LowN_par,
                                counts_top3000$LowPi_par,
                                counts_top3000$HAP4_par),
                  "hyc" = rbind(counts_top3000$CC_hyc,
                                counts_top3000$LowN_hyc,
                                counts_top3000$LowPi_hyc,
                                counts_top3000$HAP4_hyc),
                  "hyp" = rbind(counts_top3000$CC_hyp,
                                counts_top3000$LowN_hyp,
                                counts_top3000$LowPi_hyp,
                                counts_top3000$HAP4_hyp))

MEs <- map2(counts_mv, mvColors, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y)
  MEs <- MEoutput$eigengenes
  return(MEs)
})
MEDiss <- lapply(MEs, function(x) {
  return(1-cor(x))
})
METree <- lapply(MEDiss, function(x) {
  return(hclust(as.dist(x), method = "average"))
})

# plotting (relatively low effort)
# our somewhat arbitrary threshold of 0.75 correlation where we merge clusters
# (everything with a branch point at a lower height than the abline gets merged)
TestThresh <- 0.25
# all 12
par(mfcol = c(2,2))
for (i in 1:4) {
  plot(METree[[i]], main = names(counts_mv)[i], xlab = "", sub = "")
  abline(h=TestThresh, col = "red")
}

# That correlation is super high for all modules, doesn't look like merging is a good idea.
# This will be good to keep in mind as a downside to the minorityView approach

# matching module colors across networks
# matchLabels looks for modules with significantly overlapping gene content and makes them the same color (not as sophisticated as the permutation tests in modulePreservation, it's just a fisher exact test)
# Edit: I've decided I don't trust matchLabels, but here's how you'd do it:
# consensusColorsUnmatched <- consensusColors
# colorsUnmatched_all12 <- colors_all12
# colorsUnmatched_all4 <- colors_all4
# consensusColors <- lapply(consensusColors, function(x) {
#   return(matchLabels(source = x, reference = consensusColors$cer))
# }) # note that using cer as the reference is arbitrary and I can try changing that and see if it affects downstream results -- it shouldn't but hey
# colors_all12 <- lapply(colors_all12, function(x) {
#   return(matchLabels(source = x, reference = consensusColors$cer))
# }) # note that I'm still using the cer consensus as the reference so that the within-experiment module preservation below is also color-matched
# colors_all4 <- lapply(colors_all4, function(x) {
#   return(matchLabels(source = x, reference = consensusColors$cer))
# })

# final network visualizations
for (i in 1:4) {
  idx_16 <- 4*(i-1)
  plotDendroAndColors(mvTree[[i]], cbind(mvColors[[i]],colors_all16[[idx_16+1]], colors_all16[[idx_16+2]], colors_all16[[idx_16+3]], colors_all16[[idx_16+4]]),
                      c(names(mvTree)[i], names(colors_all16)[c(idx_16+1, idx_16+2, idx_16+3, idx_16+4)]), main = names(mvTree)[i],
                      dendroLabels = FALSE)
}
#TOMplots
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

# TOMplots
# minority vote
for (i in 1:4) {
  TOMplotSubset(cts = counts_mv[[i]], dissTOM = mvDissTOM[[i]], 
                moduleColors = mvColors[[i]], plotTitle = names(counts_mv)[i], sampleSize = 400)
} # wow hyp has something weird going on for sure. Which is also odd because it had the
# best range of correlation between module eigengenes
# individual experiments
for (i in 1:16) {
  TOMplotSubset(cts = counts_top3000[[i]], dissTOM = dissTOMs[[i]], 
                moduleColors = colors_all16[[i]], plotTitle = names(counts_top3000)[i], sampleSize = 400)
}

# We're done building our networks! Let's save all the relevant files
save(file = "Networks_minorityVote.RData", counts_mv, counts_top3000, adjs, TOMs, dissTOMs,
     geneTrees, colors_all16, MEs, mvTOM, mvTree, mvColors, mvDissTOM, softPower)


# TODO: I stopped modifying this code to be minorityView-specific at this point
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
