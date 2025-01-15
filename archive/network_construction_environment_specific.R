#### Environment-specific network construction ####
# following WGCNA tutorials: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Following Jeremy Miller's tutorial on Meta-analysis with WGCNA: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/
sapply(c("tidyverse", "edgeR", "WGCNA", "matrixStats"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Barkai_data_analysis")
options(stringsAsFactors = FALSE)

#### Reading & Cleaning Data ####
# read in count data (un-normalized) and sample info
load("data_files/Cleaned_Barkai_Data_Env_Specific.RData")

# edgeR's calculate norm factors might be necessary to convert raw library sizes to effective library sizes (skipping for now b/c as of 4/3/23 I don't remember what this was about)
# normFacts <- lapply(es, calcNormFactors) # TODO: applying it to all 50 lists automatically prevents it from using info across experimental conditions or species---make sure calcNormFactors is only using the DGEList groupings to figure out what blocks to apply the function to together---it's not comparing groups to groups to derive the norm factors

# Full names for plotting
PrettyNames <- character(0)
for (o in c("S. cerevisiae", "S. paradoxus", "F1 hybrid")) {
  for (e in c("TF Deletion", "Low Nitrogen", "TF Deletion x Low Nitrogen", "Cell Cycle", "Growth Curve", "Low Phosphorus")) {
    PrettyNames <- c(PrettyNames, paste(e, "\n", o))
  }
}
for (pn in PrettyNames) {
  hasReplicates <- c("TF Deletion", "Low Nitrogen", "TF Deletion x Low Nitrogen", "Cell Cycle") %>% 
    lapply(FUN = grepl, x = pn) %>% unlist() %>% any()
  if (hasReplicates) {
    PrettyNames <- c(PrettyNames, paste(pn, "\nReplicate 1"), paste(pn, "\nReplicate 2"))
  }
}
print(tibble(pretty = PrettyNames, standard = names(es)), n = 42) # making sure that they all line up

# normalize counts to cpm (don't want rpkm b/c it's tag-seq data---1 read = 1 3' UTR) (And also b/c I'm comparing the same set of genes between species ya dingus)
es <- lapply(es, cpm)

# for some unhinged reason, WGCNA has rows be samples and columns be genes
es <- lapply(es, t)

# checking that there aren't any genes with lots of missing values or outliers
gsg <- lapply(es, goodSamplesGenes, verbose = 3)
lapply(gsg, function(x) {return(x$allOK)}) 
# we'll check this again after our data cleaning. Manual inspection indicates that it's flagging 0-variance genes, not outlier samples or missing values:
sum(gsg$TFdel_cer_rep1$goodGenes)
sum(gsg$TFdel_cer_rep1$goodSamples)
dim(es$TFdel_cer_rep1)

# Removing zero-variance genes
# Because we want the same set of genes in the end, we need to eliminate samples that have zero variance in ANY of the 12 networks
geneVar <- sapply(es, colVars)
geneVarIsZero <- geneVar == 0
dropGenes <- apply(geneVarIsZero, 1, any)
sum(dropGenes)
es <- lapply(es, function(x) {
  x <- x[,!dropGenes]
  return(x)
})

# checking goodSamplesGenes again
gsg <- lapply(es, goodSamplesGenes, verbose = 3)
lapply(gsg, function(x) {return(x$allOK)}) %>% unlist() %>% all() # nice

# removing outlier samples by Euclidean distance clustering
par(cex = 0.35);
par(mar = c(0,4,2,0))
# we're looking for samples that are all off on their own---good clustering has two initial branches both leading to a decent total number of samples
dummy <- lapply(es, function(x) {
  sampleTree <- hclust(dist(x), method = "average")
  plot_name_idx <- lapply(es, function(ntwk) {
    return(identical(ntwk, x))
  }) %>% unlist %>% which()
  plot(sampleTree, main = names(es)[plot_name_idx])
}) # get ready for 42 trees
rm(dummy)
# outliers:
# WT 37 CC hyp/hyc/hyb rep2 (rep2 refers to its sample name, it's not necessarily in rep2 of each dataset)
# removing above samples by hand:
to_remove <- si$CC_hyb_rep1 %>% filter(grepl(37, sample_name)) %>% select(sample_name) %>% as.character()
to_remove %in% si$CC_hyb$sample_name
to_remove %in% si$CC_hyb_rep1$sample_name
to_remove %in% si$CC_hyb_rep2$sample_name # in this case it's in rep1, but edit if it's in rep2 for a different random subsetting round
# remove from full dataset:
si$CC_hyb <- filter(si$CC_hyb, sample_name != to_remove)
es$CC_hyb <- es$CC_hyb[rownames(es$CC_hyb) != to_remove,]
# remove from rep1:
si$CC_hyb_rep1 <- filter(si$CC_hyb_rep1, sample_name != to_remove)
es$CC_hyb_rep1 <- es$CC_hyb_rep1[rownames(es$CC_hyb_rep1) != to_remove,]

# filtering on top 3000 most connected (highest degree) genes based on average connectivity.
# (we decide on softPower in the parameter fitting section, but it's needed for 
# soft connectivity, so we'll just do the default value of 6 for now)
softPower <- 6
avgConn <- map_dfr(es, softConnectivity, type = "unsigned", power = softPower) %>% 
  rowMeans() # I have no idea why it's rowMeans not colMeans... By my calculation map_dfr should get a 16 x nGenes matrix, so colMeans should return the mean connectivity of each gene but instead it gives a vector of length 16...
hist(avgConn)
keepGenes <- rank(-avgConn) <= 3000
es <- lapply(es, function(x) {
  return(x[,keepGenes])
})

#### Parameter Fitting ####
# pick power parameter ß for determining soft threshold of edge presence in our network
powers <- c(1:20)

# helper function for pickSoftThresholdRedgrave
testScaleFreeFit <- function(countsMatrix, countsMatrixName, pwr, nBreaks = 10) {
  k <- adjacency(countsMatrix, type = "unsigned", power = pwr) %>% rowSums()
  bins <- cut(k, nBreaks)
  dk <- tapply(k, bins, mean) # d for discretized or like dx/dy I think
  p_dk <- as.vector(tapply(k, bins, length)/length(k))
  dk <- dk[!is.na(dk)] # have to remove NAs manually b/c tapply is annoying
  p_dk <- p_dk[!is.na(p_dk)]
  log_dk <- as.vector(log10(dk))
  log_p_dk <- as.numeric(log10(p_dk))
  lm1 <- lm(log_p_dk ~ log_dk)
  plot(x = log_dk, y = log_p_dk, main = paste0(countsMatrixName, "\npower = ", pwr))
  lines(x = log_dk, y = predict(lm1))
  meank <- mean(k)
  slp <- round(lm1$coefficients[[2]], 2)
  return(list(Rsq = summary(lm1)$r.squared, MeanK = meank, slope = slp))
}
  
# Correcting a behavior found in WGCNA::pickSoftThreshold:
# To estimate the R^2 of the scale free plot (log-log plot of degree distribution), 
# we must bin the actual degree distribution and quantify how many genes are in each bin.
# BUT WGCNA::pickSoftThreshold doesn't throw out bins with 0 genes, it simply corrects for log(0)
# bins by adding a minute value (log(0 = 1e-9)), meaning that all the empty bins will all be in
# a row arbitrarily at -9 and throw off the R^2 of the whole thing!
# I have no idea why it does this
pickSoftThresholdRedgrave <- function(cts, cts_name, powerVector, nBreaks = 10) {
  output <- map_dfr(powerVector, testScaleFreeFit, nBreaks = nBreaks, countsMatrix = cts, countsMatrixName = cts_name)
  output$power <- powerVector
  return(output)
}
# This takes like 10 minutes! The WGCNA version is definitely faster. Sure wish it didn't retain 0 bins...
scaleFreeFits <- map2(.x = es, .y = names(es), .f = pickSoftThresholdRedgrave, powerVector = powers)

# degree distribution plots using different ß-based soft thresholds

# wrapping this WGCNA code in a function for better repeatability
plotSoftPower <- function(sft, sft_name) {
  # plotting results to pick power
  par(mfrow=c(1,2))
  cex1 = 0.9
  # Scale-free topology fit index vs soft-thresholding power
  plot(sft$power, -sign(sft$slope)*sft$Rsq, # I think this is just to make sure that the slope is negative? Otherwise you could just use the 2nd column without the 3rd
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence\n", sft_name))
  text(sft$power, -sign(sft$slope)*sft$Rsq,
       labels=sft$power,cex=cex1,col="red")
  abline(h=0.9, col="red") # this line corresponds to using an R^2 cut-off of h apparently
  # Mean connectivity vs soft-thresholding power
  plot(sft$power, sft$MeanK,
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$power, sft$MeanK, labels=sft$power, cex=cex1,col="red")
}
dummy <- map2(scaleFreeFits, names(scaleFreeFits), function(x, y) {
  return(plotSoftPower(x, y))
}) # get ready for 50 plots once again

# LowPi hybrid seems the weirdest, so let's look at that specifically
pickSoftThresholdRedgrave(es$LowPi_hyb, "LowPi_hyb", powerVector = powers) # definitely something weird happening for the highly connected genes

# 6 seems like a good tradeoff between the LowPi samples that would prefer 3 and the rest of them that would prefer 6-8
softPower <- 6

##### Data exploration ####
# see how well average expression level of each gene and overall connectivity of the inferred networks correlate between data sets
correlateExpressionAndConnectivity <- function(cts1, cts2, cts_name1, cts_name2) {
  rankExpr1 <- rank(rowMeans(t(cts1)))
  rankExpr2 <- rank(rowMeans(t(cts2)))
  rankConn1 <- rank(softConnectivity(cts1, type="unsigned", power=softPower))
  rankConn2 <- rank(softConnectivity(cts2, type="unsigned", power=softPower))
  plotTibble <- tibble("rankExpr1" = rankExpr1,
                       "rankExpr2" = rankExpr2,
                       "rankConn1" = rankConn1,
                       "rankConn2" = rankConn2)
  exprCorr <- cor(rankExpr1, rankExpr2, method = "pearson")
  exprPvalue <- cor.test(rankExpr1, rankExpr2, method = "pearson")$p.value
  connCorr <- cor(rankConn1, rankConn2, method = "pearson")
  connPvalue <- cor.test(rankConn1, rankConn2, method = "pearson")$p.value
  p_expr <- ggplot(data = plotTibble, aes(x = rankExpr1, y = rankExpr2)) + 
    xlab(paste("Ranked Expression", cts_name1)) + 
    ylab(paste("Ranked Expression", cts_name2)) + 
    ggtitle(paste("Correlation = ", exprCorr, "\n", "p-value = ", exprPvalue)) +
    geom_point()
  p_conn <- ggplot(data = plotTibble, aes(x = rankConn1, y = rankConn2)) + 
    xlab(paste("Ranked Connectivity", cts_name1)) + 
    ylab(paste("Ranked Conenctivity", cts_name2)) + 
    ggtitle(paste("Correlation = ", connCorr, "\n", "p-value = ", connPvalue)) +
    geom_point()
  return(list(expression = p_expr, connectivity = p_conn))
}
# comparisons: 50 networks would be insane to look at systematically. We'll pick two random ones to look at
# run these few lines to your heart's desire
library(gridExtra)
ntwk_idxs <- sample(c(1:50), 2, replace = FALSE)
ExprConn <- correlateExpressionAndConnectivity(es[[ntwk_idxs[1]]], es[[ntwk_idxs[2]]],
                                   names(es)[ntwk_idxs[1]], names(es)[ntwk_idxs[2]])
grid.arrange(grobs = ExprConn, ncol=2)

#### Building networks ####

### Tree building: adjacency matrix, topical overlap (TOM), and heirachical clustering

# set soft power (decided in parameter fitting section)
softPower <- 6

nNets <- length(es)
nConditions <- 14

# getting adjacency matrix from count matrix (this takes a bit)
adjFromCounts <- function(cts) {
  adj <- adjacency(cts, type = "unsigned", power = softPower)
  diag(adj) <- 0
  return(adj)
}

# creating adjacency matrix and TOM
adjs <- lapply(es, adjFromCounts)
TOMs <- lapply(adjs, TOMsimilarity, TOMType = "unsigned")
dissTOMs <- lapply(TOMs, function(x) {
  return(1-x)
})
geneTrees <- lapply(dissTOMs, function(x) {
  return(hclust(as.dist(x), method = "average"))
})

sizeGrWindow(9,6)
#### Fig 1C: trees showing lack of modularity in CC and TFdel ####
for (i in 1:nNets) {
  plot(geneTrees[[i]], xlab="", sub="", main=PrettyNames[i],
       labels = FALSE, hang=0.04) 
}

# Verifying that this lack of modularity is because genes in TFdel and CC have less variation across samples in that experimental class than the other experimental classes
library(ggbeeswarm)
gene_vars_df <- tibble("gene_var" = numeric(0), "experiment" = character(0))
for (i in 1:nNets) {
  gene_vars_df <- bind_rows(gene_vars_df, tibble("gene_var" = colVars(es[[i]]),
                                                 "experiment" = gsub("_.*", "", names(es)[i])))
}
gene_vars_df$log_var <- log(gene_vars_df$gene_var)
ggplot(gene_vars_df, aes(x = experiment, y = log_var), alpha = 0.5) + 
  geom_quasirandom(aes(color = experiment)) + ylab("Expression \nVariance") + xlab("Experiment") 
map(es, function(x) return(mean(colVars(x)))) %>% unlist() %>% sort() # CC and TFdel do indeed have the lowest average variance across genes, but it's not the world's most obvious difference

# Due to the lack of variation exhibited by TFdel and CC, we will eliminate these from the module analyses (also eliminate TFdelxLowN b/c it's just lowN with more messiness)
CC_names <- names(es)[grepl("CC", names(es))]
TFdel_names <- names(es)[grepl("TFdel", names(es))]
keep <- !(names(es) %in% c(CC_names, TFdel_names))

esfull <- es
sifull <- si
adjsfull <- adjs
TOMsfull <- TOMs
dissTOMsfull <- dissTOMs
geneTreesfull <- geneTrees

es <- es[keep]
si <- si[keep]
adjs <- adjs[keep]
TOMs <- TOMs[keep]
dissTOMs <- dissTOMs[keep]
geneTrees <- geneTrees[keep]

nNets <- length(es)

# export to non-environment-specific network construction stage
counts <- list(cer = rbind(es$LowN_cer,
                           es$HAP4_cer,
                           es$LowPi_cer),
               par = rbind(es$LowN_par,
                           es$HAP4_par,
                           es$LowPi_par),
               hyb = rbind(es$LowN_hyb,
                           es$HAP4_hyb,
                           es$LowPi_hyb),
               cer_rep1 = rbind(es$LowN_cer_rep1,
                                es$HAP4_cer_rep1,
                                es$LowPi_cer_rep1),
               par_rep1 = rbind(es$LowN_par_rep1,
                                es$HAP4_par_rep1,
                                es$LowPi_par_rep1),
               hyb_rep1 = rbind(es$LowN_hyb_rep1,
                                es$HAP4_hyb_rep1,
                                es$LowPi_hyb_rep1),
               cer_rep2 = rbind(es$LowN_cer_rep2,
                                es$HAP4_cer_rep2,
                                es$LowPi_cer_rep2),
               par_rep2 = rbind(es$LowN_par_rep2,
                                es$HAP4_par_rep2,
                                es$LowPi_par_rep2),
               hyb_rep2 = rbind(es$LowN_hyb_rep2,
                                es$HAP4_hyb_rep2,
                                es$LowPi_hyb_rep2))

info <- list(cer = rbind(si$LowN_cer,
                         si$HAP4_cer,
                         si$LowPi_cer),
             par = rbind(si$LowN_par,
                         si$HAP4_par,
                         si$LowPi_par),
             hyb = rbind(si$LowN_hyb,
                         si$HAP4_hyb,
                         si$LowPi_hyb),
             cer_rep1 = rbind(si$LowN_cer_rep1,
                              si$HAP4_cer_rep1,
                              si$LowPi_cer_rep1),
             par_rep1 = rbind(si$LowN_par_rep1,
                              si$HAP4_par_rep1,
                              si$LowPi_par_rep1),
             hyb_rep1 = rbind(si$LowN_hyb_rep1,
                              si$HAP4_hyb_rep1,
                              si$LowPi_hyb_rep1),
             cer_rep2 = rbind(si$LowN_cer_rep2,
                              si$HAP4_cer_rep2,
                              si$LowPi_cer_rep2),
             par_rep2 = rbind(si$LowN_par_rep2,
                              si$HAP4_par_rep2,
                              si$LowPi_par_rep2),
             hyb_rep2 = rbind(si$LowN_hyb_rep2,
                              si$HAP4_hyb_rep2,
                              si$LowPi_hyb_rep2))

# export
save(counts, info, file = "data_files/Cleaned_Data_StarvationsOnly.RData")

# TODO: Everything that follows is duplciated from the non-environment-specific network
# construction, which the above dataset is being imported into anyway,
# so it's probably all worth deleting, once it's clear it's not needed (once we have final networks from the non-env-specific version)

### Module selection: tree cutting
minModuleSize <- 15

#### Clustering exploration ####
# We're gonna do a bit of a March Madness situation: first find the best parameters 
# for each of the three tree-cutting algorithms: static, dynamic tree (def needs a better name), or hybrid, 
# then find which algorithm with its best parameters is best over all

# dynamic tree with 2 values of deepSplit (T or F) and three different max heights
colorList_Dynamic <- vector(mode = "list", length = nNets)
for (i in 1:nNets) {
  cat("currently processing ", i, "/", nNets, "\n")
  for (h in c(0.85, 0.9, 0.99, 0.995)) {
    deepSplitTrue <- cutreeDynamicTree(dendro = geneTrees[[i]], 
                             deepSplit = TRUE, maxTreeHeight = h,
                             minModuleSize = minModuleSize) %>% labels2colors() # Note: cutreeDynamic defaults to cutreeHybrid! Very confusing!! cutreeDynamicTree is the one that actually does the dynamic algorithm talked about it Langfelder et al. 2007
    colorList_Dynamic[[i]] <- cbind(colorList_Dynamic[[i]], deepSplitTrue)
  }
  for (h in c(0.85, 0.9, 0.99, 0.995)) {
    deepSplitFalse <- cutreeDynamicTree(dendro = geneTrees[[i]], 
                                        deepSplit = FALSE, maxTreeHeight = h,
                                        minModuleSize = minModuleSize) %>% labels2colors()
    colorList_Dynamic[[i]] <- cbind(colorList_Dynamic[[i]], deepSplitFalse)
  }
}
# plotting
for (i in 1:nNets) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(unlist(colorList_Dynamic[[i]])),
                      groupLabels = c("deepSplit: True, height: 0.85",
                                      "deepSplit: True, height: 0.9",
                                      "deepSplit: True, height: 0.99",
                                      "deepSplit: True, height: 0.995",
                                      "deepSplit: False, height: 0.85",
                                      "deepSplit: False, height: 0.9",
                                      "deepSplit: False, height: 0.99",
                                      "deepSplit: False, height: 0.995"), 
                      main = names(geneTrees)[i],
                      dendroLabels = FALSE)
} # deepSplit = False, height = 0.995 seems best

# hybrid with different deepSplits
colorList_Hybrid <- vector(mode = "list", length = nNets)
for (i in 1:nNets) {
  for (ds in 0:3) {
    chopped_tree <- cutreeHybrid(dendro = geneTrees[[i]], 
                                 distM = dissTOMs[[i]],
                                 deepSplit = ds, cutHeight = 0.99,
                                 minClusterSize = minModuleSize,
                                 pamRespectsDendro = FALSE)
    hybrid <- chopped_tree$labels %>% labels2colors()
    colorList_Hybrid[[i]] <- cbind(colorList_Hybrid[[i]], hybrid)
  }
  # same as above without the PAM stage (Partitioning Around Medoids)
  for (ds in 0:3) {
    chopped_tree <- cutreeHybrid(dendro = geneTrees[[i]], 
                                 distM = dissTOMs[[i]],
                                 deepSplit = ds, cutHeight = 0.99,
                                 minClusterSize = minModuleSize,
                                 pamStage = FALSE)
    hybrid <- chopped_tree$labels %>% labels2colors()
    colorList_Hybrid[[i]] <- cbind(colorList_Hybrid[[i]], hybrid)
  }
}
# compare dynamic and hybrid tree cutting (module selection) methods with different deepSplit values
for (i in 1:nNets) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(unlist(colorList_Hybrid[[i]])),
                      groupLabels = c("deepSplit: 0, pamStage: True",
                                      "deepSplit: 1, pamStage: True",
                                      "deepSplit: 2, pamStage: True",
                                      "deepSplit: 3, pamStage: True",
                                      "deepSplit: 0, pamStage: False",
                                      "deepSplit: 1, pamStage: False",
                                      "deepSplit: 2, pamStage: False",
                                      "deepSplit: 3, pamStage: False"), 
                      main = names(geneTrees)[i],
                      dendroLabels = FALSE)
} # Yes we definitely want something between all or no PAM step to avoid the turquoise megamodule absorbing all grey genes. deepSplit = 1 looks best (0 also works)

# hybrid with different maxPamDists
colorList_HybridPAMvar <- vector(mode = "list", length = nNets)
for (i in 1:nNets) {
  # first no pam stage
  chopped_tree <- cutreeHybrid(dendro = geneTrees[[i]], 
                               distM = dissTOMs[[i]],
                               deepSplit = 1, cutHeight = 0.99,
                               minClusterSize = minModuleSize,
                               pamStage = FALSE)
  hybrid <- chopped_tree$labels %>% labels2colors()
  colorList_HybridPAMvar[[i]] <- cbind(colorList_HybridPAMvar[[i]], hybrid)
  # now, different max dists
  for (ds in c(0, 0.5, 0.75, 0.9)) {
    chopped_tree <- cutreeHybrid(dendro = geneTrees[[i]], 
                                 distM = dissTOMs[[i]],
                                 deepSplit = 1, cutHeight = 0.99,
                                 minClusterSize = minModuleSize,
                                 pamRespectsDendro = FALSE,
                                 maxPamDist = ds)
    hybrid <- chopped_tree$labels %>% labels2colors()
    colorList_HybridPAMvar[[i]] <- cbind(colorList_HybridPAMvar[[i]], hybrid)
  }
}
# compare dynamic and hybrid tree cutting (module selection) methods with different deepSplit values
for (i in 1:nNets) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(unlist(colorList_HybridPAMvar[[i]])),
                      groupLabels = c("no PAM",
                                      "PAMmaxdist = 0",
                                      "PAMmaxdist = 0.5",
                                      "PAMmaxdist = 0.75",
                                      "PAMmaxdist = 0.9"), 
                      main = names(geneTrees)[i],
                      dendroLabels = FALSE)
} 
# When maxPamDist is any value less than the cluster radius, PAM only adds genes to a cluster
# when the average distance between that gene and each object in the cluster is less than the cluster's radius.
# Cluster radius is defined as the maximum of the average distance between each object
# in the cluster and every other object (so each object has one average distance to all 
# other objects, then the radius is the max across all objects)
# That is all to say that it appears we can't change this parameter...
# still feels weird that it always assigns stray genes to the turquoise module though. I wish
# it would assign them to the nearest core instead

# let's also try some static cuts with different heights
colorList_Static <- vector(mode = "list", length = nNets)
for (i in 1:nNets) {
  for (h in c(0.95, 0.925, 0.9)) {
    static <- cutreeStaticColor(dendro = geneTrees[[i]],
                                cutHeight = h,
                                minSize = minModuleSize)
    colorList_Static[[i]] <- cbind(colorList_Static[[i]], static)
  }
}

for (i in 1:nNets) {
  plotDendroAndColors(dendro = geneTrees[[i]], 
                      colors = cbind(unlist(colorList_Static[[i]])),
                      groupLabels = c("0.95", "0.925", "0.9"), 
                      main = names(geneTrees)[i],
                      dendroLabels = FALSE)
} # RIP. Doesn't work at all on nested trees like HAP4_cer (which was exactly why the dynamic methods were developed in the first place)

# Anna's version: hybrid to ID cluster cores
# then kMeans with these cores as the seeds
# @intput: geneTree and corresponding dissTOM for cutreeHybrid
# @output: vector of indices for one random gene in each module
findSeeds <- function(tree, dists) {
  chopped_tree <- cutreeHybrid(dendro = tree,
                               distM = dists,
                               deepSplit = 1, cutHeight = 0.99,
                               minClusterSize = minModuleSize,
                               pamStage = FALSE)
  clrs <- chopped_tree$labels %>% labels2colors()
  # visualizing 
  plotDendroAndColors(dendro = tree, colors = clrs, dendroLabels = FALSE, main = paste(length(unique(clrs)), "seeds will be generated"))
  seeds <- vector(mode = "numeric", length = 0L)
  for(clr in setdiff(unique(clrs), "grey")) {
    seeds <- c(seeds, sample(which(clrs == clr), 1))
  }
  return(seeds)
}

# tests for findSeeds
test <- findSeeds(geneTrees$LowN_rep1_cer, dissTOMs$LowN_rep1_cer)
test

# clusters all genes by 
kMeansByDist <- function(distM, clusters, clusterDists, distThresh = 0.99) {
  tripsThresh <- lapply(clusterDists, function(x) {
    return(min(x[-unlist(clusters)]) > distThresh)
  }) %>% unlist()
  if(all(tripsThresh)) {
    return(clusters)
  }
  for(i in which(!tripsThresh)) {
    newIdx <- which(clusterDists[[i]] == min(clusterDists[[i]][-unlist(clusters)])) # not allowing the new index to be one of the ones already assigned to a cluster
    clusterDists[[i]] <- pmean(clusterDists[[i]], distM[, newIdx])
    clusters[[i]] <- c(clusters[[i]], newIdx)
  }
  return(kMeansByDist(distM, clusters, clusterDists, distThresh))
}

clusterListToColorVector <- function(clusters, nGenes = 3000) {
  colorVector <- rep("grey", nGenes)
  colorPalette <- c("magenta","black","red","turquoise","purple","green","cyan","yellow","brown","blue",      
                    "lightgreen","tan","lightcyan","greenyellow","lightyellow,midnightblue","pink","salmon")
  if(length(colorPalette) < length(clusters)) {
    stop(error = "make up some new colors dingus")
  }
  for(i in 1:length(clusters)) {
    colorVector[clusters[[i]]] <- colorPalette[i]
  }
  return(colorVector)
}
# testing with my favorite nested tree, HAP4_cer
seeds_cer_HAP4 <- findSeeds(geneTrees$HAP4_cer, dissTOMs$HAP4_cer)
test_clusters <- as.list(seeds_cer_HAP4)
test_clusterDists <- lapply(test_clusters, function(x) {
  return(dissTOMs$HAP4_cer[, x]) # it's symmetrical so arbitrary that I chose column x not row x. Just makes more sense to me I guess
})
test <- kMeansByDist(dissTOMs$HAP4_cer, test_clusters, test_clusterDists)
test_colors <- clusterListToColorVector(test)

plotDendroAndColors(geneTrees$HAP4_cer, test_colors, dendroLabels = FALSE) # wow coming up with a good clustering algorithm is harder than it looks
# TODO: replace random selection of seed gene with minimum height member of each initial hybrid cluster
table(test_colors) # TODO: figure out why the clusters are so evenly partitioned

#### Clustering exploration over ####
# We're going with dynamic tree deepSplit = FALSE, maxTreeHeight = 0.995
colors <- lapply(geneTrees, function(x) {
  chopped_tree <- cutreeDynamicTree(dendro = x,
                                    maxTreeHeight = 0.995,
                                    deepSplit = FALSE,
                                    minModuleSize = minModuleSize) %>% labels2colors()
  return(chopped_tree)
})

# Within-species module preservation between experimental conditions
for (i in seq(from = 1, to = 20, by = 4)) { # the lowN trees have the nicest separation, so they seem like the best tree reference to use
  plotDendroAndColors(dendro = geneTrees[[i]], colors =  do.call(cbind, colors[i:(i+3)]),
                      groupLabels = names(colors[i:(i+3)]), main = names(geneTrees)[i],
                      dendroLabels = FALSE)
} # We can revisit this after merging close modules and synching colors

### Merging close modules
# initial METree construction
MEs <- map2(es, colors, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y)
  MEs <- MEoutput$eigengenes
  return(MEs)
})
MEDiss <- lapply(MEs, function(x) {
  return(1-abs(cor(x))) # Anna made this absolute value of cor to keep heights between 0 and 1 again (plus significant negative correlations are also meaningful, and we already used unsigned network in the first place, so it's more consistent)
})
METree <- lapply(MEDiss, function(x) {
  return(hclust(as.dist(x), method = "average"))
})

# Taking advantage of replicate networks to determine a threshold for module eigengene similarity:
# we should merge any modules that have about the same ME difference as homologous modules from replicate networks

# Function to relate the total number of clusters generated vs mod overlap 
# (% genes in same module, including grey)
# Designed to be used for replicate pair, where hopefully there will be significant 
# module overlap at some number of clusters
relateNclustAndModuleOverlap <- function(es1, es2, colors1, colors2) {
  # checking which set of modules is bigger then matching colors to bigger set
  if(length(unique(colors1)) >= length(unique(colors2))) {
    colors2 <- matchLabels(source = colors2, reference = colors1) # Note: this gives colors 1 the dimensions 3000 1 when colors1 has dimensions NULL length 3000, so maybe keep note of that in case it matters
  }
  if(length(unique(colors1)) < length(unique(colors2))) {
    colors1 <- matchLabels(source = colors1, reference = colors2)
  }
  # making our happy little trees
  MEoutput1 <- moduleEigengenes(es1, colors = colors1, excludeGrey = TRUE)
  MEoutput2 <- moduleEigengenes(es2, colors = colors2, excludeGrey = TRUE)
  tree1 <- (1-abs(cor(MEoutput1$eigengenes))) %>% as.dist() %>% hclust(method = "average")
  tree2 <- (1-abs(cor(MEoutput2$eigengenes))) %>% as.dist() %>% hclust(method = "average")
  
  # collecting mod overlap for each module
  output <- data.frame(matrix(nrow = 0, ncol = 4))
  names(output) <- c("n_clust", "branch_height1", "branch_height2", "mod_overlap")
  
  hs1 <- sort(tree1$height, decreasing = TRUE) # cutting at the branch height given in the hclust object does NOT break that branch in 2 (I know it's confusing, I need to use the edge case to have a one:to:one correspondence between nClust and height)
  hs2 <- sort(tree2$height, decreasing = TRUE)
  for (i in 1:max(length(hs1), length(hs2))) {
    
    h1 <- hs1[i]
    h2 <- hs2[i]
    
    if(is.na(h1)) {h1 <- 0}
    if(is.na(h2)) {h2 <- 0}

    newMods1 <- mergeCloseModules(es1, colors1, cutHeight = h1, useAbs = TRUE, iterate = FALSE) # use Abs takes absolute value of ME correlation, iterate = FALSE just does one merge instead of continuing to iterate at a given cut height until modules no longer change (which could be useful I guess, but it's hard to standardize between networks)
    newMods2 <- mergeCloseModules(es2, colors2, cutHeight = h2, useAbs = TRUE, iterate = FALSE)
    newColors1 <- newMods1$colors
    newColors2 <- newMods2$colors %>% matchLabels(reference = newColors1)
    
    
    output <- rbind(output, data.frame("n_clust1" = (length(unique(newColors1)) - 1),
                                       "n_clust2" = (length(unique(newColors2)) - 1),
                                       "branch_height1" = h1,
                                       "branch_height2" = h2,
                                       "mod_overlap" = sum(newColors1 == newColors2)/length(newColors1)))
    
  }
  return(output)
}


# tests for relateNclustAndModuleOverlap
# test1 is hyc, where rep2 has more initial modules than rep1
test1 <- relateNclustAndModuleOverlap(es$LowN_rep1_hyc, 
                                     es$LowN_rep2_hyc, 
                                     colors$LowN_rep1_hyc, 
                                     colors$LowN_rep2_hyc)
plot(METree$LowN_rep1_hyc, main = "rep1")
abline(h = test1$branch_height1, col = "red")
test1$n_clust1 # visually verify that that's the heights are at the branches and the number of clusters from each cut correspond to the values in n_clust columns
plot(METree$LowN_rep2_hyc, main = "rep2")
abline(h = test1$branch_height2, col = "red")
test1$n_clust2

# test2 is cer, where rep1 has more initial modules than rep2
test2 <- relateNclustAndModuleOverlap(es$LowN_rep1_cer, 
                                      es$LowN_rep2_cer, 
                                      colors$LowN_rep1_cer, 
                                      colors$LowN_rep2_cer)
plot(METree$LowN_rep1_cer, main = "rep1")
abline(h = test2$branch_height1, col = "red")
test2$n_clust1 # good luck visually verifying this one
plot(METree$LowN_rep2_cer, main = "rep2")
abline(h = test2$branch_height2, col = "red")
test2$n_clust2

# ok now actually collecting comparisons for each pair of replicates
rep1_idxs <- seq(from = 1, to = 20, by = 4)
branch_vs_modOverlap_df <- vector(mode = "list", length = length(rep1_idxs))
for(i in 1:length(rep1_idxs)) { 
  idx <- rep1_idxs[i]
  cat(i, "/", length(rep1_idxs), ":", gsub("_rep1_", "", names(es)[idx]), "\n")
  branch_vs_modOverlap_df[[i]] <- relateNclustAndModuleOverlap(es[[idx]], es[[idx+1]], colors[[idx]], colors[[idx+1]])
}

# function to convert between a given gene tree cut height and 
# the number of clusters it would generate from the tree
# @input: query cut height h, vector of branch heights (like from <hclust object>$height),
# corresponding vector containing the number of clusters resulting from a cut at 
# each branch height
# @output: number of clusters generated (scalar integer), not including grey
height2NClust <- function(h, heights, nclusts) {
  nHeights <- length(heights)
  whichVec <- rep(FALSE, nHeights)
  if (h >= heights[1]) {
    return(nclusts[1])
  }
  if (h < heights[nHeights]) {
    return(nclusts[nHeights])
  }
  for (i in 2:nHeights) {
    whichVec[i] <- (h >= heights[i]) & (h < heights[i-1])
  }
  if (sum(whichVec) != 1) {
    stop(gettext(paste("no definitive bounding heights found", sum(whichVec))))
  }
  return(nclusts[whichVec])
}

# tests for height2NClust
# pick whatever pair of replicates you like, plot them and visually count clusters, then verify that height2NClust reports the same number of clusters (remembering that grey genes don't count as a cluster)
test <- relateNclustAndModuleOverlap(es$LowN_rep1_cer, es$LowN_rep2_cer, colors$LowN_rep1_cer, colors$LowN_rep2_cer)
plot(METree$LowN_rep1_cer, main = "test for height2NClust")
abline(h = 0.2, col = "red")
height2NClust(0.2, test$branch_height1, test$n_clust1) == 4
plot(METree$LowN_rep2_cer, main = "test for height2NClust")
abline(h = 0.2, col = "red")
height2NClust(0.2, test$branch_height2, test$n_clust2) == 3

# given a branch vs mod overlap dataset, produces a ggplot of cut height vs % module overlap between replicates
# draws a vertical line at the query height
# TODO: not sure what the best way to plot both cut heights vs mod overlap. We want to say "if you cut at this height, you will get this much mod overlap" and we're not quite there yet. I think the issue is that I didn't cut at the same height in the relateNclustAndModuleOverlap function and I probably should've. But I went with min for now, becasue cutting at min doesn't cut either one at that height 
createCutHeightVsModOverlapPlots <- function(d, name, queryHeight = 0.2) {
  p <- ggplot(d, aes(x = pmin(branch_height1, branch_height2), y = mod_overlap)) + 
    geom_line() + 
    geom_point() +
    ggtitle(name) +
    geom_vline(xintercept = queryHeight, color = "red") +
    xlab("cut height") +
    ylab("% genes in same module")
  return(p)
}

# same as above, but x axis is number of clusters instead of cut height
# draws one vertical line for each replicate at the number of clusters a 
# cut at the query height would generate
# TODO: this isn't totally accurate yet b/c it relies on average nclusters vs mod overlap and the vertical lines aren't actually at the correct mod overlap (I think it's the midpoint between the lines?)
createNClustVsModOverlapPlots <- function(d, name, queryHeight = 0.2) {
  nClust1 <- height2NClust(queryHeight, d$branch_height1, d$n_clust1)
  nClust2 <- height2NClust(queryHeight, d$branch_height2, d$n_clust2)
  p <- ggplot(d, aes(x = pmean(n_clust1, n_clust2), y = mod_overlap)) + 
    geom_line() + 
    geom_point() +
    ggtitle(name) +
    geom_vline(xintercept = nClust1, color = "red") +
    geom_vline(xintercept = nClust2, color = "blue") +
    xlab("number of clusters") +
    ylab("% genes in same module")
  return(p)
}

# visualizing module preservation by cut height
pHeight <- map2(branch_vs_modOverlap_df, gsub("_rep1_", "", names(es[rep1_idxs])), createCutHeightVsModOverlapPlots, queryHeight = 0.05) # try out different query heights here
pNClust <- map2(branch_vs_modOverlap_df, gsub("_rep1_", "", names(es[rep1_idxs])), createNClustVsModOverlapPlots, queryHeight = 0.05)
library(gridExtra)
grid.arrange(grobs = pHeight, nrow = 3)
grid.arrange(grobs = pNClust, nrow = 3)

# 0.05 looks good. Gets rid of crazy modules in cer but not much on any other one
MEDissThresh <- 0.05

# Final check of proposed cut height
# (everything with a branch point at a lower height than the abline gets merged, except grey)
par(mfrow = c(1,1))
for (i in 1:nNets) {
  plot(METree[[i]], main = names(es)[i], xlab = "", sub = "")
  abline(h=MEDissThresh, col = "red")
} # grey won't be merged by mergeCloseModules (Note: for some reason mergeCloseModules doesn't include MEgrey in its dendro even though it doesn't merge grey genes)

# actually merging
unmergedColors <- colors
unmergedMEs <- MEs
merge <- map2(es, colors, function(x,y) {
  return(mergeCloseModules(x, y, cutHeight = MEDissThresh, verbose = 3))
})
colors <- lapply(merge, function(x) {
  return(x$colors)
})


# matching module color labels
# first match replicate pairs within species
for (i in rep1_idxs) {
  colors[[i+1]] <- matchLabels(source = colors[[i+1]], reference = colors[[i]])
}
for (i in rep1_idxs) {
  plotDendroAndColors(geneTrees[[i]], cbind(colors[[i]], colors[[i+1]]),
                      main = names(geneTrees)[i],
                      groupLabels = names(es)[c(i, i+1)],
                      dendroLabels = FALSE)
} # looking good!



# before matching across species
species_idxs <- seq(from = 0, to = nNets-nConditions, by = nConditions)
for (i in 1:nConditions) {
  plotDendroAndColors(geneTrees[[i]], do.call(cbind, colors[i + species_idxs]),
                      main = names(geneTrees)[i],
                      groupLabels = c(names(es)[i + species_idxs], "colorRef"),
                      dendroLabels = FALSE)
}

# matching module colors across networks
# matchLabels looks for modules with significantly overlapping gene content and makes them the same color (not as sophisticated as the permutation tests in modulePreservation, it's just a fisher exact test)
colorsUnmatched <- colors

# TODO: create a match labels that can take n color vectors at once
colors <- colorsUnmatched

# helper function for pmacthLabels
# given a set of vectors of 1s and 0s indicating whether each gene belongs to the current module being compared or not,
# returns a pvalue obtained from a Fisher's Exact test of the following 2 x 2 contingency table


# function to match color labels across n color vectors at once
pmatchLabels <- function(cs) {
  # first figure out what genes are in a group together the most and identify consensus groups
  
  # then give these consensus groups colors
  
  # then do groups only in a subset of networks
  
  # then do groups only found in lone networks
}

# recalculate MEs post merge/color match
MEs <- map2(es, colors, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y)
  MEs <- MEoutput$eigengenes
  return(MEs)
})

# final network visualizations
for (i in 1:nNets) {
  plotDendroAndColors(geneTrees[[i]], colors[[i]],
                      names(geneTrees)[i], main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}

#TOMplots

# oh my goodness what a ride it has been
# this is far from the ideal heatmap/dendrogram but at least it allows for multiple colors
TOMplotSubsetMultiColor <- function(dissTOM, colorList, sampleSize = 1000, plotTitle = "TOM plot") {
  nGenes <- ncol(dissTOM)
  randomSample <- sample(nGenes, sampleSize)
  selectTOM <- dissTOM[randomSample, randomSample]
  plotTOM <- selectTOM^7 # I guess 7 is just a power they thought looked pretty
  diag(plotTOM) <- NA
  plotTree <-  hclust(as.dist(plotTOM), method = "average")
  layout_matrix <- cbind(c(2,3,0), c(0,0,5), c(6,1,4))
  layout(layout_matrix, widths = c(1,1,4), heights = c(1,2,4))
  colorMat <- NULL
  for (i in 1:length(colorList)) {
    colorMat <- cbind(colorMat, colorList[[i]][randomSample])
  }
  plotOrderedColors(plotTree$order, colorMat, rowLabels = NA, align = "edge")
  TOMplot(-plotTOM, plotTree, Colors = NULL, main = plotTitle, setLayout = FALSE, margins = c(2,2))
}

# TOMplots just for RStudio
TOMplotTitles <- c("S. cerevisiae", "S. paradoxus", "F1 hybrid, cerevisiae allele", "F1 hybrid, paradoxus allele")
par(mar = c(2,2,2,2)) # don't ask me why you need this
TOMplotSubsetMultiColor(dissTOM = dissTOMs[[1]], 
                        colorList = colors, plotTitle = TOMplotTitles[1], sampleSize = 100) # sacrificial plot --- for some delicious reason, the first plot has its module definitions lined up incorrectly, but it's fine when you re-run it

# for all4
for(i in 1:4) {
  TOMplotSubsetMultiColor(dissTOM = dissTOMs[[i]], 
                          colorList = colors, plotTitle = TOMplotTitles[i], sampleSize = 1000)
} 

# for just cer and par
for(i in 1:2) {
  TOMplotSubsetMultiColor(dissTOM = dissTOMs[[i]], 
                          colorList = list(colors[[1]], colors[[2]]), plotTitle = TOMplotTitles[i], sampleSize = 1000)
} 

# believe it or not, this is all for the legend
library(circlize)
library(ComplexHeatmap)
colorFunction <- colorRamp2(breaks = seq(from = 0, to = 1, length.out = 12), colors = hcl.colors(12, "YlOrRd", rev = TRUE))
lgd <- Legend(col_fun = colorFunction, title = "Topological Overlap")
pushViewport(viewport(width = 0.9, height = 0.9))
grid.rect()  # border
draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
popViewport()

# it has come to my attention that the default WGCNA colors are all in RGB and look pretty muted in CMYK
# sooo for posters and printing, here's how you can convert colors if you so choose
# library(scales)
# AnnasHandpickedPalette <- c("aquamarine", "tomato", "yellow", "seagreen", "darkorchid", "deepskyblue","lawngreen", "midnightblue", "deeppink", "thistle", "saddlebrown")
# show_col(AnnasHandpickedPalette)
# ColorConversionTable <- rbind(c("grey", "grey"), cbind(setdiff(unique(c(colors$cer, colors$par, colors$hyc, colors$hyp)), "grey"), AnnasHandpickedPalette))
# changeColor <- function(oldColors, conversionTable) {
#   newColorIdx <- which(conversionTable[, 1] == oldColors)
#   return(conversionTable[newColorIdx, 2])
# }
# newColors <- lapply(colors, function(x) {
#   return(sapply(x, changeColor, conversionTable = ColorConversionTable))
# })

# outputting TOMplots (Note this DOES NOT WORK with the multicolor! It misaligns the 4 color rows)
TOMplotFileNames <- c("TOMplot_cer.pdf", "TOMplot_par.pdf", "TOMplot_hyc.pdf", "TOMplot_hyp.pdf")
for (i in 1:4) {
  pdf(file = paste0("../../../posters/", TOMplotFileNames[i]), colormodel = "cmyk")
  TOMplotSubset(dissTOM = dissTOMs[[i]], 
                colorList = colors, plotTitle = TOMplotTitles[i], sampleSize = 1000)
  dev.off()
} 

# We're done building our networks! Let's save all the relevant files
save(file = "Networks.RData", counts_all4, counts_top3000, adjs, TOMs, dissTOMs,
     geneTrees, MEs, TOMs, colors, dissTOMs, softPower)