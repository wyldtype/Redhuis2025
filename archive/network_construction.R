#### Network construction from Gene Expression Data ####
# Following WGCNA tutorials: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Following Jeremy Miller's tutorial on Meta-analysis with WGCNA: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/
sapply(c("dplyr", "tidyr", "ggplot2", "purrr", "edgeR", "WGCNA", "matrixStats"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")
options(stringsAsFactors = FALSE)

# TODO: How did YPL274W end up in a CCM? It's not expressed in paradoxus,
# So it shouldn't have been grouped into any module

# TODO: YHR035W has a stronger correlation to the brown module in cerevisiae
# and the black module in paradoxus, even though it was placed into
# the black module in cerevisiae and the red module in paradoxus
# how is this possible if we're grouping genes based on their correlation?
# topological overlap?

#### Reading & Cleaning Data ####
# read in count data (normalized) and sample info
# cer, par
load("data_files/Cleaned_Redhuis_Data.RData")
# hyc, hyp
load("data_files/Cleaned_Redhuis_Data_AlleleSpecific.RData")

#### optional: including Fay et al. 2023 temp data ####
# cer, par
load("data_files/Cleaned_Fay_Counts.RData")
fay <- fay[rownames(fay) %in% rownames(counts),] # ~500 mostly lowly expressed genes filtered out
# hyc hyp
load("data_files/Cleaned_Fay_Counts_allele.RData")
fay_allele <- fay_allele[rownames(fay_allele) %in% rownames(counts),]
# there's no easy way to left_join matrices unfortunately
missing_fay_genes <- rownames(counts)[!(rownames(counts) %in% rownames(fay))]
missing_fay <- matrix(NA, nrow = length(missing_fay_genes), ncol = ncol(fay))
rownames(missing_fay) <- missing_fay_genes
colnames(missing_fay) <- colnames(fay)
missing_fay_allele <- matrix(NA, nrow = length(missing_fay_genes), ncol = ncol(fay_allele))
rownames(missing_fay_allele) <- missing_fay_genes
colnames(missing_fay_allele) <- colnames(fay_allele)
fay <- rbind(fay, missing_fay)
fay_allele <- rbind(fay_allele, missing_fay_allele)
fay <- fay[rownames(counts),]
fay_allele <- fay_allele[rownames(counts),]
sum(rownames(fay) == rownames(counts))/nrow(counts)
sum(rownames(fay_allele) == rownames(counts_allele))/nrow(counts_allele) # cbind doesn't check that rownames match, so make sure they match beforehand
# joining
counts <- cbind(counts, fay)
counts_allele <- cbind(counts_allele, fay_allele)
sample_info <- bind_rows(sample_info, sample_info_fay) # bind rows thankfully does check that the colnames are the same
sample_info_allele <- bind_rows(sample_info_allele, sample_info_fay_allele)

#### partitioning and cleaning counts ####
# checking that there aren't any genes with lots of missing values or outliers
gsg <- goodSamplesGenes(t(cbind(counts, counts_allele)), verbose = 3)
gsg$allOK

# filtering down to WT, removing TF deletions
# Rationale: we want the modules to reflect WT relationships between genes.
# Very few genes respond to the TF deletions, so including them just
# inflates the influence of the LowN portion of the dataset
counts_all4 <- list(cer = t(counts[, sample_info$organism == "cer" & sample_info$genotype == "WT"]), # unlike DESeq2, WGCNA expects ROWS to be samples and COLUMNS to be genes, likely because that's what R's cor function expects to get a nGene x nGene cor matrix
                    par = t(counts[, sample_info$organism == "par" & sample_info$genotype == "WT"]),
                    hyc = t(counts_allele[, sample_info_allele$allele == "cer" & sample_info_allele$genotype == "WT"]),
                    hyp = t(counts_allele[, sample_info_allele$allele == "par" & sample_info_allele$genotype == "WT"]))

# filtering for WT samples with same set of conditions in all 4 species
sample_info <- filter(sample_info, genotype == "WT")
sample_info_allele <- filter(sample_info_allele, genotype == "WT")
symdiff(x = unique(sample_info$condition),
        y = unique(sample_info_allele$condition)) # conditions we're removing
common_conditions <- intersect(unique(sample_info$condition),
                               unique(sample_info_allele$condition))
sample_info <- filter(sample_info, condition %in% common_conditions)
sample_info_allele <- filter(sample_info_allele, condition %in% common_conditions)
counts_all4$cer <- counts_all4$cer[sample_info$sample_name[sample_info$organism == "cer"],]
counts_all4$par <- counts_all4$par[sample_info$sample_name[sample_info$organism == "par"],]
counts_all4$hyc <- counts_all4$hyc[sample_info_allele$sample_name[sample_info_allele$allele == "cer"],]
counts_all4$hyp <- counts_all4$hyp[sample_info_allele$sample_name[sample_info_allele$allele == "par"],]

# checking for outlier samples by Euclidean distance clustering
par(cex = 0.2);
par(mar = c(0,4,2,0))
for (i in 1:4) {
  sampleTree <- hclust(dist(counts_all4[[i]]), method = "average")
  plot(sampleTree, main = names(counts_all4)[i])
  abline(h = 50000, col = "red")
} # we're looking for samples that are way off on their own---good clustering has two initial branches both leading to a decent total number of samples (like at least 3, but it's a judgement call)
# # (note that the y axis is relative to the largest euclidean distance exhibited between any two samples, so it's just good luck that it's so consistent)
# pruneEucDistTree <- function(cts, cuttingHeight = 50000) {
#   tree <- hclust(dist(cts), method = "average")
#   clust <- cutreeStatic(tree, cutHeight = cuttingHeight, minSize = 10)
#   cat("removing", sum(clust == 2), "samples\n")
#   cat("keeping", sum(clust == 1), "samples\n")
#   keepSamples <- (clust == 1)
#   pruned_cts <- cts[keepSamples,]
#   return(pruned_cts)
# } 

# # Skipping because there doesn't seem to be a reason to filter for connected genes in our dataset
# # filtering on top 3000 most connected (highest degree) genes based on average connectivity.
# # (we decide on softPower in the parameter fitting section, but it's needed for 
# # soft connectivity, so we'll just do the default value of 6 for now)
# softPower <- 6
# avgConn <- map_dfr(counts_all4, softConnectivity, type = "unsigned", power = softPower) %>% 
#   rowMeans() # I have no idea why it's rowMeans not colMeans... By my calculation map_dfr should get a 16 x nGenes matrix, so colMeans should return the mean connectivity of each gene but it gives a vector of length 16...
# hist(avgConn, breaks = 50)
# keepGenes <- rank(-avgConn) <= 3000
# counts_top3000 <- lapply(counts_all4, function(x) {
#   return(x[,keepGenes])
# })

#### dcor similarity matrices ####
adjFromCounts <- function(.cts) {
  # adj <- cor(.cts, use = "pairwise.complete.obs")
  adj <- matrix(NA, nrow = ncol(.cts), ncol = ncol(.cts))
  for (i in c(1:(ncol(adj) - 1))) {
    cat(i, "/", ncol(adj), "\n")
    for (j in c((i+1):ncol(adj))) {
      non_na <- cbind(.cts[,i], .cts[,j])[!is.na(.cts[,i]) & !is.na(.cts[,j]),]
      dcor_ij <- dcor(non_na[,1], non_na[,2])
      adj[i, j] <- dcor_ij
      adj[j, i] <- dcor_ij
    }
  }
  diag(adj) <- 0
  return(adj)
}
test <- adjFromCounts(counts_all2$cer[,1:100])
max(test)
min(test[test > 0])
test[1:3,1:3]
test <- adjFromCounts_WGCNA(counts_all2$cer[1:10,1:100])
test[1:3,1:3]
max(test)
min(test)

# generating adj/similarity matrices (note this takes overnight b/c dcor is O(n^2) whereas cor is O(n))
# adjs <- lapply(list(counts_all2$cer, counts_all2$par), adjFromCounts)
# save(adjs, file = "data_files/Adjacency_dcor.RData")
load("data_files/Adjacency_dcor.RData")
names(adjs) <- c("cer", "par")


################# Pearson Correlation Unsigned Version ############ 

#### Parameter Fitting ####

# pick power parameter ß for determining soft threshold of edge presence in our network
powers <- seq(from = 2, to = 20, by = 1)
# sfts <- lapply(counts_all2, pickSoftThreshold, 
#                powerVector = powers, verbose = 5, networkType = "signed") # It really just calculates R^2 for different log(k) vs log(p(k)) (aka the classic scale free network) plots
sfts <- lapply(adjs, pickSoftThreshold.fromSimilarity,
               powerVector = powers, verbose = 5)

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
  abline(h=0.9, col="red") # this line corresponds to using an R^2 cut-off of h, fairly arbitrary threshold
  # Mean connectivity vs soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}
map2(sfts, names(sfts), plotSoftPower)

# setting softPower
softPower <- 5

#### Building networks ####
# Module selection: Adjacency matrix, topical overlap (TOM), and heirachical clustering

# soft power (decided in parameter fitting section)
softPower

# getting adjacency matrix from count matrix (this takes a bit)
adjFromCounts <- function(cts) {
  adj <- adjacency(cts, type = "signed", power = softPower,
                   corOptions = "use = 'pairwise.complete.obs'") # Genes with NA for some samples have those samples excluded
  diag(adj) <- 0
  return(adj)
}

# creating adjacency matrix and TOM
# adjs <- lapply(counts_all2, adjFromCounts)
# TODO: I do not understand why unsigned and 
# signed produce different results for dcor similarity matrics,
# which are all non-negative entries. But unsigned look much better
adjs <- lapply(adjs, adjacency.fromSimilarity, type = "unsigned",
               power = softPower)
TOMs <- lapply(adjs, TOMsimilarity, TOMType = "unsigned")
dissTOMs <- lapply(TOMs, function(x) {
  return(1-x)
})
geneTrees <- lapply(dissTOMs, function(x) {
  return(hclust(as.dist(x), method = "average"))
})

par(mfrow = c(length(geneTrees)/2, 2))
for (i in 1:length(geneTrees)) {
  plot(geneTrees[[i]], xlab="", sub="", main=names(geneTrees)[i],
       labels = FALSE, hang=0.04)
}

################## module selection ##################
minModuleSize <- 10

colors <- lapply(geneTrees, function(x) {
  chopped_tree <- cutreeDynamicTree(dendro = x,
                                    maxTreeHeight = 0.995,
                                    deepSplit = FALSE,
                                    minModuleSize = minModuleSize) %>% labels2colors()
  return(chopped_tree)
})
for (i in 1:length(geneTrees)) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = colors[[i]],
                      main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}

# merging close modules
# calculate eigengenes (eigenvector with largest eigenvalue aka 1st PC of expression matrix (woah way back to that thing) for each module)
# Note: colors is not an aesthetic argument, it's literally the only information on module membership
# moduleEigengenes for a p x n module X is equivalent to svd(t(scale(x)))$v[,1], the first principle component (length p) of X. If align = "along average", the sign is also potentially flipped so that the correlation between the average gene expression of each module across environments and its eigengene is non-negative (cor(A.E., M.E.) >= 0)
MEs <- map2(counts_all2, colors, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y, excludeGrey = TRUE)
  MEs <- MEoutput$eigengenes
  return(MEs)
})
MEDiss <- lapply(MEs, function(x) {
  return(1-cor(x))
})
METree <- lapply(MEDiss, function(x) {
  return(hclust(as.dist(x), method = "average"))
})

# merging close modules
# using multiple threshold to test how robust conclusions are to module definition
thresh35 <- 0.35
thresh25 <- 0.25 
thresh10 <- 0.10
# thresh00 <- 0.00 # not a real threshold, just to point out that the third condition is no merging

# visualizing potential new modules
# .35 thresh, merging at 65% similarity (extreme example of merging)
par(mfrow = c(length(MEs)/2, 2))
for (i in 1:length(MEs)) {
  plot(METree[[i]], main = names(counts_all4)[i], xlab = "", sub = "")
  abline(h=thresh35, col = "red")
} # grey won't be merged by mergeCloseModules

# .25 thresh, merging at 75% similarity
dev.off()
par(mfrow = c(length(MEs)/2, 2))
for (i in 1:length(MEs)) {
  plot(METree[[i]], main = names(counts_all4)[i], xlab = "", sub = "")
  abline(h=thresh25, col = "red")
}

# .10 thresh, merging at 90% similarity
dev.off()
par(mfrow = c(length(MEs)/2, 2))
for (i in 1:length(MEs)) {
  plot(METree[[i]], main = names(counts_all4)[i], xlab = "", sub = "")
  abline(h=thresh10, col = "red")
}

# actually merging
colors00 <- colors
unmergedMEs <- MEs
merge10 <- map2(counts_all2, colors, function(x,y) {
  return(mergeCloseModules(x, y, cutHeight = thresh10, verbose = 3))
})
colors10 <- lapply(merge10, function(x) {
  return(x$colors)
})
merge25 <- map2(counts_all2, colors, function(x,y) {
  return(mergeCloseModules(x, y, cutHeight = thresh25, verbose = 3))
})
colors25 <- lapply(merge25, function(x) {
  return(x$colors)
})
merge35 <- map2(counts_all2, colors, function(x,y) {
  return(mergeCloseModules(x, y, cutHeight = thresh35, verbose = 3))
})
colors35 <- lapply(merge35, function(x) {
  return(x$colors)
})


# visualizing new modules

# unmerged
for (i in 1:2) {
  plotDendroAndColors(geneTrees[[i]], colors00[[i]],
                      names(geneTrees)[i], main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}
# merged at 90% similarity
for (i in 1:2) {
  plotDendroAndColors(geneTrees[[i]], colors10[[i]],
                      names(geneTrees)[i], main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}
# merged at 75% similarity
for (i in 1:2) {
  plotDendroAndColors(geneTrees[[i]], colors25[[i]],
                      names(geneTrees)[i], main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}
# merged at 65% similarity
for (i in 1:2) {
  plotDendroAndColors(geneTrees[[i]], colors35[[i]],
                      names(geneTrees)[i], main = names(geneTrees)[i],
                      dendroLabels = FALSE)
} # ok 65% similarity really isn't necessary

# matching module color labels
# before matching
for (i in 1:2) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors00[[1]], colors00[[2]]),
                      groupLabels = names(colors), main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}
for (i in 1:2) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors10[[1]], colors10[[2]]),
                      groupLabels = names(colors), main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}
for (i in 1:2) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors25[[1]], colors25[[2]]),
                      groupLabels = names(colors), main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}
for (i in 1:2) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors35[[1]], colors35[[2]]),
                      groupLabels = names(colors), main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}

# matching module colors across networks
# matchLabels looks for modules with significantly overlapping gene content and makes them the same color (not as sophisticated as the permutation tests in modulePreservation, it's just a fisher exact test)
colorsUnmatched00 <- colors00
colorsUnmatched10 <- colors10
colorsUnmatched25 <- colors25
colorsUnmatched35 <- colors35
colors00 <- lapply(colors00, function(x) {
  return(matchLabels(source = x, reference = colors00$cer)) 
})
colors10 <- lapply(colors10, function(x) {
  return(matchLabels(source = x, reference = colors10$cer)) 
})
colors25 <- lapply(colors25, function(x) {
  return(matchLabels(source = x, reference = colors25$cer)) 
})
colors35 <- lapply(colors35, function(x) {
  return(matchLabels(source = x, reference = colors35$cer)) 
})

# after matching
for (i in 1:2) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors00[[1]], colors00[[2]]),
                      groupLabels = names(colors), main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}
for (i in 1:2) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors10[[1]], colors10[[2]]),
                      groupLabels = names(colors), main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}
for (i in 1:2) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors25[[1]], colors25[[2]]),
                      groupLabels = names(colors), main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}
for (i in 1:2) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors35[[1]], colors35[[2]]),
                      groupLabels = names(colors), main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}

# recalculate MEs post merge/color match
MEs00 <- map2(counts_all2, colors00, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y)
  MEs <- MEoutput$eigengenes
  return(MEs)
})
MEs10 <- map2(counts_all2, colors10, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y)
  MEs <- MEoutput$eigengenes
  return(MEs)
})
MEs25 <- map2(counts_all2, colors25, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y)
  MEs <- MEoutput$eigengenes
  return(MEs)
})
MEs35 <- map2(counts_all2, colors35, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y)
  MEs <- MEoutput$eigengenes
  return(MEs)
})

# final network visualizations
for (i in 1:2) {
  plotDendroAndColors(geneTrees[[i]], colors00[[i]],
                      names(geneTrees)[i], main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}
for (i in 1:2) {
  plotDendroAndColors(geneTrees[[i]], colors10[[i]],
                      names(geneTrees)[i], main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}
for (i in 1:2) {
  plotDendroAndColors(geneTrees[[i]], colors25[[i]],
                      names(geneTrees)[i], main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}
for (i in 1:2) {
  plotDendroAndColors(geneTrees[[i]], colors35[[i]],
                      names(geneTrees)[i], main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}

# We're done building our networks! Let's save all the relevant files
save(file = "data_files/Networks.RData", geneTrees, colors00, colors10, colors25, colors35)

################# Supplement: Leave one out ################# 

# How does the selection of environments present influence the module definitions?
# I predict it's similar to changing merge threshold, certain modules
# are merged with each other b/c they only differ in expression in an
# environment that isn't present in the dataset

# function to run a full network construction for given partitions of counts
# @input: unfiltered counts and counts_allele, 4 logical vectors for which
#         samples to keep in each of the 4 species/alleles
# @output: list of length 4 giving color definitions
# hardcoded parameters: 
# SoftPower = 5
# no merging of modules
# minModuleSize = 10
leaveOneOut <- function(.cts, .cts_allele, .keep_cer, .keep_par) {
  counts_all2 <- list(cer = t(.cts[, .keep_cer]), # unlike DESeq2, WGCNA expects ROWS to be samples and COLUMNS to be genes, likely because that's what R's cor function expects for meaningful pairwise correlations
                      par = t(.cts[, .keep_par]))

  # setting softPower based on above plots
  softPower <- 5
  
  # getting adjacency matrix from count matrix (this takes a bit)
  adjFromCounts <- function(.cts) {
    adj <- matrix(NA, nrow = ncol(.cts), ncol = ncol(.cts))
    for (i in c(1:(ncol(adj) - 1))) {
      cat(i, "/", ncol(adj), "\n")
      for (j in c((i+1):ncol(adj))) {
        # non_na <- cbind(.cts[,i], .cts[,j])[!is.na(.cts[,i]) & !is.na(.cts[,j]),]
        # dcor_ij <- dcor(non_na[,1], non_na[,2])
        # adj[i, j] <- dcor_ij
        # adj[j, i] <- dcor_ij
        cor_ij <- cor(.cts[,i], .cts[,j], use = "pairwise.complete.obs")
        adj[i, j] <- cor_ij
        adj[j, i] <- cor_ij
      }
    }
    diag(adj) <- 0
    return(adj)
  }
  
  # creating adjacency matrix and TOM
  adjs <- lapply(counts_all2, adjFromCounts)
  TOMs <- lapply(adjs, TOMsimilarity, TOMType = "signed")
  dissTOMs <- lapply(TOMs, function(x) {
    return(1-x)
  })
  geneTrees <- lapply(dissTOMs, function(x) {
    return(hclust(as.dist(x), method = "average"))
  })
  
  # module selection
  minModuleSize <- 10
  
  colors <- lapply(geneTrees, function(x) {
    chopped_tree <- cutreeDynamicTree(dendro = x,
                                      maxTreeHeight = 0.995,
                                      deepSplit = FALSE,
                                      minModuleSize = minModuleSize) %>% labels2colors()
    return(chopped_tree)
  })
  
  return(list(geneTrees = geneTrees, colors = colors))
}
# exclude LowN
keep_cer <- sample_info$organism == "cer" & 
  sample_info$genotype == "WT" &
  sample_info$experiment != "LowN"
keep_par <- sample_info$organism == "par" &
  sample_info$genotype == "WT" &
  sample_info$experiment != "LowN"
# constructing network
networks_noLowN <- leaveOneOut(.cts = counts, 
                               .cts_allele = counts_allele,
                               .keep_cer = keep_cer,
                               .keep_par = keep_par)
for (i in 1:2) {
  plotDendroAndColors(networks_noLowN$geneTrees[[i]], networks_noLowN$colors[[i]],
                      names(networks_noLowN$geneTrees)[i], main = paste(names(networks_noLowN$geneTrees)[i], "no LowN"),
                      dendroLabels = FALSE)
}
# exclude HAP4
keep_cer <- sample_info$organism == "cer" & 
  sample_info$genotype == "WT" &
  sample_info$experiment != "HAP4"
keep_par <- sample_info$organism == "par" &
  sample_info$genotype == "WT" &
  sample_info$experiment != "HAP4"

# constructing network
networks_noHAP4 <- leaveOneOut(.cts = counts, 
                               .cts_allele = counts_allele,
                               .keep_cer = keep_cer,
                               .keep_par = keep_par)
for (i in 1:2) {
  plotDendroAndColors(networks_noHAP4$geneTrees[[i]], networks_noHAP4$colors[[i]],
                      names(networks_noHAP4$geneTrees)[i], main = paste(names(networks_noHAP4$geneTrees)[i], "no HAP4"),
                      dendroLabels = FALSE)
}
# exclude CC
keep_cer <- sample_info$organism == "cer" & 
  sample_info$genotype == "WT" &
  sample_info$experiment != "CC"
keep_par <- sample_info$organism == "par" &
  sample_info$genotype == "WT" &
  sample_info$experiment != "CC"
# constructing network
networks_noCC <- leaveOneOut(.cts = counts, 
                               .cts_allele = counts_allele,
                               .keep_cer = keep_cer,
                               .keep_par = keep_par)
for (i in 1:2) {
  plotDendroAndColors(networks_noCC$geneTrees[[i]], networks_noCC$colors[[i]],
                      names(networks_noCC$geneTrees)[i], main = paste(names(networks_noCC$geneTrees)[i], "no CC"),
                      dendroLabels = FALSE)
}
# exclude LowPi
keep_cer <- sample_info$organism == "cer" & 
  sample_info$genotype == "WT" &
  sample_info$experiment != "LowPi"
keep_par <- sample_info$organism == "par" &
  sample_info$genotype == "WT" &
  sample_info$experiment != "LowPi"
# constructing network
networks_noLowPi <- leaveOneOut(.cts = counts, 
                               .cts_allele = counts_allele,
                               .keep_cer = keep_cer,
                               .keep_par = keep_par)
for (i in 1:2) {
  plotDendroAndColors(networks_noLowPi$geneTrees[[i]], networks_noLowPi$colors[[i]],
                      names(networks_noLowPi$geneTrees)[i], main = paste(names(networks_noLowPi$geneTrees)[i], "no LowPi"),
                      dendroLabels = FALSE)
}

# exclude Temp
keep_cer <- sample_info$organism == "cer" & 
  sample_info$genotype == "WT" &
  sample_info$experiment != "Temp"
keep_par <- sample_info$organism == "par" &
  sample_info$genotype == "WT" &
  sample_info$experiment != "Temp"
# constructing network
networks_noTemp <- leaveOneOut(.cts = counts, 
                                .cts_allele = counts_allele,
                                .keep_cer = keep_cer,
                                .keep_par = keep_par,
                                .keep_hyc = keep_hyc,
                                .keep_hyp = keep_hyp)
for (i in 1:2) {
  plotDendroAndColors(networks_noTemp$geneTrees[[i]], networks_noTemp$colors[[i]],
                      names(networks_noTemp$geneTrees)[i], main = paste(names(networks_noTemp$geneTrees)[i], "no Temp"),
                      dendroLabels = FALSE)
}

# saving
save(networks_noLowN, networks_noHAP4, networks_noCC, networks_noLowPi, networks_noTemp,
     file = "data_files/LeaveOneOut.RData")

########################### Archive #######################

# ##### Data exploration ####
# # see how well average expression level of each gene and overall connectivity of the inferred networks correlate between data sets
# # comparisons: each of the 4 experiments within one organism (ex: HAP4 vs LowPi in hyc), (4 choose 2) * 4 = 24 total comparisons
# correlateExpressionAndConnectivity <- function(cts_list, ntwks_titles) {
#   n_ntwks <- length(cts_list)
#   n_comparisons <- choose(n_ntwks, 2)
#   plotList <- vector(mode = "list", length = n_comparisons*2)
#   plotCounter <- 1
#   
#   for (i in 1:(n_ntwks-1)) {
#     for (j in (i+1):n_ntwks) {
#       rankExpr1 <- rank(rowMeans(t(cts_list[[i]])))
#       rankExpr2 <- rank(rowMeans(t(cts_list[[j]])))
#       rankConn1 <- rank(softConnectivity(cts_list[[i]], type="unsigned", power=softPower))
#       rankConn2 <- rank(softConnectivity(cts_list[[j]], type="unsigned", power=softPower))
#       plotTibble <- tibble("rankExpr1" = rankExpr1,
#                            "rankExpr2" = rankExpr2,
#                            "rankConn1" = rankConn1,
#                            "rankConn2" = rankConn2)
#       exprCorr <- cor(rankExpr1, rankExpr2, method = "pearson")
#       exprPvalue <- cor.test(rankExpr1, rankExpr2, method = "pearson")$p.value
#       connCorr <- cor(rankConn1, rankConn2, method = "pearson")
#       connPvalue <- cor.test(rankConn1, rankConn2, method = "pearson")$p.value
#       p_expr <- ggplot(data = plotTibble, aes(x = rankExpr1, y = rankExpr2)) + 
#         xlab(paste("Ranked Expression", ntwks_titles[[i]])) + 
#         ylab(paste("Ranked Expression", ntwks_titles[[j]])) + 
#         ggtitle(paste("Correlation = ", exprCorr, "\n", "p-value = ", exprPvalue)) +
#         geom_point()
#       p_conn <- ggplot(data = plotTibble, aes(x = rankConn1, y = rankConn2)) + 
#         xlab(paste("Ranked Connectivity", ntwks_titles[[i]])) + 
#         ylab(paste("Ranked Conenctivity", ntwks_titles[[j]])) + 
#         ggtitle(paste("Correlation = ", connCorr, "\n", "p-value = ", connPvalue)) +
#         geom_point()
#       plotList[[plotCounter]] <- p_expr
#       plotCounter <- plotCounter + 1
#       plotList[[plotCounter]] <- p_conn
#       plotCounter <- plotCounter + 1
#     }
#   }
#   return(plotList)
# }
# 
# # not good for QC, as there is real evolution captured here when comparing between species,
# # but still interesting
# ExprConn <- correlateExpressionAndConnectivity(cts_list = list(counts_top3000$cer,
#                                                                counts_top3000$par,
#                                                                counts_top3000$hyc,
#                                                                counts_top3000$hyp), # anna... why didn't you just give it counts_top3000? That's a list... Is this the dumbest thing I've ever done or is this the dumbest function that's ever existed
#                                                ntwks_titles = names(counts_top3000))
# library(gridExtra)
# grid.arrange(grobs = ExprConn, nrow = 3) # it checks out! Expression correlation is high across the board, but highest for cer-cer/par-par comparisons. Connectivity is more variable but highest for cer-cer/par-par
# 
# #### TOMplots ####
# # oh my goodness what a ride it has been
# # this is far from the ideal heatmap/dendrogram but at least it allows for multiple colors
# TOMplotSubsetMultiColor <- function(dissTOM, colorList, sampleSize = 1000, plotTitle = "TOM plot") {
#   nGenes <- ncol(dissTOM)
#   randomSample <- sample(nGenes, sampleSize)
#   selectTOM <- dissTOM[randomSample, randomSample]
#   plotTOM <- selectTOM^7 # I guess 7 is just a power they thought looked pretty
#   diag(plotTOM) <- NA
#   plotTree <-  hclust(as.dist(plotTOM), method = "average")
#   layout_matrix <- cbind(c(2,3,0), c(0,0,5), c(6,1,4))
#   layout(layout_matrix, widths = c(1,1,4), heights = c(1,2,4))
#   colorMat <- NULL
#   for (i in 1:length(colorList)) {
#     colorMat <- cbind(colorMat, colorList[[i]][randomSample])
#   }
#   plotOrderedColors(plotTree$order, colorMat, rowLabels = NA, align = "edge")
#   TOMplot(-plotTOM, plotTree, Colors = NULL, main = plotTitle, setLayout = FALSE, margins = c(2,2))
# }
# 
# # TOMplots just for RStudio
# TOMplotTitles <- c("S. cerevisiae", "S. paradoxus", "F1 hybrid - summed alleles", "F1 hybrid - cerevisiae allele", "F1 hybrid - paradoxus allele")
# par(mar = c(2,2,2,2)) # don't ask me why you need this
# 
# # sacrificial plot --- for some delicious reason, the first plot has its module definitions lined up incorrectly, but it's fine when you re-run it
# TOMplotSubsetMultiColor(dissTOM = dissTOMs[[1]], 
#                        colorList = colors, plotTitle = TOMplotTitles[1], sampleSize = 100) 
# 
# # for all4
# for(i in 1:4) {
#   TOMplotSubsetMultiColor(dissTOM = dissTOMs[[i]], 
#                           colorList = colors, plotTitle = TOMplotTitles[i], sampleSize = 1000)
# } 
# 
# # for just cer and par
# for(i in 1:2) {
#   TOMplotSubsetMultiColor(dissTOM = dissTOMs[[i]], 
#                 colorList = list(colors[[1]], colors[[2]]), plotTitle = TOMplotTitles[i], sampleSize = 1000)
# } 
# 
# # believe it or not, this is all for the legend
# library(circlize)
# library(ComplexHeatmap)
# colorFunction <- colorRamp2(breaks = seq(from = 0, to = 1, length.out = 12), colors = hcl.colors(12, "YlOrRd", rev = TRUE))
# lgd <- Legend(col_fun = colorFunction, title = "Topological Overlap")
# pushViewport(viewport(width = 0.9, height = 0.9))
# grid.rect()  # border
# draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
# popViewport()
# 
# # it has come to my attention that the default WGCNA colors are all in RGB and look pretty muted in CMYK
# # sooo for posters and printing, here's how you can convert colors if you so choose
# # library(scales)
# # AnnasHandpickedPalette <- c("aquamarine", "tomato", "yellow", "seagreen", "darkorchid", "deepskyblue","lawngreen", "midnightblue", "deeppink", "thistle", "saddlebrown")
# # show_col(AnnasHandpickedPalette)
# # ColorConversionTable <- rbind(c("grey", "grey"), cbind(setdiff(unique(c(colors$cer, colors$par, colors$hyc, colors$hyp)), "grey"), AnnasHandpickedPalette))
# # changeColor <- function(oldColors, conversionTable) {
# #   newColorIdx <- which(conversionTable[, 1] == oldColors)
# #   return(conversionTable[newColorIdx, 2])
# # }
# # newColors <- lapply(colors, function(x) {
# #   return(sapply(x, changeColor, conversionTable = ColorConversionTable))
# # })
# 
# # outputting TOMplots (Note this DOES NOT WORK with the multicolor! It misaligns the 4 color rows)
# TOMplotFileNames <- c("TOMplot_cer.pdf", "TOMplot_par.pdf", "TOMplot_hyc.pdf", "TOMplot_hyp.pdf")
# for (i in 1:4) {
#   pdf(file = paste0("../../../posters/", TOMplotFileNames[i]), colormodel = "cmyk")
#   TOMplotSubset(dissTOM = dissTOMs[[i]], 
#                           colorList = colors, plotTitle = TOMplotTitles[i], sampleSize = 1000)
#   dev.off()
# } 
# 
# # We're done building our networks! Let's save all the relevant files
# save(file = "data_files/Networks_mergemodules25.RData", counts_all5, counts_top3000, adjs, TOMs, dissTOMs,
#      geneTrees, MEs, TOMs, colors, dissTOMs, softPower)




# ############### Distance Correlation Version ###############
# # following Hou et al. 2021 K-module algorithm's suggestion
# # Distance Correlation originally described in Szekely et al. 2007
# # very helpful article comparing Distance Correlation to Pearson: https://blogs.sas.com/content/iml/2018/04/04/distance-correlation.html
# library(energy)
# test_counts <- counts_all5$cer %>% t()
# # TDH3 and TDH2, strongly correlated
# plot(test_counts["YGR192C", ], test_counts["YJR009C", ])
# dcor(test_counts["YGR192C", ], test_counts["YJR009C", ])
# cor(test_counts["YGR192C", ], test_counts["YJR009C", ])
# # random pair, shouldn't be correlated
# gene_idxs <- sample(c(1:nrow(test_counts)), 2, replace = FALSE)
# plot(test_counts[gene_idxs[1], ], test_counts[gene_idxs[2], ])
# dcor(test_counts[gene_idxs[1], ], test_counts[gene_idxs[2], ])
# cor(test_counts[gene_idxs[1], ], test_counts[gene_idxs[2], ])
# 
# # Creating similarity matrix
# test_list <- t(test_counts) %>% as.data.frame() %>% as.list()
# dcor_on_lists <- function(.x, .y) {
#   # defining helper function inside other function cause I feel fancy
#   helper <- function(.a, .b, .a_name, .b_name) {
#     cat("currently correlating", .a_name, "and", .b_name, "\n")
#     if (length(.a) != length(.b)) {
#       stop("vectors", .a_name, "and", .b_name, "are not the same length")
#     }
#     return(dcor(.a, .b)) # TODO: this is horrendously slow. There must be a trick. First of all not calculating every value twice as it's symmetric...
#   }
#   dcor_vec <- map(c(1:length(.x)), \(i) {
#     output <- helper(.x[[i]], .y[[i]], names(.x)[i], names(.y)[i])
#     return(output)
#   })
#   return(dcor_vec)
# }
# test_sim <- outer(test_list, test_list, FUN = dcor_on_lists) # waaaay too slow
# # 
# # parallel version
# library(doParallel)
# registerDoParallel(cores = 4)
# getDoParWorkers()
# test_dcor <- foreach(i = test_list, .combine = "cbind", .packages = "energy") %:%
#   foreach(j = test_list, .combine = c, .packages = "energy") %dopar% { # TODO: change j to always be the rows below i (j > which(as.list(as.data.frame(test_counts))) == i or something like that should work I'm guessing)
#     dcor(i, j)
#   }
# 
# # function taken from dyerlab github (Rodney J Dyer)
# cov2dist <- function( C ) {
#   if( dim(C)[1] != dim(C)[2] )
#     stop("Cannot use non-symmetric matrices for this...")
#   K <- dim(C)[1]
#   
#   D <- matrix(0,nrow=K,ncol=K)
#   for( i in 1:K){
#     for( j in (i+1):K) {
#       if( j <= K )
#         D[i,j] <- C[i,i] + C[j,j] - 2.0*C[i,j]
#     }
#   }
#   return( D + t(D) )
# }
# test_cov <- cov(t(test_counts))
# test_dist <- test_counts[1,] %>% calc_dist()
# 
# AnnasDistanceCor <- function(.x, .y) {
#   nObs <- length(.x)
#   if (length(.y) != nObs) {
#     stop("my version only works with vectors of the same length")
#   }
#   A <- 0
#   for (i in 1:nObs) {
#     for (j in 1:nObs) {
#       A <- A + abs(.x[i] - .x[j])*abs(.y[i] - .y[j])
#     }
#   }
#   A <- A/(nObs^2)
#   return(A)
# }
# random_idx1 <- sample(c(1:nrow(test_counts)), 1)
# random_idx2 <- sample(c(1:nrow(test_counts)), 1)
# # cov2dist version:
# test_dist[random_idx1, random_idx2]
# # my understanding of SAS article version:
# AnnasDistanceCor(test_counts[random_idx1,], test_counts[random_idx2,])
# 
# 
################# Pearson Correlation Signed Version ############ 
# Archived b/c this should be exactly what we do to construct main networks now, but kept for reference in case it isn't
# # Parameter Fitting 
# # pick power parameter ß for determining soft threshold of edge presence in our network
# powers <- c(1:20)
# sfts <- lapply(counts_all4, pickSoftThreshold, powerVector = powers, 
#                verbose = 5, networkType = "signed") # It really just calculates R^2 for different log(k) vs log(p(k)) (aka the classic scale free network) plots
# 
# # degree distribution plots using different ß-based soft thresholds
# 
# # wrapping this WGCNA code in a function for better repeatability
# plotSoftPower <- function(sft, sft_name) {
#   # plotting results to pick power
#   par(mfrow=c(1,2))
#   cex1 = 0.9
#   # Scale-free topology fit index vs soft-thresholding power
#   plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], # I think this is just to make sure that the slope is negative? Otherwise you could just use the 2nd column without the 3rd
#        xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#        main = paste("Scale independence\n", sft_name))
#   text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#        labels=powers,cex=cex1,col="red")
#   abline(h=0.9, col="red") # this line corresponds to using an R^2 cut-off of h apparently
#   # Mean connectivity vs soft-thresholding power
#   plot(sft$fitIndices[,1], sft$fitIndices[,5],
#        xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#        main = paste("Mean connectivity"))
#   text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# }
# map2(sfts, names(sfts), function(x, y) {
#   return(plotSoftPower(x, y))
# })
# 
# # setting softPower based on above plots
# softPower <- 14
# 
# # Building networks
# # Module selection: Adjacency matrix, topical overlap (TOM), and heirachical clustering
# 
# # set soft power (decided in parameter fitting section)
# softPower <- 14
# 
# # getting adjacency matrix from count matrix (this takes a bit)
# adjFromCounts <- function(cts) {
#   adj <- adjacency(cts, type = "signed", power = softPower)
#   diag(adj) <- 0
#   return(adj)
# }
# 
# # creating adjacency matrix and TOM
# adjs <- lapply(counts_all4, adjFromCounts)
# TOMs <- lapply(adjs, TOMsimilarity, TOMType = "signed")
# dissTOMs <- lapply(TOMs, function(x) {
#   return(1-x)
# })
# geneTrees <- lapply(dissTOMs, function(x) {
#   return(hclust(as.dist(x), method = "average"))
# })
# 
# par(mfrow = c(2,3))
# for (i in 1:4) {
#   plot(geneTrees[[i]], xlab="", sub="", main=names(counts_all4)[i],
#        labels = FALSE, hang=0.04)
# } 
# 
# # module selection
# minModuleSize <- 30
# 
# colors <- lapply(geneTrees, function(x) {
#   chopped_tree <- cutreeDynamicTree(dendro = x,
#                                     maxTreeHeight = 0.995,
#                                     deepSplit = FALSE,
#                                     minModuleSize = minModuleSize) %>% labels2colors()
#   return(chopped_tree)
# })
# for (i in 1:4) {
#   plotDendroAndColors(dendro = geneTrees[[i]], colors = colors[[i]],
#                       main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# 
# # merging close modules
# # calculate eigengenes (eigenvector with largest eigenvalue aka 1st PC of expression matrix (woah way back to that thing) for each module)
# # Note: colors is not an aesthetic argument, it's literally the only information on module membership
# # moduleEigengenes for a p x n module X is equivalent to svd(t(scale(x)))$v[,1], the first principle component (length p) of X. If align = "along average", the sign is also potentially flipped so that the correlation between the average gene expression of each module across environments and its eigengene is non-negative (cor(A.E., M.E.) >= 0)
# MEs <- map2(counts_all4, colors, function(x, y) {
#   MEoutput <- moduleEigengenes(x, colors = y)
#   MEs <- MEoutput$eigengenes
#   return(MEs)
# })
# MEDiss <- lapply(MEs, function(x) {
#   return(1-cor(x))
# })
# METree <- lapply(MEDiss, function(x) {
#   return(hclust(as.dist(x), method = "average"))
# })
# 
# # merging close modules
# # using multiple threshold to test how robust conclusions are to module definition
# thresh35 <- 0.35
# thresh25 <- 0.25 
# thresh10 <- 0.10
# # thresh00 <- 0.00 # not a real threshold, just to point out that the fourth condition is no merging
# 
# # visualizing potential new modules
# # .35 thresh, merging at 65% similarity
# par(mfrow = c(2,3))
# for (i in 1:4) {
#   plot(METree[[i]], main = names(counts_all4)[i], xlab = "", sub = "")
#   abline(h=thresh35, col = "red")
# } # grey won't be merged by mergeCloseModules
# 
# # .25 thresh, merging at 75% similarity
# dev.off()
# par(mfrow = c(2,3))
# for (i in 1:4) {
#   plot(METree[[i]], main = names(counts_all4)[i], xlab = "", sub = "")
#   abline(h=thresh25, col = "red")
# }
# 
# # .10 thresh, merging at 90% similarity (extreme)
# dev.off()
# par(mfrow = c(2,3))
# for (i in 1:4) {
#   plot(METree[[i]], main = names(counts_all4)[i], xlab = "", sub = "")
#   abline(h=thresh10, col = "red")
# }
# 
# # actually merging
# colors00 <- colors
# unmergedMEs <- MEs
# merge10 <- map2(counts_all4, colors, function(x,y) {
#   return(mergeCloseModules(x, y, cutHeight = thresh10, verbose = 3))
# })
# colors10 <- lapply(merge10, function(x) {
#   return(x$colors)
# })
# merge25 <- map2(counts_all4, colors, function(x,y) {
#   return(mergeCloseModules(x, y, cutHeight = thresh25, verbose = 3))
# })
# colors25 <- lapply(merge25, function(x) {
#   return(x$colors)
# })
# merge35 <- map2(counts_all4, colors, function(x,y) {
#   return(mergeCloseModules(x, y, cutHeight = thresh35, verbose = 3))
# })
# colors35 <- lapply(merge35, function(x) {
#   return(x$colors)
# })
# 
# 
# # visualizing new modules
# 
# # unmerged
# for (i in 1:2) {
#   plotDendroAndColors(geneTrees[[i]], colors00[[i]],
#                       names(geneTrees)[i], main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# # merged at 90% similarity
# for (i in 1:2) {
#   plotDendroAndColors(geneTrees[[i]], colors10[[i]],
#                       names(geneTrees)[i], main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# # merged at 75% similarity
# for (i in 1:2) {
#   plotDendroAndColors(geneTrees[[i]], colors25[[i]],
#                       names(geneTrees)[i], main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# # merged at 65% similarity
# for (i in 1:2) {
#   plotDendroAndColors(geneTrees[[i]], colors35[[i]],
#                       names(geneTrees)[i], main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# 
# # matching module color labels
# # before matching
# for (i in 1:2) {
#   plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors00[[1]], colors00[[2]], colors00[[3]], colors00[[4]]),
#                       groupLabels = names(colors), main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# for (i in 1:2) {
#   plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors10[[1]], colors10[[2]], colors10[[3]], colors10[[4]]),
#                       groupLabels = names(colors), main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# for (i in 1:2) {
#   plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors25[[1]], colors25[[2]], colors25[[3]], colors25[[4]]),
#                       groupLabels = names(colors), main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# 
# # matching module colors across networks
# # matchLabels looks for modules with significantly overlapping gene content and makes them the same color (not as sophisticated as the permutation tests in modulePreservation, it's just a fisher exact test)
# colorsUnmatched00 <- colors00
# colorsUnmatched10 <- colors10
# colorsUnmatched25 <- colors25
# colorsUnmatched35 <- colors35
# colors00 <- lapply(colors00, function(x) {
#   return(matchLabels(source = x, reference = colors00$cer)) 
# })
# colors10 <- lapply(colors10, function(x) {
#   return(matchLabels(source = x, reference = colors10$cer)) 
# })
# colors25 <- lapply(colors25, function(x) {
#   return(matchLabels(source = x, reference = colors25$cer)) 
# })
# colors35 <- lapply(colors35, function(x) {
#   return(matchLabels(source = x, reference = colors35$cer)) 
# })
# 
# # after matching
# for (i in 1:2) {
#   plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors00[[1]], colors00[[2]]),
#                       groupLabels = names(colors), main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# for (i in 1:2) {
#   plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors10[[1]], colors10[[2]]),
#                       groupLabels = names(colors), main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# for (i in 1:2) {
#   plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors25[[1]], colors25[[2]]),
#                       groupLabels = names(colors), main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# for (i in 1:2) {
#   plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors35[[1]], colors35[[2]]),
#                       groupLabels = names(colors), main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# 
# # recalculate MEs post merge/color match
# MEs00 <- map2(counts_all4, colors00, function(x, y) {
#   MEoutput <- moduleEigengenes(x, colors = y)
#   MEs <- MEoutput$eigengenes
#   return(MEs)
# })
# MEs10 <- map2(counts_all4, colors10, function(x, y) {
#   MEoutput <- moduleEigengenes(x, colors = y)
#   MEs <- MEoutput$eigengenes
#   return(MEs)
# })
# MEs25 <- map2(counts_all4, colors25, function(x, y) {
#   MEoutput <- moduleEigengenes(x, colors = y)
#   MEs <- MEoutput$eigengenes
#   return(MEs)
# })
# MEs35 <- map2(counts_all4, colors35, function(x, y) {
#   MEoutput <- moduleEigengenes(x, colors = y)
#   MEs <- MEoutput$eigengenes
#   return(MEs)
# })
# 
# # final network visualizations
# for (i in 1:2) {
#   plotDendroAndColors(geneTrees[[i]], colors00[[i]],
#                       names(geneTrees)[i], main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# for (i in 1:2) {
#   plotDendroAndColors(geneTrees[[i]], colors10[[i]],
#                       names(geneTrees)[i], main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# for (i in 1:2) {
#   plotDendroAndColors(geneTrees[[i]], colors25[[i]],
#                       names(geneTrees)[i], main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# for (i in 1:2) {
#   plotDendroAndColors(geneTrees[[i]], colors35[[i]],
#                       names(geneTrees)[i], main = names(geneTrees)[i],
#                       dendroLabels = FALSE)
# }
# 
# # We're done building our networks! Let's save all the relevant files
# save(file = "data_files/NetworksSigned.RData", geneTrees, colors00, colors10, colors25, colors35)

# #### Parameter Fitting ####
# # pick power parameter ß for determining soft threshold of edge presence in our network
# powers <- c(1:20)
# sfts <- lapply(counts_top3000, pickSoftThreshold, powerVector = powers, verbose = 5) # It really just calculates R^2 for different log(k) vs log(p(k)) (aka the classic scale free network) plots
# 
# # degree distribution plots using different ß-based soft thresholds
# 
# # wrapping this WGCNA code in a function for better repeatability
# plotSoftPower <- function(sft, sft_name) {
#   # plotting results to pick power
#   par(mfrow=c(1,2))
#   cex1 = 0.9
#   # Scale-free topology fit index vs soft-thresholding power
#   plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], # I think this is just to make sure that the slope is negative? Otherwise you could just use the 2nd column without the 3rd
#        xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#        main = paste("Scale independence\n", sft_name))
#   text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#        labels=powers,cex=cex1,col="red")
#   abline(h=0.9, col="red") # this line corresponds to using an R^2 cut-off of h apparently
#   # Mean connectivity vs soft-thresholding power
#   plot(sft$fitIndices[,1], sft$fitIndices[,5],
#        xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#        main = paste("Mean connectivity"))
#   text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# }
# map2(sfts, names(sfts), function(x, y) {
#   return(plotSoftPower(x, y))
# })
# 
# # 6 looks good as it turns out
# softPower <- 6
# 
# # getting adjacency matrix from count matrix (this takes a bit)
# adjFromCounts <- function(cts) {
#   adj <- adjacency(cts, type = "unsigned", power = softPower)
#   diag(adj) <- 0
#   return(adj)
# }
# 
# # creating adjacency matrix and TOM
# adjs <- lapply(counts_top3000, adjFromCounts)
# TOMs <- lapply(adjs, TOMsimilarity, TOMType = "unsigned")
# dissTOMs <- lapply(TOMs, function(x) {
#   return(1-x)
# })
# geneTrees <- lapply(dissTOMs, function(x) {
#   return(hclust(as.dist(x), method = "average"))
# })
# 
# par(mfrow = c(2,3))
# for (i in 1:5) {
#   plot(geneTrees[[i]], xlab="", sub="", main=names(counts_all5)[i],
#        labels = FALSE, hang=0.04)
# } 
# 

#### Exploratory Section: Alternative ways of creating similarity matrix ####

# Rationale: Based on feedback at 2023 Gordon Conference, a Pearson correlation
# of gene expression across conditions is going to tend to be fairly high if there
# are many conditions in the dataset where most genes aren't changing expression all that much (there really are)
# So instead of correlating each gene pair's expression vectors, correlate *some sort* of eigenvalue decomposition
# (so like correlate the vector orthogonal to the gene expression, the principle eigenvector)
# Issue: How do you get one principle eigenvector for each gene? Is this even the right question to be asking?

# Archived PCA dimension reduction, as it doesn't seem to work
# for clustering groups of genes by expression across environments 
# (only works for clustering cell types by which genes they express)
# but might be interesting to explore exactly why this is the case someday
# # testing out PCA
# test_YtY <- t(counts_all5$cer) %*% counts_all5$cer
# test_eigen <- eigen(test_YtY)
# 
# # first getting colors for standard, non-reduced WGCNA as our "goldish standard"
# nonred_counts <- counts_all5$cer
# # visualizing network connectivity
# softPower <- 6
# softConn <- softConnectivity(nonred_counts, type = "unsigned", power = softPower)
# hist(softConn, breaks = 50)
# 
# # soft power parameter fitting
# powers <- c(1:20)
# sfts <- pickSoftThreshold(nonred_counts, powerVector = powers, verbose = 5)
# plotSoftPower(sfts, "cer not reduced")
# # picking soft power
# softPower <- 6
# 
# # adjacency and TO
# adj <- adjFromCounts(nonred_counts)
# TOM <- TOMsimilarity(adj, TOMType = "unsigned")
# dissTOM <- 1 - TOM
# geneTree <- dissTOM %>% as.dist() %>% hclust(method = "average")
# plot(geneTree, xlab="", sub="", main="cer not reduced",
#      labels = FALSE, hang=0.04)
# 
# # module definitions
# nonred_colors <- cutreeDynamicTree(dendro = geneTree,
#                                    maxTreeHeight = 0.995,
#                                    deepSplit = FALSE,
#                                    minModuleSize = 15) %>% labels2colors()
# 
# plotDendroAndColors(dendro = geneTree, colors = nonred_colors, dendroLabels = FALSE)
# 
# # plotting Principal Components to see how well they recapitulate module colors
# pcdf <- bind_cols(test_eigen$vectors[, c(1:14)], nonred_colors)
# names(pcdf) <- c(paste0("pc", c(1:14)), "colors")
# for (i in c(1:13)) {
#   plot(x = pcdf[,i,drop=TRUE], y = pcdf[,i + 1,drop=TRUE], col = nonred_colors, pch = 16,
#        xlab = paste0("PC", i), ylab = paste0("PC", i+1))
# } # There seems like there's something going on there... it's not just a multicolor blob, different PCs coax out different modules (to some extent, it's not crazy dramatic)
# 
# # Very exploratory: How does kmeans clustering do versus WGCNA?
# kmeans_output <- kmeans(pcdf[,1:14], centers = 6)
# pcdf$kmeans <- kmeans_output$cluster %>% labels2colors()
# pcdf$kmeans%>% table() # seems like no matter how many centers I specify, over 4000 genes always end up in one module. This is consistent with the PCs individually only separating a few genes from 4000 central genes
# for (i in c(1:13)) {
#   plot(x = pcdf[,i,drop=TRUE], y = pcdf[,i + 1,drop=TRUE], col = labels2colors(pcdf$kmeans), pch = 16,
#        xlab = paste0("PC", i), ylab = paste0("PC", i+1))
# } 
# # kmeans as it turns out works great with PCs. The problem is the huge blob around 0,0 when plotting any two PCs
# 
# # If you ever feel like you've reduced dimensionality in a useful way, you can proceed:
# 
# # believe it or not, the values are literally the proportion of variance explained by each component
# # (because the sum of all the eigen values is the total variance equal to sum of the diagonal in covariance matrix!)
# plot(y = test_eigen$values[1:100]/test_eigen$values[1], x = c(1:100), ylab = "% variance explained", xlab = "PC")
# # looks like 16 components captures the vast majority of the variability, but we can test with any number
# test_red <- test_eigen$vectors[, c(1:3)] # eigen() outputs a matrix where columns are eigenvectors
# 
# # TODO: troubleshoot by finding a pair of genes that are very related (maybe a couple of the sporulation genes?)
# 
# # creating similarity matrix - pairwise L2 norms of each gene's vector of PC1-16
# # TODO: a distance matrix of 16 PCs still weighs those 16 equally. Should I scale them by eigenvalue or something?
# red_dist <- dist(test_red) %>% as.matrix()
# red_norm <- red_dist/max(red_dist) # normalize to 0-1 scale
# red_sim <- 1 - red_norm # the smaller the distance the more similar
# 
# # visualizing network connectivity
# softConn <- adjacency.fromSimilarity(red_sim, power = 6) %>% colSums() # play around with different powers to see how it affects degree distribution
# hist(softConn, breaks = 50)
# 
# # soft power parameter fitting
# powers <- c(1:20)
# red_sfts <- pickSoftThreshold.fromSimilarity(red_sim, powerVector = powers, verbose = 5)
# plotSoftPower(red_sfts, "cer reduced")
# # picking soft power
# softPower <- 6
# 
# # adjacency and TO
# red_adj <- adjacency.fromSimilarity(red_sim, power = softPower)
# red_TOM <- TOMsimilarity(red_adj, TOMType = "unsigned")
# red_dissTOM <- 1 - red_TOM
# red_geneTree <- red_dissTOM %>% as.dist() %>% hclust(method = "average")
# plot(red_geneTree, xlab="", sub="", main="cer reduced",
#      labels = FALSE, hang=0.04)
# 
# # TODO: hmmm. Something is still weird about this. The individual PCs seem to be
# # picking up on some real aspects of clustering, but that's not reflected in the
# # gene trees that are just continuously increasing/decreasing distances with each gene
# # a few options to try:
# # 1) Read https://365datascience.com/tutorials/python-tutorials/pca-k-means/ to see
# #    if it has any ideas for using PCs as "orthogonal covariates" for clustering
# # 2) Make sure the L2 Norm isn't inappropriate to use on coordinates from 16 eigenvectors
# #    that each explain different proportions of variance. Like maybe they should be scaled
# #    by their eigenvalue?
# # 3) If all else fails, abandon PCA and try distance correlation instead
