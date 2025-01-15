#### Network Probing via Iterative Gene Swaps ####
# in which we identify the genes most responsible for differences in network structure between cerevisiae and paradoxus by iteratively altering the expression vectors of each species
sapply(c("tidyverse", "edgeR", "WGCNA", "matrixStats"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Barkai_data_analysis")
options(stringsAsFactors = FALSE)

# our starting networks and count data (counts_top3000)
load("data_files/Networks.RData")

# swapHighestDegree relies on genes vectors having the same length, not necessarily the case in our data
# So we need to filter using condition and well_flask_ID to get a standard set of samples comparable between species
load("data_files/Cleaned_Barkai_Data.RData") # just need this to know what samples to keep
cer_samples <- filter(sample_info, organism == "cer")
par_samples <- filter(sample_info, organism == "par")
paired_samples <- filter(sample_info, paste(condition, well_flask_ID) %in% paste(cer_samples$condition, cer_samples$well_flask_ID) &
                           paste(condition, well_flask_ID) %in% paste(par_samples$condition, par_samples$well_flask_ID) &
                           organism %in% c("cer", "par"))
table(paste(paired_samples$condition, paired_samples$well_flask_ID)) %>% table()
counts_cer <- counts_top3000$cer[paired_samples[paired_samples$organism == "cer",]$sample_name,]
counts_par <- counts_top3000$par[paired_samples[paired_samples$organism == "par",]$sample_name,]

############################ Function Definitions ############################ 

# Comparing Networks

# Calculate percent of genes in same module between species
# @input: two species color vectors of module labels, both length nGenes
# @output: a percent (double) of genes with same label
calculateModuleOverlap <- function(c1, c2) {
  c1 <- matchLabels(c1, c2)
  return(sum(c1 == c2)/length(c1))
}

# Calculate the absolute value of the differences in adjacencies at each edge in two networks then return the average difference
# @input: Two adjacency matrices to compare, A1 and A2
# @output: A double giving the average value of matrix abs(A1-A2)
calculateAvgAdjacencyDiff <- function(A1, A2) {
  return(mean(abs(A1-A2)))
}

# Choosing genes to swap

# Rank every gene in order of highest to lowest degree (highest degree gene gets the value 1)
# @input: adjacency matrix of the acceptor species
# @ouput: ranking of each gene (vector of length nGenes with values 1:nGenes corresponding to each gene's rank)
rankHighestDegree <- function(A) {
  degrees <- rowSums(A)
  return(order(degrees, decreasing = TRUE))
}

# Swaps the gene with the highest ranking that hasn't been swapped yet from acceptor's expression vector to donor's expression vector
# @input: ranking of each gene (vector of length nGenes with values 1:nGenes corresponding to each gene's rank), boolean vector of length nGenes that indicates if a gene has already been swapped, donor and acceptor's gene count matrices
# @ouput: acceptor's gene count matrix with one gene swapped to donor's vector, updated swap vector with gene that was just swapped set to TRUE
swapHighestDegree <- function(r, s, gcma, gcmd) {
  swap_idx <- which(r == min(r[!s]))
  new_gcm <- gcma
  new_gcm[, swap_idx] <- gcmd[, swap_idx]
  s[swap_idx] <- TRUE
  return(list(updated_gcm = new_gcm, updated_swap = s))
}

swapRandomGene <- function(s, gcma, gcmd) {
  swap_idx <- which(!s) %>% sample(size = 1)
  new_gcm <- gcma
  new_gcm[, swap_idx] <- gcmd[, swap_idx]
  s[swap_idx] <- TRUE
  return(list(updated_gcm = new_gcm, updated_swap = s))
}

# Re-Constructing Networks

# creates adjacency matrix from gene-count matrix with power 6 and diagonal set to 0
adjFromCounts <- function(gcm) {
  adj <- adjacency(gcm, type = "unsigned", power = 6)
  diag(adj) <- 0
  return(adj)
}

# given a gene count matrix, recalculates the adjacency matrix, topological overlap matrix, heiarchichal clustering, and tree cutting to define new modules
# @input: a gene count matrix
# @output: a new vector of module names of length nGenes corresponding to each gene's new module
recalculateModules <- function(gcm) {
  adj <- adjacency(gcm, power = 6) # I think for module definitions we don't want the diagonal to be set to 0
  TOM <- TOMsimilarity(adj, TOMType = "unsigned")
  dissTOM <- 1-TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  clrs <- cutreeDynamicTree(dendro = geneTree,
                            maxTreeHeight = 0.995,
                            deepSplit = FALSE,
                            minModuleSize = 15) %>% labels2colors()
  MEoutput <- moduleEigengenes(gcm, colors = clrs)
  MEs <- MEoutput$eigengenes
  MEDiss <- 1-cor(MEs)

  METree <- hclust(as.dist(MEDiss), method = "average")
  merge <- mergeCloseModules(gcm, clrs, cutHeight = 0.1, verbose = 3)
  clrs <- merge$colors
  return(clrs)
}

# Recursive version of Iteration (doesn't work seems to overload memory)
# 
# @input:
# TODO: right now it just outputs the vector of adj differences each iteration. For convenience we probably also want to know what genes (although ranking should tell us)
# @output: 
# swapGenesWholeNetwork <- function(gcma, gcmd, adjd, diffs = vector(mode = "double", length = 0), s, r, cutoff) {
#   adja <- adjacency(gcma, power = 6)
#   diff <- calculateAvgAdjacencyDiff(adja, adjd)
#   cat("diff = ", diff, "\n")
#   if (cutoff > diff) {
#     return(diffs)
#   }
#   diffs <- c(diffs, diff)
#   swap_results <- swapHighestDegree(r, s, gcma, gcmd)
#   gcma <- swap_results$updated_gcm
#   s <- swap_results$updated_swap
#   return(swapGenesWholeNetwork(gcma, gcmd, adjd, diffs, s, r, cutoff))
# }

################################## Iteration ################################## 

# module overlap check, random gene swap
gcma <- counts_cer
gcmd <- counts_par
adja <- adjFromCounts(gcma)
adjd <- adjFromCounts(gcmd)
clrsa <- recalculateModules(gcma)
clrsd <- recalculateModules(gcmd)
overlaps <- vector(mode = "double", length = 0)
OverlapCutoff <- 0.75
swapped <- rep(FALSE, ncol(adjd))
overlap <- calculateModuleOverlap(clrsa, clrsd)
overlaps <- c(overlaps, overlap)
while(overlap < OverlapCutoff & sum(swapped) < length(swapped)) {
  swap_results <- swapRandomGene(swapped, gcma, gcmd)
  swapped <- swap_results$updated_swap
  gcma <- swap_results$updated_gcm
  clrsa <- recalculateModules(gcma)
  overlap <- calculateModuleOverlap(clrsa, clrsd)
  cat("overlap = ", overlap, "\n")
  overlaps <- c(overlaps, overlap)
}
ggplot(data = tibble(iteration = c(1:length(overlaps)), overlap =  overlaps), aes(x = iteration, y = overlap)) + ylab("% Module Overlap") + xlab("Iteration") + ggtitle("Effect of Gene Swap Iterations on Module Differences") + geom_line() + theme_classic()


# module overlap check, highest degree swap
gcma <- counts_cer
gcmd <- counts_par
adja <- adjFromCounts(gcma)
adjd <- adjFromCounts(gcmd)
rankings <- rankHighestDegree(adjd)
clrsa <- recalculateModules(gcma)
clrsd <- recalculateModules(gcmd)
overlaps <- vector(mode = "double", length = 0)
OverlapCutoff <- 0.9
swapped <- rep(FALSE, ncol(adjd))
overlap <- calculateModuleOverlap(clrsa, clrsd)
overlaps <- c(overlaps, overlap)
while(overlap < OverlapCutoff) {
  swap_results <- swapHighestDegree(rankings, swapped, gcma, gcmd)
  swapped <- swap_results$updated_swap
  gcma <- swap_results$updated_gcm
  clrsa <- recalculateModules(gcma)
  overlap <- calculateModuleOverlap(clrsa, clrsd)
  cat("overlap = ", overlap, "\n")
  overlaps <- c(overlaps, overlap)
}
ggplot(data = tibble(iteration = c(1:length(overlaps)), overlap =  overlaps), aes(x = iteration, y = overlap)) + ylab("% Module Overlap") + xlab("Iteration") + ggtitle("Effect of Gene Swap Iterations on Module Differences") + geom_line()

# swap highest degree gene, adjacency network check
# (I'm not loving the adjacency check cause that'll obviously decrease linearlly as genes swap)
# gcma <- counts_cer
# gcmd <- counts_par
# adja <- adjFromCounts(counts_cer)
# adjd <- adjFromCounts(counts_par)
# diffs <- vector(mode = "double", length = 0)
# AdjDiffCutoff <- 0.001 # replicates have 0.006 - 0.009 adjacency diff
# rankings <- rankHighestDegree(adjd)
# swapped <- rep(FALSE, ncol(adjd))
# diff <- calculateAvgAdjacencyDiff(adja, adjd)
# diffs <- c(diffs, diff)
# while(diff > AdjDiffCutoff) {
#   swap_results <- swapHighestDegree(rankings, swapped, gcma, gcmd)
#   swapped <- swap_results$updated_swap
#   gcma <- swap_results$updated_gcm
#   adja <- adjFromCounts(gcma)
#   diff <- calculateAvgAdjacencyDiff(adja, adjd)
#   cat("diff = ", diff, "\n")
#   diffs <- c(diffs, diff)
# }
# plot(c(1:length(diffs)), diffs, ylab = "% Adjacency difference", xlab = "Iteration", main = "Effect of Gene Swap Iterations on Adjacency Matrix Differences")
# ggplot(data = tibble(iteration = c(1:length(diffs)), adj_diffs =  diffs), aes(x = iteration, y = adj_diffs)) + ylab("Average Difference in Edge Weight") + xlab("Iteration") + ggtitle("Effect of Gene Swap Iterations on Network Edge Differences") + geom_line() + theme_classic()
