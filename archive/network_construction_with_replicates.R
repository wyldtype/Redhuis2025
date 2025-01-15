#### Network construction from Barkai Data Using only the 143 conditions that had 2 replicates per species ####
# Following WGCNA tutorials: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Following Jeremy Miller's tutorial on Meta-analysis with WGCNA: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/
sapply(c("tidyverse", "edgeR", "WGCNA", "matrixStats"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/paralogs/DDivergence/Barkai_data_analysis")
options(stringsAsFactors = FALSE)

#### Reading & Cleaning Data ####
# read in count data (un-normalized) and sample info
load("data_files/Cleaned_Barkai_Data_RepSplit.RData")

# Note that hyc1 and hyp1 came from the same cell (and so did hyc2 and hyp2)
counts <- list(cer1 = cer_counts_rep1, cer2 = cer_counts_rep2, par1 = par_counts_rep1, par2 = par_counts_rep2,
               hyc1 = hyc_counts_rep1, hyc2 = hyc_counts_rep2, hyp1 = hyp_counts_rep1, hyp2 = hyp_counts_rep2)

nNetworks <- length(counts)

# normalize counts to cpm (don't want rpkm b/c it's tag-seq data---1 read = 1 3' UTR)
counts <- lapply(counts, cpm)

# for some unhinged reason, WGCNA has rows be samples and columns be genes from here out
counts <- lapply(counts, function(x) {return(t(x))})

# Verifying that the conditions are in the same order in all samples
getConditionOrder <- function(gcm) {
  conditions <- sample_info[match(rownames(gcm), sample_info$sample_name),]$condition
  return(conditions)
}
consensusOrder <- lapply(counts, getConditionOrder) %>% reduce(intersect)
cat("Percent samples in same order: ", length(consensusOrder)/nrow(counts[[1]]))

# Oh great, for downstream applications we need to make sure that none of the genes have zero variance
# but because we want the same set of genes in the end, we need to eliminate samples that have zero variance in ANY of the 12 networks
geneVar <- sapply(counts, colVars) # creates a nGene x length(counts) dimension matrix of variances
geneVarIsZero <- geneVar == 0
dropGenes <- apply(geneVarIsZero, 1, any)
sum(dropGenes)
counts <- lapply(counts, function(x) {
  x <- x[,!dropGenes]
  return(x)
})

# checking that there aren't any genes with lots of missing values or zero variance
gsg <- lapply(counts, function(x) {
  good <- goodSamplesGenes(x)
  cat("All OK: ", good$allOK, "\n")
  return(good)
})

# removing outlier samples by Euclidean distance clustering
par(cex = 0.4);
par(mar = c(0,4,2,0))
topOutliers <- vector(mode = "character", length = length(counts))
for (i in 1:nNetworks) {
  sampleTree <- hclust(dist(counts[[i]]), method = "average")
  topOutliers[i] <- sampleTree$labels[which(sampleTree$height == max(sampleTree$height))]
  plot(sampleTree, main = names(counts)[i])
  abline(h = 50000, col = "red")
} # we're looking for samples that are all off on their own---good clustering has two initial branches both leading to a decent total number of samples

topOutliers # TODO: WT_32 isn't the tallest sample shown in any of the hclust trees. Why is that? Is it cut off for some reason?
sample_info[match(topOutliers, sample_info$sample_name),] # all same condition, let's remove it
for (i in 1:nNetworks) {
  counts[[i]] <- counts[[i]][-which(rownames(counts[[i]]) == topOutliers[i]),]
}

# again check that sample conditions are still in same order
consensusOrder <- lapply(counts, getConditionOrder) %>% reduce(intersect)
cat("Percent samples in same order: ", length(consensusOrder)/nrow(counts[[1]]))

# filtering on top 3000 most connected (highest degree) genes based on average connectivity.
# (we decide on softPower in the parameter fitting section, but it's needed for 
# soft connectivity, so we'll just do the default value of 6 for now)
softPower <- 6
avgConn <- map_dfr(counts, softConnectivity, type = "unsigned", power = softPower) %>% 
  rowMeans() # I have no idea why it's rowMeans not colMeans... By my calculation map_dfr should get a 16 x nGenes matrix, so colMeans should return the mean connectivity of each gene but it gives a vector of length 16...
hist(avgConn)
keepGenes <- rank(-avgConn) <= 3000
counts_top3000 <- lapply(counts, function(x) {
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

# 6 looks good as it turns out 
# TODO: paradoxus is volatile as expected. But par2 much more so than par1. Repeat with a new replicate permutation (re-run clean_barkai_data.R and create a new Cleaned_Barkai_Data_RepSplit.RData file) to see if par2 is still the weird one. If it is, then there may be something not-random about the way I randomly partition replicates
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
      rankConn1 <- rank(softConnectivity(cts_list[[i]], type="unsigned", power=softPower))
      rankConn2 <- rank(softConnectivity(cts_list[[j]], type="unsigned", power=softPower))
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

# not good for QC, as there is real evolution captured here when comparing between species,
# but still interesting
ExprConn <- correlateExpressionAndConnectivity(cts_list = counts_top3000,
                                               ntwks_titles = names(counts_top3000))
library(gridExtra)
grid.arrange(grobs = ExprConn[1:8], ncol = 2) # ExprConn has 56 elements, 7 sets of 8, where each of the first 7 species are compared to the other8 (hyp2 is excluded b/c it's already compared in the first 7 sets)

#### Building networks ####
# Module selection: Adjacency matrix, topical overlap (TOM), and heirachical clustering

# set soft power (decided in parameter fitting section)
softPower <- 6

# getting adjacency matrix from count matrix (this takes a bit)
adjFromCounts <- function(cts) {
  adj <- adjacency(cts, type = "unsigned", power = softPower)
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

dev.off()
par(mfrow = c(2,4))
for (i in 1:nNetworks) {
  plot(geneTrees[[i]], xlab="", sub="", main=names(counts_top3000)[i],
       labels = FALSE, hang=0.04)
} # honestly these look better than anything I was getting with my fancy meta analysis strategies

# module selection
minModuleSize <- 15

# trying out two different tree cutting algorithms, dynamic and hybrid, and 3 values of deepSplit
colorList <- vector(mode = "list", length = nNetworks)
for (i in 1:nNetworks) {
  for (ds in 0:2) {
    dynamic <- cutreeDynamic(dendro = geneTrees[[i]], 
                             distM = dissTOMs[[i]],
                             deepSplit = ds, cutHeight = 0.995,
                             minClusterSize = minModuleSize,
                             pamRespectsDendro = FALSE) %>% labels2colors()
    
    chopped_tree <- cutreeHybrid(dendro = geneTrees[[i]], 
                                 distM = dissTOMs[[i]],
                                 deepSplit = ds, cutHeight = 0.995,
                                 minClusterSize = minModuleSize,
                                 pamRespectsDendro = FALSE)
    hybrid <- chopped_tree$labels %>% labels2colors()
    colorList[[i]] <- cbind(colorList[[i]], dynamic, hybrid)
  }
}
# compare tree cutting (module selection) methods
for (i in 1:nNetworks) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(unlist(colorList[[i]])),
                      groupLabels = c("dynamic deepSplit = 0", "hybrid deepSplit = 0",
                                      "dynamic deepSplit = 1", "hybrid deepSplit = 1",
                                      "dynamic deepSplit = 2", "hybrid deepSplit = 2"), 
                      main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}

# let's also try some static cuts with different heights
colorList <- vector(mode = "list", length = nNetworks)
for (i in 1:nNetworks) {
  for (h in c(0.99, 0.95, 0.9)) {
    static <- cutreeStaticColor(dendro = geneTrees[[i]],
                                cutHeight = h,
                                minSize = minModuleSize)
    colorList[[i]] <- cbind(colorList[[i]], static)
  }
}

for (i in 1:nNetworks) {
  plotDendroAndColors(dendro = geneTrees[[i]], 
                      colors = cbind(unlist(colorList[[i]])),
                      groupLabels = c("0.99", "0.95", "0.9"), 
                      main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}

# dynamic deepSplit=1 looks best, let's try some different cut heights
colorList <- vector(mode = "list", length = nNetworks)
for (i in 1:nNetworks) {
  for (ch in c(0.8, 0.85, 0.9, 0.95, 0.995)) {
    dynamic <- cutreeDynamic(dendro = geneTrees[[i]], 
                             distM = dissTOMs[[i]],
                             deepSplit = 1, cutHeight = ch,
                             minClusterSize = minModuleSize,
                             pamRespectsDendro = FALSE) %>% labels2colors()
    colorList[[i]] <- cbind(colorList[[i]], dynamic)
  }
}

for (i in 1:nNetworks) {
  plotDendroAndColors(dendro = geneTrees[[i]], 
                      colors = cbind(unlist(colorList[[i]])),
                      groupLabels = c("0.8", "0.85", "0.9", "0.95", "0.995"), 
                      main = names(geneTrees)[i],
                      dendroLabels = FALSE)
} # dynamic, deepsplit = 1, cutheight = 0.9 looks best

# final color setting
colors <- map2(geneTrees, dissTOMs, function(x, y) {
  dynamic <- cutreeDynamic(dendro = x,
                           distM = y,
                           deepSplit = 1, cutHeight = 0.9,
                           minClusterSize = minModuleSize,
                           pamRespectsDendro = FALSE) %>% labels2colors()
})

# merging close modules
# calculate eigengenes (eigenvector with largest eigenvalue aka 1st PC of expression matrix (woah way back to that thing) for each module)
# Note: colors is not an aesthetic argument, it's literally the only information on module membership
# moduleEigengenes for a p x n module X is equivalent to svd(t(scale(x)))$v[,1], the first principle component (length p) of X. If align = "along average", the sign is also potentially flipped so that the correlation between the average gene expression of each module across environments and its eigengene is non-negative (cor(A.E., M.E.) >= 0)
MEs <- map2(counts_top3000, colors, function(x, y) {
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
par(mfcol = c(2,4))
for (i in 1:nNetworks) {
  plot(METree[[i]], main = names(counts_top3000)[i], xlab = "", sub = "")
  abline(h=TestThresh, col = "red")
} # grey won't be merged by mergeCloseModules

# trying another threshold, this one looks better (less aggressive):
TestThresh <- 0.1

MEDissThresh <- 0.1
dev.off()
par(mfrow = c(2,4))
for (i in 1:nNetworks) {
  plot(METree[[i]], main = names(counts_top3000)[i], xlab = "", sub = "")
  abline(h=MEDissThresh, col = "red")
} 

# actually merging
unmergedColors <- colors
unmergedMEs <- MEs
merge <- map2(counts_top3000, colors, function(x,y) {
  return(mergeCloseModules(x, y, cutHeight = MEDissThresh, verbose = 3))
})
colors <- lapply(merge, function(x) {
  return(x$colors)
})

# visualizing final modules
for (i in 1:nNetworks) {
  plotDendroAndColors(geneTrees[[i]], colors[[i]],
                      names(geneTrees)[i], main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}

# matching module color labels
# before matching
for (i in 1:nNetworks) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors[[1]], colors[[2]], colors[[3]], colors[[4]], colors[[5]], colors[[6]], colors[[7]], colors[[8]]),
                      groupLabels = names(colors), main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}

# matching module colors across networks
# matchLabels looks for modules with significantly overlapping gene content and makes them the same color (not as sophisticated as the permutation tests in modulePreservation, it's just a fisher exact test)
colorsUnmatched <- colors

# if you need to re-set and try color matching again:
colors <- colorsUnmatched

# TODO: I don't think matchLabels takes into account majority gene membership enough when assigning labels (ex) I think paradoxus' largest module should be turquoise, based on how the genes line up but matchLables makes it green)
# it shouldn't be too hard to write my own function that just looks at whether each module in each species shares the majority of genes with a module of another color and then change colors to match... hmmm maybe it's a little tricky haha
colors <- lapply(colors, function(x) {
  return(matchLabels(source = x, reference = colors$cer1)) # through nothing better than trial, visualization, and error, I determind that hyp was the best reference. I tried par and hyc first and it totally screwed things up (maybe because they had more modules than cer or hyp)
})

# after matching
for (i in 1:nNetworks) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors[[1]], colors[[2]], colors[[3]], colors[[4]], colors[[5]], colors[[6]], colors[[7]], colors[[8]]),
                      groupLabels = names(colors), main = names(geneTrees)[i],
                      dendroLabels = FALSE)
} # par1's pink module really ought to be blue but otherwise looks alright

# recalculate MEs post merge/color match
MEs <- map2(counts_top3000, colors, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y)
  MEs <- MEoutput$eigengenes
  return(MEs)
})

# final network visualizations
for (i in 1:nNetworks) {
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

# outputting TOMplots (Note this DOES NOT WORK with more than one species' colors! It misaligns the 4 color rows)
TOMplotFileNames <- c("TOMplot_cer.pdf", "TOMplot_par.pdf", "TOMplot_hyc.pdf", "TOMplot_hyp.pdf")
for (i in 1:4) {
  pdf(file = paste0("../../../posters/", TOMplotFileNames[i]), colormodel = "cmyk")
  TOMplotSubset(dissTOM = dissTOMs[[i]], 
                colorList = colors, plotTitle = TOMplotTitles[i], sampleSize = 1000)
  dev.off()
} 

# We're done building our networks! Let's save all the relevant files
save(file = "data_files/Networks_RepSplit.RData", counts, counts_top3000, adjs, TOMs, dissTOMs,
     geneTrees, MEs, TOMs, colors, dissTOMs, softPower)
