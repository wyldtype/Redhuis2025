#### Heirarchical clustering of expression data
# Observation from combined-experiment clustering: each experiment
# only has a small number of distinct "shapes" to the expression
# a happy medium between experiment-specific clustering and combined-experiment
# clustering is to describe each gene by its combination of 
# expression "shapes" in each experiment
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")
options(stringsAsFactors = FALSE)
source("functions_for_figure_scripts.R")
# load("data_files/Clustering_Counts.RData")
load("data_files/Clustering_Counts_Allele.RData")
info <- info_allele
rm(info_allele)
set.seed(23)

############### Defining Clustering Functions ############ 
# Method:
# 1) Filter out genes with low var/mean ratio (low dispersion)
# 2) Correlation cluster remaining genes

# initial values for low var filter applied in corCluster (decided in following section)
var_thresh <- 3

# helper function for corCluster
# runs nIter clusterings on random subsets of genes 
# returns the majority vote for each gene's cluster
bootstrapClustering <- function(.filtered_cts, .nClust,
                                .nIter = 100, .frac = 0.75) {
  clusterSubset <- function(.cts) {
    gene_idxs <- sample(c(1:nrow(.cts)), size = nrow(.cts)*.frac, replace = FALSE)
    frac_cts <- .cts[gene_idxs,]
    cor_mat <- frac_cts |> t() |> cor(use = "pairwise.complete.obs")
    tree <- hclust(as.dist(-cor_mat), method = "average")
    topBranchHeight <- sort(tree$height, decreasing = TRUE)[.nClust]
    tree_labels <- cutreeStatic(tree, cutHeight = topBranchHeight,
                                minSize = 1)
    output <- rep(NA, nrow(.cts))
    output[gene_idxs] <- tree_labels
    return(output)
  }
  label_list <- vector(mode = "list", length = .nIter)
  for (i in c(1:.nIter)) {
    label_list[[i]] <- clusterSubset(.filtered_cts)
  }
  labeldf <- purrr::reduce(label_list, .f = cbind)
  cat("finished iterations, nGenes x nIter:", dim(labeldf), "\n")
  rownames(labeldf) <- rownames(.filtered_cts)
  # aligning labels between iterations
  matched <- apply(labeldf, 2, matchLabels, reference = labeldf[,1],
                   ignoreLabels = NA)
  matched[matched > .nClust] <- NA # matchLabels adds new labels to any gene over pThreshold
  labels <- apply(matched, 1, \(x) {
    label <- names(sort(table(x, useNA = "no"), decreasing = TRUE))[1]
    if (length(label) == 0) {
      return(NA)
    }
    return(as.numeric(label))
  })
  return(labels)
}
# given a count matrix where rows are genes and columns are samples/conditions,
# returns heirarchical cluster assignments as vector of length nGenes
# labels follow convention:
#   1) NA = too lowly expressed (mean expression < .min_expr)
#   2) 0 = too low variance (mean expression >= .min_expr & var expression < .min_var)
#   3) 1 - .nClust = clusters of co-varying genes
corCluster <- function(.cts, .nClust, .min_var = var_thresh,
                       .tree_too = FALSE, .gap_stat = FALSE, 
                       .just_counts = FALSE, .bootstrapIter = 100,
                       .bootstrapFrac = 0.75) {
  gene_idxs <- rownames(.cts)
  gene_means <- apply(.cts, 1, mean, na.rm = TRUE)
  gene_disps <- apply(.cts, 1, \(x) {var(x, na.rm = TRUE)/abs(mean(x, na.rm = TRUE))})
  gene_disps[gene_means == 0] <- 0 # avoiding divide by 0 NaN if gene is fully not expressed
  # identifying low var and low expr genes
  low_var_genes <- gene_idxs[gene_disps <= .min_var]
  cat(length(low_var_genes), "low variance genes, assigned to 0 label\n")
  good_clustering_genes <- setdiff(x = gene_idxs, y = low_var_genes)
  filtered_cts <- .cts[good_clustering_genes,]
  if (.gap_stat) {
    cluster_fun <- function(x, k) {
      list(cluster = cutree(hclust(d = as.dist(-cor(t(x), use = "pairwise.complete.obs")), 
                                   method = "average"), 
                            k = k))}
    gap_stat <- clusGap(filtered_cts, FUNcluster = cluster_fun, K.max = 10, B = 5, method = "firstmax")
    return(gap_stat)
  }
  if (.just_counts) {
    return(filtered_cts)
  }
  # clustering
  if (sum(apply(filtered_cts, 1, \(x) {all(is.na(x))})) != 0) { # if any rows (genes) have all NA values, they will cause cor to fail below
    cat("NA genes in counts matrix, returning counts matrix only\n")
    return(cor_mat)
  }
  if (.bootstrapIter > 1) {
    tree_labels <- bootstrapClustering(filtered_cts, .nClust = .nClust,
                        .nIter = .bootstrapIter, .frac = .bootstrapFrac)
    labeldf <- tibble(gene_name = rownames(filtered_cts),
                      label = tree_labels) |> 
      bind_rows(tibble(gene_name = low_var_genes,
                       label = 0))
    # making single tree for .treeToo = TRUE
    cor_mat <- filtered_cts |> t() |> cor(use = "pairwise.complete.obs")
    tree <- hclust(as.dist(-cor_mat), method = "average")
  }
  if (.bootstrapIter <= 1) {
    cor_mat <- filtered_cts |> t() |> cor(use = "pairwise.complete.obs")
    tree <- hclust(as.dist(-cor_mat), method = "average") # negative because hclust expects distance mat --- the higher the ij value the LESS similar genes i and j are
    topBranchHeight <- sort(tree$height, decreasing = TRUE)[.nClust]
    tree_labels <- cutreeStatic(tree, cutHeight = topBranchHeight,
                                minSize = 1) # guaranteeing no 0 class
    cat("cutting tree at height", topBranchHeight, "into", length(unique(tree_labels)), "clusters\n")
    labeldf <- tibble(gene_name = tree$labels,
                      label = tree_labels) |> 
      bind_rows(tibble(gene_name = low_var_genes,
                       label = 0))
  }
  # matching each label to each gene
  # rearranging genes into their original order
  outdf <- left_join(tibble(gene_name = rownames(.cts)),
                     labeldf, by = "gene_name")
  if (!.tree_too) {
    return(outdf)
  }
  if (.tree_too) {
    return(list(tree = tree, df = outdf))
  }
}
### tests for corCluster/bootstrapClustering
# toy genes
# known co-expressed genes in LowN
# should split these into genes highest at TP2 (^ genes)
# and genes lowest at TP2 (v genes)
toy_idxs <- c("YBR083W", "YBR172C", "YML015C", # v genes
              "YBR162W-A", "YKL196C", "YBR171W") # ^ genes (the excess of YBRs are coincidental, as I was just scrolling through that part of the module --- although the 171 172 W/C gene pair is probably overlapping)
toy_mat <- counts_list$par_LowN[toy_idxs,]
toydf <- makeDf(toy_mat, info, .join_by = "condition")
toyoutput <- corCluster(toy_mat, .nClust = 2, .tree_too = TRUE, 
                        .bootstrapIter = 0)
toydf <- left_join(toydf, toyoutput$df, by = "gene_name")
ggplot(toydf, aes(x = time_point_str, y = log2(expr + 1))) + 
  geom_line(aes(group = gene_name,
                color = gene_name)) +
  facet_wrap(~label)

# full dataset
toy_mat <- counts_list$par_LowN
toydf <- makeDf(toy_mat, info, .join_by = "condition")
# no bootstrap
test_labels <- corCluster(toy_mat, .nClust = 4, .min_var = 2,
                          .bootstrapIter = 0)
# yes bootstrap
test_labels <- corCluster(toy_mat, .nClust = 4, .min_var = 2,
                          .bootstrapIter = 10)
plotdf <- toydf |> 
  group_by(gene_name, time_point_str) |> 
  summarise(expr = mean(expr)) |> 
  ungroup() |> 
  reframe(mean_expr = mean(expr), 
          sd_expr = sd(expr),
          expr = expr,
          time_point_str = time_point_str,
          .by = "gene_name") |> 
  mutate(scaled_expr = (expr - mean_expr)/sd_expr) |> 
  left_join(y = test_labels, by = "gene_name")
ggplot(plotdf, aes(x = time_point_str, y = scaled_expr)) + 
  geom_line(aes(group = gene_name)) +
  facet_wrap(~label)

# # additional bootstrap comparisons
# toy_mat <- counts_list$par_LowN
# toydf <- makeDf(toy_mat, info, .join_by = "condition")
# labels_noboot <- corCluster(toy_mat, .nClust = 4, .min_var = 0.25)
# filtered_mat <- corCluster(toy_mat, .nClust = 4, .min_var = 0.25,
#                            .just_counts = TRUE)
# labels_yesboot <- bootstrapClustering(filtered_mat, .nClust = 4)
# labels_100boot <- bootstrapClustering(filtered_mat, .nClust = 4, .nIter = 100)
# labels_1000boot <- bootstrapClustering(filtered_mat, .nClust = 4, .nIter = 1000)
# label_compare <- labels_noboot |> filter(label != 0) |> select(label) |> cbind(labels_yesboot)
# sum(label_compare[,1] == label_compare[,2])
# sum(label_compare[,1] != label_compare[,2])
# label_compare <- cbind(labels_1000boot, labels_100boot)
# sum(label_compare[,1] == label_compare[,2])
# sum(label_compare[,1] != label_compare[,2]) # nice

# Wrapper function that converts named list of counts (from 1 experiment) into
# tidy counts data frame with nGenes * nConditions * length(counts list) number of rows
clusterCountsList <- function(.cts_list, .nClust = 4,
                              .bootstrapIter = 100,
                              .bootstrapFrac = 0.75,
                              .min_var = var_thresh,
                              .tree_too = FALSE,
                              .gap_stat = FALSE,
                              .just_counts = FALSE) {
  # giving each ortholog a unique name: YGR192C becomes YGR192C_cer
  cts <- map2(.cts_list, names(.cts_list), .f = \(x, nm) {
    nm <- gsub("_.*", "", nm)
    rownames(x) <- paste(nm, rownames(x), sep = "_")
    return(x)
  })
  dfs <- map(cts, makeDf, .info = info, .join_by = "condition")
  outdf <- purrr::reduce(dfs, bind_rows)
  if (all(unlist(map(cts, \(x) {all(colnames(x) == colnames(cts[[1]]))})))) {
    cat("counts are in same order, joining counts into single matrix\n")
    cts <- purrr::reduce(cts, .f = rbind)
  }
  else {
    cat("counts colnames don't match, returning counts list \n")
    return(cts)
  }
  # after collapsing replicates, samples should be in the same condition order
  # and able to be rbound
  if (.just_counts) {
    filtered_cts <- corCluster(cts, .nClust = .nClust, .min_var = .min_var,
                               .just_counts = TRUE, .bootstrapIter = .bootstrapIter,
                               .bootstrapFrac = .bootstrapFrac)
    return(filtered_cts)
  }
  if (!.tree_too) {
    labeldf <- corCluster(cts, .nClust = .nClust, .min_var = .min_var,
                          .bootstrapIter = .bootstrapIter,
                          .bootstrapFrac = .bootstrapFrac)
    outdf <- left_join(outdf, labeldf, by = "gene_name")
    return(outdf)
  }
  if (.tree_too) {
    output <- corCluster(cts, .nClust = .nClust, .min_var = .min_var,
                         .bootstrapIter = .bootstrapIter,
                         .bootstrapFrac = .bootstrapFrac,
                         .tree_too = TRUE)
    output$df <- left_join(outdf, output$df, by = "gene_name")
    return(output)
  }
}
# tests for clusterCountsList
# par is missing 1 condition versus cer
testout <- clusterCountsList(list("cer_LowN" = counts_list$cer_LowN,
                                  "par_LowN" = counts_list$par_LowN),
                             .bootstrapIter = 10, .min_var = 3)
testout$gene_name[1:10]

# given a dataframe with labels for each gene, returns a facet plot
# showing expression of each gene in each cluster (usually randomly 
# downsample option to save plotting computation)
plotClusters <- function(.df, .nDownsample = 0, .normalization = "scale",
                         .showProblem = FALSE) {
  if (.nDownsample != 0) {
    gene_idxs <- .df$gene_name |> sample(size = .nDownsample, replace = FALSE)
    .df <- .df |> filter(gene_name %in% gene_idxs)
  }
  if (.normalization == "scale") {
    plotdf <- .df |>
      reframe(mean_expr = mean(expr, na.rm = TRUE), 
              sd_expr = sd(expr, na.rm = TRUE),
              expr = expr,
              time_point_num = time_point_num,
              label = label,
              .by = "gene_name") |> 
      mutate(plot_expr = (expr - mean_expr)/sd_expr)
  }
  if (.normalization == "log2") {
    plotdf <- .df |> 
      mutate(plot_expr = log2(expr + 1))
  }
  p <- ggplot(plotdf, aes(x = time_point_num, y = plot_expr)) + 
    geom_line(aes(group = gene_name)) +
    # adding mean expr line for each cluster:
    geom_line(data = summarise(group_by(plotdf, time_point_num, label),
                                mean_expr = mean(plot_expr, na.rm = TRUE)),
              aes(x = time_point_num, y = mean_expr),
              color = "gold") +
    # adding mean expr for all low expressed genes (group_by won't create an NA group):
    geom_line(data = summarise(group_by(filter(plotdf, is.na(label)),
                                        time_point_num),
                               mean_expr = mean(plot_expr, na.rm = TRUE),
                               label = NA),
              aes(x = time_point_num, y = mean_expr),
              color = "gold") +
    facet_wrap(~ label)
  return(p)
}

#### Picking low expression and low variance thresholds ####
var_thresh <- 3

getCollapsedCountsByExperiment <- function(.experiment) {
  name_cer <- paste("cer", .experiment, sep = "_")
  name_par <- paste("par", .experiment, sep = "_")
  .cts_list <- list(name_cer = counts_list[[name_cer]],
                    name_par = counts_list[[name_par]])
  cts <- map2(.cts_list, names(.cts_list), .f = \(x, nm) {
    nm <- gsub("_.*", "", nm)
    rownames(x) <- paste(nm, rownames(x), sep = "_")
    return(x)
  }) 
  dfs <- map(cts, makeDf, .info = info, .join_by = "condition")
  outdf <- purrr::reduce(dfs, bind_rows)
  # after collapsing replicates, samples should be in the same condition order
  # and able to be rbound
  cts <- cts |> 
    purrr::reduce(rbind)
  return(cts)
}
# HAP4
toy_mat_collapsed <- getCollapsedCountsByExperiment("HAP4")
# note that it is important to specify whether we're thresholding based on log2(mean(expr))
# versus mean(log2(expr)):
plot(rowMeans(log2(toy_mat_collapsed)), 
     log2(rowMeans(toy_mat_collapsed)))
abline(a = 0, b = 1, col = "gold") 
# they're mostly a similar measure, 
# but for genes with high variance,
# log2(mean(expr)) tends to be higher 
# (hence those spikes of genes well above the y=x)
# we are thresholding by log2(mean(expr)) because it makes more intuitive sense
plot(log2(rowMeans(toy_mat_collapsed)), log2(rowVars(toy_mat_collapsed)))
mod <- lm(log2(rowVars(toy_mat_collapsed) + 1) ~ log2(rowMeans(toy_mat_collapsed) + 1))
abline(a = mod$coefficients[1], 
       b = mod$coefficients[2], col = "gold")
# the point where the low expr genes "lift off" this mean-var relationship line is where 
# lowly expressed genes have an abnormally high variance, variance that is more likely
# due to noise than response to the environment, and should therefore be in the low expr group
# low var filtering:
# var thresh is actually a dispersion thresh --- var/mean
# to filter out lowly varying genes, we apply disp threshold after the expr threshold
# the higher you're expressed, the more you need to vary in order to not be put in low var category
abline(a = 0, b = log2(var_thresh), col = "blue") 
# LowN
toy_mat_collapsed <- getCollapsedCountsByExperiment("LowN")
plot(log2(rowMeans(toy_mat_collapsed)), log2(rowVars(toy_mat_collapsed)))
mod <- lm(log2(rowVars(toy_mat_collapsed) + 1) ~ log2(rowMeans(toy_mat_collapsed) + 1))
abline(a = mod$coefficients[1], 
       b = mod$coefficients[2], col = "gold") 
abline(a = 0, b = log2(var_thresh), col = "blue")
# LowPi
toy_mat_collapsed <- getCollapsedCountsByExperiment("LowPi")
plot(log2(rowMeans(toy_mat_collapsed)), log2(rowVars(toy_mat_collapsed)))
mod <- lm(log2(rowVars(toy_mat_collapsed) + 1) ~ log2(rowMeans(toy_mat_collapsed) + 1))
abline(a = mod$coefficients[1], 
       b = mod$coefficients[2], col = "gold") 
abline(a = 0, b = log2(var_thresh), col = "blue")
# Heat
toy_mat_collapsed <- getCollapsedCountsByExperiment("Heat")
plot(log2(rowMeans(toy_mat_collapsed)), log2(rowVars(toy_mat_collapsed)))
mod <- lm(log2(rowVars(toy_mat_collapsed) + 1) ~ log2(rowMeans(toy_mat_collapsed) + 1))
abline(a = mod$coefficients[1], 
       b = mod$coefficients[2], col = "gold") 
abline(a = 0, b = log2(var_thresh), col = "blue")
# Cold
toy_mat_collapsed <- getCollapsedCountsByExperiment("Cold")
plot(log2(rowMeans(toy_mat_collapsed)), log2(rowVars(toy_mat_collapsed)))
mod <- lm(log2(rowVars(toy_mat_collapsed) + 1) ~ log2(rowMeans(toy_mat_collapsed) + 1))
abline(a = mod$coefficients[1], 
       b = mod$coefficients[2], col = "gold") 
abline(a = 0, b = log2(var_thresh), col = "blue")

#### Clustering Tests ####
# 1) Do the same genes cluster when using scaled counts or unscaled?
# unscaled, reference labels
toy_idxs <- sample(rownames(counts_list$cer_CC), 1000, replace = FALSE)
toy_mat <- counts_list$cer_CC[toy_idxs,]
toydf <- makeDf(toy_mat, .info = info, .join_by = "condition")
labelsdf <- corCluster(toy_mat, .nClust = 4, 
                       .min_var = 0, .bootstrapIter = 0) # no 0 or NA genes, for comparability, as those thresholds will be different
table(labelsdf$label, useNA = "always")
labelsdf <- rename(labelsdf, "unscaled"="label")

# scaled just prior to clustering
scaleGene <- function(.cts) {
    out_cts <- apply(.cts, 1, \(x) {
      output <- (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
      if (sd(x, na.rm = TRUE) == 0) {
        output <- x - mean(x, na.rm = TRUE)
      }
      return(output)
    }) |> t()
    return(out_cts)
  }
toy_mat <- counts_list$cer_CC[toy_idxs,]
toydf <- makeDf(toy_mat, .info = info, .join_by = "condition")
toy_mat_scaled <- scaleGene(toy_mat)
labelsdf_after <- corCluster(toy_mat_scaled, .nClust = 4, 
                             .min_var = 0, .bootstrapIter = 0)
table(labelsdf_after$label, useNA = "always")
labelsdf <- left_join(labelsdf, rename(labelsdf_after, "after"="label"),
                      by = "gene_name")

# scaled prior to mean across replicates
scaled_counts <- map(counts_list, scaleGene)
toy_mat <- scaled_counts$cer_CC[toy_idxs,]
toydf <- makeDf(toy_mat, .info = info, .join_by = "condition")
labelsdf_before <- corCluster(toy_mat, .nClust = 4,
                              .min_var = 0, .bootstrapIter = 0)
table(labelsdf_before$label, useNA = "always")
labelsdf <- left_join(labelsdf, rename(labelsdf_before, "before"="label"),
                      by = "gene_name")
# View(labelsdf)
apply(select(labelsdf, -one_of("gene_name")), 1, n_distinct) |> 
  table() 
# conclusion: scaling doesn't affect clustering 
# (but it will affect which genes are flagged as low var or 
# low expression if you don't change thresholds)

#### Choosing number of clusters per experiment ####
findTreeDropoff <- function(.tree_heights, .kmax) {
  heightdiffs <- sort(.tree_heights, decreasing = TRUE)[1:(.kmax - 1)] - 
    sort(.tree_heights, decreasing = TRUE)[2:.kmax]
  branchdiffs <- heightdiffs[1:(.kmax - 1)] - heightdiffs[2:.kmax]
  df <- tibble(branchdiffs = branchdiffs, nclust = c(2:.kmax))
  return(df)
}
# Saturated Growth, 2 clusters
tree_test <- clusterCountsList(counts_list[grepl("HAP4", names(counts_list))],
                               .tree_too = TRUE, .bootstrapIter = 0)
plot(tree_test$tree, labels = FALSE)
treedf <- findTreeDropoff(tree_test$tree$height, .kmax = 20)
plot(treedf$nclust, treedf$branchdiffs, pch = 20)
text(treedf$nclust, treedf$branchdiffs, labels = treedf$nclust, pos = 4)
treedf$nclust[which.max(treedf$branchdiffs)] # 2 clusters

# Low Nitrogen, 2 clusters
tree_test <- clusterCountsList(counts_list[grepl("LowN", names(counts_list))],
                               .tree_too = TRUE, .bootstrapIter = 0)
plot(tree_test$tree, labels = FALSE)
treedf <- findTreeDropoff(tree_test$tree$height, .kmax = 20)
plot(treedf$nclust, treedf$branchdiffs, pch = 20)
text(treedf$nclust, treedf$branchdiffs, labels = treedf$nclust, pos = 4) # 2 clusters
treedf$nclust[which.max(treedf$branchdiffs)]
# note if you think that tree looks concerningly flat:
tree_test$tree$height |> round(digits = 3) |> table() # A ton of correlation captured in the first few branches, most genes have a height close to -1, the limit
sum(tree_test$tree$height == -1)
# LowN_output <- clusterCountsList(counts_list[grepl("LowN", names(counts_list))],
#                                  .nClust = 4, .bootstrapIter = 100, .bootstrapFrac = 0.90)
# plotClusters(LowN_output, .nDownsample = 0)

# Heat Stress, 4 clusters
tree_test <- clusterCountsList(counts_list[grepl("Heat", names(counts_list))],
                               .tree_too = TRUE, .bootstrapIter = 0)
plot(tree_test$tree, labels = FALSE)
treedf <- findTreeDropoff(tree_test$tree$height, .kmax = 20)
plot(treedf$nclust, treedf$branchdiffs, pch = 20)
text(treedf$nclust, treedf$branchdiffs, labels = treedf$nclust, pos = 4)
treedf$nclust[which.max(treedf$branchdiffs)] # 4 clusters
# Heat_output <- clusterCountsList(counts_list[grepl("Heat", names(counts_list))],
#                                  .nClust = 5, .bootstrapIter = 100)
# plotClusters(Heat_output, .nDownsample = 0)

# Cold Stress, 3 clusters
tree_test <- clusterCountsList(counts_list[grepl("Cold", names(counts_list))],
                               .tree_too = TRUE, .bootstrapIter = 0)
plot(tree_test$tree, labels = FALSE)
treedf <- findTreeDropoff(tree_test$tree$height, .kmax = 20)
plot(treedf$nclust, treedf$branchdiffs, pch = 20)
text(treedf$nclust, treedf$branchdiffs, labels = treedf$nclust, pos = 4)
treedf$nclust[which.max(treedf$branchdiffs)] # 3 clusters
# Cold_output <- clusterCountsList(counts_list[grepl("Cold", names(counts_list))],
#                                  .nClust = 3, .bootstrapIter = 100)
# plotClusters(Cold_output, .nDownsample = 0)

# Cell Cycle, 2 clusters
tree_test <- clusterCountsList(counts_list[grepl("CC", names(counts_list))],
                               .tree_too = TRUE, .bootstrapIter = 0)
plot(tree_test$tree, labels = FALSE)
treedf <- findTreeDropoff(tree_test$tree$height, .kmax = 20)
plot(treedf$nclust, treedf$branchdiffs, pch = 20)
text(treedf$nclust, treedf$branchdiffs, labels = treedf$nclust, pos = 4)
treedf$nclust[which.max(treedf$branchdiffs)] # 2 clusters
# CC_output <- clusterCountsList(counts_list[grepl("CC", names(counts_list))],
#                                .nClust = 2, .bootstrapIter = 100)
# plotClusters(CC_output, .nDownsample = 0)

# Low Phosphorus, 2 clusters 
tree_test <- clusterCountsList(counts_list[grepl("LowPi", names(counts_list))],
                               .tree_too = TRUE, .bootstrapIter = 0)
plot(tree_test$tree, labels = FALSE)
treedf <- findTreeDropoff(tree_test$tree$height, .kmax = 20)
plot(treedf$nclust, treedf$branchdiffs, pch = 20)
text(treedf$nclust, treedf$branchdiffs, labels = treedf$nclust, pos = 4)
treedf$nclust[which.max(treedf$branchdiffs)] # 3 clusters
# LowPi_output <- clusterCountsList(counts_list[grepl("LowPi", names(counts_list))],
#                                .nClust = 3, .bootstrapIter = 100)
# plotClusters(LowPi_output, .nDownsample = 0)

#### Clustering ####
# clustering all-4 species per experiment, with number of clusters decided in previous section

# change for different parameter values we'll use
# var_thresh <- 1
var_thresh <- 3
# var_thresh <- 5

clusterdf_list <- vector(mode = "list", length = 0)
nclust_lookup <- tibble(experiment = c("HAP4", "LowPi", "CC", "LowN", "Cold", "Heat"),
                        nclust = c(2, 2, 2, 2, 2, 2))
nclust_lookup
for (e in nclust_lookup$experiment) {
  nclust <- nclust_lookup |> filter(experiment == e) |> 
    select(nclust) |> as.numeric()
  cat("*********** working on", nclust, "clusters in", e, "*********** \n")
  output <- clusterCountsList(counts_list[grepl(e, names(counts_list))], 
                          .nClust = nclust, .tree_too = TRUE,
                          .min_var = var_thresh)
  clusterdf_list[[paste(e, nclust, sep = "_")]] <- output
}

# getting gene clusters
getClusterCombination <- function(.clust_list) {
  cluster_comb <- tibble()
  for (nm in names(.clust_list)) {
    e <- gsub("_.*", "", nm)
    nclust <- gsub(".*_", "", nm)
    cat("working on", paste(e, nclust, sep = "_"), "\n")
    e_clust <- clusterdf_list[[paste(e, nclust, sep = "_")]]$df |> 
      select(label, gene_name) |> 
      unique()
    e_clust$gene_ID <- map(e_clust$gene_name, \(.g) {
      return(gsub(".*_", "", .g))
    }) |> unlist()
    e_clust$species <- map(e_clust$gene_name, \(.g) {
      return(gsub("_.*", "", .g))
    }) |> unlist()
    e_clust <- e_clust |> 
      select(gene_ID, species, label) |> 
      pivot_wider(id_cols = gene_ID, names_from = species, values_from = label) |> 
      mutate(experiment = e)
    cluster_comb <- bind_rows(cluster_comb, e_clust)
  }
  return(cluster_comb)
}
clusterdf <- getClusterCombination(clusterdf_list)

# saving
# save(clusterdf, clusterdf_list, file = "data_files/CorrelationClustering1Disp.RData")
# save(clusterdf, clusterdf_list, file = "data_files/CorrelationClustering3Disp.RData")
# save(clusterdf, clusterdf_list, file = "data_files/CorrelationClustering5Disp.RData")
save(clusterdf, clusterdf_list, file = "data_files/CorrelationClustering3Disp_Allele.RData")

### Visualizing cluster expression patterns per experiment
load("data_files/CorrelationClustering3Disp.RData")
plotlist <- vector(mode = "list", length = 0)
for (nm in names(clusterdf_list)) {
  plotlist[[nm]] <- plotClusters(clusterdf_list[[nm]]$df, .nDownsample = 4000) + 
    ggtitle(nm)
}
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_clusters.pdf",
    width = 7, height = 16)
ggarrange(plotlist = plotlist, nrow = 4, ncol = 2)
dev.off()

### Supplementary figure: heiarchical clustering gene trees
for (nm in names(clusterdf_list)) {
  pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/tree_",
                    nm, ".pdf"), width = 9, height = 3)
  plot(clusterdf_list[[nm]]$tree, labels = FALSE, xlab = "", 
       ylab = "distance at merge",
       main = nm)
  dev.off()
}

########################## Archive ########################## 
#### centering & scaling ####
# # test count
# test_cts <- counts_list[[sample(names(counts_list), 1)]]
# test_row_idx <- sample(c(1:nrow(test_cts)), 1)
# test_col_idx <- sample(c(1:ncol(test_cts)), 1)
# rownames(test_cts)[test_row_idx]
# # what it should be centered and scaled:
# (test_cts[test_row_idx, test_col_idx] - mean(test_cts[test_row_idx,], na.rm = TRUE))/sd(test_cts[test_row_idx,], na.rm = TRUE)
# # what is is using our function (that converts sd=0 to vecs of 0 instead of NaN)
# test_scaled <- scaleGene(test_cts)
# test_scaled[test_row_idx, test_col_idx]
# # what it is centered and scaled (gold standard function):
# test_scaled <- scale(t(test_cts)) |> t()
# test_scaled[test_row_idx, test_col_idx]
# attr(test_scaled,"scaled:scale")[test_row_idx]
# attr(test_scaled,"scaled:center")[test_row_idx]
# # applying to all counts
# scaled_counts <- map(counts_list, .f = scaleGene)

#### single species corCluster ####
# archived b/c the all 4 clustering works so well, and has a built in
# decision criteria for conservation of plasticity (whether orthologs are clustered together)
# corClusterBySpecies <- function(.cts, .species = "all", .experiment, .plotToo = FALSE) {
#   if (.species %in% c("cer", "par")) {
#     info <- sample_info
#   }
#   if (.species %in% c("hyc", "hyp")) {
#     info <- sample_info_allele
#   }
#   outputdf <- as_tibble(t(.cts)) |> 
#     bind_cols(tibble(gene_name = colnames(.cts))) |> 
#     pivot_longer(cols = rownames(.cts),
#                  names_to = "sample_name", values_to = "expr") |> 
#     left_join(info, by = "sample_name") |> 
#     filter(experiment %in% .experiment) 
#   outputdf <- left_join(x = outputdf, 
#                         y = group_by(outputdf, gene_name) |> 
#                           summarise(mean_expr = mean(expr),
#                                     sd_expr = sd(expr)),
#                         by = "gene_name", relationship = "many-to-one") |> 
#     mutate(scaled_expr = (expr - mean_expr)/sd_expr)
#   filtered_cts <- .cts[outputdf$sample_name, unique(outputdf$gene_name)]
#   # setting genes with no var to 0 cor b/c hclust can't handle NA
#   # and removing genes not present in Hot/cold
#   if (.experiment %in% c("Heat", "Cold")) {
#     all_na_gene_idxs <- apply(filtered_cts, 2, \(x) {all(is.na(x))})
#     filtered_cts <- filtered_cts[,!all_na_gene_idxs]
#   }
#   no_var_genes <- colnames(filtered_cts)[which(colVars(filtered_cts) == 0)]
#   cor_mat <- filtered_cts |> cor(use = "pairwise.complete.obs")
#   cor_mat[rownames(cor_mat) %in% no_var_genes, ] <- 0
#   cor_mat[, colnames(cor_mat) %in% no_var_genes] <- 0
#   if (sum(is.na(cor_mat)) != 0) {
#     cat("NAs in cor_mat, returning cor_mat only\n")
#     return(cor_mat)
#   }
#   tree <- hclust(as.dist(-cor_mat), method = "average") # negative because hclust expects distance mat --- the higher the ij value the LESS similar genes i and j are
#   topBranchHeight <- sort(tree$height, decreasing = TRUE)[2]
#   tree_labels <- cutreeDynamicTree(tree, maxTreeHeight = 0.999,
#                                    minModuleSize = max(c(0.1*nrow(cor_mat))), 10)
#   labeldf <- tibble(gene_name = tree$labels, label = tree_labels)
#   outputdf <- left_join(outputdf, labeldf, by = "gene_name",
#                         relationship = "many-to-one")
#   if (.plotToo) {
#     return(list(df = outputdf, tree = tree))
#   }
#   return(tree_labels)
# }
#### Saturated Growth vertical slice ####
# 
# # Answering two questions QC questions:
# # First, are single species clusters capturing more real variance
# # than clustering all 4 species? 
# # Second, is using 2 clusters best?
# all4_output <- clusterCountsList(counts_list[grepl("HAP4", names(counts_list))],
#                                  .nClust = 2)
# # number of genes in each cluster:
# all4_output |> select(gene_name, label) |> unique() |> select(label) |> table()
# plotClusters(all4_output, .nDownsample = 100) # messy, but you can see 1 increases and 2 decreases
# # individual species clusters
# cer_output <- clusterCountsList(counts_list["cer_HAP4"], .nClust = 2)
# plotClusters(cer_output, .nDownsample = 100)
# par_output <- clusterCountsList(counts_list["par_HAP4"], .nClust = 2)
# plotClusters(par_output, .nDownsample = 100)
# # hyc_output <- clusterCountsList(counts_list["hyc_HAP4"], .nClust = 2)
# # plotClusters(hyc_output, .nDownsample = 100)
# # hyp_output <- clusterCountsList(counts_list["hyp_HAP4"], .nClust = 2)
# # plotClusters(hyp_output, .nDownsample = 100)
# # conclusion: no difference between all 4 and single species
# # how many clusters? What does the tree look like?
# test <- clusterCountsList(counts_list[grepl("HAP4", names(counts_list))],
#                           .nClust = 2, .tree_too = TRUE)
# plot(test$tree, labels = FALSE) # 2 groups, only obvious height separation at 2
# # Now asking experimental questions:
# # Did any genes switch from cluster 1, 2, NA, or 0? in how many species?
# plotdf <- test$df |> select(gene_name, label) |> 
#   mutate(species = gsub("_.*", "", gene_name),
#          gene_name = gsub(".*_", "", gene_name)) |> 
#   unique() |> 
#   pivot_wider(id_cols = gene_name, names_from = species,
#               values_from = label)
# plotdf$code <- apply(plotdf, 1, \(x) {
#   paste(as.numeric(x["cer"]), 
#         as.numeric(x["par"]),
#         sep = " ")
# })
# length(unique(plotdf$code)) # way less than 256 possible codes
# sort(table(plotdf$code), decreasing = TRUE)[1:20] 
# code_order <- names(sort(table(plotdf$code), decreasing = TRUE))
# # upset plot with categories 1111, 2222, 1211, NANANANA, 0202, etc.
# #       probably with minimum count so that there are ~ 20 categories
# ggplot(filter(plotdf, code %in% code_order), aes(x = code)) + 
#   geom_bar() +
#   scale_x_discrete(breaks = code_order, limits = code_order) +
#   theme(axis.text.x = element_text(angle = 90))
#  
# # 3) Low N maybe did seem to have a difference between all 4 and 
# # single species (par looks especially nice)
# # Full species cluster output to group into - v ^ / and \ like the
# all4_output <- clusterCountsList(counts_list[grepl("LowN", names(counts_list))])
# plotClusters(all4_output, .nDownsample = 4000)
# 
# # par
# par_output <- clusterCountsList(counts_list["par_LowN"])
# plotClusters(par_output, .nDownsample = 4000)
# # in case you're curious, 
# # this is how we generated the par LowN plot in the corCluster test
# # (it looks nicer, but that's mostly the categorical timepoint)
# par_output <- clusterCountsList(counts_list["par_LowN"],
#                                 .min_var = 0.1)
# plotdf <- par_output |> 
#   group_by(gene_name, time_point_str, label) |> 
#   summarise(expr = mean(expr)) |> 
#   ungroup() |> 
#   reframe(mean_expr = mean(expr), 
#           sd_expr = sd(expr),
#           expr = expr,
#           time_point_str = time_point_str,
#           label = label,
#           .by = "gene_name") |> 
#   mutate(scaled_expr = (expr - mean_expr)/sd_expr)
# ggplot(plotdf, aes(x = time_point_str, y = scaled_expr)) + 
#   geom_line(aes(group = gene_name)) +
#   facet_wrap(~label)
# 
# # TODO: continue adapting the individual species to include - genes
# # whole experiment LowN in single species: should split into 5 categories: - v ^ \ and / genes
# # cer
# testoutput <- clusterCountsList(counts_list["cer_LowN"])
# plotClusters(testoutput, .nDownsample = 500)
# # hardcoding -, v, ^, /, and \ labels
# testoutput$label <- map(testoutput$label, \(x) {
#   if (is.na(x)) {
#     return("low")
#   }
#   if (x == 0) {
#     return("-")
#   }
#   if (x == 1) {
#     return("/")
#   }
#   if (x == 2) {
#     return("\\") # first backslash is to escape backslash
#   }
#   if (x == 3) {
#     return("^") 
#   }
#   if (x == 4) {
#     return("v")
#   }
# }) |> unlist()
# testdf_all4 <- rename(testoutput, "label_cer"="label")
# testdf_all4$gene_name <- gsub("cer_", "", testdf_all4$gene_name) 
# testdf_all4 <- testdf_all4 |> 
#   select(gene_name, label_cer) |> 
#   unique()
# 
# # par
# testoutput <- clusterCountsList(counts_list["par_LowN"])
# plotClusters(testoutput, .nDownsample = 500)
# testoutput$gene_name <- gsub("par_", "", testoutput$gene_name) 
# testoutput <- rename(testoutput, "label_par"="label") |> 
#   select(gene_name, label_par) |> 
#   unique()
# testdf_all4 <- left_join(testdf_all4, testoutput, by = "gene_name")
# table(paste(testdf_all4$label_cer, testdf_all4$label_par)) |> sort(decreasing = TRUE) # looks like the only real shapes in par are / and \
# 
# # Where I left off: After excluding low expr and low var genes, single-species
# # clustering isn't looking any better than all-4 clustering.
# # Since all-4 is easier and simpler, I'll continue forward with that but if
# # I want to continue this troubleshooting, the key is to label genes based on their expression shape
# # so that labels are comparable between species. If it isn't clear which shape is applicable,
# # the tie-breaker is which cluster shares more of those genes from other species
# 
# # how many genes are in the same shape group in both species?
# View(testdf_all4)
# sum(testdf_all4$label_cer == testdf_all4$label_par, na.rm = TRUE)
# 
# # Next find this set of conserved genes then run below all-4 clustering
# # to see if these conserved genes are still together
# # (spoiler: I don't think they will be because the clusters are not
# # visibly /, \, v, and ^)

#### Spectral clustering ####
# archived b/c heiarchical works so well for correlation clustering
# adjFromCounts <- function(.cts) {
#   # adj <- cor(.cts, use = "pairwise.complete.obs")
#   adj <- matrix(NA, nrow = ncol(.cts), ncol = ncol(.cts))
#   # TODO: make i > j, run in half the time
#   for (i in c(1:(ncol(adj) - 1))) {
#     cat(i, "/", ncol(adj), "\n")
#     for (j in c((i+1):ncol(adj))) {
#       non_na <- cbind(.cts[,i], .cts[,j])[!is.na(.cts[,i]) & !is.na(.cts[,j]),]
#       dcor_ij <- dcor(non_na[,1], non_na[,2])
#       adj[i, j] <- dcor_ij
#       adj[j, i] <- dcor_ij
#     }
#   }
#   diag(adj) <- 0
#   return(adj)
# }
# test <- adjFromCounts(counts_all4$cer[,1:100])
# max(test)
# min(test[test > 0])
# test[1:3,1:3]
# test <- adjFromCounts_WGCNA(counts_all4$cer[1:10,1:100])
# test[1:3,1:3]
# max(test)
# min(test)
# 
# # generating adj matrices
# adjs <- lapply(list(counts_all4$cer, counts_all4$par), adjFromCounts)
# save(adjs, file = "data_files/Adjacency_dcor.RData")
