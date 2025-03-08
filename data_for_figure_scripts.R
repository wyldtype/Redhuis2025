sapply(c("dplyr", "purrr", "tidyr", "ggpubr", "readr",
         "data.table", "ggplot2", "data.table", "msir", 
         "WGCNA", "energy", "matrixStats"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2025/")
source("functions_for_figure_scripts.R")

#### Reading in Count Data ####
load("data_files/Cleaned_Count_Data.RData")
load("data_files/Cleaned_Count_Data_AlleleSpecific.RData")

#### filtering low expr ####
# changing to non-log2 scale for actual filtering (log2 scale is for visualizing)
cutoffExpr <- 30

# Criteria: mean expr less than threshold (30 cpm, not log scale) in
# cer, par, hyc, and hyp in all experiments
getGoodExprGeneNames <- function(.organism, .allele, .experiment, .expr_thresh) {
  if (.organism != "hyb") {
    sample_names <- sample_info |> filter(allele == .allele &
                                            experiment == .experiment) |> 
      select(sample_name) |> pull()
    cts <- counts[, sample_names]
    means <- rowMeans(cts)
    return(rownames(cts)[means >= .expr_thresh])
  }
  if (.organism == "hyb") {
    sample_names <- sample_info_allele |> filter(allele == .allele &
                                                   experiment == .experiment) |> 
      select(sample_name) |> pull()
    cts <- counts_allele[, colnames(counts_allele) %in% sample_names]
    means <- rowMeans(cts)
    return(rownames(cts)[means >= .expr_thresh])
  }
}
# # tests for getGoodExprGeneNames
# # in Scer
# test <- getGoodExprGeneNames(.organism = "cer", .allele = "cer",
#                              .experiment = "HAP4", .expr_thresh = cutoffExpr)
# "YPR199C" %in% test # should be
# "YFL051C" %in% test # shouldn't be
# # in Spar
# test <- getGoodExprGeneNames(.organism = "par", .allele = "par",
#                              .experiment = "HAP4", .expr_thresh = cutoffExpr)
# "YPR199C" %in% test # shouldn't be
# "YFL051C" %in% test # should be

# applying to all genes/groups/experiments
griddf <- expand_grid(organism = c("cer", "par", "hyb"),
                      allele = c("cer", "par"),
                      experiment = unique(sample_info$experiment))
keep <- map(c(1:nrow(griddf)), \(i) getGoodExprGeneNames(griddf$organism[i],
                                                         griddf$allele[i],
                                                         griddf$experiment[i],
                                                         .expr_thresh = cutoffExpr)) |> 
  unlist() |> unique()
"YPR199C" %in% keep
"YFL051C" %in% keep

# filtering lowly expressed genes
# our 46 TF deletions
TFdel_lookup <- read_delim("data_files/downloaded_genomes_and_features/yeastract_46TFs.csv", col_names = FALSE, col_select = c(1,2), delim = ";") # gets some warnings, but so far has been fine
colnames(TFdel_lookup) <- c("common", "systematic")

length(keep) # number of genes we're keeping
TFdel_lookup$common[!(TFdel_lookup$systematic %in% keep)] # all the TFdel genes we'd be removing based on this expression criteria
TFdel_lookup$common[which(!(TFdel_lookup$systematic %in% rownames(counts)))] # but some aren't annotated in unfiltered dataset either
# preserving TFs
keep <- c(keep, TFdel_lookup$systematic) |> unique()
keep <- keep[keep %in% rownames(counts)]
length(keep)
# filtering
counts <- counts[keep,]
counts_allele <- counts_allele[keep,]
#### Filtering genotypes ####
full_counts <- counts
full_sample_info <- sample_info
full_counts_allele <- counts_allele
full_sample_info_allele <- sample_info_allele
sample_info <- sample_info |> filter(genotype == "WT")
sample_info_allele <- sample_info_allele |> filter(genotype == "WT")
counts <- counts[,sample_info$sample_name]
counts_allele <- counts_allele[,sample_info_allele$sample_name]
# TFdel counts
sample_info_tfdel <- full_sample_info |> filter(experiment == "LowN")
sample_info_tfdel_allele <- full_sample_info_allele |> filter(experiment == "LowN")
counts_tfdel <- full_counts[, sample_info_tfdel$sample_name]
counts_tfdel_allele <- full_counts_allele[, sample_info_tfdel_allele$sample_name]

rm(full_sample_info, full_sample_info_allele, 
   full_counts, full_counts_allele)

# Saving TFdel counts for analysis and Fig script
save(sample_info_tfdel, sample_info_tfdel_allele,
     counts_tfdel, counts_tfdel_allele, file = "data_files/TFdel.RData")

#### taking mean count across replicates (aka collapsing) ####

# collapsing
collapsed <- collapseReplicates(sample_info, counts)
collapsed_allele <- collapseReplicates(sample_info_allele, counts_allele)

# also filtering sample infos down to single info and 
# info_allele with all conditions
info <- sample_info |> select("condition", "experiment", 
                              "time_point_num", "time_point_str") |> 
  unique() |> arrange(experiment, time_point_num)
info_allele <- sample_info_allele |> select("condition", "experiment",
                                            "time_point_num", "time_point_str") |> 
  unique() |> arrange(experiment, time_point_num)
# ordering collapsed counts to same order as info
collapsed$cer <- collapsed$cer[,info$condition]
collapsed$par <- collapsed$par[,info$condition]
collapsed_allele$cer <- collapsed_allele$cer[,info_allele$condition]
collapsed_allele$par <- collapsed_allele$par[,info_allele$condition]

# pre-collapse
dim(counts[,sample_info$experiment == "LowN" & sample_info$allele == "cer"])
sample_info |> filter(experiment == "LowN" & allele == "cer") |> 
  select(time_point_num) |>
  table()
dim(counts[,sample_info$experiment == "HAP4" & sample_info$allele == "par"])
sample_info |> filter(experiment == "HAP4" & allele == "par") |> 
  select(time_point_num) |>
  table()
# post-collapse
dim(collapsed$cer[,info$experiment == "LowN"]) # LowN has been reduced (~30 replicates)
info |> filter(experiment == "LowN") |> select(time_point_num) |> table()
dim(collapsed$par[,info$experiment == "HAP4"]) # HAP4 has not been reduced (no replicates)
info |> filter(experiment == "HAP4") |> select(time_point_num) |> table()

#### Tests for collapsed/movavg counts ####
makeDf <- function(.cts, .info, .join_by = "sample_name") {
  outdf <- as_tibble(.cts) |>
    bind_cols(tibble(gene_name = rownames(.cts))) |>
    pivot_longer(cols = colnames(.cts),
                 names_to = .join_by, values_to = "expr") |>
    left_join(.info,
              by = .join_by)
  return(outdf)
}

# additional tests for collapsing replicates
# toy count matrix using random sample of genes (to test ability to separate low var genes out)
toy_mat <- counts[, sample_info$organism == "par" & sample_info$experiment == "LowN"]
toydf <- makeDf(toy_mat, sample_info)
# assessing distribution of variance, before collapsing replicates
plotdf0 <- toydf |> group_by(gene_name) |>
  summarise(var_expr = var(expr),
            mean_expr = mean(expr))
p0 <- ggplot(plotdf0, aes(x = log2(var_expr/mean_expr))) + geom_density() +
  geom_vline(xintercept = 0, color = "red") +
  ggtitle("before collapsing replicates")
# after collapsing replicates
toy_mat_collapsed <- collapsed$par[, info$experiment == "LowN"]
plotdf1 <- tibble(var_expr = apply(toy_mat_collapsed, 1, var),
                  mean_expr = apply(toy_mat_collapsed, 1, mean),
                  gene_name = rownames(toy_mat_collapsed))
p1 <- ggplot(plotdf1, aes(x = log2(var_expr/mean_expr))) + geom_density() +
  geom_vline(xintercept = 0, color = "red") +
  ggtitle("after collapsing replicates")
ggarrange(p0, p1, nrow = 1, ncol = 2)
# mean expression pre and post (shouldn't have changed)
plotdf <- left_join(plotdf0, plotdf1, by = "gene_name", suffix = c("_pre", "_post"))
ggplot(plotdf, aes(x = log2(mean_expr_pre), y = log2(mean_expr_post))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "gold")
# difference in var/mean pre and post versus mean expr
plotdf <- select(plotdf, var_expr_pre, var_expr_post, mean_expr_pre, gene_name) |>
  rename("mean_expr"="mean_expr_pre")
ggplot(plotdf, aes(x = log2(mean_expr + 1),
                   y = log2(var_expr_pre) - log2(var_expr_post))) +
  geom_hex() +
  geom_hline(yintercept = 0, color = "gold")
# as expected, there is lower variance post collapse (very few genes below y = 0)
# most genes have fairly little difference in var when collapsing replicates,
# but the difference can be extreme for a few genes, of middling expression level
# For genes with a difference, do we trust the collapsed variance to be more accurate to expression shape?
# first random gene that tend to have very little difference in variance
gene_idx <- plotdf |> select(gene_name) |> pull() |> sample(1)
ggplot(filter(toydf, gene_name == gene_idx), aes(x = time_point_str, y = log2(expr + 1))) +
  geom_jitter() +
  geom_line(data = summarise(group_by(filter(toydf, gene_name == gene_idx), time_point_str, gene_name), mean_expr = mean(expr)),
            aes(x = time_point_str, y = log2(mean_expr + 1), group = gene_name),
            color = "red")
# second a gene with variance substantially reduced in collapsed form
gene_idx <- plotdf |> filter(log2(var_expr_pre) - log2(var_expr_post) > 7) |> select(gene_name) |> pull() |> sample(1)
ggplot(filter(toydf, gene_name == gene_idx), aes(x = time_point_str, y = log2(expr + 1))) +
  geom_jitter() +
  geom_line(data = summarise(group_by(filter(toydf, gene_name == gene_idx), time_point_str, gene_name), mean_expr = mean(expr)),
            aes(x = time_point_str, y = log2(mean_expr + 1), group = gene_name),
            color = "red") # these are genes with very little expression change attributable to timepoint,
# the exact kind we want to have significantly reduced variance in collapsed counts

#### taking moving average of LowPi and HAP4 ####
movavg_LowPi_cer <- getMovingAverage(collapsed$cer[,info$experiment == "LowPi"])
movavg_LowPi_par <- getMovingAverage(collapsed$par[,info$experiment == "LowPi"])
movavg_LowPi_hyc <- getMovingAverage(collapsed_allele$cer[,info_allele$experiment == "LowPi"])
movavg_LowPi_hyp <- getMovingAverage(collapsed_allele$par[,info_allele$experiment == "LowPi"])
movavg_HAP4_cer <- getMovingAverage(collapsed$cer[,info$experiment == "HAP4"])
movavg_HAP4_par <- getMovingAverage(collapsed$par[,info$experiment == "HAP4"])
movavg_HAP4_hyc <- getMovingAverage(collapsed_allele$cer[,info_allele$experiment == "HAP4"])
movavg_HAP4_hyp <- getMovingAverage(collapsed_allele$par[,info_allele$experiment == "HAP4"])

# adding moving average counts to counts
# cer
collapsed$cer[, c(colnames(movavg_HAP4_cer),
                  colnames(movavg_LowPi_cer))] <- cbind(movavg_HAP4_cer, 
                                                        movavg_LowPi_cer)
# par
collapsed$par[, c(colnames(movavg_HAP4_par),
                  colnames(movavg_LowPi_par))] <- cbind(movavg_HAP4_par, 
                                                        movavg_LowPi_par)
# hyc
collapsed_allele$cer[, c(colnames(movavg_HAP4_hyc),
                         colnames(movavg_LowPi_hyc))] <- cbind(movavg_HAP4_hyc, 
                                                               movavg_LowPi_hyc)
# hyp
collapsed_allele$par[, c(colnames(movavg_HAP4_hyp),
                         colnames(movavg_LowPi_hyp))] <- cbind(movavg_HAP4_hyp, 
                                                               movavg_LowPi_hyp)

# tests for getMovingAverage
# plotting before and after moving average
gene_idx <- sample(rownames(counts), 1)
# change to try different experiments (you have to keep the same experiment/allele/species):
test_counts <- collapseReplicates(sample_info, counts)
test_collapsed <- test_counts$cer[,info$condition]
test_movavg <- collapsed$cer
test_info <- info
test_experiment <- "LowPi"

plotdf <- tibble(expr = test_collapsed[gene_idx, test_info$experiment == test_experiment],
                 condition = colnames(test_collapsed[, test_info$experiment == test_experiment])) |>
  left_join(test_info, by = "condition")
ggplot(plotdf, aes(x = time_point_num, y = expr)) + geom_line()

# same random gene after moving average
plotdf2 <- tibble(expr = test_movavg[gene_idx, info$experiment == test_experiment],
                  status = "after",
                  condition = colnames(test_movavg[,info$experiment == test_experiment])) |>
  left_join(test_info, by = "condition")
plotdf$status <- "before"
plotdf <- bind_rows(plotdf, plotdf2)
ggplot(plotdf, aes(x = time_point_num, y = expr)) +
  geom_line(aes(group = status, color = status))

# checking conditions order still matches info
sum(colnames(collapsed$cer) == info$condition)/nrow(info)
sum(colnames(collapsed$par) == info$condition)/nrow(info)
sum(colnames(collapsed_allele$cer) == info_allele$condition)/nrow(info_allele)
sum(colnames(collapsed_allele$par) == info_allele$condition)/nrow(info_allele)

# done with data wrangling, cleaning up extra variables
rm(movavg_LowPi_cer,
   movavg_LowPi_par,
   movavg_LowPi_hyc,
   movavg_LowPi_hyp,
   movavg_HAP4_cer,
   movavg_HAP4_par,
   movavg_HAP4_hyc,
   movavg_HAP4_hyp,
   plotdf, plotdf2,
   test_collapsed,
   test_info, test_movavg)

# saving
save(counts, sample_info, collapsed, info, 
     file = "data_files/Cleaned_Counts.RData")
save(counts_allele, sample_info_allele, collapsed_allele, info_allele,
     file = "data_files/Cleaned_Counts_Allele.RData")

#### Data for single gene GLMs ####
# using the cleaned counts, we now create data structures for heiarchical
# clustering and for generalized linear modeling. Both need lists of
# counts data with a separate entry for each experiment. Both need pre-filtering
# of genes that are lowly expressed in both species in that environment.
# But the GLM data needs to include replicates while the heiarchical clustering
# used collapsed counts, averaged among replicates or moving average if there are not replicates
load("data_files/Cleaned_Counts.RData")
load("data_files/Cleaned_Counts_Allele.RData")
# not using collapsed counts (replicates are vital for the glm)
rm(collapsed, info, collapsed_allele, info_allele)

# combining two species' counts for each experiment into the same dataset
# filtering down to conditions/replicates (ex. LowN_TP1_rep2) where both species have data
ExperimentNames <- c("LowN", "CC", "HAP4", "LowPi", "Heat", "Cold")
spcts <- vector(mode = "list", length = length(ExperimentNames))
names(spcts) <- ExperimentNames
spinfo <- vector(mode = "list", length = length(ExperimentNames))
names(spinfo) <- ExperimentNames
alcts <- vector(mode = "list", length = length(ExperimentNames))
names(alcts) <- ExperimentNames
alinfo <- vector(mode = "list", length = length(ExperimentNames))
names(alinfo) <- ExperimentNames
for(e in ExperimentNames) {
  cer_cts <- counts[, sample_info$experiment == e & 
                      sample_info$allele == "cer"]
  par_cts <- counts[, sample_info$experiment == e & 
                      sample_info$allele == "par"]
  hyc_cts <- counts_allele[, sample_info_allele$experiment == e & 
                             sample_info_allele$allele == "cer"]
  hyp_cts <- counts_allele[, sample_info_allele$experiment == e & 
                             sample_info_allele$allele == "par"]
  cer_info <- sample_info[sample_info$experiment == e & 
                            sample_info$organism == "cer",]
  par_info <- sample_info[sample_info$experiment == e & 
                            sample_info$organism == "par",]
  hyc_info <- sample_info_allele[sample_info_allele$experiment == e & 
                                   sample_info_allele$allele == "cer",]
  hyp_info <- sample_info_allele[sample_info_allele$experiment == e & 
                                   sample_info_allele$allele == "par",]
  # combining
  spcts[[e]] <- cbind(cer_cts, par_cts)
  spinfo[[e]] <- bind_rows(cer_info, par_info)
  alcts[[e]] <- cbind(hyc_cts, hyp_cts)
  alinfo[[e]] <- bind_rows(hyc_info, hyp_info)
  cat("number of times each condition is represented in", e, ":", 
      names(table(table(spinfo[[e]]$condition))), "\n") # number of replicates*2 species in each environment (LowPi and HAP4 don't have replicates)
}
# factorize genotype, experiment, and time_point_str (no specific reference level for the last two)
spinfo <- lapply(spinfo, mutate, genotype = as.factor(genotype) %>% relevel(ref = "WT"))
alinfo <- lapply(alinfo, mutate, genotype = as.factor(genotype) %>% relevel(ref = "WT"))
# releveling time_point_str is a tad more tricky b/c the reference is named differently in each experiment (and isn't 0 in LowPi)
spinfo <- lapply(spinfo, mutate, time_point_str = as.factor(time_point_str) %>% relevel(ref = time_point_str[which.min(parse_number(time_point_str))]))
spinfo <- lapply(spinfo, mutate, allele = as.factor(allele) %>% relevel(ref = "par")) # we always do cer/par for expression ratios, so we want a big estimate to also mean cer is more strongly expressed
alinfo <- lapply(alinfo, mutate, time_point_str = as.factor(time_point_str) %>% relevel(ref = time_point_str[which.min(parse_number(time_point_str))]))
alinfo <- lapply(alinfo, mutate, allele = as.factor(allele) %>% relevel(ref = "par"))

# removing genes that are lowly expressed in both parental species
expr_thresh <- 30
for(e in ExperimentNames) {
  spcts[[e]] <- spcts[[e]][apply(spcts[[e]], 1, \(x) !all(is.na(x))),]
  alcts[[e]] <- alcts[[e]][apply(alcts[[e]], 1, \(x) !all(is.na(x))),]
  good_idxs <- rowMeans(spcts[[e]]) > expr_thresh
  spcts[[e]] <- spcts[[e]][good_idxs,]
  alcts[[e]] <- alcts[[e]][good_idxs,]
}
# now datasets from different environments have different nGenes
dim(spcts$LowN)
dim(alcts$LowN)
dim(spcts$HAP4)
dim(alcts$HAP4)
dim(spcts$Cold)
dim(alcts$Cold)

# saving
save(spcts, spinfo, alcts, alinfo, file = "data_files/GLM_Counts.RData")

#### Data for Heiarchical Clustering ####
load("data_files/Cleaned_Counts.RData")
load("data_files/Cleaned_Counts_Allele.RData")
ExperimentNames <- c("LowN", "CC", "HAP4", "LowPi", "Heat", "Cold")

# parents only for clustering
# using info and collapsed counts from parents only, all WT samples, replicates averaged
dim(info)
dim(collapsed$cer)
dim(collapsed$par)
sum(grepl("WT", info$condition))/nrow(info) # should have all WT, ratio is 1
sum(grepl("WT", colnames(collapsed$cer)))/ncol(collapsed$cer) # should have all WT, ratio is 1
table(colnames(collapsed$cer)) |> table() # should all have one entry, 1 for each condition

### Splitting counts into species/experiments
counts_list <- vector(mode = "list", length = 0)
for (e in unique(info$experiment)) {
  for (a in c("cer", "par")) {
    cts <- collapsed[[a]][, info$experiment == e]
    cts <- cts[apply(cts, 1, \(x) !all(is.na(x))),] # removing na genes from Heat/Cold
    counts_list[[paste(a, e, sep = "_")]] <- cts
  }
}
rm(cts)
# removing genes that are lowly expressed in BOTH species
expr_thresh <- 30
for (e in ExperimentNames) {
  cer_name <- paste("cer", e, sep = "_")
  par_name <- paste("par", e, sep = "_")
  low_cer <- rowMeans(counts_list[[cer_name]]) <= expr_thresh
  low_par <- rowMeans(counts_list[[par_name]]) <= expr_thresh
  good_idxs <- !(low_cer & low_par)
  counts_list[[cer_name]] <- counts_list[[cer_name]][good_idxs,]
  counts_list[[par_name]] <- counts_list[[par_name]][good_idxs,]
}
# now experiments have different nGenes (but species from the same experiment will have the same nGenes)
dim(counts_list$cer_CC)
dim(counts_list$par_CC)
dim(counts_list$cer_Cold)
dim(counts_list$par_Cold)
dim(info)
# saving
save(counts_list, info, file = "data_files/Clustering_Counts.RData")

#### Categorizing level and dynamics of all genes ####
# colors used throughout paper
levdyn_colordf <- tibble(type = c("salmon", "aquamarine", "gold", "greenyellow"),
                         labels = c("conserved level \nand dynamics", 
                                    "conserved level, \ndiverged dynamics", 
                                    "diverged level, \nconserved dynamics", 
                                    "diverged level \nand dynamics"),
                         limits = c("conserved level and dynamics", 
                                    "conserved level, diverged dynamics", 
                                    "diverged level, conserved dynamics", 
                                    "diverged level and dynamics"))

# Running the following section requires running the clustering.R 
# and single_gene_model_construction_and_QC>R scripts first
load("data_files/CorrelationClustering3Disp.RData")
load("data_files/single_gene_models.RData")
load("data_files/Cleaned_Count_Data.RData")
rm(counts, sample_info)
ExperimentNames <- c("CC", "HAP4", "LowN", "LowPi", "Heat", "Cold")
p_thresh <- 0.05
eff_thresh <- log2(1.5) # gene needs to be 1.5x higher in one species than the other to be considered DE

# Benjamini-Hochberg FDR correction of glm p-values
spaldf$padj <- p.adjust(spaldf$pvalue, method = "BH")
sum(spaldf$pvalue < p_thresh)
sum(spaldf$padj < p_thresh)

# Creating mutually exclusive sets of genes based on
# level and dynamics divergence in parents
# each gene is a combination of its level category 
# and its dynamics category
### level categories:
# 1) conserved: pval >= p_thresh
# 2) diverged: parent pval < p_thresh, not taking into account hybrid yet
### dynamics categories (just looking at parents):
# 1) conserved: 11, 22, NANA
# 2) diverged: 12, 21
# 3) diverged in level: NA for one value
finaldf <- pivot_wider(spaldf, id_cols = c("gene_name", "experiment"), 
                       values_from = c("effect_size", "padj"),
                       names_from = "coefficient") |> 
  drop_na() |> # 1 gene missing from species, YMR107W
  inner_join(y =  clusterdf, 
             by = join_by("gene_name"=="gene_ID",
                          "experiment"))

finaldf$dynamics <- map2(finaldf$cer, finaldf$par, \(x, y) {
  if (x == y) {
    return("conserved")
  }
  if (x != y) {
    return("diverged")
  }
}) |> unlist()

finaldf$level <- map(c(1:nrow(finaldf)), \(i) {
  x <- finaldf$padj_species[i] |> as.numeric()
  y <- finaldf$effect_size_species[i] |> as.numeric()
  output <- if_else(x < p_thresh & abs(y) > eff_thresh,
                    true = "diverged",
                    false = "conserved")
  return(output)
}) |> unlist()

# for testing how eff_threshold affects level divergence
# finaldf <- finaldf |> pivot_longer(cols = c("level0", "level15", "level2"),
#                                    names_to = "eff_threshold",
#                                    values_to = "level")

table(finaldf$dynamics, finaldf$level, useNA = "always") # checking that all observations are in the 2x2 box
finaldf <- drop_na(finaldf)
sum(table(finaldf$level, finaldf$dynamics))
nrow(finaldf) # mutually exclusive group categories
finaldf |> filter(dynamics == "diverged") |> select(experiment) |> table()
finaldf |> filter(level == "diverged") |> select(experiment) |> table()

### Removing hybrid par allele-deleted genes from CC experiment (where the deletion was present)
omit_list <- c("YLR078C", "YLR077W", "YLR074C", "YLR072W", "YLR075W", "YLR073C", # large hybrid CC paradoxus haplotype deletion
               "YNL247W", "YNL244C")
finaldf |> filter(experiment == "CC" & gene_name %in% omit_list)
dim(finaldf)
finaldf <- finaldf |> filter(!(experiment == "CC" & gene_name %in% omit_list))
dim(finaldf)

### Checking allele mapping bias in individual genes
# adding bias column to finaldf
finaldf <- finaldf |> mutate(bias = if_else(gene_name %in% cer_biased_genes,
                                           true = "cer",
                                           false = if_else(gene_name %in% par_biased_genes,
                                                           true = "par",
                                                           false = if_else(gene_name %in% both_biased_genes,
                                                                           true = "both",
                                                                           false = "none")))) 
plotdf <- finaldf |> 
  filter(level == "diverged")
# effect of bias on level divergence
ggplot(plotdf, aes(x = effect_size_species)) +
  geom_density(aes(fill = bias), alpha = 0.5) 
# pretty clear how allele mapping bias affects lfc
### amending so that biased genes are not considered diverged in level
finaldf$level <- if_else(finaldf$bias == "none",
                         true = finaldf$level,
                         false = "biased")

# effect of bias on dynamics divergence
p_none <- plotdf |> filter(bias == "none") |> 
  select(cer, par) |> table()
p_cer <- plotdf |> filter(bias == "cer") |> 
  select(cer, par) |> table()
p_par <- plotdf |> filter(bias == "par") |> 
  select(cer, par) |> table()

round(p_none/sum(p_none), digits = 2)
round(p_cer/sum(p_cer), digits = 2)
round(p_par/sum(p_par), digits = 2)
# proportions are the same between the three groups

### numbers for workflow figure
nrow(finaldf)
n_per_experiment <- table(finaldf$experiment)
# min filtered out for low expression:
max(n_per_experiment)
(nrow(genedf) - max(n_per_experiment))/nrow(genedf)
# max filtered out for low expression:
min(n_per_experiment)
(nrow(genedf) - min(n_per_experiment))/nrow(genedf)
# percent low variance per experiment:
finaldf |> group_by(experiment) |> 
  summarise(n_cer_lowvar = sum(cer == 0),
            n_par_lowvar = sum(par == 0),
            n_genes = n()*2) |> 
  mutate(pct_low_var = (n_cer_lowvar + n_par_lowvar)/n_genes)

### Barplot of 4 divergence categories
finaldf$group4 <- map2(finaldf$level, finaldf$dynamics, \(l, d) {
  if (l != "diverged" & d == "conserved") {
    return("conserved level and dynamics")
  }
  if (l != "diverged" & d == "diverged") {
    return("conserved level, diverged dynamics")
  }
  if (l == "diverged" & d == "conserved") {
    return("diverged level, conserved dynamics")
  }
  if (l == "diverged" & d == "diverged") {
    return("diverged level and dynamics")
  }
}) |> unlist()

# # for testing how eff_threshold affects level divergence
# # bar plot of total number of level and dynamics divergers in each environment
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/bar.pdf",
#     width = 11, height = 5)
# ggplot(mutate(finaldf,
#               experiment = factor(experiment, 
#                                   levels = c("HAP4", "CC", "LowN", "LowPi", "Heat", "Cold"))), 
#        aes(x = group4)) +
#   geom_bar(aes(fill = group4)) +
#   geom_text(stat='count', aes(label = after_stat(count)), vjust=-1) +
#   scale_fill_discrete(limits =  levdyn_colordf$limits,
#                       type = levdyn_colordf$type) +
#   scale_x_discrete(limits = levdyn_colordf$limits,
#                    breaks = levdyn_colordf$limits,
#                    labels = levdyn_colordf$labels) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#         legend.position = "none") +
#   xlab("") +
#   ylim(c(0, 3500)) +
#   facet_grid(eff_threshold ~ experiment)
# dev.off()

### group codes
# gene groups factor in a) whether level is conserved, b) whether dynamics are conserved,
# and c) if level or dynamics are diverged, which direction that is happening in (i.e. up in Scer, or 2 in Scer and 1 in Spar)
# group code format: type_direction_cluster(s)
# type: <cons, lev, dyn, or levdyn>
# direction (of level divergence): <uppar, upcer, or empty>
# clusters: <0, 1, 2, 01, 02, etc.>
finaldf$group <- apply(finaldf, 1, \(x) {
  x_cer <- x["cer"] |> as.numeric()
  x_par <- x["par"] |> as.numeric()
  x_effect_size <- x["effect_size_species"] |> as.numeric()
  x_level <- x["level"] |> as.character()
  x_dynamics <- x["dynamics"] |> as.character()
  if (x_level != "diverged" & x_dynamics == "conserved") {
    type_str <- "cons"
  }
  if (x_level == "diverged" & x_dynamics == "conserved") {
    type_str <- "lev"
  }
  if (x_level != "diverged" & x_dynamics == "diverged") {
    type_str <- "dyn"
  }
  if (x_level == "diverged" & x_dynamics == "diverged") {
    type_str <- "levdyn"
  }
  if (x_level != "diverged") {
    dir_str <- NULL
  }
  if (x_level == "diverged") {
    if (x_effect_size > 0) {
      dir_str <- "upcer"
    }
    if (x_effect_size < 0) {
      dir_str <- "uppar"
    }
  }
  if (identical(x_cer, x_par)) {
    full_str <- paste0(type_str, dir_str, x_cer)
  }
  if (!identical(x_cer, x_par)) {
    full_str <- paste0(type_str, dir_str, x_cer, x_par)
  }
  return(full_str)
}) |> unlist()

# saving
save(finaldf, levdyn_colordf, file = "data_files/FinalDataframe3Disp.RData")

#### Gene Ontology Enrichment ####
# table of gene ontology terms associated with each gene
# each row is a term, so this is a very long table
goslim <- read.table("data_files/downloaded_genomes_and_features/go_slim_mapping.tab", header = FALSE, sep = "\t") |>
  as_tibble()
colnames(goslim) <- c("ORF", # (mandatory) - Systematic name of the gene (or gene complex if it starts with CPX-)
                      "gene", # (optional) - Gene name, if one exists
                      "SGDID", # SGDID (mandatory) - the SGDID, unique database identifier for the gene
                      "GO_aspect", # (mandatory) - which ontology: P=Process, F=Function, C=Component
                      "GOslim_term", # (mandatory) - the name of the GO term that was selected as a GO Slim term
                      "GOID", # (optional) - the unique numerical identifier of the GO term
                      "feature_type") # (mandatory) - a description of the sequence feature, such as ORF or tRNA
too_vague_terms <- c("cellular process", "molecular function", "biological process")
goslim <- filter(goslim, !(GOslim_term %in% too_vague_terms) & 
                   (ORF %in% finaldf$gene_name))
save(goslim, file = "data_files/GO_Slim.RData")

################################ Archive ############################
# #### Reading & Cleaning Data - Lab-specific datasets ####
# # read in count data (normalized) and sample info
# # cer, par
# load("data_files/Cleaned_Redhuis_Data.RData")
# # hyc, hyp
# load("data_files/Cleaned_Redhuis_Data_AlleleSpecific.RData")
# 
# #### optional: including Fay et al. 2023 temp data ####
# # cer, par
# load("data_files/Cleaned_Fay_Counts.RData")
# fay <- fay[rownames(fay) %in% rownames(counts),] # ~500 mostly lowly expressed genes filtered out
# # hyc hyp
# load("data_files/Cleaned_Fay_Counts_allele.RData")
# fay_allele <- fay_allele[rownames(fay_allele) %in% rownames(counts),]
# # there's no easy way to left_join matrices unfortunately
# missing_fay_genes <- rownames(counts)[!(rownames(counts) %in% rownames(fay))]
# missing_fay <- matrix(NA, nrow = length(missing_fay_genes), ncol = ncol(fay))
# rownames(missing_fay) <- missing_fay_genes
# colnames(missing_fay) <- colnames(fay)
# missing_fay_allele <- matrix(NA, nrow = length(missing_fay_genes), ncol = ncol(fay_allele))
# rownames(missing_fay_allele) <- missing_fay_genes
# colnames(missing_fay_allele) <- colnames(fay_allele)
# fay <- rbind(fay, missing_fay)
# fay_allele <- rbind(fay_allele, missing_fay_allele)
# fay <- fay[rownames(counts),]
# fay_allele <- fay_allele[rownames(counts),]
# sum(rownames(fay) == rownames(counts))/nrow(counts)
# sum(rownames(fay_allele) == rownames(counts_allele))/nrow(counts_allele) # cbind doesn't check that rownames match, so make sure they match beforehand
# # joining
# counts <- cbind(counts, fay)
# counts_allele <- cbind(counts_allele, fay_allele)
# sample_info <- bind_rows(sample_info, sample_info_fay) # bind rows thankfully does check that the colnames are the same
# sample_info_allele <- bind_rows(sample_info_allele, sample_info_fay_allele)
# # cleaning variables
# rm(fay, fay_allele, missing_fay, missing_fay_allele,
#    missing_fay_genes)
# #### TF deletion analysis ####
# if (!file.exists("data_files/TFdel.RData")) {
#   #load("data_files/Cleaned_Barkai_Data_Env_Specific.RData")
#   # library(MASS, include.only = "glm.nb")
#   # 
#   factorizeGenotypeAndTimepoint <- function(.info) {
#     if (length(unique(.info$genotype)) > 1) {
#       .info$genotype <- as.factor(.info$genotype) %>% relevel(ref = "WT")
#     }
#     if (length(unique(.info$time_point_str)) > 1) {
#       timepoints <- .info$time_point_str %>% unique()
#       reftime <- timepoints[which.min(parse_number(timepoints))] # can't just do 0 because LowPi has a -5 fml
#       .info$time_point_str <- as.factor(.info$time_point_str) %>% relevel(ref = reftime)
#     }
#     return(.info)
#   }
#   infos_TFdel <- lapply(infos_TFdel, factorizeGenotypeAndTimepoint)
#   # 
#   # verify that counts are counts per million
#   test_cts <- counts_all2$cer # change to par to change which set you're testing
#   test_rowIdx <- sample(c(1:nrow(test_cts)), 1)
#   test_colIdx <- sample(c(1:ncol(test_cts)), 1)
#   test_count <- test_cts[test_rowIdx, test_colIdx]
#   ((test_count/rowSums(test_cts, na.rm = TRUE)[test_rowIdx])*1e6) %>% round() # what it should be
#   test_cts[test_rowIdx, test_colIdx] # what it is currently, should be the same, +/- 1 for rounding
#   
#   # calcualting LFC, mean, and sd for one gene upon each TF deletion
#   # separately in each of 3 timepoints (3 * nTFs LFCs per gene)
#   # log2 fold change is calculated as:
#   # log2(mean(tfdel)/mean(wt)) = log2(m) - log2(m1)
#   # where m1 and sd1 are the mean and sd of counts from WT 
#   # and m2 and sd2 are the same for the TF deletion
#   # lfc > 0 means the gene increased expr upon deletion, lfc < 0 means it decreased)
#   getLFCandSDByGene <- function(.gene_idx, .cts, .info) {
#     gdf <- bind_cols(tibble(expr = .cts[,.gene_idx]), .info)
#     gdf$genotype <- factor(gdf$genotype) |> relevel(ref = "WT")
#     output <- expand_grid(genotype = setdiff(unique(gdf$genotype), "WT"), 
#                           time_point_str = unique(gdf$time_point_str))
#     # calculating LFC
#     m_wt <- gdf |> filter(genotype == "WT") |> select(expr) |> pull() |> mean()
#     lfcs <- map(unique(output$genotype), \(.g) {
#       m_del <- gdf |> filter(genotype == .g) |> select(expr) |> pull() |> mean()
#       if (m_wt < 100 & m_del < 100) {
#         lfc <- NA # excluding LFCs of lowly expressed genes because they can be quite high for a change in only ~30 raw counts
#       }
#       else {
#         lfc <- log2(m_del + 1)-log2(m_wt + 1) # adding a count of 1 prevents -Inf values
#       }
#       return(tibble(deletion = gsub("delete", "", .g),
#                     lfc = lfc))
#     }) |> reduce(.f = bind_rows)
#     output <- map2(output$genotype, output$time_point_str, \(.g, .t) {
#       gdf2gen1tp <- gdf |> filter(genotype %in% c("WT", .g) &
#                                     time_point_str == .t)
#       if (length(unique(gdf2gen1tp$genotype)) != 2) {
#         return(tibble(sd_wt = NA,
#                       mean_wt = NA,
#                       sd_del = NA,
#                       mean_del = NA,
#                       sig = NA,
#                       deletion = gsub("delete", "", .g),
#                       time_point_str = .t))
#       }
#       # calculating confidence intervals
#       wt_vec <- filter(gdf2gen1tp, genotype == "WT") |> select(expr) |> pull()
#       del_vec <- filter(gdf2gen1tp, genotype != "WT") |> select(expr) |> pull()
#       m1 <- mean(wt_vec)
#       sd1 <- sd(wt_vec)
#       m2 <- mean(del_vec)
#       sd2 <- sd(del_vec)
#       CI_wt <- c(m1 - sd1, m1 + sd1)
#       CI_del <- c(m2 - sd2, m2 + sd2)
#       sig <- all((all(CI_wt < m2) | all(CI_wt > m2)) &
#                    (all(CI_del < m1) | all(CI_del > m1)))
#       return(tibble(sd_wt = sd1,
#                     mean_wt = m1,
#                     sd_del = sd2,
#                     mean_del = m2,
#                     sig = sig,
#                     deletion = gsub("delete", "", .g),
#                     time_point_str = .t))
#     }) |> reduce(.f = bind_rows)
#     output$gene_name <- .gene_idx
#     output <- left_join(output, lfcs, by = "deletion", relationship = "many-to-one")
#     return(output)
#   }
#   
#   # tests for getLFCandSDByGene
#   # # first exploring detecting which TF deletions a gene is DE in with a positive control
#   # # example gene: YGR192C (TDH3) in GCR2 delete (GCR2 is direct positive regulator of TDH3)
#   # coef_thresh <- 0.5
#   # gene_idx <- "YGR192C"
#   # test <- getLFCandSDByGene(gene_idx, counts_TFdel$cer, infos_TFdel$cer)
#   # test <- bind_rows(test, tibble(gene_name = gene_idx,
#   #                                deletion = "WT",
#   #                                smd = NA,
#   #                                sd_wt = NA,
#   #                                mean_wt = NA,
#   #                                sd_del = NA,
#   #                                mean_del = NA,
#   #                                sig = NA))
#   # checkIfDE <- function(.lfcs, .sigs) {
#   #   output <- tibble(lfc = .lfcs, sig = .sigs) |> 
#   #     filter(sig)
#   #   # checking if all significant timepoints are sig in the same direction and of sufficient magnitude:
#   #   if (length(table(sign(output$lfc))) != 1) {
#   #     return(FALSE)
#   #   }
#   #   if (!all(abs(output$lfc) > coef_thresh)) {
#   #     return(FALSE)
#   #   }
#   #   return(sum(output$sig) >= 2)
#   # }
#   # test <- test |> group_by(gene_name, deletion) |> 
#   #   summarise(DE = checkIfDE(.lfcs = lfc, .sigs = sig))
#   # genedf <- bind_cols(tibble(expr = counts_TFdel$cer[,gene_idx]), infos_TFdel$cer) |>
#   #   mutate(gene_name = gene_idx,
#   #          deletion = gsub("delete", "", genotype))
#   # plotdf <- genedf |>
#   #   select(expr, gene_name, deletion, time_point_str) |>
#   #   left_join(test, by = c("gene_name", "deletion"),
#   #             relationship = "many-to-one")
#   # ggplot(plotdf, aes(x = deletion, y = expr)) +
#   #   geom_point(aes(color = DE, shape = time_point_str)) + theme_classic() +
#   #   theme(axis.text.x = element_text(angle = 90)) # PHD1 looks at first glance to be comparable to the 5 DE deletions, but it has high variation between 0h and 1h replicates
#   # 
#   # # test 2: example of why you can't take the mean LFC for each timepoint individually:
#   # # YLR053C in NRG1 delete in hyc has LFC 2.5 for the mean of each timepoint's LFC,
#   # # but LFC of 0.79 overall
#   # gene_idx <- "YLR053C"
#   # # overall lfc:
#   # counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "NRG1delete",
#   #                         "YLR053C"] |> mean()
#   # counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "WT",
#   #                         "YLR053C"] |> mean()
#   # log2(187) - log2(108)
#   # # just TP1:
#   # counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "NRG1delete" &
#   #                           infos_TFdel_allele$cer$time_point_num == 0,
#   #                         "YLR053C"] |> mean()
#   # counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "NRG1delete" &
#   #                           infos_TFdel_allele$cer$time_point_num == 0,
#   #                         "YLR053C"] |> sd()
#   # counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "WT" &
#   #                           infos_TFdel_allele$cer$time_point_num == 0,
#   #                         "YLR053C"] |> mean()
#   # counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "WT" &
#   #                           infos_TFdel_allele$cer$time_point_num == 0,
#   #                         "YLR053C"] |> sd()
#   # log2(33.5) - log2(3.5)
#   # 
#   # # just TP2:
#   # counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "NRG1delete" &
#   #                           infos_TFdel_allele$cer$time_point_num == 60,
#   #                         "YLR053C"] |> mean()
#   # counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "WT" &
#   #                           infos_TFdel_allele$cer$time_point_num == 60,
#   #                         "YLR053C"] |> mean()
#   # log2(32.5) - log2(2.7)
#   # # just TP3:
#   # counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "NRG1delete" &
#   #                           infos_TFdel_allele$cer$time_point_num == 960,
#   #                         "YLR053C"] |> mean()
#   # counts_TFdel_allele$cer[infos_TFdel_allele$cer$genotype == "WT" &
#   #                           infos_TFdel_allele$cer$time_point_num == 960,
#   #                         "YLR053C"] |> mean()
#   # log2(496) - log2(318.3)
#   # # mean of each single-TP's LFC:
#   # mean(c(3.25, 3.5, .64))
#   # # what getLFCandSDByGene says it is:
#   # test <- getLFCandSDByGene(gene_idx, counts_TFdel_allele$cer, infos_TFdel_allele$cer)
#   # test |> filter(deletion == "NRG1") |> select(lfc) |> pull() |> mean()
#   # 
#   
#   # applying to all genes in both species
#   TFdeldf_cer <- map_dfr(colnames(counts_TFdel$cer), \(g) {
#     cat("working on cer", g, which(colnames(counts_TFdel$cer) == g), "/", ncol(counts_TFdel$cer), "\n")
#     return(getLFCandSDByGene(g, counts_TFdel$cer, infos_TFdel$cer))
#   })
#   
#   TFdeldf_par <- map_dfr(colnames(counts_TFdel$par), \(g) {
#     cat("working on par", g, which(colnames(counts_TFdel$par) == g), "/", ncol(counts_TFdel$par), "\n")
#     return(getLFCandSDByGene(g, counts_TFdel$par, infos_TFdel$par))
#   })
#   
#   # also making them for hybrid alleles
#   TFdeldf_hyc <- map_dfr(colnames(counts_TFdel_allele$cer), \(g) {
#     cat("working on hyc", g, which(colnames(counts_TFdel_allele$cer) == g), "/", ncol(counts_TFdel_allele$cer), "\n")
#     return(getLFCandSDByGene(g, counts_TFdel_allele$cer, infos_TFdel_allele$cer))
#   })
#   
#   TFdeldf_hyp <- map_dfr(colnames(counts_TFdel_allele$par), \(g) {
#     cat("working on hyp", g, which(colnames(counts_TFdel_allele$par) == g), "/", ncol(counts_TFdel_allele$par), "\n")
#     return(getLFCandSDByGene(g, counts_TFdel_allele$par, infos_TFdel_allele$par))
#   })
#   
#   ### Collapsing TFdeldfs into one row per gene/TF combination
#   
#   # by default, a gene is considered sig DE upon deletion of a certain TF at a certain timepoint
#   # if it's mean expr in deletion is at least 1 sd away from the mean expr of the WT,
#   # but means and sds are included in un-collapsed TFdeldfs in case
#   # we want to change that definition to be more/less stringent
#   
#   coef_thresh <- 0.5
#   
#   # to collapse, we summarise each gene/TF combination as DE
#   # if they have all 3 timepoints in the same direction
#   # with non-overlapping confidence intervals (1 sd on either side of mean)
#   # and if all significant timepoints have LFCs of sufficient magnitude
#   # and in the same direction
#   checkIfDE <- function(.g, .lfc, .mean_wt_vec, .sd_wt_vec,
#                         .mean_tfdel_vec, .sd_tfdel_vec) {
#     # helper function that checks if the WT and TFdel means
#     # are separated by at least one standard deviation
#     # and at least one mean is greater than 50 cpm (same threshold as lfc)
#     checkCI <- function(m1, sd1, m2, sd2) {
#       CI_wt <- c(m1 - sd1, m1 + sd1)
#       CI_del <- c(m2 - sd2, m2 + sd2)
#       if (m1 < 50 & m2 < 50) {
#         sig <- FALSE
#       }
#       else {
#         sig <- all((all(CI_wt < m2) | all(CI_wt > m2)) &
#                      (all(CI_del < m1) | all(CI_del > m1)))
#       }
#       return(sig)
#     }
#     .g <- unique(.g)
#     cat("working on", .g, which(colnames(counts_TFdel$cer) == .g), "/", ncol(counts_TFdel$cer), "\n")
#     output <- tibble(mean_wt = .mean_wt_vec,
#                      sd_wt = .sd_wt_vec,
#                      mean_tfdel = .mean_tfdel_vec,
#                      sd_tfdel = .sd_tfdel_vec,
#                      meandiff = .mean_tfdel_vec - .mean_wt_vec)
#     output$sig <- apply(output, 1, \(x) {
#       checkCI(m1 = x["mean_wt"],
#               sd1 = x["sd_wt"],
#               m2 = x["mean_tfdel"],
#               sd2 = x["sd_tfdel"])
#     })
#     output <- output |>
#       filter(abs(.lfc) > coef_thresh & sig)
#     return(output)
#     # checking if all significant timepoints are sig in the same direction:
#     if (length(table(sign(output$meandiff))) != 1) {
#       return(FALSE)
#     }
#     return(nrow(output) == 3)
#   }
#   # tests for checkIfDE
#   test <- checkIfDE(.g = "YLR053C",
#                     .lfc = as.numeric(unique(TFdeldf_cer[TFdeldf_cer$gene_name == "YLR053C" &
#                                                            TFdeldf_cer$deletion == "GLN3", 
#                                                          "lfc"])),
#                     .mean_wt_vec = as.numeric(unlist(TFdeldf_cer[TFdeldf_cer$gene_name == "YLR053C" &
#                                                                    TFdeldf_cer$deletion == "GLN3", 
#                                                                  "mean_wt"])), 
#                     .sd_wt_vec = as.numeric(unlist(TFdeldf_cer[TFdeldf_cer$gene_name == "YLR053C" &
#                                                                  TFdeldf_cer$deletion == "GLN3", 
#                                                                "sd_wt"])),
#                     .mean_tfdel_vec = as.numeric(unlist(TFdeldf_cer[TFdeldf_cer$gene_name == "YLR053C" &
#                                                                       TFdeldf_cer$deletion == "GLN3", 
#                                                                     "mean_del"])), 
#                     .sd_tfdel_vec = as.numeric(unlist(TFdeldf_cer[TFdeldf_cer$gene_name == "YLR053C" &
#                                                                     TFdeldf_cer$deletion == "GLN3", 
#                                                                   "sd_del"])))
#   
#   # collapsing
#   # cer
#   TFdeldf_cer_collapsed <- TFdeldf_cer |> group_by(gene_name, deletion) |>
#     summarise(DE = checkIfDE(.g = gene_name,
#                              .mean_wt_vec = mean_wt, .sd_wt_vec = sd_wt,
#                              .mean_tfdel_vec = mean_del, .sd_tfdel_vec = sd_del),
#               lfc = unique(lfc),
#               mean_sd_wt = mean(sd_wt),
#               mean_sd_del = mean(sd_del))
#   # par
#   TFdeldf_par_collapsed <- TFdeldf_par |> group_by(gene_name, deletion) |>
#     summarise(DE = checkIfDE(.g = gene_name,
#                              .mean_wt_vec = mean_wt, .sd_wt_vec = sd_wt,
#                              .mean_tfdel_vec = mean_del, .sd_tfdel_vec = sd_del),
#               lfc = unique(lfc),
#               mean_sd_wt = mean(sd_wt),
#               mean_sd_del = mean(sd_del))
#   
#   # hyc
#   TFdeldf_hyc_collapsed <- TFdeldf_hyc |> group_by(gene_name, deletion) |>
#     summarise(DE = checkIfDE(.g = gene_name,
#                              .mean_wt_vec = mean_wt, .sd_wt_vec = sd_wt,
#                              .mean_tfdel_vec = mean_del, .sd_tfdel_vec = sd_del),
#               lfc = unique(lfc),
#               mean_sd_wt = mean(sd_wt),
#               mean_sd_del = mean(sd_del))
#   
#   # hyp
#   TFdeldf_hyp_collapsed <- TFdeldf_hyp |> group_by(gene_name, deletion) |>
#     summarise(DE = checkIfDE(.g = gene_name,
#                              .mean_wt_vec = mean_wt, .sd_wt_vec = sd_wt,
#                              .mean_tfdel_vec = mean_del, .sd_tfdel_vec = sd_del),
#               lfc = unique(lfc),
#               mean_sd_wt = mean(sd_wt),
#               mean_sd_del = mean(sd_del))
#   
#   # Looking at the TFs themselves
#   TFdel_lookup <- read_delim(file = "data_files/yeastract_46TFs.csv", delim = ";", col_names = FALSE) |> select(X1, X2)
#   colnames(TFdel_lookup) <- c("common", "systematic")
#   # dataframe with just information on TFs, not how their deletion affects other genes
#   TFdf <- TFdel_lookup |> rename("gene_name"="systematic") |> 
#     left_join(y = select(module_genedf25, gene_name,
#                          effect_size_CC_species, effect_size_LowN_species,
#                          effect_size_HAP4_species, effect_size_LowPi_species,
#                          pvalue_CC_species, pvalue_LowN_species,
#                          pvalue_HAP4_species, pvalue_LowPi_species), by = "gene_name") |> 
#     drop_na()
#   names(TFdf) <- gsub("_species", "", names(TFdf))
#   
#   # visualizing each TF by its significant effect sizes in each environment
#   TFbyexpdf <- TFdf |> pivot_longer(cols = c("effect_size_CC", "effect_size_LowN",
#                                              "effect_size_HAP4", "effect_size_LowPi"), names_to = "experiment",
#                                     values_to = "effect_size", names_prefix = "effect_size_")
#   TFbyexpdf$pvalue <- apply(TFbyexpdf, 1, \(x) {
#     if (x["experiment"] == "CC") {
#       return(as.numeric(x["pvalue_CC"]))
#     }
#     if (x["experiment"] == "LowN") {
#       return(as.numeric(x["pvalue_LowN"]))
#     }
#     if (x["experiment"] == "HAP4") {
#       return(as.numeric(x["pvalue_HAP4"]))
#     }
#     if (x["experiment"] == "LowPi") {
#       return(as.numeric(x["pvalue_LowPi"]))
#     }
#   })
#   TFbyexpdf <- select(TFbyexpdf, -c("pvalue_CC", "pvalue_LowN",
#                                     "pvalue_HAP4", "pvalue_LowPi")) |> 
#     mutate(adjusted_effect_size = if_else(pvalue < 1e-5,
#                                           true = effect_size,
#                                           false = 0))
#   TFdf$adj_mean_effect_size <- TFbyexpdf |> group_by(common) |> summarise(adj_mean_effect_size = mean(adjusted_effect_size)) |> pull()
#   coef_thresh <- 0.25
#   TFdf$direction <- map(TFdf$adj_mean_effect_size, \(x) {
#     if (abs(x) < coef_thresh) {
#       return("conserved")
#     }
#     if (x < -coef_thresh) {
#       return("up_par")
#     }
#     if (x > coef_thresh) {
#       return("up_cer")
#     }
#   }) |> unlist()
#   table(TFdf$direction)
#   
#   # # most genes are not affected by any TF deletion:
#   # sum(rowSums(TFdel_pval_cer[,-1] < 1e-5) != 0)
#   # sum(rowSums(TFdel_pval_par[,-1] < 1e-5) != 0)
#   # # but all 42 TF deletions affect at least one gene:
#   # sum(colSums(TFdel_pval_cer[,-1] < 1e-5) != 0)
#   # sum(colSums(TFdel_pval_cer[,-1] < 1e-5) != 0)
#   # # many (but not most) genes that are affected by one TFdel are affected by multiple:
#   # table(rowSums(TFdel_pval_cer[,-1] < 1e-5))
#   # table(rowSums(TFdel_pval_par[,-1] < 1e-5))
#   save(TFdel_lookup, TFdeldf_cer, TFdeldf_par,
#        TFdeldf_hyc, TFdeldf_hyp, 
#        TFdeldf_cer_collapsed, TFdeldf_par_collapsed, TFdeldf_hyc_collapsed, TFdeldf_hyp_collapsed,
#        TFdf, TFbyexpdf, file = "data_files/TFdel.RData")
# }
# load("data_files/TFdel.RData")
# # based on above pvals and coefficients, determining which genes are regulated 
# # (directly or indirectly, can't distinguish here) by which TFs
# 
# #### TF Regulation Enrichment Tests ####
# # common_TFs <- intersect(setdiff(infos_TFdel$cer$genotype, "WT"), 
# #                         setdiff(infos_TFdel$par$genotype, "WT")) |> unique() |> 
# #   gsub(pattern = "delete", replacement = "")
# p_thresh <- 1e-5
# coef_thresh <- 2
# if (!file.exists("data_files/TFreg.RData")) {
#   # construct regdf where each row is one gene that is regulated by one TF
#   # in 1 of 4 experiments (so max 4 entries per gene, maybe filtering to experiments where TF is DE)
#   # takes a little while for reasons that elude me, as we're really just compiling data that's already been computed elsewhere
#   regdf <- map_dfr(common_TFs, \(tf) {
#     cat("starting on", tf, "\n")
#     gene_idxs_cer <- TFdeldf_cer |> 
#       filter(deletion == tf & abs(coef) > coef_thresh & pval < p_thresh) |> 
#       select(gene_name) |> pull()
#     gene_idxs_par <- TFdeldf_par |> 
#       filter(deletion == tf & abs(coef) > coef_thresh & pval < p_thresh) |> 
#       select(gene_name) |> pull()
#     gene_idxs <- union(gene_idxs_cer, gene_idxs_par)
#     output <- tibble(gene_name = character(0),
#                      TF = character(0),
#                      experiment = character(0),
#                      effect_size = numeric(0),
#                      TF_relationship = character(0))
#     for (g in gene_idxs) {
#       for (e in unique(info$experiment)) {
#         if (e == "LowN") {
#           e <- "LowN"
#         }
#         mod <- spaldf |> 
#           filter(gene_name == g & experiment == e & coefficient == "species") |> 
#           select(effect_size, pvalue)
#         effsize <- if_else(mod$pvalue < p_thresh, true = mod$effect_size, false = 0)
#         # TF relationship, one of 8 catagories:
#         #   activated in both cer and par (gene is DE with negative effect size upon TF deletion)
#         #   repressed in both cer and par (gene is DE with positive effect size upon TF deletion)
#         #   activated only in cer
#         #   activated only in par
#         #   repressed only in cer
#         #   repressed only in par
#         #   activated in cer, repressed in par (unlikely)
#         #   repressed in cer, activated in par (unlikely)
#         #   9th, missing category: not affected in either species, which shouldn't be possible b/c we filtered these genes out earlier
#         tfmod_cer <- TFdeldf_cer |> filter(gene_name == g & deletion == tf) |> 
#           select(coef, pval)
#         tfmod_par <- TFdeldf_par |> filter(gene_name == g & deletion == tf) |> 
#           select(coef, pval)
#         if (tfmod_cer$pval >= p_thresh | tfmod_cer$coef == 0) {
#           tfrel_cer <- "noRelationship"
#         }
#         if (tfmod_cer$pval < p_thresh & tfmod_cer$coef > 0) {
#           tfrel_cer <- "repressed"
#         }
#         if (tfmod_cer$pval < p_thresh & tfmod_cer$coef < 0) {
#           tfrel_cer <- "activated"
#         }
#         if (tfmod_par$pval >= p_thresh | tfmod_par$coef == 0) {
#           tfrel_par <- "noRelationship"
#         }
#         if (tfmod_par$pval < p_thresh & tfmod_par$coef > 0) {
#           tfrel_par <- "repressed"
#         }
#         if (tfmod_par$pval < p_thresh & tfmod_par$coef < 0) {
#           tfrel_par <- "activated"
#         }
#         tfrel <- paste(tfrel_cer, tfrel_par, sep = "_")
#         output <- bind_rows(output, 
#                             tibble(gene_name = g,
#                                    TF = tf,
#                                    experiment = e,
#                                    effect_size = effsize,
#                                    TF_relationship = tfrel))
#       }
#     }
#     return(output)
#   })
#   # Adding TF expression divergence direction (up in cer/up in par)
#   regdf <- TFdf |> select(common, direction) |> rename("TF"="common") |> 
#     right_join(y = regdf, by = "TF")
#   
#   # checking for enriched TF-regulated gene groups
#   makeEnrichmentDataFrame <- function(.m_genedf) {
#     tfenrich_genedf <- TFdeldf_cer |> mutate(DE_cer = abs(coef) > coef_thresh & pval < p_thresh) |> 
#       select(gene_name, deletion, DE_cer)
#     tfenrich_genedf <- TFdeldf_par |> mutate(DE_par = abs(coef) > coef_thresh & pval < p_thresh) |> 
#       select(gene_name, deletion, DE_par) |> 
#       full_join(y = tfenrich_genedf, by = c("gene_name", "deletion"))
#     tfenrich_genedf <- tfenrich_genedf |> filter(DE_cer | DE_par)
#     tfenrich_genedf <- left_join(x = tfenrich_genedf, 
#                                  y = select(.m_genedf,
#                                             gene_name, cer_color, par_color),
#                                  by = "gene_name")
#     # collapsing genedf version to sum number of genes per module (cer color/par color combination)
#     tfenrichdf <- tfenrich_genedf |> group_by(deletion, cer_color, par_color) |> 
#       summarise(nDE_cer = sum(DE_cer, na.rm = TRUE),
#                 nDE_par = sum(DE_par, na.rm = TRUE))
#     # for one subset of one module that's regulated by one TF,
#     # conducts a Fisher Exact Test to determine if this is a significant enrichment
#     calculateEnrichment <- function(.enrich_row, .species = c("cer", "par"), .m_genedf) {
#       del <- .enrich_row["deletion"] |> as.character()
#       if (.species == "cer") {
#         a <- .enrich_row["nDE_cer"] |> as.numeric()
#         nDE <- ungroup(tfenrichdf) |> filter(deletion == del) |> select(nDE_cer) |> pull() |> as.numeric() |> sum(na.rm = TRUE)
#       }
#       if (.species == "par") {
#         a <- .enrich_row["nDE_par"] |> as.numeric()
#         nDE <- ungroup(tfenrichdf) |> filter(deletion == del) |> select(nDE_par) |> pull() |> as.numeric() |> sum(na.rm = TRUE)
#       }
#       if (a < 5) {
#         return(FALSE)
#       }
#       nGenes <- nrow(.m_genedf)
#       cercol <- .enrich_row["cer_color"] |> as.character()
#       parcol <- .enrich_row["par_color"] |> as.character()
#       modsize <- .m_genedf |> filter(cer_color == cercol & par_color == parcol) |> nrow()
#       contingencytab <- rbind(c(a, nDE - a),
#                               c(modsize - a, nGenes - nDE - modsize + a))
#       pval <- fisher.test(contingencytab, alternative = "greater")$p.value
#       return(pval < 0.05)
#     }
#     tfenrichdf$enriched_cer <- apply(tfenrichdf, 1, calculateEnrichment, .species = "cer", .m_genedf = .m_genedf)
#     tfenrichdf$enriched_par <- apply(tfenrichdf, 1, calculateEnrichment, .species = "par", .m_genedf = .m_genedf)
#     
#     tfenrichdf <- select(.m_genedf, cer_color, par_color, is_CCM) |> unique() |> 
#       right_join(y = tfenrichdf,
#                  by = c("cer_color", "par_color"))
#     block_size_table <- table(.m_genedf$cer_color, .m_genedf$par_color)
#     tfenrichdf$block_size <- map2(tfenrichdf$cer_color, tfenrichdf$par_color, \(x, y) {
#       return(block_size_table[x, y])
#     }) |> unlist()
#     return(tfenrichdf)
#   }
#   TFenrichdf00 <- makeEnrichmentDataFrame(module_genedf00)
#   TFenrichdf10 <- makeEnrichmentDataFrame(module_genedf10)
#   TFenrichdf25 <- makeEnrichmentDataFrame(module_genedf25)
#   TFenrichdf35 <- makeEnrichmentDataFrame(module_genedf35)
#   random_TFenrichdf00 <- makeEnrichmentDataFrame(random_module_genedf00)
#   random_TFenrichdf10 <- makeEnrichmentDataFrame(random_module_genedf10)
#   random_TFenrichdf25 <- makeEnrichmentDataFrame(random_module_genedf25)
#   random_TFenrichdf35 <- makeEnrichmentDataFrame(random_module_genedf35)
#   
#   save(regdf, TFenrichdf00, TFenrichdf10, TFenrichdf25, TFenrichdf35,
#        random_TFenrichdf00, random_TFenrichdf10, random_TFenrichdf25, random_TFenrichdf35,
#        file = "data_files/TFreg.RData")
# }
# load("data_files/TFreg.RData")
