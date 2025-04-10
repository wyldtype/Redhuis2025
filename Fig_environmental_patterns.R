sapply(c("dplyr", "tidyr", "purrr", "WGCNA", "stringr", "ggplot2", "RColorBrewer", "ggpubr", "grid", "ComplexHeatmap", "circlize", "huxtable"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/aligning_the_molecular_phenotype/Redhuis2025/")
source("functions_for_figure_scripts.R")
load("data_files/FinalDataframe3Disp.RData")
load("data_files/Cleaned_Counts.RData")
load("data_files/Cleaned_Counts_Allele.RData")
load(file = "data_files/GO_Slim.RData")

# Goal: Here is the place to link findings from each individual environment
# to identify general patterns. The main way to convince people they are general
# patterns is to categorize every gene into one group (conserved, cis/level divergent,
# trans/dynamics divergent, or other) and show how "other" doesn't encompass many genes

#### Table of cluster counts separately for each species ####
ht <- finaldf |> select(cer, par, experiment) |> 
  pivot_longer(cols = c("cer", "par"),
               names_to = "organism", values_to = "cluster") |> 
  mutate(experiment = factor(experiment, levels = c("HAP4", "CC", "LowN", "LowPi", "Heat", "Cold"))) |> 
  arrange(experiment, cluster, organism) |> 
  group_by(experiment, cluster, organism) |> 
  summarise(n = n()) |> 
  mutate(title = paste0(experiment, cluster)) |> 
  select(title, organism, n) |> 
  pivot_wider(id_cols = "title", 
              names_from = "organism",
              values_from = "n") |> 
    hux(add_colnames = FALSE) |> 
  t() |> 
  theme_article()

# bold(ht)[1,]           <- TRUE
# bottom_border(ht)[1,]  <- 0.4
# align(ht)[,2]          <- "right"
# right_padding(ht)      <- 10
# left_padding(ht)       <- 10
# width(ht)              <- 1
# number_format(ht)      <- 0

quick_html(ht, file = "../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/cluster_table.html")

#### Visualize cluster expression patterns in each experiment ####
display.brewer.all()
load("data_files/CorrelationClustering3Disp.RData")
# color palettes we'll use for the 6 environments
palettedf <- tibble(experiment = ExperimentNames,
                    palette = c("Greys", "YlGn", "Greens", "Purples", "BuPu", "YlOrBr"),
                    long_name = LongExperimentNames)

plotClusterPatternByExperiment <- function(.df, .experiment, .title = NULL) {
  plotdf <- summarise(group_by(.df, time_point_num, label),
                      mean_expr = mean(expr, na.rm = TRUE))
  plotdf$label <- as.factor(plotdf$label)
  color_plt <- palettedf |> filter(experiment == .experiment) |> select(palette) |> pull()
  if (is.null(.title)) {
    .title <- palettedf |> filter(experiment == .experiment) |> select(long_name) |> pull()
  }
  ggplot(plotdf, aes(x = time_point_num, y = log2(mean_expr + 1))) +
    geom_line(aes(group = label, color = label), linewidth = 4) +
    geom_point(color = "black", alpha = 0.4) +
    xlab("timepoint (min)") +
    ylab("expression (log2)") +
    scale_color_brewer(palette = color_plt, name = "cluster",
                       direction = -1) +
    theme_classic() +
    theme(legend.position = "right") +
    ggtitle(.title)
}
pdf("../../aligning_the_molecular_phenotype/paper_figures/ExperimentOverview/cluster_ref.pdf",
    width = 8, height = 4)
ggarrange(plotClusterPatternByExperiment(clusterdf_list$HAP4_2$df, .experiment = "HAP4"),
          plotClusterPatternByExperiment(clusterdf_list$CC_2$df, .experiment = "CC"),
          plotClusterPatternByExperiment(clusterdf_list$LowN_2$df, .experiment = "LowN"),
          plotClusterPatternByExperiment(clusterdf_list$LowPi_2$df, .experiment = "LowPi"),
          plotClusterPatternByExperiment(clusterdf_list$Heat_2$df, .experiment = "Heat"),
          plotClusterPatternByExperiment(clusterdf_list$Cold_2$df, .experiment = "Cold"),
          nrow = 2, ncol = 3, common.legend = FALSE)
dev.off()

# pdf("../../aligning_the_molecular_phenotype/paper_figures/ExperimentOverview/cluster_ref.pdf",
#     width = 15, height = 2)
# ggarrange(plotClusterPatternByExperiment(clusterdf_list$HAP4_2$df, .experiment = "HAP4"),
#           plotClusterPatternByExperiment(clusterdf_list$CC_2$df, .experiment = "CC"),
#           plotClusterPatternByExperiment(clusterdf_list$LowN_2$df, .experiment = "LowN"),
#           plotClusterPatternByExperiment(clusterdf_list$LowPi_2$df, .experiment = "LowPi"),
#           plotClusterPatternByExperiment(clusterdf_list$Heat_2$df, .experiment = "Heat"),
#           plotClusterPatternByExperiment(clusterdf_list$Cold_2$df, .experiment = "Cold"),
#           nrow = 1, ncol = 6, common.legend = FALSE)
# dev.off()

# most common Scer shape
plotdf <- expand_grid(experiment = unique(finaldf$experiment),
            gene_name = unique(finaldf$gene_name)) |>
  left_join(y = finaldf, by = c("gene_name","experiment")) |>
  arrange(experiment) |>
  mutate(cer_clust = paste0(experiment, cer),
         par_clust = paste0(experiment, par)) |>
  ungroup() |>
  group_by(gene_name) |>
  summarise(code_cer = reduce(cer_clust, .f = paste, sep = "_"),
            code_par = reduce(par_clust, .f = paste, sep = "_"))
sort(table(plotdf$code_cer), decreasing = TRUE)[1:30]
sort(table(plotdf$code_par), decreasing = TRUE)[1:30]

#### Dynamics-divergers have higher expression in decreasing species at TP0 ####
# Rationale: In order for genes to have switched from increasing to decreasing,
# they must have higher or smaller expression difference at TP0

# example gene first, 2-1 HAP4
gene_idx <- finaldf |> filter(experiment == "HAP4" & 
                                level == "conserved" &
                                cer == 2 & par == 1) |> 
  select(gene_name) |> pull() |> sample(1)
plotEnvironments(.gene_idxs = "YMR221C") # first example we found, strong CC level divergence
plotEnvironments(.gene_idxs = gene_idx)

# Figure: dynamics-diverging genes have higher expression at TP0 in decreasing species
plotdf <- finaldf |> filter(dynamics == "diverged" &
                              level == "conserved" & 
                              cer %in% c(1, 2) & 
                              par %in% c(1, 2)) |> 
  mutate("decreasing_species" = if_else(cer == 1,
                                        true = "par",
                                        false = "cer")) |> 
  select(gene_name, experiment, cer, par, decreasing_species)
# ...in environment they were diverging in dynamics in:
getTP0AvgExpr <- function(.gene_name, .experiment, .organism) {
  if (.organism == "cer") {
    cts_mat <- collapsed$cer
    info_df <- info
  }
  if (.organism == "par") {
    cts_mat <- collapsed$par
    info_df <- info
  }
  tp0_condition <- info_df |> filter(experiment == .experiment) |> 
    filter(time_point_num == min(time_point_num)) |> 
    select(condition) |> pull()
  mean_expr <- cts_mat[.gene_name, tp0_condition] |> mean()
  return(mean_expr)
}
plotdf$tp0_avgExpr_decreasing_species <- map(c(1:nrow(plotdf)), .f = \(i) {
  getTP0AvgExpr(plotdf$gene_name[i],
                .experiment = plotdf$experiment[i],
                .organism = plotdf$decreasing_species[i])
}) |> unlist()
plotdf$tp0_avgExpr_increasing_species <- map(c(1:nrow(plotdf)), .f = \(i) {
  getTP0AvgExpr(plotdf$gene_name[i],
                .experiment = plotdf$experiment[i],
                .organism = setdiff(c("cer", "par"), plotdf$decreasing_species[i]))
}) |> unlist()
plotdf$lower_in_increasing_species <- plotdf$tp0_avgExpr_increasing_species < plotdf$tp0_avgExpr_decreasing_species
# ...in all environments:
p_tp0 <- ggplot(pivot_longer(plotdf, cols = c("tp0_avgExpr_decreasing_species",
                                     "tp0_avgExpr_increasing_species"),
                    names_to = "decreasing_or_increasing",
                    values_to = "tp0_avgExpr"), 
       aes(x = decreasing_or_increasing, 
           y = log2(tp0_avgExpr))) +
  geom_line(aes(group = gene_name,
                color = lower_in_increasing_species)) +
  facet_grid(~experiment + lower_in_increasing_species) +
  theme(axis.text.x = element_blank()) +
  ggtitle("decreasing species has higher expression at timepoint 0\n(in all environments except Cold and Cell Cycle)")
# True except for CC and Cold (where clusters 1 and 2 don't separate increasing/decreasing genes)

# ...level divergence in CC:
plotdf <- left_join(plotdf, y = finaldf |> 
                      filter(experiment == "CC") |> 
                      select(gene_name, effect_size_species, level),
                    by = "gene_name", relationship = "many-to-one")
p_cc <- ggplot(plotdf, aes(x = lower_in_increasing_species, 
                   y = effect_size_species)) +
  geom_boxplot(aes(color = decreasing_species)) +
  facet_wrap(~experiment) +
  ggtitle("decreasing species has higher expression level across Cell Cycle \n(in all environments except Cold and Cell Cycle)")
### Supplementary figure: decreasing species has higher expression in standard lab conditions
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/decreasing_vs_increasing_dynamics_at_tp0.pdf",
    width = 12, height = 8)
ggarrange(p_tp0, p_cc, nrow = 2, ncol = 1)
dev.off()
# TODO: if you end up using this, add counts above each barplot and line group

#### Similar proportions of level/dynamics divergence in all 6 environments ####
# Stacked bars of 4 groups in each environment
plotdf <- finaldf |> select(gene_name, experiment, group4) |> 
  pivot_wider(id_cols = c("gene_name"), names_from = "experiment",
              values_from = "group4") |> 
  pivot_longer(cols = ExperimentNames, names_to = "experiment",
               values_to = "group5") # to add NA values for low expression
plotdf$group5 <- if_else(is.na(plotdf$group5), true = "lowly expressed",
                         false = plotdf$group5) |> 
  factor(levels = c("conserved level and dynamics",
                    "conserved level, diverged dynamics",
                    "diverged level, conserved dynamics",
                    "diverged level and dynamics",
                    "lowly expressed"))
plotdf$experiment <- factor(plotdf$experiment, levels = ExperimentNames)
pdf("../paper_figures/EnvironmentalPatterns/stacked_bars.pdf",
    width = 6, height = 4)
ggplot(plotdf, aes(x = experiment)) + 
  geom_bar(aes(fill = group5)) +
  scale_fill_discrete(type = levdyn_colordf$type,
                      limits = levdyn_colordf$limits) +
  scale_x_discrete(limits = ExperimentNames, labels = LongExperimentNames) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("number of genes") +
  xlab("")
dev.off()
# Ridgeline plot of log2 fold change in each environment
library(ggridges)
plotdf <- filter(finaldf, level == "diverged")
filter(finaldf, dynamics == "conserved" & level == "conserved") |> group_by(experiment) |> 
  summarise(n = n())
# accompanying percentages of upScer and upSpar genes
pctdf <- plotdf |> group_by(experiment) |> 
  summarise(pct_upcer = round(sum(effect_size_species > 0)/length(effect_size_species), digits = 2)*100,
            pct_uppar = round(sum(effect_size_species < 0)/length(effect_size_species), digits = 2)*100,
            n = length(effect_size_species)) |> 
  pivot_longer(cols = c("pct_upcer", "pct_uppar"),
               names_to = "direction", values_to = "pct") |> 
  mutate(x_pos = if_else(direction == "pct_upcer",
                         true = 5, false = -5))
pctdf$y_pos <- sapply(pctdf$experiment, \(e) which(ExperimentNames == e) + 0.5)
pdf("../paper_figures/EnvironmentalPatterns/ridgeline.pdf",
    width = 4, height = 3)
ggplot(data = plotdf, aes(x = effect_size_species, y = experiment, fill = experiment)) +
  geom_density_ridges(data = plotdf) +
  geom_text(data = pctdf, aes(x = x_pos, y = y_pos, label = paste0(pct, "%"))) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("log2 fold change") +
  scale_y_discrete(limits = ExperimentNames, labels = LongExperimentNames)
dev.off()

# Heatmap of species-specific and reversals of expression plasticity
plotdf <- finaldf |> # filter(dynamics == "diverged") |> 
  group_by(experiment, cer, par) |> 
  summarise(n = n())
plotdf$long_experiment <- map(plotdf$experiment, \(e) {
  LongExperimentNames[which(ExperimentNames == e)]
}) |> factor(levels = LongExperimentNames)
plotdf$cer <- factor(plotdf$cer, levels = c("0", "1", "2"))
plotdf$par <- factor(plotdf$par, levels = c("2", "1", "0"))
plotdf$type <- map2(plotdf$cer, plotdf$par, \(x, y) {
  if (x == y) {
    return("conserved")
  }
  if (x == 0) {
    return("Spar-unique")
  }
  if (y == 0) {
    return("Scer-unique")
  }
  else {
    return("reversal")
  }
}) |> unlist()
pdf("../paper_figures/EnvironmentalPatterns/dynamics_counts_heatmap.pdf",
    width = 4, height = 3)
ggplot(plotdf, aes(x = cer, y = par, fill = type)) +
  geom_tile() +
  geom_text(aes(label = n,
                color = type == "conserved")) +
  scale_color_discrete(limits = c(TRUE, FALSE),
                       type = c("grey40", "white")) +
  scale_fill_discrete(limits = c("conserved", "Scer-unique", "Spar-unique", "reversal"),
                      type = c("grey60", "orange1", "blue2", "purple")) +
   facet_wrap(~long_experiment) +
  theme_classic() +
  xlab("Scer dynamics cluster") +
  ylab("Spar dynamics cluster") +
  guides(fill=guide_legend(title="number of genes")) +
  theme(legend.position = "none")
dev.off()

#### Dynamics divergence is unique to each environment #### 

# Y axis: environment those 2-1 or 1-2 divergers were ID'd in
# X axis: how their dynamics divergence looks in other environments

# How do we measure dynamics divergence? Correlation of avg expr between species

# Given gene idxs and experiment, 
# returns between-species correlation
# of average expression across timepoints
getDynamicsCorr <- function(.gene_idxs, .experiment_name) {
  condition_vec <- info |> filter(experiment == .experiment_name) |> 
    select(condition) |> pull()
  cer_vec <- collapsed$cer[.gene_idxs, condition_vec] |> 
    colMeans(na.rm = TRUE)
  par_vec <- collapsed$par[.gene_idxs, condition_vec] |> 
    colMeans(na.rm = TRUE)
  return(cor(cer_vec, par_vec))
}
# tests for getDynamics Corr
# HAP4 2-1, should be very uncorrelated in Sat Growth
gene_idxs <- finaldf |> filter(experiment == "HAP4" &
                                 cer == 2 & par == 1) |> 
  select(gene_name) |> pull()
getDynamicsCorr(gene_idxs, .experiment_name = "HAP4")
getDynamicsCorr(gene_idxs, .experiment_name = "LowN")
getDynamicsCorr(gene_idxs, .experiment_name = "LowPi")
getDynamicsCorr(gene_idxs, .experiment_name = "CC")
getDynamicsCorr(gene_idxs, .experiment_name = "Heat")
getDynamicsCorr(gene_idxs, .experiment_name = "Cold")

# 2-1 divergers
plot_mat <- matrix(nrow = length(ExperimentNames),
                   ncol = length(ExperimentNames))
for (e_row in ExperimentNames) {
  for (e_col in ExperimentNames) {
    e_gene_idxs <- finaldf |> filter(experiment == e_row &
                                       cer == 2 & par == 1) |> 
      select(gene_name) |> pull()
    e_cor <- getDynamicsCorr(e_gene_idxs, .experiment_name = e_col)
    plot_mat[which(ExperimentNames == e_row),
             which(ExperimentNames == e_col)] <- e_cor
  }
}
colnames(plot_mat) <- LongExperimentNames
rownames(plot_mat) <- LongExperimentNames

# plotting
col_fun = colorRamp2(c(-1, 0, 1), c("red2", "white", "skyblue"))
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/heatmap_21.pdf",
    width = 5, height = 3)
Heatmap(plot_mat, col = col_fun,
        row_order = LongExperimentNames, column_order = LongExperimentNames,
        row_names_side = "left", heatmap_legend_param = list(title = ""),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", plot_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

# 1-2 divergers
plot_mat <- matrix(nrow = length(ExperimentNames),
                   ncol = length(ExperimentNames))
for (e_row in ExperimentNames) {
  for (e_col in ExperimentNames) {
    e_gene_idxs <- finaldf |> filter(experiment == e_row &
                                       cer == 1 & par == 2) |> 
      select(gene_name) |> pull()
    e_cor <- getDynamicsCorr(e_gene_idxs, .experiment_name = e_col)
    plot_mat[which(ExperimentNames == e_row),
             which(ExperimentNames == e_col)] <- e_cor
  }
}
colnames(plot_mat) <- LongExperimentNames
rownames(plot_mat) <- LongExperimentNames

# plotting
col_fun = colorRamp2(c(-1, 0, 1), c("red2", "white", "skyblue"))
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/heatmap_12.pdf",
    width = 5, height = 3)
Heatmap(plot_mat, col = col_fun,
        row_order = LongExperimentNames, column_order = LongExperimentNames,
        row_names_side = "left", heatmap_legend_param = list(title = ""),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", plot_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

### All dynamics divergers
getDynamicsCorrAllClusters <- function(.gene_idxs, .experiment_name) {
  condition_vec <- info |> filter(experiment == .experiment_name) |> 
    select(condition) |> pull()
  
  clust_pairs <- finaldf |> filter(gene_name %in% .gene_idxs &
                                     experiment == .experiment_name) |> 
    group_by(cer, par) |> summarise(n = n()) |> ungroup()
  cors <- map2(clust_pairs$cer, clust_pairs$par, \(x, y) {
    clust_idxs <- finaldf |> filter(gene_name %in% .gene_idxs &
                                      experiment == .experiment_name &
                                      cer == x & par == y) |> 
      select(gene_name) |> pull()
    cer_vec <- collapsed$cer[clust_idxs, condition_vec, drop = FALSE] |> 
      colMeans(na.rm = TRUE)
    par_vec <- collapsed$par[clust_idxs, condition_vec, drop = FALSE] |> 
      colMeans(na.rm = TRUE)
    return(cor(cer_vec, par_vec))
  }) |> unlist()
  return(as.numeric(sum(cors*clust_pairs$n)/sum(clust_pairs$n))) # weighted average
}

# tests for getDynamicsCorrBoth
# LowPi 1-2 alone:
gene_idxs <- finaldf |> filter(experiment == "LowPi" &
                                 dynamics == "diverged" &
                                 cer == 1 & par == 2) |> 
  select(gene_name) |> pull()
testcor12 <- getDynamicsCorrAllClusters(gene_idxs, .experiment_name = "LowPi")
testcor12
getDynamicsCorr(gene_idxs, .experiment_name = "LowPi") # should be same number
# LowPi 2-1 alone:
gene_idxs <- finaldf |> filter(experiment == "LowPi" &
                                 dynamics == "diverged" &
                                 cer == 2 & par == 1) |> 
  select(gene_name) |> pull()
testcor21 <- getDynamicsCorrAllClusters(gene_idxs, .experiment_name = "LowPi")
testcor21
getDynamicsCorr(gene_idxs, .experiment_name = "LowPi") # should be same number
# weights
finaldf |> filter(experiment == "LowPi" &
                    dynamics == "diverged" &
                    cer != 0 & par != 0) |> 
  select(cer, par) |> table()
# what both should be:
sum(testcor12*372, testcor21*795)/(372 + 795) # with weighted average
mean(c(testcor12, testcor21)) # without weighted average
# Both 1-2 and 2-1:
gene_idxs <- finaldf |> filter(experiment == "LowPi" &
                                 dynamics == "diverged" &
                                 cer != 0 & par != 0) |> 
  select(gene_name) |> pull()
getDynamicsCorrAllClusters(gene_idxs, .experiment_name = "LowPi") # should be same as manual calculation
# test2: HAP4/CC failed until we added drop = FALSE
gene_idxs <- finaldf |> filter(experiment == "HAP4" &
                                 dynamics == "diverged" &
                                 cer != 0 & par != 0) |> 
  select(gene_name) |> pull()
finaldf |> filter(experiment == "CC" & gene_name %in% gene_idxs) |> 
  select(cer, par) |> table() # b/c there's only one 1-0 gene
getDynamicsCorrAllClusters(gene_idxs, .experiment_name = "CC")

# plotting
plot_mat <- matrix(nrow = length(ExperimentNames),
                   ncol = length(ExperimentNames))
for (e_row in ExperimentNames) {
  for (e_col in ExperimentNames) {
    cat(e_row, e_col, "\n")
    e_gene_idxs <- finaldf |> filter(experiment == e_row &
                                         cer != 0 & par != 0 &
                                       dynamics == "diverged") |> 
      select(gene_name) |> pull()
    e_cor <- getDynamicsCorrAllClusters(.gene_idxs = e_gene_idxs,
                                        .experiment_name = e_col)
    plot_mat[which(ExperimentNames == e_row),
             which(ExperimentNames == e_col)] <- e_cor
  }
}
colnames(plot_mat) <- LongExperimentNames
rownames(plot_mat) <- LongExperimentNames

# plotting
col_fun = colorRamp2(c(-1, 0, 1), c("red", "white", "skyblue"))
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/heatmap_dynamics_all_clusters.pdf",
    width = 5, height = 4)
Heatmap(plot_mat, col = col_fun,
        row_order = LongExperimentNames, column_order = LongExperimentNames,
        row_names_side = "left", heatmap_legend_param = list(title = ""),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", plot_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

### Supplementary figure: Diagonal is still the most uncorrelated when taking average cor of each ortholog pair
# (as opposed to one cor of average expr of all divergent genes)
# But now the correlations are much weaker magnitude b/c most single genes are weakly correlated moving the average towards 0
getDynamicsCorrPerGene <- function(.gene_idxs, .experiment_name) {
  # filter for genes expressed high enough in this environment
  good_expr_genes <- finaldf |> filter(experiment == .experiment_name) |> 
    select(gene_name) |> pull()
  .gene_idxs <- .gene_idxs[.gene_idxs %in% good_expr_genes]
  # calculating correlation of each gene btwn parental orthologs
  condition_vec <- info |> filter(experiment == .experiment_name) |> 
    select(condition) |> pull()
  cer_mat <- collapsed$cer[.gene_idxs, condition_vec]
  par_mat <- collapsed$par[.gene_idxs, condition_vec]
  cor_vec <- map(.gene_idxs, \(g) {
    return(cor(cer_mat[g,], par_mat[g,]))
  }) |> unlist()
  return(mean(cor_vec, na.rm = TRUE)) # NAs can arise from one species' counts being all 0s
}
# tests for getDynamicsCorrPerGene
# HAP4 all diverged dynamics
gene_idxs <- finaldf |> filter(experiment == "HAP4" &
                                 dynamics == "diverged" & 
                                 cer != 0 & par != 0) |> 
  select(gene_name) |> pull()
getDynamicsCorrPerGene(gene_idxs, .experiment_name = "HAP4")
getDynamicsCorrPerGene(gene_idxs, .experiment_name = "LowN")
getDynamicsCorrPerGene(gene_idxs, .experiment_name = "LowPi")
getDynamicsCorrPerGene(gene_idxs, .experiment_name = "CC")
getDynamicsCorrPerGene(gene_idxs, .experiment_name = "Heat")
getDynamicsCorrPerGene(gene_idxs, .experiment_name = "Cold")

# all dynamics divergers
plot_mat <- matrix(nrow = length(ExperimentNames),
                   ncol = length(ExperimentNames))
for (e_row in ExperimentNames) {
  for (e_col in ExperimentNames) {
    e_gene_idxs <- finaldf |> filter(experiment == e_row &
                                       dynamics == "diverged" & 
                                       cer != 0 & par != 0) |> 
      select(gene_name) |> pull()
    e_cor <- getDynamicsCorrPerGene(e_gene_idxs, .experiment_name = e_col)
    plot_mat[which(ExperimentNames == e_row),
             which(ExperimentNames == e_col)] <- e_cor
  }
}
colnames(plot_mat) <- LongExperimentNames
rownames(plot_mat) <- LongExperimentNames

# plotting
col_fun = colorRamp2(c(0, 0.25, 0.5), c("red2", "white", "skyblue"))
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/dynamics_divergence_heatmap_per_gene.pdf",
    width = 5, height = 3)
Heatmap(plot_mat, col = col_fun,
        row_order = LongExperimentNames, column_order = LongExperimentNames,
        row_names_side = "left", heatmap_legend_param = list(title = ""),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", plot_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

#### Supplemental figure: 6x6 plots of all expression divergence ####
# instead or in addition to the heatmap, have a 6x6 grid of 
# expression profiles from
# each divergence group in all 6 environments
plotlist21 <- vector(mode = "list", length = length(ExperimentNames))
names(plotlist21) <- ExperimentNames
plotlist12 <- vector(mode = "list", length = length(ExperimentNames))
names(plotlist12) <- ExperimentNames
for (e in ExperimentNames) {
  gene_idxs21 <- finaldf |> filter(experiment == e &
                                     cer == 2 & par == 1) |> 
    select(gene_name) |> pull()
  gene_idxs12 <- finaldf |> filter(experiment == e &
                                     cer == 1 & par == 2) |> 
    select(gene_name) |> pull()
  plotlist21[[e]] <- plotExpressionProfilePair(.cts1 = collapsed$cer[gene_idxs21,],
                                               .cts2 = collapsed$par[gene_idxs21,],
                                               .info1 = info,
                                               .info2 = info,
                                               .method = "line",
                                               .show_points = FALSE,
                                               .show_confidence_intervals = TRUE,
                                               .normalization = "centered log2",
                                               .plot_titles = "ngenes")
  plotlist12[[e]] <- plotExpressionProfilePair(.cts1 = collapsed$cer[gene_idxs12,],
                                               .cts2 = collapsed$par[gene_idxs12,],
                                               .info1 = info,
                                               .info2 = info,
                                               .method = "line",
                                               .show_points = FALSE,
                                               .show_confidence_intervals = TRUE,
                                               .normalization = "centered log2",
                                               .plot_titles = "ngenes")
}
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_env_patterns21.pdf",
    width = 13, height = 10)
ggarrange(plotlist = plotlist21, nrow = 6, ncol = 1)
dev.off()
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_env_patterns12.pdf",
    width = 13, height = 10)
ggarrange(plotlist = plotlist12, nrow = 6, ncol = 1)
dev.off()

#### Heatmaps of level divergers, Upcer Uppar separated ####
plot_mat <- matrix(nrow = length(ExperimentNames),
                   ncol = length(ExperimentNames))
# up in Scer
for (e_row in ExperimentNames) {
  for (e_col in ExperimentNames) {
    e_gene_idxs <- finaldf |> filter(experiment == e_row &
                                       level == "diverged" &
                                       sign(effect_size_species) == 1) |> 
      select(gene_name) |> pull()
    avg_lfc <- finaldf |> filter(experiment == e_col & 
                                   gene_name %in% e_gene_idxs) |> 
      select(effect_size_species) |> pull() |> mean()
    plot_mat[which(ExperimentNames == e_row),
             which(ExperimentNames == e_col)] <- avg_lfc
  }
}
colnames(plot_mat) <- LongExperimentNames
rownames(plot_mat) <- LongExperimentNames

col_fun = colorRamp2(c(-1, 0, 1), c("blue2", "white", "orange1"))
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/level_divergence_heatmap_up_scer.pdf",
    width = 5, height = 3)
Heatmap(plot_mat, col = col_fun,
        row_order = LongExperimentNames, column_order = LongExperimentNames,
        row_names_side = "left", heatmap_legend_param = list(title = ""),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", plot_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()
# up in Spar
plot_mat <- matrix(nrow = length(ExperimentNames),
                   ncol = length(ExperimentNames))
for (e_row in ExperimentNames) {
  for (e_col in ExperimentNames) {
    e_gene_idxs <- finaldf |> filter(experiment == e_row &
                                       level == "diverged" & 
                                       sign(effect_size_species) == -1) |> 
      select(gene_name) |> pull()
    avg_lfc <- finaldf |> filter(experiment == e_col & 
                                   gene_name %in% e_gene_idxs) |> 
      select(effect_size_species) |> pull() |> mean()
    plot_mat[which(ExperimentNames == e_row),
             which(ExperimentNames == e_col)] <- avg_lfc
  }
}
colnames(plot_mat) <- LongExperimentNames
rownames(plot_mat) <- LongExperimentNames

col_fun = colorRamp2(c(-1, 0, 1), c("blue2", "white", "orange1"))
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/level_divergence_heatmap_up_spar.pdf",
    width = 5, height = 3)
Heatmap(plot_mat, col = col_fun,
        row_order = LongExperimentNames, column_order = LongExperimentNames,
        row_names_side = "left", heatmap_legend_param = list(title = ""),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", plot_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

#### 10% Highest cor vs 10% lowest cor in each environment ####

# we've confirmed that dynamics divergence in one environment doesn't predict
# dynamics divergence in others. But is the opposite true?
# Do genes with strong correlation in one environment have divergent dynamics
# in another?
# The Pho4 regulon in LowPi with divergent dynamics in HAP4 is our best example of this

# Top n most correlated
nGenes <- 500
plotlist <- vector(mode = "list", length = length(ExperimentNames))
names(plotlist) <- ExperimentNames
for (e in ExperimentNames) {
  top_table <- finaldf |> filter(experiment == e & cer != 0 & par != 0) |> 
    mutate(cor_rank = rank(-cor_parents)) |> 
    filter(cor_rank <= nGenes) |> select(cer, par) |> table()
  max_table <- which(top_table == max(top_table), arr.ind = TRUE)
  cer_clust <- colnames(top_table)[as.numeric(max_table[1, 1])]
  par_clust <- colnames(top_table)[as.numeric(max_table[1, 2])]
  top_gene_idxs <- finaldf |> filter(experiment == e & cer != 0 & par != 0) |> 
    mutate(cor_rank = rank(-cor_parents)) |> 
    filter(cor_rank <= nGenes & cer == cer_clust & par == par_clust) |> # just plotting the subset from the most represented cluster
    select(gene_name) |> pull()
  plotlist[[e]] <- annotate_figure(plotEnvironments(.gene_idxs = top_gene_idxs,
                                                    .normalization = "scale"),
                                   top = paste(e, length(top_gene_idxs), "genes"))
}
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/top", nGenes, ".pdf"),
    width = 10, height = 10)
ggarrange(plotlist = plotlist, 
          nrow = length(ExperimentNames), ncol = 1,
          common.legend = TRUE)
dev.off()

# TODO: Bottom nGenes least correlated
plotlist <- vector(mode = "list", length = length(ExperimentNames))
names(plotlist) <- ExperimentNames
for (e in ExperimentNames) {
  bottom_table <- finaldf |> filter(experiment == e & cer != 0 & par != 0) |> 
    mutate(cor_rank = rank(cor_parents)) |> 
    filter(cor_rank <= nGenes) |> select(cer, par) |> table()
  max_table <- which(bottom_table == max(bottom_table), arr.ind = TRUE)
  cer_clust <- colnames(bottom_table)[as.numeric(max_table[1, 1])]
  par_clust <- colnames(bottom_table)[as.numeric(max_table[1, 2])]
  bottom_gene_idxs <- finaldf |> filter(experiment == e & cer != 0 & par != 0) |> 
    mutate(cor_rank = rank(cor_parents)) |> 
    filter(cor_rank <= nGenes & cer == cer_clust & par == par_clust) |> # just plotting the subset from the most represented cluster
    select(gene_name) |> pull()
  plotlist[[e]] <- annotate_figure(plotEnvironments(.gene_idxs = bottom_gene_idxs,
                                                    .normalization = "scale"),
                                   top = paste(e, length(bottom_gene_idxs), "genes"))
}
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/bottom", nGenes, ".pdf"),
    width = 10, height = 10)
ggarrange(plotlist = plotlist, 
          nrow = length(ExperimentNames), ncol = 1,
          common.legend = TRUE)
dev.off()

# TODO: these plots still suffer from the fact that most genes aren't still in the same cluster
# together in a different environment. So divergent dynamics in other environments is potentially
# masked by divergence in different directions. Couple things to try:
# 1) Screen each pair of environments for a paucity or enrichment of dynamics divergers:
#    First environment is "identified in" environment, where you select top nGenes most correlated btwn species genes
#    Second environment is "measured in" environment, where you check what percent of those genes are 
#    dynamics divergers versus the remaining genes in the genome to see  which pairs of environments
#    have the most different percents, paucity or enrichment
# 2) Plot these paired sets of environments. The idea is that if genes are related in two environments,
#    they're probably related enough in the remaining environments and further dissecting can actually
#    obscure patterns
# 3) Instead of having the "measured in" be % dynamics divergers, it could be parental cor

#### Probably Archive: Heatmaps of pct dynamics and level divergers ####
# probably archive b/c
# Pros: intuitive --- "what percent of level or dynamics divergers in one environment are also diverging in the other environment?"
# Cons: diagonal guaranteed to be 1.0 (highest possible value)

# dynamics divergence
getPercentDynamicsDiv <- function(.gene_idxs, .experiment_name) {
  ndiv <- finaldf |> filter(gene_name %in% .gene_idxs &
                              experiment == .experiment_name &
                              dynamics == "diverged") |> 
    nrow()
  return(ndiv/length(.gene_idxs))
}
# plotting
plot_mat <- matrix(nrow = length(ExperimentNames),
                   ncol = length(ExperimentNames))
for (e_row in ExperimentNames) {
  for (e_col in ExperimentNames) {
    cat(e_row, e_col, "\n")
    e_gene_idxs <- finaldf |> filter(experiment == e_row &
                                       dynamics == "diverged") |> 
      select(gene_name) |> pull()
    e_mean <- getPercentDynamicsDiv(.gene_idxs = e_gene_idxs,
                                    .experiment_name = e_col)
    plot_mat[which(ExperimentNames == e_row),
             which(ExperimentNames == e_col)] <- e_mean
  }
}
colnames(plot_mat) <- LongExperimentNames
rownames(plot_mat) <- LongExperimentNames

# plotting
col_fun = colorRamp2(c(0, 1), c("white", "red2"))
Heatmap(plot_mat, col = col_fun,
        row_order = LongExperimentNames, column_order = LongExperimentNames,
        row_names_side = "left", heatmap_legend_param = list(title = ""))

# level divergence
getPercentLevelDiv <- function(.gene_idxs, .experiment_name) {
  ndiv <- finaldf |> filter(gene_name %in% .gene_idxs &
                              experiment == .experiment_name &
                              level == "diverged") |> 
    nrow()
  return(ndiv/length(.gene_idxs))
}
# plotting
plot_mat <- matrix(nrow = length(ExperimentNames),
                   ncol = length(ExperimentNames))
for (e_row in ExperimentNames) {
  for (e_col in ExperimentNames) {
    cat(e_row, e_col, "\n")
    e_gene_idxs <- finaldf |> filter(experiment == e_row &
                                       level == "diverged") |> 
      select(gene_name) |> pull()
    e_mean <- getPercentLevelDiv(.gene_idxs = e_gene_idxs,
                                 .experiment_name = e_col)
    plot_mat[which(ExperimentNames == e_row),
             which(ExperimentNames == e_col)] <- e_mean
  }
}
colnames(plot_mat) <- LongExperimentNames
rownames(plot_mat) <- LongExperimentNames

# plotting
col_fun = colorRamp2(c(0, 1), c("white", "red2"))
Heatmap(plot_mat, col = col_fun,
        row_order = LongExperimentNames, column_order = LongExperimentNames,
        row_names_side = "left", heatmap_legend_param = list(title = ""))

#### Level divergence is consistent across environments ####

# Heat map ordering genes by lfc in one environment per row (cols are genes)
# Probably filter down to genes that are DE in at least one experiment
levmatdf <- finaldf |> filter(experiment %in% ExperimentNames) |> 
  mutate(effect_size_species = if_else((pvalue_species > 1e-5 & abs(effect_size_species) > 1),
                                       true = 0,
                                       false = effect_size_species)) |> 
  select(gene_name, experiment, effect_size_species) |> 
  pivot_wider(id_cols = experiment, names_from = gene_name, values_from = effect_size_species)
levmat <- as.matrix(levmatdf[,-1])
rownames(levmat) <- levmatdf$experiment
levmat[is.na(levmat)] <- 0
at_least_1_idxs <- finaldf |> filter(level == "diverged") |>
  select(gene_name) |> pull() |> unique() |> intersect(y = colnames(levmat))

# plotting
col_fun = colorRamp2(c(-2, 0, 2), c("blue2", "lightyellow", "orange1"))
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/lev_heatmap.pdf",
    width = 7, height = 2.5)
Heatmap(levmat[, at_least_1_idxs],
        column_labels = rep("", length(at_least_1_idxs)),
        col = col_fun)
dev.off()

# example level-diverging gene groups across all environments
# Diauxic Shift increasing cluster, up in Scer
gene_idxs <- finaldf |> filter(experiment == "HAP4" & group == "levupcer1") |> 
  select(gene_name) |> pull()
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/avg_expr_level.pdf",
    width = 7, height = 3)
plotExpressionProfilePair(.cts1 = collapsed$cer[gene_idxs, info$experiment %in% ExperimentNames],
                          .cts2 = collapsed$par[gene_idxs, info$experiment %in% ExperimentNames],
                          .info1 = info[info$experiment %in% ExperimentNames,],
                          .info2 = info[info$experiment %in% ExperimentNames,],
                          .method = "line",
                          .show_points = FALSE,
                          .show_confidence_intervals = TRUE,
                          .normalization = "log2",
                          .plotlims = c(4, 9))
dev.off()

### Upset plot of strongest level divergers
# TODO: Also collecting gene idxs to see if environment specificity affects cis/trans
plotdf <- finaldf |> select(gene_name, experiment, effect_size_species) |> 
  pivot_wider(id_cols = gene_name, values_from = effect_size_species,
              names_from = experiment) |> 
  drop_na() |> # removing genes that were lowly expressed in any environments
  pivot_longer(cols = c("HAP4", "CC", "LowN", "LowPi", "Heat", "Cold"),
               names_to = "experiment", values_to = "effect_size_species") |> 
  filter(abs(effect_size_species) > 1) |>
  mutate(up_cer = sign(effect_size_species) == 1,
         up_par = sign(effect_size_species) == -1)
# counts of how many environments the divergence was detected in
# tagseq
plotdf_tagseq <- plotdf |> filter(!(experiment %in% c("Heat", "Cold"))) |> 
  group_by(gene_name) |> summarise(n_upcer = sum(effect_size_species > 1),
                                   n_uppar = sum(effect_size_species < -1))
# heat/cold
plotdf_heatcold <- plotdf |> filter(experiment %in% c("Heat", "Cold")) |> 
  group_by(gene_name) |> summarise(n_upcer = sum(effect_size_species > 1),
                                   n_uppar = sum(effect_size_species < -1))

# up_cer, tagseq
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/level_upsets_upcertagseq.pdf",
    width = 5, height = 3)
makeUpsetPlot(filter(plotdf, up_cer), 
              .group_name = "experiment",
              .group_members = c("HAP4", "CC", "LowN", "LowPi"),
              .item_names = "gene_name")
dev.off()
# constitutive level-divergers
constitutive_upcer_tagseq_idxs <- plotdf_tagseq |> filter(n_upcer == 4) |>
  select(gene_name) |> 
  pull()
# single environment level-divergers
hap4_upcer_idxs <- plotdf_tagseq |> filter(n_upcer == 1) |>
  select(gene_name) |> 
  pull() |> 
  intersect(y = filter(plotdf, up_cer & experiment == "HAP4") |> 
              select(gene_name) |> pull())
cc_upcer_idxs <- plotdf_tagseq |> filter(n_upcer == 1) |>
  select(gene_name) |> 
  pull() |> 
  intersect(y = filter(plotdf, up_cer & experiment == "CC") |> 
              select(gene_name) |> pull())
lown_upcer_idxs <- plotdf_tagseq |> filter(n_upcer == 1) |>
  select(gene_name) |> 
  pull() |> 
  intersect(y = filter(plotdf, up_cer & experiment == "LowN") |> 
              select(gene_name) |> pull())
lowpi_upcer_idxs <- plotdf_tagseq |> filter(n_upcer == 1) |>
  select(gene_name) |> 
  pull() |> 
  intersect(y = filter(plotdf, up_cer & experiment == "LowPi") |> 
              select(gene_name) |> pull())

# up_par, tagseq
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/level_upsets_uppartagseq.pdf",
    width = 5, height = 3)
makeUpsetPlot(filter(plotdf, up_par), 
              .group_name = "experiment",
              .group_members = c("HAP4", "CC", "LowN", "LowPi"),
              .item_names = "gene_name")
dev.off()
# constitutive level-divergers
constitutive_uppar_tagseq_idxs <- plotdf_tagseq |> filter(n_uppar == 4) |>
  select(gene_name) |> 
  pull()
# single environment level-divergers
hap4_uppar_idxs <- plotdf_tagseq |> filter(n_uppar == 1) |>
  select(gene_name) |> 
  pull() |> 
  intersect(y = filter(plotdf, up_par & experiment == "HAP4") |> 
              select(gene_name) |> pull())
cc_uppar_idxs <- plotdf_tagseq |> filter(n_uppar == 1) |>
  select(gene_name) |> 
  pull() |> 
  intersect(y = filter(plotdf, up_par & experiment == "CC") |> 
              select(gene_name) |> pull())
lown_uppar_idxs <- plotdf_tagseq |> filter(n_uppar == 1) |>
  select(gene_name) |> 
  pull() |> 
  intersect(y = filter(plotdf, up_par & experiment == "LowN") |> 
              select(gene_name) |> pull())
lowpi_uppar_idxs <- plotdf_tagseq |> filter(n_uppar == 1) |>
  select(gene_name) |> 
  pull() |> 
  intersect(y = filter(plotdf, up_par & experiment == "LowPi") |> 
              select(gene_name) |> pull())

# up_cer, heat/cold
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/level_upsets_upcerheatcold.pdf",
    width = 3, height = 3)
makeUpsetPlot(filter(plotdf, up_cer), 
              .group_name = "experiment",
              .group_members = c("Heat", "Cold"),
              .item_names = "gene_name")
dev.off()
constitutive_upcer_heatcold_idxs <- plotdf_heatcold |> 
  filter(n_upcer == 2) |>
  select(gene_name) |> 
  pull()
heat_upcer_idxs <- plotdf_heatcold |> filter(n_upcer == 1) |>
  select(gene_name) |> 
  pull() |> 
  intersect(y = filter(plotdf, up_cer & experiment == "Heat") |> 
              select(gene_name) |> pull())
cold_upcer_idxs <- plotdf_heatcold |> filter(n_upcer == 1) |>
  select(gene_name) |> 
  pull() |> 
  intersect(y = filter(plotdf, up_cer & experiment == "Cold") |> 
              select(gene_name) |> pull())

# up_par, heat/cold
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/level_upsets_upparheatcold.pdf",
    width = 3, height = 3)
makeUpsetPlot(filter(plotdf, up_par), 
              .group_name = "experiment",
              .group_members = c("Heat", "Cold"),
              .item_names = "gene_name")
dev.off()
constitutive_uppar_heatcold_idxs <- plotdf_heatcold |> 
  filter(n_uppar == 2) |>
  select(gene_name) |> 
  pull()
heat_uppar_idxs <- plotdf_heatcold |> filter(n_uppar == 1) |>
  select(gene_name) |> 
  pull() |> 
  intersect(y = filter(plotdf, up_par & experiment == "Heat") |> 
              select(gene_name) |> pull())
cold_uppar_idxs <- plotdf_heatcold |> filter(n_uppar == 1) |>
  select(gene_name) |> 
  pull() |> 
  intersect(y = filter(plotdf, up_par & experiment == "Cold") |> 
              select(gene_name) |> pull())

### Average expression and Gene Ontology of select level-divergers
plotlims <- c(4, 10)

# 1) Constitutive level-divergers upcer tagseq
plotEnvironments(constitutive_upcer_tagseq_idxs)
plotEnvironments(constitutive_upcer_tagseq_idxs, .quartet = TRUE)
getGOSlimDf(.idxs = constitutive_upcer_tagseq_idxs, 
            .group_name = "upcer_tagseq",
            .min_hits = 3) |> View() # extracellular/cell wall proteins

# 2) Consistent level-divergers uppar tagseq
plotEnvironments(constitutive_uppar_tagseq_idxs)
plotEnvironments(constitutive_uppar_tagseq_idxs, .quartet = TRUE)
getGOSlimDf(.idxs = gene_idxs, .group_name = "uppar_tagseq",
            .min_hits = 3) |> View() # mitochondrial translation/ribosomes

# 3) Single environment upcer tagseq
# HAP4
plotEnvironments(hap4_upcer_idxs)
plotEnvironments(hap4_upcer_idxs, .quartet = TRUE)
# CC
plotEnvironments(cc_upcer_idxs)
plotEnvironments(cc_upcer_idxs, .quartet = TRUE)
finaldf |> filter(gene_name %in% cc_upcer_idxs) |> 
  select(dynamics, experiment) |> table() # these guys have a lot of dynamics divergence,
# seems like they're expressed higher in Scer in lab conditions 
# (hence why they're detected as level divergers in CC and dynamics
# divergers in other environments)
# LowN
plotEnvironments(lown_upcer_idxs)
plotEnvironments(lown_upcer_idxs, .quartet = TRUE)
# LowPi
plotEnvironments(lowpi_upcer_idxs)
plotEnvironments(lowpi_upcer_idxs, .quartet = TRUE) # pretty trans!
finaldf |> filter(gene_name %in% lowpi_upcer_idxs) |> 
  select(dynamics, experiment) |> table() # genes aren't any more likely to diverge in dynamics in LowPi

# 4) Single environment uppar tagseq
# HAP4
plotEnvironments(hap4_uppar_idxs) # fairly trans
plotEnvironments(hap4_uppar_idxs, .quartet = TRUE)
finaldf |> filter(gene_name %in% hap4_uppar_idxs) |> 
  select(dynamics, experiment) |> table() 
# a lot of dynamics divergence, but not significantly more in HAP4
# CC
plotEnvironments(cc_uppar_idxs)
plotEnvironments(cc_uppar_idxs, .quartet = TRUE)
finaldf |> filter(gene_name %in% cc_uppar_idxs) |> 
  select(dynamics, experiment) |> table() # not as much dynamics divergence as upcer cc

# LowN
plotEnvironments(lown_uppar_idxs)
plotEnvironments(lown_uppar_idxs, .quartet = TRUE) # these just look like weak level-divergers in every environment
# LowPi
plotEnvironments(lowpi_uppar_idxs)
plotEnvironments(lowpi_uppar_idxs, .quartet = TRUE) # these just look like weak level-divergers in every environment

# 5) Constitutive level-divergers upcer heat/cold
plotEnvironments(constitutive_upcer_heatcold_idxs)
plotEnvironments(constitutive_upcer_heatcold_idxs, .quartet = TRUE)
getGOSlimDf(.idxs = constitutive_upcer_heatcold_idxs, .group_name = "upcer_heatcold",
            .min_hits = 3) |> View() # mitochondrial

# 6) Consistent level-divergers uppar heat/cold
plotEnvironments(constitutive_uppar_heatcold_idxs)
plotEnvironments(constitutive_uppar_heatcold_idxs, .quartet = TRUE)
getGOSlimDf(.idxs = constitutive_uppar_heatcold_idxs, .group_name = "uppar_heatcold",
            .min_hits = 3) |> View() # ribosome

# 7) Upcer/Uppar unique to Heat/Cold
plotEnvironments(heat_upcer_idxs)
plotEnvironments(heat_upcer_idxs, .quartet = TRUE)
plotEnvironments(heat_uppar_idxs)
plotEnvironments(heat_uppar_idxs, .quartet = TRUE)
plotEnvironments(cold_upcer_idxs)
plotEnvironments(cold_upcer_idxs, .quartet = TRUE)
plotEnvironments(cold_uppar_idxs)
plotEnvironments(cold_uppar_idxs, .quartet = TRUE)

# Takeaways:
# 1) level divergence is consistent between environments for the same strains
# 2) level divergence is 50-50 conserved-not conserved between different strains
# 3) environment-specific level divergence has a moderate association with dynamics divergence

#### Level divergence is independent from dynamics ####
# We know from Krieger et al. 2020 that divergence in level does not predict dynamics
# divergence in dynamics and vice-versa, but it's worth illustrating how this looks
# among genes in different dynamics categories. How it's simply a raising or lowering of the
# same dynamic pattern

# For each divergence category (cons1, cons2, dyn12, dyn21) in each experiment, 
# highlight portions of the l2fc density plot showing how the curves stay the same
# and the height difference between species is what is changing

griddf1 <- expand_grid(experiment = c("CC", "HAP4", "LowN", "LowPi"),
                      cer = c(1, 2),
                      par = c(1, 2))
griddf2 <- expand_grid(experiment = c("Cold"),
                       cer = c(1, 2, 3),
                       par = c(1, 2, 3))
griddf3 <- expand_grid(experiment = c("Heat"),
                       cer = c(1, 2, 3, 4),
                       par = c(1, 2, 3, 4))
griddf <- bind_rows(griddf1, griddf2, griddf3)

# plotting
for (i in 1:nrow(griddf)) {
  ex <- griddf$experiment[i]
  cerclust <- griddf$cer[i]
  parclust <- griddf$par[i]
  plotdf <- filter(finaldf, experiment == ex & cer == cerclust & par == parclust) |> 
    select(effect_size_species, gene_name)
  if (nrow(plotdf) < 10) {
    next
  }
  distbins <- quantile(plotdf$effect_size_species, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
  pdens <- ggplot(plotdf, aes(x = effect_size_species)) + geom_density() + theme_classic() +
    xlab("log2 fold change") + 
    geom_vline(xintercept = distbins)
  plotlist <- vector(mode = "list", length = length(distbins) - 1)
  for (i in 1:(length(distbins) - 1)) {
    gidxs <- filter(plotdf, effect_size_species > distbins[i] & 
                      effect_size_species <= distbins[i + 1]) |> 
      select(gene_name) |> pull()
    plotlist[[i]] <- annotate_figure(plotGenes(gidxs, .plotlims = c(0, 1200), .experiment_name = ex,
                               .normalization = "none"), top = paste0(length(gidxs), " genes"))
  }
  pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/level_independence/", 
                    ex, cerclust, parclust, ".pdf"), width = 9, height = 5)
  print(annotate_figure(ggarrange(pdens,
            ggarrange(plotlist = plotlist, ncol = length(plotlist), nrow = 1),
            nrow = 2), top = paste(nrow(plotdf), "genes", ex, "\nScer cluster:", cerclust, 
                                    "Spar cluster:", parclust)))
  dev.off()
}

#### discrete heatmap of level and dynamics divergence across experiments ####

plotdf <- finaldf
plot_mat <- plotdf |> 
  select(gene_name, experiment, group4) |>
  pivot_wider(id_cols = "gene_name",
              names_from = "experiment", 
              values_from = "group4")
colnames_plotmat <- plot_mat$gene_name
rownames_plotmat <- colnames(plot_mat)
plot_mat <- data.frame(plot_mat) |> t()
colnames(plot_mat) <- colnames_plotmat
rownames(plot_mat) <- rownames_plotmat
plot_mat <- plot_mat[rownames(plot_mat) != "gene_name",]
plot_mat <- plot_mat |> as.matrix()
clrs <- structure(levdyn_colordf$type, 
                  names = levdyn_colordf$limits)
majority_NA <- apply(plot_mat, 2, \(x) {return(sum(is.na(x)) > 2)})
sum(majority_NA)
plot_mat <- plot_mat[,!majority_NA]
# sorting cols by HAP4's order
col_order_vec <- order(factor(plot_mat["HAP4",], levels = levdyn_colordf$limits))
# recursively order each subset of previous row
orderGenesByGroup <- function(.mat, .row_idx = 1,
                              .labels = levdyn_colordf$limits) {
  if (ncol(.mat) == 1) {
    return(.mat)
  }
  vec <- .mat[.row_idx,]
  labels_present <- .labels[.labels %in% unique(vec)]
  breaks <- factor(vec, levels = labels_present) |> 
    table(useNA = "ifany") |> as.numeric()
  cumulative_breaks <- breaks
  for(i in 1:length(breaks)) {
    cumulative_breaks[i] <- sum(breaks[1:i])
  }
  vec_order <- factor(vec, levels = .labels) |> order()
  reordered_mat <- .mat[,vec_order]
  if (.row_idx == nrow(.mat)) {
    return(reordered_mat)
  }
  out_vec <- map2(.x = c(0, cumulative_breaks[-length(cumulative_breaks)]) + 1, 
                  .y = cumulative_breaks, 
                  .f = \(col_start, col_end) {
                    orderGenesByGroup(.mat = reordered_mat[, col_start:col_end, drop = FALSE],
                                      .row_idx = .row_idx + 1,
                                      .labels = .labels)
                  }) |> purrr::reduce(.f = cbind)
  return(out_vec)
}
# # tests for orderGenesByGroup
# orderGenesByGroup(.mat = rbind(c("hai", "der", "der", "hai"),
#                                c("der", "der", "der", "hai")),
#                   .labels = c("hai", "der"))
# # subset of plotmat
# orderGenesByGroup(.mat = plot_mat[, 1:100]) |> dim()

# ordering full plotmat
ordered_plot_mat <- orderGenesByGroup(.mat = plot_mat[ExperimentNames,])
rownames(ordered_plot_mat) <- LongExperimentNames
# plotting
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/discrete_heatmap.pdf",
    width = 9, height = 1.5)
Heatmap(ordered_plot_mat, 
        col = clrs, show_column_names = FALSE,
        row_order = rownames(ordered_plot_mat), 
        column_order = colnames(ordered_plot_mat),
        row_names_side = "left", heatmap_legend_param = list(title = ""))
dev.off()

#### Example expression profiles for figure ####
# Low Phosphorus 1-2, dynamics divergers
gene_idxs <- finaldf |> filter(experiment == "LowPi" & group == "dyn12") |> 
  select(gene_name) |> pull()
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/LowPi12.pdf",
    width = 12, height = 2)
plotExpressionProfilePair(.cts1 = collapsed$cer[gene_idxs,],
                          .cts2 = collapsed$par[gene_idxs,],
                          .info1 = info,
                          .info2 = info,,
                          .method = "line",
                          .show_points = FALSE,
                          .show_confidence_intervals = TRUE,
                          .normalization = "log2")
dev.off()

# upCer LowN, level divergers
gene_idxs <- finaldf |> filter(experiment == "LowN" & level == "diverged" &
                                 sign(effect_size_species) == 1) |> 
  select(gene_name) |> pull()
pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/LowN_upcer.pdf",
    width = 12, height = 2)
plotExpressionProfilePair(.cts1 = collapsed$cer[gene_idxs,],
                          .cts2 = collapsed$par[gene_idxs,],
                          .info1 = info,
                          .info2 = info,,
                          .method = "line",
                          .show_points = FALSE,
                          .show_confidence_intervals = TRUE,
                          .normalization = "log2")
dev.off()

#### Does the most conserved/correlated/variable subset of each environment have more divergent expression in other environments? ####
# TODO: see if we can sidestep the yeastract and just take the most strongly
# co-regulated subset of each environmental cluster and look at their dynamics
# in other environments compared to the rest of that cluster
plotCoRegGroup <- function(.environment, .clust, .pct = 0.1) {
  
}

############################## Archive ######################################
# #### Probably Archive from here on down: Level and dynamics divergers only retain level divergence across environments ####
# plotCode <- function(.code, .experiment = "LowN", .normalization = "log2",
#                      .plotlims = NULL) {
#   plotdf <- clusterdf |> filter(experiment == .experiment)
#   plotdf$code <- apply(plotdf, 1, \(x) {
#     paste(as.numeric(x["cer"]), 
#           as.numeric(x["par"]), 
#           as.numeric(x["hyc"]), 
#           as.numeric(x["hyp"]), sep = " ")
#   })
#   gene_idxs <- plotdf |> filter(code == .code) |> select(gene_ID) |> pull()
#   plotExpressionProfileQuartet(.cts1 = counts_all2$cer[info$experiment == .experiment, gene_idxs],
#                                .cts2 = counts_all2$par[info$experiment == .experiment, gene_idxs],
#                                .cts3 = counts_all2_allele$cer[info$experiment == .experiment, gene_idxs],
#                                .cts4 = counts_all2_allele$par[info$experiment == .experiment, gene_idxs],
#                                .info1 = info[info$experiment == .experiment,],
#                                .info2 = info[info$experiment == .experiment,],
#                                .info3 = info[info$experiment == .experiment,],
#                                .info4 = info[info$experiment == .experiment,],
#                                .name1 = "S. cer",
#                                .name2 = "S. par",
#                                .name3 = "F1, cer allele",
#                                .name4 = "F1, par allele",
#                                .color1 = "orange1",
#                                .color2 = "blue2",
#                                .color3 = "orange4",
#                                .color4 = "blue4",
#                                .normalization = .normalization,
#                                .method = "line",
#                                .show_points = FALSE,
#                                .show_confidence_intervals = TRUE,
#                                .plotlims = .plotlims)
# }
# #### Look at most interesting divergence examples in each environment ####
# 
# # cer dominance:
# # Diauxic Shift: par up cer/hyc/hyp down (2122 in 2 clust)
# plotCode("2 1 2 2", .experiment = "HAP4")
# # Low Pi: par down cer/hyc/hyp up (1211 in 3 clust)
# plotCode("1 2 1 1", .experiment = "LowPi")
# 
# # LowN:
# # Low N: cer/hyc/hyp up TP2 and 3, par only up TP3 (only visible in 6 clusters)
# #        par/hyc/hyp up TP2 and 3, cer down TP3 (3111 in 4 clust)
# plotCode("3 1 1 1", .experiment = "LowN")
# #        par up TP3, cer/hyc/hyp down TP2 and 3 (2422 in 4 clust)
# plotCode("2 4 2 2", .experiment = "LowN")
# #        par down TP2, sometimes also hyp (1411 and 1414 in 4 clust)
# plotCode("1 4 1 1", .experiment = "LowN", .plotlims = c(5, 8))
# plotCode("1 4 1 4", .experiment = "LowN", .plotlims = c(5, 8))
# # ^ Look at LowN ones in TF deletions especially
# 
# # par dominance:
# # Cold: par/hyc/hyp up, cer down
# plotCode("2 1 1 1", .experiment = "Cold", .plotlims = c(6, 8))
# 
# # wild hybrid:
# # Heat: hybrid more consistent patterns overall than either parent
# plotCode("2 1 1 1", .experiment = "Heat")
# plotCode("1 2 1 1", .experiment = "Heat")
# plotCode("1 2 1 2", .experiment = "Heat")
# plotCode("1 3 3 3", .experiment = "Heat", .plotlims = c(4,9)) # less-wild hybrid
# plotCode("3 1 1 1", .experiment = "Heat", .plotlims = c(5.5,8.5)) # very similar to 1 3 3 3, but now Spar increases
# 
# # Leaving out CC/Urea until I can figure out if it's worth using for clusters
# # or if it's better as a cross reference for how CC dynamics might
# # be related to the modules we find
# 
# # First off, are these all the same genes divergent in every experiment?
# # or are different subsets divergent in different experiments?
# clusterdf$code <- apply(clusterdf, 1, \(x) {
#   paste(as.numeric(x["cer"]), 
#         as.numeric(x["par"]), 
#         as.numeric(x["hyc"]), 
#         as.numeric(x["hyp"]), sep = " ")
# })
# HAP4_2122 <- clusterdf |> filter(experiment == "HAP4" &
#                                    code == "2 1 2 2") |> 
#   select(gene_ID) |> pull()
# LowPi_1211 <- clusterdf |> filter(experiment == "LowPi" &
#                                    code == "1 2 1 1") |> 
#   select(gene_ID) |> pull()
# LowN_3111 <- clusterdf |> filter(experiment == "LowN" &
#                                     code == "3 1 1 1") |> 
#   select(gene_ID) |> pull()
# LowN_2422 <- clusterdf |> filter(experiment == "LowN" &
#                                    code == "2 4 2 2") |> 
#   select(gene_ID) |> pull()
# LowN_1411 <- clusterdf |> filter(experiment == "LowN" &
#                                    code == "1 4 1 1") |> 
#   select(gene_ID) |> pull()
# LowN_1414 <- clusterdf |> filter(experiment == "LowN" &
#                                    code == "1 4 1 4") |> 
#   select(gene_ID) |> pull()
# Cold_2111 <- clusterdf |> filter(experiment == "Cold" &
#                                    code == "2 1 1 1") |> 
#   select(gene_ID) |> pull()
# Heat_2111 <- clusterdf |> filter(experiment == "Heat" &
#                                    code == "2 1 1 1") |> 
#   select(gene_ID) |> pull()
# table(c(HAP4_2122, LowPi_1211, LowN_3111,
#         LowN_2422, LowN_1411, LowN_1414,
#         Cold_2111, Heat_2111)) |> table() # most only appear once
# 
# table(c(HAP4_2122, Heat_2111)) |> table() # these two have a decent number in common
# table(c(HAP4_2122, LowN_3111)) |> table() # so do these
# 
# # Level ones
# # LowN
# LowN_1NA1NA <- clusterdf |> filter(experiment == "LowN" &
#                                    code == "1 NA 1 NA") |> 
#   select(gene_ID) |> pull()
# LowN_1NA11 <- clusterdf |> filter(experiment == "LowN" &
#                                      code == "1 NA 1 1") |> 
#   select(gene_ID) |> pull()
# # Heat
# Heat_1NANANA <- clusterdf |> filter(experiment == "Heat" &
#                                      code == "1 NA NA NA") |> 
#   select(gene_ID) |> pull()
# # Cold
# Cold_NA1NA1 <- clusterdf |> filter(experiment == "Cold" &
#                                       code == "NA 1 NA 1") |> 
#   select(gene_ID) |> pull()
# Cold_1NANANA <- clusterdf |> filter(experiment == "Cold" &
#                                      code == "1 NA NA NA") |> 
#   select(gene_ID) |> pull()
# # SatGrowth
# HAP4_1NA11 <- clusterdf |> filter(experiment == "HAP4" &
#                                       code == "1 NA 1 1") |> 
#   select(gene_ID) |> pull()
# HAP4_1NA1NA <- clusterdf |> filter(experiment == "HAP4" &
#                                     code == "1 NA 1 NA") |> 
#   select(gene_ID) |> pull()
# HAP4_NA1NA1 <- clusterdf |> filter(experiment == "HAP4" &
#                                     code == "NA 1 NA 1") |> 
#   select(gene_ID) |> pull()
# 
# # not as much overlap as I was expecting
# table(c(LowN_1NA1NA, HAP4_1NA1NA)) |> table()
# 
# # How are these divergent gene groups behaving in all environments?
# annotate_figure(plotExpressionProfileQuartet(.cts1 = counts_all2$cer[, HAP4_1NA1NA],
#                              .cts2 = counts_all2$par[, HAP4_1NA1NA],
#                              .cts3 = counts_all2_allele$cer[, HAP4_1NA1NA],
#                              .cts4 = counts_all2_allele$par[, HAP4_1NA1NA],
#                              .info1 = info,
#                              .info2 = info,
#                              .info3 = info,
#                              .info4 = info,
#                              .name1 = "S. cer",
#                              .name2 = "S. par",
#                              .name3 = "F1, cer allele",
#                              .name4 = "F1, par allele",
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .normalization = "log2",
#                              .method = "line",
#                              .show_points = FALSE,
#                              .show_confidence_intervals = TRUE),
#                 top = "cis level divergence")
# 
# # SatGrowth genes are interestingly lower level in Cell Cycle
# # LowPi genes only have divergence in Spar
# # LowN 3111 interestingly have way stronger expression in SatGrowth in par, also lower level in par CC
# # LowN 2422 seems to be just slightly higher level in par by the end of SatGrowth, LowN, and LowPi
# # LowN 1411 higher in par CC, 1414 is just too noisy
# # Cold 2111 only visible in Cold
# # Heat 2111 also visible in SatGrowth
# 

# #### Cluster conservation data wrangling ####
# 
# # given at least 2 vectors, checks if all of them are equal (order matters, NAs allowed)
# vecsEqual <- function(.x, .y, ...) {
#   checkVecs <- function(.a, .b) {
#     if (identical(.a, .b)) {
#       return(.a)
#     }
#     if (!identical(.a, .b)) {
#       return(FALSE)
#     }
#   }
#   consensus <- reduce(list(.x, .y, ...), .f = checkVecs)
#   return(!isFALSE(consensus))
# }
# # tests for vecsEqual
# vecsEqual(c(1,1,1), c(1,1,1)) # should be TRUE
# vecsEqual(c(1,1,1), c(1,1,1), c(1,1,2)) # should be FALSE
# vecsEqual(c(1,2,3), c(1,2,3)) # should be TRUE
# vecsEqual(c(1,2,3), c(1,2,3), c(1,2,3)) # should be TRUE
# vecsEqual(c(1,2,3), c(1,2,4), c(1,2,3)) # should be FALSE
# vecsEqual(c(3,2,1), c(1,2,3), c(1,2,3)) # should be FALSE
# vecsEqual(c(3,2,1), c(NA,2,3), c(1,2,3)) # should be FALSE
# vecsEqual(c(1,2,NA), c(1,2,NA), c(1,2,NA)) # should be TRUE
# testvec_cer <- clusterdf |> filter(gene_ID == "YAL013W") |> 
#   select(cer) |> pull()
# testvec_par <- clusterdf |> filter(gene_ID == "YAL013W") |> 
#   select(par) |> pull()
# testvec_hyc <- clusterdf |> filter(gene_ID == "YAL013W") |> 
#   select(hyc) |> pull()
# testvec_hyp <- clusterdf |> filter(gene_ID == "YAL013W") |> 
#   select(hyp) |> pull()
# vecsEqual(testvec_cer, testvec_par, testvec_hyc, testvec_hyp) # should be false (and not throw a warning about vectors being different lengths)
# 
# getConsGroup <- function(.vec_list) {
#   # removing NAs, if any vector has entirely NA, returns "low"
#   vecdf <- bind_cols(.vec_list) |> drop_na()
#   if (nrow(vecdf) == 0) {
#     return("low")
#   }
#   .vec_list <- select(vecdf, colnames(vecdf)) |> lapply(unlist)
#   if (vecsEqual(.vec_list[[1]], .vec_list[[2]], .vec_list[[3]], .vec_list[[4]])) {
#     return(reduce(names(.vec_list), .f = paste, sep = "_"))
#   }
#   for (idx in c(1:length(.vec_list))) {
#     s <- .vec_list[idx]
#     vec_list3 <- .vec_list[setdiff(c(1:length(.vec_list)), idx)]
#     if (vecsEqual(vec_list3[[1]], vec_list3[[2]], vec_list3[[3]])) {
#       return(reduce(names(vec_list3), .f = paste, sep = "_"))
#     }
#   }
#   for (idx1 in c(1:length(.vec_list))) {
#     for (idx2 in setdiff(c(1:length(.vec_list)), idx1)) {
#       s <- .vec_list[c(idx1, idx2)]
#       vec_list2 <- .vec_list[setdiff(c(1:length(.vec_list)), c(idx1, idx2))]
#       if (vecsEqual(vec_list2[[1]], vec_list2[[2]])) {
#         if (vecsEqual(s[[1]], s[[2]])) { # if one pair is equal, you have to check if the other pair is equal too
#           group1 <- reduce(names(vec_list2), .f = paste, sep = "_")
#           group2 <- reduce(names(s), .f = paste, sep = "_")
#           return(reduce(c(group1, group2), .f = paste, sep = "_and_"))
#         }
#         if (!vecsEqual(s[[1]], s[[2]])) {
#           return(reduce(names(vec_list2), .f = paste, sep = "_"))
#         }
#       }
#     }
#   }
#   return("none")
# }
# # tests for getConsGroup
# gene_idx <- c("YMR086W") # cons cer-par and hyc-hyp
# gene_idx <- c("YAL003W") # fully conserved
# gene_idx <- c("YOR126C") # fully diverged
# gene_idx <- c("YOL125W") # cer-hyc conserved
# gene_idx <- c("YAL033W") # 1 1 1 3 in LowN, should have 3 species conserved but it currently evaluates to "none"
# testvec_cer <- clusterdf |> filter(gene_ID == gene_idx) |> 
#   select(cer) |> pull()
# testvec_par <- clusterdf |> filter(gene_ID == gene_idx) |> 
#   select(par) |> pull()
# testvec_hyc <- clusterdf |> filter(gene_ID == gene_idx) |> 
#   select(hyc) |> pull()
# testvec_hyp <- clusterdf |> filter(gene_ID == gene_idx) |> 
#   select(hyp) |> pull()
# cbind(testvec_cer, testvec_par, testvec_hyc, testvec_hyp) # check which species should be conserved
# getConsGroup(list("cer" = testvec_cer, "par" = testvec_par,
#                   "hyc" = testvec_hyc, "hyp" = testvec_hyp))
# getConsGroup(list("cer" = 1, "par" = 1,
#                   "hyc" = 1, "hyp" = 3))
# getConsGroup(list("cer" = NA, "par" = 1,
#                   "hyc" = 1, "hyp" = 3))
# 
# assignStatusDiscrete <- function(g) {
#   if (g == "cer_par_hyc_hyp") {
#     return("fully conserved")
#   }
#   if (g == "none") {
#     return("fully diverged")
#   }
#   if (g == "hyc_hyp") {
#     return("trans-diverging")
#   }
#   if (g == "cer_hyc_hyp") {
#     return("cer dominant")
#   }
#   if (g == "par_hyc_hyp") {
#     return("par dominant")
#   }
#   if (g %in% c("hyc_hyp_and_cer_par", "cer_par")) {
#     return("compensatory, no allelic effect")
#   }
#   if (g %in% c("cer_par_hyp", "cer_par_hyc")) {
#     return("compensatory, allelic effect")
#   }
#   if (g %in% c("par_hyp", "cer_hyc", "par_hyp_and_cer_hyc")) {
#     return("cis-diverging")
#   }
#   if (g %in% c("cer_hyp", "par_hyc", "par_hyc_and_cer_hyp")) {
#     return("cross-allele")
#   }
#   if (g == "low") {
#     return("lowly expressed in at least 1/4 groups")
#   }
#   if (g == "undetermined") {
#     return("undetermined")
#   }
#   return("you havent thought of everything you dingus")
# }
# 
# # combine these mutually exclusive categories 
# # into fewer categories that make more sense
# # then make barplots to display them
# 
# # applying per experiment
# clusterdf$cons_group <- apply(clusterdf, MARGIN = 1, FUN = \(x) {
#   grp <- getConsGroup(list("cer" = as.numeric(x["cer"]), "par" = as.numeric(x["par"]), 
#                            "hyc" = as.numeric(x["hyc"]), "hyp" = as.numeric(x["hyp"])))
# })
# clusterdf$status <- map(clusterdf$cons_group, assignStatusDiscrete) |> unlist()
# 
# # How many would we expect if we randomly gave genes labels (at same frequency as existing labels)
# random_clusterdf <- map(unique(clusterdf$experiment), \(e) {
#   gene_names <- clusterdf |> drop_na() |>
#     filter(experiment == e) |> 
#     select(gene_ID) |> pull()
#   e_mat <- clusterdf |> drop_na() |>
#     filter(experiment == e) |> 
#     select(cer, par, hyc, hyp) |> as.matrix()
#   n_mat <- nrow(e_mat)*ncol(e_mat)
#   e_mat_shuffled <- e_mat[sample(c(1:n_mat),
#                                  size = n_mat,
#                                  replace = FALSE)] |> 
#     matrix(ncol = 4)
#   colnames(e_mat_shuffled) <- c("cer", "par", "hyc", "hyp")
#   e_mat_shuffled <- as_tibble(e_mat_shuffled)
#   return(bind_cols(e_mat_shuffled, tibble(gene_ID = gene_names, experiment = e)))
# }) |> reduce(bind_rows)
# # adding back in NA genes (scrambling NAs propagates the number of low genes to be like 75% of them)
# na_idxs <- apply(clusterdf, 1, \(x) {any(is.na(c(x["cer"], x["par"], x["hyc"], x["hyp"])))})
# random_clusterdf <- bind_rows(random_clusterdf, clusterdf[na_idxs,])
# 
# # adding divergence group info
# random_clusterdf$cons_group <- apply(random_clusterdf, MARGIN = 1, FUN = \(x) {
#   grp <- getConsGroup(list("cer" = as.numeric(x["cer"]), "par" = as.numeric(x["par"]), 
#                            "hyc" = as.numeric(x["hyc"]), "hyp" = as.numeric(x["hyp"])))
# })
# random_clusterdf$status <- map(random_clusterdf$cons_group, assignStatusDiscrete) |> unlist()
# 
# #### stacked barplots for cluster conservation per experiment ####
# # Two questions: 
# # 1) are real clusters different than randomly assigning 
# #    genes labels in each experiment?
# # 2) are any environments enriched for certain cluster
# #    conservation patterns?
# 
# catagory_order <- c("fully conserved", 
#                     "trans-diverging", 
#                     "cis-diverging", 
#                     "cross-allele",
#                     "compensatory, allelic effect", 
#                     "compensatory, no allelic effect",
#                     "par dominant", 
#                     "cer dominant",
#                     "fully diverged",
#                     "undetermined",
#                     "lowly expressed in at least 1/4 groups")
# # real clusters
# plotdf <- data.frame(table(clusterdf$status, clusterdf$experiment))
# colnames(plotdf) <- c("status", "experiment", "count")
# plotdf$experiment <- factor(plotdf$experiment, 
#                             levels = c("CC", "HAP4", "LowN", "LowPi", "Heat", "Cold"),
#                             labels = c("Cell\nCycle", "Saturated\nGrowth", "Low\nNitrogen",
#                                        "Low\nPhosphorus", "Heat\nStress", "Cold\nStress"))
# plotdf$status <- factor(plotdf$status, 
#                         levels = catagory_order)
# p_cons <- ggplot(plotdf, aes(x = experiment, y = count)) +
#   geom_bar(stat = "identity", position = "stack", aes(fill = status)) +
#   scale_fill_brewer(palette = "Spectral", name = "status",
#                     direction = -1) + ggtitle("Actual Clusters")
# 
# # random
# plotdf <- data.frame(table(random_clusterdf$status, random_clusterdf$experiment))
# colnames(plotdf) <- c("status", "experiment", "count")
# plotdf$experiment <- factor(plotdf$experiment, 
#                             levels = c("CC", "HAP4", "LowN", "LowPi", "Heat", "Cold"),
#                             labels = c("Cell\nCycle", "Saturated\nGrowth", "Low\nNitrogen",
#                                        "Low\nPhosphorus", "Heat\nStress", "Cold\nStress"))
# plotdf$status <- factor(plotdf$status, 
#                         levels = catagory_order)
# p_rand <- ggplot(plotdf, aes(x = experiment, y = count)) +
#   geom_bar(stat = "identity", position = "stack", aes(fill = status)) +
#   scale_fill_brewer(palette = "Spectral", name = "status",
#                     direction = -1) + ggtitle("Random Clusters")
# 
# # plotting
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/conservation_group_barplots.pdf",
#     width = 10, height = 4)
# ggarrange(p_cons, p_rand, ncol = 2, nrow = 1, common.legend = TRUE)
# dev.off()
# 
# #### Combining conservation groups between experiments ####
# 
# # Let's look at random genes to see what we would consider to be 
# # each of the major conservation categories:
# # 1) conserved
# # 2) trans-diverging
# # 3) cis-diverging
# # 4) cross-allele (negative control of sorts for cis-diverging)
# # 5) compensatory
# # 6) really not conserved at all
# 
# # YCL040W is poster child for should be conserved
# gene_idx <- sample(unique(clusterdf$gene_ID), 1)
# clusterdf |> filter(gene_ID == gene_idx)
# test_vec <- clusterdf |> filter(gene_ID == gene_idx) |> 
#   select(cons_group) |> pull()
# plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,gene_idx, drop = FALSE],
#                              .cts2 = counts_all2$par[,gene_idx, drop = FALSE],
#                              .cts3 = counts_all2_allele$cer[,gene_idx, drop = FALSE],
#                              .cts4 = counts_all2_allele$par[,gene_idx, drop = FALSE],
#                              .info1 = info,
#                              .info2 = info,
#                              .info3 = info,
#                              .info4 = info,
#                              .name1 = "cer",
#                              .name2 = "par",
#                              .name3 = "hyc",
#                              .name4 = "hyp",
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .normalization = "log2",
#                              .show_points = FALSE,
#                              .show_confidence_intervals = FALSE,
#                              .method = "loess")
# 
# assignStatusScore <- function(.g_vec, .verbose = FALSE) {
#   # splitting double pairs into separate elements
#   .g_vec <- map(.g_vec, \(g) {str_split(g, pattern = "_and_")}) |> unlist()
#   # quantifying every pair's occurrence
#   pairdf <- tibble(species1 = c("cer", "cer", "cer", "par", "par", "hyc"), 
#                    species2 = c("par", "hyc", "hyp", "hyc", "hyp", "hyp"))
#   pairdf$n <- map2(pairdf$species1, pairdf$species2, \(s1, s2) {
#     return(sum(grepl(s1, .g_vec) & grepl(s2, .g_vec)))
#   }) |> unlist()
#   if (.verbose) {
#     print("pair frequencies:")
#     print(pairdf)
#   }
#   # return whichever is the most common pair. If there are multiple,
#   # returns all members of the pairs
#   max_n <- max(pairdf$n)
#   output <- filter(pairdf, n > 2) |> 
#     select(species1, species2) |> 
#     unlist() |> unique()
#   if (length(output) == 0) {
#     output <- "undetermined"
#   }
#   output <- reduce(output, paste, sep = "_") |> assignStatusDiscrete()
#   return(output)
# }
# test <- assignStatusScore(test_vec, .verbose = TRUE)
# test # TODO: for now I just thresholded to take pairs that show up at least 3 times,
# # but I could consider doing something more mathematically pure
# 
# getConsensusStatus <- function(.stat_vec) {
#   outtab <- sort(table(.stat_vec), decreasing = TRUE)[1]
#   return(names(outtab))
# }
# clusterdf |> filter(gene_ID == "YCL040W") |> 
#   select(status) |> 
#   pull() |> 
#   getConsensusStatus()
# # applying to all genes and see breakdown overall
# df1 <- group_by(clusterdf, gene_ID) |> 
#   summarise(pairs = assignStatusScore(cons_group),
#             consensus = getConsensusStatus(status))
# bind_rows(table(df1$pairs), table(df1$consensus)) |> t()
# # repeat for random
# random_df1 <- group_by(random_clusterdf, gene_ID) |> 
#   summarise(pairs = assignStatusScore(cons_group),
#             consensus = getConsensusStatus(status))
# bind_rows(table(random_df1$pairs), table(random_df1$consensus)) |> t()
# 
# # note that assigning cis/trans across all experiments first
# # results in a bias of trans
# df0 <- group_by(clusterdf, gene_ID) |> 
#   summarise(cons_group = getConsGroup(list(cer = cer,
#                                            par = par,
#                                            hyc = hyc,
#                                            hyp = hyp)))
# table(df0$cons_group)
# df0$status <- map(df0$cons_group, assignStatusDiscrete) |> unlist()
# table(df0$status)
# # it's because hybrids are most likely to have all clusters in common:
# clusterdf |> filter(gene_ID == gene_idx)
# 
# ### Barplot of cross-experiment groups
# # real clusters
# plotdf <- tibble(count = as.numeric(table(df1$consensus)),
#                  consensus = names(table(df1$consensus)))
# plotdf$type = "all experiments"
# plotdf$consensus <- factor(plotdf$consensus, levels = catagory_order)
# p_all <- ggplot(plotdf, aes(x = type, y = count)) +
#   geom_bar(stat = "identity", position = "stack", aes(fill = consensus)) +
#   scale_fill_brewer(palette = "Spectral", name = "consensus",
#                     direction = -1) + ggtitle("All Experiments") +
#   xlab("")
# # random clusters
# plotdf <- tibble(count = as.numeric(table(random_df1$consensus)),
#                  consensus = names(table(random_df1$consensus)))
# plotdf$type = "all experiments"
# plotdf$consensus <- factor(plotdf$consensus, levels = catagory_order)
# p_rand_all <- ggplot(plotdf, aes(x = type, y = count)) +
#   geom_bar(stat = "identity", position = "stack", aes(fill = consensus)) +
#   scale_fill_brewer(palette = "Spectral", name = "consensus",
#                     direction = -1) + ggtitle("All Experiments, Random") +
#   xlab("")
# pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/barplots.pdf",
#     width = 20, height = 5)
# ggarrange(p_cons, p_all, p_rand, p_rand_all, 
#           nrow = 1, ncol = 4, common.legend = TRUE)
# dev.off()
# 
# 
# #### Which experiments reveal plasticity divergence the most? ####
# # TODO: which experiments are orthologs the most variable in?
# #       which ones are they the most conserved in? Controlling for how common 
# #       that cluster is in the general population of genes?
# 
# #### Conservation scores ####
# 
# #### Supplemental Figure: all expression profiles ####
# 
# # merge 10
# plotdf <- moduledf10 |> filter(block_size >= 10)
# plotdf$modtype <- if_else(plotdf$coexpressed == "conserved co-expressed", true = "conserved", false = "non-conserved")
# plotlist_ccm <- vector(mode = "list", length = sum(plotdf$is_CCM))
# plotlist_div <- vector(mode = "list", length = sum(!plotdf$is_CCM))
# plotlist <- list("conserved" = plotlist_ccm,
#                  "non-conserved" = plotlist_div)
# 
# for (idx in 1:nrow(plotdf)) {
#   m <- plotdf$module_name[idx]
#   cat("working on module", m, "\n")
#   gene_idxs <- module_genedf10 |> filter(module_name == m) |> select(gene_name) |> pull()
#   modtype <- if_else(plotdf$is_CCM[idx], true = "conserved", false = "non-conserved")
#   p <- plotExpressionProfilePair(counts_all2$cer[,gene_idxs],
#                                  counts_all2$par[,gene_idxs],
#                                  info,
#                                  info,
#                                  .method = "line", .show_points = FALSE,
#                                  .normalization = "log2")
#   p <- annotate_figure(p, top = paste0(m, 
#                                        " cor_CC=", round(plotdf$avgCor_CC[idx], 3),
#                                        " cor_HAP4=", round(plotdf$avgCor_HAP4[idx], 3),
#                                        " cor_LowN=", round(plotdf$avgCor_LowN[idx], 3),
#                                        " cor_LowPi=", round(plotdf$avgCor_LowPi[idx], 3),
#                                        " cor_Heat=", round(plotdf$avgCor_Heat[idx], 3),
#                                        " cor_Cold=", round(plotdf$avgCor_Cold[idx], 3)))
#   # arranging plots in order of decreasing AvgCor
#   p_idx <- which(sort(plotdf$avgCor[plotdf$modtype == modtype], decreasing = TRUE) == plotdf$avgCor[idx])
#   plotlist[[modtype]][[p_idx]] <- p
# }
# # random
# plotdf <- random_moduledf10 |> filter(block_size >= 30)
# plotlist_Random <- vector(mode = "list", length = nrow(plotdf))
# for (idx in 1:nrow(plotdf)) {
#   m <- plotdf$module_name[idx]
#   gene_idxs <- random_module_genedf10 |> filter(module_name == m) |> select(gene_name) |> pull()
#   p <- plotExpressionProfilePair(counts_all2$cer[,gene_idxs],
#                                  counts_all2$par[,gene_idxs],
#                                  info,
#                                  info,
#                                  .method = "line", .show_points = FALSE,
#                                  .normalization = "log2")
#   p <- annotate_figure(p, top = paste0(m, 
#                                        " cor_CC=", round(plotdf$avgCor_CC[idx], 3),
#                                        " cor_HAP4=", round(plotdf$avgCor_HAP4[idx], 3),
#                                        " cor_LowN=", round(plotdf$avgCor_LowN[idx], 3),
#                                        " cor_LowPi=", round(plotdf$avgCor_LowPi[idx], 3),
#                                        " cor_Heat=", round(plotdf$avgCor_Heat[idx], 3),
#                                        " cor_Cold=", round(plotdf$avgCor_Cold[idx], 3)))
#   # arranging plots in order of decreasing AvgCor
#   p_idx <- which(sort(plotdf$avgCor, decreasing = TRUE) == plotdf$avgCor[idx])
#   plotlist_Random[[p_idx]] <- p
# }
# # making pdfs
# # CCMs
# p_CCM <- ggarrange(plotlist = plotlist[["conserved"]], 
#                    nrow = length(plotlist[["conserved"]]), 
#                    ncol = 1, 
#                    common.legend = TRUE, 
#                    legend = "right")
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_expr_profiles_CCM10.pdf",
#     width = 12, height = length(plotlist_ccm)*2)
# p_CCM
# dev.off()
# # non-conserved modules
# p_Divergent <- ggarrange(plotlist = plotlist[["non-conserved"]], 
#                    nrow = length(plotlist[["non-conserved"]]), 
#                    ncol = 1, 
#                    common.legend = TRUE, 
#                    legend = "right")
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_expr_profiles_Divergent10.pdf",
#     width = 12, height = length(plotlist[["non-conserved"]])*2)
# p_Divergent
# dev.off()
# 
# # Random
# p_Random <- ggarrange(plotlist = plotlist_Random, 
#                       nrow = length(plotlist_Random), 
#                       ncol = 1, 
#                       common.legend = TRUE, 
#                       legend = "right")
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_expr_profiles_Random10.pdf",
#     width = 12, height = length(plotlist_Random)*4)
# p_Random
# dev.off()
# 
# rm(plotdf)
# 
# #### Distribution of module expression conservation ####
# library(ggExtra, ggpubr)
# plotModuleLineAndDensityByBlockSize <- function(.mdf, .rmdf, .show_labels = FALSE, .just_density = FALSE) {
#   plotdf <- .mdf
#   plotdf$modtype <- if_else(plotdf$is_CCM, true = "conserved", false = "non-conserved")
#   random_plotdf <- .rmdf |> mutate(modtype = "random",
#                                    CCM_color = "grey20")
#   plotdf <- bind_rows(plotdf, random_plotdf) |> 
#     mutate(block_rank = rank(block_size, ties.method = "first")) |>
#     select(c("block_rank", setdiff(colnames(plotdf), "avgCor"))) |> 
#     pivot_longer(cols = c("avgCor_CC",
#                           "avgCor_LowN",
#                           "avgCor_LowPi",
#                           "avgCor_HAP4",
#                           "avgCor_Heat",
#                           "avgCor_Cold"),
#                  names_to = "experiment",
#                  values_to = "avgCor",
#                  names_prefix = "avgCor_")
#   plotdf$CCM_color <- gsub("none", "grey80", plotdf$CCM_color)
#   plotdf <- plotdf |> filter(block_size >= 30)
#   # only density plot
#   if (.just_density) {
#     p_density <- ggplot(plotdf, aes(x = avgCor)) +
#       geom_density(aes(fill = modtype), alpha = 0.7) +
#       theme_classic() +
#       theme(legend.title = element_blank()) +
#       xlab("Average expression correlation\n between species") +
#       scale_fill_discrete(limits = c("conserved",
#                                      "non-conserved",
#                                      "random"),
#                           type = c("gold", "grey80", "black"),
#                           labels = c("conserved module cores",
#                                      "non-conserved module portions",
#                                      "random modules")) +
#       theme(axis.ticks.y = element_blank(),
#             axis.text.y = element_blank(),
#             legend.position = "right")
#     return(p_density)
#   }
#   # lineplot
#   p_line <- ggplot(plotdf, aes(y = block_rank, x = avgCor)) + 
#     geom_line(aes(color = CCM_color, group = module_name)) +
#     geom_point(aes(color = CCM_color, shape = experiment)) +
#     #geom_text(aes(label = module_name), color = "black", size = 3, check_overlap = TRUE) + # <- add this line to see which modules are which
#     scale_color_discrete(limits = plotdf$CCM_color,
#                          type = plotdf$CCM_color) +
#     ylab("<-smallest          Module size          largest ->") +
#     xlab("Average expression correlation\n between species") +
#     theme_classic() +
#     theme(axis.ticks.y = element_blank(),
#           axis.text.y = element_blank(),
#           legend.position = "none")
#   if (.show_labels) {
#     p_line <- p_line + geom_text(aes(label = module_name, y = block_rank, x = 1.1),
#                                  size = 1)
#   }
#   # corresponding density plot
#   if (!.just_density) {
#     p_density <- ggplot(plotdf, aes(x = avgCor)) +
#       geom_density(aes(fill = modtype), alpha = 0.7) +
#       theme_classic() +
#       theme(legend.title = element_blank()) +
#       xlab("") +
#       scale_fill_discrete(limits = c("conserved",
#                                      "non-conserved",
#                                      "random"),
#                           type = c("gold", "grey80", "black"),
#                           labels = c("conserved module cores",
#                                      "non-conserved module portions",
#                                      "random modules")) +
#       theme(axis.ticks.y = element_blank(),
#             axis.text.y = element_blank(),
#             axis.ticks.x = element_blank(),
#             axis.text.x = element_blank(),
#             legend.position = "top")
#     return(list("line" = p_line, "density" = p_density))
#   }
# }
# ### Density plot for main paper, showing CCMs have more conserved 
# # expression between species versus random or divergent module portions
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/density.pdf",
#     height = 3, width = 5)
# plotModuleLineAndDensityByBlockSize(moduledf25, random_moduledf25, .just_density = TRUE)
# dev.off()
# 
# ### Supplementary figure: conserved modules having more conserved expression plasticity
# # isn't an artifact of conserved modules having larger block sizes
# # merge 00
# p_00 <- plotModuleLineAndDensityByBlockSize(moduledf00, random_moduledf00)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/modDist00.pdf",
#     width = 5, height = 12)
# ggarrange(p_00$density, p_00$line, nrow = 2, ncol = 1, widths = c(7, 7),
#           heights = c(3, 10))
# dev.off()
# 
# # merge 10
# p_10 <- plotModuleLineAndDensityByBlockSize(moduledf10, random_moduledf10)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/modDist10.pdf",
#     width = 3, height = 12)
# ggarrange(p_10$density, p_10$line, nrow = 2, ncol = 1, widths = c(7, 7),
#           heights = c(3, 10))
# dev.off()
# 
# # merge 25
# p_25 <- plotModuleLineAndDensityByBlockSize(moduledf25, random_moduledf25)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/modDist25.pdf",
#     width = 3, height = 12)
# ggarrange(p_25$density, p_25$line, nrow = 2, ncol = 1, widths = c(7, 7),
#           heights = c(3, 10))
# dev.off()
# 
# # merge 35
# p_35 <- plotModuleLineAndDensityByBlockSize(moduledf35, random_moduledf35)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/modDist35.pdf",
#     width = 3, height = 12)
# ggarrange(p_35$density, p_35$line, nrow = 2, ncol = 1, widths = c(7, 7),
#           heights = c(3, 10))
# dev.off()
# 
# #### Case example most divergent: module a in Low Pi ####
# # module 
# # selecting modules to plot
# p <- plotModuleLineAndDensityByBlockSize(filter(moduledf25, is_CCM), 
#                                          filter(random_moduledf25), .show_labels = TRUE)
# p$line
# 
# # a, purple, up in par LowPi but not in cer LowPi
# gene_idxs <- module_genedf25 |> filter(module_name == "a") |> select(gene_name) |> pull()
# conditions <- info |> filter(experiment == "LowPi") |> select(condition) |> pull()
# # TODO: make this plot two centered and scaled lineplots,
# # one for cer one for par, where each gene is a line
# getGeneCount <- function(.g, .condition, .species) {
#   if (.species == "cer") {
#     output <- counts_all2$cer[.condition, .g] |> as.numeric()
#   }
#   if (.species == "par") {
#     output <- counts_all2$par[.condition, .g] |> as.numeric()
#   }
#   if (length(output) != 1) {
#     stop("output isn't a single count")
#   }
#   return(output)
# }
# # tests for getGeneCount
# getGeneCount(.g = "YGR192C", .condition = "WT_CC_0", .species = "cer")
# 
# plotdf <- expand_grid(gene_name = gene_idxs,
#                       condition = conditions)
# plotdf$expr_cer <- map2(plotdf$gene_name, plotdf$condition, getGeneCount, .species = "cer") |> unlist()
# plotdf$expr_par <- map2(plotdf$gene_name, plotdf$condition, getGeneCount, .species = "par") |> unlist()
# plotdf_means <- plotdf |> group_by(gene_name) |> summarise(mean_expr_cer = mean(expr_cer),
#                                                            mean_expr_par = mean(expr_par),
#                                                            var_expr_cer = var(expr_cer),
#                                                            var_expr_par = var(expr_par),
#                                                            max_expr_cer = max(expr_cer),
#                                                            max_expr_par = max(expr_par))
# plotdf <- left_join(plotdf, info, by = "condition", relationship = "many-to-one") |> 
#   left_join(plotdf_means, by = "gene_name", relationship = "many-to-one")
# 
# # adding columns to say if each gene increased or decreased it's expression
# # at the final timepoint, relative to its mean expression
# plotdf <- left_join(x = plotdf,
#                     y = map(unique(plotdf$gene_name), \(g) {
#                       end_expr_cer <- plotdf |> filter(gene_name == g &
#                                                          time_point_num == 360) |> 
#                         select(expr_cer) |> pull()
#                       mean_expr_cer <- plotdf |> filter(gene_name == g) |> 
#                         select(mean_expr_cer) |> 
#                         unique() |> pull()
#                       end_expr_par <- plotdf |> filter(gene_name == g &
#                                                          time_point_num == 360) |> 
#                         select(expr_par) |> pull()
#                       mean_expr_par <- plotdf |> filter(gene_name == g) |> 
#                         select(mean_expr_par) |> 
#                         unique() |> pull()
#                       output <- tibble(gene_name = g,
#                                        upAtEnd_cer = end_expr_cer > mean_expr_cer,
#                                        upAtEnd_par = end_expr_par > mean_expr_par)
#                       return(output)
#                     }) |> bind_rows(),
#                     by = "gene_name", relationship = "many-to-one")
# 
# # removing lowly expressed genes
# plotdf <- filter(plotdf, max_expr_cer > 100 | max_expr_par > 100)
# # # centering and scaling
# # plotdf$expr_cer <- if_else(plotdf$var_expr_cer < 1,
# #                            true = 0, false = (plotdf$expr_cer - plotdf$mean_expr_cer)/plotdf$var_expr_cer)
# # plotdf$expr_par <- if_else(plotdf$var_expr_par < 1,
# #                            true = 0, false = (plotdf$expr_par - plotdf$mean_expr_par)/plotdf$var_expr_par)
# # log2 transforming
# plotdf$expr_cer <- log2(plotdf$expr_cer + 1)
# plotdf$expr_par <- log2(plotdf$expr_par + 1)
# arrange(plotdf, desc(expr_par))
# 
# p_cer_up <- ggplot(filter(plotdf, upAtEnd_cer), aes(x = time_point_num, y = expr_cer)) + 
#   geom_line(aes(group = gene_name), color = "black") + 
#   geom_line(data = group_by(filter(plotdf, upAtEnd_cer), time_point_num) |> 
#               summarise(mean_expr = mean(expr_cer)), aes(x = time_point_num,
#                                                          y = mean_expr),
#             color = "gold") +
#   theme_classic() + 
#   theme(legend.position = "none") +
#   ggtitle(paste("Up in Scer, nGenes =", sum(plotdf$upAtEnd_cer))) +
#   xlab("timepoint (min)") +
#   ylab("expression \n(log2)") + 
#   ylim(c(0, 13))
# p_cer_down <- ggplot(filter(plotdf, !upAtEnd_cer), aes(x = time_point_num, y = expr_cer)) + 
#   geom_line(aes(group = gene_name), color = "black") + 
#   geom_line(data = group_by(filter(plotdf, !upAtEnd_cer), time_point_num) |> 
#               summarise(mean_expr = mean(expr_cer)), aes(x = time_point_num,
#                                                          y = mean_expr),
#             color = "gold") +
#   theme_classic() +
#   theme(legend.position = "none") +
#   ggtitle(paste("Down in Scer, nGenes =", sum(!plotdf$upAtEnd_cer))) +
#   xlab("timepoint (min)") +
#   ylab("expression \n(log2)") + 
#   ylim(c(0, 13))
# p_par_up <- ggplot(filter(plotdf, upAtEnd_par), aes(x = time_point_num, y = expr_par)) + 
#   geom_line(aes(group = gene_name), color = "black") +
#   geom_line(data = group_by(filter(plotdf, upAtEnd_par), time_point_num) |> 
#               summarise(mean_expr = mean(expr_par)), aes(x = time_point_num,
#                                                          y = mean_expr),
#             color = "gold") +
#   theme_classic() +
#   theme(legend.position = "none") +
#   ggtitle(paste("Up in Spar, nGenes =", sum(plotdf$upAtEnd_par))) +
#   xlab("timepoint (min)") +
#   ylab("expression \n(log2)") + 
#   ylim(c(0, 13))
# p_par_down <- ggplot(filter(plotdf, !upAtEnd_par), aes(x = time_point_num, y = expr_par)) + 
#   geom_line(aes(group = gene_name), color = "black") +
#   geom_line(data = group_by(filter(plotdf, !upAtEnd_par), time_point_num) |> 
#               summarise(mean_expr = mean(expr_par)), aes(x = time_point_num,
#                                                          y = mean_expr),
#             color = "gold") +
#   theme_classic() +
#   theme(legend.position = "none") +
#   ggtitle(paste("Down in Spar, nGenes =", sum(!plotdf$upAtEnd_par))) +
#   xlab("timepoint (min)") +
#   ylab("expression \n(log2)") + 
#   ylim(c(0, 13))
# 
# pdf("../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/module_a_LowPi.pdf",
#     height = 4, width = 6)
# ggarrange(p_cer_up, p_cer_down, p_par_up, p_par_down, nrow = 2, ncol = 2)
# dev.off()
# 
# #### Expression signatures ####
# # TODO: plot every gene's closest "expression signature",
# # aka CCM average expression pattern
# # probably work best to
# # 1) decide closest signature by area between curves (could also be correlation probably)
# # 2) bipartite where each gene is connected to
# # its closest signature in each species, genes colored as conserved or divergent
# 
# # dataframe recording the average expression of each CCM
# # across all conditions in cer, par, hyc, and hyp
# signaturesdf <- expand_grid(info, module_name = filter(moduledf25, is_CCM) |> 
#                               select(module_name))
# signaturesdf$expr_cer <- apply(signaturesdf, 1, \(x) {
#   m <- x["module_name"] |> as.character()
#   cond <- x["condition"] |> as.character()
#   g_idxs <- module_genedf25 |> filter(module_name == m) |> 
#     select(gene_name) |> pull()
#   mean_expr <- counts_all2$cer[info$condition == cond, g_idxs] |> mean()
#   return(mean_expr)
# })
# signaturesdf$expr_par <- apply(signaturesdf, 1, \(x) {
#   m <- x["module_name"] |> as.character()
#   cond <- x["condition"] |> as.character()
#   g_idxs <- module_genedf25 |> filter(module_name == m) |> 
#     select(gene_name) |> pull()
#   mean_expr <- counts_all2$par[info$condition == cond, g_idxs] |> mean()
#   return(mean_expr)
# })
# signaturesdf$expr_hyc <- apply(signaturesdf, 1, \(x) {
#   m <- x["module_name"] |> as.character()
#   cond <- x["condition"] |> as.character()
#   g_idxs <- module_genedf25 |> filter(module_name == m) |> 
#     select(gene_name) |> pull()
#   mean_expr <- counts_all2_allele$cer[info$condition == cond, g_idxs] |> mean()
#   return(mean_expr)
# })
# signaturesdf$expr_hyp <- apply(signaturesdf, 1, \(x) {
#   m <- x["module_name"] |> as.character()
#   cond <- x["condition"] |> as.character()
#   g_idxs <- module_genedf25 |> filter(module_name == m) |> 
#     select(gene_name) |> pull()
#   mean_expr <- counts_all2_allele$par[info$condition == cond, g_idxs] |> mean()
#   return(mean_expr)
# })
# 
# # assigning each gene its closest signature, 
# # whichever its average expression is most correlated with
# getClosestSignatureName <- function(.g, .sigdf = signaturesdf, 
#                                     .species, .allele) {
#   if (.species == "cer") {
#     signaturesdf$gene_expr <- apply(signaturesdf, 1, \(x) {
#       cond <- x["condition"] |> as.character()
#       mean_expr <- counts_all2$cer[info$condition == cond, .g] |> mean(na.rm = TRUE)
#       return(mean_expr)
#     })
#     modcors <- signaturesdf |> group_by(module_name) |> 
#       summarise(gene_cor = cor(expr_cer, gene_expr, use = "pairwise.complete.obs"))
#   }
#   if (.species == "par") {
#     signaturesdf$gene_expr <- apply(signaturesdf, 1, \(x) {
#       cond <- x["condition"] |> as.character()
#       mean_expr <- counts_all2$par[info$condition == cond, .g] |> mean(na.rm = TRUE)
#       return(mean_expr)
#     })
#     modcors <- signaturesdf |> group_by(module_name) |> 
#       summarise(gene_cor = cor(expr_par, gene_expr, use = "pairwise.complete.obs"))
#   }
#   if (.species == "hyb" & .allele == "cer") {
#     signaturesdf$gene_expr <- apply(signaturesdf, 1, \(x) {
#       cond <- x["condition"] |> as.character()
#       mean_expr <- counts_all2_allele$cer[info$condition == cond, .g] |> mean(na.rm = TRUE)
#       return(mean_expr)
#     })
#     modcors <- signaturesdf |> group_by(module_name) |> 
#       summarise(gene_cor = cor(expr_hyc, gene_expr, use = "pairwise.complete.obs"))
#   }
#   if (.species == "hyb" & .allele == "par") {
#     signaturesdf$gene_expr <- apply(signaturesdf, 1, \(x) {
#       cond <- x["condition"] |> as.character()
#       mean_expr <- counts_all2_allele$par[info$condition == cond, .g] |> mean(na.rm = TRUE)
#       return(mean_expr)
#     })
#     modcors <- signaturesdf |> group_by(module_name) |> 
#       summarise(gene_cor = cor(expr_hyp, gene_expr, use = "pairwise.complete.obs"))
#   }
#   closest_idx <- which.max(modcors$gene_cor)
#   if (length(closest_idx) != 1) {
#     return(NA)
#   }
#   closest_modname <- unlist(modcors$module_name)[closest_idx] |> 
#     as.character()
#   return(closest_modname)
# }
# # CCM gene
# getClosestSignatureName(.g = "YDR078C", 
#                         .species = "cer", .allele = "cer")
# getClosestSignatureName(.g = "YDR078C", 
#                         .species = "par", .allele = "par")
# getClosestSignatureName(.g = "YDR078C", 
#                         .species = "hyb", .allele = "cer")
# getClosestSignatureName(.g = "YDR078C", 
#                         .species = "hyb", .allele = "par")
# module_genedf25 |> filter(gene_name == "YDR078C") |> 
#   select(module_name)
# # divergent gene
# getClosestSignatureName(.g = "YHR035W", 
#                         .species = "cer", .allele = "cer")
# getClosestSignatureName(.g = "YHR035W", 
#                         .species = "par", .allele = "par")
# module_genedf25 |> filter(gene_name == "YHR035W") |> select(CCM_color, cer_color, par_color)
# moduledf25 |> filter(module_name == "f") |> select(CCM_color, cer_color, par_color)
# moduledf25 |> filter(module_name == "e") |> select(CCM_color, cer_color, par_color)
# # gene with no expression in paradoxus
# getClosestSignatureName(.g = "YPR199C", 
#                         .species = "par", .allele = "par")
# 
# # getting closest signature for each gene in cer/par/hyc/hyp
# module_genedf25$signature_cer <- map(module_genedf25$gene_name,
#                                      getClosestSignatureName,
#                                      .species = "cer",
#                                      .allele = "cer") |> unlist()
# module_genedf25$signature_par <- map(module_genedf25$gene_name,
#                                      getClosestSignatureName,
#                                      .species = "par",
#                                      .allele = "par") |> unlist()
# module_genedf25$signature_hyc <- map(module_genedf25$gene_name,
#                                      getClosestSignatureName,
#                                      .species = "hyb",
#                                      .allele = "cer") |> unlist()
# module_genedf25$signature_hyp <- map(module_genedf25$gene_name,
#                                      getClosestSignatureName,
#                                      .species = "hyb",
#                                      .allele = "par") |> unlist()
# 
# # What proportion of CCM genes are closest to their own CCM?
# sum(module_genedf25$signature_cer[module_genedf25$is_CCM] == 
#       module_genedf25$module_name[module_genedf25$is_CCM])
# sum(module_genedf25$signature_cer[module_genedf25$is_CCM] != 
#       module_genedf25$module_name[module_genedf25$is_CCM])
# sum(module_genedf25$signature_par[module_genedf25$is_CCM] == 
#       module_genedf25$module_name[module_genedf25$is_CCM])
# sum(module_genedf25$signature_par[module_genedf25$is_CCM] != 
#       module_genedf25$module_name[module_genedf25$is_CCM]) # about 2/3 for parents
# sum(module_genedf25$signature_hyc[module_genedf25$is_CCM] == 
#       module_genedf25$module_name[module_genedf25$is_CCM])
# sum(module_genedf25$signature_hyc[module_genedf25$is_CCM] != 
#       module_genedf25$module_name[module_genedf25$is_CCM])
# sum(module_genedf25$signature_hyp[module_genedf25$is_CCM] == 
#       module_genedf25$module_name[module_genedf25$is_CCM])
# sum(module_genedf25$signature_hyp[module_genedf25$is_CCM] != 
#       module_genedf25$module_name[module_genedf25$is_CCM]) # much more even for hybrids
# 
# # what proportion of divergent genes were given different signatures
# # in cer and par?
# sum(module_genedf25$signature_cer[!module_genedf25$is_CCM] != 
#       module_genedf25$signature_par[!module_genedf25$is_CCM], na.rm = TRUE)
# sum(module_genedf25$signature_cer[!module_genedf25$is_CCM] == 
#       module_genedf25$signature_par[!module_genedf25$is_CCM], na.rm = TRUE) # about 3/4
# # in hyc and hyp?
# sum(module_genedf25$signature_hyc[!module_genedf25$is_CCM] != 
#       module_genedf25$signature_hyp[!module_genedf25$is_CCM], na.rm = TRUE)
# sum(module_genedf25$signature_hyc[!module_genedf25$is_CCM] == 
#       module_genedf25$signature_hyp[!module_genedf25$is_CCM], na.rm = TRUE) # close to 50:50
# # what proportion of signatures are simply the same for parent and hybrid alleles?
# sum(module_genedf25$signature_cer == module_genedf25$signature_par, na.rm = TRUE)
# sum(module_genedf25$signature_cer != module_genedf25$signature_par, na.rm = TRUE)
# sum(module_genedf25$signature_hyc == module_genedf25$signature_hyp)
# sum(module_genedf25$signature_hyc != module_genedf25$signature_hyp) # hybrids have more in common than parents
# 
# # TODO: where I left off: what does this mean? Parent results largely
# # make sense, most CCM genes are most correlated with their CCM
# # hybrids often are not, why would this be?
# # TODO: when hybrid CCM genes aren't most correlated with their CCM,
# # what are they correlated with? The same CCM for both alleles?
# 
# 
# # TODO: plotting, if necessary for the supplement
# # mostly we want to do this to talk about the hybrid results
# plotdf <- module_genedf25 |> 
#   select(gene_name, signature_cer, signature_par, signature_hyc, signature_hyp,
#          is_CCM, module_name, cer_color, par_color)
# plotdf$sigcolor_cer <- map(plotdf$signature_cer, \(m) {
#   moduledf25 |> filter(module_name == m) |> select(CCM_color) |> 
#     pull()
# }) |> unlist()
# plotdf$sigcolor_par <- map(plotdf$signature_par, \(m) {
#   color <- moduledf25 |> filter(module_name == m) |> select(CCM_color) |> 
#     pull()
#   if (length(color) != 1) {return(NA)}
#   return(color)
# }) |> unlist()
# plotdf$sigcolor_hyc <- map(plotdf$signature_hyc, \(m) {
#   moduledf25 |> filter(module_name == m) |> select(CCM_color) |> 
#     pull()
# }) |> unlist()
# plotdf$sigcolor_hyp <- map(plotdf$signature_hyp, \(m) {
#   moduledf25 |> filter(module_name == m) |> select(CCM_color) |> 
#     pull()
# }) |> unlist()
# # parents plot
# plotdf <- pivot_longer(plotdf, cols = c("sigcolor_cer", "sigcolor_par"),
#                        names_to = "allele", values_to = "sigcolor")
# ggplot(plotdf, aes(x = allele, y = sigcolor)) + 
#   geom_point(color = plotdf$sigcolor) +
#   geom_line(aes(color = is_CCM, group = gene_name))
# # hybrid plot
# plotdf <- pivot_longer(plotdf, cols = c("sigcolor_hyc", "sigcolor_hyp"),
#                        names_to = "allele", values_to = "sigcolor")
# 
# ggplot(plotdf, aes(x = allele, y = sigcolor)) + 
#   geom_point(color = plotdf$sigcolor) +
#   geom_line(aes(color = is_CCM, group = gene_name))
# 
# 
# #### Supplement: Visualizing module extreme subsets (divergence in expression level) ####
# # Thus far we have found that conserved co-expressed module membership
# # does not predict expression level divergence of individual genes
# # but it does tend to increase the likelihood of expression pattern 
# # conservation, how is this possible?
# 
# # First thing to check: do the most divergent genes in each CCM still
# # adhere to the expression pattern? Is there anything special about 
# # them that might explain why there's this discrepancy between
# # expression level and pattern?
# 
# coef_thresh <- 0.25
# p_thresh <- 1e-5
# 
# # Random, example CCM first (can try any merge level here, we just need the gene idxs)
# random_CCM_name <- moduledf25 |> 
#   filter(is_CCM) |>
#   select(module_name) |> 
#   pull() |> unique() |> sample(1)
# gene_idxs <- module_genedf25 |> 
#   filter(module_name == random_CCM_name) |> 
#   select(gene_name) |> pull()
# # visualizing effect sizes (average of all 4 experiments)
# spaldf |> filter(gene_name %in% gene_idxs & 
#                    coefficient == "species" &
#                    experiment != "all") |> 
#   mutate(sig = abs(effect_size) > coef_thresh & pvalue < p_thresh) |> 
#   mutate(adj_effect_size = if_else(sig, true = effect_size, false = 0)) |> 
#   group_by(gene_name) |> summarise(mean_effect = mean(adj_effect_size)) |> 
#   select(mean_effect) |> round(1) |> table()
# # seems like there's a manageable number of divergent genes on both
# # ends that we can each put them in three groups: up_cer, up_par, and conserved
# getDivergentAndConservedGeneIdxs <- function(.ccm_name) {
#   gene_idxs <- module_genedf10 |> 
#     filter(module_name == .ccm_name) |> 
#     select(gene_name) |> pull()
#   moddf <- spaldf |> filter(gene_name %in% gene_idxs & 
#                               coefficient == "species" &
#                               experiment != "all") |> 
#     mutate(sig = abs(effect_size) > coef_thresh & pvalue < p_thresh) |> 
#     mutate(adj_effect_size = if_else(sig, true = effect_size, false = 0)) |> 
#     group_by(gene_name) |> summarise(mean_effect = mean(adj_effect_size)) |> 
#     select(mean_effect, gene_name) |> 
#     mutate(direction = if_else(mean_effect < -coef_thresh,
#                                true = "up_par", 
#                                false = if_else(mean_effect > coef_thresh,
#                                                true = "up_cer", false = "conserved")))
#   table(moddf$direction)
#   up_par_idxs <- moddf |> filter(direction == "up_par") |> 
#     select(gene_name) |> pull()
#   up_cer_idxs <- moddf |> filter(direction == "up_cer") |> 
#     select(gene_name) |> pull()
#   conserved_idxs <- moddf |> filter(direction == "conserved") |> 
#     select(gene_name) |> pull()
#   return(list(conserved = conserved_idxs,
#               up_cer = up_cer_idxs,
#               up_par = up_par_idxs))
# }
# idx_list <- getDivergentAndConservedGeneIdxs(random_CCM_name)
# # up cer and conserved, comparing between species
# p <- plotExpressionProfileQuartetStack(.cts1 = counts_all2$cer[,idx_list$conserved],
#                              .cts2 = counts_all2$par[,idx_list$conserved],
#                              .cts3 = counts_all2$cer[,idx_list$up_cer],
#                              .cts4 = counts_all2$par[,idx_list$up_cer],
#                              .info1 = info, .info2 = info, .info3 = info, .info4 = info,
#                              .name1 = "S. cerevisiae, conserved",
#                              .name2 = "S. paradoxus, conserved",
#                              .name3 = "S. cerevisiae, up cer",
#                              .name4 = "S. paradoxus, up cer",
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .normalization = "log2") 
# annotate_figure(p, top = paste("Conserved and up cer genes of\n", random_CCM_name, "CCM"))
# 
# # up par and conserved, comparing between species
# p <- plotExpressionProfileQuartetStack(.cts1 = counts_all2$cer[,idx_list$conserved],
#                              .cts2 = counts_all2$par[,idx_list$conserved],
#                              .cts3 = counts_all2$cer[,idx_list$up_par],
#                              .cts4 = counts_all2$par[,idx_list$up_par],
#                              .info1 = info, .info2 = info, .info3 = info, .info4 = info,
#                              .name1 = "S. cerevisiae, conserved",
#                              .name2 = "S. paradoxus, conserved",
#                              .name3 = "S. cerevisiae, up par",
#                              .name4 = "S. paradoxus, up par",
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .normalization = "log2") 
# annotate_figure(p, top = paste("Conserved and up par genes of\n", random_CCM_name, "CCM"))
# 
# # two extremes versus each other within one species
# # cer
# plotExpressionProfilePairStack(.cts1 = counts_all2$cer[,idx_list$up_cer,drop=FALSE],
#                           .cts2 = counts_all2$cer[,idx_list$up_par,drop=FALSE],
#                           info, info,
#                           .name1 = "up cer", .name2 = "up par",
#                           .color1 = "orange1", .color2 = "orange4",
#                           .method = "line", .normalization = "log2",
#                           .show_points = TRUE) + 
#   ggtitle(paste("Divergent genes of CCM", random_CCM_name, "\nin S. cerevisiae"))
# 
# # par
# plotExpressionProfilePairStack(.cts1 = counts_all2$par[,idx_list$up_cer,drop=FALSE],
#                           .cts2 = counts_all2$par[,idx_list$up_par,drop=FALSE],
#                           info, info,
#                           .name1 = "up cer", .name2 = "up par",
#                           .color1 = "blue2", .color2 = "blue4",
#                           .method = "line", .normalization = "log2",
#                           .show_points = TRUE) + 
#   ggtitle(paste("Divergent genes of CCM", random_CCM_name, "\nin S. paradoxus"))
# 
# ### A simple story emerges: The divergent subsets are indeed expressed
# # higher or lower than the other genes, but they're clearly still
# # correlated (sometimes strongly negatively correlated, like merge25 yellow)
# # example for main figure plot
# random_CCM_name <- moduledf10 |> filter(is_CCM) |> select(module_name) |> pull() |> sample(1)
# idx_list <- getDivergentAndConservedGeneIdxs(random_CCM_name)
# p_upcer <- plotExpressionProfileQuartetStack(.cts1 = counts_all2$cer[,idx_list$conserved,drop=FALSE],
#                                         .cts2 = counts_all2$par[,idx_list$conserved,drop=FALSE],
#                                         .cts3 = counts_all2$cer[,idx_list$up_cer,drop=FALSE],
#                                         .cts4 = counts_all2$par[,idx_list$up_cer,drop=FALSE],
#                                         .info1 = info, .info2 = info, .info3 = info, .info4 = info,
#                                         .name1 = "S. cerevisiae, conserved",
#                                         .name2 = "S. paradoxus, conserved",
#                                         .name3 = "S. cerevisiae, up cer",
#                                         .name4 = "S. paradoxus, up cer",
#                                         .color1 = "orange1",
#                                         .color2 = "blue2",
#                                         .color3 = "orange4",
#                                         .color4 = "blue4",
#                                         .normalization = "log2")
# p_uppar <- plotExpressionProfileQuartetStack(.cts1 = counts_all2$cer[,idx_list$conserved,drop=FALSE],
#                                         .cts2 = counts_all2$par[,idx_list$conserved,drop=FALSE],
#                                         .cts3 = counts_all2$cer[,idx_list$up_par,drop=FALSE],
#                                         .cts4 = counts_all2$par[,idx_list$up_par,drop=FALSE],
#                                         .info1 = info, .info2 = info, .info3 = info, .info4 = info,
#                                         .name1 = "S. cerevisiae, conserved",
#                                         .name2 = "S. paradoxus, conserved",
#                                         .name3 = "S. cerevisiae, up par",
#                                         .name4 = "S. paradoxus, up par",
#                                         .color1 = "orange1",
#                                         .color2 = "blue2",
#                                         .color3 = "orange4",
#                                         .color4 = "blue4",
#                                         .normalization = "log2")
# 
# p <- ggarrange(annotate_figure(p_upcer, top = "upregulated in S. cerevisiae"),
#                annotate_figure(p_uppar, top = "upregulated in S. paradoxus"), 
#                common.legend = TRUE, nrow = 1, ncol = 2, legend = "right")
# pdf("../../aligning_the_molecular_phenotype/paper_figures/ExtremeSubsets/extreme_subset_example.pdf",
#     width = 7, height = 14)
# annotate_figure(p, top = paste("extreme subsets of module", random_CCM_name))
# dev.off()
# ### Supplemental plots: all CCMs with extreme subsets indicated
# # Just merge 10 for now because this quickly gets out of hand
# CCM_names10 <- moduledf10 |> filter(is_CCM) |> select(module_name) |> pull() |> unique()
# plotlist_upcer <- vector(mode = "list", length = length(CCM_names10))
# names(plotlist_upcer) <- CCM_names10
# plotlist_uppar <- vector(mode = "list", length = length(CCM_names10))
# names(plotlist_uppar) <- CCM_names10
# for (ccm in CCM_names10) {
#   cat("starting on", ccm, "\n")
#   idx_list <- getDivergentAndConservedGeneIdxs(ccm)
#   if (length(idx_list$up_cer > 0)) {
#     p <- plotExpressionProfileQuartetStack(.cts1 = counts_all2$cer[,idx_list$conserved,drop=FALSE],
#                                                           .cts2 = counts_all2$par[,idx_list$conserved,drop=FALSE],
#                                                           .cts3 = counts_all2$cer[,idx_list$up_cer,drop=FALSE],
#                                                           .cts4 = counts_all2$par[,idx_list$up_cer,drop=FALSE],
#                                                           .info1 = info, .info2 = info, .info3 = info, .info4 = info,
#                                                           .name1 = "S. cerevisiae, conserved",
#                                                           .name2 = "S. paradoxus, conserved",
#                                                           .name3 = "S. cerevisiae, up cer",
#                                                           .name4 = "S. paradoxus, up cer",
#                                                           .color1 = "orange1",
#                                                           .color2 = "blue2",
#                                                           .color3 = "orange4",
#                                                           .color4 = "blue4",
#                                                           .normalization = "log2")
#     plotlist_upcer[[ccm]] <- annotate_figure(p, top = ccm)
#   }
#   if (length(idx_list$up_par) > 0) {
#     p <- plotExpressionProfileQuartetStack(.cts1 = counts_all2$cer[,idx_list$conserved,drop=FALSE],
#                                                           .cts2 = counts_all2$par[,idx_list$conserved,drop=FALSE],
#                                                           .cts3 = counts_all2$cer[,idx_list$up_par,drop=FALSE],
#                                                           .cts4 = counts_all2$par[,idx_list$up_par,drop=FALSE],
#                                                           .info1 = info, .info2 = info, .info3 = info, .info4 = info,
#                                                           .name1 = "S. cerevisiae, conserved",
#                                                           .name2 = "S. paradoxus, conserved",
#                                                           .name3 = "S. cerevisiae, up par",
#                                                           .name4 = "S. paradoxus, up par",
#                                                           .color1 = "orange1",
#                                                           .color2 = "blue2",
#                                                           .color3 = "orange4",
#                                                           .color4 = "blue4",
#                                                           .normalization = "log2")
#     plotlist_uppar[[ccm]] <- annotate_figure(p, top = ccm)
#   }
# }
# # upregulated in cer
# p_upcer <- ggarrange(plotlist = plotlist_upcer, 
#                    nrow = ceiling(length(plotlist_upcer)/3), 
#                    ncol = 3, 
#                    common.legend = TRUE, 
#                    legend = "right")
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/supp_upcer_profiles.pdf",
#     width = 12, height = 4*length(plotlist_upcer))
# annotate_figure(p_upcer, left = textGrob("Expression (log scale)", rot = 90, vjust = 1, gp = gpar(cex = 3)),
#                 bottom = textGrob("Timepoint", gp = gpar(cex = 3)),
#                 top = textGrob("Conserved vs. upregulated \nin S. cerevisiae", gp = gpar(cex = 3)))
# dev.off()
# # upregulated in par
# p_uppar <- ggarrange(plotlist = plotlist_uppar, 
#                      nrow = ceiling(length(plotlist_uppar)/3), 
#                      ncol = 3, 
#                      common.legend = TRUE, 
#                      legend = "right")
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/supp_uppar_profiles.pdf",
#     width = 12, height = 4*length(plotlist_uppar))
# # TODO: add expression correlations for comparison to quantification below
# annotate_figure(p_uppar, left = textGrob("Expression (log scale)", rot = 90, vjust = 1, gp = gpar(cex = 3)),
#                 bottom = textGrob("Timepoint", gp = gpar(cex = 3)),
#                 top = textGrob("Conserved vs. upregulated \nin S. paradoxus", gp = gpar(cex = 3)))
# dev.off()
# 
# # Jitterish plot to show that extreme subset is
# # "still quite correlated" with rest of module
# # (no comparison group, but we're expecting them all to be quite high aka > 0.8)
# # use getExpressionSimilarity to get the avg expr correlation
# # for each extreme subset with the conserved portion of its module
# # y axis: avg expression correlation between extreme subset and rest of module
# # x axis: 4 groups: up cer, in cer, up par in cer, up cer in par, up par in par
# # and indicate where module b is
# plotdf <- tibble(module_name = rep(CCM_names10, 2),
#                  direction = c(rep("up_cer", length(CCM_names10)),
#                                rep("up_par", length(CCM_names10))))
# plotdf$avg_cor_cer <- map2(plotdf$module_name, plotdf$direction, \(ccm, d) {
#   idx_list <- getDivergentAndConservedGeneIdxs(ccm)
#   if (d == "up_cer") {
#     extremes <- idx_list$up_cer
#   }
#   if (d == "up_par") {
#     extremes <- idx_list$up_par
#   }
#   if (length(extremes) == 0) {
#     return(NA)
#   }
#   avgcors_cer <- getExpressionSimilarity(counts_all2$cer[,idx_list$conserved],
#                                          counts_all2$cer[,extremes])
#   return(cor(avgcors_cer$avg1, avgcors_cer$avg2))
# }) |> unlist()
# plotdf$avg_cor_par <- map2(plotdf$module_name, plotdf$direction, \(ccm, d) {
#   idx_list <- getDivergentAndConservedGeneIdxs(ccm)
#   if (d == "up_cer") {
#     extremes <- idx_list$up_cer
#   }
#   if (d == "up_par") {
#     extremes <- idx_list$up_par
#   }
#   if (length(extremes) == 0) {
#     return(NA)
#   }
#   avgcors_par <- getExpressionSimilarity(counts_all2$par[,idx_list$conserved],
#                                          counts_all2$par[,extremes])
#   return(cor(avgcors_par$avg1, avgcors_par$avg2))
# }) |> unlist()
# plotdf <- plotdf |> pivot_longer(cols = c("avg_cor_cer", "avg_cor_par"),
#                                  names_to = "species", names_prefix = "avg_cor_",
#                                  values_to = "avg_cor")
# plotdf <- left_join(plotdf, select(moduledf10, module_name, CCM_color),
#                     by = "module_name")
# plotdf$label = map2(plotdf$direction, plotdf$species, \(x, y) {
#   if (x == "up_cer" & y == "cer") {
#     return("upregulated in S. cerevisiae\ncorrelation in S. cerevisiae")
#   }
#   if (x == "up_cer" & y == "par") {
#     return("upregulated in S. cerevisiae\ncorrelation in S. paradoxus")
#   }
#   if (x == "up_par" & y == "cer") {
#     return("upregulated in S. paradoxus\ncorrelation in S. cerevisiae")
#   }
#   if (x == "up_par" & y == "par") {
#     return("upregulated in S. paradoxus\ncorrelation in S. paradoxus")
#   }
# }) |> unlist()
# plotdf$module_name <- gsub("ccm", "", plotdf$module_name)
# pdf("../../aligning_the_molecular_phenotype/paper_figures/ExtremeSubsets/avg_cor_extremes.pdf",
#     width = 9, height = 4)
# ggplot(plotdf, aes(x = label, y = avg_cor)) +
#   geom_jitter(aes(color = CCM_color), size = 7,
#               position = position_jitter(seed = 1),
#               alpha = 0.75) +
#   scale_color_discrete(limits = plotdf$CCM_color,
#                        type = plotdf$CCM_color) +
#   geom_text(aes(label = module_name),
#             check_overlap = TRUE, size = 5,
#              position = position_jitter(seed = 1)) +
#   theme_classic() +
#   theme(legend.position = "none") +
#   ylim(c(0, 1)) +
#   xlab("") +
#   ylab("average expression correlation\nbetween extreme subset and rest of module")
# dev.off()
# #### Pattern divergence - CCM vs rest of cluster ####
# 
# # TODO: This is a clever idea, but currently it's not clear how much pattern divergence
# # would be considered significantly different
# # example of why this is a problem:
# # the rest of module sets of genes tend to have quite similar patterns between species,
# # why were they only grouped in the same module in one species?
# # it *appears* it's because the ccm is actually diverging between species, but
# # because we don't have a cutoff for when an expression pattern is significantly different,
# # we can't say this is the case for sure
# # maybe we could randomly divide the ccm in half a few times to create an expectation
# # for how much genes within a highly correlated block are allowed to differ?
# 
# # Extreme subsets have given us an idea of how expression level can
# # diverge with relatively little of an effect on expression pattern
# 
# # But is the opposite ever true? Can pattern diverge while expression level
# # is relatively unaffected?
# 
# # We have a way to look at this---the co-expressed modules have portions
# # that are conserved between species and portions that are unique
# # to one species
# 
# # For each CCM, we'll compare the expression profiles of the conserved co-
# # expressed genes in each species with the genes that were only
# # grouped into that same module in one species or the other.
# # Presumably the species-unique portions will have slightly
# # different expression profiles
# 
# # let's start with a random CCM
# random_CCM_name <- moduledf25 |> filter(is_CCM) |> select(module_name) |> pull() |> sample(1)
# random_cer_color <- moduledf25 |> filter(module_name == random_CCM_name) |> select(cer_color) |> pull()
# random_par_color <- moduledf25 |> filter(module_name == random_CCM_name) |> select(par_color) |> pull()
# 
# ccm_idxs <- module_genedf25 |> filter(module_name == random_CCM_name) |> select(gene_name) |> pull()
# restof_idxs_cer <- module_genedf25 |> filter(cer_color == random_cer_color) |> select(gene_name) |> pull()
# restof_idxs_par <- module_genedf25 |> filter(par_color == random_par_color) |> select(gene_name) |> pull()
# 
# # CCMs and cer-unique genes
# plotExpressionProfileQuartetStack(.cts1 = counts_all2$cer[,ccm_idxs],
#                              .cts2 = counts_all2$par[,ccm_idxs],
#                              .cts3 = counts_all2$cer[,restof_idxs_cer, drop = FALSE],
#                              .cts4 = counts_all2$par[,restof_idxs_cer, drop = FALSE],
#                              .info1 = info,
#                              .info2 = info,
#                              .info3 = info,
#                              .info4 = info,
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .name1 = "S. cerevisiae, CCM",
#                              .name2 = "S. paradoxus, CCM",
#                              .name3 = "S. cerevisiae, rest of module",
#                              .name4 = "S. cerevisiae, rest of module",
#                              .normalization = "scale")
# # CCMs and par-unique genes
# plotExpressionProfileQuartetStack(.cts1 = counts_all2$cer[,ccm_idxs],
#                                   .cts2 = counts_all2$par[,ccm_idxs],
#                                   .cts3 = counts_all2$cer[,restof_idxs_par, drop = FALSE],
#                                   .cts4 = counts_all2$par[,restof_idxs_par, drop = FALSE],
#                                   .info1 = info,
#                                   .info2 = info,
#                                   .info3 = info,
#                                   .info4 = info,
#                                   .color1 = "orange1",
#                                   .color2 = "blue2",
#                                   .color3 = "orange4",
#                                   .color4 = "blue4",
#                                   .name1 = "S. cerevisiae, CCM",
#                                   .name2 = "S. paradoxus, CCM",
#                                   .name3 = "S. cerevisiae, rest of module",
#                                   .name4 = "S. cerevisiae, rest of module",
#                                   .normalization = "scale")
# 
# # Same thing with hybrid alleles
# plotExpressionProfileQuartetStack(.cts1 = counts_all2_allele$cer[,ccm_idxs],
#                              .cts2 = counts_all2_allele$par[,ccm_idxs],
#                              .cts3 = counts_all2_allele$cer[,restof_idxs_cer, drop = FALSE],
#                              .cts4 = counts_all2_allele$par[,restof_idxs_cer, drop = FALSE],
#                              .info1 = info,
#                              .info2 = info,
#                              .info3 = info,
#                              .info4 = info,
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .name1 = "Hybrid cerevisiae allele, CCM",
#                              .name2 = "Hybrid paradoxus allele, CCM",
#                              .name3 = "Hybrid cerevisiae allele, rest of cerevisiae module",
#                              .name4 = "Hybrid paradoxus allele, rest of cerevisiae module",
#                              .normalization = "scale")
# plotExpressionProfileQuartetStack(.cts1 = counts_all2_allele$cer[,ccm_idxs],
#                                   .cts2 = counts_all2_allele$par[,ccm_idxs],
#                                   .cts3 = counts_all2_allele$cer[,restof_idxs_par, drop = FALSE],
#                                   .cts4 = counts_all2_allele$par[,restof_idxs_par, drop = FALSE],
#                                   .info1 = info,
#                                   .info2 = info,
#                                   .info3 = info,
#                                   .info4 = info,
#                                   .color1 = "orange1",
#                                   .color2 = "blue2",
#                                   .color3 = "orange4",
#                                   .color4 = "blue4",
#                                   .name1 = "Hybrid cerevisiae allele, CCM",
#                                   .name2 = "Hybrid paradoxus allele, CCM",
#                                   .name3 = "Hybrid cerevisiae allele, rest of paradoxus module",
#                                   .name4 = "Hybrid paradoxus allele, rest of paradoxus module",
#                                   .normalization = "scale")
# # hybrid alleles should be more similar than species orthologs were,
# # but do they revert to be more similar to the parental expression 
# # or not is the question
# 
# # Does this affect how genes are grouped in the hybrid?
# # In cerevisiae, these divergent idxs are grouped with the ccm:
# colors25$cer[colnames(counts_all2$cer) %in% c(divergent_idxs, ccm_idxs)] |> table()
# # In paradoxus, they're grouped separate:
# colors25$par[colnames(counts_all2$par) %in% ccm_idxs] |> table()
# colors25$par[colnames(counts_all2$par) %in% divergent_idxs] |> table()
# # In hyc, are they grouped the same?
# colors25$hyc[colnames(counts_all2$cer) %in% divergent_idxs] |> table()
# # how about hyp?
# colors25$hyp[colnames(counts_all2$cer) %in% divergent_idxs] |> table()
# 
# #### Expression pattern vs level ####
# 
# # The main plot we want to get from this is an illustration of how
# # level and pattern are independent properties and that
# # pattern is more correlated between species/alleles in conserved, co-expressed genes
# 
# # (But also, correlation is a finicky measure for individual genes,
# # one of the reasons WGCNA builds networks on topological overlap rather 
# # than pairwise correlations alone, so it might be best to only
# # compare correlations of module avg expression like we do in the cis_trans script)
# 
# # Before we conclude that it's easier to evolve changes in expression
# # level than pattern, we need to look at what's happening to non-CCM
# # genes. How does their degree of expression level and pattern divergence
# # (as measured by shape correlation and effect size respectively)
# # compare to what's happening with the CCMs?
# 
# # plotting mean effect size (maybe magnitude so divergent are
# # all on one side) versus shape correlation for each individual
# # gene and indicate which genes are part of CCMs and which are
# # divergent (and maybe which are the grey CCM? If they seem special)
# 
# # first step is to get a divergentdf where each gene has 4 rows, 1 for each experiment
# # with 3 columns: color in each species, mean effect size, and shape correlation
# 
# # making divergentdf where each gene is one row and it has both the parent and hybrid cor and effect sizes,
# # averaged across all 4 experiments as 4 separate columns
# divergentdf <- spaldf |> filter(experiment != "all") |> 
#   mutate(lfc = if_else(pvalue < p_thresh & 
#                                      abs(effect_size) > coef_thresh,
#                                    true = effect_size,
#                                    false = 0)) |> 
#   select(lfc, gene_name, experiment, coefficient) |> 
#   right_join(y = select(module_genedf10, gene_name, cer_color, par_color, coexpressed),
#              by = "gene_name", relationship = "many-to-one")
# divergentdf$experiment <- gsub("TFdelx", "", divergentdf$experiment) # I hate that I did this at some point
# divergentdf$shape_cor <- apply(divergentdf, 1, \(x) {
#   e <- x["experiment"] |> as.character()
#   g <- x["gene_name"] |> as.character()
#   if (as.character(x["coefficient"]) == "species") {
#     exprSim <- getExpressionSimilarity(counts_all2$cer[info$experiment == e,g,drop=FALSE],
#                                        counts_all2$par[info$experiment == e,g,drop=FALSE])
#   }
#   if (as.character(x["coefficient"]) == "allele") {
#     exprSim <- getExpressionSimilarity(counts_all2_allele$cer[info$experiment == e,g,drop=FALSE],
#                                        counts_all2_allele$par[info$experiment == e,g,drop=FALSE])
#   }
#   return(cor(exprSim$avg1, exprSim$avg2))
# }) |> unlist()
# divergentdf <- pivot_wider(divergentdf, id_cols = c("gene_name", "cer_color", "par_color", "experiment", "coexpressed"),
#                            names_from = c("coefficient"), values_from = c("lfc", "shape_cor")) |> 
#   group_by(gene_name, cer_color, par_color, coexpressed) |> 
#   summarise(mean_lfc_parent = mean(lfc_species),
#             mean_shape_cor_parent = mean(shape_cor_species),
#             mean_lfc_hybrid = mean(lfc_allele),
#             mean_shape_cor_hybrid = mean(shape_cor_allele))
# library(ggExtra)
# max_lfc <- max(c(abs(divergentdf$mean_lfc_parent),
#                  abs(divergentdf$mean_lfc_hybrid)), na.rm = TRUE)
# p_avgcor_effsize_parents <- ggplot(divergentdf, aes(x = mean_shape_cor_parent, y = mean_lfc_parent)) +
#   geom_point(aes(color = coexpressed)) + 
#   scale_color_discrete(type = c("gold", "mediumseagreen", "orange1", "blue2", "grey"),
#                        limits = c("conserved co-expressed",
#                                   "diverged co-expressed",
#                                   "S. cerevisiae co-expressed",
#                                   "S. paradoxus co-expressed",
#                                   "never co-expressed")) +
#   theme_classic() +
#   theme(legend.title = element_blank()) +
#   xlab("expression pattern correlation between species") +
#   ylab("log fold change between species") +
#   xlim(c(-0.5, 1)) +
#   ylim(c(-max_lfc, max_lfc))
# p_parents <- ggMarginal(p_avgcor_effsize_parents, groupColour = TRUE, groupFill = TRUE)
# p_avgcor_effsize_hyb <- ggplot(divergentdf, aes(x = mean_shape_cor_hybrid, y = mean_lfc_hybrid)) +
#   geom_point(aes(color = coexpressed)) + 
#   scale_color_discrete(type = c("gold", "mediumseagreen", "orange1", "blue2", "grey"),
#                        limits = c("conserved co-expressed",
#                                   "diverged co-expressed",
#                                   "S. cerevisiae co-expressed",
#                                   "S. paradoxus co-expressed",
#                                   "never co-expressed")) +
#   theme_classic() +
#   theme(legend.title = element_blank()) +
#   xlab("expression pattern correlation between hybrid alleles") +
#   ylab("log fold change") +
#   xlim(c(-0.5, 1)) +
#   ylim(c(-max_lfc, max_lfc))
# p_hyb <- ggMarginal(p_avgcor_effsize_hyb, groupColour = TRUE, groupFill = TRUE)
# pdf("../../aligning_the_molecular_phenotype/paper_figures/PatternVsLevel/single_gene_avg_cor_effsizes_parents.pdf",
#     width = 6, height = 6)
# p_parents
# dev.off()
# pdf("../../aligning_the_molecular_phenotype/paper_figures/PatternVsLevel/single_gene_avg_cor_effsizes_hybrid.pdf",
#     width = 6, height = 6)
# p_hyb
# dev.off()
# 
# #### Hybrid expression similarity improvement ####
# 
# # What does the distribution of correlation IMPROVEMENT look like for genes of different coexpression groups?
# # Does every CCM gene tend to have slightly more improvement than every non-CCM gene?
# # Or is there a subset of CCM genes that have REALLY strong improvement (like jump from 0.5 to 1)
# # and there isn't really the same sort of subset of divergent genes?
# 
# # To address this, we will look at the improvement in correlation for each gene,
# # paying attention to if the gene is a CCM gene or divergent, and if CCM, which CCM it's part of
# 
# # Pairing every gene's parental and hybrid correlation to see how many are 
# # improving versus not improving
# p_improve_shapecor <- ggplot(divergentdf, aes(x = mean_shape_cor_parent, y = mean_shape_cor_hybrid)) + 
#   geom_point(aes(color = coexpressed)) +
#   scale_color_discrete(type = c("gold", "mediumseagreen", "orange1", "blue2", "grey"),
#                        limits = c("conserved co-expressed",
#                                   "diverged co-expressed",
#                                   "S. cerevisiae co-expressed",
#                                   "S. paradoxus co-expressed",
#                                   "never co-expressed")) +
#   theme_classic() +
#   xlim(c(-0.5, 1)) +
#   ylim(c(-0.5, 1)) +
#   geom_abline(slope = 1, intercept = 0) +
#   xlab("scaled correlation - parents") +
#   ylab("scaled correlation - hybrid alleles")
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/scaled_cor_improvement.pdf",
#     width = 6, height = 5)
# ggMarginal(p_improve_shapecor, groupColour = TRUE, groupFill = TRUE)
# dev.off()
# 
# # repeat for lfc
# p_improve_lfc <- ggplot(divergentdf, aes(x = mean_lfc_parent, y = mean_lfc_hybrid)) + 
#   geom_point(aes(color = coexpressed)) +
#   geom_text(data = filter(divergentdf, abs(mean_lfc_parent) > 2.5 | abs(mean_lfc_hybrid) > 2.5),
#             aes(label = gene_name), check_overlap = TRUE, nudge_x = 1.2) +
#   scale_color_discrete(type = c("gold", "mediumseagreen", "orange1", "blue2", "grey"),
#                        limits = c("conserved co-expressed",
#                                   "diverged co-expressed",
#                                   "S. cerevisiae co-expressed",
#                                   "S. paradoxus co-expressed",
#                                   "never co-expressed")) +
#   theme_classic() +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_hline(yintercept = 0) +
#   xlab("log fold change - parents") +
#   ylab("log fold change - hybrid alleles")
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/lfc_improvement.pdf",
#     width = 6, height = 5)
# ggMarginal(p_improve_lfc, groupColour = TRUE, groupFill = TRUE)
# dev.off()
# 
# # Is level and pattern improvement related?
# plotdf <- divergentdf |> mutate()

#### Separating out module shadows ####
# 
# # Archived because we switched to just straight-up correlation clustering using
# # non-abs value Pearson cor (but this code was instrumental in figuring out
# # how to do that and why it might be useful!)
# 
# # Issue: clustering is more informative using correlation magnitude, ignoring direction
# # but if I want to see the average expression of a module and some genes are negatively correlated
# # while others are positive correlated, and average can suppress this.
# # So I first need to check how many, if any, genes are negatively correlated with most of their module
# 
# # helper function for getModuleShadowGenes
# # iterates to minimize 
# splitCorMat <- function(.cor_mat) {
#   shadow_tree <- hclust(as.dist(-.cor_mat), method = "average") # negative because hclust expects distance mat --- the higher the ij value the LESS similar genes i and j are
#   topBranchHeight <- sort(shadow_tree$height, decreasing = TRUE)[2]
#   shadow_labels <- cutreeStatic(shadow_tree, cutHeight = topBranchHeight, 
#                                 minSize = 0.05*nrow(.cor_mat))
#   if (length(unique(shadow_labels)) != 2) {
#     stop("more than two groups identified!\n")
#   }
#   genes1 <- rownames(.cor_mat)[shadow_labels == 1]
#   genes2 <- rownames(.cor_mat)[shadow_labels != 1]
#   if (length(genes1) > length(genes2)) {
#     shadow_genes <- genes2
#   }
#   if (length(genes1) <= length(genes2)) {
#     shadow_genes <- genes1
#   }
#   return(shadow_genes)
# }
# getModuleShadowGenes <- function(.gene_idxs, .species, 
#                                  .experiment = unique(info$experiment), 
#                                  .plotToo = FALSE) {
#   # TODO: drop na genes in Heat/Cold
#   shadow_genes <- cor(counts_all2[[.species]][info$experiment %in% .experiment, .gene_idxs],
#                       use = "pairwise.complete.obs") |> splitCorMat()
#   if (.plotToo) {
#     outputdf <- as_tibble(t(counts_all2[[.species]][, .gene_idxs])) |> 
#       bind_cols(tibble(gene_name = .gene_idxs)) |> 
#       pivot_longer(cols = rownames(counts_all2[[.species]]),
#                    names_to = "condition", values_to = "expr") |> 
#       left_join(info, by = "condition") |> 
#       filter(experiment %in% .experiment) 
#     outputdf <- left_join(x = outputdf, 
#                           y = group_by(outputdf, gene_name) |> 
#                             summarise(mean_expr = mean(expr),
#                                       sd_expr = sd(expr)),
#                           by = "gene_name", relationship = "many-to-one") |> 
#       mutate(scaled_expr = (expr - mean_expr)/sd_expr)
#     outputdf$isShadow <- outputdf$gene_name %in% shadow_genes
#     return(outputdf)
#   } 
#   if (!.plotToo) {
#     return(shadow_genes)
#   }
# }
# ### tests for splitCorMat and getModuleShadowGenes
# # tests for splitCorMat
# # module m in LowN, should split these into genes highest at TP2 (^ genes) and genes lowest at TP2 (v genes)
# toy_idxs <- c("YBR083W", "YBR172C", "YML015C", # v genes
#               "YBR162W-A", "YKL196C", "YBR171W") # ^ genes (the YBRs are coincidental, as I was just scrolling through that part of the module --- although the 171 172 W/C genes are probably overlapping)
# toy_mat <- cor(counts_all2[["par"]][info$experiment == "LowN", toy_idxs])
# toy_tree <- hclust(-as.dist(toy_mat))
# plot(toy_tree) # v genes should be separated from ^ genes
# toy_shadow_genes <- splitCorMat(toy_mat)
# toy_shadow_genes # should be all on one branch of the tree
# # now full module test
# random_module <- module_genedf10 |> filter(is_CCM) |> select(module_name) |> pull() |> unique() |> sample(1)
# random_experiment <- sample(unique(info$experiment), 1)
# random_species <- sample(c("cer", "par"), 1)
# gene_idxs <- module_genedf10 |> filter(module_name == random_module) |> select(gene_name) |> pull()
# test_mat <-  cor(counts_all2[[random_species]][info$experiment == random_experiment, gene_idxs])
# sum(is.na(test_mat))
# test_tree <- hclust(-as.dist(test_mat))
# test_height <- sort(test_tree$height, decreasing = TRUE)[2]
# test_labels <- cutreeDynamicTree(test_tree, 
#                                  maxTreeHeight = 0.999,
#                                  deepSplit = FALSE)
# plot(test_tree, labels = FALSE)
# labelsdf <- tibble(gene_name = rownames(test_mat),
#                    label = test_labels)
# testdf <- getModuleShadowGenes(gene_idxs, .species = "par", .experiment = "LowN", .plotToo = TRUE)
# testdf <- left_join(testdf, labelsdf, by = "gene_name", relationship = "many-to-one")
# ggplot(testdf, aes(x = time_point_str, y = scaled_expr)) + 
#   geom_line(aes(group = gene_name)) +
#   facet_wrap(~label)
# 
# # tests for getModuleShadowGenes
# # module m in LowN has most genes peaking at TP2 (^ shape) and a subset peaking at TP1 or TP3 (V shape, or \ / shapes)
# gene_idxs <- module_genedf10 |> filter(module_name == "m") |> select(gene_name) |> pull()
# # cer
# test <- getModuleShadowGenes(.gene_idxs = gene_idxs, .species = "cer",
#                              .plotToo = TRUE)
# ggplot(test, aes(x = time_point_num, y = scaled_expr)) + 
#   geom_line(aes(group = gene_name, color = isShadow)) +
#   facet_wrap(~experiment, scales = "free")
# ggplot(filter(test, experiment == "LowN"),
#        aes(x = time_point_num, y = scaled_expr)) + 
#   geom_line(aes(group = gene_name, color = isShadow)) +
#   facet_wrap(~isShadow)
# 
# # par
# test <- getModuleShadowGenes(.gene_idxs = gene_idxs, .species = "par", 
#                              .experiment = "LowN", .plotToo = TRUE)
# ggplot(test, aes(x = time_point_str, y = scaled_expr)) + 
#   geom_line(aes(group = gene_name)) +
#   facet_wrap(~isShadow)
# 
# # TODO: applying to all modules
# shadowdf00 <- 
#   for (m in module_genedf00$module_name) {
#     gene_idxs <- module_genedf00 |> filter(module_name == m) |> select(gene_name) |> pull()
#     isShadow_cer <- getModuleShadowGenes(.gene_idxs = gene_idxs, .species = "cer", 
#                                          .plotToo = FALSE)
#   }
# 
# 
# # TODO: where I left off on this: this isn't clearly separating + and - correlated genes,
# # which may mean there aren't many instances where it matters, so we can safely just use
# # module average expression as the signature.
# # We could adapt this into a QC measure to create module avgExpr based on shadow vs not
# # shadow then plot just to see which modules are that different
# 

#### Module line plots ####
# # goal was to show how much ccm expression conservation varied between
# # environments
# # archived b/c no one liked how the line plots looked, too hard to follow one module
# plotModuleLines <- function(.mdf, .none_color = "black", .show_labels = FALSE) {
#   plotdf <- .mdf |> 
#     select(setdiff(colnames(.mdf), "avgCor")) |> 
#     pivot_longer(cols = c("avgCor_CC",
#                           "avgCor_LowN",
#                           "avgCor_LowPi",
#                           "avgCor_HAP4",
#                           "avgCor_Heat",
#                           "avgCor_Cold"),
#                  names_to = "experiment",
#                  values_to = "avgCor",
#                  names_prefix = "avgCor_") |> 
#     group_by(module_name) |> 
#     mutate(avgCor_rank = rank(-avgCor))
#   plotdf <- plotdf |> filter(block_size >= 30)
#   plotdf$CCM_color <- gsub("none", .none_color, plotdf$CCM_color)
#   # lineplot
#   p_line <- ggplot(plotdf, aes(y = avgCor, x = experiment)) + 
#     geom_line(aes(color = CCM_color, group = module_name)) +
#     geom_point(aes(color = CCM_color)) +
#     scale_color_discrete(type = plotdf$CCM_color, 
#                          limits = plotdf$CCM_color) +
#     scale_x_discrete(breaks = c("CC", "HAP4", "LowN", "LowPi", "Heat", "Cold"),
#                      limits = c("CC", "HAP4", "LowN", "LowPi", "Heat", "Cold"),
#                      labels = c("Urea Shock", "Diauxic Shift", "Low Nitrogen",
#                                 "Low Phosphorus", "Heat Shock", "Cold Shock")) +
#     ylab("Module similarity between species\n (average expression correlation)") +
#     xlab("") +
#     ylim(c(-1, 1)) +
#     theme_classic() +
#     theme(legend.position = "none")
#   if (.show_labels) {
#     p_line <- p_line + geom_text(aes(label = module_name, y = avgCor, x = 1.1),
#                                  size = 1)
#   }
#   return(p_line)
# }
# 
# # Supplementary figure: other module levels also have same avgCor null patterns
# # merge 00
# p_00 <- plotModuleLines(filter(moduledf00, is_CCM))
# p_n00 <- plotModuleLines(filter(moduledf00, !is_CCM), .none_color = "grey80")
# p_r00 <- plotModuleLines(random_moduledf00, .none_color = "black")
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/modLines00.pdf",
#     height = 10, width = 10)
# ggarrange(p_00, p_n00, p_r00, nrow = 3, ncol = 1)
# dev.off()
# 
# # merge 10
# p_10 <- plotModuleLines(filter(moduledf10, is_CCM))
# p_n10 <- plotModuleLines(filter(moduledf10, !is_CCM), .none_color = "grey80")
# p_r10 <- plotModuleLines(random_moduledf10, .none_color = "black")
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/modLines10.pdf",
#     height = 10, width = 10)
# ggarrange(p_10, p_n10, p_r10, nrow = 3, ncol = 1)
# dev.off()
# 
# # merge 25
# p_25 <- plotModuleLines(filter(moduledf25, is_CCM))
# p_n25 <- plotModuleLines(filter(moduledf25, !is_CCM), .none_color = "grey80")
# p_r25 <- plotModuleLines(random_moduledf25, .none_color = "black")
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/modLines25.pdf",
#     height = 10, width = 10)
# ggarrange(p_25, p_n25, p_r25, nrow = 3, ncol = 1)
# dev.off()
# 
# # merge 35
# p_35 <- plotModuleLines(filter(moduledf35, is_CCM))
# p_n35 <- plotModuleLines(filter(moduledf35, !is_CCM), .none_color = "grey80")
# p_r35 <- plotModuleLines(random_moduledf35, .none_color = "black")
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/modLines35.pdf",
#     height = 10, width = 10)
# ggarrange(p_35, p_n35, p_r35, nrow = 3, ncol = 1)
# dev.off()

#### Comparing expression divergence quantification methods ####
# # QC: How our different measures of expression similarity relate to each other
# # 1) AvgCor: correlation between each species' module's average expression vector---the closer to 1, the more similar
# # 2) EigenDist: euclidean distance between each species' module's first eigenvector---the closer to 0, the more similar
# # 3) Area Between Curves: summed absolute difference between each member of the average expression vectors, the closer to 0, the more similar
# # 4) Distance Correlation: should be more robust to lowly expressed modules
# 
# # AvgCor vs dcor
# # overall
# plotdf <- moduledf10 |> 
#   filter(is_CCM) |> 
#   select(CCM_color, dcor, avgCor)
# ggplot(plotdf, aes(x = avgCor, y = dcor)) + 
#   geom_point(color = plotdf$CCM_color) + 
#   xlab("Pearson correlation") + ylab("distance correlation") + 
#   ggtitle("comparison of two measures of module\n average expression correlation between species") +
#   scale_shape_discrete(limits = c("CC", "HAP4", "LowN", "LowPi"),
#                        labels = c("Urea Shock", "Diauxic Shift", "Low Nitrogen", "Low Phosphorus")) +
#   theme_classic() + xlim(c(-1, 1)) + ylim(c(0, 1)) + geom_abline(slope = 0.5, intercept = 0.5, color = "gold")
# 
# # experiment-specific
# plotdf_avgCor <- moduledf10 |> 
#   filter(is_CCM) |> 
#   select(avgCor_CC, avgCor_HAP4, avgCor_LowN, avgCor_LowPi, CCM_color) |> 
#   pivot_longer(cols = c("avgCor_CC", "avgCor_HAP4", "avgCor_LowN", "avgCor_LowPi"),
#                names_to = "experiment", values_to = "avgCor", names_prefix = "avgCor_")
# plotdf_dcor <- moduledf10 |> 
#   filter(is_CCM) |> 
#   select(dcor_CC, dcor_HAP4, dcor_LowN, dcor_LowPi, CCM_color) |> 
#   pivot_longer(cols = c("dcor_CC", "dcor_HAP4", "dcor_LowN", "dcor_LowPi"),
#                names_to = "experiment", values_to = "dcor", names_prefix = "dcor_")
# plotdf <- left_join(plotdf_avgCor, plotdf_dcor, by = c("experiment", "CCM_color"))
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/avgCor_vs_dcor.pdf",
#     width = 5, height = 3)
# ggplot(plotdf, aes(x = avgCor, y = dcor)) + 
#   geom_point(aes(shape = experiment), color = plotdf$CCM_color) + 
#   xlab("Pearson correlation") + ylab("distance correlation") + 
#   ggtitle("comparison of two measures of module\n average expression correlation between species") +
#   scale_shape_discrete(limits = c("CC", "HAP4", "LowN", "LowPi"),
#                        labels = c("Urea Shock", "Diauxic Shift", "Low Nitrogen", "Low Phosphorus")) +
#   theme_classic() + xlim(c(-1, 1)) + ylim(c(0, 1)) + geom_abline(slope = 0.5, intercept = 0.5, color = "gold")
# dev.off()
# 
# plot(abs(moduledf10$avgCor[moduledf10$is_CCM]), moduledf10$dcor[moduledf10$is_CCM], type = "n",
#      xlab = "Average Expression Correlation", ylab = "Distance Correlation",
#      main = "Two measures of expression similarity,\n weighted evenly between all 4 experiments")
# text(abs(moduledf10$avgCor[moduledf10$is_CCM]), moduledf10$dcor[moduledf10$is_CCM], moduledf10$CCM_color[moduledf10$is_CCM])
# 
# # AvgCor vs Eigengene distance
# # Just CCMs
# plot(abs(moduledf10$avgCor[moduledf10$is_CCM]), moduledf10$eigenDist[moduledf10$is_CCM], type = "n",
#      xlab = "Average Expression Correlation", ylab = "Eigengene Distance",
#      main = "Two measures of expression similarity,\n weighted evenly between all 4 experiments")
# text(abs(moduledf10$avgCor[moduledf10$is_CCM]), moduledf10$eigenDist[moduledf10$is_CCM], moduledf10$CCM_color[moduledf10$is_CCM]) # mostly negatively correlated
# # all modules (of size 10 or greater)
# plot(abs(moduledf10$avgCor[moduledf10$block_size >= 10]), moduledf10$eigenDist[moduledf10$block_size >= 10], type = "n",
#      xlab = "Average Expression Correlation", ylab = "Eigengene Distance",
#      main = "Two measures of expression similarity,\n weighted evenly between all 4 experiments")
# text(abs(moduledf10$avgCor[moduledf10$block_size >= 10]), moduledf10$eigenDist[moduledf10$block_size >= 10], 
#      moduledf10$module_name[moduledf10$block_size >= 10]) # mostly negatively correlated
# # correlation vs area between curves
# # Just CCMs
# plot(abs(moduledf10$avgCor[moduledf10$is_CCM]), moduledf10$abc[moduledf10$is_CCM], type = "n",
#      xlab = "Average Expression Correlation", ylab = "Area Between Curves",
#      main = "Two measures of expression similarity,\n weighted evenly between all 4 experiments")
# text(abs(moduledf10$avgCor[moduledf10$is_CCM]), moduledf10$abc[moduledf10$is_CCM], moduledf10$CCM_color[moduledf10$is_CCM])
# # all modules (of size 10 or greater)
# plot(abs(moduledf10$avgCor[moduledf10$block_size >= 10]), moduledf10$abc[moduledf10$block_size >= 10], type = "n",
#      xlab = "Average Expression Correlation", ylab = "Area Between Curves",
#      main = "Two measures of expression similarity,\n weighted evenly between all 4 experiments")
# text(abs(moduledf10$avgCor[moduledf10$block_size >= 10]), moduledf10$abc[moduledf10$block_size >= 10], 
#      moduledf10$module_name[moduledf10$block_size >= 10])
# # eigengene dist vs area between curves
# plot(abs(moduledf10$eigenDist[moduledf10$is_CCM]), moduledf10$abc[moduledf10$is_CCM], type = "n",
#      xlab = "Eigengene Distance", ylab = "Area Between Curves",
#      main = "Two measures of expression similarity,\n weighted evenly between all 4 experiments")
# text(abs(moduledf10$eigenDist[moduledf10$is_CCM]), moduledf10$abc[moduledf10$is_CCM], moduledf10$CCM_color[moduledf10$is_CCM])
# # all modules (of size 10 or greater)
# plot(abs(moduledf10$eigenDist[moduledf10$block_size >= 10]), moduledf10$abc[moduledf10$block_size >= 10], type = "n",
#      xlab = "Eigengene Distance", ylab = "Area Between Curves",
#      main = "Two measures of expression similarity,\n weighted evenly between all 4 experiments")
# text(abs(moduledf10$eigenDist[moduledf10$block_size >= 10]), moduledf10$abc[moduledf10$block_size >= 10], 
#      moduledf10$module_name[moduledf10$block_size >= 10])
# # Avgcor and area between curves seem the least similar, 
# # makes sense b/c ABC takes into account level differences that
# # correlation would ignore
# # abc and eigendist are fairly related, although there's a group with high
# # eigendist and not that high abc
# 
# 

# #### Identifying true plasticity divergence using ABC ####
# 
# # Comparing plasticity between species in each environment in the dataset
# 
# # plot of all modules x each environment
# plotdf0 <- moduledf25 |> 
#   filter(block_size >= 10) |> 
#   select(module_name, coexpressed, 
#          avgCor_CC, avgCor_LowN,
#          avgCor_LowPi, avgCor_HAP4) |> 
#   pivot_longer(cols = c(avgCor_CC, avgCor_LowN,
#                         avgCor_LowPi, avgCor_HAP4),
#                names_to = "experiment", values_to = "avgCor", names_prefix = "avgCor_")
# plotdf1 <- moduledf10 |> 
#   filter(block_size >= 10) |> 
#   select(module_name, coexpressed, 
#          mean_CC, mean_LowN, mean_LowPi, mean_HAP4) |> 
#   pivot_longer(cols = c(mean_CC, mean_LowN, mean_LowPi, mean_HAP4),
#                names_to = "experiment", values_to = "mean", names_prefix = "mean_")
# plotdf <- left_join(plotdf0, plotdf1, by = c("module_name", "coexpressed", "experiment"))
# plotdf$module_number <- parse_number(plotdf$module_name)
# 
# ggplot(plotdf, aes(x = avgCor, y = log2(mean))) + 
#   geom_text(aes(label = module_number)) + 
#   geom_point(aes(shape = experiment,
#                  color = coexpressed)) +
#   scale_color_discrete(limits = c("conserved co-expressed",
#                                   "diverged co-expressed",
#                                   "S. cerevisiae co-expressed",
#                                   "S. paradoxus co-expressed",
#                                   "never co-expressed"),
#                        type = c("gold","mediumseagreen","orange1", "blue2", "grey")) +
#   theme_classic() +
#   geom_vline(xintercept = 0.5, color = "red")
# 
# # TODO: consider replacing avgCor vs block size plot in figure with something like the above
# # that better highlights which experiments are exposing the most module divergence and how this
# # divergence relates to expression level
# # (also consider using residual variance instead of mean?)
# # dcor vs abc
# plotdf <- moduledf10 |> filter(is_CCM) |> select(dcor_CC,
#                                                  dcor_LowN,
#                                                  dcor_HAP4,
#                                                  dcor_LowPi,
#                                                  module_name,
#                                                  CCM_color,
#                                                  block_size) |> 
#   pivot_longer(cols = c("dcor_CC", 
#                         "dcor_LowN",
#                         "dcor_HAP4",
#                         "dcor_LowPi"), names_to = "experiment",
#                values_to = "dcor", names_prefix = "dcor_") |> 
#   left_join(y = pivot_longer(select(moduledf10, 
#                                     abc_CC,
#                                     abc_LowN,
#                                     abc_HAP4,
#                                     abc_LowPi, 
#                                     module_name),
#                              cols = c("abc_CC",
#                                       "abc_LowN",
#                                       "abc_HAP4",
#                                       "abc_LowPi"),
#                              names_to = "experiment",
#                              values_to = "abc", 
#                              names_prefix = "abc_"),
#             by = c("module_name", "experiment"))
# 
# ggplot(plotdf, aes(x = abc, y = dcor)) + 
#   geom_point(aes(color = CCM_color,
#                  shape = experiment)) +
#   geom_text(aes(label = module_name,
#                 color = CCM_color),
#             nudge_x = 0.05) +
#   scale_color_discrete(type = plotdf$CCM_color,
#                        limits = plotdf$CCM_color) +
#   theme_classic() +
#   xlab("Area Between Curves") +
#   ylab("Module correlation")
#### Old distribution plot avgCor w/ ranked block size ####
# # line and density plots looking at distribution of
# # expression variability and divergence between the 4 experiments for each module
# # (and how module size doesn't predict divergence)
# library(ggExtra, ggpubr)
# plotModuleLineAndDensity <- function(.mdf, .rmdf) {
#   plotdf <- .mdf
#   random_plotdf <- .rmdf |> mutate(coexpressed = "random")
#   plotdf <- bind_rows(plotdf, random_plotdf) |> 
#     mutate(block_rank = rank(block_size, ties.method = "first")) |> 
#     pivot_longer(cols = c("avgCor_CC",
#                           "avgCor_LowN",
#                           "avgCor_LowPi",
#                           "avgCor_HAP4"),
#                  names_to = "experiment",
#                  values_to = "avg_expression")
#   plotdf <- plotdf |> filter(block_size >= 10)
#   # lineplot
#   p_line <- ggplot(plotdf, aes(y = block_rank, x = avg_expression)) + 
#     geom_line(aes(color = coexpressed, group = paste(module_name, coexpressed))) +
#     geom_point(aes(color = coexpressed)) +
#     #geom_text(aes(label = module_name), color = "black", size = 3, check_overlap = TRUE) + # <- add this line to see which modules are which
#     scale_color_discrete(limits = c("conserved co-expressed",
#                                     "diverged co-expressed",
#                                     "S. cerevisiae co-expressed",
#                                     "S. paradoxus co-expressed",
#                                     "never co-expressed",
#                                     "random"),
#                          type = c("gold","mediumseagreen","orange1", "blue2", "grey", "purple1")) +
#     ylab("Module size") +
#     xlab("Module similarity between species\n (average expression correlation)") +
#     theme_classic() +
#     theme(axis.ticks.y = element_blank(),
#           axis.text.y = element_blank())
#   
#   # corresponding density plot
#   # p_density <- ggplot(plotdf, aes(x = avg_expression)) + 
#   #   geom_density(aes(fill = coexpressed), alpha = 0.7) +
#   #   theme_classic() +
#   #   theme(legend.title = element_blank()) +
#   #   xlab("") +
#   #   scale_fill_discrete(limits = c("conserved co-expressed",
#   #                                    "diverged co-expressed",
#   #                                    "S. cerevisiae co-expressed",
#   #                                    "S. paradoxus co-expressed",
#   #                                    "never co-expressed",
#   #                                    "random"),
#   #                       type = c("gold","mediumseagreen","orange1", "blue2", "grey", "purple1")) +
#   #   theme(axis.ticks.y = element_blank(),
#   #         axis.text.y = element_blank(),
#   #         axis.ticks.x = element_blank(),
#   #         axis.text.x = element_blank())
#   return(ggMarginal(p_line, groupFill = TRUE, groupColour = TRUE, margins = "x"))
# }
# 
# # merge 10
# p_10 <- plotModuleLineAndDensity(filter(moduledf10, is_CCM), random_moduledf10)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/EnvironmentalPatterns/avgCor10.pdf",
#     width = 5, height = 10)
# p_10
# dev.off()
# 
# # supplementary line/density plots at other merge levels
# # merge 00
# p_00 <- plotModuleLineAndDensity(moduledf00, random_moduledf00)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/avgCor00.pdf",
#     width = 5, height = 10)
# p_00
# dev.off()
# 
# # merge 25
# p_25 <- plotModuleLineAndDensity(moduledf25, random_moduledf25)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/avgCor25.pdf",
#     width = 5, height = 10)
# p_25
# dev.off()
# 
# # merge 35
# p_35 <- plotModuleLineAndDensity(moduledf35, random_moduledf35)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/avgCor35.pdf",
#     width = 5, height = 10)
# p_35
# dev.off()
# 
# # # Supplemental lineplot of avg cor variation by experiment
# # so you know which experiments exposed the most amount of divergence
# # # TODO: it's a little too messy, might need to separate the 3 classes into 3 plots with the same y axis
# # pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/AvgCor_lineplots.pdf",
# #     width = 7, height = 3)
# # ggplot(data = filter(plotdf, cer_color != "grey" & par_color != "grey")) + 
# #   geom_line(data = filter(plotdf, module_name %in% c("h Conserved", "f Conserved")),
# #               aes(x = experiment, y = avg_expression,
# #                   group = module_name),
# #               color = "blue", alpha = 0.5, linewidth = 3) +
# #   geom_line(data = filter(plotdf, module_name %in% c("18 Divergent", "8 Divergent")),
# #             aes(x = experiment, y = avg_expression,
# #                 group = module_name),
# #             color = "red", alpha = 0.5, linewidth = 3) +
# #   geom_line(aes(y = avg_expression, x = experiment, group = module_name, color = type)) +
# #   geom_text(aes(y = avg_expression, x = experiment, label = module_name), check_overlap = TRUE) + # add this line for checking which module is which
# #   theme_classic() +
# #   xlab("Experiment") +
# #   scale_color_discrete(labels = c("Random", "Divergent", "Conserved"),
# #                        type = c("mediumpurple1", "grey80", "gold")) +
# #   scale_x_discrete(labels = c("Cell Cycle", "Growth Curve", "Low Nitrogen", "Low Phosphorus")) +
# #   ylab("Correlation in module expression \npattern between species") +
# #   theme(legend.title = element_blank())
# # dev.off()
# 
# #### Old Supplemental Figure: All expression profiles ####
# 
# # merge 10
# plotdf <- moduledf10 |> filter(block_size >= 10)
# plotlist_ccm <- vector(mode = "list", length = sum(plotdf$coexpressed == "conserved co-expressed"))
# plotlist_div <- vector(mode = "list", length = sum(plotdf$coexpressed == "diverged co-expressed"))
# plotlist_cer <- vector(mode = "list", length = sum(plotdf$coexpressed == "S. cerevisiae co-expressed"))
# plotlist_par <- vector(mode = "list", length = sum(plotdf$coexpressed == "S. paradoxus co-expressed"))
# plotlist_nev <- vector(mode = "list", length = sum(plotdf$coexpressed == "never co-expressed"))
# for (idx in 1:nrow(plotdf)) {
#   m <- plotdf$module_name[idx]
#   gene_idxs <- module_genedf10 |> filter(module_name == m) |> select(gene_name) |> pull()
#   if (plotdf$coexpressed[idx] == "conserved co-expressed") {
#     p <- plotExpressionProfilePairStack(counts_all2$cer[,gene_idxs],
#                                         counts_all2$par[,gene_idxs],
#                                         info,
#                                         info,
#                                         .method = "line", .show_points = TRUE,
#                                         .normalization = "scale")
#     p <- annotate_figure(p, top = m)
#     p_idx <- which(plotdf$module_name[plotdf$coexpressed == "conserved co-expressed"] == m)
#     plotlist_ccm[[p_idx]] <- p
#   }
#   if (plotdf$coexpressed[idx] == "diverged co-expressed") {
#     p <- plotExpressionProfilePairStack(counts_all2$cer[,gene_idxs],
#                                         counts_all2$par[,gene_idxs],
#                                         info,
#                                         info,
#                                         .method = "line", .show_points = TRUE,
#                                         .normalization = "scale")
#     p <- annotate_figure(p, top = m)
#     p_idx <- which(plotdf$module_name[plotdf$coexpressed == "diverged co-expressed"] == m)
#     plotlist_div[[p_idx]] <- p
#   }
#   if (plotdf$coexpressed[idx] == "S. cerevisiae co-expressed") {
#     p <- plotExpressionProfilePairStack(counts_all2$cer[,gene_idxs],
#                                         counts_all2$par[,gene_idxs],
#                                         info,
#                                         info,
#                                         .method = "line", .show_points = TRUE,
#                                         .normalization = "scale")
#     p <- annotate_figure(p, top = m)
#     p_idx <- which(plotdf$module_name[plotdf$coexpressed == "S. cerevisiae co-expressed"] == m)
#     plotlist_cer[[p_idx]] <- p
#   }
#   if (plotdf$coexpressed[idx] == "S. paradoxus co-expressed") {
#     p <- plotExpressionProfilePairStack(counts_all2$cer[,gene_idxs],
#                                         counts_all2$par[,gene_idxs],
#                                         info,
#                                         info,
#                                         .method = "line", .show_points = TRUE,
#                                         .normalization = "scale")
#     p <- annotate_figure(p, top = m)
#     p_idx <- which(plotdf$module_name[plotdf$coexpressed == "S. paradoxus co-expressed"] == m)
#     plotlist_par[[p_idx]] <- p
#   }
#   if (plotdf$coexpressed[idx] == "never co-expressed") {
#     p <- plotExpressionProfilePairStack(counts_all2$cer[,gene_idxs],
#                                         counts_all2$par[,gene_idxs],
#                                         info,
#                                         info,
#                                         .method = "line", .show_points = TRUE,
#                                         .normalization = "scale")
#     p <- annotate_figure(p, top = m)
#     p_idx <- which(plotdf$module_name[plotdf$coexpressed == "never co-expressed"] == m)
#     plotlist_nev[[p_idx]] <- p
#   }
# }
# # random
# plotdf <- random_moduledf10 |> filter(block_size >= 10)
# plotlist_Random <- vector(mode = "list", length = nrow(plotdf))
# for (idx in 1:nrow(plotdf)) {
#   m <- plotdf$module_name[idx]
#   gene_idxs <- random_module_genedf10 |> filter(module_name == m) |> select(gene_name) |> pull()
#   p <- plotExpressionProfilePairStack(counts_all2$cer[,gene_idxs],
#                                       counts_all2$par[,gene_idxs],
#                                       info,
#                                       info,
#                                       .method = "line", .show_points = TRUE,
#                                       .normalization = "scale")
#   p <- annotate_figure(p, top = m)
#   plotlist_Random[[idx]] <- p
# }
# # making pdfs
# # CCMs
# p_CCM <- ggarrange(plotlist = plotlist_ccm, 
#                    nrow = ceiling(length(plotlist_ccm)/3), 
#                    ncol = 3, 
#                    common.legend = TRUE, 
#                    legend = "right")
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_expr_profiles_CCM10.pdf",
#     width = 12, height = length(plotlist_ccm)*4)
# annotate_figure(p_CCM, bottom = textGrob("Timepoint", gp = gpar(cex = 3)))
# dev.off()
# # Divergent
# p_Divergent <- ggarrange(plotlist = plotlist_div, 
#                          nrow = ceiling(length(plotlist_div)/3), 
#                          ncol = 3, 
#                          common.legend = TRUE, 
#                          legend = "right")
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_expr_profiles_Divergent10.pdf",
#     width = 12, height = length(plotlist_div)*4)
# annotate_figure(p_Divergent, bottom = textGrob("Timepoint", gp = gpar(cex = 3)))
# dev.off()
# # S. cerevisiae co-expressed
# p_Scer <- ggarrange(plotlist = plotlist_cer, 
#                     nrow = ceiling(length(plotlist_cer)/3), 
#                     ncol = 3, 
#                     common.legend = TRUE, 
#                     legend = "right")
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_expr_profiles_Scer10.pdf",
#     width = 12, height = length(plotlist_cer)*4)
# annotate_figure(p_Scer, bottom = textGrob("Timepoint", gp = gpar(cex = 3)))
# dev.off()
# # S. paradoxus co-expressed
# p_Spar <- ggarrange(plotlist = plotlist_par, 
#                     nrow = ceiling(length(plotlist_par)/3), 
#                     ncol = 3, 
#                     common.legend = TRUE, 
#                     legend = "right")
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_expr_profiles_Spar10.pdf",
#     width = 12, height = length(plotlist_par)*4)
# annotate_figure(p_Spar, bottom = textGrob("Timepoint", gp = gpar(cex = 3)))
# dev.off()
# # Never co-expressed
# p_Never <- ggarrange(plotlist = plotlist_nev, 
#                      nrow = ceiling(length(plotlist_nev)), 
#                      ncol = 1, 
#                      common.legend = TRUE, 
#                      legend = "right")
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_expr_profiles_Never10.pdf",
#     width = 4, height = 16)
# annotate_figure(p_Never, bottom = textGrob("Timepoint", gp = gpar(cex = 3)))
# dev.off()
# # Random
# p_Random <- ggarrange(plotlist = plotlist_Random, 
#                       nrow = ceiling(length(plotlist_Random)/3), 
#                       ncol = 3, 
#                       common.legend = TRUE, 
#                       legend = "right")
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_expr_profiles_Random10.pdf",
#     width = 12, height = length(plotlist_Random)*4)
# annotate_figure(p_Random,
#                 bottom = textGrob("Timepoint", gp = gpar(cex = 3)))
# dev.off()
# 

#### Hybrid improvement exploration ####
# ultimately led me back to separating into 5 co-expression groups
# and acknowledging that there are differences in environmental variability 
# in each group that in turn affect differences in correlation
## just Conserved
# p_improve_conserved <- ggplot(filter(divergentdf, type == "Conserved"), aes(x = parental_cor, y = hybrid_cor)) + 
#   geom_point(aes(color = type)) +
#   scale_color_discrete(type = c("gold", "mediumseagreen", "orange1", "blue2", "grey"),
#                        limits = c("conserved co-expressed",
#                                   "diverged co-expressed",
#                                   "S. cerevisiae co-expressed",
#                                   "S. paradoxus co-expressed",
#                                   "never co-expressed")) +
#   theme_classic() +
#   xlim(c(-0.5, 1)) +
#   ylim(c(-0.5, 1)) +
#   geom_abline(slope = 1, intercept = 0)
# ggMarginal(p_improve_conserved, groupColour = TRUE, groupFill = TRUE)
# 
# p_improve_diverged <- ggplot(filter(improvedf, type == "Divergent"), aes(x = parental_cor, y = hybrid_cor)) + 
#   geom_point(aes(color = type)) +
#   scale_color_discrete(type = c("gold", "grey80"),
#                        limits = c("Conserved", "Divergent")) +
#   theme_classic() +
#   xlim(c(-0.5, 1)) +
#   ylim(c(-0.5, 1)) +
#   geom_abline(slope = 1, intercept = 0)
# ggMarginal(p_improve_diverged, groupColour = TRUE, groupFill = TRUE)
# # there's a real glut of genes that are uncorrelated in both parents and hybrids
# improvedf$const_uncor <- improvedf$parental_cor < 0.5 & improvedf$hybrid_cor < 0.5
# improvedf |> select(type, const_uncor) |> table() # 6x as many constitutively uncorrelated are in divergent versus CCM, despite the fact that only 3x as many genes are in divergent versus CCM
# # what modules are they in?
# const_uncor_idxs <- improvedf |> filter(const_uncor) |> select(gene_name) |> pull()
# modcoltab <- module_genedf25[, c("cer_color", "par_color")] |> table()
# constuncor_modcoltab <- module_genedf25[module_genedf25$gene_name %in% const_uncor_idxs, c("cer_color", "par_color")] |> table()
# constuncor_modcoltab # counts
# round(constuncor_modcoltab/modcoltab, 2) # fractions of modules
# # bingo-bango our greys are well represented
# 
# # what proportion of const_uncor genes are grey in either or both species?
# improvedf <- left_join(improvedf, select(module_genedf25, gene_name, cer_color, par_color),
#                        by = "gene_name")
# improvedf$evergrey <- improvedf$cer_color == "grey" | improvedf$par_color == "grey"
# p_evergrey <- ggplot(improvedf, aes(x = parental_cor, y = hybrid_cor)) + 
#   geom_point(aes(color = evergrey)) +
#   theme_classic() +
#   xlim(c(-0.5, 1)) +
#   ylim(c(-0.5, 1)) +
#   geom_abline(slope = 1, intercept = 0) +
#   scale_color_discrete(limits = c(FALSE, TRUE),
#                        labels = c("co-expressed", "not co-expressed"),
#                        type = c("gold", "grey80"),
#                        name = "") +
#   xlab("parental species") +
#   ylab("hybrid alleles") +
#   ggtitle("Expression pattern correlation\nbetween S. cerevisiae and S. paradoxus")
# pdf("../../aligning_the_molecular_phenotype/paper_figures/PatternVslevel/evergrey.pdf",
#     width = 7, height = 5)
# ggMarginal(p_evergrey, groupColour = TRUE, groupFill = TRUE) # a really high proportion as it turns out
# dev.off()
# improvedf |> select(evergrey, const_uncor) |> table()
# 
# 
# # TODO: show how this same pattern isn't true for effect size
# # scale is going to be tricky, "improve" means decrease in magnitude,
# # so maybe plot abs value? 
# # But there will still be a couple that are like lfc of 6 who will
# # throw off the axis limits
# p_evergreyLFC <- ggplot(improvedf, aes(x = abs(parental_lfc), y = abs(hybrid_lfc))) + 
#   geom_point(aes(color = evergrey)) +
#   geom_text(data = filter(improvedf,
#                           abs(parental_lfc) > 3 |
#                             abs(hybrid_lfc) > 3),
#             aes(x = abs(parental_lfc), y = abs(hybrid_lfc),
#                 label = gene_name), check_overlap = TRUE,
#             size = 2.5) +
#   theme_classic() +
#   xlim(c(0, max_lfc)) +
#   ylim(c(0, max_lfc)) +
#   geom_abline(slope = 1, intercept = 0) +
#   scale_color_discrete(limits = c(FALSE, TRUE),
#                        labels = c("co-expressed", "not co-expressed"),
#                        type = c("gold", "grey80"),
#                        name = "") +
#   xlab("parental species") +
#   ylab("hybrid alleles") +
#   ggtitle("Expression level change (LFC)\nbetween S. cerevisiae and S. paradoxus")
# pdf("../../aligning_the_molecular_phenotype/paper_figures/PatternVslevel/evergrey_lfc.pdf",
#     width = 7, height = 5)
# ggMarginal(p_evergreyLFC, groupColour = TRUE, groupFill = TRUE)
# dev.off()
# 
# # TODO: illustrate how evergrey genes are much less likely to improve
# # their correlation in the hybrid than co-expressed genes with line plots
# # of all the improving, worsening, and staying-the-same genes separated
# # for both groups (6 line plots total)
# improvedf$cor_improvement <- improvedf$hybrid_cor - improvedf$parental_cor
# # co-expressed, positive cor improvement
# plotdf <- filter(improvedf, cor_improvement > 0.1 & !evergrey) |> 
#   pivot_longer(cols = c("parental_cor", "hybrid_cor"),
#                names_to = "parent_or_hybrid",
#                values_to = "cor") |> 
#   mutate(parent_or_hybrid = gsub("_cor", "", parent_or_hybrid))
# ngenes <- plotdf$gene_name |> unique() |> length()
# p_coex_improve <- ggplot(plotdf, aes(x = parent_or_hybrid, y = cor)) + 
#   geom_line(aes(group = gene_name), color = "gold") +
#   theme_classic() +
#   xlab("") + 
#   ylab("correlation") +
#   annotate("text", x = 1.5, y = 1.3, label = paste(ngenes, "genes")) +
#   scale_x_discrete(breaks = c("parental", "hybrid"),
#                    limits = c("parental", "hybrid")) +
#   ggtitle("Co-expressed genes with \ncorrelation improvement in hybrid")
# # co-expressed, no change in cor
# plotdf <- filter(improvedf, cor_improvement < 0.1 & 
#                    cor_improvement > -0.1 & !evergrey) |> 
#   pivot_longer(cols = c("parental_cor", "hybrid_cor"),
#                names_to = "parent_or_hybrid",
#                values_to = "cor") |> 
#   mutate(parent_or_hybrid = gsub("_cor", "", parent_or_hybrid))
# ngenes <- plotdf$gene_name |> unique() |> length()
# p_coex_nochange <- ggplot(plotdf, aes(x = parent_or_hybrid, y = cor)) + 
#   geom_line(aes(group = gene_name), color = "gold") +
#   theme_classic() +
#   xlab("") + 
#   ylab("correlation") +
#   annotate("text", x = 1.5, y = 1.3, label = paste(ngenes, "genes")) +
#   scale_x_discrete(breaks = c("parental", "hybrid"),
#                    limits = c("parental", "hybrid")) +
#   ggtitle("Co-expressed genes with \nlittle change in correlation in hybrid")
# # co-expressed, cor worsens in hybrid
# plotdf <- filter(improvedf, cor_improvement < -0.1 & !evergrey) |> 
#   pivot_longer(cols = c("parental_cor", "hybrid_cor"),
#                names_to = "parent_or_hybrid",
#                values_to = "cor") |> 
#   mutate(parent_or_hybrid = gsub("_cor", "", parent_or_hybrid))
# ngenes <- plotdf$gene_name |> unique() |> length()
# p_coex_worse <- ggplot(plotdf, aes(x = parent_or_hybrid, y = cor)) + 
#   geom_line(aes(group = gene_name), color = "gold") +
#   theme_classic() +
#   xlab("") + 
#   ylab("correlation") +
#   annotate("text", x = 1.5, y = 1.3, label = paste(ngenes, "genes")) +
#   scale_x_discrete(breaks = c("parental", "hybrid"),
#                    limits = c("parental", "hybrid")) +
#   ggtitle("Co-expressed genes with \nworse correlation in hybrid")
# # evergrey, positive cor improvement
# plotdf <- filter(improvedf, cor_improvement > 0.1 & evergrey) |> 
#   pivot_longer(cols = c("parental_cor", "hybrid_cor"),
#                names_to = "parent_or_hybrid",
#                values_to = "cor") |> 
#   mutate(parent_or_hybrid = gsub("_cor", "", parent_or_hybrid))
# ngenes <- plotdf$gene_name |> unique() |> length()
# p_evergrey_improve <- ggplot(plotdf, aes(x = parent_or_hybrid, y = cor)) + 
#   geom_line(aes(group = gene_name), color = "grey80") +
#   theme_classic() +
#   xlab("") + 
#   ylab("correlation") +
#   annotate("text", x = 1.5, y = 1.3, label = paste(ngenes, "genes")) +
#   scale_x_discrete(breaks = c("parental", "hybrid"),
#                    limits = c("parental", "hybrid")) +
#   ggtitle("Not co-expressed genes with \ncorrelation improvement in hybrid")
# # evergrey, no change in cor
# plotdf <- filter(improvedf, cor_improvement < 0.1 & 
#                    cor_improvement > -0.1 & evergrey) |> 
#   pivot_longer(cols = c("parental_cor", "hybrid_cor"),
#                names_to = "parent_or_hybrid",
#                values_to = "cor") |> 
#   mutate(parent_or_hybrid = gsub("_cor", "", parent_or_hybrid))
# ngenes <- plotdf$gene_name |> unique() |> length()
# p_evergrey_nochange <- ggplot(plotdf, aes(x = parent_or_hybrid, y = cor)) + 
#   geom_line(aes(group = gene_name), color = "grey80") +
#   theme_classic() +
#   xlab("") + 
#   ylab("correlation") +
#   annotate("text", x = 1.5, y = 1.3, label = paste(ngenes, "genes")) +
#   scale_x_discrete(breaks = c("parental", "hybrid"),
#                    limits = c("parental", "hybrid")) +
#   ggtitle("Not co-expressed genes with \nlittle change in correlation in hybrid")
# # evergrey, cor worsens in hybrid
# plotdf <- filter(improvedf, cor_improvement < -0.1 & evergrey) |> 
#   pivot_longer(cols = c("parental_cor", "hybrid_cor"),
#                names_to = "parent_or_hybrid",
#                values_to = "cor") |> 
#   mutate(parent_or_hybrid = gsub("_cor", "", parent_or_hybrid))
# ngenes <- plotdf$gene_name |> unique() |> length()
# p_evergrey_worse <- ggplot(plotdf, aes(x = parent_or_hybrid, y = cor)) + 
#   geom_line(aes(group = gene_name), color = "grey80") +
#   theme_classic() +
#   xlab("") + 
#   ylab("correlation") +
#   annotate("text", x = 1.5, y = 1.3, label = paste(ngenes, "genes")) +
#   scale_x_discrete(breaks = c("parental", "hybrid"),
#                    limits = c("parental", "hybrid")) +
#   ggtitle("Not co-expressed genes with \nworse correlation in hybrid")
# 
# ggarrange(p_coex_improve, p_coex_nochange, p_coex_worse,
#           p_evergrey_improve, p_evergrey_nochange, p_evergrey_worse,
#           nrow = 2, ncol = 3) # yeahhhh there's too many to see what's going on
# # but based on counts we can conclude that there's not really an enrichment for 
# # a certain cor change direction in either group---seems to be more that
# # the subset of genes with the highest hybrid correlation is all co-expressed,
# # and the vast majority of evergrey (88%) have really low hybrid correlation 
# # (strong cis effects on their expression pattern)
# improvedf |> mutate("hybrid_diff" = hybrid_cor < 0.5) |> select(evergrey, hybrid_diff) |> table()
# 
# # plotExpressionProfilePair one random evergrey gene in parent and hybrid
# # and one random co-expressed gene in parent and hybrid as a companion to the above point plots
# random_evergrey_idx <- improvedf |> filter(evergrey) |> select(gene_name) |> pull() |> sample(1)
# improvedf |> filter(gene_name == random_evergrey_idx) |> t()
# plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,random_evergrey_idx, drop = FALSE],
#                              .cts2 = counts_all2$par[,random_evergrey_idx, drop = FALSE],
#                              .cts3 = counts_all2_allele$cer[,random_evergrey_idx, drop = FALSE],
#                              .cts4 = counts_all2_allele$par[,random_evergrey_idx, drop = FALSE],
#                              .info1 = info,
#                              .info2 = info,
#                              .info3 = info,
#                              .info4 = info,
#                              .name1 = "S. cerevisiae parent",
#                              .name2 = "S. cerevisiae hybrid",
#                              .name3 = "S. paradoxus parent",
#                              .name4 = "S. paradoxus hybrid",
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .show_points = FALSE, .method = "line",
#                              .normalization = "none")
# # random co-expressed idx
# random_coex_idx <- improvedf |> filter(!evergrey) |> select(gene_name) |> pull() |> sample(1)
# improvedf |> filter(gene_name == random_coex_idx) |> t()
# plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,random_coex_idx, drop = FALSE],
#                              .cts2 = counts_all2$par[,random_coex_idx, drop = FALSE],
#                              .cts3 = counts_all2_allele$cer[,random_coex_idx, drop = FALSE],
#                              .cts4 = counts_all2_allele$par[,random_coex_idx, drop = FALSE],
#                              .info1 = info,
#                              .info2 = info,
#                              .info3 = info,
#                              .info4 = info,
#                              .name1 = "S. cerevisiae parent",
#                              .name2 = "S. cerevisiae hybrid",
#                              .name3 = "S. paradoxus parent",
#                              .name4 = "S. paradoxus hybrid",
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .show_points = FALSE, .method = "line",
#                              .normalization = "none")
# 
# # are evergrey genes more noisy? Or lower expressed?
# # mean and variance of each gene's counts in both parents combined
# improvedf$mean_expr_parents <- colMeans(rbind(counts_all2$cer, counts_all2$par))
# improvedf$var_parents <- apply(rbind(counts_all2$cer, counts_all2$par), 2, var)
# p_meanvar_parents <- ggplot(improvedf, aes(x = log2(mean_expr_parents), y = log2(var_parents))) +
#   geom_point(aes(color = evergrey)) +
#   scale_color_discrete(limits = c(TRUE, FALSE),
#                        type = c("grey80", "gold")) +
#   theme_classic()
# ggMarginal(p_meanvar_parents, groupColour = TRUE, groupFill = TRUE)
# # mean and var in hybrids, both alleles
# improvedf$mean_expr_hyb <- colMeans(rbind(counts_all2_allele$cer, counts_all2_allele$par))
# improvedf$var_hyb <- apply(rbind(counts_all2_allele$cer, counts_all2_allele$par), 2, var)
# p_meanvar_hyb <- ggplot(improvedf, aes(x = log2(mean_expr_hyb), y = log2(var_hyb))) +
#   geom_point(aes(color = evergrey)) +
#   scale_color_discrete(limits = c(TRUE, FALSE),
#                        type = c("grey80", "gold")) +
#   theme_classic()
# ggMarginal(p_meanvar_hyb, groupColour = TRUE, groupFill = TRUE)
# # not really more noisy, both groups have the same mean-var relationship,
# # but a little lower expressed on average
# 
# # Are the genes with the higher correlation the ones that are higher expressed?
# p_mean_cor_parents <- ggplot(improvedf, aes(x = log2(mean_expr_parents), y = parental_cor)) +
#   geom_point(aes(color = evergrey)) +
#   scale_color_discrete(limits = c(TRUE, FALSE),
#                        type = c("grey80", "gold"))
# ggMarginal(p_mean_cor_parents, groupColour = TRUE, groupFill = TRUE) # not really?
# p_var_cor_parents <- ggplot(improvedf, aes(x = log2(var_parents/mean_expr_parents), y = parental_cor)) +
#   geom_point(aes(color = evergrey)) +
#   scale_color_discrete(limits = c(TRUE, FALSE),
#                        type = c("grey80", "gold"))
# ggMarginal(p_var_cor_parents, groupColour = TRUE, groupFill = TRUE) 
# # slightly more of a relationship, but in the opposite way that I was
# # expecting. The genes with higher scaled variance also have higher correlation
# # TODO: look at Tau for mean-independent variability and 
# # variance among LowN 0h, YPD WT samples to gauge noise
# 
# # Tau measure of expression variability independent of mean
# # rows are genes, cols are samples
# calculateTau <- function(cts) {
#   m <- ncol(cts)
#   max_per_gene <- apply(cts, 1, max)
#   unscaled_tau <- (1 - cts/max_per_gene) %>% rowSums(na.rm = TRUE)
#   return(unscaled_tau/(m - 1))
# }
# improvedf$tau_parents <- calculateTau(t(rbind(counts_all2$cer, counts_all2$par)))
# p_tau_cor_parents <- ggplot(improvedf, aes(x = tau_parents, y = parental_cor)) +
#   geom_point(aes(color = evergrey)) +
#   scale_color_discrete(limits = c(TRUE, FALSE),
#                        type = c("grey80", "gold"))
# ggMarginal(p_tau_cor_parents, groupColour = TRUE, groupFill = TRUE) 
# # The highest cor genes have the highest tau, but plenty have high cor and low tau
# 
# # What about the genes with good hybrid correlation?
# p_mean_cor_hyb <- ggplot(improvedf, aes(x = log2(mean_expr_hyb), y = hybrid_cor)) +
#   geom_point(aes(color = evergrey)) +
#   scale_color_discrete(limits = c(TRUE, FALSE),
#                        type = c("grey80", "gold"))
# ggMarginal(p_mean_cor_hyb, groupColour = TRUE, groupFill = TRUE) # not really?
# p_var_cor_hyb <- ggplot(improvedf, aes(x = log2(var_hyb/mean_expr_hyb), y = hybrid_cor)) +
#   geom_point(aes(color = evergrey)) +
#   scale_color_discrete(limits = c(TRUE, FALSE),
#                        type = c("grey80", "gold"))
# ggMarginal(p_var_cor_hyb, groupColour = TRUE, groupFill = TRUE) 
# # slightly more of a relationship, but in the opposite way that I was
# # expecting. The genes with higher scaled variance also have higher correlation
# 
# # Noise, var/mean in WT lowN 0h, YPD - have to go to TFdel to get replicate info
# improvedf$noise_parents <- apply(rbind(counts_TFdel$cer[infos_TFdel$cer$time_point_str == "0 h, YPD" &
#                                                           infos_TFdel$cer$genotype == "WT",],
#                                        counts_TFdel$par[infos_TFdel$par$time_point_str == "0 h, YPD" &
#                                                           infos_TFdel$par$genotype == "WT",]),
#                                  2,
#                                  \(x) {
#                                    var_x <- var(x)
#                                    mean_x <- mean(x)
#                                    return(var_x/mean_x)
#                                  })
# 
# p_noise_cor_parents <- ggplot(improvedf, aes(x = log2(noise_parents), y = parental_cor)) +
#   geom_point(aes(color = evergrey)) +
#   scale_color_discrete(limits = c(TRUE, FALSE),
#                        type = c("grey80", "gold"))
# ggMarginal(p_noise_cor_parents, groupColour = TRUE, groupFill = TRUE) 
# # Nope noise doesn't really seem to explain it. Variability across environments
# # is more of an explainer, but Tau can't be trusted
# 
# # Variability across environments, var/mean
# improvedf$scaled_var_parents <- improvedf$var_parents/improvedf$mean_expr_parents
# improvedf$scaled_var_hyb <- improvedf$var_hyb/improvedf$mean_expr_hyb
# # high scaled var, coexpressed
# high_var_coex_parents <- improvedf |> filter(!evergrey & scaled_var_parents > 2^9) |> 
#   select(gene_name) |> pull() |> sample(1)
# plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,high_var_coex_parents, drop = FALSE],
#                              .cts2 = counts_all2$par[,high_var_coex_parents, drop = FALSE],
#                              .cts3 = counts_all2_allele$cer[,high_var_coex_parents, drop = FALSE],
#                              .cts4 = counts_all2_allele$par[,high_var_coex_parents, drop = FALSE],
#                              .info1 = info,
#                              .info2 = info,
#                              .info3 = info,
#                              .info4 = info,
#                              .name1 = "S. cerevisiae parent",
#                              .name2 = "S. cerevisiae hybrid",
#                              .name3 = "S. paradoxus parent",
#                              .name4 = "S. paradoxus hybrid",
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .show_points = FALSE, .method = "loess",
#                              .normalization = "log2")
# # high var, evergrey
# high_var_evergrey_parents <- improvedf |> filter(evergrey & scaled_var_parents > 2^9) |> 
#   select(gene_name) |> pull() |> sample(1)
# plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,high_var_evergrey_parents, drop = FALSE],
#                              .cts2 = counts_all2$par[,high_var_evergrey_parents, drop = FALSE],
#                              .cts3 = counts_all2_allele$cer[,high_var_evergrey_parents, drop = FALSE],
#                              .cts4 = counts_all2_allele$par[,high_var_evergrey_parents, drop = FALSE],
#                              .info1 = info,
#                              .info2 = info,
#                              .info3 = info,
#                              .info4 = info,
#                              .name1 = "S. cerevisiae parent",
#                              .name2 = "S. cerevisiae hybrid",
#                              .name3 = "S. paradoxus parent",
#                              .name4 = "S. paradoxus hybrid",
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .show_points = FALSE, .method = "loess",
#                              .normalization = "log2")
# # low scaled var, coexpressed
# low_var_coex_parents <- improvedf |> filter(!evergrey & scaled_var_parents < 2^6) |> 
#   select(gene_name) |> pull() |> sample(1)
# plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,low_var_coex_parents, drop = FALSE],
#                              .cts2 = counts_all2$par[,low_var_coex_parents, drop = FALSE],
#                              .cts3 = counts_all2_allele$cer[,low_var_coex_parents, drop = FALSE],
#                              .cts4 = counts_all2_allele$par[,low_var_coex_parents, drop = FALSE],
#                              .info1 = info,
#                              .info2 = info,
#                              .info3 = info,
#                              .info4 = info,
#                              .name1 = "S. cerevisiae parent",
#                              .name2 = "S. cerevisiae hybrid",
#                              .name3 = "S. paradoxus parent",
#                              .name4 = "S. paradoxus hybrid",
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .show_points = FALSE, .method = "loess",
#                              .normalization = "log2")
# # low var, evergrey
# low_var_evergrey_parents <- improvedf |> filter(evergrey & scaled_var_parents < 2^6) |> 
#   select(gene_name) |> pull() |> sample(1)
# plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,low_var_evergrey_parents, drop = FALSE],
#                              .cts2 = counts_all2$par[,low_var_evergrey_parents, drop = FALSE],
#                              .cts3 = counts_all2_allele$cer[,low_var_evergrey_parents, drop = FALSE],
#                              .cts4 = counts_all2_allele$par[,low_var_evergrey_parents, drop = FALSE],
#                              .info1 = info,
#                              .info2 = info,
#                              .info3 = info,
#                              .info4 = info,
#                              .name1 = "S. cerevisiae parent",
#                              .name2 = "S. cerevisiae hybrid",
#                              .name3 = "S. paradoxus parent",
#                              .name4 = "S. paradoxus hybrid",
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .show_points = FALSE, .method = "loess",
#                              .normalization = "log2")
# # which genes improved their correlation the most in the hybrid?
# improvedf |> arrange(desc(cor_improvement)) |> select(gene_name, parental_cor, hybrid_cor)
# plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,"YAL003W", drop = FALSE],
#                              .cts2 = counts_all2$par[,"YAL003W", drop = FALSE],
#                              .cts3 = counts_all2_allele$cer[,"YAL003W", drop = FALSE],
#                              .cts4 = counts_all2_allele$par[,"YAL003W", drop = FALSE],
#                              .info1 = info,
#                              .info2 = info,
#                              .info3 = info,
#                              .info4 = info,
#                              .name1 = "S. cerevisiae parent",
#                              .name2 = "S. cerevisiae hybrid",
#                              .name3 = "S. paradoxus parent",
#                              .name4 = "S. paradoxus hybrid",
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .show_points = FALSE, .method = "line",
#                              .normalization = "none")
# # well I mean, that looks pretty correlated to begin with
# # TODO: how is the parental cor for YAL003W so low? It seems super correlated except for CC
# test <- getExpressionSimilarity(counts_all2$cer[,"YAL003W"],
#                                 counts_all2$par[,"YAL003W"])
# cor_overall <- cor(test$avg1, test$avg2)
# test <- getExpressionSimilarity(counts_all2$cer[info$experiment == "HAP4","YAL003W"],
#                                 counts_all2$par[info$experiment == "HAP4","YAL003W"])
# cor_HAP4 <- cor(test$avg1, test$avg2)
# test <- getExpressionSimilarity(counts_all2$cer[info$experiment == "CC","YAL003W"],
#                                 counts_all2$par[info$experiment == "CC","YAL003W"])
# cor_CC <- cor(test$avg1, test$avg2)
# test <- getExpressionSimilarity(counts_all2$cer[info$experiment == "LowN","YAL003W"],
#                                 counts_all2$par[info$experiment == "LowN","YAL003W"])
# cor_LowN <- cor(test$avg1, test$avg2)
# test <- getExpressionSimilarity(counts_all2$cer[info$experiment == "LowPi","YAL003W"],
#                                 counts_all2$par[info$experiment == "LowPi","YAL003W"])
# cor_LowPi <- cor(test$avg1, test$avg2)
# mean(c(cor_HAP4, cor_CC, cor_LowN, cor_LowPi))
# cor_overall # lettttts separate by experiment k thx
