sapply(c("plyr", "dplyr", "purrr", "tidyr", "ggpubr", "readr", "igraph", "Matrix", "waffle",
         "data.table", "ggplot2", "data.table", "stringr", "ComplexHeatmap"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/aligning_the_molecular_phenotype/Redhuis2025/")
load("data_files/FinalDataframe3Disp.RData")
load("data_files/Cleaned_Count_Data.RData")
load("data_files/Cleaned_Counts.RData")
load("data_files/Cleaned_Counts_Allele.RData")
source("functions_for_figure_scripts.R")

# Regulators are rows
# Targets are columns
# i,j = 1 means i regulates j (could be positive or negative or both depending on regulatory matrix used)

#### Data Organizing - reading in regulatory matrices ####
# Making one shared gene idx list
idxs <- rownames(counts)
grep("/", idxs, value = TRUE) # paralog pairs with too similar sequence 
# need to have separate rows for yeastract:
idxs <- map(idxs, \(g) {
  if (grepl("/", g)) {
    g <- strsplit(g, split = "/") |> unlist()
  }
  return(g)
}) |> unlist() |> sort()
write_delim(tibble(gene_name = idxs),
            delim = "\n", file = "yeastract_genes.txt",
            col_names = FALSE)

# reading in table to convert from common to systematic names
yeastract_lookup <- read_delim("data_files/downloaded_genomes_and_features/yeastract_orftogene.csv",
                               delim = ";", col_names = c("common_protein", 
                                                          "systematic", 
                                                          "common"))
yeastract_lookup$description <- select(yeastract_lookup, 
                                       setdiff(colnames(yeastract_lookup),
                                               c("common_protein", 
                                                 "systematic", 
                                                 "common"))) |> 
  apply(MARGIN = 1, \(x) reduce(x, .f = paste0)) |> as.character()
yeastract_lookup <- select(yeastract_lookup,
                           common_protein,
                           systematic,
                           common,
                           description)
# removing duplicate rows
yeastract_lookup[duplicated(yeastract_lookup$common_protein),] # 2 duplicates
yeastract_lookup <- yeastract_lookup[!duplicated(yeastract_lookup$common_protein),]

# reading in matrices
readInRegulatoryMatrix <- function(.file_name) {
  matdf <- read_delim(.file_name,
                      delim = ";")
  mat <- matdf[, -1] |> as.matrix()
  colnames(mat) <- left_join(tibble(common = colnames(mat)),
                             yeastract_lookup, by = "common") |> 
    select(systematic) |> pull()
  rownames(mat) <- left_join(tibble(common_protein = unlist(matdf[,1])),
                             yeastract_lookup, by = "common_protein") |> 
    select(systematic) |> pull()
  return(mat)
}
# setting file suffix to determine which matrix type
# file_name_suffix <- "_directOrExpression.csv" # Binding OR expression regulatory matrices
file_name_suffix <- "_bindingANDexpression.csv" # Binding AND expression regulatory matrices
# file_name_suffix <- "_Activators_directAndExpression.csv" # Binding AND expression regulatory matrices, positive results on expression only
# No-stress control matrix
ypd <- readInRegulatoryMatrix(.file_name = paste0("../../aligning_the_molecular_phenotype/yeastract/Control",
                                                  file_name_suffix))
unique(as.vector(ypd)) # binary matrix
# Diauxic Shift
hap4 <- readInRegulatoryMatrix(.file_name = paste0("../../aligning_the_molecular_phenotype/yeastract/HAP4",
                                                   file_name_suffix))
hap4_activators <- readInRegulatoryMatrix(.file_name = paste0("../../aligning_the_molecular_phenotype/yeastract/HAP4_Activators",
                                                              file_name_suffix))
hap4_inhibitors <- readInRegulatoryMatrix(.file_name = paste0("../../aligning_the_molecular_phenotype/yeastract/HAP4_Inhibitors",
                                                              file_name_suffix))
# Cell Cycle
cc <- readInRegulatoryMatrix(.file_name = paste0("../../aligning_the_molecular_phenotype/yeastract/CC",
                                                 file_name_suffix))
# Low Nitrogen
lown <- readInRegulatoryMatrix(.file_name = paste0("../../aligning_the_molecular_phenotype/yeastract/LowN",
                                                   file_name_suffix))
# Low Phosphate
lowpi <- readInRegulatoryMatrix(.file_name = paste0("../../aligning_the_molecular_phenotype/yeastract/LowPi",
                                                    file_name_suffix))
# Heat Shock
heat <- readInRegulatoryMatrix(.file_name = paste0("../../aligning_the_molecular_phenotype/yeastract/Heat",
                                                   file_name_suffix))
# Cold Shock
cold <- readInRegulatoryMatrix(.file_name = paste0("../../aligning_the_molecular_phenotype/yeastract/Cold",
                                                   file_name_suffix))

# Are activators/inhibitors complete subsets of full matrix?
# round(abs(hap4 - (hap4_activators + hap4_inhibitors))) |> table() # no

regmats <- list("YPD" = ypd, "HAP4" = hap4,
                # "HAP4_up" = hap4_activators,
                # "HAP4_down" = hap4_inhibitors,
                "CC" = cc, "LowN" = lown, "LowPi" = lowpi, 
                "Heat" = heat, "Cold" = cold)

#### Filtering out unused regulators ####
# filtering out any regulators (rows) with zero counts in all
# environments and not included in the targets (columns)
zero_count_elements <- rownames(ypd)[rowSums(purrr::reduce(regmats, .f = cbind)) == 0]
sum(zero_count_elements %in% colnames(ypd))
sum(!(zero_count_elements %in% colnames(ypd)))
toRemove <- setdiff(zero_count_elements, colnames(ypd))
regmats <- map(regmats, .f = \(x) {x[setdiff(rownames(x), toRemove),]})

#### Making regmats square ####
# filtering out zero count regulators made rows/columns nearly identical
# but there are still a few non-zero regulators that aren't found in columns

# checking for any genes missing from idxs that have network connections
setdiff(rownames(regmats$YPD)[rowSums(regmats$YPD) > 0], idxs) # YER109C = FLO8, YLR256W = HAP1, YOL028C = YAP7 all unannotated b/c of large sequence differences between Scer and other saccharomyces
# how many connections are we missing?
rowSums(regmats$YPD[!(rownames(regmats$YPD) %in% colnames(regmats$YPD)),]) # 1% is acceptable

### making regulator matrices square
regmats <- map(regmats, \(x) {
  common_idxs <- colnames(x)
  out_mat <- x[common_idxs, common_idxs]
  return(out_mat)
})

### removing unconnected nodes
# genes with no in or out edges in any environment
zero_count_rows <- rownames(regmats$YPD)[rowSums(purrr::reduce(regmats, .f = cbind)) == 0]
zero_count_cols <- colnames(regmats$YPD)[colSums(purrr::reduce(regmats, .f = rbind)) == 0]
zero_count_genes <- intersect(zero_count_rows, zero_count_cols)
length(zero_count_genes)
regmats <- map(regmats, \(x) {
  out_mat <- x[setdiff(rownames(x), zero_count_genes), 
               setdiff(colnames(x), zero_count_genes)]
  return(out_mat)
})

#### Finding network overlap ####
# how many regulatory connections are found in both each environment and YPD? 

# creating regdf, each regulator (non-zero rowSum in at least one environment),
# sum and list of its targets, environments it's active in
regdf <- map2(regmats, names(regmats), .f = \(x, y) {
  x_regulators <- rownames(x)[rowSums(x) != 0]
  outdf <- map(x_regulators, \(r) {
    targets <- colnames(x)[x[r,] != 0]
    return(tibble(regulator = r,
                  target = targets,
                  environment = y))
  }) |> purrr::reduce(.f = bind_rows)
  return(outdf)
}) |> purrr::reduce(.f = bind_rows) |> 
  # adding regulator descriptions
  left_join(select(yeastract_lookup, systematic, common, description),
                     by = c("target"="systematic")) |> 
  rename(c("target_description"="description",
           "target_common"="common")) |> 
  left_join(select(yeastract_lookup, systematic, common, description),
            by = c("regulator"="systematic")) |> 
  rename(c("regulator_description"="description",
           "regulator_common"="common"))
regdf
# How many regulators/targets are unique to each environment?
makeUpsetPlot(regdf, .group_name = "environment", 
              .group_members = names(regmats),
              .item_names = "regulator", .min_comb_size = 0)
# all environment-specific regulators are shared with YPD
makeUpsetPlot(regdf, .group_name = "environment", 
              .group_members = names(regmats),
              .item_names = "target", .min_comb_size = 0)
# most but not all targets are shared with YPD
# regulator-target pairs
makeUpsetPlot(mutate(regdf, pair = paste(regulator, target, sep = "_")), 
              .group_name = "environment", 
              .group_members = names(regmats),
              .item_names = "pair", .min_comb_size = 0)
# fewer, but still the vast majority of, regulator-target pairs 
# are shared with YPD
# all environments except Cold have at least one unique regulator-target pair
# (Cold only has 2 regulators and 3 targets total)

# TODO: remove env-specific regulatory connections from YPD

# adding regulatory info to finaldf (for this script)
finaldf_full <- finaldf
finaldf <- finaldf |> 
  select(gene_name, experiment, cer, par, dynamics,
         cor_hybrid, cor_parents, cor_scer, cor_spar) |> 
  left_join(y = regdf |> select(target, regulator, environment) |> 
              group_by(target, environment) |> # combining rows of targets with multiple regulators
              summarise(regulator = list(regulator)) |> 
              filter(environment != "YPD"),
                     by = c("gene_name"="target", 
                            "experiment"="environment"))

#### Visualizing networks of environment-specific regulators ####
getConnectedGraph <- function(.adj) {
  colnames(.adj) <- left_join(tibble(systematic = colnames(.adj)),
                             yeastract_lookup, by = "systematic") |> 
    select(common) |> pull()
  rownames(.adj) <- left_join(tibble(systematic = rownames(.adj)),
                             yeastract_lookup, by = "systematic") |> 
    select(common) |> pull()
  regnames <- rownames(.adj)[rowSums(.adj) > 0]
  g <- graph_from_adjacency_matrix(.adj)
  sub_g <- subgraph(g, which(degree(g) > 0))
  V(sub_g)$isHub <- names(V(sub_g)) %in% regnames
  V(sub_g)$vertex_size <- if_else(V(sub_g)$isHub, true = 7, false = 1)
  V(sub_g)$color <- if_else(V(sub_g)$isHub, true = "steelblue", false = "orange")
  V(sub_g)$vertex_label <- if_else(V(sub_g)$isHub, true = names(V(sub_g)), false = "")
  return(sub_g)
}
### Comented out b/c every time you re-run it, iGraph uses a different random seed

# # plotting each environment
# # downsampling to 3000 (there are 3500ish features) because igraph seems to have a problem plotting
# # plot(getConnectedGraph(regmats$YPD)) # better as a heatmap
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Yeastract/yeastract_networks.pdf",
#     width = 15, height = 20)
# par(mfrow = c(3, 2))
# for (r in setdiff(names(regmats), "YPD")) {
#   mat <- regmats[[r]]
#   if (r == "HAP4") {
#     downsample_idxs <- sample(c(1:nrow(mat)), 3000, replace = FALSE)
#     g <- getConnectedGraph(mat[downsample_idxs, downsample_idxs])
#   }
#   else {
#     g <- getConnectedGraph(mat)
#   }
#   plot(g, vertex.size = V(g)$vertex_size,
#        vertex.label = V(g)$vertex_label,
#        edge.arrow.size = 0.1,
#        main = r, layout = layout_nicely(g),
#        frame = TRUE)
# }
# dev.off()

#### Number of regulators and targets in each environmental network ####
# getting numbers for figure
for (e in ExperimentNames) {
  cat(e, "nTargets:", sum(colSums(regmats[[e]]) != 0),
      "nRegulators:", sum(rowSums(regmats[[e]]) != 0), "\n")
}

#### proportion of regulatory connections in each plasticity group per environment ####
# barplots with mutually exclusive plasticity groups:
# 1) conserved plastic, 2) conserved static, 3) Scer-unique plastic,
# 4) Spar-unique plastic, 5) plasticity reversal
# (YPD and env-specific connections in each)

# helper for plotting function
getPlasticityGroup <- function(.cer, .par) {
  if (.cer == 0 & .par == 0) {
    return("conserved static")
  }
  if (.cer == .par) {
    return("conserved plastic")
  }
  if (.cer == 0) {
    return("Spar-unique")
  }
  if (.par == 0) {
    return("Scer-unique")
  }
  else {
    return("reversal")
    }
}
plotdf <- map(ExperimentNames, \(.e) {
  expr_genes <- filter(finaldf, experiment == .e) |> 
    select(gene_name) |> pull() |> 
    intersect(colnames(regmats$YPD))
  tibble(gene_name = expr_genes,
         in_degree_e = colSums(regmats[[.e]][expr_genes, expr_genes]),
         in_degree_YPD = colSums(regmats$YPD[expr_genes, expr_genes]),
         out_degree_e = rowSums(regmats[[.e]][expr_genes, expr_genes]),
         out_degree_YPD = rowSums(regmats$YPD[expr_genes, expr_genes])) |> 
    right_join(filter(finaldf, experiment == .e),
              by = "gene_name")
}) |> purrr::reduce(.f = bind_rows)
# NA degree means gene isn't in regmat 
# (likely never had any connections in any environment),
# so setting these to a degree of 0
plotdf$in_degree_e[is.na(plotdf$in_degree_e)] <- 0
plotdf$in_degree_YPD[is.na(plotdf$in_degree_YPD)] <- 0
plotdf$out_degree_e[is.na(plotdf$out_degree_e)] <- 0
plotdf$out_degree_YPD[is.na(plotdf$out_degree_YPD)] <- 0
# getting plasticity group by gene/environment
plotdf$plasticity <- map2(plotdf$cer, plotdf$par, getPlasticityGroup) |> unlist()
totalsdf <- plotdf |> group_by(experiment) |> 
  summarise(total_in_e = sum(in_degree_e),
            total_in_YPD = sum(in_degree_YPD),
            total_out_e = sum(out_degree_e),
            total_out_YPD = sum(out_degree_YPD), 
            total_genes = n())
plotdf <- plotdf |> group_by(experiment, plasticity) |> 
  summarise(n_in_e = sum(in_degree_e),
            n_in_YPD = sum(in_degree_YPD),
            n_out_e = sum(out_degree_e),
            n_out_YPD = sum(out_degree_YPD),
            n_genes = n()) |> 
  left_join(y =  totalsdf,
            by = "experiment") |> 
  mutate(prop_in_e = n_in_e/total_in_e,
         prop_in_YPD = n_in_YPD/total_in_YPD,,
         prop_out_e = n_out_e/total_out_e,
         prop_out_YPD = n_out_YPD/total_out_YPD,
         prop_genes = n_genes/total_genes) |>
  pivot_longer(cols = c("prop_in_e", "prop_in_YPD", "prop_genes",
                        "prop_out_e", "prop_out_YPD"),
               names_to = "prop_type", values_to = "proportion") |> 
  mutate(degree = if_else(grepl("in", prop_type),
                          true = "in", 
                          false = if_else(grepl("out", prop_type),
                                          true = "out", 
                                          false = "genes")),
         network = if_else(grepl("YPD", prop_type),
                           true = "YPD", 
                           false = if_else(grepl("genes", prop_type),
                                           true = "genes",
                                           false = "env")))
totalsdf <- pivot_longer(totalsdf,
                         cols = c("total_in_e", "total_in_YPD",
                                  "total_out_e", "total_out_YPD",
                                  "total_genes"),
                         names_to = "prop_type", values_to = "n") |> 
  mutate(degree = if_else(grepl("in", prop_type),
                          true = "in", 
                          false = if_else(grepl("out", prop_type),
                                          true = "out", 
                                          false = "genes")),
         network = if_else(grepl("YPD", prop_type),
                           true = "YPD", 
                           false = if_else(grepl("genes", prop_type),
                                           true = "genes",
                                           false = "env")))
plotdf$network <- factor(plotdf$network, levels = c("genes", "YPD", "env"),
                         labels = c("number of genes", "un-stressed network\nconnections",
                                    "environment-specific network\nconnections"))
totalsdf$network <- factor(totalsdf$network, levels = c("genes", "YPD", "env"),
                           labels = c("number of genes", "un-stressed network\nconnections",
                                      "environment-specific network\nconnections"))
# in degree
p_in <- ggplot(filter(plotdf, degree %in% c("in", "genes")), aes(x = network, y = proportion)) +
  geom_bar(aes(fill = plasticity), position = "stack", stat = "identity") +
  geom_text(data = filter(totalsdf, degree %in% c("in", "genes")), aes(x = network, label = n), y = 1.1) +
  ylim(c(0, 1.1)) +
  scale_fill_discrete(limits = colordf[colordf$scheme == "dynamics",]$limits,
                      type = colordf[colordf$scheme == "dynamics",]$type) +
  ylab("%") +
  xlab("") +
  ggtitle("in degree") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2)) +
  facet_wrap(~factor(experiment, levels = ExperimentNames, labels = LongExperimentNames))

# out degree
p_out <- ggplot(filter(plotdf, degree %in% c("out", "genes")), aes(x = network, y = proportion)) +
  geom_bar(aes(fill = plasticity), position = "stack", stat = "identity") +
  geom_text(data = filter(totalsdf, degree %in% c("out", "genes")), aes(x = network, label = n), y = 1.1) +
  ylim(c(0, 1.1)) +
  scale_fill_discrete(limits = colordf[colordf$scheme == "dynamics",]$limits,
                      type = colordf[colordf$scheme == "dynamics",]$type) +
  ylab("%") +
  xlab("") +
  ggtitle("out degree") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2)) +
  facet_wrap(~factor(experiment, levels = ExperimentNames, labels = LongExperimentNames))

# plotting
pdf("../paper_figures/Yeastract/propConnectionsBars.pdf",
    width = 9, height = 7)
ggarrange(p_in, p_out, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()

#### Reversals have similar expression variance to conserved dynamics genes ####
# might be needed as a supplemental figure to point out that plasticity reveresals
# and conserved plasticity have the highest expression variation 
# across both species and environments

load("data_files/Cleaned_Counts.RData")
cer_counts <- counts[, sample_info$organism == "cer" &
                       sample_info$experiment == "LowN"]
par_counts <- counts[, sample_info$organism == "par" &
                       sample_info$experiment == "LowN"]
plotdf <- bind_rows(x = tibble(gene_name = rownames(cer_counts),
                               mean_expr = rowMeans(cer_counts),
                               var_expr = rowVars(cer_counts),
                               organism = "cer"),
                    y = tibble(gene_name = rownames(par_counts),
                               mean_expr = rowMeans(par_counts),
                               var_expr = rowVars(par_counts),
                               organism = "par")) |> 
  pivot_wider(id_cols = "gene_name",
              names_from = "organism",
              values_from = c("mean_expr", "var_expr")) |> 
  right_join(finaldf, by = c("gene_name"))
plotdf$plasticityCategory <- map2(plotdf$cer, plotdf$par, \(x, y) {
  if ((x == "1" & y == "1") |
      (x == "2" & y == "2")) {
    return("Conserved plastic")
  }
  if (x == "0" & y == "0") {
    return("Conserved static")
  }
  if ((x %in% c("1", "2") & y == "0") |
      (y %in% c("1", "2") & x == "0")) {
    return("Single-species plasticity")
  }
  if ((x == "1" & y == "2") |
      (x == "2" & y == "1")) {
    return("Plasticity reversal")
  }
}) |> unlist()
# Var/Mean by plasticity category, Scer
p1 <- ggplot(plotdf, aes(x = log2(mean_expr_cer))) + 
  geom_density(aes(fill = plasticityCategory), alpha = 0.5) +
  theme_classic() +
  xlab("Mean Expression (log2)")
p2 <- ggplot(plotdf, aes(x = log2(var_expr_cer))) + 
  geom_density(aes(fill = plasticityCategory), alpha = 0.5) +
  theme_classic() +
  xlab("Var Expression (log2)")
annotate_figure(ggarrange(p1, p2, common.legend = TRUE,
                          legend = "right", ncol = 1, nrow = 2), top = "Scer")
# Var/Mean by plasticity category, Spar
p1 <- ggplot(plotdf, aes(x = log2(mean_expr_par))) + 
  geom_density(aes(fill = plasticityCategory), alpha = 0.5) +
  theme_classic() +
  xlab("Mean Expression (log2)")
p2 <- ggplot(plotdf, aes(x = log2(var_expr_par))) + 
  geom_density(aes(fill = plasticityCategory), alpha = 0.5) +
  theme_classic() +
  xlab("Var Expression (log2)")
annotate_figure(ggarrange(p1, p2, common.legend = TRUE,
                          legend = "right", ncol = 1, nrow = 2), top = "Spar")

# TODO: could change this from n_tfs to be degree in yeastract data
# Mean
ggplot(plotdf, aes(x = log2(mean_expr), y = n_tfs)) + 
  geom_point(aes(color = dynamics))
ggplot(plotdf, aes(x = log2(mean_expr))) + 
  geom_density(aes(fill = dynamics), alpha = 0.5)
ggplot(plotdf, aes(x = log2(mean_expr))) + 
  geom_density(aes(fill = is_highly_connected), alpha = 0.5)
# Variance
ggplot(plotdf, aes(x = log2(var_expr), y = n_tfs)) + 
  geom_point(aes(color = dynamics)) +
  geom_smooth(aes(color = dynamics), method = "lm")
ggplot(plotdf, aes(x = log2(var_expr))) + 
  geom_density(aes(fill = dynamics), alpha = 0.5)
# the most highly connected genes do have much higher variance:
pdf("../../aligning_the_molecular_phenotype/paper_figures/TFdel/var_density.pdf",
    width = 3.5, height = 2)
ggplot(plotdf, aes(x = log2(var_expr))) + 
  geom_density(aes(fill = is_highly_connected), alpha = 0.5) +
  scale_fill_discrete(type = c("grey", "purple"),
                      limits = c(FALSE, TRUE),
                      labels = c(paste0("< ", ntf_cutoff, " TFs affecting gene"),
                                 paste0(ntf_cutoff, " + TFs affecting gene"))) + 
  xlab("expression variance across \nenvironments (log2 scale)") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right")
dev.off()

# Is there an excess of variance among highly connected genes, controlling for mean?
# not highly connected
ggplot(plotdf, aes(x = log2(mean_expr), y = log2(var_expr))) +
  geom_point(data = filter(plotdf, !is_highly_connected),
             aes(color = is_highly_connected), alpha = 0.1) +
  geom_smooth(method = "lm")
# highly connected
ggplot(plotdf, aes(x = log2(mean_expr), y = log2(var_expr))) +
  geom_point(data = filter(plotdf, is_highly_connected),
             aes(color = is_highly_connected), alpha = 0.1) +
  geom_smooth(method = "lm") # doesn't look like it

#### Waffle plots of environment-specific regulators ####
# in the stacked bars, we saw that regulatory out degree had very unequal representation
# of plasticity groups, particularly an over-representation of single-species plasticiy connections
# Here we want to compare regulators and their regulons from each plasticity group
# by visualizing first the regulator expression across environments 
# and then their regulon (or at least the most common plasticity pattern)

# first how many regulators from each environment are from each plasticity group?
# and how many connections do each of them have?
tfdf <- map(ExperimentNames, \(.e) {
  expr_genes <- filter(finaldf, experiment == .e) |> 
    select(gene_name) |> pull() |> 
    intersect(colnames(regmats$YPD))
  tibble(gene_name = expr_genes,
         out_degree_e = rowSums(regmats[[.e]][expr_genes, expr_genes])) |> 
    right_join(filter(finaldf, experiment == .e),
               by = "gene_name")
}) |> purrr::reduce(.f = bind_rows)
tfdf$out_degree_e[is.na(tfdf$out_degree_e)] <- 0
tfdf <- tfdf |> filter(out_degree_e != 0)
tfdf$plasticity <- map2(tfdf$cer, tfdf$par, getPlasticityGroup) |> unlist()
plotdf <- tfdf |> 
  dplyr::count(plasticity, experiment, out_degree_e, gene_name) |> 
  arrange(experiment, desc(out_degree_e))
  #arrange(experiment, plasticity)

# manually setting x and y of out degree labels
# cause waffle plot doesn't currently allow labels
# NOTE: Flip=TRUE needs to be set on geom_waffle for this to work
nrows <- 5
ncols <- max(table(plotdf$experiment))/nrows
plotdf$x <- ((c(1:nrow(plotdf)) - rank(plotdf$experiment, ties.method = "min")) %% nrows) + 1
plotdf$y <- sapply(as.numeric(table(plotdf$experiment)), \(x) {
                   full_vec <- sapply(c(1:ncols), \(i) {rep(i, times = nrows)}) |> c()
                   return(full_vec[c(1:x)])
                 }) |> unlist()
               
pdf("../paper_figures/Yeastract/regulator_waffle.pdf",
    width = 10.5, height = 2.5)  
ggplot(plotdf, aes(fill = plasticity, values = n)) +
  geom_waffle(n_rows = nrows,
    size = 1, 
    colour = "white",
    flip = TRUE) +
  geom_text(aes(x = x, y = y, label = out_degree_e), color = "white") +
  scale_fill_discrete(
    type = colordf[colordf$scheme == "dynamics",]$type,
    limits = colordf[colordf$scheme == "dynamics",]$limits) +
  facet_wrap(~factor(experiment, levels = ExperimentNames, labels = LongExperimentNames),
             nrow = 1) +
  theme_classic() +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  xlab("") +
  ylab("") +
  ggtitle("targets per regulator")
dev.off()

### Also creating YPD waffle plot where genes are colored by their plasticity group in each environment
nTargets_thresh <- 20
ypd_genes <- intersect(rownames(regmats$YPD), unique(finaldf$gene_name))
ypddf <- tibble(gene_name = ypd_genes,
                out_degree_ypd = rowSums(regmats[["YPD"]][ypd_genes, 
                                                          ypd_genes]))

ypddf$out_degree_ypd[is.na(ypddf$out_degree_ypd)] <- 0
ypddf <- ypddf |> filter(out_degree_ypd > nTargets_thresh)
ypddf$plasticity <- map2(ypddf$cer, ypddf$par, getPlasticityGroup) |> unlist()
plotdf <- tfdf |> 
  dplyr::count(plasticity, experiment, out_degree_e, gene_name) |> 
  arrange(experiment, desc(out_degree_e))
#arrange(experiment, plasticity)

# manually setting x and y of out degree labels
# cause waffle plot doesn't currently allow labels
# NOTE: Flip=TRUE needs to be set on geom_waffle for this to work
nrows <- 5
ncols <- max(table(plotdf$experiment))/nrows
plotdf$x <- ((c(1:nrow(plotdf)) - rank(plotdf$experiment, ties.method = "min")) %% nrows) + 1
plotdf$y <- sapply(as.numeric(table(plotdf$experiment)), \(x) {
  full_vec <- sapply(c(1:ncols), \(i) {rep(i, times = nrows)}) |> c()
  return(full_vec[c(1:x)])
}) |> unlist()

#### visualizing example regulator and regulon expression ####

# main points:
# 1) environment-specific regulons have conserved expression dynamics
#    in their home enivonment
# 2) Many environment-specific regulons have divergent expression dynamics
#    in other environments. These seem to follow a pattern: 
#         a) LowPi regulon has HAP4 divergence
#         b) HAP4 regulon has Heat/LowPi divergence
#         c) LowN regulon has Heat divergence (maybe a little LowPi)
#         a) Heat regulon has LowPi/HAP4/maybe some LowN divergence
# ok maybe it's just that LowPi-HAP4-LowN-Heat are the environments

#### TODO: for Fig 8D 4/23/25  ####
# TODO: consolidate this and the next two sections into one final
# section. This is where we focus on specific examples:
# 1) Pho4 has divergent dynamics in LowPi but not divergent regulon in LowPi,
#    its regulon is divergent in HAP4
# 2) GCN4 and CBF1 (big regulatory players in general) have conserved 
#    HAP4 regulon expr but diverged Heat and (GCN4 more) LowPi
# 3) TEC1, no dynamics divergence anywhere but very high LowN expr
# 4) MSN2 LowN regulon has LowPi and Heat divergence, but Heat
#    regulon only has LowPi divergence and strong Heat conserved plasticity
# 5) MSN4 Heat regulon has LowPi divergence
# 6) HSF1 Heat regulon has LowPi divergence
# 7) RAP1 has Heat plasticity reversal, and had an opposite pattern of LowPi
#    divergence in its heat regulon, maybe LowN/HAP4 divergence too
# 8) MSN2/4 have divergent trans-based LowPi and HAP4 dynamics in 5/5 regulon members

# HAP4, conserved plasticity
e <- "HAP4"
# regulator
reg_idx <- tfdf |> filter(experiment == e) |> 
  filter(out_degree_e == 383) |> 
  select(gene_name) |> pull()
filter(yeastract_lookup, systematic == reg_idx) # GCN4
plotEnvironments(.gene_idxs = reg_idx, .quartet = TRUE,
                 .normalization = "scale") # very highly expressed, slight level divergence in opposite directions for different strains. Interestingly plummets in Scer Heat

# TODO: something's wrong with centering, b/c Scer Cold can't have consistently above 0 scaled expression

# regulon
regulon_idxs <- filter(regdf, environment == e & regulator == reg_idx) |> 
  select(target) |> pull() |> 
  intersect(y = finaldf[finaldf$experiment == e,]$gene_name)
plotEnvironments(regulon_idxs, .normalization = "scale")

# HAP4, Spar unique
# regulator
reg_idx <- tfdf |> filter(experiment == e) |> 
  filter(out_degree_e == 377) |> 
  select(gene_name) |> pull()
filter(yeastract_lookup, systematic == reg_idx) # CBF1
plotEnvironments(.gene_idxs = reg_idx, .quartet = TRUE,
                 .normalization = "log2")
# regulon
regulon_idxs <- filter(regdf, environment == e & regulator == reg_idx) |> 
  select(target) |> pull() |> 
  intersect(y = finaldf[finaldf$experiment == e,]$gene_name)
plotEnvironments(.gene_idxs = regulon_idxs, .normalization = "scale")

# Conclusion: no affect on HAP4 regulon, but
# Both regulons have the same LowPi (GCN4 more) 
# and Heat dynamics divergence

# LowPi, Scer unique
e <- "LowPi"
# regulator
reg_idx <- tfdf |> filter(experiment == e) |> 
  filter(out_degree_e == 24) |> 
  select(gene_name) |> pull()
filter(yeastract_lookup, systematic == reg_idx) # Pho4
plotEnvironments(.gene_idxs = reg_idx, .quartet = TRUE,
                 .normalization = "log2")
# regulon
regulon_idxs <- filter(regdf, environment == e & regulator == reg_idx) |> 
  select(target) |> pull() |> 
  intersect(y = finaldf[finaldf$experiment == e,]$gene_name)
plotEnvironments(.gene_idxs = regulon_idxs, .normalization = "scale")

# LowN, Scer unique
e <- "LowN"
# regulator
reg_idx <- tfdf |> filter(experiment == e) |> 
  filter(out_degree_e == 373) |> 
  select(gene_name) |> pull()
filter(yeastract_lookup, systematic == reg_idx) # CBF1 again
plotEnvironments(.gene_idxs = reg_idx, .quartet = TRUE,
                 .normalization = "log2")
# regulon
regulon_idxs <- filter(regdf, environment == e & regulator == reg_idx) |> 
  select(target) |> pull() |> 
  intersect(y = finaldf[finaldf$experiment == e,]$gene_name)
plotEnvironments(.gene_idxs = regulon_idxs, .normalization = "scale")

# LowN, Spar unique
e <- "LowN"
# regulator
reg_idx <- tfdf |> filter(experiment == e) |> 
  filter(out_degree_e == 73) |> 
  select(gene_name) |> pull()
filter(yeastract_lookup, systematic == reg_idx) # TEC1
plotEnvironments(.gene_idxs = reg_idx, .quartet = TRUE,
                 .normalization = "log2") # This one just seems too lowly expressed in Scer to have plasticity
# regulon
regulon_idxs <- filter(regdf, environment == e & regulator == reg_idx) |> 
  select(target) |> pull() |> 
  intersect(y = finaldf[finaldf$experiment == e,]$gene_name)
plotEnvironments(.gene_idxs = regulon_idxs, .normalization = "scale")

# LowN, conserved plasticity
e <- "LowN"
# regulator
reg_idx <- tfdf |> filter(experiment == e) |> 
  filter(out_degree_e == 78) |> 
  select(gene_name) |> pull()
filter(yeastract_lookup, systematic == reg_idx) # MSN2
plotEnvironments(.gene_idxs = reg_idx, .quartet = TRUE,
                 .normalization = "log2")
# regulon
regulon_idxs <- filter(regdf, environment == e & regulator == reg_idx) |> 
  select(target) |> pull() |> 
  intersect(y = finaldf[finaldf$experiment == e,]$gene_name)
plotEnvironments(.gene_idxs = regulon_idxs, .normalization = "scale") # LowPi and Heat diverged

# Heat, conserved static
e <- "Heat"
# regulator
reg_idx <- tfdf |> filter(experiment == e) |> 
  filter(out_degree_e == 138) |> 
  select(gene_name) |> pull()
filter(yeastract_lookup, systematic == reg_idx) # also MSN2
plotEnvironments(.gene_idxs = reg_idx, .quartet = TRUE,
                 .normalization = "log2")
# regulon
regulon_idxs <- filter(regdf, environment == e & regulator == reg_idx) |> 
  select(target) |> pull() |> 
  intersect(y = finaldf[finaldf$experiment == e,]$gene_name)
plotEnvironments(.gene_idxs = regulon_idxs, .normalization = "scale") # LowPi but not Heat diverged. 
# Nice b/c it demonstrates that the enivronment they regulate in 
# doesn't have the divergence, it's other environments

# Heat, conserved plastic
e <- "Heat"
# regulator
reg_idx <- tfdf |> filter(experiment == e) |> 
  filter(out_degree_e == 91) |> 
  select(gene_name) |> pull()
filter(yeastract_lookup, systematic == reg_idx) # MSN4
plotEnvironments(.gene_idxs = reg_idx, .quartet = TRUE,
                 .normalization = "log2")
# regulon
regulon_idxs <- filter(regdf, environment == e & regulator == reg_idx) |> 
  select(target) |> pull() |> 
  intersect(y = finaldf[finaldf$experiment == e,]$gene_name)
plotEnvironments(.gene_idxs = regulon_idxs, .normalization = "scale") # also LowPi

# Heat, Spar unique
e <- "Heat"
# regulator
reg_idx <- tfdf |> filter(experiment == e) |> 
  filter(out_degree_e == 89) |> 
  select(gene_name) |> pull()
filter(yeastract_lookup, systematic == reg_idx) # HSF1
plotEnvironments(.gene_idxs = reg_idx, .quartet = TRUE,
                 .normalization = "log2")
# regulon
regulon_idxs <- filter(regdf, environment == e & regulator == reg_idx) |> 
  select(target) |> pull() |> 
  intersect(y = finaldf[finaldf$experiment == e,]$gene_name)
plotEnvironments(.gene_idxs = regulon_idxs, .normalization = "scale") # also LowPi

# Heat, Reversal
e <- "Heat"
# regulator
reg_idx <- tfdf |> filter(experiment == e) |> 
  filter(out_degree_e == 74) |> 
  select(gene_name) |> pull()
filter(yeastract_lookup, systematic == reg_idx) # RAP1
plotEnvironments(.gene_idxs = reg_idx, .quartet = TRUE,
                 .normalization = "log2")
# regulon
regulon_idxs <- filter(regdf, environment == e & regulator == reg_idx) |> 
  select(target) |> pull() |> 
  intersect(y = finaldf[finaldf$experiment == e,]$gene_name)
plotEnvironments(.gene_idxs = regulon_idxs, .normalization = "scale") # reversal is seen in LowPi and a bit in HAP4

# Cold, msn2
e <- "Cold"
# regulator
reg_idx <- tfdf |> filter(experiment == e) |> 
  filter(out_degree_e == 3) |> 
  select(gene_name) |> pull()
filter(yeastract_lookup, systematic == reg_idx) # MSN2
plotEnvironments(.gene_idxs = reg_idx, .quartet = TRUE,
                 .normalization = "log2")
# regulon, as individual genes
regulon_idxs <- filter(regdf, environment == e & regulator == reg_idx) |> 
  select(target) |> pull() |> 
  intersect(y = finaldf[finaldf$experiment == e,]$gene_name)
plotEnvironments(.gene_idxs = regulon_idxs[1], .normalization = "scale", .quartet = TRUE)
plotEnvironments(.gene_idxs = regulon_idxs[2], .normalization = "scale", .quartet = TRUE)
plotEnvironments(.gene_idxs = regulon_idxs[3], .normalization = "scale", .quartet = TRUE) # lowPi and (kind of) HAP4 dynamics divergence, Scer dominant

# Cold, msn4
e <- "Cold"
# regulator
reg_idx <- tfdf |> filter(experiment == e) |> 
  filter(out_degree_e == 2) |> 
  select(gene_name) |> pull()
filter(yeastract_lookup, systematic == reg_idx) # MSN4
plotEnvironments(.gene_idxs = reg_idx, .quartet = TRUE,
                 .normalization = "log2")
# regulon, as individual genes
regulon_idxs <- filter(regdf, environment == e & regulator == reg_idx) |> 
  select(target) |> pull() |> 
  intersect(y = finaldf[finaldf$experiment == e,]$gene_name)
plotEnvironments(.gene_idxs = regulon_idxs[1], .normalization = "scale", .quartet = TRUE)
plotEnvironments(.gene_idxs = regulon_idxs[2], .normalization = "scale", .quartet = TRUE) # same pattern as MSN2

#### Exploring examples of environmental sensors and their targets ####
# In each environment, plot expression profiles of selected environmental sensors followed by
# Expression of the group of genes each sensor targets

# First example: Gcr1 in Saturated Growth
r <- "HAP4"
sensor_idx <- yeastract_lookup |> filter(common == "GCR1") |> 
  select(systematic) |> pull()
target_idxs <- colnames(regmats[[r]])[regmats[[r]][sensor_idx, ] != 0]
plotGenes(sensor_idx, .quartet = TRUE, .experiment_name = r)
target_idx <- target_idxs[sample(c(1:length(target_idxs)), 1)]
plotGenes(target_idx, .quartet = TRUE, 
          .experiment_name = r, .plot_titles = target_idx) 
# random gene to gauge how responsive these GCR1 targets are:
random_idx <- sample(rownames(regmats[[r]]), 1)
plotGenes(random_idx, .quartet = TRUE, 
          .experiment_name = r, .plot_titles = random_idx) # random genes seem to be lower expressed
# very clearly a mix of some up some down:
finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  select(cer, par) |> 
  table()
# %conserved dynamics and mean parental cor
finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  select(dynamics, cor_parents) |> 
  summarise(pct_consdyn = sum(dynamics == "conserved")/length(dynamics),
            mean_cor = mean(cor_parents))
# same numbers for entire genome in HAP4
finaldf |> filter(experiment == r) |> 
  select(dynamics, cor_parents) |> 
  summarise(pct_consdyn = sum(dynamics == "conserved")/length(dynamics),
            mean_cor = mean(cor_parents, na.rm = TRUE)) # much lower
# grouping by dynamics
target_idxs_11 <- finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  filter(cer == 1 & par == 1) |> select(gene_name) |> pull()
plotGenes(target_idxs_11, .quartet = TRUE, 
          .experiment_name = r, .plot_titles = paste(length(target_idxs_11), "genes")) 
target_idxs_22 <- finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  filter(cer == 2 & par == 2) |> select(gene_name) |> pull()
plotGenes(target_idxs_22, .quartet = TRUE, 
          .experiment_name = r, .plot_titles = paste(length(target_idxs_22), "genes")) 
# diverging dynamics
target_idxs_12 <- finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  filter(cer == 1 & par == 2) |> select(gene_name) |> pull()
plotGenes(target_idxs_12[1], .quartet = TRUE, 
          .experiment_name = r, .plot_titles = paste(length(target_idxs_12), "genes")) 
plotGenes(target_idxs_12[2], .quartet = TRUE, 
          .experiment_name = r, .plot_titles = paste(length(target_idxs_12), "genes")) 
target_idxs_21 <- finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  filter(cer == 2 & par == 1) |> select(gene_name) |> pull()
plotGenes(target_idxs_21[1], .quartet = TRUE, 
          .experiment_name = r, .plot_titles = paste(length(target_idxs_21), "genes")) 
plotGenes(target_idxs_21[2], .quartet = TRUE, 
          .experiment_name = r, .plot_titles = paste(length(target_idxs_21), "genes")) 
# Conclusion: While GCR1 is fairly lowly expressed and doesn't change much at the diauxic shift,
# The majority of its targets are very responsive to the diauxic shift
# (and have quite conserved dynamics)
# checking other environments
plotGenes(target_idxs_11, .quartet = TRUE, 
          .experiment_name = "CC", 
          .plot_titles = "GCR1 regulon in CC") # Scer-specific spike in CC Urea shock
plotGenes(target_idxs_11, .quartet = TRUE, 
          .experiment_name = "LowN", 
          .plot_titles = "GCR1 regulon in LowN")
plotGenes(target_idxs_11, .quartet = TRUE, 
          .experiment_name = "LowPi", 
          .plot_titles = "GCR1 regulon in LowPi")
plotGenes(target_idxs_11, .quartet = TRUE, 
          .experiment_name = "Heat", 
          .plot_titles = "GCR1 regulon in Heat")
plotGenes(target_idxs_11, .quartet = TRUE, 
          .experiment_name = "Cold", 
          .plot_titles = "GCR1 regulon in Cold")

# Second example: PHO4 in LowPi
r <- "LowPi"
sensor_idx <- yeastract_lookup |> filter(common == "PHO4") |> 
  select(systematic) |> pull()
target_idxs <- colnames(regmats[[r]])[regmats[[r]][sensor_idx, ] != 0]
plotGenes(sensor_idx, .quartet = TRUE, .experiment_name = r)
target_idx <- target_idxs[sample(c(1:length(target_idxs)), 1)]
plotGenes(target_idx, .quartet = TRUE, 
          .experiment_name = r, .plot_titles = target_idx)
# random gene to gauge how responsive these GCR1 targets are:
random_idx <- sample(rownames(regmats[[r]]), 1)
plotGenes(random_idx, .quartet = TRUE, 
          .experiment_name = r, .plot_titles = random_idx) # random genes seem to be lower expressed
# mostly conserved up:
finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  select(cer, par) |> 
  table()
# %conserved dynamics and mean parental cor
finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  select(dynamics, cor_parents) |> 
  summarise(pct_consdyn = sum(dynamics == "conserved")/length(dynamics),
            mean_cor = mean(cor_parents))
# same numbers for entire genome in this environment
finaldf |> filter(experiment == r) |> 
  select(dynamics, cor_parents) |> 
  summarise(pct_consdyn = sum(dynamics == "conserved")/length(dynamics),
            mean_cor = mean(cor_parents, na.rm = TRUE)) # much lower
# grouping by dynamics
target_idxs_11 <- finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  filter(cer == 1 & par == 1) |> select(gene_name) |> pull()
plotGenes(target_idxs_11, .quartet = TRUE, 
          .experiment_name = r, .plot_titles = paste(length(target_idxs_11), "genes")) 
# diverging dynamics
target_idxs_21 <- finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  filter(cer == 2 & par == 1) |> select(gene_name) |> pull()
plotGenes(target_idxs_21[1], .quartet = TRUE, 
          .experiment_name = r, .plot_titles = target_idxs_21[1]) 
plotGenes(target_idxs_21[2], .quartet = TRUE, 
          .experiment_name = r, .plot_titles = target_idxs_21[2]) 
plotGenes(target_idxs_21[3], .quartet = TRUE, 
          .experiment_name = r, .plot_titles = target_idxs_21[3])

# This was an accident, but these PHO4 connected genes have Spar-dominant divergent
# dynamics in Saturated Growth (but not LowPi)
finaldf |> filter(experiment == "HAP4" & gene_name %in% target_idxs) |> 
  select(cer, par) |> 
  table() # a lot of 2-1 divergers (LowPi had nearly all 1-1 conserved)
plotGenes(target_idxs, .quartet = TRUE, .experiment_name = "HAP4", 
          .plot_titles = paste(length(target_idxs), "genes in Saturated Growth"))

# Would this be more or less true for all PHO4-regulated genes across environments?
target_idxs <- colnames(regmats$YPD)[regmats$YPD[sensor_idx, ] != 0]
finaldf |> filter(experiment == "HAP4" & gene_name %in% target_idxs) |> 
  select(cer, par) |> 
  table()
target_idxs_11 <- finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  filter(cer == 1 & par == 1) |> select(gene_name) |> pull()
plotGenes(target_idxs_11, .quartet = TRUE, .experiment_name = "LowPi", 
          .plot_titles = paste(length(target_idxs), "genes in Low Phosphate"))
plotGenes(target_idxs_11, .quartet = TRUE, .experiment_name = "HAP4", 
          .plot_titles = paste(length(target_idxs), "genes in Saturated Growth"))
# less true, but still there
# is PHO4 differentially expressed btwn Scer and Spar in HAP4?
plotEnvironments(sensor_idx, .normalization = "log2")
# PHO4 has higher level in Spar:
finaldf_full |> filter(gene_name == sensor_idx) |> select(level, effect_size_species)

# Third example: HSF1 in Heat
r <- "Heat"
sensor_idx <- yeastract_lookup |> filter(common == "HSF1") |> 
  select(systematic) |> pull()
target_idxs <- colnames(regmats[[r]])[regmats[[r]][sensor_idx, ] != 0]
plotGenes(sensor_idx, .quartet = TRUE, .experiment_name = r)
target_idx <- target_idxs[sample(c(1:length(target_idxs)), 1)]
plotGenes(target_idx, .quartet = TRUE, 
          .experiment_name = r, .plot_titles = target_idx)
# random gene to gauge how responsive these GCR1 targets are:
random_idx <- sample(rownames(regmats[[r]]), 1)
plotGenes(random_idx, .quartet = TRUE, 
          .experiment_name = r, .plot_titles = random_idx) # random genes seem to be lower expressed
# mostly 3-3s:
finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  select(cer, par) |> 
  table()
# %conserved dynamics and mean parental cor
finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  select(dynamics, cor_parents) |> 
  summarise(pct_consdyn = sum(dynamics == "conserved")/length(dynamics),
            mean_cor = mean(cor_parents))
# same numbers for entire genome in this environment
finaldf |> filter(experiment == r) |> 
  select(dynamics, cor_parents) |> 
  summarise(pct_consdyn = sum(dynamics == "conserved")/length(dynamics),
            mean_cor = mean(cor_parents, na.rm = TRUE)) # much lower
# grouping by dynamics
target_idxs_33 <- finaldf |> filter(experiment == r & gene_name %in% target_idxs) |> 
  filter(cer == 3 & par == 3) |> select(gene_name) |> pull()
plotGenes(target_idxs_33, .quartet = TRUE, 
          .experiment_name = r, .plot_titles = paste(length(target_idxs_33), "genes"))
# Conclusion: Hsf1 spikes early in response to heat
# What about other environments?
plotGenes(target_idxs_33, .quartet = TRUE, 
          .experiment_name = "HAP4", .plot_titles = "HAP4")
plotGenes(target_idxs_33, .quartet = TRUE, 
          .experiment_name = "CC", .plot_titles = "CC") # level change seen in HAP4 btwn Scer/Spar isn't constitutive
plotGenes(target_idxs_33, .quartet = TRUE, 
          .experiment_name = "LowN", .plot_titles = "LowN")
plotGenes(target_idxs_33, .quartet = TRUE, 
          .experiment_name = "LowPi", .plot_titles = "LowPi")
plotGenes(target_idxs_33, .quartet = TRUE, 
          .experiment_name = "Cold", .plot_titles = "Cold")

#### Regulon expression in home and other environments ####
visualizeRegulon <- function(.regulator, .environment, .clust, .df = regdf,
                             .normalization = "log2", .plotlims = NULL) {
  regulator_systematic <- yeastract_lookup |> filter(common == .regulator) |> 
    select(systematic) |> pull()
  target_idxs <- .df |> filter(regulator == regulator_systematic & 
                                 environment == .environment) |> 
    select(target) |> pull()
  idxs <- finaldf |> filter(gene_name %in% target_idxs & experiment == .environment &
                              cer == .clust & par == .clust) |> 
    select(gene_name) |> pull()
  null_idxs <- finaldf |> filter(experiment == .environment &
                                   cer == .clust & par == .clust) |> 
    slice_sample(n = length(idxs)) |> select(gene_name) |> pull()
  p <- plotExpressionProfileQuartet(.cts1 = collapsed$cer[idxs, , drop = FALSE],
                               .cts2 = collapsed$par[idxs, , drop = FALSE],
                               .cts3 = collapsed_allele$cer[idxs, , drop = FALSE],
                               .cts4 = collapsed_allele$par[idxs, , drop = FALSE],
                               .info1 = info,
                               .info2 = info,
                               .info3 = info_allele,
                               .info4 = info_allele,
                               .name1 = "S. cer",
                               .name2 = "S. par",
                               .name3 = "F1, cer allele",
                               .name4 = "F1, par allele",
                               .color1 = "orange1",
                               .color2 = "blue2",
                               .color3 = "orange4",
                               .color4 = "blue4",
                               .normalization = .normalization,
                               .method = "line",
                               .show_points = FALSE,
                               .show_confidence_intervals = TRUE,
                               .plotlims = .plotlims,
                               .plot_titles = "experiment")
  p_null <- plotExpressionProfileQuartet(.cts1 = collapsed$cer[null_idxs, , drop = FALSE],
                                    .cts2 = collapsed$par[null_idxs, , drop = FALSE],
                                    .cts3 = collapsed_allele$cer[null_idxs, , drop = FALSE],
                                    .cts4 = collapsed_allele$par[null_idxs, , drop = FALSE],
                                    .info1 = info,
                                    .info2 = info,
                                    .info3 = info_allele,
                                    .info4 = info_allele,
                                    .name1 = "S. cer",
                                    .name2 = "S. par",
                                    .name3 = "F1, cer allele",
                                    .name4 = "F1, par allele",
                                    .color1 = "orange1",
                                    .color2 = "blue2",
                                    .color3 = "orange4",
                                    .color4 = "blue4",
                                    .normalization = .normalization,
                                    .method = "line",
                                    .show_points = FALSE,
                                    .show_confidence_intervals = TRUE,
                                    .plotlims = .plotlims,
                                    .plot_titles = "experiment")
  return(annotate_figure(ggarrange(p, p_null, nrow = 2, ncol = 1),
                         top = paste(length(idxs), "genes, top = regulon, bottom = random cluster sample")))
}
# examples from each environment
visualizeRegulon(.regulator = "HSF1", .environment = "Heat", .clust = 3)
visualizeRegulon(.regulator = "PHO4", .environment = "LowPi", .clust = 1)
visualizeRegulon(.regulator = "STE12", .environment = "CC", .clust = 2)
visualizeRegulon(.regulator = "MCM1", .environment = "CC", .clust = 1)
visualizeRegulon(.regulator = "GCN4", .environment = "HAP4", .clust = 2)
visualizeRegulon(.regulator = "PHD1", .environment = "LowN", .clust = 2)
visualizeRegulon(.regulator = "MSN2", .environment = "Cold", .clust = 1)

#### Divergent dynamics clusters have overrepresentation of other environment regulons ####

# TODO: same as above, just adapt this code so that
# it calculates for setdiff(all environments, current environment)
plotdf <- regdf |> filter(environment != "YPD") |> 
  group_by(environment, regulator) |> 
  summarise(genome = calculateConsDyn(.env = unique(environment)),
            regulon = calculateConsDyn(.env = unique(environment),
                                       .idxs = target),
            regulon_size = sum(target %in% unlist(finaldf[finaldf$experiment == unique(environment), "gene_name"]))) |> 
  filter(regulon_size > 1)

########################### Archive ##############################
# #### Verifying environment sensing in Saturated Growth ####
# archived b/c this was mainly exploratory, but might be useful if we do it for all enivonments
# not very precise though, b/c without using CHiP data like yeastract does,
# it relies on different paper authors having the same definiton of env sensing
# # verifying that downsampling HAP4 matrix makes for a more organized graph
# downsample_idxs <- sample(c(1:nrow(regmats$HAP4)), 3000, replace = FALSE)
# missing_idxs <- setdiff(1:nrow(regmats$HAP4), downsample_idxs)
# rowSums(regmats$HAP4[missing_idxs, missing_idxs]) |> table()
# colSums(regmats$HAP4[missing_idxs, missing_idxs]) |> table() # if no connections are removed, downsampling won't help (we checked)
# hap4_g <- getConnectedGraph(regmats$HAP4[downsample_idxs, downsample_idxs])
# plot(hap4_g, vertex.size = V(hap4_g)$vertex_size,
#      vertex.label = V(hap4_g)$vertex_label,
#      edge.arrow.size = 0.1,
#      main = "Saturated Growth", layout = layout_nicely(hap4_g),
#      frame = TRUE)
# # number of regulators
# yeastract_regulators <- rownames(regmats$HAP4)[rowSums(regmats$HAP4) > 0]
# length(yeastract_regulators)
# # range of number of genes each regulates
# quantile(rowSums(regmats$HAP4)[yeastract_regulators])
# hist(rowSums(regmats$HAP4)[yeastract_regulators])
# rownames(regmats$HAP4)[which(rowSums(regmats$HAP4) == 444)] # GCN4
# 
# # Yeastract environment: Carbon source quality/availability > Non-fermentable carbon source
# # (diauxic shift environment in yeastract has very few connections)
# 
# # environment sensors from the literature:
# yeastract_lookup |> filter(common %in% c("TPK1", "TPK2", "TPK3", "TOR1", "TOR2", "RIM15", "SNF1", "MIG1", "GCR1", "ADR1"))
# literature_sensors <- yeastract_lookup |> filter(common %in% c("TPK1", "TPK2", "TPK3", "TOR1", "TOR2", "RIM15", "SNF1", "MIG1", "GCR1", "ADR1")) |> 
#   select(systematic) |> pull()
# setdiff(literature_sensors, colnames(regmats$HAP4))
# literature_sensors <- intersect(literature_sensors, colnames(regmats$HAP4))
# intersect(yeastract_regulators, literature_sensors)
# yeastract_lookup |> filter(systematic %in% literature_sensors) |> View()
# # Note: even when not requiring there to be DNA binding evidence, there is a bias 
# # for interactions where the regulator is a transcription factor that regulates expression via DNA binding---
# # the TORs, PKAs, and RIM15 are all kinases that don't bind DNA and don't have any yeastract binding interactions
# regmats$HAP4[literature_sensors,] |> sum() # nGenes sensors regulate
# regmats$HAP4[,literature_sensors] |> sum() # nSensors regulated
# regmats$HAP4[yeastract_regulators,] |> sum() # nGenes yeastract hubs regulate
# regmats$HAP4[,yeastract_regulators] |> sum() # nHubs regulated
# 
# # sensor targets:
# targets <- yeastract_lookup |> filter(systematic %in% rownames(regmats$HAP4)[colSums(regmats$HAP4[literature_sensors,]) != 0]) |> 
#   select(common) |> pull()
# # sensor regulators:
# yeastract_lookup |> filter(systematic %in% colnames(regmats$HAP4)[rowSums(regmats$HAP4[,literature_sensors]) != 0]) |> 
#   View()
# # targets per sensor:
# rowSums(regmats$HAP4[literature_sensors,])
# yeastract_lookup |> filter(systematic == "YPL075W")
#### Separating Activators and Inhibitors ####
# Archived b/c A) Activator/Inhibitor matrices are not perfect subsets of combination,
# and B) Activator/Inhibitor distinction did not relate to direction of expression change
# # Do inhibitory interactions go down and activating interactions go up?
# r <- "HAP4_down"
# sensor_idx <- yeastract_lookup |> filter(common == "GCR1") |> 
#   select(systematic) |> pull()
# target_idxs <- colnames(regmats[[r]])[regmats[[r]][sensor_idx, ] != 0]
# plotGenes(sensor_idx, .quartet = TRUE, .experiment_name = "HAP4")
# target_idx <- target_idxs[sample(c(1:length(target_idxs)), 1)]
# plotGenes(target_idx, .quartet = TRUE, 
#           .experiment_name = "HAP4", .plot_titles = target_idx) 
# # still a mix of some up some down:
# finaldf |> filter(experiment == "HAP4" & gene_name %in% target_idxs) |> 
#   select(cer, par) |> 
#   table()
# # repeat with activators
# r <- "HAP4_up"
# sensor_idx <- yeastract_lookup |> filter(common == "GCR1") |> 
#   select(systematic) |> pull()
# target_idxs <- colnames(regmats[[r]])[regmats[[r]][sensor_idx, ] != 0]
# plotGenes(sensor_idx, .quartet = TRUE, .experiment_name = "HAP4")
# target_idx <- target_idxs[sample(c(1:length(target_idxs)), 1)]
# plotGenes(target_idx, .quartet = TRUE, 
#           .experiment_name = "HAP4", .plot_titles = target_idx) 
# # still a mix of some up some down:
# finaldf |> filter(experiment == "HAP4" & gene_name %in% target_idxs) |> 
#   select(cer, par) |> 
#   table()

#### Probably Archive: Proportions of regulatory connections in conserved vs diverged dynamics genes ####
# # archived b/c the stacked bars with more than just conserved/diverged
# # dynamics are more informative
# 
# # stacked barplots from each environment with the following columns:
# #    1) fraction of genes with conserved vs diverged dynamics
# #    2-5) fraction of regulatory connections pointing towards 
# #         genes with conserved vs diverged dynamics
# #         with the following filters:
# #   2) YPD and specific environment connections combined
# #   3) just YPD connections
# #   4) just specific environment connections
# #   5) diverged dynamics are just the reversals
# # Then repeat for second plot but instead of conserved vs diverged, 
# # it's low varying in either species vs not low varying in either species
# 
# plotConnectionsBars <- function(.e, .col = "dynamics", .df = finaldf,
#                                 .cons_fill = colordf[colordf$limits == "conserved level and dynamics",]$type,
#                                 .div_fill = colordf[colordf$limits == "conserved level, diverged dynamics",]$type) {
#   low_expr <- setdiff(colnames(regmats[[.e]]), .df[.df$experiment == .e, ]$gene_name)
#   plotdf <- tibble(gene_name = colnames(regmats[[.e]])[!(colnames(regmats[[.e]]) %in% low_expr)],
#                    in_degree_e = colSums(regmats[[.e]])[!(colnames(regmats[[.e]]) %in% low_expr)],
#                    in_degree_YPD = colSums(regmats$YPD)[!(colnames(regmats[[.e]]) %in% low_expr)],
#                    in_degree_full = colSums(regmats[[.e]] + regmats$YPD)[!(colnames(regmats[[.e]]) %in% low_expr)],
#                    n_genes = 1) |> 
#     left_join(filter(.df, experiment == .e),
#               by = "gene_name") |> 
#     mutate(fill_col = .data[[.col]])
#   plotdf_consdiv <- plotdf |> 
#     pivot_longer(cols = c("in_degree_e", "in_degree_YPD", "in_degree_full", "n_genes"),
#                  names_to = "type", values_to = "count") |> group_by(fill_col, type) |> 
#     summarise(in_degree = sum(count)) |> 
#     group_by(type) |> 
#     reframe(n = sum(in_degree), in_degree = in_degree,
#             fill_col = fill_col, type = type) |> 
#     mutate(frac_in_degree = in_degree/n)
#   frac_fill <- plotdf_consdiv |> filter(type == "n_genes") |> select(in_degree) |> pull()
#   frac_fill <- min(frac_fill)/sum(frac_fill)
#   p <- ggplot(plotdf_consdiv, aes(x = type, y = frac_in_degree)) +
#     geom_bar(aes(fill = fill_col), stat = "identity") +
#     geom_hline(yintercept = frac_fill) +
#     scale_fill_discrete(limits = c("conserved", "diverged"),
#                         type = c(.cons_fill, .div_fill))
#   return(p)
# }
# # conserved dynamics have more regulatory connections
# plotConnectionsBars(.e = "HAP4")
# plotConnectionsBars(.e = "CC")
# plotConnectionsBars(.e = "LowN")
# plotConnectionsBars(.e = "LowPi")
# plotConnectionsBars(.e = "Heat")
# plotConnectionsBars(.e = "Cold")
# # is this more or less pronounced for high varying vs low varying?
# plotConnectionsBars(.e = "HAP4", 
#                     .df = mutate(finaldf, isLow = (cer == 0 | par == 0)),
#                     .col = "isLow")
# plotConnectionsBars(.e = "CC", 
#                     .df = mutate(finaldf, isLow = (cer == 0 | par == 0)),
#                     .col = "isLow")
# plotConnectionsBars(.e = "LowN", 
#                     .df = mutate(finaldf, isLow = (cer == 0 | par == 0)),
#                     .col = "isLow")
# plotConnectionsBars(.e = "LowPi", 
#                     .df = mutate(finaldf, isLow = (cer == 0 | par == 0)),
#                     .col = "isLow")
# plotConnectionsBars(.e = "Heat", 
#                     .df = mutate(finaldf, isLow = (cer == 0 | par == 0)),
#                     .col = "isLow")
# plotConnectionsBars(.e = "Cold", 
#                     .df = mutate(finaldf, isLow = (cer == 0 | par == 0)),
#                     .col = "isLow")
# # more pronounced except for Cold YPD regulon
# # is the pattern still there if we filter out lowvar genes?
# plotConnectionsBars(.e = "HAP4", 
#                     .df = filter(finaldf, !(cer == 0 | par == 0)),
#                     .col = "dynamics")
# plotConnectionsBars(.e = "CC", 
#                     .df = filter(finaldf, !(cer == 0 | par == 0)),
#                     .col = "dynamics")
# plotConnectionsBars(.e = "LowN", 
#                     .df = filter(finaldf, !(cer == 0 | par == 0)),
#                     .col = "dynamics")
# plotConnectionsBars(.e = "LowPi", 
#                     .df = filter(finaldf, !(cer == 0 | par == 0)),
#                     .col = "dynamics")
# plotConnectionsBars(.e = "Heat", 
#                     .df = filter(finaldf, !(cer == 0 | par == 0)),
#                     .col = "dynamics")
# plotConnectionsBars(.e = "Cold", 
#                     .df = filter(finaldf, !(cer == 0 | par == 0)),
#                     .col = "dynamics")
# # YPD connections are very even now
# # env-specific connections depend on environment:
# # HAP4, LowN, no distribution difference
# # CC, looks the same as without filtering, b/c there are very few 0 genes
# # LowPi, Heat, have way fewer connections in divergent dynamics genes
# # Cold has way more connections in divergent dynamics genes (but the only env sensors are MSN2/4)
# 
# # What if we only look at lowVar vs single species plasiticty?
# plotConnectionsBars(.e = "HAP4", 
#                     .df = filter(finaldf, (cer == 0 | par == 0)),
#                     .col = "dynamics")
# plotConnectionsBars(.e = "CC", 
#                     .df = filter(finaldf, (cer == 0 | par == 0)),
#                     .col = "dynamics")
# plotConnectionsBars(.e = "LowN", 
#                     .df = filter(finaldf, (cer == 0 | par == 0)),
#                     .col = "dynamics")
# plotConnectionsBars(.e = "LowPi", 
#                     .df = filter(finaldf, (cer == 0 | par == 0)),
#                     .col = "dynamics")
# plotConnectionsBars(.e = "Heat", 
#                     .df = filter(finaldf, (cer == 0 | par == 0)),
#                     .col = "dynamics")
# plotConnectionsBars(.e = "Cold", 
#                     .df = filter(finaldf, (cer == 0 | par == 0)),
#                     .col = "dynamics")
# # Single-species plasticity have way more connections than cons0, but LowPi and Cold
# # have no environment-specific connections in 0 genes

# #### Archived: Environmental sensor targets (regulons) tend to have conserved dynamics in their home environment ####
# # Archived b/c this is a less-intuitive version of the
# # stacked barplots that have even more information about regulatory dynamics
# 
# # This seems to be true from observation but we have to verify
# # by checking if the targets of each environmental regulator have
# # a higher chance of being conserved dynamics (i.e. > genome avg of 70%)
# 
# # for each regulator in regdf, calculating conserved dynamics
# # for regulon versus entire genome 
# calculateConsDyn <- function(.env, .idxs = NULL) {
#   if (!is.null(.idxs)) {
#     output <- finaldf[finaldf$experiment == .env & 
#                         finaldf$gene_name %in% .idxs,] |> 
#       select(dynamics, cor_parents) |> 
#       summarise(pct_consdyn = sum(dynamics == "conserved")/length(dynamics),
#                 mean_cor = mean(cor_parents, na.rm = TRUE))
#   }
#   else {
#     output <- finaldf[finaldf$experiment == .env,] |> 
#       select(dynamics, cor_parents) |> 
#       summarise(pct_consdyn = sum(dynamics == "conserved")/length(dynamics),
#                 mean_cor = mean(cor_parents, na.rm = TRUE))
#   }
#   return(output)
# }
# # tests for 
# calculateConsDyn(.env = "CC")
# # regulon
# test_regulon <- regdf |> filter(regulator == "YDR207C" &
#                                   environment == "CC") |> 
#   select(target) |> pull()
# finaldf[finaldf$experiment == "CC" & 
#           finaldf$gene_name %in% test_regulon,] |> 
#   select(dynamics, cor_parents) # no genes high enough expressed in this regulon
# calculateConsDyn(.env = "CC", .idxs = test_regulon)
# 
# plotdf <- regdf |> filter(environment != "YPD") |> 
#   group_by(environment, regulator) |> 
#   summarise(genome = calculateConsDyn(.env = unique(environment)),
#             regulon = calculateConsDyn(.env = unique(environment),
#                                        .idxs = target),
#             regulon_size = sum(target %in% unlist(finaldf[finaldf$experiment == unique(environment), "gene_name"]))) |> 
#   filter(regulon_size > 1)
# # splitting results
# plotdf$genome_pct_consdyn <- plotdf$genome$pct_consdyn
# plotdf$regulon_pct_consdyn <- plotdf$regulon$pct_consdyn
# plotdf$genome_mean_cor <- plotdf$genome$mean_cor
# plotdf$regulon_mean_cor <- plotdf$regulon$mean_cor
# plotdf <- select(plotdf, -genome, -regulon)
# 
# # paired t-tests of each property
# # % conserved dynamics
# tdf <- plotdf |> pivot_longer(cols = c("genome_pct_consdyn", 
#                                        "regulon_pct_consdyn")) |> 
#   mutate(name = factor(name, levels = c("regulon_pct_consdyn",
#                                         "genome_pct_consdyn")))
# t_pct_consdyn <- t.test(value ~ name,
#        data = tdf, alternative = "greater") # greater checks if first factor level is greater than second
# 
# # mean correlation
# tdf <- plotdf |> pivot_longer(cols = c("genome_mean_cor", 
#                                        "regulon_mean_cor")) |> 
#   mutate(name = factor(name, levels = c("regulon_mean_cor",
#                                         "genome_mean_cor")))
# t_mean_cor <- t.test(value ~ name,
#        data = tdf, alternative = "greater") # greater checks if first factor level is greater than second
# 
# # Supplementary figure: percent conserved dynamics
# # in regulon vs genome in each environment
# pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/regulon_consdyn.pdf",
#     width = 7, height = 4)
# p_pct <- ggplot(plotdf, aes(x = genome_pct_consdyn,
#                    y = regulon_pct_consdyn)) + 
#   geom_point(aes(color = environment)) +
#   xlim(c(0, 1)) + ylim(c(0, 1)) +
#   xlab("all genes") +
#   ylab("regulon") +
#   geom_abline(slope = 1, intercept = 0) +
#   theme_classic() + 
#   ggtitle(paste("% conserved dynamics\npaired t.test p <=", round_any(t_pct_consdyn$p.value, 
#                                                        accuracy = 0.000001, f = ceiling)))
# p_avgCor <- ggplot(plotdf, aes(x = genome_mean_cor,
#                             y = regulon_mean_cor)) + 
#   geom_point(aes(color = environment)) +
#   xlim(c(-1, 1)) + ylim(c(-1, 1)) +
#   xlab("all genes") +
#   ylab("regulon") +
#   geom_abline(slope = 1, intercept = 0) +
#   theme_classic() + 
#   ggtitle(paste("mean ortholog correlation\npaired t.test p <=", round_any(t_mean_cor$p.value, 
#                                                             accuracy = 0.000001, f = ceiling)))
# ggarrange(p_pct, p_avgCor, common.legend = TRUE,
#           nrow = 1, ncol = 2)
# dev.off()
#### Probably Archive: Separate matrices for each cluster ####
# 
# # possibly nice to keep as a control though, to check that these are the
# # same counts as the all-genes matrices, just with only specific rows
# 
# # Making single environment/cluster idx lists
# idxs <- finaldf |> filter(experiment == "LowN" & cer == 2 & par == 0) |> select(gene_name)
# write_delim(idxs, delim = "\n", file = "gene_ontology/20_LowN.txt")
# 
# # conserved plastic
# mat11 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_11.csv",
#                     delim = ";")
# mat22 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_22.csv",
#                     delim = ";")
# # conserved static
# mat00 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_00.csv",
#                     delim = ";")
# # diverged plastic
# mat12 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_12.csv",
#                     delim = ";")
# mat21 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_21.csv",
#                     delim = ";")
# # diverged, plastic in one species
# mat10 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_10.csv",
#                     delim = ";")
# mat01 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_01.csv",
#                     delim = ";")
# mat20 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_20.csv",
#                     delim = ";")
# mat02 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_02.csv",
#                     delim = ";")
# 
# ### Connectivity
# plotdf <- map2(list(mat11, mat22, mat12, mat21, mat10, mat01, mat20, mat02, mat00),
#                list("11", "22", "12", "21", "10", "01", "20", "02", "00"), \(x, y) {
#   output <- tibble(gene_name = x[1,-1],
#                    nConnections = colSums(x[,-1]),
#                    cluster = y)
#   return(output)
# }) |> purrr::reduce(.f = bind_rows)
# ggplot(plotdf, 
#        aes(x = cluster, y = nConnections)) +
#   geom_violin(aes(color = cluster))


