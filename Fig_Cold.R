setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2025/")
source("functions_for_figure_scripts.R")
load("data_files/FinalDataframe3Disp.RData")
load("data_files/Cleaned_Counts.RData")
load("data_files/Cleaned_Counts_Allele.RData")

experiment_name <- "Cold" # change these for different environment-specific scripts
clusterdf_name <- "Cold_3"
expdf <- filter(finaldf, experiment == "Cold")
folder_name <- "ColdStress"

# summary of groups in this experiment
table(expdf$group) |> sort(decreasing = TRUE)
sum(table(expdf$group))
length(unique(expdf$group))
code_order <- names(sort(table(expdf$group), decreasing = TRUE))

# fisher exact test for un-enrichment of dynamics divergence 
# among level divergent set
plotdf <- expdf |> filter(level != "low_expr")
n_leveldiv <- sum(plotdf$level == "diverged")
n_levelcons <- sum(plotdf$level == "conserved")
n_dyndiv <- sum(plotdf$dynamics == "diverged")
n_dyncons <- sum(plotdf$dynamics == "conserved")
n_dynandlevdiv <- sum(plotdf$dynamics == "diverged" &
                        plotdf$level == "diverged")
fisher.test(matrix(c(n_dynandlevdiv, n_dyndiv - n_dynandlevdiv,
                     n_leveldiv - n_dynandlevdiv, 
                     n_dyncons - (n_leveldiv - n_dynandlevdiv)), 
                   nrow = 2))
# barplot of all divergence groups
ggplot(expdf, aes(x = group)) + 
  geom_bar(aes(fill = group4)) +
  scale_x_discrete(breaks = code_order, limits = code_order) +
  scale_fill_discrete(type = levdyn_colordf$type,
                      labels = levdyn_colordf$labels,
                      limits = levdyn_colordf$limits) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  ylim(c(0, 2750))

#### Gene count tables ####
# Biased genes, can't determine level divergence:
expdf |> filter(level == "biased") |> nrow()

# Lowly expressed genes in this environment:
setdiff(unique(finaldf$gene_name), expdf$gene_name) |> length()

# Scer gene counts:
expdf |> select(cer) |> table()
expdf |> select(par) |> table()

# Conserved genes:
expdf |> filter(dynamics != "diverged" & 
                  level != "diverged") |> 
  select(cer, par) |> table()

# Diverged level:
expdf |> filter(dynamics != "diverged" & 
                  level == "diverged") |> 
  mutate(direction = if_else(effect_size_species > 0,
                             true = "up in Scer",
                             false = "up in Spar")) |> 
  select(cer, par, direction) |> table()

# Diverged dynamics:
expdf |> filter(dynamics == "diverged" & 
                  level != "diverged") |> 
  mutate(isZero = if_else(cer == 0 | par == 0,
                          true = if_else(cer == 0,
                                         true = "Scer0",
                                         false = "Spar0"),
                          false = "neither0")) |> 
  select(cer, par, isZero) |> table()

# Diverged level and dynamics:
expdf |> filter(dynamics == "diverged" & 
                  level == "diverged") |> 
  mutate(direction = if_else(effect_size_species > 0,
                             true = "up in Scer",
                             false = "up in Spar")) |> 
  select(cer, par, direction) |> table()

#### Average expression of genes in each group ####
# Conserved genes:
ylims <- c(3, 9)
cluster_names <- expdf |> filter(dynamics != "diverged" & 
                                   level != "diverged") |> 
  select(cer, par) |> pull() |> unique()
plotlist <- vector(mode = "list", length = length(cluster_names))
names(plotlist) <- cluster_names
for (clust in cluster_names) {
  gene_idxs <- expdf |> filter(dynamics != "diverged" & 
                                 level != "diverged" &
                                 cer == clust &
                                 par == clust) |> 
    select(gene_name) |> pull()
  plotlist[[as.character(clust)]] <- annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = FALSE,
                                                               .experiment_name = experiment_name,
                                                               .plotlims = ylims) + 
                                                       theme(plot.margin = unit(c(0,0.25,0,0), "cm")),
                                                     top = paste0(length(gene_idxs), " genes"))
}
### Supplement: top 3 gene groups
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/", 
                  folder_name, "/conserved_top3.pdf"), width = 2, height = 5)
ggarrange(plotlist = plotlist[1:3], nrow = 3, ncol = 1)
dev.off()
### Supplement: all gene groups
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/", 
                  folder_name, "/conserved_all.pdf"), width = 2, height = (5/3)*length(plotlist))
ggarrange(plotlist = plotlist, nrow = length(plotlist), ncol = 1)
dev.off()

# Diverged level:
griddf <- expdf |> filter(dynamics != "diverged" & 
                            level == "diverged") |> 
  mutate(lfc_sign = sign(effect_size_species)) |> 
  group_by(cer, par, lfc_sign) |> summarise(n = n()) |> 
  arrange(desc(n))
plotlist <- vector(mode = "list", length = 0)
for (i in c(1:nrow(griddf))) {
  clust <- griddf$cer[i]
  lfcsign <- griddf$lfc_sign[i]
  gene_idxs <- expdf |> filter(dynamics != "diverged" & 
                                 level == "diverged" &
                                 cer == clust &
                                 par == clust &
                                 sign(effect_size_species) == lfcsign) |> 
    select(gene_name) |> pull()
  plotlist[[paste0(clust, lfcsign)]] <- annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = FALSE,
                                                                  .experiment_name = experiment_name,
                                                                  .plotlims = ylims) + 
                                                          theme(plot.margin = unit(c(0,0.25,0,0), "cm")),
                                                        top = paste0(length(gene_idxs), " genes"))
}

### Supplement: top 3 gene groups
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/", 
                  folder_name, "/level_top3.pdf"), width = 2, height = 5)
ggarrange(plotlist = plotlist[1:3], nrow = 3, ncol = 1)
dev.off()
### Supplement: all gene groups
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/", 
                  folder_name, "/level_all.pdf"), width = 2, height = (5/3)*length(plotlist))
ggarrange(plotlist = plotlist, nrow = length(plotlist), ncol = 1)
dev.off()

# Diverged dynamics:
griddf <- expdf |> filter(dynamics == "diverged" & 
                            level != "diverged") |> 
  group_by(cer, par) |> 
  summarise(n = n()) |> 
  arrange(desc(n))
plotlist <- vector(mode = "list", length = 0)
for (i in c(1:nrow(griddf))) {
  clust_cer <- griddf$cer[i]
  clust_par <- griddf$par[i]
  gene_idxs <- expdf |> filter(dynamics == "diverged" & 
                                 level != "diverged" &
                                 cer == clust_cer &
                                 par == clust_par) |> 
    select(gene_name) |> pull()
  plotlist[[paste0(clust_cer, clust_par)]] <- annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = FALSE,
                                                                        .experiment_name = experiment_name,
                                                                        .plotlims = ylims) + 
                                                                theme(plot.margin = unit(c(0,0.25,0,0), "cm")),
                                                              top = paste0(length(gene_idxs), " genes"))
}

### Supplement: top 3 gene groups
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/", 
                  folder_name, "/dynamics_top3.pdf"), width = 2, height = 5)
ggarrange(plotlist = plotlist[1:3], nrow = 3, ncol = 1)
dev.off()
### Supplement: all gene groups
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/", 
                  folder_name, "/dynamics_all.pdf"), width = 2, height = (5/3)*length(plotlist))
ggarrange(plotlist = plotlist, nrow = length(plotlist), ncol = 1)
dev.off()

# Diverged level and dynamics:
griddf <- expdf |> filter(dynamics == "diverged" & 
                            level == "diverged") |> 
  mutate(lfc_sign = sign(effect_size_species)) |> 
  group_by(cer, par, lfc_sign) |> 
  summarise(n = n()) |> 
  arrange(desc(n))
plotlist <- vector(mode = "list", length = 0)
for (i in c(1:nrow(griddf))) {
  clust_cer <- griddf$cer[i]
  clust_par <- griddf$par[i]
  lfcsign <- griddf$lfc_sign[i]
  gene_idxs <- expdf |> filter(dynamics == "diverged" & 
                                 level == "diverged" &
                                 cer == clust_cer &
                                 par == clust_par &
                                 sign(effect_size_species) == lfcsign) |> 
    select(gene_name) |> pull()
  plotlist[[paste0(clust_cer, clust_par, lfcsign)]] <- annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = FALSE,
                                                                                 .experiment_name = experiment_name,
                                                                                 .plotlims = ylims) + 
                                                                         theme(plot.margin = unit(c(0,0.25,0,0), "cm")),
                                                                       top = paste0(length(gene_idxs), " genes"))
}

### Supplement: top 3 gene groups
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/", 
                  folder_name, "/level_and_dynamics_top3.pdf"), width = 2, height = 5)
ggarrange(plotlist = plotlist[1:3], nrow = 3, ncol = 1)
dev.off()
### Supplement: all gene groups
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/", 
                  folder_name, "/level_and_dynamics_all.pdf"), width = 2, height = (5/3)*length(plotlist))
ggarrange(plotlist = plotlist, nrow = length(plotlist), ncol = 1)
dev.off()

### Proportional area plot for Main figure
prop_table <- expdf |> mutate(dynamics_diverged = dynamics == "diverged",
                              level_diverged = level == "diverged") |> 
  select(dynamics_diverged, level_diverged) |> table()
prop_table
sum(prop_table)
prop_table/sum(prop_table)
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/", 
                  folder_name, "/prop_area.pdf"), width = 4, height = 4)
plotPropAreaSingle(.counts = c(prop_table[1, 2], prop_table[1, 1],
                               prop_table[2, 1], prop_table[2, 2]),
                   levdyn_colordf$type[c(3, 1, 2, 4)],
                   .text_size = 10, .text_bounds = 20)
dev.off()
