setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")
source("functions_for_figure_scripts.R")
load("data_files/FinalDataframe.RData")
load("data_files/Cleaned_Counts.RData")
load("data_files/Cleaned_Counts_Allele.RData")

experiment_name <- "LowPi" # change these for different environment-specific scripts
clusterdf_name <- "LowPi_2"
expdf <- filter(finaldf, experiment == "LowPi")
folder_name <- "LowPhosphorus"

# summary of groups in this experiment
table(expdf$group4)
table(expdf$group) |> sort(decreasing = TRUE)
sum(table(expdf$group))
length(unique(expdf$group))
code_order <- names(sort(table(expdf$group), decreasing = TRUE))

# upset plot of all divergence groups
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/upset", experiment_name, ".pdf"),
    width = 10, height = 3.25)
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
dev.off()

#### Average expression of genes in each group ####

#### average expression of all gene groups ####
ylims <- c(3.5, 9.5)
plotlist <- vector(mode = "list", length = length(code_order))
names(plotlist) <- code_order
for (g in code_order) {
  gene_idxs <- expdf |> filter(group == g) |> 
    select(gene_name) |> pull()
  plotlist[[g]] <- annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .plotlims = ylims,
                                             .experiment_name = experiment_name),
                                   top = paste0(length(gene_idxs), " / ", nrow(expdf), " genes"))
}
### Supplement: all gene groups
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/avg_expr_",
                  experiment_name, ".pdf"), width = 9, height = ceiling(length(plotlist)/3)*3)
ggarrange(plotlist = plotlist, ncol = 3, nrow = ceiling(length(plotlist)/3))
dev.off()

### main figure: top 3 groups in each category
# collecting sizes of all 16 groups (top 3 + other from each of the 4 groups)
sizedf <- map2(.x = rep(c("conserved", "diverged"), each = 2),
               .y = rep(c("conserved", "diverged"), times = 2), \(l, d) {
                 groups <- filter(expdf, level == l & dynamics == d) |> 
                   select(group) |> table() |> sort(decreasing = TRUE)
                 group_names <- groups |> names()
                 # data for prop area plots
                 output <- tibble(level = l, dynamics = d,
                                  n = sum(groups), n1 = as.numeric(groups[1]),
                                  n2 = as.numeric(groups[2]), n3 = as.numeric(groups[3]),
                                  groups = list(group_names[1:3]))
                 output$n_other <- output$n - (output$n1 + output$n2 + output$n3)
                 return(output)
               }) |> bind_rows()

sizedf
p_blank <- ggplot() + theme_void()

# line plots
# conserved
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", 
                  folder_name, "/conserved.pdf"), width = 5, height = 5)
ggarrange(plotlist = plotlist[unlist(sizedf$groups[sizedf$level == "conserved" &
                                                     sizedf$dynamics == "conserved"])], nrow = 2, ncol = 2)
dev.off()
# diverged level
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", 
                  folder_name, "/level.pdf"), width = 5, height = 5)
groups <- unlist(sizedf$groups[sizedf$level == "diverged" &
                                 sizedf$dynamics == "conserved"])
ggarrange(plotlist[[groups[1]]], plotlist[[groups[2]]],
          p_blank, plotlist[[groups[3]]], nrow = 2, ncol = 2)
dev.off()
# diverged dynamics
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", 
                  folder_name, "/dynamics.pdf"), width = 5, height = 5)
groups <- unlist(sizedf$groups[sizedf$level == "conserved" &
                                 sizedf$dynamics == "diverged"])
ggarrange(plotlist[[groups[1]]], p_blank,
          plotlist[[groups[2]]], plotlist[[groups[3]]], nrow = 2, ncol = 2)
dev.off()
# diverged level and dynamics
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", 
                  folder_name, "/level_and_dynamics.pdf"), width = 5, height = 5)
groups <- unlist(sizedf$groups[sizedf$level == "diverged" &
                                 sizedf$dynamics == "diverged"])
ggarrange(p_blank, plotlist[[groups[1]]],
          plotlist[[groups[2]]], plotlist[[groups[3]]], nrow = 2, ncol = 2)
dev.off()

sizedf[is.na(sizedf)] <- 0
# proportional area plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", 
                  folder_name, "/proportional_area.pdf"), width = 5, height = 5)
# size ordered by left to right, top to bottom + proportional area ordered CCW around unit circle = chaos
plotPropArea(x1 = c(sizedf$n2[sizedf$level == "diverged" & sizedf$dynamics == "conserved"],
                    sizedf$n1[sizedf$level == "diverged" & sizedf$dynamics == "conserved"],
                    sizedf$n_other[sizedf$level == "diverged" & sizedf$dynamics == "conserved"],
                    sizedf$n3[sizedf$level == "diverged" & sizedf$dynamics == "conserved"]),
             x2 = c(sizedf$n2[sizedf$level == "conserved" & sizedf$dynamics == "conserved"],
                    sizedf$n1[sizedf$level == "conserved" & sizedf$dynamics == "conserved"],
                    sizedf$n3[sizedf$level == "conserved" & sizedf$dynamics == "conserved"],
                    sizedf$n_other[sizedf$level == "conserved" & sizedf$dynamics == "conserved"]),
             x3 = c(sizedf$n_other[sizedf$level == "conserved" & sizedf$dynamics == "diverged"],
                    sizedf$n1[sizedf$level == "conserved" & sizedf$dynamics == "diverged"],
                    sizedf$n2[sizedf$level == "conserved" & sizedf$dynamics == "diverged"],
                    sizedf$n3[sizedf$level == "conserved" & sizedf$dynamics == "diverged"]),
             x4 = c(sizedf$n1[sizedf$level == "diverged" & sizedf$dynamics == "diverged"],
                    sizedf$n_other[sizedf$level == "diverged" & sizedf$dynamics == "diverged"],
                    sizedf$n2[sizedf$level == "diverged" & sizedf$dynamics == "diverged"],
                    sizedf$n3[sizedf$level == "diverged" & sizedf$dynamics == "diverged"]),
             .colors = levdyn_colordf$type[c(3, 1, 2, 4)])
dev.off()
# plot headings
sizedf

#### Divergence in level is preserved in the hybrid ####
plotdf <- filter(expdf, level != "low_expr")
# scatter plot of LFCs in parent and hybrid with dynamics divergence colored
# R squared for all (non-lowly expressed) genes
mod <- lm(effect_size_allele ~ effect_size_species, data = plotdf)
summary(mod)
# R squared for level-diverged genes
mod <- lm(effect_size_allele ~ effect_size_species, data = filter(plotdf, level == "diverged"))
summary(mod)
# plotting
p_scatter <- ggplot(plotdf,
                    aes(x = effect_size_species,
                        y = effect_size_allele)) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_vline(xintercept = c(-1, 1), color = "black", linetype = "dashed") +
  geom_point(aes(color = group4)) +
  scale_color_discrete(type = levdyn_colordf$type,
                       labels = levdyn_colordf$labels,
                       limits = levdyn_colordf$limits) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(color = "gene group") +
  xlab("log2 fold change, parents") +
  ylab("log2 fold change, hybrid") +
  ggtitle("") +
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  stat_cor(cor.coef.name = "R")

# accompanying barplot to show that most of these genes have conserved dynamics
p_bar <- ggplot(plotdf, aes(x = group4)) +
  geom_bar(aes(fill = group4)) +
  geom_text(stat='count', aes(label = after_stat(count)), vjust=-1) +
  scale_fill_discrete(type = levdyn_colordf$type,
                      labels = levdyn_colordf$labels,
                      limits = levdyn_colordf$limits) +
  scale_x_discrete(limits = levdyn_colordf$limits) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "right",
        plot.margin = margin(t = 1, r = 0.1, b = 0.1, l = 0.1, "cm")) +
  ggtitle("genes in each group") +
  ylim(c(0, 4100))

# level and dynamics summary
library(ggExtra)
ggarrange(p_scatter, p_bar, nrow = 1, ncol = 2,
          common.legend = TRUE, widths = c(2, 1), hjust = 1,
          legend = "bottom")

pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/l2fc.pdf"),
    width = 4, height = 4)
ggMarginal(p_scatter, groupColour = TRUE, groupFill = FALSE)
dev.off()

pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/bar.pdf"),
    width = 6, height = 4)
p_bar
dev.off()

#### Divergence in dynamics is not preserved in the hybrid ####
# collecting metrics
common_conditions <- intersect(info[info$experiment == experiment_name, "condition"],
                               info_allele[info_allele$experiment == experiment_name, "condition"])
parent_cor <- map(c(1:nrow(collapsed$cer)), \(i) {
  cor(collapsed$cer[i, common_conditions],
      collapsed$par[i, common_conditions], use = "pairwise.complete.obs")
}) |> unlist()
hybrid_cor <- map(c(1:nrow(collapsed_allele$cer)), \(i) {
  cor(collapsed_allele$cer[i, common_conditions],
      collapsed_allele$par[i, common_conditions], use = "pairwise.complete.obs")
}) |> unlist()

cordf <- tibble(gene_name = rownames(collapsed$cer),
                parent_cor,
                hybrid_cor)
plotdf <- left_join(select(expdf,
                           gene_name, effect_size_species,
                           effect_size_allele, pvalue_species,
                           pvalue_allele, cer, par, experiment,
                           dynamics, level, group4), cordf, by = "gene_name")

# testing metrics
# random conserved gene, increasing group
gene_idx <- plotdf |> 
  filter(level == "conserved" & cer == 1 & par == 1) |> 
  select(gene_name) |> pull() |> sample(1)
plotdf |> filter(gene_name == gene_idx) |> 
  select(effect_size_species, pvalue_species,
         effect_size_allele, pvalue_allele)
plotGenes(.gene_idxs = gene_idx, .quartet = TRUE, .experiment_name = experiment_name) 

# random level gene
gene_idx <- "YGR192C"
gene_idx <- plotdf |> 
  filter(level == "diverged" & dynamics == "conserved") |> 
  select(gene_name) |> pull() |> sample(1)
gene_idx
plotGenes(.gene_idxs = gene_idx, .quartet = TRUE, .experiment_name = experiment_name) 
# YIL136W YKL216W, YGL028C, YNL195C, pure trans
# YFL030W, YGR192C, YDR111C, YLR179C, YDR007W, YIL107C, YHR071W, YJL217W, cis
# YDR461C-A, YMR107W, cis then trans
# YDR530C, cis, then trans, then cis

# random dynamics gene
gene_idx <- plotdf |> 
  filter(level == "conserved" & dynamics == "diverged") |> 
  select(gene_name) |> pull() |> sample(1)
plotGenes(.gene_idxs = gene_idx, .quartet = TRUE, .experiment_name = experiment_name)

# quantifying all genes
# cor
p <- ggplot(plotdf, aes(x = parent_cor, y = hybrid_cor)) + 
  geom_point(aes(color = group4)) +
  scale_color_discrete(type = levdyn_colordf$type,
                       labels = levdyn_colordf$labels,
                       limits = levdyn_colordf$limits) +
  theme_classic() + 
  xlim(c(-1, 1)) + 
  ylim(c(-1, 1)) +
  theme(legend.position = "none") +
  xlab("parental correlation") +
  ylab("hybrid correlation")
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor.pdf"),
    width = 4, height = 4)
ggMarginal(p, groupColour = TRUE, groupFill = FALSE)
dev.off()

#### Separating cis and trans diverging level and dynamics genes ####

# Based on correlation and log2 fold change metrics,
# we separate out the genes
# with cis divergence from trans

### trans/level
# level divergers, cluster 1, cer up
# main plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/lfc_level_trans.pdf"),
    width = 3, height = 3)
ggplot(filter(plotdf, group4 == "diverged level, conserved dynamics"),
       aes(x = effect_size_species, y = effect_size_allele)) + 
  geom_point(color = levdyn_colordf$type[levdyn_colordf$limits == "diverged level, conserved dynamics"]) +
  theme_classic() + 
  xlim(c(-8, 8)) +
  ylim(c(-8, 8)) +
  theme(legend.position = "none") +
  xlab("log2 fold change, parents") +
  ylab("log2 fold change, hybrid") +
  geom_rect(xmin = 0, xmax = 5, ymin = -0.5, ymax = 0.5, color = "black", alpha = 0)
dev.off()
# corresponding avgExpr plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/lfc_level_trans_avgExpr.pdf"),
    width = 3, height = 3)
gene_idxs <- plotdf |> filter(cer == 1 & par == 1 &
                                group4 == "diverged level, conserved dynamics" &
                                effect_size_allele > -0.5 &
                                effect_size_allele < 0.5 &
                                effect_size_species > 1) |>
  select(gene_name) |> pull() |> sample(1)
annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste(gene_idxs, "example gene,\n level diverging in trans"))
dev.off()

### cis/level
# main plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/lfc_level_cis.pdf"),
    width = 3, height = 3)
ggplot(filter(plotdf, group4 == "diverged level, conserved dynamics"),
       aes(x = effect_size_species, y = effect_size_allele)) + 
  geom_point(color = levdyn_colordf$type[levdyn_colordf$limits == "diverged level, conserved dynamics"]) +
  theme_classic() + 
  theme(legend.position = "none") +
  xlab("log2 fold change, parents") +
  ylab("log2 fold change, hybrid") +
  geom_rect(xmin = 0, xmax = 10.5, ymin = 1, ymax = 10.5, color = "black", alpha = 0)
dev.off()
# corresponding avgExpr plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/lfc_level_cis_avgExpr.pdf"),
    width = 3, height = 3)
gene_idxs <- plotdf |> filter(cer == 1 & par == 1 &
                                group4 == "diverged level, conserved dynamics" &
                                effect_size_species > 0 &
                                effect_size_allele > 1) |>
  select(gene_name) |> pull()
annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste0(length(gene_idxs), " genes,\nlevel diverging in cis"))
dev.off()

### dynamics/cis
# dynamics divergers, 2-1
# main plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_dynamics_cis.pdf"),
    width = 3, height = 3)
ggplot(filter(plotdf, group4 == "conserved level, diverged dynamics"),
       aes(x = parent_cor, y = hybrid_cor)) + 
  geom_point(color = levdyn_colordf$type[levdyn_colordf$limits == "conserved level, diverged dynamics"]) +
  theme_classic() + 
  xlim(c(-1, 1)) + 
  ylim(c(-1, 1)) + 
  theme(legend.position = "none") +
  xlab("correlation, parents") +
  ylab("correlation, hybrid") +
  geom_rect(xmin = -1, xmax = 1, ymin = -1, ymax = -0.75, color = "black", alpha = 0)
dev.off()
# corresponding avgExpr plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_dynamics_cis_avgExpr.pdf"),
    width = 3, height = 3)
gene_idxs <- plotdf |> filter(cer == 1 & par == 2 &
                                group4 == "conserved level, diverged dynamics" &
                                hybrid_cor < -0.5) |>
  select(gene_name) |> pull()
length(gene_idxs)
gene_idx <- sample(gene_idxs, 1)
annotate_figure(plotGenes(.gene_idxs = gene_idx, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste(gene_idx, "example gene,\n dynamics diverging in cis"))
dev.off()

### dynamics/trans
# dynamics divergers, 2-1
# main plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_dynamics_trans.pdf"),
    width = 3, height = 3)
ggplot(filter(plotdf, group4 == "conserved level, diverged dynamics"),
       aes(x = parent_cor, y = hybrid_cor)) + 
  geom_point(color = levdyn_colordf$type[levdyn_colordf$limits == "conserved level, diverged dynamics"]) +
  theme_classic() + 
  xlim(c(-1, 1)) + 
  ylim(c(-1, 1)) + 
  theme(legend.position = "none") +
  xlab("correlation, parents") +
  ylab("correlation, hybrid") +
  geom_rect(xmin = -1, xmax = 1, ymin = 0, ymax = 1, color = "black", alpha = 0)
dev.off()
# corresponding avgExpr plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_dynamics_trans_avgExpr.pdf"),
    width = 3, height = 3)
gene_idxs <- plotdf |> filter(cer == 1 & par == 2 &
                                group4 == "conserved level, diverged dynamics" &
                                hybrid_cor > 0) |>
  select(gene_name) |> pull()
annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste0(length(gene_idxs), " genes,\n dynamics diverging in trans"))
dev.off()

### level and dynamics divergers, cis x trans
# main plot, correlation, trans
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_leveldyn_trans.pdf"),
    width = 3, height = 3)
ggplot(filter(plotdf, group4 == "diverged level and dynamics"),
       aes(x = parent_cor, y = hybrid_cor)) + 
  geom_point(color = levdyn_colordf$type[levdyn_colordf$limits == "diverged level and dynamics"]) +
  theme_classic() + 
  theme(legend.position = "none") +
  xlab("correlation, parents") +
  ylab("correlation, hybrid") +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  geom_rect(xmin = -1, xmax = 1, ymin = 0, ymax = 1, color = "black", alpha = 0)
dev.off()
# corresponding avgExpr plot
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_leveldyn_cistrans_avgExpr.pdf"),
    width = 3, height = 3)
# these guys are doing too many different things for avgExpr to make sense
gene_idxs <- plotdf |> filter(group4 == "diverged level and dynamics" &
                                hybrid_cor > 0) |>
  select(gene_name) |> pull() |> sample(1)
annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste(gene_idxs, "example gene,\nlevel diverging in cis,\ndynamics diverging in trans"))
dev.off()
# also looking for a pure trans lev dyn gene
gene_idxs <- plotdf |> filter(group4 == "diverged level and dynamics" &
                                hybrid_cor > 0) |>
  select(gene_name) |> pull() |> sample(1)
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_leveldyn_trans_avgExpr.pdf"),
    width = 3, height = 3)
annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste(gene_idxs, "example gene,\nlevel and dynamics \ndiverging in trans"))
dev.off()
# main plot, correlation, outliers
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_leveldyn_outliers.pdf"),
    width = 3, height = 3)
ggplot(filter(plotdf, group4 == "diverged level and dynamics"),
       aes(x = parent_cor, y = hybrid_cor)) + 
  geom_point(color = levdyn_colordf$type[levdyn_colordf$limits == "diverged level and dynamics"]) +
  theme_classic() + 
  theme(legend.position = "none") +
  xlab("correlation, parents") +
  ylab("correlation, hybrid") +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  geom_rect(xmin = -1, xmax = 1, ymin = -1, ymax = -0.75, color = "black", alpha = 0)
dev.off()
# corresponding avgExpr plot
gene_idxs <- plotdf |> filter(group4 == "diverged level and dynamics" &
                                hybrid_cor < -0.75) |>
  select(gene_name) |> pull() |> sample(1)
pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/", folder_name, "/cor_leveldyn_outlier1_avgExpr.pdf"),
    width = 3, height = 3)
annotate_figure(plotGenes(.gene_idxs = gene_idxs, .quartet = TRUE, .experiment_name = experiment_name),
                top = paste(gene_idxs, "example gene\nlevel and dynamics\ndiverging in cis"))
dev.off()
