sapply(c("tidyr", "dplyr", "purrr", "ggplot2", "bipartite", "ggpubr", "ggExtra"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")
source("functions_for_figure_scripts.R")

### Investigating the genetic basis of module expression 
# pattern conservation and divergence

# Thus far we have seen that co-expression conservation
# does not predict conservation of expression level (what Krieger 2020 saw)
# But it DOES predict conservation of expression pattern (Krieger 2020 didn't look at this)

# Additionally, we have seen that expression pattern correlation
# tends to improve MORE in CCM genes than divergent genes in the hybrid,
# indicating that CCM genes tend to be diverging more in trans-regulation than divergent genes

#### Supplemental figure: all hybrid expression profiles ####

# merge 10
plotdf <- moduledf_allele10 |> filter(block_size >= 10)
plotlist_ccm <- vector(mode = "list", length = sum(plotdf$coexpressed == "conserved co-expressed"))
plotlist_div <- vector(mode = "list", length = sum(plotdf$coexpressed == "diverged co-expressed"))
plotlist_cer <- vector(mode = "list", length = sum(plotdf$coexpressed == "S. cerevisiae co-expressed"))
plotlist_par <- vector(mode = "list", length = sum(plotdf$coexpressed == "S. paradoxus co-expressed"))
plotlist_nev <- vector(mode = "list", length = sum(plotdf$coexpressed == "never co-expressed"))
plotlist <- list("conserved co-expressed" = plotlist_ccm,
                 "diverged co-expressed" = plotlist_div,
                 "S. cerevisiae co-expressed" = plotlist_cer,
                 "S. paradoxus co-expressed" = plotlist_par,
                 "never co-expressed" = plotlist_nev)
for (idx in 1:nrow(plotdf)) {
  m <- plotdf$module_name[idx]
  gene_idxs <- module_genedf10 |> filter(module_name == m) |> select(gene_name) |> pull()
  modtype <- plotdf$coexpressed[idx]
  p <- plotExpressionProfilePair(counts_all2_allele$cer[,gene_idxs],
                                 counts_all2_allele$par[,gene_idxs],
                                 info,
                                 info,
                                 .method = "line", .show_points = TRUE,
                                 .normalization = "log2")
  p <- annotate_figure(p, top = paste0(m, 
                                       " dcor_CC=", round(plotdf$dcor_CC[idx], 3),
                                       " dcor_HAP4=", round(plotdf$dcor_HAP4[idx], 3),
                                       " dcor_LowN=", round(plotdf$dcor_LowN[idx], 3),
                                       " dcor_LowPi=", round(plotdf$dcor_LowPi[idx], 3)))
  # arranging plots in order of decreasing AvgCor
  p_idx <- which(sort(plotdf$dcor[plotdf$coexpressed == modtype], decreasing = TRUE) == plotdf$dcor[idx])
  plotlist[[modtype]][[p_idx]] <- p
}
# random
plotdf <- random_moduledf_allele10 |> filter(block_size >= 10)
plotlist_Random <- vector(mode = "list", length = nrow(plotdf))
for (idx in 1:nrow(plotdf)) {
  m <- plotdf$module_name[idx]
  gene_idxs <- random_module_genedf10 |> filter(module_name == m) |> select(gene_name) |> pull()
  p <- plotExpressionProfilePair(counts_all2_allele$cer[,gene_idxs],
                                 counts_all2_allele$par[,gene_idxs],
                                 info,
                                 info,
                                 .method = "line", .show_points = TRUE,
                                 .normalization = "log2")
  p <- annotate_figure(p, top = paste0(m, 
                                       " dcor_CC=", round(plotdf$dcor_CC[idx], 3),
                                       " dcor_HAP4=", round(plotdf$dcor_HAP4[idx], 3),
                                       " dcor_LowN=", round(plotdf$dcor_LowN[idx], 3),
                                       " dcor_LowPi=", round(plotdf$dcor_LowPi[idx], 3)))
  p_idx <- which(sort(plotdf$dcor, decreasing = TRUE) == plotdf$dcor[idx])
  plotlist_Random[[p_idx]] <- p
}
# making pdfs
# CCMs
p_CCM <- ggarrange(plotlist = plotlist[["conserved co-expressed"]], 
                   nrow = length(plotlist[["conserved co-expressed"]]), 
                   ncol = 1, 
                   common.legend = TRUE, 
                   legend = "right")
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_hybrid_expr_profiles_CCM10.pdf",
    width = 12, height = length(plotlist_ccm)*4)
p_CCM
dev.off()
# Divergent
p_Divergent <- ggarrange(plotlist = plotlist[["diverged co-expressed"]], 
                         nrow = length(plotlist[["diverged co-expressed"]]), 
                         ncol = 1, 
                         common.legend = TRUE, 
                         legend = "right")
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_hybrid_expr_profiles_Divergent10.pdf",
    width = 12, height = length(plotlist[["diverged co-expressed"]])*4)
p_Divergent
dev.off()

# S. cerevisiae co-expressed
p_Scer <- ggarrange(plotlist = plotlist[["S. cerevisiae co-expressed"]], 
                    nrow = length(plotlist[["S. cerevisiae co-expressed"]]), 
                    ncol = 1, 
                    common.legend = TRUE, 
                    legend = "right")
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_hybrid_expr_profiles_Scer10.pdf",
    width = 12, height = length(plotlist[["S. cerevisiae co-expressed"]])*4)
p_Scer
dev.off()

# S. paradoxus co-expressed
p_Spar <- ggarrange(plotlist = plotlist[["S. paradoxus co-expressed"]], 
                    nrow = length(plotlist[["S. paradoxus co-expressed"]]), 
                    ncol = 1, 
                    common.legend = TRUE, 
                    legend = "right")
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_hybrid_expr_profiles_Spar10.pdf",
    width = 12, height = length(plotlist[["S. paradoxus co-expressed"]])*4)
p_Scer
dev.off()

# Never co-expressed
p_Never <- ggarrange(plotlist = plotlist[["never co-expressed"]], 
                     nrow = length(plotlist[["never co-expressed"]]), 
                     ncol = 1, 
                     common.legend = TRUE, 
                     legend = "right")
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_hybrid_expr_profiles_Never10.pdf",
    width = 12, height = length(plotlist[["never co-expressed"]])*4)
p_Never
dev.off()

# Random
p_Random <- ggarrange(plotlist = plotlist_Random, 
                      nrow = length(plotlist_Random), 
                      ncol = 1, 
                      common.legend = TRUE, 
                      legend = "right")
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/all_hybrid_expr_profiles_Random10.pdf",
    width = 12, height = length(plotlist_Random)*4)
p_Random
dev.off()

#### Cis/trans basis of expression divergence - distribution of all modules ####

# We have seen that expression pattern correlation
# tends to improve MORE in CCM genes than divergent genes in the hybrid
# is this generally true (like every CCM improves correlation in
# the hybrid more than every divergent module does), 
# or are a few CCMs really jumping up their correlation?

# Species vs Allele dcor for each module

# regressing avg correlation between species with avg correlation between
# alleles where each point is a module
plotdf1 <- moduledf25 %>% 
  filter(block_size >= 30 & is_CCM) |> 
  select(module_name, avgCor_CC,
         avgCor_LowN, avgCor_LowPi,
         avgCor_HAP4, avgCor_Heat,
         avgCor_Cold, is_CCM,
         cer_color, par_color,
         coexpressed, CCM_color) %>% 
  pivot_longer(cols = c("avgCor_CC", 
                        "avgCor_LowN", 
                        "avgCor_LowPi", 
                        "avgCor_HAP4",
                        "avgCor_Heat",
                        "avgCor_Cold"), 
               names_to = "experiment",
               values_to = "avgCor")
plotdf1$allele_or_species <- "species"
plotdf2 <- moduledf_allele25  %>% 
  filter(block_size >= 25 & is_CCM) |>
  select(module_name, avgCor_CC,
         avgCor_LowN, avgCor_LowPi,
         avgCor_HAP4, avgCor_Heat,
         avgCor_Cold, is_CCM,
         cer_color, par_color,
         coexpressed, CCM_color) %>% 
  pivot_longer(cols = c("avgCor_CC", 
                        "avgCor_LowN", 
                        "avgCor_LowPi", 
                        "avgCor_HAP4",
                        "avgCor_Heat",
                        "avgCor_Cold"), 
               names_to = "experiment",
               values_to = "avgCor")
plotdf2$allele_or_species <- "allele"

plotdf <- bind_rows(plotdf1, plotdf2) %>% pivot_wider(id_cols = c("module_name", "experiment", "is_CCM", "cer_color", "par_color", "coexpressed", "CCM_color"),
                                                      names_from = allele_or_species,
                                                      values_from = avgCor)
table(plotdf$module_name) %>% table() # should all have 6 entries, 1 for each experiment
# adding random
rplotdf1 <- random_moduledf25 %>%
  filter(block_size >= 30) |> 
  mutate(CCM_color = "black") |> 
  select(module_name, avgCor_CC,
         avgCor_LowN, avgCor_LowPi,
         avgCor_HAP4, avgCor_Heat,
         avgCor_Cold, is_CCM,
         cer_color, par_color,
         coexpressed, CCM_color) %>% 
  pivot_longer(cols = c("avgCor_CC", 
                        "avgCor_LowN", 
                        "avgCor_LowPi", 
                        "avgCor_HAP4",
                        "avgCor_Heat",
                        "avgCor_Cold"), 
               names_to = "experiment",
               values_to = "avgCor")
rplotdf1$allele_or_species <- "species"
rplotdf2 <- random_moduledf_allele25  %>% 
  filter(block_size >= 30) |>  
  mutate(CCM_color = "black") |> 
  select(module_name, avgCor_CC,
         avgCor_LowN, avgCor_LowPi,
         avgCor_HAP4, avgCor_Heat,
         avgCor_Cold, is_CCM,
         cer_color, par_color,
         CCM_color,
         coexpressed) %>% 
  pivot_longer(cols = c("avgCor_CC", 
                        "avgCor_LowN", 
                        "avgCor_LowPi", 
                        "avgCor_HAP4",
                        "avgCor_Heat",
                        "avgCor_Cold"), 
               names_to = "experiment",
               values_to = "avgCor")
rplotdf2$allele_or_species <- "allele"

rplotdf <- bind_rows(rplotdf1, rplotdf2) %>% pivot_wider(id_cols = c("module_name", "experiment", "is_CCM", "cer_color", "par_color", "coexpressed", "CCM_color"),
                                                         names_from = allele_or_species,
                                                         values_from = avgCor)
table(rplotdf$module_name) %>% table()
plotdf <- bind_rows(plotdf, rplotdf)
plotdf$experiment <- gsub("avgCor_", "", plotdf$experiment)
p <- plotdf |> 
  arrange(desc(coexpressed)) |> 
  ggplot(aes(x = species, y = allele)) + 
  # to add CCM colors (breaks ggmarginal):
  geom_point(data = filter(plotdf, is_CCM),
             aes(x = species, y = allele),
             color = plotdf$CCM_color[plotdf$is_CCM],
             size = 3) +
  geom_point(aes(color = coexpressed, shape = experiment)) +
  # # if you wanted to highlight one particular module:
  # geom_point(data = filter(plotdf, experiment == "LowPi" & module_name == "j"),
  #            aes(x = species, y = allele), color = "red", shape = "+", size = 4) +
  # add this line to see which modules are which (but too ugly for final plot):
  #geom_text(data = plotdf, aes(x = species, y = allele, label = module_name), nudge_x = .05) +
  scale_color_discrete(type = c("gold",
                                "black"),
                       limits = c("conserved co-expressed",
                                  "random"),
                       labels = c("conserved modules",
                                  "random modules")) +
  scale_shape_discrete(limits = c("CC", "HAP4", "LowN", "LowPi", "Heat", "Cold"),
                       labels = c("Urea Shock", "Saturated Growth",
                                  "Low Nitrogen", "Low Phosphorus",
                                  "Heat Shock", "Cold Shock")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 15),
        legend.position = "left") +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  xlab("Module correlation \nbetween parental species") +
  ylab("Module correlation \nbetween hybrid alleles") +
  xlim(c(0, 1)) +
  ylim(c(0, 1))
pdf("../../aligning_the_molecular_phenotype/paper_figures/CisTrans/allale_vs_species.pdf",
    width = 7, height = 5)
p
# ggMarginal(p, groupColour = TRUE, groupFill = TRUE)
dev.off()

#### Supplemental figure: allele vs species for all 5 module co-expression classes and all merge levels ####
makeCisTransPlot <- function(.modsize = "25") {
  if (.modsize == "00") {
    .mdf <- moduledf00
    .mdf_allele <- moduledf_allele00
    .rmdf <- random_moduledf00
    .rmdf_allele <- random_moduledf_allele00
  }
  if (.modsize == "10") {
    .mdf <- moduledf10
    .mdf_allele <- moduledf_allele10
    .rmdf <- random_moduledf10
    .rmdf_allele <- random_moduledf_allele10
  }
  if (.modsize == "25") {
    .mdf <- moduledf25
    .mdf_allele <- moduledf_allele25
    .rmdf <- random_moduledf25
    .rmdf_allele <- random_moduledf_allele25
  }
  if (.modsize == "35") {
    .mdf <- moduledf35
    .mdf_allele <- moduledf_allele35
    .rmdf <- random_moduledf35
    .rmdf_allele <- random_moduledf_allele35
  }
  plotdf1 <- .mdf %>% 
    filter(block_size >= 25) |> 
    select(module_name, avgCor_CC,
           avgCor_LowN, avgCor_LowPi,
           avgCor_HAP4, is_CCM,
           cer_color, par_color,
           coexpressed, CCM_color) %>% 
    pivot_longer(cols = c("avgCor_CC", 
                          "avgCor_LowN", 
                          "avgCor_LowPi", 
                          "avgCor_HAP4"), 
                 names_to = "experiment",
                 values_to = "avgCor")
  plotdf1$allele_or_species <- "species"
  plotdf2 <- .mdf_allele  %>% 
    filter(block_size >= 25) |>
    select(module_name, avgCor_CC,
           avgCor_LowN, avgCor_LowPi,
           avgCor_HAP4, is_CCM,
           cer_color, par_color,
           coexpressed, CCM_color) %>% 
    pivot_longer(cols = c("avgCor_CC", 
                          "avgCor_LowN", 
                          "avgCor_LowPi", 
                          "avgCor_HAP4"), 
                 names_to = "experiment",
                 values_to = "avgCor")
  plotdf2$allele_or_species <- "allele"
  
  plotdf <- bind_rows(plotdf1, plotdf2) %>% pivot_wider(id_cols = c("module_name", "experiment", "is_CCM", "cer_color", "par_color", "coexpressed", "CCM_color"),
                                                        names_from = allele_or_species,
                                                        values_from = avgCor)
  # adding random
  rplotdf1 <- .rmdf %>%
    filter(block_size >= 25) |> 
    mutate(CCM_color = "purple1") |> 
    select(module_name, avgCor_CC,
           avgCor_LowN, avgCor_LowPi,
           avgCor_HAP4, is_CCM,
           cer_color, par_color,
           coexpressed, CCM_color) %>% 
    pivot_longer(cols = c("avgCor_CC", 
                          "avgCor_LowN", 
                          "avgCor_LowPi", 
                          "avgCor_HAP4"), 
                 names_to = "experiment",
                 values_to = "avgCor")
  rplotdf1$allele_or_species <- "species"
  rplotdf2 <- .rmdf_allele  %>% 
    filter(block_size >= 25) |>  
    mutate(CCM_color = "purple1") |> 
    select(module_name, avgCor_CC,
           avgCor_LowN, avgCor_LowPi,
           avgCor_HAP4, is_CCM,
           cer_color, par_color,
           CCM_color,
           coexpressed) %>% 
    pivot_longer(cols = c("avgCor_CC", 
                          "avgCor_LowN", 
                          "avgCor_LowPi", 
                          "avgCor_HAP4"), 
                 names_to = "experiment",
                 values_to = "avgCor")
  rplotdf2$allele_or_species <- "allele"
  
  rplotdf <- bind_rows(rplotdf1, rplotdf2) %>% pivot_wider(id_cols = c("module_name", "experiment", "is_CCM", "cer_color", "par_color", "coexpressed", "CCM_color"),
                                                           names_from = allele_or_species,
                                                           values_from = avgCor)
  plotdf <- bind_rows(plotdf, rplotdf)
  plotdf$experiment <- gsub("avgCor_", "", plotdf$experiment)
  p <- plotdf |> 
    arrange(desc(coexpressed)) |> 
    ggplot(aes(x = species, y = allele)) + 
    # to add CCM colors (breaks ggmarginal):
    # geom_point(data = filter(plotdf, is_CCM), 
    #            aes(x = species, y = allele), 
    #            color = plotdf$CCM_color[plotdf$is_CCM],
    #            size = 3) +
    geom_point(aes(color = coexpressed, shape = experiment)) +
    # # if you wanted to highlight one particular module:
    # geom_point(data = filter(plotdf, experiment == "LowPi" & module_name == "j"),
    #            aes(x = species, y = allele), color = "red", shape = "+", size = 4) +
    # add this line to see which modules are which (but too ugly for final plot):
    #geom_text(data = plotdf, aes(x = species, y = allele, label = module_name), nudge_x = .05) +
    scale_color_discrete(type = c("gold",
                                  "mediumseagreen",
                                  "orange1",
                                  "blue2",
                                  "grey",
                                  "purple1"),
                         limits = c("conserved co-expressed",
                                    "diverged co-expressed",
                                    "S. cerevisiae co-expressed",
                                    "S. paradoxus co-expressed",
                                    "never co-expressed",
                                    "random")) +
    scale_shape_discrete(limits = c("CC", "HAP4", "LowN", "LowPi"),
                         labels = c("Urea Shock", "Saturated Growth",
                                    "Low Nitrogen", "Low Phosphorus")) +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = "left") +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    xlab("Module correlation \nbetween parental species") +
    ylab("Module correlation \nbetween hybrid alleles") +
    xlim(c(0, 1)) +
    ylim(c(0, 1))
  
  return(p)
}

# merge 00
p <- makeCisTransPlot("00")
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/allele_vs_species_00.pdf",
    width = 7, height = 5)
ggMarginal(p, groupColour = TRUE, groupFill = TRUE)
dev.off()

# merge 10
p <- makeCisTransPlot("10")
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/allele_vs_species_10.pdf",
    width = 7, height = 5)
ggMarginal(p, groupColour = TRUE, groupFill = TRUE)
dev.off()

# merge 25
p <- makeCisTransPlot("25")
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/allele_vs_species_25.pdf",
    width = 7, height = 5)
ggMarginal(p, groupColour = TRUE, groupFill = TRUE)
dev.off()

# merge 35
p <- makeCisTransPlot("35")
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/allele_vs_species_35.pdf",
    width = 7, height = 5)
ggMarginal(p, groupColour = TRUE, groupFill = TRUE)
dev.off()

#### Possibly supplemental figure: overall dcor ####

# motivation: dcor of modules in environments they don't vary in may
# be low despite the module having similar expression between species
# looking at overall dcor may reduce this possibility, as all modules
# vary in their expression in at least one environment

plotdf1 <- moduledf10 %>% 
  filter(block_size >= 25) |> 
  select(module_name, dcor, is_CCM,
         cer_color, par_color,
         coexpressed, CCM_color)
plotdf1$allele_or_species <- "species"
plotdf2 <- moduledf_allele10  %>% 
  filter(block_size >= 25) |>
  select(module_name, dcor, is_CCM,
         cer_color, par_color,
         coexpressed, CCM_color)
plotdf2$allele_or_species <- "allele"

plotdf <- bind_rows(plotdf1, plotdf2) %>% pivot_wider(id_cols = c("module_name", "is_CCM", "cer_color", "par_color", "coexpressed", "CCM_color"),
                                                      names_from = allele_or_species,
                                                      values_from = dcor)
table(plotdf$module_name) %>% table() # should all have 1 entrie, 1 per module
# adding random
rplotdf1 <- random_moduledf10 %>%
  filter(block_size >= 25) |> 
  mutate(CCM_color = "purple1") |> 
  select(module_name, dcor, is_CCM,
         cer_color, par_color,
         coexpressed, CCM_color)
rplotdf1$allele_or_species <- "species"
rplotdf2 <- random_moduledf_allele10  %>% 
  filter(block_size >= 25) |>  
  mutate(CCM_color = "purple1") |> 
  select(module_name, dcor, is_CCM,
         cer_color, par_color,
         CCM_color,
         coexpressed)
rplotdf2$allele_or_species <- "allele"

rplotdf <- bind_rows(rplotdf1, rplotdf2) %>% pivot_wider(id_cols = c("module_name", "is_CCM", "cer_color", "par_color", "coexpressed", "CCM_color"),
                                                         names_from = allele_or_species,
                                                         values_from = dcor)
table(rplotdf$module_name) %>% table()
plotdf <- bind_rows(plotdf, rplotdf)
p <- plotdf |> 
  arrange(desc(coexpressed)) |> 
  ggplot(aes(x = species, y = allele)) + 
  # to add CCM colors (breaks ggmarginal):
  geom_point(data = filter(plotdf, is_CCM),
             aes(x = species, y = allele),
             color = plotdf$CCM_color[plotdf$is_CCM],
             size = 3) +
  geom_point(aes(color = coexpressed)) +
  # # if you wanted to highlight one particular module:
  # geom_point(data = filter(plotdf, experiment == "LowPi" & module_name == "j"),
  #            aes(x = species, y = allele), color = "red", shape = "+", size = 4) +
  # add this line to see which modules are which (but too ugly for final plot):
  #geom_text(data = plotdf, aes(x = species, y = allele, label = module_name), nudge_x = .05) +
  scale_color_discrete(type = c("gold",
                                "mediumseagreen",
                                "orange1",
                                "blue2",
                                "grey",
                                "purple1"),
                       limits = c("conserved co-expressed",
                                  "diverged co-expressed",
                                  "S. cerevisiae co-expressed",
                                  "S. paradoxus co-expressed",
                                  "never co-expressed",
                                  "random")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "left") +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  xlab("Module correlation \nbetween parental species") +
  ylab("Module correlation \nbetween hybrid alleles") +
  xlim(c(0, 1)) +
  ylim(c(0, 1))
p
# ggMarginal(p, groupColour = TRUE, groupFill = TRUE)

# conclusions: Using overall dcor doesn't change conclusions.
# tan, pink, and black are still the most divergent modules in the hybrid

#### Hybrid module correlation improvement example ####

# module h, black, cis regulation (as much as we can find, which isn't much)
gene_idxs <- module_genedf10 |> filter(module_name == "h") |> select(gene_name) |> pull()
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTrans/h_parents.pdf",
        width = 12, height = 3)
plotExpressionProfilePair(counts_all2$cer[,gene_idxs],
                          counts_all2$par[,gene_idxs],
                          info,
                          info,
                          .method = "line",
                          .normalization = "log2",
                          .show_points = TRUE,
                          .show_confidence_intervals = TRUE)
dev.off()
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTrans/h_hybrid.pdf",
    width = 12, height = 3)
plotExpressionProfilePair(counts_all2_allele$cer[,gene_idxs],
                          counts_all2_allele$par[,gene_idxs],
                          info,
                          info,
                          .method = "line",
                          .normalization = "log2",
                          .show_points = TRUE,
                          .show_confidence_intervals = TRUE)
dev.off()
# module i, trans regulation
gene_idxs <- module_genedf10 |> filter(module_name == "i") |> select(gene_name) |> pull()
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTrans/i_parents.pdf",
    width = 12, height = 3)
plotExpressionProfilePair(counts_all2$cer[,gene_idxs],
                          counts_all2$par[,gene_idxs],
                          info,
                          info,
                          .method = "line",
                          .normalization = "log2",
                          .show_points = TRUE,
                          .show_confidence_intervals = TRUE)
dev.off()
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTrans/i_hybrid.pdf",
    width = 12, height = 3)
plotExpressionProfilePair(counts_all2_allele$cer[,gene_idxs],
                          counts_all2_allele$par[,gene_idxs],
                          info,
                          info,
                          .method = "line",
                          .normalization = "log2",
                          .show_points = TRUE,
                          .show_confidence_intervals = TRUE)
dev.off()

# random module 73
gene_idxs <- random_module_genedf10 |> filter(module_name == "ran73") |> select(gene_name) |> pull()
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTrans/ran73_parents.pdf",
    width = 12, height = 3)
plotExpressionProfilePair(counts_all2$cer[,gene_idxs],
                          counts_all2$par[,gene_idxs],
                          info,
                          info,
                          .method = "line",
                          .normalization = "log2",
                          .show_points = TRUE,
                          .show_confidence_intervals = TRUE)
dev.off()
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTrans/ran73_hybrid.pdf",
    width = 12, height = 3)
plotExpressionProfilePair(counts_all2_allele$cer[,gene_idxs],
                          counts_all2_allele$par[,gene_idxs],
                          info,
                          info,
                          .method = "line",
                          .normalization = "log2",
                          .show_points = TRUE,
                          .show_confidence_intervals = TRUE)
dev.off()



#### Which parent is dominant? ####
# Hybrid alleles of modules tend to be more correlated to each other
# But which parent do they resemble more?
# Module h shows cerevisiae dominance, is this a common pattern?

### Example of how Area Between Curves is calculated for Module J (ccm45)
gene_idxs_j <- module_genedf10 %>% filter(module_name == "ccm45") %>% 
  select(gene_name) %>% pull()
expr_j_cer <- map_dfc(gene_idxs_j, \(i) counts_all2$cer[,i]) %>% 
  scale() %>% as_tibble()
colnames(expr_j_cer) <- gene_idxs_j
gdf_j_cer <- expr_j_cer %>% 
  bind_cols(info) %>% 
  pivot_longer(cols = gene_idxs_j, names_to = "gene_name", values_to = "expr")
expr_j_par <- map_dfc(gene_idxs_j, \(i) counts_all2$par[,i]) %>% 
  scale()
colnames(expr_j_par) <- gene_idxs_j
gdf_j_par <- expr_j_par %>% 
  bind_cols(info) %>% 
  pivot_longer(cols = gene_idxs_j, names_to = "gene_name", values_to = "expr")
gdf_j <- inner_join(gdf_j_cer, gdf_j_par, by = setdiff(names(gdf_j_cer), "expr"),
                    suffix = c("_cer", "_par"))
# hybrid module j
expr_j_hyc <- map_dfc(gene_idxs_j, \(i) counts_all2_allele$cer[,i]) %>% 
  scale()
colnames(expr_j_hyc) <- gene_idxs_j
gdf_j_hyc <- expr_j_hyc %>% 
  bind_cols(info) %>% 
  pivot_longer(cols = gene_idxs_j, names_to = "gene_name", values_to = "expr")
expr_j_hyp <- map_dfc(gene_idxs_j, \(i) counts_all2_allele$par[,i]) %>% 
  scale()
colnames(expr_j_hyp) <- gene_idxs_j
gdf_j_hyp <- expr_j_hyp %>% 
  bind_cols(info) %>% 
  pivot_longer(cols = gene_idxs_j, names_to = "gene_name", values_to = "expr")
gdf_j_hyb <- inner_join(gdf_j_hyc, gdf_j_hyp, by = setdiff(names(gdf_j_hyc), "expr"),
                        suffix = c("_cer", "_par"))
# averages
avg_gdf_j <- filter(gdf_j, experiment == "LowPi") %>% 
  group_by(time_point_num) %>% 
  summarise(avgExpr_cer = mean(expr_cer),
            avgExpr_par = mean(expr_par))
# now hybrid average
avg_gdf_j_hyb <- filter(gdf_j_hyb, experiment == "LowPi") %>% 
  group_by(time_point_num) %>% 
  summarise(avgExpr_cer = mean(expr_cer),
            avgExpr_par = mean(expr_par))
# Making sure samples are not out of order
avg_gdf_j_sampleMatched <- tibble(cer = avg_gdf_j$avgExpr_cer,
                                  par = avg_gdf_j$avgExpr_par,
                                  time_point_num = avg_gdf_j$time_point_num) %>% 
  inner_join(avg_gdf_j_hyb, by = "time_point_num")
avg_gdf_j_sampleMatched$hyb <- avg_gdf_j_sampleMatched$avgExpr_cer + avg_gdf_j_sampleMatched$avgExpr_par # don't need to re-scale, I checked
plotdf <- avg_gdf_j_sampleMatched %>% select(time_point_num, cer, par, hyb) %>% 
  pivot_longer(cols = c("cer", "par", "hyb"), names_to = "type", values_to = "expr")
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/moduleJ_AreaBetweenCurves.pdf",
    width = 4, height = 3)
ggplot() + 
  geom_line(data = plotdf, aes(x = time_point_num, y = expr, color = type)) + theme_classic() +
  geom_ribbon(data = avg_gdf_j_sampleMatched, 
              aes(ymax = cer,
                  ymin = hyb,
                  x = time_point_num), alpha = 0.5) +
  xlab("Timepoint (min)") +
  ylab("Average expression") +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = c("S. cerevisiae", "F1 hybrid", "S. paradoxus")) +
  ggtitle("Area Between Curves\n for Module 45")
dev.off()

### Now calculating hybrid-parent similarity for all modules
# helper for below
# long story short, this is a more reliable way to determine which parental average 
# expression (cer or par) the hybrid is more similar to than correlation
# because correlation of two flat-ish lines is low but if hyb and one parent are both flat lines, we'd like to call them the most similar)
calculateAreaBetweenCurves <- function(.vec1, .vec2) {
  stopifnot(length(.vec1) == length(.vec2))
  output <- abs(.vec1 - .vec2) %>% sum()
  return(output)
}

# average cor values for the hybrid (hyc + hyp) versus each parent
modnames <- moduledf10 |> 
  filter(block_size >= 10) |> 
  select(module_name) |> 
  pull()
module_dominance_results_cer <- map(modnames, \(m) {
  cat("working on module", m, "\n")
  gene_idxs <- module_genedf10[module_genedf10$module_name == m,] %>% select(gene_name) %>% pull()
  counts_cer <- counts_all2$cer[, gene_idxs, drop = FALSE]
  counts_hyc <- counts_all2_allele$cer[, gene_idxs, drop = FALSE]
  counts_hyp <- counts_all2_allele$par[, gene_idxs, drop = FALSE]
  return(getExpressionSimilarity(counts_cer, counts_hyc + counts_hyp))
})
names(module_dominance_results_cer) <- modnames
module_dominance_results_par <- map(modnames, \(m) {
  cat("working on module", m, "\n")
  gene_idxs <- module_genedf10[module_genedf10$module_name == m,] %>% select(gene_name) %>% pull()
  counts_par <- counts_all2$par[, gene_idxs, drop = FALSE]
  counts_hyc <- counts_all2_allele$cer[, gene_idxs, drop = FALSE]
  counts_hyp <- counts_all2_allele$par[, gene_idxs, drop = FALSE]
  return(getExpressionSimilarity(counts_par, counts_hyc + counts_hyp))
})
names(module_dominance_results_par) <- modnames

module_dominance_by_experiment_cer <- map_dfr(module_dominance_results_cer, \(x) {
  ABC_CC <- calculateAreaBetweenCurves(x$avg1[info$experiment == "CC"], 
                                       x$avg2[info$experiment == "CC"]) # Area Between Curves, get it? Also I'm writing this at Arbor Brewing Company right now so it just fits
  ABC_LowN <- calculateAreaBetweenCurves(x$avg1[info$experiment == "LowN"], 
                                         x$avg2[info$experiment == "LowN"])
  ABC_LowPi <- calculateAreaBetweenCurves(x$avg1[info$experiment == "LowPi"], 
                                          x$avg2[info$experiment == "LowPi"])
  ABC_HAP4 <- calculateAreaBetweenCurves(x$avg1[info$experiment == "HAP4"], 
                                         x$avg2[info$experiment == "HAP4"])
  return(list(ABC_CC = ABC_CC,
              ABC_LowN = ABC_LowN,
              ABC_LowPi = ABC_LowPi,
              ABC_HAP4 = ABC_HAP4))
})
module_dominance_by_experiment_par <- map_dfr(module_dominance_results_par, \(x) {
  ABC_CC <- calculateAreaBetweenCurves(x$avg1[info$experiment == "CC"], 
                                       x$avg2[info$experiment == "CC"]) 
  ABC_LowN <- calculateAreaBetweenCurves(x$avg1[info$experiment == "LowN"], 
                                         x$avg2[info$experiment == "LowN"])
  ABC_LowPi <- calculateAreaBetweenCurves(x$avg1[info$experiment == "LowPi"], 
                                          x$avg2[info$experiment == "LowPi"])
  ABC_HAP4 <- calculateAreaBetweenCurves(x$avg1[info$experiment == "HAP4"], 
                                         x$avg2[info$experiment == "HAP4"])
  return(list(ABC_CC = ABC_CC,
              ABC_LowN = ABC_LowN,
              ABC_LowPi = ABC_LowPi,
              ABC_HAP4 = ABC_HAP4))
})

# It is more intuitive to think about how a module has higher similarity with one parent than the other,
# than that a module has higher difference with one parent than the other.
# Area Between Curves measures the difference between each parent and the hybrid (the higher the ABC, the higher the difference)
# that the smaller the ABC, the more similar the expression
module_dominance_by_experiment_cer$species <- "cer"
module_dominance_by_experiment_cer$module_name <- modnames
module_dominance_by_experiment_par$species <- "par"
module_dominance_by_experiment_par$module_name <- modnames
# Area Between Curves
plotdf <- bind_rows(module_dominance_by_experiment_cer, 
                    module_dominance_by_experiment_par) %>% 
  pivot_longer(cols = c("ABC_CC",
                        "ABC_LowN",
                        "ABC_LowPi",
                        "ABC_HAP4"),
               names_to = "experiment", values_to = "ABC") %>% 
  pivot_wider(id_cols = c("module_name", "experiment"), 
              names_from = "species", 
              values_from = "ABC")
plotdf <- left_join(plotdf, select(moduledf10, c("module_name", "coexpressed", "cer_color", "par_color")),
                    by = "module_name")

table(plotdf$module_name) %>% table() # should all have 4 entries, 1 for each experiment
plotdf$experiment <- gsub("ABC_", "", plotdf$experiment)

# adding random modules
random_modnames <- random_moduledf10 |> 
  filter(block_size >= 10) |> 
  select(module_name) |> 
  pull()
random_module_dominance_results_cer <- map(random_modnames, \(m) {
  cat("working on module", m, "\n")
  gene_idxs <- random_module_genedf10[random_module_genedf10$module_name == m,] %>% 
    select(gene_name) %>% pull()
  counts_cer <- counts_all2$cer[, gene_idxs, drop = FALSE]
  counts_hyc <- counts_all2_allele$cer[, gene_idxs, drop = FALSE]
  counts_hyp <- counts_all2_allele$par[, gene_idxs, drop = FALSE]
  return(getExpressionSimilarity(counts_cer, counts_hyc + counts_hyp))
})
names(random_module_dominance_results_cer) <- random_modnames
random_module_dominance_results_par <- map(random_modnames, \(m) {
  cat("working on module", m, "\n")
  gene_idxs <- random_module_genedf10[random_module_genedf10$module_name == m,] %>% 
    select(gene_name) %>% pull()
  counts_par <- counts_all2$par[, gene_idxs, drop = FALSE]
  counts_hyc <- counts_all2_allele$cer[, gene_idxs, drop = FALSE]
  counts_hyp <- counts_all2_allele$par[, gene_idxs, drop = FALSE]
  return(getExpressionSimilarity(counts_par, counts_hyc + counts_hyp))
})
names(random_module_dominance_results_par) <- random_modnames

random_module_dominance_by_experiment_cer <- map_dfr(random_module_dominance_results_cer, \(x) {
  ABC_CC <- calculateAreaBetweenCurves(x$avg1[info$experiment == "CC"], 
                                       x$avg2[info$experiment == "CC"]) # Area Between Curves, get it? Also I'm writing this at Arbor Brewing Company right now so it just fits
  ABC_LowN <- calculateAreaBetweenCurves(x$avg1[info$experiment == "LowN"], 
                                         x$avg2[info$experiment == "LowN"])
  ABC_LowPi <- calculateAreaBetweenCurves(x$avg1[info$experiment == "LowPi"], 
                                          x$avg2[info$experiment == "LowPi"])
  ABC_HAP4 <- calculateAreaBetweenCurves(x$avg1[info$experiment == "HAP4"], 
                                         x$avg2[info$experiment == "HAP4"])
  return(list(ABC_CC = ABC_CC,
              ABC_LowN = ABC_LowN,
              ABC_LowPi = ABC_LowPi,
              ABC_HAP4 = ABC_HAP4))
})
random_module_dominance_by_experiment_par <- map_dfr(random_module_dominance_results_par, \(x) {
  ABC_CC <- calculateAreaBetweenCurves(x$avg1[info$experiment == "CC"], 
                                       x$avg2[info$experiment == "CC"]) 
  ABC_LowN <- calculateAreaBetweenCurves(x$avg1[info$experiment == "LowN"], 
                                         x$avg2[info$experiment == "LowN"])
  ABC_LowPi <- calculateAreaBetweenCurves(x$avg1[info$experiment == "LowPi"], 
                                          x$avg2[info$experiment == "LowPi"])
  ABC_HAP4 <- calculateAreaBetweenCurves(x$avg1[info$experiment == "HAP4"], 
                                         x$avg2[info$experiment == "HAP4"])
  return(list(ABC_CC = ABC_CC,
              ABC_LowN = ABC_LowN,
              ABC_LowPi = ABC_LowPi,
              ABC_HAP4 = ABC_HAP4))
})
random_module_dominance_by_experiment_cer$species <- "cer"
random_module_dominance_by_experiment_cer$module_name <- random_modnames
random_module_dominance_by_experiment_par$species <- "par"
random_module_dominance_by_experiment_par$module_name <- random_modnames
# Area Between Curves
rplotdf <- bind_rows(random_module_dominance_by_experiment_cer, 
                     random_module_dominance_by_experiment_par) %>% 
  pivot_longer(cols = c("ABC_CC",
                        "ABC_LowN",
                        "ABC_LowPi",
                        "ABC_HAP4"),
               names_to = "experiment", values_to = "ABC") %>% 
  pivot_wider(id_cols = c("module_name", "experiment"), 
              names_from = "species", 
              values_from = "ABC")
rplotdf <- left_join(rplotdf, select(moduledf10, c("module_name", "coexpressed", "cer_color", "par_color")),
                     by = "module_name")

table(rplotdf$module_name) %>% table() # should all have 4 entries, 1 for each experiment
rplotdf$experiment <- gsub("ABC_", "", rplotdf$experiment)
rplotdf$coexpressed <- "random"
plotdf <- bind_rows(plotdf, rplotdf)

# Which Random is the farthest from y=x?
random_bounds <- plotdf %>% filter(coexpressed == "random") %>% select(cer, par) |> 
  mutate(distance_up_from_yx = abs(cer - par)) %>% select(distance_up_from_yx) %>% pull() %>% max(na.rm = TRUE) # farthest distance out from y=x, either along rise or run

### Hybrid parent dominance for all modules
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/hybrid_parent_dominance.pdf",
    width = 5, height = 4)
plotdf %>% 
  ggplot(aes(x = cer, y = par)) + 
  geom_point(aes(color = coexpressed, shape = experiment)) +
  geom_point(data = filter(plotdf, experiment == "LowPi" & module_name == "j"),
             aes(x = cer, y = par), color = "red", shape = "+", size = 4) +
  # purple dotted lines bounding random modules
  geom_abline(slope = 1, intercept = random_bounds, color = "purple1", linetype="dashed") +
  geom_abline(slope = 1, intercept = -random_bounds, color = "purple1", linetype="dashed") +
  #geom_text(data = plotdf, aes(x = cer, y = par, label = module_name), nudge_x = 1) +
  scale_color_discrete(type = c("gold",
                                "mediumseagreen",
                                "orange1",
                                "blue2",
                                "grey",
                                "purple1"),
                       limits = c("conserved co-expressed",
                                  "diverged co-expressed",
                                  "S. cerevisiae co-expressed",
                                  "S. paradoxus co-expressed",
                                  "never co-expressed",
                                  "random")) +
  scale_shape_discrete(limits = c("CC", "HAP4", "LowN", "LowPi"),
                       labels = c("Urea Shock", "Saturated Growth", "Low Nitrogen", "Low Phosphorus")) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  geom_abline(slope = 1, intercept = 0, color = "midnightblue") +
  ggtitle("Expression similarity of \nhybrid with each parent") +
  xlab("Area Between Curves\nhybrid - S. cerevisiae") +
  ylab("Area Between Curves\nhybrid - S. paradoxus")
dev.off()
# neither appear to be dominant very frequently

#### Extreme subsets in hybrid ####

# We've seen that the extreme subsets of parental ccms are still strongly correlated with the rest of the parental module
# And we've seen that parental ccms become much more correlated with each other in the hybrid
# so... do the extreme subsets also become more correlated? Do they remain extreme?
# I'm expecting they will become more correlated but remain extreme
# indicating that this cis regulatory tuning is affecting level while trans affects shape

coef_thresh <- 0.25
p_thresh <- 1e-5

# same function in environmental_patterns.R
getDivergentAndConservedGeneIdxs <- function(.ccm_name, .module_gdf = module_genedf10, .allele_or_species = "species") {
  gene_idxs <- .module_gdf |> 
    filter(module_name == .ccm_name) |> 
    select(gene_name) |> pull()
  moddf <- spaldf |> filter(gene_name %in% gene_idxs & 
                              coefficient == .allele_or_species &
                              experiment != "all") |> 
    mutate(sig = abs(effect_size) > coef_thresh & pvalue < p_thresh) |> 
    mutate(adj_effect_size = if_else(sig, true = effect_size, false = 0)) |> 
    group_by(gene_name) |> summarise(mean_effect = mean(adj_effect_size)) |> 
    select(mean_effect, gene_name) |> 
    mutate(direction = if_else(mean_effect < -coef_thresh,
                               true = "up_par", 
                               false = if_else(mean_effect > coef_thresh,
                                               true = "up_cer", false = "conserved")))
  table(moddf$direction)
  up_par_idxs <- moddf |> filter(direction == "up_par") |> 
    select(gene_name) |> pull()
  up_cer_idxs <- moddf |> filter(direction == "up_cer") |> 
    select(gene_name) |> pull()
  conserved_idxs <- moddf |> filter(direction == "conserved") |> 
    select(gene_name) |> pull()
  return(list(conserved = conserved_idxs,
              up_cer = up_cer_idxs,
              up_par = up_par_idxs))
}

# divergent CCM
divergent_CCM_name <- "ccm43"
idx_list <- getDivergentAndConservedGeneIdxs(divergent_CCM_name, 
                                             .module_gdf = module_genedf10, 
                                             .allele_or_species = "species")

# up cer in parents
plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,idx_list$conserved],
                             .cts2 = counts_all2$par[,idx_list$conserved],
                             .cts3 = counts_all2$cer[,idx_list$up_cer],
                             .cts4 = counts_all2$par[,idx_list$up_cer],
                             .info1 = info,
                             .info2 = info,
                             .info3 = info,
                             .info4 = info,
                             .color1 = "orange1",
                             .color2 = "blue2",
                             .color3 = "orange4",
                             .color4 = "blue4",
                             .name1 = "conserved genes, S. cerevisiae",
                             .name2 = "conserved genes, S. paradoxus",
                             .name3 = "up cer genes, S. cerevisiae",
                             .name4 = "up cer genes, S. paradoxus",
                             .normalization = "log2",
                             .show_points = FALSE,
                             .show_confidence_intervals = FALSE,
                             .method = "line") |> 
  annotate_figure(top = paste(divergent_CCM_name, "up cer parental expression"))

# up cer in hybrid
plotExpressionProfileQuartet(.cts1 = counts_all2_allele$cer[,idx_list$conserved],
                             .cts2 = counts_all2_allele$par[,idx_list$conserved],
                             .cts3 = counts_all2_allele$cer[,idx_list$up_cer],
                             .cts4 = counts_all2_allele$par[,idx_list$up_cer],
                             .info1 = info,
                             .info2 = info,
                             .info3 = info,
                             .info4 = info,
                             .color1 = "orange1",
                             .color2 = "blue2",
                             .color3 = "orange4",
                             .color4 = "blue4",
                             .name1 = "conserved genes, cer allele",
                             .name2 = "conserved genes, par allele",
                             .name3 = "up cer genes, cer allele",
                             .name4 = "up cer genes, par allele",
                             .normalization = "log2",
                             .show_points = FALSE,
                             .show_confidence_intervals = FALSE,
                             .method = "line") |> 
  annotate_figure(top = paste(divergent_CCM_name, "up cer hybrid expression"))

# up par in parents
plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,idx_list$conserved],
                             .cts2 = counts_all2$par[,idx_list$conserved],
                             .cts3 = counts_all2$cer[,idx_list$up_par],
                             .cts4 = counts_all2$par[,idx_list$up_par],
                             .info1 = info,
                             .info2 = info,
                             .info3 = info,
                             .info4 = info,
                             .color1 = "orange1",
                             .color2 = "blue2",
                             .color3 = "orange4",
                             .color4 = "blue4",
                             .name1 = "conserved genes, S. cerevisiae",
                             .name2 = "conserved genes, S. paradoxus",
                             .name3 = "up par genes, S. cerevisiae",
                             .name4 = "up par genes, S. paradoxus",
                             .normalization = "log2",
                             .show_points = FALSE,
                             .show_confidence_intervals = FALSE,
                             .method = "line") |> 
  annotate_figure(top = paste(divergent_CCM_name, "up par parental expression"))

# up par in hybrid
plotExpressionProfileQuartet(.cts1 = counts_all2_allele$cer[,idx_list$conserved],
                             .cts2 = counts_all2_allele$par[,idx_list$conserved],
                             .cts3 = counts_all2_allele$cer[,idx_list$up_par],
                             .cts4 = counts_all2_allele$par[,idx_list$up_par],
                             .info1 = info,
                             .info2 = info,
                             .info3 = info,
                             .info4 = info,
                             .color1 = "orange1",
                             .color2 = "blue2",
                             .color3 = "orange4",
                             .color4 = "blue4",
                             .name1 = "conserved genes, cer allele",
                             .name2 = "conserved genes, par allele",
                             .name3 = "up par genes, cer allele",
                             .name4 = "up par genes, par allele",
                             .normalization = "log2",
                             .show_points = FALSE,
                             .show_confidence_intervals = FALSE,
                             .method = "line") |> 
  annotate_figure(top = paste(divergent_CCM_name, "up par hybrid expression"))

# Are the same genes up cer and up par in the hybrid?
for (ccm in moduledf10$module_name[moduledf10$is_CCM]) {
  idx_list_parents <- getDivergentAndConservedGeneIdxs(ccm, 
                                                       .module_gdf = module_genedf10, 
                                                       .allele_or_species = "species")
  idx_list_hybrid <- getDivergentAndConservedGeneIdxs(ccm, 
                                                      .module_gdf = module_genedf10, 
                                                      .allele_or_species = "allele")
  ccm_idxs <- module_genedf10$gene_name[module_genedf10$module_name == ccm]
  # guess what's an excellent way to visualize this conservation?
  extremeSubsetMatrix <- makeBipartiteMatrix(.top_colors = if_else(ccm_idxs %in% idx_list_parents$conserved,
                                                                   true = "conserved", 
                                                                   false = if_else(ccm_idxs %in% idx_list_parents$up_cer,
                                                                                   true = "upregulated cerevisiae", 
                                                                                   false = "upregulated paradoxus")),
                                             .bottom_colors = if_else(ccm_idxs %in% idx_list_hybrid$conserved,
                                                                      true = "conserved", 
                                                                      false = if_else(ccm_idxs %in% idx_list_hybrid$up_cer,
                                                                                      true = "upregulated cerevisiae", 
                                                                                      false = "upregulated paradoxus")))
  color_table <- data.frame(colors =  c("gold", "orange1", "blue2"))
  rownames(color_table) <- c("conserved", "upregulated cerevisiae", "upregulated paradoxus")
  plotweb(extremeSubsetMatrix$matrix,
          method = "normal",
          col.high = color_table[colnames(extremeSubsetMatrix$matrix), "colors"],
          col.low = color_table[rownames(extremeSubsetMatrix$matrix), "colors"],
          sequence = list(seq.high = c("conserved", setdiff(colnames(extremeSubsetMatrix$matrix), "conserved")),
                          seq.low = c("conserved", setdiff(rownames(extremeSubsetMatrix$matrix), "conserved"))))
  title(main = ccm)
}



#### Identifying hybrid CCMs ####

# TODO: we do see that several modules get globbed together
# in the hybrid networks, so it might be useful to at least
# know what they are and see if there's a reason (i.e. are they
# the trans-varying modules whose correlation improves the best
# in the hybrid?)

### cer-hyc
cer_colors <- rep(unique(colors10$cer), length(unique(colors10$hyc))) %>% as.character()
hyc_colors <- sapply(unique(colors10$hyc), rep, times = length(unique(colors10$cer))) %>% as.character()
CCMdf_cerhyc <- tibble(cer_color = cer_colors,
                       hyc_color = hyc_colors)
CCMdf_cerhyc$p_value <- map2(CCMdf_cerhyc$cer_color, CCMdf_cerhyc$hyc_color, \(x, y) {
  cat("starting on", x, "and", y, "\n")
  return(callConservedConnection(.color_name1 = x,
                                 .color_name2 = y,
                                 .colors1 = colors10$cer,
                                 .colors2 = colors10$hyc))
}) %>% unlist()
# CCMdf_cerhyc$block_size <- map2(CCMdf_cerhyc$cer_color, CCMdf_cerhyc$hyc_color, \(x, y) {
#   return(sum(colors10$cer == x & colors10$hyc == y))
# }) %>% unlist()

# par-hyp, identify cis-varying modules
par_colors <- rep(unique(colors10$par), length(unique(colors10$hyp))) %>% as.character()
hyp_colors <- sapply(unique(colors10$hyp), rep, times = length(unique(colors10$par))) %>% as.character()
CCMdf_parhyp <- tibble(par_color = par_colors,
                       hyp_color = hyp_colors)
CCMdf_parhyp$p_value <- map2(CCMdf_parhyp$par_color, CCMdf_parhyp$hyp_color, \(x, y) {
  cat("starting on", x, "and", y, "\n")
  return(callConservedConnection(.color_name1 = x,
                                 .color_name2 = y,
                                 .colors1 = colors10$par,
                                 .colors2 = colors10$hyp))
}) %>% unlist()
# CCMdf_parhyp$block_size <- map2(CCMdf_parhyp$par_color, CCMdf_parhyp$hyp_color, \(x, y) {
#   return(sum(colors10$par == x & colors10$hyp == y))
# }) %>% unlist()

# hyc-hyp, , identify trans-varying modules
hyc_colors <- rep(unique(colors10$hyc), length(unique(colors10$hyp))) %>% as.character()
hyp_colors <- sapply(unique(colors10$hyp), rep, times = length(unique(colors10$hyc))) %>% as.character()
CCMdf_hychyp <- tibble(hyc_color = hyc_colors,
                       hyp_color = hyp_colors)
CCMdf_hychyp$p_value <- map2(CCMdf_hychyp$hyc_color, CCMdf_hychyp$hyp_color, \(x, y) {
  cat("starting on", x, "and", y, "\n")
  return(callConservedConnection(.color_name1 = x,
                                 .color_name2 = y,
                                 .colors1 = colors10$hyc,
                                 .colors2 = colors10$hyp))
}) %>% unlist()
# CCMdf_hychyp$block_size <- map2(CCMdf_hychyp$hyc_color, CCMdf_hychyp$hyp_color, \(x, y) {
#   return(sum(colors10$hyc == x & colors10$hyp == y))
# }) %>% unlist()

save(CCMdf_cerhyc, CCMdf_parhyp, CCMdf_hychyp, file = "data_files/CCMs_hybrid.RData")
load(file = "data_files/CCMs_hybrid.RData")
load(file = "data_files/CCMs.RData")

# which groups are more connected than expected by chance?
CCMdf_cerhyc %>% filter(p_value < 1/10000)
CCMdf_parhyp %>% filter(p_value < 1/10000)
CCMdf_hychyp %>% filter(p_value < 1/10000) # looks like a lot

# are CCMs in one category more or less likely to be a CCM in another category?
# for this we will construct an upset plot (Venn diagram with 4+ categories)
CCM_genedf <- rename(CCM_genedf, "is_CCM_cerpar"="is_conserved_block", "p_value_cerpar"="p_value")
CCMdf_cerhyc <- rename(CCMdf_cerhyc, "p_value_cerhyc"="p_value")
CCMdf_parhyp <- rename(CCMdf_parhyp, "p_value_parhyp"="p_value")
CCMdf_hychyp <- rename(CCMdf_hychyp, "p_value_hychyp"="p_value")

# Removing grey CCMs
# setting grey to a pvalue of 1 to effectively remove them (without throwing a dplyr error for not having a category of gene color present)
CCMdf[CCMdf$cer_color == "grey" | CCMdf$par_color == "grey",]$p_value <- 1
CCMdf_cerhyc[CCMdf_cerhyc$hyc_color == "grey" | CCMdf_cerhyc$cer_color == "grey",]$p_value_cerhyc <- 1
CCMdf_parhyp[CCMdf_parhyp$par_color == "grey" | CCMdf_parhyp$hyp_color == "grey",]$p_value_parhyp <- 1
CCMdf_hychyp[CCMdf_hychyp$hyc_color == "grey" | CCMdf_hychyp$hyp_color == "grey",]$p_value_hychyp <- 1

CCM_genedf$hyc_color <- colors10$hyc
CCM_genedf$hyp_color <- colors10$hyp

CCM_genedf <- left_join(CCM_genedf, CCMdf_cerhyc, by = c("cer_color", "hyc_color")) %>% 
  left_join(CCMdf_parhyp, by = c("par_color", "hyp_color")) %>% 
  left_join(CCMdf_hychyp, by = c("hyc_color", "hyp_color"))

CCM_genedf$is_CCM_cerhyc <- CCM_genedf$p_value_cerhyc < 0.05/(length(unique(colors10$cer))*length(unique(colors10$hyc)))
CCM_genedf$is_CCM_parhyp <- CCM_genedf$p_value_parhyp < 0.05/(length(unique(colors10$par))*length(unique(colors10$hyp)))
CCM_genedf$is_CCM_hychyp <- CCM_genedf$p_value_hychyp < 0.05/(length(unique(colors10$hyc))*length(unique(colors10$hyp)))
  
### Supplemental Figure: Upset plot of if any of those CCM groups co-occur more than others
# Adapted from FigNo_singleGene...etc...
library(ComplexHeatmap)
# species
plotdf <- CCM_genedf[, c("is_CCM_cerpar", "is_CCM_cerhyc",
                     "is_CCM_parhyp", "is_CCM_hychyp")] %>% 
  apply(MARGIN = 2, FUN = as.integer) %>% as.matrix()
colnames(plotdf) <- c( "cerevisiae - paradoxus", "cerevisiae - hybrid cer allele",
                      "paradoxus - hybrid par allele", "hybrid alleles")
# replacing NA with 0 because make_comb_mat cannot handle NA for some reason
plotdf[is.na(plotdf)] <- 0
plotdf <- make_comb_mat(plotdf)
p <- UpSet(plotdf, set_order = c("cerevisiae - paradoxus", "cerevisiae - hybrid cer allele",
                                 "paradoxus - hybrid par allele", "hybrid alleles"), 
           comb_order = order(comb_size(plotdf)), # this, I will say, makes ggplot look like a freakin' dream. I literally just wanted to rename the histogram y axis and add counts to the top of the bars
           top_annotation = HeatmapAnnotation( 
             "number of genes" = anno_barplot(comb_size(plotdf), 
                                              ylim = c(0, max(comb_size(plotdf))*1.1),
                                              border = FALSE, 
                                              gp = gpar(fill = "black"), 
                                              height = unit(4, "cm")), 
             annotation_name_side = "left", 
             annotation_name_rot = 90))

pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/supp4_upsetOfHybridCCMs.pdf",
    width = 5, height = 3)
draw(p)
decorate_annotation("number of genes", {
  grid.text(comb_size(plotdf)[column_order(p)], x = seq_along(comb_size(plotdf)), y = unit(comb_size(plotdf)[column_order(p)], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})
dev.off()

#### CCM case example ####
# Are expression-divergent CCMs more likely to have trans-reg divergence than expression-conserved CCMs?
# to address this, I'll start out with investigating a case example
random_ccm_name <- moduledf10 |> filter(is_CCM) |> select(module_name) |> pull() |> sample(1)

# TODO: adapt this for a random/generic ccm, if it seems useful to
# have a case example, otherwise archive

# first, are j genes really diverging in expression between cer and par (parental species)? We'll compare to h, which seem very conserved

# comparing avg expression in cer/par in j vs h
gene_idxs_j <- module_genedf %>% filter(module_name_cer == "j") %>% 
  select(gene_name) %>% pull()
gene_idxs_h <- module_genedf %>% filter(module_name_cer == "h") %>% 
  select(gene_name) %>% pull()
# centered and scaled average expression of genes in each module
expr_j_cer <- counts_all2$cer[,gene_idxs_j] %>% 
  scale() %>% 
  rowMeans(na.rm = TRUE) # each observation is a condition
expr_j_par <- counts_all2$par[,gene_idxs_j] %>% 
  scale() %>% 
  rowMeans(na.rm = TRUE)
gdf_j <- bind_cols(tibble(cer = expr_j_cer,
                          par = expr_j_par),
                   info) 
### 3D: Module J expression in cer vs par
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTransBasis/moduleJ/3A.pdf",
    width = 4, height = 3)
ggplot(gdf_j, aes(x = cer, y = par)) + 
  geom_point(aes(color = experiment)) + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_classic() +
  xlab("S. cerevisiae") +
  ylab("S. paradoxus") +
  ggtitle("Expression of Module J \n(condition-matched and\n averaged across all genes)") # definitely something going on in LowPi
dev.off()
# repeat for h
expr_h_cer <- counts_all2$cer[,gene_idxs_h] %>% 
  scale() %>% 
  rowMeans(na.rm = TRUE) # each observation is a condition
expr_h_par <- counts_all2$par[,gene_idxs_h] %>% 
  scale() %>% 
  rowMeans(na.rm = TRUE)
gdf_h <- bind_cols(tibble(cer = expr_h_cer,
                          par = expr_h_par),
                   info) 
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/supp8_moduleJ/module_H_average.pdf",
    width = 4, height = 3)
ggplot(gdf_h, aes(x = cer, y = par)) + 
  geom_point(aes(color = experiment)) + 
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() +
  xlab("S. cerevisiae") +
  ylab("S. paradoxus") +
  ggtitle("Expression of Module H \n(condition-matched and\n averaged across all genes)") # no divergence to observe
dev.off()
# centered and scaled BUT NOT average expression of genes in each module
expr_j_cer <- map_dfc(gene_idxs_j, \(i) counts_all2$cer[,i]) %>% 
  scale() %>% as_tibble()
colnames(expr_j_cer) <- gene_idxs_j
gdf_j_cer <- expr_j_cer %>% 
  bind_cols(info) %>% 
  pivot_longer(cols = c(gene_idxs_j, ), names_to = "gene_name", values_to = "expr")
expr_j_par <- map_dfc(gene_idxs_j, \(i) counts_all2$par[,i]) %>% 
  scale()
colnames(expr_j_par) <- gene_idxs_j
gdf_j_par <- expr_j_par %>% 
  bind_cols(info) %>% 
  pivot_longer(cols = gene_idxs_j, names_to = "gene_name", values_to = "expr")
gdf_j <- inner_join(gdf_j_cer, gdf_j_par, by = setdiff(names(gdf_j_cer), "expr"),
                    suffix = c("_cer", "_par"))
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/supp8_moduleJ/module_J_allGenes.pdf",
    width = 4, height = 3)
ggplot(filter(gdf_j, experiment == "LowPi"), aes(x = expr_cer, y = expr_par)) +
  geom_point(aes(color = time_point_num)) +
  scale_color_continuous(name = "Timepoint (min)") +
  xlab("S. cerevisiae") +
  ylab("S. paradoxus") +
  ggtitle("Expression of Module J \n(condition-matched, normalized \ncounts of all genes)") +
  theme_classic() # trend is still present. At later time points expr spikes in paradoxus
dev.off()
# (I also checked, there are 22 genes represented in the paradoxus peak---this isn't being driven by one gene)
high_J_genes <- gdf_j %>% filter(expr_par > 4 & experiment == "LowPi") %>% select(gene_name) %>% unique()
high_J_genes %>% print(n = 22)
write.table(high_J_genes, file = "gene_ontology/module_j_high_genes.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

# is it present in h?
expr_h_cer <- map_dfc(gene_idxs_h, \(i) counts_all2$cer[,i]) %>% 
  scale()
colnames(expr_h_cer) <- gene_idxs_h
gdf_h_cer <- expr_h_cer %>% 
  bind_cols(info) %>% 
  pivot_longer(cols = c(gene_idxs_h, ), names_to = "gene_name", values_to = "expr")
expr_h_par <- map_dfc(gene_idxs_h, \(i) counts_all2$par[,i]) %>% 
  scale()
colnames(expr_h_par) <- gene_idxs_h
gdf_h_par <- expr_h_par %>% 
  bind_cols(info) %>% 
  pivot_longer(cols = gene_idxs_h, names_to = "gene_name", values_to = "expr")
gdf_h <- inner_join(gdf_h_cer, gdf_h_par, by = setdiff(names(gdf_h_cer), "expr"),
                    suffix = c("_cer", "_par"))
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/supp8_moduleJ/module_H_allGenes.pdf",
    width = 4, height = 3)
ggplot(filter(gdf_h, experiment == "LowPi"), aes(x = expr_cer, y = expr_par)) +
  geom_point(aes(color = time_point_num)) +
  xlab("S. cerevisiae") +
  ylab("S. paradoxus") +
  ggtitle("Expression of Module H \n(condition-matched, normalized \ncounts of all genes)") +
theme_classic() # nope
dev.off()

# 3A-B
# plot mostly useful for comparing with hybrid expr below
# how expression of every gene in module J increases over time
# in par:
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/supp8_moduleJ/module_J_timepoint_par.pdf",
    width = 4, height = 3)
ggplot(filter(gdf_j, experiment == "LowPi"), aes(x = time_point_num, y = expr_par)) + 
  geom_point(aes(color = time_point_num)) + 
  geom_line(data = filter(gdf_j, experiment == "LowPi") %>% group_by(time_point_num) %>% summarise(avgExpr = mean(expr_par)),
            aes(x = time_point_num, y = avgExpr), color = "#619CFF") +
  xlab("Timepoint (min)") +
  ylab("Normalized expression level") +
  ggtitle("Module J in S. paradoxus") +
  ylim(c(-3, 10)) +
  theme_classic()
dev.off()
# not in cer:
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/supp8_moduleJ/module_J_timepoint_cer.pdf",
    width = 4, height = 3)
ggplot(filter(gdf_j, experiment == "LowPi"), aes(x = time_point_num, y = expr_cer)) + 
  geom_point(aes(color = time_point_num)) + 
  geom_line(data = filter(gdf_j, experiment == "LowPi") %>% group_by(time_point_num) %>% summarise(avgExpr = mean(expr_cer)),
            aes(x = time_point_num, y = avgExpr), color = "#F8766D") +
  xlab("Timepoint (min)") +
  ylab("Normalized expression level") +
  ggtitle("Module J in S. cerevisiae") +
  ylim(c(-3, 10)) +
  theme_classic()
dev.off()
# now that we've established the divergence is real and specific to j, is this divergence in cis or trans?

# does hybrid module j recapitulate parental finding?
expr_j_cer <- counts_all2_allele$cer[,gene_idxs_j] %>% 
  scale() %>% 
  rowMeans(na.rm = TRUE) # each observation is a condition
expr_j_par <- counts_all2_allele$par[,gene_idxs_j] %>% 
  scale() %>% 
  rowMeans(na.rm = TRUE)
gdf_j <- bind_cols(tibble(cer = expr_j_cer,
                          par = expr_j_par),
                   info_allele) 
# 3E
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/supp8_moduleJ/module_J_average_hybrid.pdf",
    width = 4, height = 3)
ggplot(gdf_j, aes(x = cer, y = par)) + 
  geom_point(aes(color = experiment)) + 
  geom_abline(slope = 1, intercept = 0) +
  xlab("F1 hybrid, cerevisiae allele") +
  ylab("F1 hybrid, paradoxus allele") +
  ggtitle("Expression of Module J\n(condition-matched and \n averaged across all genes)") +
  theme_classic() # seems pretty trans? There's like a slight remnant of allele specificity
dev.off()
# just LowPi
ggplot(filter(gdf_j, experiment == "LowPi"), aes(x = cer, y = par)) + 
  geom_point(aes(color = time_point_num)) + 
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() # you can see the remnant in the sense that the later timepoints tend to be on the par side of y=x line

# same thing but not average
expr_j_cer <- map_dfc(gene_idxs_j, \(i) counts_all2_allele$cer[,i]) %>% 
  scale()
colnames(expr_j_cer) <- gene_idxs_j
gdf_j_cer <- expr_j_cer %>% 
  bind_cols(info_allele) %>% 
  pivot_longer(cols = c(gene_idxs_j, ), names_to = "gene_name", values_to = "expr")
expr_j_par <- map_dfc(gene_idxs_j, \(i) counts_all2_allele$par[,i]) %>% 
  scale()
colnames(expr_j_par) <- gene_idxs_j
gdf_j_par <- expr_j_par %>% 
  bind_cols(info_allele) %>% 
  pivot_longer(cols = gene_idxs_j, names_to = "gene_name", values_to = "expr")
gdf_j <- inner_join(gdf_j_cer, gdf_j_par, by = setdiff(names(gdf_j_cer), "expr"),
                    suffix = c("_cer", "_par"))
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/supp8_moduleJ/module_J_allGenes_hybrid.pdf",
    width = 4, height = 3)
ggplot(filter(gdf_j, experiment == "LowPi"), aes(x = expr_cer, y = expr_par)) +
  geom_point(aes(color = time_point_num)) + 
  xlab("F1 hybrid, cerevisiae allele") +
  ylab("F1 hybrid, paradoxus allele") +
  ggtitle("Expression of Module J\n") +
  theme_classic() # yeah nothing
dev.off()
# 3C
gdf_j$expr_hyb <- gdf_j$expr_cer + gdf_j$expr_par
# is the module still upregulated in later timepoints? Doesn't really look like it
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/supp8_moduleJ/module_J_timepoint_hyb.pdf",
    width = 4, height = 3)
ggplot(filter(gdf_j, experiment == "LowPi"), aes(x = time_point_num, y = expr_hyb)) + 
  geom_point(aes(color = time_point_num)) +
  geom_line(data = filter(gdf_j, experiment == "LowPi") %>% group_by(time_point_num) %>% summarise(avgExpr = mean(expr_hyb)),
            aes(x = time_point_num, y = avgExpr), color = "#00BA38") +
  xlab("Timepoint (min)") +
  ylab("Normalized expression level") +
  ggtitle("Module J in F1 Hybrid") +
  ylim(c(-3, 10)) +
  theme_classic()
dev.off()


############################### Archive ######################################## 
#### line plots of each module type in parents and hybrid, to emphasize degree of improvement ####
# pdf("../../aligning_the_molecular_phenotype/paper_figures/CisTrans/allele_vs_species_lines.pdf",
#     width = 7, height = 4)
# plotdf |> 
#   pivot_longer(cols = c("allele", "species"), 
#                names_to = "allele_or_species",
#                values_to = "dcor") |> 
#   group_by(coexpressed, allele_or_species) |> 
#   ggplot() +
#   geom_line(aes(color = coexpressed, group = paste(module_name, experiment), 
#                 x = factor(paste(coexpressed, allele_or_species),
#                            levels = c("conserved co-expressed species", "conserved co-expressed allele",
#                                       "diverged co-expressed species", "diverged co-expressed allele",
#                                       "S. cerevisiae co-expressed species", "S. cerevisiae co-expressed allele",
#                                       "S. paradoxus co-expressed species", "S. paradoxus co-expressed allele",
#                                       "never co-expressed species", "never co-expressed allele",
#                                       "random species", "random allele")),
#                 y = dcor)) + 
#   xlab("") +
#   ylab("module correlation") +
#   scale_color_discrete(type = c("gold",
#                                 "purple1"),
#                        limits = c("conserved co-expressed",
#                                   "random"),
#                        labels = c("high confidence modules",
#                                   "random modules")) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())
# dev.off()
# Bipartite-style cis-trans
# archived b/c it's even more confusing than module-level cis-trans
# which in turn is even more confusing than single gene cis-trans
## Probably won't use: Hypothetical schematics of what full-trans and full-cis would mean for hybrid networks ####
# # 3A: Full-trans expectation: 
# # no difference between hyc and hyp
# # aka create toy gene colors and bipartite it to itself
# full_trans_colors_hyb <- c(rep("blue", 300), rep("red", 100), rep("green", 150), rep("yellow", 75), rep("purple", 43))
# full_trans_colors_cer <- c(rep("blue", 100), rep("red", 160), rep("greenyellow", 150), rep("purple", 53), rep("brown", 15)) %>% sample(size = length(full_trans_colors_hyb), replace = TRUE)
# full_trans_colors_par <- c(rep("blue", 60), rep("red", 200), rep("yellow", 15), rep("turquoise", 50), rep("magenta", 140), rep("greenyellow", 150), rep("purple", 53)) %>% sample(size = length(full_trans_colors_hyb), replace = TRUE)
# # hyb vs cer
# full_trans_cerhyb_color_counts <- makeBipartiteMatrix(.top_colors = full_trans_colors_cer, .bottom_colors = full_trans_colors_hyb)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTransBasis/hypotheticals/3A_full_trans_cerhyb.pdf",
#     width = 4, height = 3)
# plotweb(full_trans_cerhyb_color_counts$matrix, high.lablength = 0, low.lablength = 0,
#         col.high = full_trans_cerhyb_color_counts$colors_high,
#         col.low = full_trans_cerhyb_color_counts$colors_low)
# dev.off()
# 
# # hyb vs par
# full_trans_parhyb_color_counts <- makeBipartiteMatrix(.top_colors = full_trans_colors_par, .bottom_colors = full_trans_colors_hyb)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTransBasis/hypotheticals/3A_full_trans_parhyb.pdf",
#     width = 4, height = 3)
# plotweb(full_trans_parhyb_color_counts$matrix, high.lablength = 0, low.lablength = 0,
#         col.high = full_trans_parhyb_color_counts$colors_high,
#         col.low = full_trans_parhyb_color_counts$colors_low)
# dev.off()
# 
# # hyb vs hyb, no difference
# full_trans_hybhyb_color_counts <- makeBipartiteMatrix(.top_colors = full_trans_colors_hyb, .bottom_colors = full_trans_colors_hyb)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTransBasis/hypotheticals/3A_full_trans_hybhyb.pdf",
#     width = 4, height = 3)
# plotweb(full_trans_hybhyb_color_counts$matrix, high.lablength = 0, low.lablength = 0,
#         col.high = full_trans_hybhyb_color_counts$colors_high,
#         col.low = full_trans_hybhyb_color_counts$colors_low)
# dev.off()
# 
# # 3B: Full-cis expectation: no difference btwn cer and hyc, and par and hyp
# full_cis_colors_cer <- c(rep("blue", 30), rep("purple", 143), rep("red", 100), rep("green", 150), rep("yellow", 175))
# full_cis_colors_par <- c(rep("blue", 60), rep("red", 200), rep("yellow", 15), rep("turquoise", 50), rep("magenta", 140), rep("greenyellow", 150), rep("purple", 53)) %>% sample(size = length(full_cis_colors_cer), replace = TRUE)
# 
# # cer vs hyc
# full_cis_cerhyc_color_counts <- makeBipartiteMatrix(.top_colors = full_cis_colors_cer, .bottom_colors = full_cis_colors_cer)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTransBasis/hypotheticals/3B_full_cis_cerhyc.pdf",
#     width = 4, height = 3)
# plotweb(full_cis_cerhyc_color_counts$matrix, high.lablength = 0, low.lablength = 0,
#         col.high = full_cis_cerhyc_color_counts$colors_high,
#         col.low = full_cis_cerhyc_color_counts$colors_low)
# dev.off()
# # par vs hyp
# full_cis_parhyp_color_counts <- makeBipartiteMatrix(.top_colors = full_cis_colors_par, .bottom_colors = full_cis_colors_par)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTransBasis/hypotheticals/3B_full_cis_parhyp.pdf",
#     width = 4, height = 3)
# plotweb(full_cis_parhyp_color_counts$matrix, high.lablength = 0, low.lablength = 0,
#         col.high = full_cis_parhyp_color_counts$colors_high,
#         col.low = full_cis_parhyp_color_counts$colors_low)
# dev.off()
# # hyc vs hyp (identical to cer/par)
# full_cis_hychyp_color_counts <- makeBipartiteMatrix(.top_colors = full_cis_colors_cer, .bottom_colors = full_cis_colors_par)
# pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/CisTransBasis/hypotheticals/3B_full_cis_hychyp.pdf",
#     width = 4, height = 3)
# plotweb(full_cis_hychyp_color_counts$matrix, high.lablength = 0, low.lablength = 0,
#         col.high = full_cis_hychyp_color_counts$colors_high,
#         col.low = full_cis_hychyp_color_counts$colors_low)
# dev.off()
# 
# # 3C:   real bipartites but it's hybrid alleles now
# #       also cer-hyc and par-hyp bipartites because 
# #       you know people will be curious
# 
# # in supplemental bipartite comparisons folder, created in Fig1
# # Quantifying % module overlap for all 3 bipartites (like those shown in toy examples in 3A-B)
# # module preservation is supposed to be an intuitive measure of how similar 
# # module definitions are without envoking permutation tests yet
# # module preservation takes the color set with MORE colors as reference
# # then takes the largest CCM for each module in reference as numerator
# # and total gene number as denominator
# calculatePercentModuleOverlap <- function(.colors1, .colors2) {
#   if (length(.colors1) != length(.colors2)) {
#     stop("color vectors are not the same length!\n")
#   }
#   n1 <- .colors1 %>% unique() %>% length()
#   n2 <- .colors2 %>% unique() %>% length()
#   if (n1 >= n2) {
#     ref_colors <- .colors1
#     test_colors <- .colors2
#   }
#   if (n1 < n2) {
#     ref_colors <- .colors2
#     test_colors <- .colors1
#   }
#   numerator <- 0
#   for (color in unique(ref_colors)) {
#     largest_preserved <- test_colors[which(ref_colors == color)] %>% table() %>% max()
#     numerator <- numerator + largest_preserved
#   }
#   denominator <- length(.colors1)
#   return(numerator/denominator)
# }
# # tests for calculatePercentModuleOverlap
# calculatePercentModuleOverlap(full_trans_colors_cer, full_trans_colors_cer) # should be 1.00
# calculatePercentModuleOverlap(full_trans_colors_hyb, full_trans_colors_cer) # should be < 1.00
# calculatePercentModuleOverlap(full_trans_colors_hyb, full_trans_colors_par) # should be < 1.00
# calculatePercentModuleOverlap(full_cis_colors_cer, full_cis_colors_par) # should be < 1.00
# 
# # actual calculating on 3 real bipartites
# # cer-hyc
# calculatePercentModuleOverlap(colors10$cer, colors10$hyc)
# # par-hyp
# calculatePercentModuleOverlap(colors10$par, colors10$hyp)
# # hyc-hyp
# calculatePercentModuleOverlap(colors10$hyc, colors10$hyp)
# # cer-par (not in Fig 3, ust for curiosity)
# calculatePercentModuleOverlap(colors10$cer, colors10$par)
# 
# # 3D:  Classifying parental CCMs as preserved or not preserved in hybrid
# 
# 
# # I wanted to do a bipartite visualization but 
# # A) we have too many and 
# # B) bipartites have thus far been ALL genes and this would filter for CCM genes which might be confusing
# 
# # Bipartite code (archived for now)
# # # Alright we def need to filter for CCMs that have at least min_CCM_size genes, in parents OR hybrids 
# # CCM_names_parents <- table(paste(colors10$cer, colors10$par))[which(table(paste(colors10$cer, colors10$par)) > min_CCM_size)] %>% names()
# # CCM_names_parents
# # gene_is_in_parental_CCM <- paste(colors10$cer, colors10$par) %in% CCM_names_parents
# # CCM_names_hybrid <- table(paste(colors10$hyc, colors10$hyp))[which(table(paste(colors10$hyc, colors10$hyp)) > min_CCM_size)] %>% names()
# # CCM_names_hybrid
# # gene_is_in_hybrid_CCM <- paste(colors10$hyc, colors10$hyp) %in% CCM_names_hybrid
# # sum(gene_is_in_parental_CCM | gene_is_in_hybrid_CCM)/length(gene_is_in_parental_CCM) # how much is this filtering anything out
# # 
# # # here we just pretend the combination of cer and par colors is the new set of 45 or so colors
# # parent_colors <- paste(colors10$cer[gene_is_in_parental_CCM | gene_is_in_hybrid_CCM], colors10$par[gene_is_in_parental_CCM | gene_is_in_hybrid_CCM])
# # hybrid_colors <- paste(colors10$hyc[gene_is_in_parental_CCM | gene_is_in_hybrid_CCM], colors10$hyp[gene_is_in_parental_CCM | gene_is_in_hybrid_CCM])
# # # Bipartite comparisons
# # CCM_color_counts <- makeBipartiteMatrix(.top_colors = parent_colors, .bottom_colors = hybrid_colors)
# # # aaaand this is a mess
# # plotweb(CCM_color_counts$matrix, high.lablength = 0, low.lablength = 0)
# # # TODO: The discrepancy appears to be here. There should only be 11 "color names" and yet there are 93...
# # # is it just from including all the hybrid options?
# # parent_colors %>% unique() %>% length()
# # table(paste(colors10$cer, colors10$par))[which(table(paste(colors10$cer, colors10$par)) > min_CCM_size)] %>% length()
# # # If a column or row has a sum less than min_CCM_size, it must be connected to (aka have a non-zero entry with) a column that does sum to > min_CCM_size
# # test <- sample(which(colSums(CCM_color_counts$matrix) < min_CCM_size), 1)
# # test_connections <- which(CCM_color_counts$matrix[,names(test)] > 0)
# # (rowSums(CCM_color_counts$matrix[names(test_connections),]) > min_CCM_size) %>% all() 
# # # in fact it appears to be that every single connection it has has to be a CCM... Why is that?
# # # If it were connected to a non_CCM, then that would mean there's at least one gene (the connection) that's not in a CCM in either species, that's why
# 
# # 3E: quantification (maybe some sort of pie chart) of 3D
# 
# # 3F: Environmental expression of CCMs that are preserved in hybrid versus ones that are not
# # Basically same as 2C but add hybrid data for genes in parental CCMs, and indicate whether hybrid is mostly conserved or not
# # IDK we don't want to plot the same data twice, but
# # we also don't want to bring in the hybrid too early or have nothing to compare it to in 3C
# 
# 
# 
# 
# 
# #### Probably won't use: Quantification of cis/trans expectation in bipartites ####
# 
# # TODO: I stopped after supplementary upset plot cause tackling a case example felt more meaningful currently,
# #       but I may come back to quantify across modules


#### Wittkopp plot ####
# # Fig4: Wittkopp plot of expression effect/mean_expr in parents versus hybrid alleles (effect of par allele versus cer, so effect of 30 means Par is on average 30 counts higher than the comparable cer sample)
# pct_cis_quants <- spalcomps %>% filter(allele_sig | species_sig) %>% select(pct_cis) %>% quantile(na.rm = TRUE, probs = seq(0, 1, 0.1))
# 
# # TODO: import common names for all genes in spalcomps and convert plot_name to common names
# plotdf <- filter(spalcomps, allele_sig | species_sig) %>% 
#   mutate(plot_name = if_else(abs(species_estimate) > 1000 | abs(allele_estimate) > 1000,
#                              true = gene_name, false = ""),
#          pct_cis_bin = cut(pct_cis, breaks = pct_cis_quants, include.lowest = TRUE))
# 
# # median of condition-matched expression ratios (cer/par)
# plotdf$parental_expr <- sapply(plotdf$gene_name, calculateMedianExprRatio, .mode = "parents")
# plotdf$hybrid_expr <- sapply(plotdf$gene_name, calculateMedianExprRatio, .mode = "hybrid")
# # full plot
# ggplot(filter(plotdf, !is.na(pct_cis)), aes(x = log(parental_expr + 1e-5), y = log(hybrid_expr + 1e-5), label = gene_name)) + 
#   geom_text(check_overlap = TRUE, hjust = 0, nudge_x = 0.1, size = 3, aes(color = pct_cis_bin)) +
#   geom_abline(color = "gold") + geom_hline(yintercept = 0, color = "skyblue") + 
#   geom_point(aes(color = pct_cis_bin)) + theme_classic() + xlim(c(-6, 6)) + ylim(-6, 6) +
#   ggtitle("parental vs hybrid median\nexpression ratios (cer/par)") +
#   xlab("log(parental)") + ylab("log(hybrid)") + theme(legend.title = element_blank())
# 

