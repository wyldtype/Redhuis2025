#### Network analysis of Barkai Data ####
# Following WGCNA tutorials: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Following Jeremy Miller's tutorial on Meta-analysis with WGCNA: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/
sapply(c("tidyverse", "edgeR", "WGCNA", "matrixStats"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Barkai_data_analysis")
options(stringsAsFactors = FALSE)
# read in network data
load("data_files/Networks.RData")

# this was the one I actually used in my ASBMB poster
# load("data_files/NetworksHybSplit_unmerged.RData")

#### Comparison to networks constructed via other methods (perturbation, PPI, etc.) ####

# TODO: Here would be a good place to add a comparison of our network connections with 
# networks constructed from other methods (Costanzo, Landry lab, etc.)

#### Module preservation between species ####
# Evaluating how strongly each module is preserved between species

# visualizing module preservation for each species network with all 4 species' module definitions
for (i in 1:3) {
  plotDendroAndColors(dendro = geneTrees[[i]], colors = cbind(colors[[1]], colors[[2]], colors[[3]]),
                      groupLabels = names(colors), main = names(geneTrees)[i],
                      dendroLabels = FALSE)
}

#### Comparing Module Eigengenes ####
# Next I'll compare modules based on eigengenes to better understand how the modules that aren't shared in all 4 relate to these "core" modules
# recalculating module eigengenes because the colors switched during matchLabels
MEs <- map2(counts_top3000, colors, function(x, y) {
  MEoutput <- moduleEigengenes(x, colors = y)
  MEs <- MEoutput$eigengenes
  return(MEs)
})
MEDiss <- lapply(MEs, function(x) {
  return(1-abs(cor(x, use = "p")))  # Anna changed this to absolute value to match our other distance function. It probably won't make much of a difference, I doubt many module have negative correlations
})
METree <- lapply(MEDiss, function(x) {
  return(hclust(as.dist(x), method = "average"))
})

# MDS plot of eigeingenes
MDSs <- lapply(MEDiss, function(x) {
  MDScoords <- cmdscale(as.dist(x), k=2)
  MDStibble <- tibble(module = gsub("ME", "", rownames(MDScoords)),
                      x = MDScoords[,1],
                      y = MDScoords[,2])
  return(MDStibble)
}) # cmdscale projects our modules onto two dimensions where their distance from each other in 2D is approximately equal to their distances in dist
library(gridExtra)
MDSplots <- vector(mode = "list", length = 3)
SpeciesNames <- c("S. cerevisiae", "S. paradoxus", "F1 hybrid,\ncerevisiae allele", "F1 hybrid,\nparadoxus allele")
for (i in 1:3) {
  MDSplots[[i]] <- ggplot(data = MDSs[[i]], aes(x = x, y = y)) + 
    geom_point(color = MDSs[[i]]$module, size = 8) + 
    ggtitle(SpeciesNames[i]) + xlab("coordinate 1") + ylab("coordinate 2") + theme_classic(base_size = 21)
}
grid.arrange(grobs = MDSplots, nrow = 1) # Note that we'll never be able to overlay all 3 species because their eigengenes are different lengths (and don't even THINK about trying to make them the same length!)
# pretty cool, the modules common to all3 are all on their own. Doesn't seem to have much to do with preservation (blue had much lower preservation for hybrids than parents versus turquoise, but that isn't depicted here) though, I could look into why that is

#### kME Data Cleaning ####
# calculate kME, correlation of each gene with each module eigengene
kMEs <- map2(counts_top3000, MEs, signedKME)

# moduleEigengenes tries to align sign (+ or -) between MEs in the most useful direction, but it uses average expression of ALL genes in that module across environments, which I feel is often a fairly flat line and can therefore be wrong sometimes
# After calculating kME, we can compare kME vectors for each module between species (can't compare MEs directly between species) to see if any signs should be flipped
ParentColors <- intersect(unique(colors$cer), unique(colors$par))
dev.off()
# cer vs par
for (clr in ParentColors) {
  plot(x = kMEs$cer[, paste0("kME", clr)], # We're wanting kMEs to be positively correlated
       y = kMEs$par[, paste0("kME", clr)], 
       main = clr, xlab = "kME cer", ylab = "kME par", xlim = c(-1.0,1.0), ylim = c(-1.0,1.0))
} # There are two groups, one negative one positive correlation with eigengene, especially for turquoise module. This isn't seen in the meta-analysis tutorial, I'm guessing b/c that's a signed (0.5*(1+corr))^ß network and we did unsigned |corr|^ß
dev.off()
# cer vs hyb
for (clr in ParentColors) {
  plot(x = kMEs$cer[, paste0("kME", clr)], # We're wanting kMEs to be positively correlated
       y = kMEs$hyb[, paste0("kME", clr)], 
       main = clr, xlab = "kME cer", ylab = "kME hyb", xlim = c(-1.0,1.0), ylim = c(-1.0,1.0))
} 
dev.off()
# par vs hyb
for (clr in ParentColors) {
  plot(x = kMEs$par[, paste0("kME", clr)], # We're wanting kMEs to be positively correlated
       y = kMEs$hyb[, paste0("kME", clr)], 
       main = clr, xlab = "kME par", ylab = "kME hyb", xlim = c(-1.0,1.0), ylim = c(-1.0,1.0))
} 
dev.off()

# if you need to flip any kMEs, use this:
# dev.off()
# plot(x = kMEs$par$kMEturquoise,
#        y = kMEs$hyb$kMEturquoise,
#        main = "turquoise", xlab = "kME par", ylab = "kME hyp", xlim = c(-1.0,1.0), ylim = c(-1.0,1.0))
# # flipping MEs
# MEs$par$MEturquoise <- -MEs$par$MEturquoise
# MEs$hyb$MEturquoise <- -MEs$hyb$MEturquoise

#recalculating kME
kMEs <- map2(counts_top3000, MEs, signedKME) # you can re-plot those ^ to verify it fixed the turquoise

#### kME Data Exploration ####
# now we can meaningfully select the top 10 (strongest positive correlation with ME) and bottom 10 (strongest negative correlation with ME) kMEs
# note that I do nothing to ensure that the top/bottom 10 are actually in the module yet
top10_kMEs <- lapply(kMEs, function(mat) {
  top10perModule <- apply(mat, 2, function(mod) {
    return(rownames(mat)[rank(-mod) <= 10])
  })
  return(as_tibble(top10perModule))
})
bottom10_kMEs <- lapply(kMEs, function(mat) {
  bottom10perModule <- apply(mat, 2, function(mod) {
    return(rownames(mat)[rank(mod) <= 10]) # one - sign of difference
  })
  return(as_tibble(bottom10perModule))
})
# checking which module each of these top10 genes is in (should be its column's module)
GeneNames <- colnames(counts_top3000$cer)
top10_colors <- lapply(c(1:3), function(i) {
  tib <- top10_kMEs[[i]]
  cat(dim(tib))
  tib <- apply(tib, c(1,2), function(x) {
    x_idx <- which(GeneNames == as.character(x))
    return(colors[[i]][x_idx])
  })
  return(tib)
})
for (i in c(1:3)) {
  print(top10_colors[[i]])
} # we'd expect top10 kMEs to be in their actual module, except grey which is no man's land

bottom10_colors <- lapply(c(1:3), function(i) {
  tib <- bottom10_kMEs[[i]]
  tib <- apply(tib, c(1,2), function(x) {
    x_idx <- which(GeneNames == as.character(x))
    return(colors[[i]][x_idx])
  })
  return(tib)
})
for (i in c(1:3)) {
  print(bottom10_colors[[i]])
} # hmm this is telling me that the negative correlations are less meaningful --- the vast majority of bottom10 kMEs are from the turquoise module of
# each network, aka the largest module and more likely to just happen to have negative correlation with the ME of interest

# so this is to say that we'll just work with top10 sets from here on

#### Introducing the Gene Df --- 1 row = 1 gene ####
# recalculating kME to just give one value for each gene, it's in-module kME
GeneNames <- colnames(counts_top3000$cer)
gene_dfs <- lapply(colors, function(x) {
  tib <- tibble(GeneName = GeneNames,
                module = x)
})
# in-module kME
for (i in 1:3) {
  kME_names <- names(kMEs[[i]])
  gene_dfs[[i]]$inModuleKME <- sapply(1:length(colors[[i]]), function(idx) {
    modcol <- colors[[i]][idx]
    whichKME <- which(kME_names == paste0("kME", modcol))
    return(kMEs[[i]][idx, whichKME])
  })
}
# recalculating top10 kMEs to be module-specific
top10InModule <- lapply(gene_dfs, function(x) {
  top10List <- vector(mode = "list", length = length(unique(x$module)))
  names(top10List) <- unique(x$module)
  for (modcol in unique(x$module)) {
    mod_df <- filter(x, module == modcol)
    top10 <- mod_df$GeneName[rank(-mod_df$inModuleKME) <= 10]
    top10List[[modcol]] <- top10
  }
  return(top10List)
})
for (i in 1:3) {
  gene_dfs[[i]]$isTop10 <- pmap_lgl(gene_dfs[[i]], function(GeneName, module, inModuleKME) {
    return(GeneName %in% top10InModule[[i]][[module]])
  })
}

# A final check of hub-gene-ness is to see if top10 genes have above-average connectivity (sum of adjacencies to all other genes) within their respective modules
# I won't bother to do this for just the intersect top10s unless the full top10s look weird
for (i in 1:3) {
  gene_dfs[[i]]$connectivity <- colSums(adjs[[i]])
}

connVskMEPlots <- vector(mode = "list", length = 3)
for (i in 1:3) {
  p <- ggplot(data = gene_dfs[[i]], aes(x = paste(module, isTop10), y = connectivity)) + 
    geom_point() +
    theme(axis.text.x = element_text(angle = 90))
  connVskMEPlots[[i]] <- p
}
grid.arrange(grobs = connVskMEPlots, nrow = 2) # def on the higher end of connectivity, but not the absolute max. Makes sense because there's no direct relationship between connectivity and correlation with ME

# TODO: in-module connectivity would probably be a more accurate estimation of hub gene-ness, but that's a lot of work for a Saturday

# Gene Ontology
# crete lists of which genes are in each module in each species
# SGD needs a .txt file where every line is a gene name, so we'll do one for each module in each species
# for just cer/par:
for (i in 1:2) {
  for (modcol in unique(gene_dfs[[i]]$module)) {
    mod_names <- filter(gene_dfs[[i]], module == modcol)$GeneName
    write.table(mod_names, file = paste0("GO_unmergedMods_CerPar/", modcol, "_", names(gene_dfs)[i], ".txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

# here you go to SGD and input .txt files and output info either by hand or as tables

# read in cer GO processes from SGD
GOprocessesCer <- vector(mode = "list", length = (length(unique(colors$par)) - 1))
names(GOprocessesCer) <- setdiff(unique(colors$cer), "pink")
for (clr in names(GOfunctionsCer)) {
  GOprocessesCer[[clr]] <- read.table(file = paste0("GO_unmergedMods_CerPar/", clr, "_cer_process_tab.txt"), sep = "\t", header = TRUE)
}

# I decided to do this recursively for some reason
# given a data frame of GO terms from SGD, returns a list of the minimum terms (starting at row 1) that reference all the genes referenced in the set
# so a lazy way to get a sense of the minimum set of terms but without trying all possible combinations of terms (I can do that eventually if I like)
# It's still too long for including in the data, but helpful to look at by eye
getMinGOTerms <- function(go_df, remaining_genes = NULL, genes_present = character(0), min_terms = character(0), idx = 1) {
  if (length(remaining_genes) == 0 & !is.null(remaining_genes)) {
    return(list(terms = min_terms, genes = genes_present, max_iter = idx))
  }
  if (is.null(remaining_genes)) {
    remaining_genes <- sapply(unlist(go_df$ANNOTATED_GENES), str_split, pattern = ", ") %>% unlist() %>% unique()
  }
  min_terms <- c(min_terms, go_df$TERM[idx])
  genes_i <- str_split(go_df$ANNOTATED_GENES[idx], pattern = ", ") %>% unlist()
  if (!all(genes_i %in% genes_present)) {
    remaining_genes <- setdiff(remaining_genes, genes_i)
    genes_present <- unique(c(genes_present, genes_i))
  }
  getMinGOTerms(go_df, remaining_genes = remaining_genes, genes_present = genes_present, min_terms = min_terms, idx = idx + 1)
}
test <- getMinGOTerms(read.table(file = paste0("GO_unmergedMods_CerPar/", "blue", "_cer_process_tab.txt"), sep = "\t", header = TRUE))

# pink isn't big enough for statistical power, but it's our spo module by manual inspection
purrr::reduce(list(GOlists$cer$gene_name[GOlists$cer$module == "pink"],
          GOlists$par$gene_name[GOlists$par$module == "pink"],
          GOlists$hyc$gene_name[GOlists$hyc$module == "pink"],
          GOlists$hyp$gene_name[GOlists$hyp$module == "pink"]), intersect)
# of the 17 shared pink genes (of 19 total pink genes): 10 spo, 3 nicotinic, 1 low-glucose conditions hexose transporter, 3 unknown

# greenyellow par not sig either but
# mostly proteins that localize to mitochondria and upregulate in response to DNA/oxidative stress

# tan par not sig either but mostly membrane proteins:
# 1 stationary phase, 2 vacuolar 1 autophagosomes 1 glutamate dehydrogenase 2 drug resistant 1 ubiquitin conjugating 1 glyoxylate aminotransferase 1 glucan synthase 1 cell wall 2 in lipid particles in ER

GOProcessSummary <- list(cer = tibble("turquoise" = "ribosome biogenesis",
                                      "blue" = "oxidative phosphorylation",
                                      "pink" = "sporulation",
                                      "brown" = "pyridine compound metabolism\nheat acclimation",
                                      "green" = "organic acid transport, \nresponse to stimulus",
                                      "red" = "amino acid biosynthesis",
                                      "yellow" = "sulfur compound metabolism and biosynthesis",
                                      "grey" = "no significant results",),
                            par = tibble("turquoise" = "ribosome biogenesis",
                                         "blue" = "oxidative phosphorylation",
                                         "pink" = "sporulation",
                                         "brown" = "carbohydrate metabolism and biosynthesis",
                                         "green" = "organic acid transport, \nresponse to stimulus",
                                         "red" = "amino acid biosynthesis",
                                         "greenyellow" = "localize to mitochondria,\nstress response",
                                         "salmon" = "translation",
                                         "black" = "glycolitic pathway",
                                         "tan" = "membrane proteins",
                                         "grey" = "no significant results"))

library(huxtable)
# cer GO process
ht_cer <- hux(GOProcessSummary$cer) %>%
  set_font_size(20)                               %>%
  set_header_rows(1, TRUE)                             %>%
  set_width(0.9)                                       %>%
  set_tb_padding(2)                                    %>%
  set_caption("Huxtable properties")                   %>%
  set_label("tab:props")                               %>%
  theme_bright(
    header_rows = TRUE,
    header_cols = FALSE,
    colors = names(GOProcessSummary$cer)
  )

font(ht_cer) <- "monospace"
caption(ht_cer) <- "<span style='font-family: monospace'>GO Enrichment Summary - S. cerevisiae</span>"

quick_html(ht_cer, file = "GOprocessCer.html") # TODO: Add summary of each module's role (ex) Blue is "Respiration Module"), and hub gene descriptions for each

# par GO process
ht_par <- hux(GOProcessSummary$par) %>%
  set_font_size(20)                               %>%
  set_header_rows(1, TRUE)                             %>%
  set_width(0.9)                                       %>%
  set_tb_padding(2)                                    %>%
  set_caption("Huxtable properties")                   %>%
  set_label("tab:props")                               %>%
  theme_bright(
    header_rows = TRUE,
    header_cols = FALSE,
    colors = names(GOProcessSummary$par)
  )

font(ht_par) <- "monospace"
caption(ht_par) <- "<span style='font-family: monospace'>GO Enrichment Summary - S. paradoxus</span>"

quick_html(ht_par, file = "GOprocessPar.html") # TODO: Add summary of each module's role (ex) Blue is "Respiration Module"), and hub gene descriptions for each

# Hub genes
HubGeneList <- lapply(gene_dfs, function(x) {
  clrs <- unique(x$module)
  hubs <- vector(mode = "list", length = length(clrs))
  for (i in 1:length(clrs)) {
    mod_df <- filter(x, module == clrs[i])
    idx <- which(mod_df$inModuleKME == max(mod_df$inModuleKME))
    hubs[[i]] <- mod_df$GeneName[idx]
  }
  names(hubs) <- clrs
  return(hubs)
})
# Description of the top hub gene (highest |kME|) for each module for cer and par
HubGenesInfo <- list(cer = tibble("turquoise" = ,
                                  "blue" = ,
                                  "grey" = ,
                                  "yellow" = ,
                                  "pink" = ,
                                  "brown" = ,
                                  "green" = ,
                                  "red" = ),
                     par = tibble("turquoise" = ,
                                  "blue" = ,
                                  "grey" = ,
                                  "greenyellow" = ,
                                  "pink" = ,
                                  "brown" = ,
                                  "green" = ,
                                  "red" = ,
                                  "salmon" = ,
                                  "black" = ,
                                  "tan" = ,))

#### Identifying genes that have switched modules ####

DidSwitchParents <- (colors$cer != colors$par) & (colors$cer != "grey") & (colors$par != "grey") & (colors$hyc != "grey") & (colors$hyp != "grey") # grey doesn't count b/c it's not really a module
DidSwitchHybrids <- (colors$hyc != colors$hyp) & (colors$cer != "grey") & (colors$par != "grey") & (colors$hyc != "grey") & (colors$hyp != "grey")
DidSwitchHybridsAndParents <- DidSwitchParents & DidSwitchHybrids # cis influence on module switching
DidSwitchParentsNotHybrids <- DidSwitchParents & !DidSwitchHybrids # trans influence on module switching
DidSwitchHybridsNotParents <- !DidSwitchParents & DidSwitchHybrids # This is weird, shouldn't be too many of these

# when a gene does switch in parents and hybrids, does it switch to same module?
sum(colors$cer[DidSwitchHybridsAndParents] == colors$hyc[DidSwitchHybridsAndParents])/sum(DidSwitchHybridsAndParents) # usually, not always

# which modules have the most switching?
# pie charts for starters (unfortunately there's some margin issue with fitting 12 on a grid...)
# par(mfrow = c(4,3))
# mtext(c("Parents and Hybrids", "Parents Not Hybrids", "Hybrids Not Parents"))
for (clrs in colors) {
  pie(table(clrs[DidSwitchHybridsAndParents]), col = rownames(table(clrs[DidSwitchHybridsAndParents])), labels = "")
  pie(table(clrs[DidSwitchParentsNotHybrids]), col = rownames(table(clrs[DidSwitchParentsNotHybrids])), labels = "")
  pie(table(clrs[DidSwitchHybridsNotParents]), col = rownames(table(clrs[DidSwitchHybridsNotParents])), labels = "")
}

# bipartite view of module switching
library(bipartite)
# parents
ModuleNames_cer <- setdiff(unique(colors$cer), "grey")
ModuleNames_par <- setdiff(unique(colors$par), "grey")
nColor_cer <- length(ModuleNames_cer)
nColor_par <- length(ModuleNames_par)
parents_color_counts <- matrix(nrow = nColor_cer, ncol = nColor_par, dimnames = list(ModuleNames_cer, ModuleNames_par))
for (col_cer in ModuleNames_cer) {
  for (col_par in ModuleNames_par) {
    parents_color_counts[col_cer, col_par] <- sum((colors$cer == col_cer) & (colors$par == col_par))
  }
}
plotweb(parents_color_counts, col.high = ModuleNames_par, col.low = ModuleNames_cer, labsize = 1, text.rot = 45) # oh my goodness would but everything be this easy <3 <3 <3

# cerevisiae parent & hybrid
ModuleNames_hyb <- setdiff(unique(colors$hyb), "grey")
nColor_hyb <- length(ModuleNames_hyb)
cerhyb_color_counts <- matrix(nrow = nColor_hyb, ncol = nColor_cer, dimnames = list(ModuleNames_hyb, ModuleNames_cer))
for (col_hyb in ModuleNames_hyb) {
  for (col_cer in ModuleNames_cer) {
    cerhyb_color_counts[col_hyb, col_cer] <- sum((colors$hyb == col_hyb) & (colors$cer == col_cer))
  }
}
plotweb(cerhyb_color_counts, col.high = ModuleNames_cer, col.low = ModuleNames_hyb, labsize = 1)

# paradoxus parent & hybrid
parhyb_color_counts <- matrix(nrow = nColor_hyb, ncol = nColor_par, dimnames = list(ModuleNames_hyb, ModuleNames_par))
for (col_hyb in ModuleNames_hyb) {
  for (col_par in ModuleNames_par) {
    parhyb_color_counts[col_hyb, col_par] <- sum((colors$hyb == col_hyb) & (colors$par == col_par))
  }
}
plotweb(parhyb_color_counts, col.high = ModuleNames_par, col.low = ModuleNames_hyb, labsize = 1)

### divercence mode
load(file = "data_files/cis_trans_df.RData")
# what are the most common switches that happen between cer and par modules?
table(paste(gene_dfs$cer$module, gene_dfs$par$module, sep = " to ")) %>% sort(decreasing = TRUE)
# TODO: join spaldf and the cer and par module definitions (and hyb?) and see if cis or trans-diverging genes are enriched in certain pairs of modules (i.e. blue to turquoise)

#### Is module switching predicted by connectivity, kME, or expression variability (tau)? ####
# module switching
gene_dfs$cer$switchedBetweenGroups <- colors$cer != colors$par # groups: parents or hybrids
gene_dfs$cer$switchedBetweenAlleles <- colors$cer != colors$hyc # alleles: cer or par
gene_dfs$par$switchedBetweenGroups <- colors$cer != colors$par
gene_dfs$par$switchedBetweenAlleles <- colors$par != colors$hyp
gene_dfs$hyc$switchedBetweenGroups <- colors$hyc != colors$hyp
gene_dfs$hyc$switchedBetweenAlleles <- colors$hyc != colors$cer
gene_dfs$hyp$switchedBetweenGroups <- colors$hyp != colors$hyc
gene_dfs$hyp$switchedBetweenAlleles <- colors$hyp != colors$par

# 1) Tau (expression specificity)

# calculates Tau
# input: counts data frame where rows are genes, columns are environments/samples
# output: vector of length = # of genes, where each gene has its corresponding Tau value
calculate_tau <- function(cts) {
  m <- ncol(cts)
  max_per_gene <- apply(cts, 1, max)
  unscaled_tau <- (1 - cts/max_per_gene) %>% rowSums(na.rm = TRUE)
  return(unscaled_tau/(m - 1))
}

for (i in 1:3) {
  gene_dfs[[i]]$tau <- calculate_tau(t(counts_top3000[[i]]))
}

# Tau by module + tau of MEs for each module
library(ggbeeswarm)
tau_by_module_plots <- vector(mode = "list", length = 3)
for (i in 1:3) {
  p <- ggplot(gene_dfs[[i]], aes(x = module, y = tau), alpha = 0.5) +
    geom_quasirandom(aes(color = module)) +
    # TODO: figure out how to use scale x discrete without hand-typing every single color to get the modules to be their actual colors
    # scale_x_discrete() + 
    ggtitle(SpeciesNames[i]) + ylab("Expression \nSpecificity") + xlab("") + theme(legend.position = "none")
  tau_by_module_plots[[i]] <- p
}
library(ggpubr)
ggarrange(tau_by_module_plots[[1]],
          tau_by_module_plots[[2]],
          tau_by_module_plots[[3]], nrow = 3)

# switching probability vs expression variability
switch_vs_tau_plots <- vector(mode = "list", length = 8)
for (i in 1:4) {
  p1 <- ggplot(gene_dfs[[i]], aes(x = switchedBetweenGroups, y = tau, color = switchedBetweenGroups), alpha = 0.5) +
    geom_quasirandom() + ggtitle("Switched Between Groups") + xlab(names(gene_dfs)[i]) + ylab("Expression Variability (Tau)")
  p2 <- ggplot(gene_dfs[[i]], aes(x = switchedBetweenAlleles, y = tau, color = switchedBetweenAlleles), alpha = 0.5) +
    geom_quasirandom() + ggtitle("Switched Between Alleles") + xlab(names(gene_dfs)[i]) + ylab("Expression Variability (Tau)")
  switch_vs_tau_plots[[2*i - 1]] <- p1
  switch_vs_tau_plots[[2*i]] <- p2
}

ggarrange(switch_vs_tau_plots[[1]], switch_vs_tau_plots[[2]], 
          switch_vs_tau_plots[[3]], switch_vs_tau_plots[[4]],
          switch_vs_tau_plots[[5]], switch_vs_tau_plots[[6]],
          switch_vs_tau_plots[[7]], switch_vs_tau_plots[[8]], ncol=2, nrow=4, common.legend = TRUE, legend="bottom")
# can't see much there ^

# again just looking at parents
parents_tau_vs_switch <- bind_rows(tibble(tau = gene_dfs$cer$tau, 
                                          switched = gene_dfs$cer$switchedBetweenGroups,
                                          species = "cer"), 
                                   tibble(tau = gene_dfs$par$tau, 
                                          switched = gene_dfs$par$switchedBetweenGroups, 
                                          species = "par"))
wilcox_cer <- wilcox.test(gene_dfs$cer$tau[gene_dfs$cer$switchedBetweenGroups], gene_dfs$cer$tau[!gene_dfs$cer$switchedBetweenGroups])
wilcox_par <- wilcox.test(gene_dfs$par$tau[gene_dfs$par$switchedBetweenGroups], gene_dfs$par$tau[!gene_dfs$par$switchedBetweenGroups])
ggplot(data = parents_tau_vs_switch) +
  geom_quasirandom(aes(x = paste(species, switched), y = tau, color = switched)) + scale_x_discrete(breaks = c("Hai", "der", "U", "cutie"))
mean(gene_dfs$cer$tau[gene_dfs$cer$switchedBetweenGroups])
mean(gene_dfs$cer$tau[!gene_dfs$cer$switchedBetweenGroups]) # you know it's significant when you can't tell which way it goes...
mean(gene_dfs$par$tau[gene_dfs$par$switchedBetweenGroups])
mean(gene_dfs$par$tau[!gene_dfs$par$switchedBetweenGroups]) 
# I'm not sure I believe that cer wilcox test, let's repeat with downsampled data
random500 <- sample(1:3000, 500)
wilcox_cer <- wilcox.test(gene_dfs$cer[random500,]$tau[gene_dfs$cer$switchedBetweenGroups], gene_dfs$cer[random500,]$tau[!gene_dfs$cer$switchedBetweenGroups])
# yeah def just the power

# let's look at the difference in tau seeing as switching is the same information twice
ggplot(data = gene_dfs$cer) + geom_quasirandom(aes(x = switchedBetweenGroups, y = abs(gene_dfs$cer$tau - gene_dfs$par$tau), color = switchedBetweenGroups))
wilcox.test(abs(gene_dfs$cer$tau - gene_dfs$par$tau)[gene_dfs$cer$switchedBetweenGroups], abs(gene_dfs$cer$tau - gene_dfs$par$tau)[!gene_dfs$cer$switchedBetweenGroups]) # nope

# parents tau boxplots
justTheTwoOfUs <- tibble("Average Expression \nSpecificity" = pmean(gene_dfs$cer$tau, gene_dfs$par$tau), 
                         "Switched Modules" = gene_dfs$cer$switchedBetweenGroups)
ggboxplot(data = justTheTwoOfUs, x = "Switched Modules", 
          y = "Average Expression \nSpecificity", color = "Switched Modules", palette = "jco") + 
  stat_compare_means(label.x = 1, label.y = 1.1) + theme(text=element_text(size=21), legend.position="none")
sum(justTheTwoOfUs$`Switched Modules`)

# Tau vs kME just for fun
tau_vs_kME_plots <- vector(mode = "list", length = 4)
for (i in 1:4) {
  tau_vs_kME_plots[[i]] <- ggplot(data = gene_dfs[[i]], aes(x = abs(inModuleKME), y = tau)) + geom_point(aes(color = switchedBetweenGroups))
}
ggarrange(tau_vs_kME_plots[[1]],
          tau_vs_kME_plots[[2]],
          tau_vs_kME_plots[[3]],
          tau_vs_kME_plots[[4]], nrow = 2, ncol = 2, common.legend = TRUE) # def no relationship there!


# connectivity by module
conn_by_module_plots <- vector(mode = "list", length = 4)
for (i in 1:4) {
  p <- ggplot(gene_dfs[[i]], aes(x = module, y = connectivity), alpha = 0.5) + # abs value kME feels the most useful, strong negative correlations with ME should still indicate that gene is important to that module
    geom_quasirandom(color = gene_dfs[[i]]$module) + ggtitle(SpeciesNames[i]) + ylab("Connectivity") + xlab("") +
    scale_x_discrete(labels = c("blue" = "blue", "brown" = "brown", "green" = "green", "turquoise" = "turquoise", "red" = "red", "pink" = "pink", "yellow" = "yellow", "black" = "black", "salmon" = "salmon", "greenyellow" = "greenyellow", "tan" = "tan", "grey" = "grey\n(genes w/o module)"),
                     limits = c("blue", "brown", "green", "turquoise", "red", "pink", "yellow", "black", "salmon", "greenyellow", "tan" = "tan", "grey"))
  conn_by_module_plots[[i]] <- p
}
ggarrange(conn_by_module_plots[[1]],
          conn_by_module_plots[[2]],
          conn_by_module_plots[[3]],
          conn_by_module_plots[[4]], nrow = 4, common.legend = TRUE)

# Connectivity vs switching
conn_vs_switch_plots <- vector(mode = "list", length = 4)
n_observations <- numeric()
for (i in 1:4) {
  p <- ggboxplot(data = gene_dfs[[i]], x = "switchedBetweenGroups", 
                 y = "connectivity", color = "switchedBetweenGroups", palette = "jco") + 
    stat_compare_means(label.x = 1, label.y = 250) + scale_x_discrete(breaks = c("hai", "der")) + xlab("")
  conn_vs_switch_plots[[i]] <- ggpar(p, legend.title = "Switched modules between groups or alleles", main = SpeciesNames[i])
  n_observations <- c(n_observations, sum(gene_dfs[[i]]$switchedBetweenGroups), sum(!gene_dfs[[i]]$switchedBetweenGroups))
}
n_observations # TODO: figure out how to add these n values to the boxplots. For now I'll just add them in Illustrator
ggarrange(conn_vs_switch_plots[[1]],
          conn_vs_switch_plots[[2]],
          conn_vs_switch_plots[[3]],
          conn_vs_switch_plots[[4]], nrow = 2, ncol = 2, common.legend = TRUE)

# one plot just for parents with average connectivity
justTheTwoOfUs <- tibble("Average Connectivity" = pmean(gene_dfs$cer$connectivity, gene_dfs$par$connectivity), 
                        "Switched Modules" = gene_dfs$cer$switchedBetweenGroups)
ggboxplot(data = justTheTwoOfUs, x = "Switched Modules", 
          y = "Average Connectivity", color = "Switched Modules", palette = "jco") + 
  stat_compare_means(label.x = 1, label.y = 250) + theme(text=element_text(size=21), legend.position="none")

# kME vs switching
switch_vs_kME_plots <- vector(mode = "list", length = 8)
for (i in 1:4) {
  p1 <- ggplot(gene_dfs[[i]], aes(x = switchedBetweenGroups, y = abs(inModulekME), color = switchedBetweenGroups), alpha = 0.5) + # abs value kME feels the most useful, strong negative correlations with ME should still indicate that gene is important to that module
    geom_quasirandom() + ggtitle("Switched Between Groups") + xlab(names(gene_dfs)[i]) + ylab("|kME|")
  p2 <- ggplot(gene_dfs[[i]], aes(x = switchedBetweenAlleles, y = abs(inModulekME), color = switchedBetweenAlleles), alpha = 0.5) +
    geom_quasirandom() + ggtitle("Switched Between Alleles") + xlab(names(gene_dfs)[i]) + ylab("|kME|")
  switch_vs_kME_plots[[2*i - 1]] <- p1
  switch_vs_kME_plots[[2*i]] <- p2
}
ggarrange(switch_vs_kME_plots[[1]], switch_vs_kME_plots[[2]], 
          switch_vs_kME_plots[[3]], switch_vs_kME_plots[[4]],
          switch_vs_kME_plots[[5]], switch_vs_kME_plots[[6]],
          switch_vs_kME_plots[[7]], switch_vs_kME_plots[[8]],ncol=2, nrow=4, common.legend = TRUE, legend="bottom")

justTheTwoOfUs <- tibble("Average Correlation\nWith Module" = pmean(abs(gene_dfs$cer$inModuleKME), abs(gene_dfs$par$inModuleKME)), 
                         "Switched Modules" = gene_dfs$cer$switchedBetweenGroups)
ggboxplot(data = justTheTwoOfUs, x = "Switched Modules", 
          y = "Average Correlation\nWith Module", color = "Switched Modules", palette = "jco") + 
  stat_compare_means(label.x = 1, label.y = 1.2) + theme(text=element_text(size=21), legend.position="none")


# let's just look at genes that have switched between parents in terms of their kME in both species
ggplot(data = bind_rows(tibble(kME = gene_dfs$cer$inModulekME, 
                               switched = gene_dfs$cer$switchedBetweenGroups,
                               species = "cer"), 
                        tibble(kME = gene_dfs$par$inModulekME, 
                               switched = gene_dfs$par$switchedBetweenGroups, 
                               species = "par"))) +
  geom_quasirandom(aes(x = paste(species, switched), y = kME, color = switched)) # this is mostly driven by the green-turquoise switch I'm guessing though. Maybe there's a more robust way to look at network rewiring besides if genes switched colors


#### Duplicates vs Singletons ####
marchant <- read.table("marchant_paralogs.txt", header = TRUE)

marchant$sameModuleCer <- sapply(c(1:nrow(marchant)), function(idx) {
  P1_idx <- which(gene_dfs$cer$gene_name == marchant[idx,]$P1)
  P2_idx <- which(gene_dfs$cer$gene_name == marchant[idx,]$P2)
  if (length(P1_idx) > 0 & length(P2_idx) > 0) {
    col_P1 <- colors$cer[P1_idx]
    col_P2 <- colors$cer[P2_idx]
    return(col_P1 == col_P2)
  }
  else {
    return(NA)
  }
}) %>% unlist()

sum(marchant[marchant$Duplication == "WGD",]$sameModuleCer, na.rm = TRUE)
sum(!marchant[marchant$Duplication == "WGD",]$sameModuleCer, na.rm = TRUE) # pretty perfect split, but I am missing like half my pairs due to one being dropped in the top3000 cutoff
sum(marchant[marchant$Duplication == "WGD" & marchant$Dupli.SSD_WGD == 0,]$sameModuleCer, na.rm = TRUE)

# Issue: many pairs had one member thrown out when we filtered for top 3000 genes, but this was fairly necessary b/c the networks
# do not look good if we don't filter. So we might need to focus more on a few pairs, or maybe it'll be enough to get the multi SSDs out of there


# TODO:
# comparing networks:
# connectivity by group: conserved, switched module in parents only, switched module in parents and hybrids
# overlay adj matrices
# ECCs from deBoer's Nature paper versus module switching (Supplementary Table 1)
# case studies for causative genetic change affecting module switching using barkai TF binding data
# Luis' idea: titrate expression between two species with different module membership until a module in one species but not the other appears ex) titrate purple module's paradoxus expression in cerevisiae until purple appears in cerevisiae
# Then to make Luis idea a little more realistic, find the minimum set of expression changes that shift one species' network to the other (ex) stopping criteria: 90% of genes are in same module in both species)
# other suggestion: switch genes in order of highest to lowest kME, so hubs switch first, to try and speed up the process
# determine better metric for divergence besides module switching

# QCish:
# for readability, create gene_dfs in this script with all necessary columns then move visualizations and analysis to new script: notwonks_genesDF.R
# label each gene as TF, enzyme, etc. and see what modules they're in
# repeat everything with different random sample
# differential expression vs switching? Or whichever other metrics for module divergence I come up with

# Hybrid stress buffering hypothesis:
# ID set of stress response genes
# Are stress genes up in par vs cer but not hyp vs hyc?
# Justin Fay looks at temp stress, his lab's work may have some ideas about mechanism/candidate genes

######################### Archive ##############################
# Skip Module Preservation, not clear it's what we want and appears to be biased by module size
# streamlining some of the annoyance of the WGCNA modulePreservation function (whyyy does it need lists of lists named data?)
modulePreservationAcrossSpecies <- function(cts_list, ref_colors, network_type) {
  nNetworks <- length(cts_list)
  cts <- list(ref = list(data = cts_list$ref))
  ref_colors <- list(ref = ref_colors)
  for (i in 2:nNetworks) {
    cts[[i]] <- list(data = cts_list[[i]])
    names(cts)[i] <- names(cts_list)[i]
  }
  mp <- modulePreservation(cts, ref_colors, 
                           referenceNetworks = 1, verbose = 3, 
                           networkType = network_type, nPermutations = 30, 
                           maxGoldModuleSize = 100, maxModuleSize = 400) # This number is 400 into demos, modules larger are randomly downsampled
  return(mp)
}

# too complicated to do this in a loop so here goes
mps <- list(cer = NULL,
            par = NULL,
            hyb = NULL)
# cer
mps$cer <- modulePreservationAcrossSpecies(cts_list = list(ref = counts_top3000$cer,
                                                           par = counts_top3000$par,
                                                           hyb = counts_top3000$hyb), 
                                           ref_colors = colors$cer,
                                           network_type = "unsigned")
# par
mps$par <- modulePreservationAcrossSpecies(cts_list = list(ref = counts_top3000$par,
                                                           cer = counts_top3000$cer,
                                                           hyb = counts_top3000$hyb), 
                                           ref_colors = colors$par,
                                           network_type = "unsigned")
# hyb
mps$hyb <- modulePreservationAcrossSpecies(cts_list = list(ref = counts_top3000$hyb,
                                                           cer = counts_top3000$cer,
                                                           par = counts_top3000$par), 
                                           ref_colors = colors$hyc,
                                           network_type = "unsigned")


# Save these files b/c their computation is hefty
save(file = "data_files/ModulePreservation.RData", mps)
load("data_files/ModulePreservation.RData")

# extracting preservation Zscores from mp object
plotPreservationZscores <- function(mp, comparisonNames) {
  mp1 <- mp$preservation$Z[[1]][[2]]
  mp2 <- mp$preservation$Z[[1]][[3]]
  mp3 <- mp$preservation$Z[[1]][[4]]
  pres_df <- tibble(module = c(rownames(mp1), rownames(mp2), rownames(mp3)),
                    Zscore = c(mp1$Zsummary.pres, mp2$Zsummary.pres, mp3$Zsummary.pres),
                    modSize = c(mp1$moduleSize, mp2$moduleSize, mp3$moduleSize),
                    testGroup = c(rep(comparisonNames[1], nrow(mp1)), rep(comparisonNames[2], nrow(mp2)), rep(comparisonNames[3], nrow(mp3))))
  p <- ggplot(data = pres_df, aes(x = module, y = Zscore)) + 
    geom_point(aes(color = testGroup, size = modSize), alpha = 0.5) + 
    geom_hline(aes(yintercept = 10, color = "significance threshold")) + 
    scale_color_manual(values = c("cer" = "red", "par" = "blue", "hyc" = "pink", "hyp" = "skyblue"))
  output <- list(df = pres_df, plot = p)
  return(output)
}

preZ <- list(cer = NULL,
             par = NULL,
             hyc = NULL,
             hyp = NULL)
for (i in 1:4) {
  nms <- setdiff(c("cer", "par", "hyc", "hyp"), names(mps)[i])
  preZ[[i]] <- plotPreservationZscores(mps[[i]], nms)
}

# visualizing module preservation for each "species" in the other three "species"
# cer
preZ$cer$plot # there's perfect overlap between hyp and hyc on the grey module
# par
preZ$par$plot # weird that this isn't reciprocal --- cerevisiae's blue module has poor preservation in paradoxus, but paradoxus' blue module has good preservation in cerevisiae. I'll need to think more about why that is the case. Ah it's b/c the blue module is smaller in par, my plot doesn't do a great job of showing that
# hyc
preZ$hyc$plot
# hyp
preZ$hyp$plot

# does the species with larger module size always have worse preservation zscore (Zscore when it's the reference)?
# first have to update module sizes with their actual size---I think modulePreservation's module size is either downsampled or upsampled to be equal between comparisons, but the doccumentation is slim
for (i in 1:4) {
  colorIdxs <- match(preZ[[i]]$df$module, names(table(colors[[i]])))
  preZ[[i]]$df$modSize <- table(colors[[i]])[colorIdxs] %>% as.integer()
  preZ[[i]]$df[preZ[[i]]$df$module == "gold",]$modSize <- 100
}

preZ$cer$df$refGroup <- rep("cer", nrow(preZ$cer$df))
preZ$par$df$refGroup <- rep("par", nrow(preZ$par$df))
preZ$hyc$df$refGroup <- rep("hyc", nrow(preZ$hyc$df))
preZ$hyp$df$refGroup <- rep("hyp", nrow(preZ$hyp$df))
preZdf <- rbind(preZ$cer$df, preZ$par$df, preZ$hyc$df, preZ$hyp$df)

sizeVsZ <- tibble_row()
for (mod in reduce(list(names(table(colors$cer)), 
                        names(table(colors$par)), 
                        names(table(colors$hyc)), 
                        names(table(colors$hyp))), intersect)) {
  modDf <- preZdf[preZdf$module == mod,]
  for (sp1 in c("cer", "par", "hyc", "hyp")) {
    for (sp2 in setdiff(c("cer", "par", "hyc", "hyp"), sp1)) {
      sizeDiff <- filter(modDf, refGroup == sp1)$modSize[1] - filter(modDf, refGroup == sp2)$modSize[1]
      ZDiff <- filter(modDf, refGroup == sp1 & testGroup == sp2)$Zscore - filter(modDf, refGroup == sp2 & testGroup == sp1)$Zscore
      newRow <- tibble("module" = mod,
                       "refGroup" = sp1,
                       "testGroup" = sp2,
                       "sizeDiff" = sizeDiff,
                       "ZDiff" = ZDiff)
      sizeVsZ <- rbind(sizeVsZ, newRow)
    }
  }
}
# plotting sizeVsZ to see if larger mod always has smaller Z
ggplot(data = sizeVsZ, aes(x = sizeDiff, y = ZDiff)) + 
  geom_point(aes(color = module)) +
  xlab("Module size difference (reference - test)") + 
  ylab("Z score difference (reference - test)") + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "salmon", alpha = 0.2) + # these salmon regions are modules where the species with the larger module also has the higher Zscore metric for preservation, exceptions to my expectation
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "salmon", alpha = 0.2) +
  labs(color="module") +
  ggtitle("Module size vs preservation in each\n reference-test pair")
# it's elegant b/c it's a mirror image! we're counting every pair twice, once a-b and once b-a, so of course every point has a corresponding point flipped over x and y axis, but what's interesting is that close to 0 size difference, most of the points are reverse-sign---the smaller module has the higher preservation
# let's see this even better by adding lines between pairs
sizeVsZ$flippedSign <- sign(sizeVsZ$sizeDiff) == -sign(sizeVsZ$ZDiff)
sizeVsZ$orderDoesntMatterPair <- map2(sizeVsZ$refGroup, sizeVsZ$testGroup, function(sp1, sp2) {
  pair <- c(sp1, sp2)
  sps <- c("cer", "par", "hyc", "hyp")
  for (sp1 in sps) {
    if (sp1 %in% pair) {
      for (sp2 in setdiff(sps, sp1)) {
        if (sp2 %in% pair) {
          return(paste(sp1, sp2, sep="-"))
        }
      }
    }
  }
}) %>% unlist()
ggplot(data = sizeVsZ, aes(x = sizeDiff, y = ZDiff)) + 
  geom_point(color = sizeVsZ$module, size = 4) +
  xlab("Module size difference (reference - test)") + 
  ylab("Z score difference (reference - test)") + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "salmon", alpha = 0.2) + # these salmon regions are modules where the species with the larger module also has the higher Zscore metric for preservation, exceptions to my expectation
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "salmon", alpha = 0.2) +
  labs(color="Larger module has\nweaker preservation") +
  ggtitle("Module size vs preservation in each\n reference-test pair") +
  geom_line(aes(group = paste(module, orderDoesntMatterPair, sep="-"), color=flippedSign))
# not as pretty as the first plot, but shows things more clearly:
# the modules with less size difference are less likely for the size difference to influence preservation

# This establishes that I should use lowest Zscore not lowest size to gauge preservation of each module between 2 groups---although size is more intuitive, it sometimes results in picking the higher Zscore, and I think minimum Zscore is more relevant to determining if a module is truly preserved (but remember this was a decision)
# I don't think I even need a new function for this, just use pmin with the preZ per module
SharedModules <- reduce(list(preZ$cer$df$module, preZ$par$df$module, preZ$hyc$df$module, preZ$hyp$df$module), intersect)
preZ_parents <- tibble(module = SharedModules, 
                       Zmin = pmin(preZ$cer$df[match(SharedModules, preZ$cer$df$module),]$Zscore,
                                   preZ$par$df[match(SharedModules, preZ$par$df$module),]$Zscore))
preZ_hybrids <- tibble(module = SharedModules, 
                       Zmin = pmin(preZ$hyc$df[match(SharedModules, preZ$hyc$df$module),]$Zscore,
                                   preZ$hyp$df[match(SharedModules, preZ$hyp$df$module),]$Zscore))

preZ_parents$parentOrHybrid <- rep("parent", nrow(preZ_parents))
preZ_hybrids$parentOrHybrid <- rep("hybrid", nrow(preZ_hybrids))
preZ_comp <- rbind(preZ_parents, preZ_hybrids)
ggplot(data = preZ_comp, aes(x = module, y = Zmin)) + 
  geom_point(color = preZ_comp$module, size = 5, aes(shape = parentOrHybrid)) + 
  geom_hline(yintercept = 10, color = "yellow") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.70)) +
  guides(shape=guide_legend(title=element_blank())) + 
  ylab("Minimum Zscore for pair") +
  scale_x_discrete(labels = c("blue" = "blue", "brown" = "brown", "turquoise" = "turquoise", "purple" = "purple", "red" = "red", "gold" = "gold\n(random module)", "grey" = "grey\n(genes w/o module)"),
                   limits = c("blue", "brown", "turquoise", "purple", "red", "grey", "gold")) # this sets the order so that grey and gold are at the end
# okay ^this is fairly interesting, most modules have the same level of preservation in hybrids, but blue and brown really don't. Wonder why!
# BUT: a major issue here is that I'm throwing out all the modules that don't have the same color in all 4 "species"---there's a lot of genes not accounted for here

# TODO: there is some major size-bias afoot in modulePreservation that I didn't see with my non-split data (in fact I saw the opposite!)
# I'm especially concerned about the red module, as visually that module shares significant genes between samples
# So for now I'm leaving modulePreservation results out of my data presenting
