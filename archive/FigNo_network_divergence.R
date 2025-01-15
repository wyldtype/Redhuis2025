sapply(c("tidyverse"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Barkai_data_analysis/")

# Script to construct Figure 2: Networks in both species, bipartite graph, and divergent gene depletion in conserved blocks (probably)

# load single gene model dataset
load(file = "data_files/single_gene_models.RData")

# re-call sig if you want to change thresholds
p_thresh <- 0.05/(length(unique(spaldf$gene_name))*5)
eff_thresh <- 2
spaldf$sig <- spaldf$pvalue < p_thresh & abs(spaldf$effect_size) > eff_thresh & spaldf$isexpressed
sum(spaldf$sig)
spaldf %>% group_by(gene_name) %>% summarise(sig_gene = any(sig)) %>% select(sig_gene) %>% pull() %>% sum()
# re-creating genedf, where each row is one gene
genedf <-  spaldf %>% 
  pivot_wider(names_from = c("experiment", "coefficient"), values_from = c("effect_size", "pvalue", "sig"), id_cols = c("gene_name"))
table(genedf$gene_name) %>% table() # checking that every gene is only represented once

# load networks
load(file = "data_files/Networks_mergemodules05.RData")

#### Fig 2a: Bipartite graphs and divergent gene enrichment ####
library(bipartite)

makeBipartiteMatrix <- function(.top_colors, .bottom_colors) {
  module_names_top <- setdiff(unique(.top_colors), "grey")
  module_names_bottom <- setdiff(unique(.bottom_colors), "grey")
  color_counts_matrix <- matrix(nrow = length(module_names_bottom), ncol = length(module_names_top), dimnames = list(module_names_bottom, module_names_top))
  for (col_top in module_names_top) {
    for (col_bottom in module_names_bottom) {
      color_counts_matrix[col_bottom, col_top] <- sum((.bottom_colors == col_bottom) & (.top_colors == col_top))
    }
  }
  return(list(matrix = color_counts_matrix, colors_high = module_names_top, colors_low = module_names_bottom))
}
# parents
parents_color_counts <- makeBipartiteMatrix(colors$cer, colors$par)
plotweb(parents_color_counts$matrix, col.high = parents_color_counts$colors_high, col.low = parents_color_counts$colors_low, labsize = 1) # oh my goodness would but everything be this easy <3 <3 <3

# Are divergent genes less common in conserved blocks?
block_table_cer_par <- table(paste(colors$cer, colors$par)) %>% sort(decreasing = TRUE)
block_table_cer_par
genedf$module_color_cer <- sapply(genedf$gene_name, \(g) {
  g_idx <- which(colnames(counts_top3000$cer) == g)
  color <- NA
  if (length(g_idx) == 1) {
    color <- colors$cer[g_idx] %>% unlist() %>% as.character()
  }
  return(color)
}) %>% as.character()
genedf$module_color_par <- sapply(genedf$gene_name, \(g) {
  g_idx <- which(colnames(counts_top3000$par) == g)
  color <- NA
  if (length(g_idx) == 1) {
    color <- colors$par[g_idx] %>% unlist() %>% as.character()
  }
  return(color)
}) %>% as.character()
table(paste(genedf$module_color_cer, genedf$module_color_par)) %>% sort(decreasing = TRUE) # ~ 2000 genes were thrown out as not in the top 3000 most connected genes (but we could change that!)
plotdf <- genedf
plotdf$block_size <- map2(plotdf$module_color_cer, plotdf$module_color_par, \(x, y) {
  block_size <- NA
  block_name <- paste(x, y)
  if (!is.na(x) & !is.na(y)) {
    block_size <- block_table_cer_par[which(names(block_table_cer_par) == block_name)] %>% as.numeric()
  }
  return(block_size)
}) %>% unlist()
plotdf$max_effect_mag_species <- pmax(abs(plotdf$effect_size_CC_species), 
                                       abs(plotdf$effect_size_HAP4_species),
                                       abs(plotdf$effect_size_LowPi_species),
                                       abs(plotdf$effect_size_TFdelxLowN_species))

plotdf$sig_any_species <- apply(plotdf, 1, \(x) {
  return(any(x["sig_CC_species"],
             x["sig_HAP4_species"],
             x["sig_TFdelxLowN_species"],
             x["sig_LowPi_species"]))
})
ggplot(plotdf, aes(x = block_size, y = max_effect_mag_species)) + geom_point(aes(color = sig_any_species)) + geom_vline(xintercept = 55, color = "gold") + theme_classic()
# looks like a block size of 55 is a good cutoff. Just visually it seems like there's something different happening in blocks that are smaller than that vs blocks that are larger
plotdf$is_large_block <- if_else(plotdf$block_size > 55, true = "large", false = "small")
plotdf %>% filter(sig_any_species) %>% ggplot(aes(x = is_large_block, y = max_effect_mag_species)) + geom_jitter()
plotdf %>% ggplot(aes(x = is_large_block, y = max_effect_mag_species)) + geom_jitter(aes(color = sig_any_species))
# effect size isn't necessarily larger, but are there more sig genes in smaller blocks?
table(plotdf$is_large_block, plotdf$sig_all_species)
fisher.test(table(plotdf$is_large_block, plotdf$sig_all_species))
# I'm fairly well convinced that it isn't this simple. Larger blocks of co-regulated genes aren't necessarily less divergent at the single-gene level
# But there are definitely a few heavy hitters in the small blocks. Maybe there's something going on here, but it's not as simple as all the divergent genes being in smaller blocks

# Proportion of sig genes vs block size per block
# TODO: make this work. plotdf is based on genedf right now, so sig isn't a column and my brain hurts trying to figure out what to do
plotdf2 <- plotdf %>% group_by(block_size) %>% 
  summarise(prop_sig = sum(sig)/length(sig))
ggplot(plotdf2, aes(x = block_size, y = prop_sig)) + geom_point()

#### Gene Ontology Enrichment ####

# What are the functions of the conserved blocks of co-expressed genes?
genedf$block_name <- paste(genedf$module_color_cer, genedf$module_color_par, sep = "_")
for (block in unique(genedf$block_name)) {
  gene_names <- genedf %>% filter(block_name == block) %>% select(gene_name) %>% pull()
  if (length(gene_names) > 62) {
    write.table(gene_names, file = paste0("gene_ontology/gene_lists/conserved/", block, ".txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}
# how about the hybrid where there's that giant module? Is that all the rRNA machinery together at last?
for (hyb_col in unique(colors$hyb)) {
  gene_names <- colnames(counts_all5$hyb)[colors$hyb == hyb_col]
  write.table(gene_names, file = paste0("gene_ontology/gene_lists/hybrid/", hyb_col, ".txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

GOProcessSummary <- list(cer = tibble("pink" = "",
                                      "blue" = "",
                                      "blue" = "",
                                      "blue" = "",
                                      "blue" = "",
                                      "brown" = "",
                                      "turquoise" = "",
                                      "turquoise" = "",
                                      "yellow" = "",
                                      "yellow" = "",
                                      "yellow" = "",
                                      "yellow" = ""), # best way I can think to have both colors in a header in huxtable is to put cer's blank then par's with the block enrichments
                         par = tibble("pink" = "oxidative phosphorylation",
                                      "blue" = "membrane proteins\n(no significant process enrichment)",
                                      "green" = "G-protein signalling",
                                      "lightcyan" = "golgi/vacuole transport",
                                      "turquoise" = "ribosome biogenesis",
                                      "green" = "filamentous growth/allantonin metabolism",
                                      "lightcyan" = "translation",
                                      "turquoise" = "ribosome biogenesis",
                                      "blue" = "vesicle-mediated golgi transport",
                                      "turquoise" = "endoplasmic reticulum proteins\n(no significant process enrichment)",
                                      "lightcyan" = "cytoskeleton",
                                      "red" = "sulfur/toxin metabolism"))

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

# how about genes in each module in each species that aren't part of conserved blocks?
genedf$block_size <- sapply(genedf$block_name, \(b) {
  return(sum(genedf$block_name == b))
})
conserved_blocks <- unique(genedf$block_name[genedf$block_size > 62])
# cer
for (cer_mod in unique(genedf$module_color_cer)) {
  regexblock <- paste0("^", cer_mod, "_")
  forbidden_par <- conserved_blocks[grep(regexblock, conserved_blocks)]
  forbidden_par <- gsub(regexblock, "", forbidden_par)
  gene_names <- genedf %>% filter(module_color_cer == cer_mod & !(module_color_par %in% c(cer_mod, forbidden_par))) %>% select(gene_name) %>% pull()
  if (length(gene_names) > 10) {
    write.table(gene_names, file = paste0("gene_ontology/gene_lists/divergent/", cer_mod, "_cer.txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

GOdf_cer <- tibble("blue" = "beta-D-glucan metabolism\n(cell wall component)",
                   "green" = "nitrogen utilization",
                   "tan" = "response to chemical/oxidative stress")

# par
for (par_mod in unique(genedf$module_color_par)) {
  regexblock <- paste0("_", par_mod, "$")
  forbidden_cer <- conserved_blocks[grep(regexblock, conserved_blocks)]
  forbidden_cer <- gsub(regexblock, "", forbidden_cer)
  gene_names <- genedf %>% filter(!(module_color_cer %in% c(forbidden_cer, par_mod)) & module_color_par == par_mod) %>% select(gene_name) %>% pull()
  if (length(gene_names) > 10) {
    write.table(gene_names, file = paste0("gene_ontology/gene_lists/divergent/", par_mod, "_par.txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

GOdf_par <- tibble("green" = "amino acid transport",
                   "greenyellow" = "polyphosphate metabolism",
                   "pink" = "mitochondrial respiratory chain",
                   "tan" = "glycogen metabolism")


#### Module expression patterns ####
getGeneDf <- function(.gene_name, .mode = c("parents", "hybrid"), .experiment = c("TFdelxLowN", "CC", "HAP4", "LowPi")) {
  if (.mode == "parents") {
    info <- spinfo[.experiment] %>% Reduce(f = bind_rows)
    cts <- spcts[.experiment] %>% Reduce(f = cbind)
  }
  if (.mode == "hybrid") {
    info <- alinfo[.experiment] %>% Reduce(f = bind_rows)
    cts <- alcts[.experiment] %>% Reduce(f = cbind)
  }
  genedf <- bind_cols(tibble(expr = cts[.gene_name,]), info) %>% 
    mutate(non_unique_sample_name = paste(condition, well_flask_ID)) %>% 
    select(expr, non_unique_sample_name, allele, experiment, time_point_num, genotype) %>% 
    pivot_wider(id_cols = c(non_unique_sample_name, experiment, time_point_num, genotype), values_from = expr, names_from = allele)
  return(genedf)
}
visualizeGenesAcrossEnvironments <- function(.gene_names) { # just cer/par for now, but can add more if it seems useful (cer/par networks seem like the only interesting ones)
  cts_cer <- counts_top3000[["cer"]]
  info_cer <- spinfo_full %>% filter(organism == "cer") %>% select(time_point_num, genotype, experiment, sample_name)
  plotdf_cer <- lapply(.gene_names, \(g) {
    gdf <- tibble(gene_name = g,
                  expr = cts_cer[,g],
                  sample_name = rownames(cts_cer),
                  organism = "cer")
    output <- inner_join(gdf, info_cer, by = "sample_name")
    return(output)
  }) %>% Reduce(f = bind_rows)
  cts_par <- counts_top3000[["par"]]
  info_par <- spinfo_full %>% filter(organism == "par") %>% select(time_point_num, genotype, experiment, sample_name)
  plotdf_par <- lapply(.gene_names, \(g) {
    gdf <- tibble(gene_name = g,
                  expr = cts_par[,g],
                  sample_name = rownames(cts_par),
                  organism = "par")
    output <- inner_join(gdf, info_par, by = "sample_name")
    return(output)
  }) %>% Reduce(f = bind_rows)
  plotdf <- bind_rows(plotdf_cer, plotdf_par)
  plotdf$environment <- order(plotdf$experiment,
                              plotdf$genotype,
                              plotdf$time_point_num) %>% order() # this as it turns out is how we have to do "rank" using multiple columns
  
  
  p <- plotdf %>% group_by(experiment, genotype, time_point_num) %>% slice_sample(n = 2500) %>% 
    ggplot(aes(x = environment, y = log(expr + 1), group = organism)) + 
    # geom_point(aes(color = paste(organism, experiment))) + 
    geom_smooth(method = "loess", aes(color = organism)) +
    xlab("Environment") + scale_x_discrete(labels = element_blank()) +
    ylab("Expression Level (log)") + theme_classic()
  return(p)
}
# sulfur/toxin module
gene_names <- colnames(counts_top3000[["cer"]])[colors[["cer"]] == "yellow" & colors[["par"]] == "red"]
sulftox <- visualizeGenesAcrossEnvironments(.gene_names = gene_names) + ggtitle("sulfur/toxin")

# cytoskeleton module
gene_names <- colnames(counts_top3000[["cer"]])[colors[["cer"]] == "yellow" & colors[["par"]] == "lightcyan"]
cyto <- visualizeGenesAcrossEnvironments(.gene_names = gene_names) + ggtitle("cytoskeleton")

# small ribosome subunit biogenesis module
gene_names <- colnames(counts_top3000[["cer"]])[colors[["cer"]] == "blue" & colors[["par"]] == "turquoise"]
visualizeGenesAcrossEnvironments(.gene_names = gene_names) + ggtitle("ribosome biogenesis block 1")

# oxphos
gene_names <- colnames(counts_top3000[["cer"]])[colors[["cer"]] == "pink" & colors[["par"]] == "pink"]
oxphos <- visualizeGenesAcrossEnvironments(.gene_names = gene_names) + ggtitle("oxidative phosphorylation")

# betaglucan
gene_names <- colnames(counts_top3000[["cer"]])[colors[["cer"]] == "blue" & !(colors[["par"]] %in% c("blue", "green", "turquoise", "lightcyan"))]
betaglu <- visualizeGenesAcrossEnvironments(.gene_names = gene_names) + ggtitle("beta-glucan")

ggarrange(oxphos, betaglu, cyto, sulftox, common.legend = TRUE)
