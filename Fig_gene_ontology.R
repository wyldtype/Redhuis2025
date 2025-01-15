sapply(c("dplyr", "purrr", "ggplot2", "ggpubr", "tidyr", "readr", "data.table"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")

# load functions and datasets
source(file = "functions_for_figure_scripts.R")
source(file = "data_for_figure_scripts.R")

# TODO: check the DJ-1 superfamily as those up-par down-cer 500 genes in Sat Growth


#### printing ccm gene lists for SGD ####
printGeneList <- function(.gene_idxs, .file = "gene_ontology/test.txt") {
  return(write.table(.gene_idxs, 
              file = .file,
              quote = FALSE, 
              row.names = FALSE, 
              col.names = FALSE))
}
# tests for printGeneList
test <- printGeneList(.gene_idxs = c("YGR192C", "YJR009C", "YJL052W"))


#### in-house GO enrichment test based on GOslim table ####
# load GO slim table
goslim <- read.table("data_files/go_slim_mapping.tab", header = FALSE, sep = "\t") |>
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
                   (ORF %in% genedf$gene_name))

# do this per group of interest
getGOSlimDf <- function(.idxs, .group_name, .file_prefix = "gene_ontology/results/",
                        .min_hits = 5) {
  test_table <- goslim |> filter(ORF %in% .idxs) |> select(GOslim_term) |> table()
  test_table <- test_table[test_table >= .min_hits]
  testdf <- tibble("term" = names(test_table),
                   "group_count" = as.numeric(test_table))
  testdf$overall_count <- map(testdf$term, \(x) {
    ct <- goslim |> filter(GOslim_term == x) |> select(ORF) |> 
      pull() |> unique() |> length()
    return(ct)
  }) |> unlist()
  testdf$exact_pval <- map2(testdf$group_count, testdf$overall_count, \(x, y) {
    n_mod <- length(.idxs)
    n_genes <- nrow(genedf)
    exact_table <- rbind(c(x, n_mod - x), c(y - x, n_genes - n_mod - y))
    exact_result <- fisher.test(exact_table, alternative = "greater")
    return(exact_result$p.value)
  }) |> unlist()
  testdf$sig <- testdf$exact_pval < 0.001
  testdf$genes <- map(testdf$term, \(x) {
    termgenes <- goslim |> filter(ORF %in% .idxs & GOslim_term == x) |> 
      select(ORF) |> pull() |> unique()
    return(reduce(termgenes, paste, sep = " "))
  }) |> unlist()
  
  write_csv(arrange(testdf, exact_pval, desc(group_count)), 
            file = paste0(.file_prefix, .group_name, ".csv"),
            quote = "none", 
            col_names = TRUE)
  return(testdf)
}

# tests for getGOSlimDf
test <- getGOSlimDf(.idxs = c("YGR192C", "YJR009C", "YJL052W"),
                    .group_name = "tdhs", .min_hits = 1)

#### Dynamics ####

# Saturated Growth 2-1
gene_idxs <- finaldf |> filter(experiment == "HAP4" & group == "dyn21") |> 
  select(gene_name) |> pull()
HAP4_21 <- getGOSlimDf(.idxs = gene_idxs,
                       .group_name = "HAP4_dyn21")
# vesicles, COPI/COPII, golgi ER transport

# LowPi 1-2
gene_idxs <- finaldf |> filter(experiment == "LowPi" & group == "dyn12") |> 
  select(gene_name) |> pull()
LowPi_12 <- getGOSlimDf(.idxs = gene_idxs,
                       .group_name = "LowPi_12")
# stress

# Heat 2-1
gene_idxs <- finaldf |> filter(experiment == "Heat" & group == "dyn21") |> 
  select(gene_name) |> pull()
Heat_21 <- getGOSlimDf(.idxs = gene_idxs,
                        .group_name = "Heat_21")
# tRNA/nucleobase catabolic process/regulation of RNA Pol II

# LowN 2-1
gene_idxs <- finaldf |> filter(experiment == "LowN" & group == "dyn21") |> 
  select(gene_name) |> pull()
LowN_21 <- getGOSlimDf(.idxs = gene_idxs,
                       .group_name = "LowN_21")

#### Level ####
# Sat Growth upcer 1
gene_idxs <- finaldf |> filter(experiment == "HAP4" & group == "levupcer1") |> 
  select(gene_name) |> pull()
HAP4_upcer_1 <- getGOSlimDf(.idxs = gene_idxs,
                            .group_name = "HAP4_upcer_1")
# plasma membrane/cell wall

# Sat Growth uppar 1
gene_idxs <- finaldf |> filter(experiment == "HAP4" & group == "levuppar1") |> 
  select(gene_name) |> pull()
HAP4_uppar_1 <- getGOSlimDf(.idxs = gene_idxs,
                            .group_name = "HAP4_uppar_1")
# ribosome

# Heat upcer 1
gene_idxs <- finaldf |> filter(experiment == "Heat" & group == "levupcer1") |> 
  select(gene_name) |> pull()
Heat_upcer_1 <- getGOSlimDf(.idxs = gene_idxs,
                            .group_name = "Heat_upcer_1")
# caaaaatabolic process





#### GO enrichment across cluster combinations ####
# Rationale: thus far we've looked at clusters ID'd in a single
# environment comprising like 100s of genes. Let's try to look
# more fine-scale at the genes that are all in the same
# combination of groups across environments

#### Archive ####
# generating the GO results for each ccm to see which terms are
# the most unique to each module

# merge 00
godf00 <- select(moduledf00, is_CCM, module_name, CCM_color, block_size)
godf00$nsigGOterms <- 0
for (m in godf00$module_name) {
  gene_idxs <- module_genedf00 |> filter(module_name == m) |> 
    select(gene_name) |> pull()
  printGeneList(.gene_idxs = gene_idxs,
                .file = paste0("gene_ontology/merge00/", m, ".txt"))
  testdf <- getGOSlimDf(m, 
                        .file_prefix = "gene_ontology/merge00/results/",
                        .idxs = gene_idxs)
  godf00$nsigGOterms[godf00$module_name == m] <- sum(testdf$sig)
}

# random 00
random_godf00 <- select(random_moduledf00, is_CCM, module_name, CCM_color, block_size)
random_godf00$nsigGOterms <- 0
for (m in random_moduledf00$module_name) {
  mod_idxs <- random_module_genedf00 |> filter(module_name == m) |> select(gene_name) |> pull()
  printGeneList(.gene_idxs = mod_idxs,
                .file = paste0("gene_ontology/random00/", m, ".txt"))
  testdf <- getGOSlimDf(.mod = m, 
                        .file_prefix = "gene_ontology/random00/results/", 
                        .idxs = mod_idxs)
  random_godf00$nsigGOterms[random_godf00$module_name == m] <- sum(testdf$sig)
}

# merge 10
godf10 <- select(moduledf10, is_CCM, module_name, CCM_color, block_size)
godf10$nsigGOterms <- 0
for (m in godf10$module_name) {
  gene_idxs <- module_genedf10 |> filter(module_name == m) |> 
    select(gene_name) |> pull()
  printGeneList(.gene_idxs = gene_idxs,
                .file = paste0("gene_ontology/merge10/", m, ".txt"))
  testdf <- getGOSlimDf(m, 
                        .file_prefix = "gene_ontology/merge10/results/",
                        .idxs = gene_idxs)
  godf10$nsigGOterms[godf10$module_name == m] <- sum(testdf$sig)
}

# random 10
random_godf10 <- select(random_moduledf10, is_CCM, module_name, CCM_color, block_size)
random_godf10$nsigGOterms <- 0
for (m in random_moduledf10$module_name) {
  mod_idxs <- random_module_genedf10 |> filter(module_name == m) |> select(gene_name) |> pull()
  printGeneList(.gene_idxs = mod_idxs,
                .file = paste0("gene_ontology/random10/", m, ".txt"))
  testdf <- getGOSlimDf(.mod = m, 
                        .file_prefix = "gene_ontology/random10/results/", 
                        .idxs = mod_idxs)
  random_godf10$nsigGOterms[random_godf10$module_name == m] <- sum(testdf$sig)
}

# merge 25
godf25 <- select(moduledf25, is_CCM, module_name, CCM_color, block_size)
godf25$nsigGOterms <- 0
for (m in godf25$module_name) {
  gene_idxs <- module_genedf25 |> filter(module_name == m) |> 
    select(gene_name) |> pull()
  printGeneList(.gene_idxs = gene_idxs,
                .file = paste0("gene_ontology/merge25/", m, ".txt"))
  testdf <- getGOSlimDf(m, 
                        .file_prefix = "gene_ontology/merge25/results/",
                        .idxs = gene_idxs)
  godf25$nsigGOterms[godf25$module_name == m] <- sum(testdf$sig)
}

# random 25
random_godf25 <- select(random_moduledf25, is_CCM, module_name, CCM_color, block_size)
random_godf25$nsigGOterms <- 0
for (m in random_moduledf25$module_name) {
  mod_idxs <- random_module_genedf25 |> filter(module_name == m) |> select(gene_name) |> pull()
  printGeneList(.gene_idxs = mod_idxs,
                .file = paste0("gene_ontology/random25/", m, ".txt"))
  testdf <- getGOSlimDf(.mod = m, 
                        .file_prefix = "gene_ontology/random25/results/", 
                        .idxs = mod_idxs)
  random_godf25$nsigGOterms[random_godf25$module_name == m] <- sum(testdf$sig)
}

# merge 35
godf35 <- select(moduledf35, is_CCM, module_name, CCM_color, block_size)
godf35$nsigGOterms <- 0
for (m in godf35$module_name) {
  gene_idxs <- module_genedf35 |> filter(module_name == m) |> 
    select(gene_name) |> pull()
  printGeneList(.gene_idxs = gene_idxs,
                .file = paste0("gene_ontology/merge35/", m, ".txt"))
  testdf <- getGOSlimDf(m, 
                        .file_prefix = "gene_ontology/merge35/results/",
                        .idxs = gene_idxs)
  godf35$nsigGOterms[godf35$module_name == m] <- sum(testdf$sig)
}

# random 35
random_godf35 <- select(random_moduledf35, is_CCM, module_name, CCM_color, block_size)
random_godf35$nsigGOterms <- 0
for (m in random_moduledf35$module_name) {
  mod_idxs <- random_module_genedf35 |> filter(module_name == m) |> select(gene_name) |> pull()
  printGeneList(.gene_idxs = mod_idxs,
                .file = paste0("gene_ontology/random35/", m, ".txt"))
  testdf <- getGOSlimDf(.mod = m, 
                        .file_prefix = "gene_ontology/random35/results/", 
                        .idxs = mod_idxs)
  random_godf35$nsigGOterms[random_godf35$module_name == m] <- sum(testdf$sig)
}

#### plotting nsig GO terms by module type ####

# merge 25
godf25$type <- if_else(godf25$is_CCM, true = "conserved", false = "non-conserved")
random_godf25$type <- "random"
# plotdf <- bind_rows(godf25, random_godf25) |> filter(block_size >= 30)
p25 <- ggplot(godf25, aes(x = block_size, y = nsigGOterms)) + 
  geom_point(aes(color = type)) +
  scale_color_discrete(limits = c("conserved", "non-conserved"),
                       type = c("gold", "grey")) +
  ylab("number of significant GO terms") +
  xlab("module size (number of genes)") +
  ggtitle("GO term enrichment for\n conserved vs non-conserved \nmodules") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "left")
library(ggExtra)
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/GOterms25.pdf",
    width = 5, height = 3)
ggMarginal(p25, groupColour = TRUE, groupFill = TRUE, margins = "y")
dev.off()
bind_rows(godf25, random_godf25) |> 
  filter(block_size >= 30) |> 
  group_by(type) |> 
  summarise(avg_nsig = mean(nsigGOterms),
            med_nsig = median(nsigGOterms),
            min_nsig = min(nsigGOterms),
            max_nsig = max(nsigGOterms),
            n = n(),
            n_0 = sum(nsigGOterms == 0))

# merge 00
godf00$type <- if_else(godf00$is_CCM, true = "conserved", false = "non-conserved")
random_godf00$type <- "random"
bind_rows(godf00, random_godf00) |> 
  filter(block_size >= 30) |> 
  group_by(type) |> 
  summarise(avg_nsig = mean(nsigGOterms),
            med_nsig = median(nsigGOterms),
            min_nsig = min(nsigGOterms),
            max_nsig = max(nsigGOterms))

# merge 10
godf10$type <- if_else(godf10$is_CCM, true = "conserved", false = "non-conserved")
random_godf10$type <- "random"
bind_rows(godf10, random_godf10) |> 
  filter(block_size >= 30) |> 
  group_by(type) |> 
  summarise(avg_nsig = mean(nsigGOterms),
            med_nsig = median(nsigGOterms),
            min_nsig = min(nsigGOterms),
            max_nsig = max(nsigGOterms))

# merge 35
godf35$type <- if_else(godf35$is_CCM, true = "conserved", false = "non-conserved")
random_godf35$type <- "random"
bind_rows(godf35, random_godf35) |> 
  filter(block_size >= 30) |> 
  group_by(type) |> 
  summarise(avg_nsig = mean(nsigGOterms),
            med_nsig = median(nsigGOterms),
            min_nsig = min(nsigGOterms),
            max_nsig = max(nsigGOterms))

#### summary of most unique aspects of each ccm ####
# a now g - RNA processing/splicing/metabolism
# b - phosphate metabolism
# c now e - DNA repair/metabolism/binding, cytoskeleton, some RNA processing, endosome
# d now f - same as c plus histone remodeling (SWI/SNF and RPD3l complexes)
#     and ubiquitin (SCF and ligase) complexes
# e now c- small/large ribosomal subunits, translation, rRNA processing
# f now b - small/large ribosomal subunits, transcription,
#     RNA methylation (box C/D RNP complex), protein folding (including chaperone complex)
#     translation (including eukaryotic translation initiation factor 3 complex)
# g now a - nitrogen and sulfur metabolism
# h - transmembrane transport, transferase, lyase, oxidoreductase, mitochondrion
# i - mitochondrial large/small subunit, mitochondrial translation
# j - ATP synthase complex coupling factor F(o), respiratory chain complexes, MICOS complex
# k now d - mitochondrial respiratory chain enzymes (Acetyl coA synthase, NADH dehydrogenase, aldehyde dehydrogenase, glutamate dehydrogenase
#     other mitochondrial proteins: (membrane protein, acetate transporter, inner membrane proteins)
# l now h - glucose metabolism/low glucose growth (5/79), mitochondrial organization (6/79), mitochondrial proteins (22/79)
# m - same as l, but adding glycosyl transferases/glycogen metabolism and removing membrane/mitochondrial organization including ARP2/3 complex

# a - sulfur metabolism
# b - small/large ribosomal subunits, RNA Pol III transcription,
#     RNA modification, protein folding chaperone complex,
#     translation, tRNA processing
# c - same as b without chaperone
# d - mitochondrial respiratory chain enzymes, other mitochondrial proteins
# e - DNA repair/metabolism/binding, cytoskeleton, some RNA processing, endosome
# f - same as e plus histone remodeling
#     and ubiquitin complexes
# g - RNA splicing, proteasome complex, protein folding
# h - glucose metabolism/low glucose growth (x/x), mitochondrial organization (x/x), mitochondrial proteins (x/x)
# i - small/large ribosomal subunits, transcription, RNA Pol I and III
# j - small/large ribosomal subunits, translation

#### Combining GO terms ####

# based on GOslim table results we summarise below
GOtable25 <- tibble("CCM_name" = moduledf25$module_name[moduledf25$is_CCM],
                    "block_size" = moduledf25$block_size[moduledf25$is_CCM],
                    "GOslim_combined_terms" = "")
# @input: GOslim enrichment results df for one ccm, character vector of GOslim terms to combine
# @output: character vector of all the genes associated with those terms in that ccm, each gene only appearing once
getCombinedTermGenes <- function(.df, .term_vec) {
  .df |> filter(term %in% .term_vec) |> 
    select(genes) |> pull() |> map(.f = strsplit, split = " ") |> 
    unlist() |> unique()
}

# manually combining terms
term_RNAsplicingAndMetabolism <- c("mRNA binding", 
                                   "mRNA processing",
                                   "nucleobase-containing compound catabolic process",
                                   "nucleobase-containing small molecule metabolic process",
                                   "RNA binding",
                                   "RNA catabolic process",
                                   "RNA splicing",
                                   "protein-RNA complex assembly",
                                   "small nuclear ribonucleoprotein complex",
                                   "spliceosomal complex",
                                   "U2-type spliceosomal complex",
                                   "U4/U6 x U5 tri-snRNP complex",
                                   "U6 snRNP",
                                   "box C/D RNP complex")

term_chaperone <- c("chaperonin-containing T-complex",
                    "prefoldin complex",
                    "protein folding",
                    "protein folding chaperone complex",
                    "protein glycosylation",
                    "protein maturation",
                    "protein modification by small protein conjugation or removal",
                    "protein targeting",
                    "unfolded protein binding")

term_proteasome <- c("proteasome core complex, alpha-subunit complex",
                     "proteasome core complex, beta-subunit complex",
                     "proteasome regulatory particle, lid subcomplex",
                     "proteasome storage granule")

term_mRNADegrading <- c("P-body", "Lsm1-7-Pat1 complex")

term_vacuoleGolgiTransport <- c("Golgi vesicle transport",
                                "transport",
                                "vacuolar transport",
                                "vacuole",
                                "vacuole organization",
                                "vesicle-mediated transport",
                                "regulation of transport",
                                "vesicle organization")

term_metabolism <- c("carbohydrate metabolic process",
                     "catabolic process",
                     "cellular nitrogen compound metabolic process",
                     "hydrolase activity",
                     "nucleobase-containing compound catabolic process",
                     "nucleobase-containing small molecule metabolic process",
                     "small molecule metabolic process")

term_cytoskeleton <- c("cytoplasmic microtubule",
                       "cytoskeletal protein binding",
                       "cytoskeleton organization",
                       "cytoskeleton-dependent intracellular transport",
                       "microtubule",
                       "microtubule organizing center",
                       "organelle assembly",
                       "organelle fusion",
                       "regulation of organelle organization",
                       "site of polarized growth",
                       "cell cortex")

term_cellularRespiration <- c("generation of precursor metabolites and energy",
                              "MICOS complex",
                              "mitochondrial proton-transporting ATP synthase complex  coupling factor F(o)",
                              "mitochondrial respiratory chain complex III",
                              "mitochondrial respiratory chain complex IV",
                              "cellular respiration")

term_DNArepairMetabolismBinding <- c("DNA binding",
                                     "DNA damage response",
                                     "DNA metabolic process",
                                     "DNA recombination",
                                     "DNA repair",
                                     "DNA-binding transcription factor activity",
                                     "DNA-templated transcription",
                                     "telomere organization")

term_transcription <- c("RNA polymerase II transcription regulator complex",
                        "DNA-templated transcription",
                        "transcription",
                        "DNA-templated transcription elongation",
                        "DNA-templated transcription termination",
                        "RNA polymerase I complex",
                        "RNA polymerase II core complex",
                        "RNA polymerase III complex",
                        "transcription by RNA polymerase I",
                        "transcription by RNA polymerase II",
                        "transcription by RNA polymerase III")

term_translation <- c("translation",
                      "regulation of translation",
                      "translational initiation",
                      "translational elongation",
                      "eukaryotic 48S preinitiation complex",
                      "eukaryotic translation initiation factor 3 complex",
                      "translation factor activity RNA binding")

term_tRNA <- c("tRNA aminoacylation for protein translation",
               "tRNA metabolic process",
               "tRNA methyltransferase complex",
               "tRNA processing")

term_TF <- c("DNA-binding transcription factor activity",
             "transcription factor binding",
             "transcription regulator complex")

term_chromatin <- c("chromatin",
                    "chromatin binding",
                    "chromatin organization",
                    "NuA4 histone acetyltransferase complex",
                    "histone acetyltransferase complex",
                    "histone deacetylase complex",
                    "rDNA heterochromatin")

term_cellCycle <- c("chromosome segregation",
                    "regulation of cell cycle")

term_ubiquitin <- c("ubiquitin ligase complex",
                    "ubiquitin-like protein binding")

term_exoendocytosis <- c("endocytosis",
                         "exocytosis",
                         "SNARE complex",
                         "endosomal transport",
                         "endosome")

term_endocytosis <- c("endocytosis",
                         "endosomal transport",
                         "endosome")

term_membrane <- c("membrane",
                   "membrane fusion",
                   "membrane organization",
                   "transmembrane transporter activity",
                   "transmembrane transport")

term_ribosome <- c("90S preribosome",
                   "cytosolic large ribosomal subunit",
                   "cytosolic small ribosomal subunit",
                   "cytosolic ribosome",
                   "large ribosomal subunit",
                   "preribosome",
                   "preribosome large subunit precursor",
                   "ribonucleoprotein complex",
                   "ribosome biogenesis",
                   "ribosomal subunit export from nucleus",
                   "rRNA processing",
                   "sno(s)RNA-containing ribonucleoprotein complex",
                   "rRNA binding",
                   "small ribosomal subunit",
                   "structural constituent of ribosome",
                   "preribosome small subunit precursor",
                   "Pwp2p-containing subcomplex of 90S preribosome",
                   "small-subunit processome",
                   "t-UTP complex")

term_mitochondrialTranslation <- c("mitochondrial large ribosomal subunit",
                                   "mitochondrial small ribosomal subunit",
                                   "mitochondrial translation")

term_mitochondria <- c("mitochondrion",
                       "mitochondrion organization")

#### hardcoding combined term fractions ####
# a, sulfer metabolism
testdf <- getGOSlimDf("a")
moduledf25 |> filter(module_name == "a") |> select(block_size)

genes_sulfur <- getCombinedTermGenes(testdf, c("sulfur compound metabolic process"))
genes_transferase <- getCombinedTermGenes(testdf, c("transferase activity",
                                                    "transferase activity, transferring alkyl or aryl (other than methyl) groups"))
genes_transmembrane <- getCombinedTermGenes(testdf, c("transmembrane transport",
                                                      "transmembrane transport activity"))

length(genes_transferase)
length(genes_sulfur)
length(genes_transmembrane)

sum(genes_transmembrane %in% genes_transferase)
sum(genes_sulfur %in% genes_transferase)
sum(genes_transmembrane %in% genes_sulfur)

# term heiarchy (most specific is leftmost):
# sulfur > transferase = transmembrane

# hardcoding a
GOtable25[GOtable25$CCM_name == "a", "GOslim_combined_terms"] <- 
  paste("sulfur metabolism (28/207)",
        "transferase (40/207)",
        "transmembrane transport (20/207)",
        sep = ",\n\n")

# b, small/large ribosomal subunits, transcription,
#     protein folding (including chaperone complex)
#     translation (including eukaryotic translation initiation factor 3 complex)
#     cytoplasmic stress granule
testdf <- getGOSlimDf("b")
moduledf25 |> filter(module_name == "b") |> select(block_size)

genes_ribosome <- getCombinedTermGenes(testdf, term_ribosome)
length(genes_ribosome) # 86/271

genes_translation <- getCombinedTermGenes(testdf, term_translation)
length(genes_translation) # 17/271

genes_chaperone <- getCombinedTermGenes(testdf, term_chaperone)
length(genes_chaperone) # 25/271

genes_tRNA <- getCombinedTermGenes(testdf, term_tRNA)
length(genes_tRNA) # 41/271

genes_stress <- getCombinedTermGenes(testdf, "cytoplasmic stress granule")
length(genes_stress) # 18/271

table(c(genes_ribosome,
        genes_translation,
        genes_tRNA,
        genes_chaperone,
        genes_stress)) |> table() # 133 nonoverlapping, 21 overlapping 2, 4 overlapping 3

overlapdf <- bind_rows(tibble("gene_name" = genes_ribosome,
                              "term_name" = "ribosome"),
                       tibble("gene_name" = genes_translation,
                              "term_name" = "translation"),
                       tibble("gene_name" = genes_tRNA,
                              "term_name" = "tRNA"),
                       tibble("gene_name" = genes_chaperone,
                              "term_name" = "chaperone"),
                       tibble("gene_name" = genes_stress,
                              "term_name" = "stress")) |> 
  pivot_wider(id_cols = "gene_name", names_from = "term_name", values_from = "term_name")
overlapdf$all_terms <- paste(overlapdf$ribosome,
                             overlapdf$translation,
                             overlapdf$tRNA,
                             overlapdf$chaperone,
                             overlapdf$stress, sep = "_")
table(overlapdf$all_terms)

sum(genes_stress %in% genes_chaperone)
# term heiarchy:
# stress = chaperone > tRNA > ribosome, translation not worth including

# hardcoding b
GOtable25[GOtable25$CCM_name == "b", "GOslim_combined_terms"] <- 
  paste("ribosome components (73/271) ",
        "tRNA processing (36/271)",
        "protein chaperone/modification (25/271)",
        "stress granule (18/271)",
        sep = ",\n\n")

# c, small/large ribosomal subunits, translation, rRNA processing
testdf <- getGOSlimDf("c")
moduledf25 |> filter(module_name == "c") |> select(block_size)

genes_ribosome <- getCombinedTermGenes(testdf, term_ribosome)
genes_translation <- getCombinedTermGenes(testdf, term_translation)
length(genes_ribosome) # 35/69
length(genes_translation) # 5/69
sum(genes_translation %in% genes_ribosome)

# hardcoding c
GOtable25[GOtable25$CCM_name == "c", "GOslim_combined_terms"] <- 
  paste("ribosome components & rRNA processing (35/69)",
        sep = ",\n\n")

# d, respiratory chain complexes, MICOS complex, ATP synthase complex coupling factor F(o)
testdf <- getGOSlimDf("d")
moduledf25 |> filter(module_name == "d") |> select(block_size)

genes_cellularRespiration <- getCombinedTermGenes(testdf, term_cellularRespiration)
length(genes_cellularRespiration) # 48/484

genes_mitochondrialTranslation <- getCombinedTermGenes(testdf, term_mitochondrialTranslation)
length(genes_mitochondrialTranslation) # 23/484

genes_mitochondria <- getCombinedTermGenes(testdf, term_mitochondria)
length(genes_mitochondria) # 243/484

genes_DNA <- getCombinedTermGenes(testdf, term_DNArepairMetabolismBinding)
length(genes_DNA)
sum(genes_DNA %in% c(genes_cellularRespiration, genes_mitochondrialTranslation, genes_mitochondria))

sum(c(genes_cellularRespiration, genes_mitochondrialTranslation) %in% genes_mitochondria) # nearly all of them
sum(genes_cellularRespiration %in% genes_mitochondrialTranslation)

# term heiarchy
# cellular respiration = mitochondrial translation > DNA = other mitochondria

# hardcoding d
GOtable25[GOtable25$CCM_name == "d", "GOslim_combined_terms"] <- 
  paste("cellular respiration (48/484)",
        "mitochondrial translation (23/484)",
        "mitochondrial organization & other mitochondrial genes (182/484)",
        "DNA binding & repair (23/484)",
        sep = ",\n\n")

# e, DNA repair/metabolism/binding, cytoskeleton, some RNA processing, exocytosis/endocytosis
testdf <- getGOSlimDf("e")
moduledf25 |> filter(module_name == "e") |> select(block_size)

genes_cytoskeleton <- getCombinedTermGenes(testdf, term_cytoskeleton)
length(genes_cytoskeleton) # 96/534

genes_cellCycle <- getCombinedTermGenes(testdf, term_cellCycle)
length(genes_cellCycle) # 45/534

genes_Endo <- getCombinedTermGenes(testdf, term_endocytosis)
length(genes_Endo) # 46/534

sum(genes_cytoskeleton %in% genes_cellCycle)
sum(genes_cytoskeleton %in% genes_Endo)
sum(genes_cellCycle %in% genes_Endo)

genes_RNAsplicingAndMetabolism <- getCombinedTermGenes(testdf, term_RNAsplicingAndMetabolism)
length(genes_RNAsplicingAndMetabolism) # 27/534

genes_chromatinDNArepair <- getCombinedTermGenes(testdf, c(term_chromatin, term_DNArepairMetabolismBinding))
length(genes_chromatinDNArepair) # 116/534

genes_transport <- getCombinedTermGenes(testdf, term_vacuoleGolgiTransport)
length(genes_transport) # 164/534

sum(genes_Endo %in% genes_transport)

# how many genes are in multiple terms?
table(c(genes_cytoskeleton,
        genes_cellCycle,
        genes_Endo,
        genes_RNAsplicingAndMetabolism, 
        genes_chromatinDNArepair,
        genes_transport)) |> table() # fairly few, but there is a decent subset with 2 terms

overlapdf <- bind_rows(tibble("gene_name" = genes_cytoskeleton,
                              "term_name" = "cytoskeleton"),
                       tibble("gene_name" = genes_cellCycle,
                              "term_name" = "cellCycle"),
                       tibble("gene_name" = genes_Endo,
                              "term_name" = "Endo"),
                       tibble("gene_name" = genes_RNAsplicingAndMetabolism,
                              "term_name" = "RNAsplicingAndMetabolism"),
                       tibble("gene_name" = genes_chromatinDNArepair,
                              "term_name" = "chromatinDNArepair"),
                       tibble("gene_name" = genes_transport,
                              "term_name" = "transport"),) |> 
  pivot_wider(id_cols = "gene_name", names_from = "term_name", values_from = "term_name")
overlapdf$all_terms <- paste(overlapdf$cytoskeleton,
                             overlapdf$cellCycle,
                             overlapdf$Endo, 
                             overlapdf$RNAsplicingAndMetabolism,
                             overlapdf$transport,
                             overlapdf$chromatinDNArepair, sep = "_")
table(overlapdf$all_terms)

# term heiarchy
# endocytosis = cell cycle > cytoskeleton = transport = RNA splicing = DNA repair/remodel

# hardcoding e
GOtable25[GOtable25$CCM_name == "e", "GOslim_combined_terms"] <- 
  paste("endocytosis (44/534)",
        "vacuole/vesicle transport (95/534)",
        "cytoskeleton (51/534)",
        "DNA repair, chromatin remodelling & regulation of cell cycle (86/534)",
        sep = ",\n\n")

# f, same as e but also including histone remodelling (SWI/SNF and RPD3l complexes)
#     and ubiquitin (SCF and ligase) complexes
testdf <- getGOSlimDf("f")
moduledf25 |> filter(module_name == "f") |> select(block_size)

genes_ExoEndo <- getCombinedTermGenes(testdf, term_exoendocytosis)
length(genes_ExoEndo) # 20/354

genes_cytoskeleton <- getCombinedTermGenes(testdf, term_cytoskeleton)
length(genes_cytoskeleton) # 65/354

genes_chromatinDNArepairCellCycle <- getCombinedTermGenes(testdf, c(term_chromatin, term_DNArepairMetabolismBinding, term_cellCycle))
length(genes_chromatinDNArepairCellCycle) # 81/354

genes_transport <- getCombinedTermGenes(testdf, term_vacuoleGolgiTransport)
length(genes_transport) # 75/354

# how many genes are in multiple terms?
table(c(genes_ExoEndo, 
        genes_cytoskeleton,
        genes_chromatinDNArepair,
        genes_transport)) |> table() 

overlapdf <- bind_rows(tibble("gene_name" = genes_cytoskeleton,
                              "term_name" = "cytoskeleton"),
                       tibble("gene_name" = genes_ExoEndo,
                              "term_name" = "ExoEndo"),
                       tibble("gene_name" = genes_chromatinDNArepairCellCycle,
                              "term_name" = "chromatinDNArepairCellCycle"),
                       tibble("gene_name" = genes_transport,
                              "term_name" = "transport")) |> 
  pivot_wider(id_cols = "gene_name", names_from = "term_name", values_from = "term_name")
overlapdf$all_terms <- paste(overlapdf$cytoskeleton,
                             overlapdf$ExoEndo,
                             overlapdf$chromatinDNArepairCellCycle,
                             overlapdf$transport, sep = "_")
table(overlapdf$all_terms)

# term heiarchy
# Exoendo > cytoskeleton > chromatin DNA Cell Cycle = transport

# hardcoding f
GOtable25[GOtable25$CCM_name == "f", "GOslim_combined_terms"] <- 
  paste("exo/endocytosis (20/354)",
        "vacuole/vesicle transport (37/354)",
        "cytoskeleton (44/534)",
        "DNA repair, chromatin remodelling & regulation of cell cycle (62/354)",
        sep = ",\n\n")

# g, proteasome, RNA splicing, 
testdf <- getGOSlimDf("g")
moduledf25 |> filter(module_name == "g") |> select(block_size)

# RNA splicing and metabolism
genes_RNAsplicingAndMetabolism <- getCombinedTermGenes(testdf, term_RNAsplicingAndMetabolism)
length(genes_RNAsplicingAndMetabolism) # 72/440

# mRNA degradation (P-body, Lsm1-7-Pat1 complex)
genes_mRNADegrading <- getCombinedTermGenes(testdf, term_mRNADegrading)
length(genes_mRNADegrading) # 17/440

# vacuole transport
genes_vacuoleGolgiTransport <- getCombinedTermGenes(testdf, term_vacuoleGolgiTransport)
length(genes_vacuoleGolgiTransport) # 114/440

# protein folding
genes_protein <- getCombinedTermGenes(testdf, term_chaperone)
length(genes_protein) # 68/440

# proteasome
genes_proteasome <- getCombinedTermGenes(testdf, term_proteasome)
length(genes_proteasome) # 16/440

sum(genes_proteasome %in% genes_mRNADegrading)

table(c(genes_RNAsplicingAndMetabolism,
        genes_mRNADegrading,
        genes_vacuoleGolgiTransport,
        genes_protein,
        genes_proteasome)) |> table() 

overlapdf <- bind_rows(tibble("gene_name" = genes_RNAsplicingAndMetabolism,
                              "term_name" = "RNAsplicing"),
                       tibble("gene_name" = genes_mRNADegrading,
                              "term_name" = "mRNAdegrade"),
                       tibble("gene_name" = genes_vacuoleGolgiTransport,
                              "term_name" = "transport"),
                       tibble("gene_name" = genes_proteasome,
                              "term_name" = "proteasome"),
                       tibble("gene_name" = genes_protein,
                              "term_name" = "protein")) |> 
  pivot_wider(id_cols = "gene_name", names_from = "term_name", values_from = "term_name")
overlapdf$all_terms <- paste(overlapdf$RNAsplicing,
                             overlapdf$mRNAdegrade,
                             overlapdf$transport,
                             overlapdf$protein,
                             overlapdf$proteasome, sep = "_")
table(overlapdf$all_terms)

# term heiarchy
# mRNA degrade = proteasome > protein folding = RNA splicing > transport

# hardcoding g
GOtable25[GOtable25$CCM_name == "g", "GOslim_combined_terms"] <- 
  paste("mRNA degredation (17/440)",
        "mRNA splicing (54/440)",
        "proteasome complex (16/440)",
        "protein folding (60/440)",
        "vacuole/Golgi transport (87/440)",
        sep = ",\n\n")

# h, glucose metabolism/low glucose growth, mitochondrial organization, mitochondrial proteins
testdf <- getGOSlimDf("h")
moduledf25 |> filter(module_name == "h") |> select(block_size)

genes_metabolism <- getCombinedTermGenes(testdf, term_metabolism)
length(genes_metabolism) # 93/166

genes_mitochondria <- getCombinedTermGenes(testdf, term_mitochondria)
length(genes_mitochondria) # 27/166

# goslim doesn't have a glycolysis category specifically,
# but these are the known glycolytic enzymes
glycolgenes <- c("YBR196C", "YGR240C", "YMR205C", "YLR377C", "YKL060C", "YDR050C", 
                "YJL052W", "YJR009C", "YGR192C", "YCR012W", "YKL152C", "YHR174W", 
                "YGR254W", "YOR347C", "YAL038W")
sum(glycolgenes %in% genes_metabolism)/length(glycolgenes)

sum(genes_mitochondria %in% genes_metabolism)

# hardcoding h
GOtable25[GOtable25$CCM_name == "h", "GOslim_combined_terms"] <- 
  paste("glycolytic enzymes (12/166)",
        "mitochondrial metabolism (27/166)",
        "other metabolism (81/166)",
        sep = ",\n\n")

# i, small/large ribosomal subunits
testdf <- getGOSlimDf("i")
moduledf25 |> filter(module_name == "i") |> select(block_size)

genes_ribosome <- getCombinedTermGenes(testdf, term_ribosome)
length(genes_ribosome) # 50/83

genes_translation <- getCombinedTermGenes(testdf, term_translation)
length(genes_translation) # 5/83

genes_transcription <- getCombinedTermGenes(testdf, term_transcription)
length(genes_transcription) # 7/83

genes_tRNA <- getCombinedTermGenes(testdf, term_tRNA)
length(genes_tRNA) # 13/83

table(c(genes_ribosome,
        genes_translation,
        genes_tRNA,
        genes_transcription)) |> table()

overlapdf <- bind_rows(tibble("gene_name" = genes_ribosome,
                              "term_name" = "ribosome"),
                       tibble("gene_name" = genes_translation,
                              "term_name" = "translation"),
                       tibble("gene_name" = genes_tRNA,
                              "term_name" = "tRNA"),
                       tibble("gene_name" = genes_transcription,
                              "term_name" = "transcription")) |> 
  pivot_wider(id_cols = "gene_name", names_from = "term_name", values_from = "term_name")
overlapdf$all_terms <- paste(overlapdf$ribosome,
                             overlapdf$translation,
                             overlapdf$tRNA,
                             overlapdf$transcription, sep = "_")
table(overlapdf$all_terms)

sum(genes_stress %in% genes_chaperone)
# term heiarchy:
# stress = chaperone > tRNA > ribosome, translation not worth including

# hardcoding i
GOtable25[GOtable25$CCM_name == "i", "GOslim_combined_terms"] <- 
  paste("ribosome components (50/83)",
        sep = ",\n\n")

# j, small/large ribosomal subunits, stress
testdf <- getGOSlimDf("j")
moduledf25 |> filter(module_name == "j") |> select(block_size)

genes_ribosome <- getCombinedTermGenes(testdf, term_ribosome)
length(genes_ribosome) # 111/154

genes_translation <- getCombinedTermGenes(testdf, term_translation)
length(genes_translation) # 99/154

sum(genes_translation %in% genes_ribosome)

genes_transcription <- getCombinedTermGenes(testdf, term_transcription)
length(genes_transcription) # 0/154

genes_tRNA <- getCombinedTermGenes(testdf, term_tRNA)
length(genes_tRNA) # 9/154

genes_stress <- getCombinedTermGenes(testdf, "cytoplasmic stress granule")
length(genes_stress)

table(c(genes_ribosome,
        genes_stress)) |> table()

# term heiarchcy
# stress granule > ribosome

# hardcoding j
GOtable25[GOtable25$CCM_name == "j", "GOslim_combined_terms"] <- 
  paste("ribosome components (106/166)",
        "stress granule (12/166)",
        sep = ",\n\n")

#### printing publication table ####
GOSummary <- GOtable25 |> 
  select(CCM_name, GOslim_combined_terms) |> 
  pivot_wider(names_from = "CCM_name",
              values_from = "GOslim_combined_terms")
GOSummary
library(huxtable)
# first half
halfway_idx <- ceiling(ncol(GOSummary)/2)
ht1 <- hux(GOSummary[,1:halfway_idx]) %>%
  set_font_size(20)                               %>%
  set_header_rows(1, TRUE)                             %>%
  set_width(0.9)                                       %>%
  set_tb_padding(2)                                    %>%
  set_caption("Huxtable properties")                   %>%
  set_label("tab:props")                               %>%
  theme_bright(
    header_rows = TRUE,
    header_cols = FALSE,
    colors = sapply(names(GOSummary)[1:halfway_idx], \(m) {
      output <- moduledf25[moduledf25$module_name == m, "CCM_color"] |> as.character()
      return(output)
    })
  )

font(ht1) <- "monospace"
caption(ht1) <- "<span style='font-family: monospace'>GO Enrichment Summary</span>"

quick_html(ht1, file = "../../aligning_the_molecular_phenotype/paper_figures/GeneOntology/GOSummary1.html")

# second half
ht2 <- hux(GOSummary[,(halfway_idx + 1):ncol(GOSummary)]) %>%
  set_font_size(20)                               %>%
  set_header_rows(1, TRUE)                             %>%
  set_width(0.9)                                       %>%
  set_tb_padding(2)                                    %>%
  set_caption("Huxtable properties")                   %>%
  set_label("tab:props")                               %>%
  theme_bright(
    header_rows = TRUE,
    header_cols = FALSE,
    colors = sapply(names(GOSummary)[(halfway_idx + 1):ncol(GOSummary)], \(m) {
      output <- moduledf25[moduledf25$module_name == m, "CCM_color"] |> as.character()
      return(output)
    })
  )

font(ht2) <- "monospace"
caption(ht2) <- "<span style='font-family: monospace'>GO Enrichment Summary</span>"

quick_html(ht2, file = "../../aligning_the_molecular_phenotype/paper_figures/GeneOntology/GOSummary2.html")

#### GO for TF-regulated ccm subsets? ####

# pro: easier than whole modules, you can literally check what every gene does
# con: no one is asking you to do this

##################### Archive ##################### 



                       