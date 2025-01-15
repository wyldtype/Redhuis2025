sapply(c("dplyr", "readr", "tidyr", "purrr", "ggplot2", "openxlsx"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")

# read in data
counts_cer <- read.xlsx("data_files/downloaded_from_Krieger2020/Supplemental_Table_S11.xlsx", sheet = 2, na.strings = "")
counts_par <- read.xlsx("data_files/downloaded_from_Krieger2020/Supplemental_Table_S12.xlsx", sheet = 2, na.strings = "")

cat("percent count data columns with same names in both species: ",
  sum(colnames(counts_cer)==colnames(counts_par))/ncol(counts_cer))

sample_info <- read.xlsx("data_files/downloaded_from_Krieger2020/bioSample1to999.xlsx", na.strings="not applicable", cols=c(1,4,9,13,14,15,17)) %>%
  bind_rows(read.xlsx("data_files/downloaded_from_Krieger2020/bioSample1000toEnd.xlsx", na.strings="not applicable", cols=c(1,4,9,13,14,15,17)))
colnames(sample_info) <- c("sample_name", "organism" , "collection_date", "genotype", "experiment","time_point", "well_flask_ID")
cat("percent sample info columns with same names in count data (before matching): ",
    sum(colnames(counts_cer)[-1]==sample_info$sample_name)/ncol(counts_cer[,-1]))
desired_sample_order <- sapply(c(1:nrow(sample_info)), function(i) {
  return(which(sample_info$sample_name==colnames(counts_cer[,-1])[i]))
})
sample_info <- sample_info[desired_sample_order,]
cat("percent sample info columns with same names in count data (after matching): ",
  sum(colnames(counts_cer)[-1]==sample_info$sample_name)/ncol(counts_cer[,-1]))

#### Calculating mapping bias as % Scer reads ####
# normalizing function for filtering out lowly expressed 
# genes prior to assessing mapping bias
# (used later to actually normalize count data)
# normalizing counts to adjust for differences in library size
# sums .cts_cer and .cts_par to get library size, only returns
# counts for specified allele
# @input: count matrix (genes are rows, columns are samples)
# @output: a count matrix normalzied for library size---integer counts in counts-per-million
countsPerMillionAllele <- function(.cts_cer, .cts_par, .allele) {
  librarySizes <- colSums(.cts_cer, na.rm = TRUE) + colSums(.cts_par, na.rm = TRUE)
  if (.allele == "cer") {
    .cts <- .cts_cer
  }
  if (.allele == "par") {
    .cts <- .cts_par
  }
  output <- apply(.cts, 1, function(x) {
    normalized <- (x/librarySizes)*1e6
    return(round(normalized))
  })
  return(t(output)) # For some unhinged reason, a vector output of apply across ROWS forms the COLUMNS of a new matrix
}
# tests for countsPerMillionAllele
test_cts <- counts_cer[,-1]
test_rowIdx <- sample(c(1:nrow(test_cts)), 1)
test_colIdx <- sample(c(1:ncol(test_cts)), 1)
test_count <- test_cts[test_rowIdx, test_colIdx]
test_cpm <- countsPerMillionAllele(.cts_cer = counts_cer[,-1],
                                   .cts_par = counts_par[,-1],
                                   .allele = "cer")
((test_count/(colSums(counts_cer[,-1], na.rm = TRUE) + 
                colSums(counts_par[,-1], na.rm = TRUE))[test_colIdx])*1e6) %>% 
  round() # what it should be
test_cts[test_rowIdx, test_colIdx] # what it is before using our function
test_cpm[test_rowIdx, test_colIdx] # what it is using our function

# normalize to counts per million based on total 
# library size: cer reads + par reads regardless of sample organism
cpm_cer <- countsPerMillionAllele(.cts_cer = counts_cer[,-1],
                                  .cts_par = counts_par[,-1],
                                  .allele = "cer")
cpm_par <- countsPerMillionAllele(.cts_cer = counts_cer[,-1],
                                  .cts_par = counts_par[,-1],
                                  .allele = "par")
# 2) filter lowly expressed: < 30 cpm
isHighExpr <- (rowMeans(cpm_cer + cpm_par) > 30) |> sapply(FUN = isTRUE)
common_genes <- counts_cer[isHighExpr, 1]
cpm_cer <- cpm_cer[isHighExpr,]
cpm_par <- cpm_par[isHighExpr,]

# 3) check % cer of all high-enough expressed genes is close to 1 for cer samples and 0 for par samples
plotdf <- bind_rows(bind_cols(tibble(gene = common_genes,
                                     allele = "cer"), cpm_cer),
                    bind_cols(tibble(gene = common_genes,
                                     allele = "par"), cpm_par)) |> 
  pivot_longer(cols = colnames(counts_cer[, -1]),
               names_to = c("sample_name"),
               values_to = "count") |> 
  pivot_wider(id_cols = c("sample_name", "gene"),
              values_from = "count", names_from = "allele",
              names_prefix = "counts_")
plotdf$organism <- if_else(grepl("_cer_", plotdf$sample_name),
                           true = "cerSample", 
                           false = if_else(grepl("_par_", plotdf$sample_name),
                                           true = "parSample",
                                           false = "hybSample"))
plotdf$pct_cer <- plotdf$counts_cer/(plotdf$counts_cer + plotdf$counts_par + 1)
sample_genes <- sample(plotdf$gene, size = 100)
ggplot(filter(plotdf, gene %in% sample_genes),
       aes(x = gene, y = pct_cer)) + 
  geom_point(aes(color = organism))

# Do any cer/par samples have many incorrect gene mappings?
sampdf <- plotdf |> filter(organism %in% c("cerSample", "parSample")) |> 
  group_by(sample_name, organism) |> 
  summarise(avg_pct_cer = mean(pct_cer, na.rm = TRUE)) |> 
  filter((organism == "cerSample" & avg_pct_cer < 0.75) |
           (organism == "parSample" & avg_pct_cer > 0.25))
sampdf |> print(n = nrow(sampdf))
bad_samples <- sampdf$sample_name

# Do any genes have consistently biased mappings? These might 
# be genes with local homology btwn cer and par and therefore
# shouldn't be used for allele-specific comparisons in the hybrid
cer_genedf <- plotdf |> filter(organism == "cerSample" & !(sample_name %in% bad_samples)) |> 
  group_by(gene) |> 
  summarise(avg_pct_cer = mean(pct_cer, na.rm = TRUE))
hist(pull(select(filter(cer_genedf, avg_pct_cer < 0.9), avg_pct_cer)))
par_genedf <- plotdf |> filter(organism == "parSample" &
                                 !(sample_name %in% bad_samples)) |> 
  group_by(gene) |> 
  summarise(avg_pct_cer = mean(pct_cer, na.rm = TRUE))
hist(pull(select(filter(par_genedf, avg_pct_cer > 0.1), avg_pct_cer)))
cer_biased_genes <- filter(par_genedf, avg_pct_cer > 0.1) |> select(gene) |> pull()
par_biased_genes <- filter(cer_genedf, avg_pct_cer < 0.9) |> select(gene) |> pull()

# 4) check if cer/par ratio of parents matches hybrids for each gene across samples
plotdf$bias <- if_else(plotdf$gene %in% cer_biased_genes,
                       true = "cer",
                       false = if_else(plotdf$gene %in% par_biased_genes,
                                       true = "par", false = "none"))
plotdf_vshyb <- plotdf |> group_by(gene, organism, bias) |> 
  summarise(mean_counts_cer = mean(counts_cer, na.rm = TRUE),
            mean_counts_par = mean(counts_par, na.rm = TRUE)) |> 
  drop_na() |> 
  pivot_wider(id_cols = c("gene", "bias"), 
              values_from = c("mean_counts_cer", "mean_counts_par"),
              names_from = organism)

# two ways to calculate read counts for each gene in the parents:
#    1) sum the hits for the cer allele and the par allele (after all we know that all these reads came from one allele or the other)
#    2) only count the hits for the correct parent's allele
# in the hybrid, we no longer know which allele is correct, because both alleles are present.
# This is generally equivalent to summing the parental alleles in the parents
# So if we see a difference in how correlated the parental vs hybrid %cer values
# between the following two graphs, there might be significant mapping bias
plotdf_vshyb$par_sum <- plotdf_vshyb$mean_counts_cer_parSample + plotdf_vshyb$mean_counts_par_parSample
plotdf_vshyb$cer_sum <- plotdf_vshyb$mean_counts_cer_cerSample + plotdf_vshyb$mean_counts_par_cerSample
plotdf_vshyb$pct_cer_parents_sum <- plotdf_vshyb$cer_sum/(plotdf_vshyb$par_sum + plotdf_vshyb$cer_sum)
plotdf_vshyb$pct_cer_parents_singleAllele <- plotdf_vshyb$mean_counts_cer_cerSample/(plotdf_vshyb$mean_counts_cer_cerSample + plotdf_vshyb$mean_counts_par_parSample)
plotdf_vshyb$pct_cer_hybrids_sum <- plotdf_vshyb$mean_counts_cer_hybSample/(plotdf_vshyb$mean_counts_cer_hybSample + plotdf_vshyb$mean_counts_par_hybSample)
# summing parental allele reads
ggplot(plotdf_vshyb, aes(x = pct_cer_parents_sum, y = pct_cer_hybrids_sum)) +
  geom_point(aes(color = bias), alpha = 0.5)
ggplot(filter(plotdf_vshyb, bias != "none"), aes(x = pct_cer_parents_sum, y = pct_cer_hybrids_sum)) +
  geom_point(aes(color = bias), alpha = 0.5)
# single parental allele reads
ggplot(plotdf_vshyb, aes(x = pct_cer_parents_singleAllele, y = pct_cer_hybrids_sum)) +
  geom_point(aes(color = bias), alpha = 0.5)
ggplot(filter(plotdf_vshyb, bias != "none"), aes(x = pct_cer_parents_singleAllele, y = pct_cer_hybrids_sum)) +
  geom_point(aes(color = bias), alpha = 0.5)

# Conclusion: very few genes with mapping bias, we can remove them

#### Comparing with --alignEndsType EndToEnd alignment ####
# do the same genes show up as biased with a different alignment
# strategy?
cer_genedf_barkai <- cer_genedf
par_genedf_barkai <- par_genedf

tagseq <- list.files("data_files/tagseq_counts/", full.names = TRUE) |> 
  map(read_table, col_names = FALSE, show_col_types = FALSE) |> 
  map(.f = select, X1, X3) |> # X1 are gene names, X3 is the sense strand read count
  purrr::reduce(.f = \(x, y) full_join(x = x, y = y, by = "X1"))

colnames(tagseq) <- c("gene", gsub("_ReadsPerGene.out.tab", "", list.files("data_files/tagseq_counts/", full.names = FALSE)))
QCdf <- tagseq[grepl("N_", tagseq$gene), ]
tagseq <- tagseq[!grepl("N_", tagseq$gene),]
tagseq_cer <- tagseq[grepl("^cer_", tagseq$gene),]
tagseq_par <- tagseq[grepl("^par_", tagseq$gene),]
common_genes <- intersect(gsub("^cer_", "", tagseq_cer$gene),
                          gsub("^par_", "", tagseq_par$gene))
tagseq_cer$gene <- gsub("^cer_", "", tagseq_cer$gene)
tagseq_par$gene <- gsub("^par_", "", tagseq_par$gene)
tagseq_cer <- tagseq_cer[sapply(common_genes, \(x) which(x == tagseq_cer$gene)),]
tagseq_par <- tagseq_par[sapply(common_genes, \(x) which(x == tagseq_par$gene)),]
sum(tagseq_cer$gene == tagseq_par$gene)
length(common_genes)

# normalize to counts per million based on total 
# library size: cer reads + par reads regardless of sample organism
cpm_cer <- countsPerMillionAllele(.cts_cer = tagseq_cer[,-1],
                                  .cts_par = tagseq_par[,-1],
                                  .allele = "cer")
cpm_par <- countsPerMillionAllele(.cts_cer = tagseq_cer[,-1],
                                  .cts_par = tagseq_par[,-1],
                                  .allele = "par")
# 2) filter lowly expressed: < 30 cpm
isHighExpr <- (rowMeans(cpm_cer + cpm_par) > 30)
keep_genes <- common_genes[isHighExpr]
cpm_cer <- cpm_cer[isHighExpr,]
cpm_par <- cpm_par[isHighExpr,]

# 3) check % cer of all high-enough expressed genes is close to 1 for cer samples and 0 for par samples
plotdf <- bind_rows(bind_cols(tibble(gene = keep_genes,
                                     allele = "cer"), cpm_cer),
                    bind_cols(tibble(gene = keep_genes,
                                     allele = "par"), cpm_par)) |> 
  pivot_longer(cols = colnames(tagseq_cer[,-c(1,2)]),
               names_to = c("sample_name"),
               values_to = "count") |> 
  pivot_wider(id_cols = c("sample_name", "gene"),
              values_from = "count", names_from = "allele",
              names_prefix = "counts_")
plotdf$organism <- if_else(grepl("_cer_", plotdf$sample_name),
                           true = "cerSample", 
                           false = if_else(grepl("_par_", plotdf$sample_name),
                                           true = "parSample",
                                           false = "hybSample"))
plotdf$pct_cer <- plotdf$counts_cer/(plotdf$counts_cer + plotdf$counts_par + 1)
sample_genes <- sample(plotdf$gene, size = 100)
ggplot(filter(plotdf, gene %in% sample_genes),
       aes(x = gene, y = pct_cer)) + 
  geom_point(aes(color = organism))

# Do any cer/par samples have many incorrect gene mappings?
sampdf <- plotdf |> filter(organism %in% c("cerSample", "parSample")) |> 
  group_by(sample_name, organism) |> 
  summarise(avg_pct_cer = mean(pct_cer, na.rm = TRUE)) |> 
  filter((organism == "cerSample" & avg_pct_cer < 0.75) |
           (organism == "parSample" & avg_pct_cer > 0.25))
sampdf |> print(n = nrow(sampdf))
bad_samples <- sampdf$sample_name

# Do any genes have consistently biased mappings? These might 
# be genes with local homology btwn cer and par and therefore
# shouldn't be used for allele-specific comparisons in the hybrid
cer_genedf <- plotdf |> filter(organism == "cerSample" & !(sample_name %in% bad_samples)) |> 
  group_by(gene) |> 
  summarise(avg_pct_cer = mean(pct_cer, na.rm = TRUE))
hist(pull(select(filter(cer_genedf, avg_pct_cer < 0.9), avg_pct_cer)))
par_genedf <- plotdf |> filter(organism == "parSample" &
                                 !(sample_name %in% bad_samples)) |> 
  group_by(gene) |> 
  summarise(avg_pct_cer = mean(pct_cer, na.rm = TRUE))
hist(pull(select(filter(par_genedf, avg_pct_cer > 0.1), avg_pct_cer)))

cer_genedf_endToEnd <- cer_genedf
par_genedf_endToEnd <- par_genedf

cer_genedf_endToEnd |> filter(avg_pct_cer < 0.9) # 37 genes
par_genedf_endToEnd |> filter(avg_pct_cer > 0.1) # 54 genes
cer_genedf_barkai |> filter(avg_pct_cer < 0.9) # 41 genes
par_genedf_barkai |> filter(avg_pct_cer > 0.1) # 54 genes

#### Comparing ####

sample_genedf <- plotdf # % cer of each gene in each sample (massive)
plotdf <- left_join(tibble(gene = cer_genedf_barkai$gene,
                           pct_cer_cerSample_barkai = cer_genedf_barkai$avg_pct_cer,
                           pct_cer_parSample_barkai = par_genedf_barkai$avg_pct_cer,
                           alignment = "barkai"),
                    tibble(gene = cer_genedf_endToEnd$gene,
                           pct_cer_cerSample_endToEnd = cer_genedf_endToEnd$avg_pct_cer,
                           pct_cer_parSample_endToEnd = par_genedf_endToEnd$avg_pct_cer,
                           alignment = "endToEnd"),
                    by = "gene")
plotdf$flagged_cerSample <- if_else(plotdf$pct_cer_cerSample_barkai < 0.9,
                                    true = if_else(plotdf$pct_cer_cerSample_endToEnd < 0.9,
                                                   true = "both",
                                                   false = "barkai"),
                                    false = if_else(plotdf$pct_cer_cerSample_endToEnd < 0.9,
                                                    true = "endToEnd",
                                                    false = "neither"))
plotdf$flagged_parSample <- if_else(plotdf$pct_cer_parSample_barkai > 0.1,
                                    true = if_else(plotdf$pct_cer_parSample_endToEnd > 0.1,
                                                   true = "both",
                                                   false = "barkai"),
                                    false = if_else(plotdf$pct_cer_parSample_endToEnd > 0.1,
                                                    true = "endToEnd",
                                                    false = "neither"))
# cer samples with par bias
p_cerSamples_newAlignment <- ggplot(plotdf, aes(x = pct_cer_cerSample_barkai,
                   y = pct_cer_cerSample_endToEnd)) +
  geom_point(aes(color = flagged_cerSample))
# par samples with cer bias
p_parSamples_newAlignment <- ggplot(plotdf, aes(x = pct_cer_parSample_barkai,
                   y = pct_cer_parSample_endToEnd)) +
  geom_point(aes(color = flagged_parSample))

# # (repeating with misisng29hybaligned counts to see how much soft-clipping affects bias)
# p_cerSamples_oldAlignment <- ggplot(plotdf, aes(x = pct_cer_cerSample_barkai, 
#                                                 y = pct_cer_cerSample_softClipped)) +
#   geom_point(aes(color = flagged_cerSample))
# # par samples with cer bias
# p_parSamples_oldAlignment <- ggplot(plotdf, aes(x = pct_cer_parSample_barkai, 
#                                                 y = pct_cer_parSample_softClipped)) +
#   geom_point(aes(color = flagged_parSample))

library(ggpubr)
pdf(file = "../../aligning_the_molecular_phenotype/paper_figures/Supplement/mapping_bias.pdf",
    width = 11, height = 7)
ggarrange(p_cerSamples_newAlignment, p_cerSamples_oldAlignment,
          p_parSamples_newAlignment, p_parSamples_oldAlignment,
          nrow = 2, ncol = 2)
dev.off()
