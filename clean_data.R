sapply(c("dplyr", "readr", "tidyr", "purrr", "ggplot2", "openxlsx", "matrixStats"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2025")

#### Checking alignment stats ####
# tagseq
tagseq_stats <- read_tsv("data_files/QC/QC_tag/star_alignment_plot.tsv")
tagseq_general <- read_tsv("data_files/QC/QC_tag/general_stats_table.tsv")
mean(tagseq_stats$`Uniquely mapped`)
tagseq_general[!grepl("clipped", tagseq_general$Sample), "% Aligned"] |> 
  unlist() |> as.numeric() |> mean(na.rm = TRUE)
# fay rnaseq
fay_stats <- read_tsv("data_files/QC/QC_fay/star_alignment_plot.tsv")
fay_general <- read_tsv("data_files/QC/QC_fay/general_stats_table.tsv")
mean(fay_stats$`Uniquely mapped`)
fay_general[!grepl("clipped", fay_general$Sample), "% Aligned"] |> 
  unlist() |> as.numeric() |> mean(na.rm = TRUE)

#### Reading in Krieger et al. 2020 and Lupo et al. 2021 tagseq data ####
# read in count data (needs to be read in, organized, and normalized)
# (copied from QC_comparingAlignmentsAndQuantifications, which compared this
# alignment to the Krieger et al. 2020 alignment)
tagseq <- list.files("data_files/tagseq_counts/", full.names = TRUE) |> 
  map(read_table, col_names = FALSE, show_col_types = FALSE) |> 
  map(.f = select, X1, X3) |> # X1 are gene names, X3 is the sense strand read count
  purrr::reduce(.f = \(x, y) full_join(x = x, y = y, by = "X1"))

colnames(tagseq) <- c("gene", gsub("_ReadsPerGene.out.tab", "", list.files("data_files/tagseq_counts/", full.names = FALSE)))
QCdf <- tagseq[grepl("N_", tagseq$gene), ]
tagseq <- tagseq[!grepl("N_", tagseq$gene),]
tagseq <- tagseq[!tagseq$gene %in% c("cer_NA", "par_NA"),]
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

#### Switching WT_39_CellCycle_rep2 cer and hyb samples ####
# Based on QC_comparingAlignmentsAndQuantification.R, 
# the cer and hyb counts for WT 39 rep2 have been mislabeled 
# and need to be switched
real_cer_39_counts_cerAllele <- tagseq_cer[,"WT_39_hyb_CellCycle_rep2"]
real_hyb_39_counts_cerAllele <- tagseq_cer[,"WT_39_cer_CellCycle_rep2"]
real_cer_39_counts_parAllele <- tagseq_par[,"WT_39_hyb_CellCycle_rep2"]
real_hyb_39_counts_parAllele <- tagseq_par[,"WT_39_cer_CellCycle_rep2"]

tagseq_cer[,"WT_39_hyb_CellCycle_rep2"] <- real_hyb_39_counts_cerAllele
tagseq_cer[,"WT_39_cer_CellCycle_rep2"] <- real_cer_39_counts_cerAllele
tagseq_par[,"WT_39_hyb_CellCycle_rep2"] <- real_hyb_39_counts_parAllele
tagseq_par[,"WT_39_cer_CellCycle_rep2"] <- real_cer_39_counts_parAllele

rm(real_hyb_39_counts_cerAllele, real_cer_39_counts_cerAllele,
   real_hyb_39_counts_parAllele, real_cer_39_counts_parAllele)

#### calculating % mapping to each allele for parents vs hybrids #### 
# Note: this can only be done on tagseq samples, 
# as Fay et al. 2023 pooled Scer and Spar parental samples

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
test_cts <- tagseq_cer[,-1]
test_rowIdx <- sample(c(1:nrow(test_cts)), 1)
test_colIdx <- sample(which(grepl("cer", colnames(test_cts))), 1)
test_count <- test_cts[test_rowIdx, test_colIdx]
test_cpm <- countsPerMillionAllele(.cts_cer = tagseq_cer[,-1], 
                                   .cts_par = tagseq_par[,-1],
                                   .allele = "cer")
((test_count/(colSums(tagseq_cer[,-1], na.rm = TRUE) + 
                colSums(tagseq_par[,-1], na.rm = TRUE))[test_colIdx])*1e6) %>% 
  round() # what it should be
test_cts[test_rowIdx, test_colIdx] # what it is before using our function
test_cpm[test_rowIdx, test_colIdx] # what it is using our function

# normalize to counts per million based on total 
# library size: cer reads + par reads regardless of sample organism
cpm_cer <- countsPerMillionAllele(.cts_cer = tagseq_cer[,-1],
                                      .cts_par = tagseq_par[,-1],
                                      .allele = "cer")
cpm_par <- countsPerMillionAllele(.cts_cer = tagseq_cer[,-1],
                                      .cts_par = tagseq_par[,-1],
                                      .allele = "par")
# 2) filter lowly expressed: < 30 cpm
sum(cpm_cer == 0 & cpm_par == 0)
isHighExpr <- (rowMeans(cpm_cer + cpm_par) > 30) |> sapply(FUN = isTRUE)
keep_genes <- common_genes[isHighExpr]
cpm_cer <- cpm_cer[isHighExpr,]
cpm_par <- cpm_par[isHighExpr,]
sum(cpm_cer == 0 & cpm_par == 0) # note there are still individual samples with zero counts

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
# Calculating % of reads mapping to the Scer allele 
# for each gene/sample
# (So % Spar is 1 - % Scer)
plotdf$pct_cer <- if_else(plotdf$counts_cer == 0 &
                            plotdf$counts_par == 0,
                          true = NA,
                          false = plotdf$counts_cer/(plotdf$counts_cer + plotdf$counts_par))
plotdf <- drop_na(plotdf)
sample_genes <- sample(plotdf$gene, size = 100)
ggplot(filter(plotdf, gene %in% sample_genes),
       aes(x = gene, y = pct_cer)) + 
  geom_point(aes(color = organism))

### Do any cer/par samples have many incorrect gene mappings?
sampdf <- plotdf |> filter(organism %in% c("cerSample", "parSample")) |> 
  group_by(sample_name, organism) |> 
  summarise(avg_pct_cer = mean(pct_cer, na.rm = TRUE))

sampdf |> filter(organism == "cerSample") |> arrange(avg_pct_cer)
sampdf |> filter(organism == "parSample") |> arrange(desc(avg_pct_cer))

ggplot(sampdf, aes(x = avg_pct_cer)) + 
  geom_density(aes(fill = organism)) + 
  geom_vline(xintercept = 0.75) +
  geom_vline(xintercept = 0.25)
  
filter(sampdf, (organism == "cerSample" & round(avg_pct_cer, digits = 1) < 0.8) |
           (organism == "parSample" & round(avg_pct_cer, digits = 1) > 0.2)) |> 
  arrange(organism, desc(avg_pct_cer))
biased_samples <- filter(sampdf, (organism == "cerSample" & round(avg_pct_cer, digits = 1) < 0.9) |
                           (organism == "parSample" & round(avg_pct_cer, digits = 1) > 0.1)) |> 
  select(sample_name) |> pull()
biased_samples
# do these samples have small library sizes?
bad_libsizes <- colSums(tagseq_cer[,biased_samples] + tagseq_par[,biased_samples])
bad_libsizes |> sort()
hist(colSums(tagseq_cer[,-1] + tagseq_par[,-1]), breaks = 50)
abline(v = bad_libsizes, col = "red")

# Do any genes have consistently biased mappings? These might 
# be genes with above-average seq conservation btwn cer and par and therefore
# shouldn't be used for allele-specific comparisons in the hybrid
cer_genedf <- plotdf |> filter(organism == "cerSample" & !(sample_name %in% biased_samples)) |> 
  group_by(gene) |> 
  summarise(avg_pct_cer = mean(pct_cer, na.rm = TRUE))
hist(cer_genedf$avg_pct_cer, breaks = 50)
abline(v = 0.9, col = "red")
par_genedf <- plotdf |> filter(organism == "parSample" &
                                 !(sample_name %in% biased_samples)) |> 
  group_by(gene) |> 
  summarise(avg_pct_cer = mean(pct_cer, na.rm = TRUE))
hist(par_genedf$avg_pct_cer, breaks = 50)
abline(v = 0.1, col = "red")
cer_biased_genes <- par_genedf |> filter(avg_pct_cer > 0.1) |> select(gene) |> pull()
par_biased_genes <- cer_genedf |> filter(avg_pct_cer < 0.9) |> select(gene) |> pull()
both_biased_genes <- intersect(cer_biased_genes, par_biased_genes)
cer_biased_genes <- setdiff(cer_biased_genes, par_biased_genes)
par_biased_genes <- setdiff(par_biased_genes, cer_biased_genes)
# any common genes of note?
sgd_lookup <- read_tsv("data_files/downloaded_genomes_and_features/SGD_features.tab",
                       col_names = FALSE) |> 
  select(X4, X5) |> unique() |> drop_na()
# cer biased
sgd_lookup[sgd_lookup$X4 %in% cer_biased_genes,] |> 
  print(n = length(cer_biased_genes))
par_genedf |> filter(gene %in% cer_biased_genes) |> 
  arrange(desc(avg_pct_cer)) |> 
  print(n = length(cer_biased_genes))
# par biased
sgd_lookup[sgd_lookup$X4 %in% par_biased_genes,] |> 
  print(n = length(par_biased_genes))
cer_genedf |> filter(gene %in% par_biased_genes) |> 
  arrange(avg_pct_cer) |> 
  print(n = length(par_biased_genes))
# both
sgd_lookup[sgd_lookup$X4 %in% both_biased_genes,] |> 
  print(n = length(both_biased_genes))
cer_genedf |> filter(gene %in% both_biased_genes) |> 
  arrange(avg_pct_cer)
par_genedf |> filter(gene %in% both_biased_genes) |> 
  arrange(desc(avg_pct_cer))
# Are these genes on the low expr end?
cer_genedf <- left_join(cer_genedf, tibble(gene = common_genes[isHighExpr],
                                           mean_expr = rowMeans(cpm_cer[,grepl("cer", colnames(cpm_cer))] +
                                                                  cpm_par[,grepl("cer", colnames(cpm_cer))])),
                        by = "gene")
p_parbias <- ggplot(cer_genedf, aes(x = log2(mean_expr), y = avg_pct_cer)) + 
  geom_point(aes(color = avg_pct_cer < 0.9)) +
  ylab("% reads mapping to Scer allele")
par_genedf <- left_join(par_genedf, tibble(gene = common_genes[isHighExpr],
                                           mean_expr = rowMeans(cpm_cer[,grepl("par", colnames(cpm_cer))] +
                                                                  cpm_par[,grepl("par", colnames(cpm_cer))])),
                        by = "gene")
p_cerbias <- ggplot(par_genedf, aes(x = log2(mean_expr), y = avg_pct_cer)) + 
  geom_point(aes(color = avg_pct_cer > 0.1)) +
  ylab("% reads mapping to Scer allele")
library(ggpubr)
pdf("../../aligning_the_molecular_phenotype/paper_figures/Supplement/mapping_bias_genes.pdf",
    width = 10, height = 5)
ggarrange(p_parbias, p_cerbias, nrow = 1, ncol = 2)
dev.off()

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
# So we can only count them one way. If we see a difference in how 
# correlated the parental vs hybrid %cer values between the 
# two methods of parental read counting, mapping bias
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

# Conclusion: slight differences in the counts of biased
# genes, no visible difference in the counts of unbiased genes
# we can count just reads mapping to correct parent allele
# and remove any biased genes after adding
# Fay et al. 2023 samples and filtering for low expression
cer_biased_genes
par_biased_genes
both_biased_genes
biased_samples

#### Reading in Fay et al. 2023 RNAseq data and R1/R2 QC ####
fay_filenames <- gsub("_ReadsPerGene.out.tab", "", 
                      list.files("data_files/fay_counts/", 
                                 full.names = FALSE))
fay <- list.files("data_files/fay_counts/", full.names = TRUE) |> 
  map(read_table, col_names = FALSE, show_col_types = FALSE) |> 
  map2(.y = fay_filenames,
       .f = \(x, y) {
         if (grepl("R1", y)) {
           output <- select(x, X1, X2)
           return(output)}
         if (grepl("R2", y)) {
           output <- select(x, X1, X3)
           return(output)}
       }) |> # X1 are gene names, X2 is sense strand read count for R1, X3 is the sense strand for R2
  purrr::reduce(.f = \(x, y) full_join(x = x, y = y, by = "X1"))

colnames(fay) <- c("gene", fay_filenames)
QCdf_fay <- fay[grepl("N_", fay$gene), ]
fay <- fay[!grepl("N_", fay$gene),]
fay <- fay[!fay$gene %in% c("cer_NA", "par_NA"),]
fay_cer <- fay[grepl("^cer_", fay$gene),]
fay_par <- fay[grepl("^par_", fay$gene),]
common_genes_fay <- intersect(gsub("^cer_", "", fay_cer$gene),
                          gsub("^par_", "", fay_par$gene))
setequal(common_genes, common_genes_fay) # should be the same set of genes as tagseq
rm(common_genes_fay)
fay_cer$gene <- gsub("^cer_", "", fay_cer$gene)
fay_par$gene <- gsub("^par_", "", fay_par$gene)
fay_cer <- fay_cer[sapply(common_genes, \(x) which(x == fay_cer$gene)),]
fay_par <- fay_par[sapply(common_genes, \(x) which(x == fay_par$gene)),]
sum(fay_cer$gene == fay_par$gene)

### comparing gene counts for R1 vs R2 of each sample
# to make sure they are equivalent and can be combined

# 1) Are library sizes for R1 and R2 of each sample in each allele
#    about equal? And similar to the Fay et al. 2023 alignment counts?
sampdf <- bind_rows(tibble(sample_name = colnames(fay_cer[,-1]),
                           lib_size = colSums(fay_cer[,-1]),
                           allele = "cerAllele"),
                    tibble(sample_name = colnames(fay_par[,-1]),
                           lib_size = colSums(fay_par[,-1]),
                           allele = "parAllele"))
sampdf$read <- gsub("S[0-9]{1,2}_", "", sampdf$sample_name)
sampdf$sample_name <- gsub("_R[12]", "", sampdf$sample_name)

ggplot(pivot_wider(sampdf, id_cols = c("sample_name", "allele"),
                   names_from = "read", values_from = "lib_size"), 
       aes(x = R1, y = R2)) +
  geom_point(aes(color = allele))
# good agreement in libsize between reads. 
# Scer samples have a narrower libsize range though
# was this also seen in Fay alignment?
fay_sampdf <- read.xlsx("data_files/downloaded_from_Fay2023/Supporting_Tables.xlsx",
                        sheet = 4, startRow = 2, colNames = TRUE)
sampdf <- fay_sampdf |> 
  select(Sample, Sc_gene_reads, Sp_gene_reads) |>
  pivot_longer(cols = c("Sc_gene_reads", "Sp_gene_reads"),
               names_to = "allele", values_to = "lib_size_fay") |> 
  mutate(allele = if_else(allele == "Sc_gene_reads",
                 true = "cerAllele",
                 false = "parAllele")) |> 
  dplyr::rename("sample_name" = "Sample") |> 
  right_join(sampdf, by = c("sample_name", "allele"))
ggplot(sampdf, aes(x = lib_size_fay, y = lib_size)) +
  geom_point(aes(color = allele, shape = read)) +
  geom_abline(slope = 1, intercept = 0)
# highly correlated with fay alignment, R^2 > 0.99
lm(lib_size ~ lib_size_fay, data = sampdf) |> summary()

# 2) Are read counts for each gene in each allele about equal between
#    R1 and R2?
R1R2df <- expand_grid(sample_name = unique(sampdf$sample_name),
                      allele = c("cerAllele", "parAllele"))
R1R2df$cor <- map2(R1R2df$sample_name, R1R2df$allele,
                   \(s, a) {
                     if (a == "cerAllele") {
                       R1_reads <- fay_cer[,paste0(s, "_R1")]
                       R2_reads <- fay_cer[,paste0(s, "_R2")]
                     }
                     if (a == "parAllele") {
                       R1_reads <- fay_par[,paste0(s, "_R1")]
                       R2_reads <- fay_par[,paste0(s, "_R2")]
                     }
                     return(cor(R1_reads, R2_reads))
                   }) |> unlist()
mean(R1R2df$cor)^2
min(R1R2df$cor)^2
### based on the strong R1-R2 correlation and the similarity between
# our alignment and the Fay et al. 2023 alignment, we can
# take the mean between R1 and R2 as the read count for each allele
# checking that samples/genes are in correct order
sum(common_genes == fay_cer[,1])
nrow(fay_cer)
fay_cer_R1R2 <- fay_cer
fay_cer <- map(unique(sampdf$sample_name), \(s) {
  return((fay_cer[,paste0(s, "_R1")] + fay_cer[,paste0(s, "_R2")])/2)
}) |> purrr::reduce(.f = cbind)
colnames(fay_cer) <- unique(sampdf$sample_name)
fay_par_R1R2 <- fay_par
fay_par <- map(unique(sampdf$sample_name), \(s) {
  return((fay_par[,paste0(s, "_R1")] + fay_par[,paste0(s, "_R2")])/2)
}) |> purrr::reduce(.f = cbind)
colnames(fay_par) <- unique(sampdf$sample_name)
# testing taking the mean
random_Sample <- sample(unique(sampdf$sample_name), 1)
random_gene <- sample(c(1:nrow(fay_cer)), 1)
# cer
# what it was
fay_cer_R1R2[random_gene, paste0(random_Sample, "_R1")]
fay_cer_R1R2[random_gene, paste0(random_Sample, "_R2")]
# what it should be
(fay_cer_R1R2[random_gene, paste0(random_Sample, "_R1")] +
  fay_cer_R1R2[random_gene, paste0(random_Sample, "_R2")])/2
# what it is
fay_cer[random_gene, random_Sample]
# par
# what it was
fay_par_R1R2[random_gene, paste0(random_Sample, "_R1")]
fay_par_R1R2[random_gene, paste0(random_Sample, "_R2")]
# what it should be
(fay_par_R1R2[random_gene, paste0(random_Sample, "_R1")] +
    fay_par_R1R2[random_gene, paste0(random_Sample, "_R2")])/2
# what it is
fay_par[random_gene, random_Sample]

rm(fay_cer_R1R2, fay_par_R1R2)

# setting rownames to gene names for each count matrix
rownames(fay_cer) <- common_genes
rownames(fay_par) <- common_genes
tagseq_cer <- tagseq_cer |> select(!gene) |> as.matrix()
tagseq_par <- tagseq_par |> select(!gene) |> as.matrix()
rownames(tagseq_cer) <- common_genes
rownames(tagseq_par) <- common_genes

#### Combining allele-specific counts into counts matrix ####
# in the parents this means limiting to those reads that mapped to the
# correct parent's allele (i.e. in cer samples, cer allele counts) 
# in the hybrid this means splitting allele reads into separate columns
# (we'll check that library sizes are about even for hybrid allele pairs later on)

# tagseq
# reading in sample info
info_tagseq <- read.xlsx("data_files/downloaded_from_Krieger2020/bioSample1to999.xlsx", na.strings="not applicable", cols=c(1,4,9,13,14,15,17)) %>%
  bind_rows(read.xlsx("data_files/downloaded_from_Krieger2020/bioSample1000toEnd.xlsx", na.strings="not applicable", cols=c(1,4,9,13,14,15,17)))
colnames(info_tagseq) <- c("sample_name", "organism" , "collection_date", "genotype", "experiment","time_point", "well_flask_ID")
sum(info_tagseq$sample_name %in% colnames(tagseq_cer))
sum(info_tagseq$sample_name %in% colnames(tagseq_par))
# creating count matrix
counts_tagseq <- apply(info_tagseq, 1, \(x) {
  sample_name <- x["sample_name"]
  org <- x["organism"]
  if (!sample_name %in% colnames(tagseq_cer)) {
    cat("missing sample", sample_name, "\n")
    output <- matrix(NA, nrow = nrow(tagseq_cer), ncol = 1)
    colnames(output) <- sample_name
    return(output)
  }
  if (org == "Saccharomyces cerevisiae") {
    return(tagseq_cer[,sample_name, drop = FALSE])
  }
  if (org == "Saccharomyces paradoxus") {
    return(tagseq_par[,sample_name, drop = FALSE])
  }
  if (org == "Saccharomyces cerevisiae x Saccharomyces paradoxus") {
    cer_countcol <- tagseq_cer[,sample_name]
    par_countcol <- tagseq_par[,sample_name]
    output <- cbind(cer_countcol, par_countcol)
    colnames(output) <- c(gsub("_hyb_", "_hyc_", sample_name),
                          gsub("_hyb_", "_hyp_", sample_name))
    return(output)
  }
}) |> Reduce(f = cbind)
sum(rownames(counts_tagseq) == rownames(tagseq_cer))
sum(rownames(counts_tagseq) == rownames(tagseq_par))
# adding second row for each hybrid allele in info df
info_tagseq <- map(c(1:nrow(info_tagseq)), \(i) {
  x <- info_tagseq[i,]
  org <- info_tagseq[i,"organism"]
  if (org == "Saccharomyces cerevisiae x Saccharomyces paradoxus") {
    x_cer <- x
    x_par <- x
    x_cer["sample_name"] <- gsub("_hyb_", "_hyc_", x_cer["sample_name"])
    x_par["sample_name"] <- gsub("_hyb_", "_hyp_", x_par["sample_name"])
    output <- bind_rows(x_cer, x_par)
    return(output)
  }
  if (org != "Saccharomyces cerevisiae x Saccharomyces paradoxus") {
    return(x)
  }
}) |> purrr::reduce(.f = bind_rows)
sum(colnames(counts_tagseq) == info_tagseq$sample_name)

# fay/rnaseq
# reading in sample info
info_fay <- read.xlsx("data_files/downloaded_from_Fay2023/Supporting_Tables.xlsx",
                      sheet = 4, startRow = 2) |> 
  select(Sample, Condition, Time, Strains) |> 
  dplyr::rename("sample_name"="Sample", "experiment"="Condition",
         "time_point_num"="Time", "parents_or_hybrid"="Strains") |> 
  mutate(experiment = gsub("cold", "Cold", experiment)) |> 
  mutate(experiment = gsub("heat", "Heat", experiment)) |> 
  mutate(parents_or_hybrid = if_else(grepl("ScxSp", parents_or_hybrid),
                       true = "hybrid", false = "parents")) |> 
  filter(sample_name %in% colnames(fay_cer))

info_fay <- left_join(info_fay, 
                      bind_rows(expand_grid(parents_or_hybrid = "parents",
                                            organism = c("cer", "par")),
                                expand_grid(parents_or_hybrid = "hybrid",
                                            organism = c("hyc", "hyp"))),
                      by = "parents_or_hybrid",
                      relationship = "many-to-many")
info_fay$allele <- sapply(info_fay$organism, \(org) {
  if_else(grepl("hy[pc]", org),
          false = org,
          true = if_else(grepl("hyc", org), 
                         true = "cer", false = "par"))
})
info_fay$sample_name <- paste(info_fay$sample_name, info_fay$organism, sep = "_")
info_fay$organism <- gsub("hy[cp]", "hyb", info_fay$organism)
table(info_fay$organism, info_fay$allele)
# combining allele-specific count matrices same as tagseq
colnames(fay_cer)
counts_fay <- apply(info_fay, 1, \(x) {
  sample_name <- x["sample_name"]
  sample_name_counts <- strsplit(sample_name, "_")[[1]][1]
  al <- x["allele"]
  if (!sample_name_counts %in% colnames(fay_cer)) {
    cat("missing sample", sample_name_counts, "\n")
    output <- matrix(NA, nrow = nrow(fay_cer), ncol = 1)
    colnames(output) <- sample_name
    return(output)
  }
  if (al == "cer") {
    return(fay_cer[,sample_name_counts, drop = FALSE])
  }
  if (al == "par") {
    return(fay_par[,sample_name_counts, drop = FALSE])
  }
}) |> purrr::reduce(.f = cbind)
colnames(counts_fay) <- info_fay$sample_name

# note there are 20 paralog pairs/trios with ambiguous mapping in Yue et al. 2017 annotation:
paralog_pairs <- common_genes[(common_genes %in% grep("/", common_genes, value = TRUE))]
paralog_pairs

# at the end of this section, we have:
# 1) count matrices for each sample with its correct allele counts
# 2) info dataframes with sample_name matching columns of counts
sum(colnames(counts_fay) == info_fay$sample_name)/ncol(counts_fay)
sum(colnames(counts_tagseq) == info_tagseq$sample_name)/ncol(counts_tagseq)

#### Combining counts and sample info from fay and tagseq ####

### cleaning up some sample info column values in tagseq info prior to combining
# shorten organism names
info_tagseq$organism <- map_chr(info_tagseq$sample_name, function(s) {
  if (grepl("_cer_", s)) {
    return("cer")
  }
  if (grepl("_par_", s)) {
    return("par")
  }
  if (grepl("_hy[pc]_", s)) {
    return("hyb")
  }
})

# add allele column
info_tagseq$allele <- map_chr(info_tagseq$sample_name, function(s) {
  if (grepl("_cer_", s) | grepl("_hyc_", s)) {
    return("cer")
  }
  if (grepl("_par_", s) | grepl("_hyp_", s)) {
    return("par")
  }
})

# removing space from genotype
info_tagseq$genotype <- gsub(" ", "", info_tagseq$genotype)

# converting timepoint to integer values of minutes
timepoint_to_int <- function(t) {
  if (grepl("[0-9] h", t)) {
    return(parse_number(t)*60)
  }
  else {
    return(parse_number(t))
  }
}
info_tagseq$time_point_num <- map_dbl(info_tagseq$time_point, timepoint_to_int)
colnames(info_tagseq) <- map_chr(colnames(info_tagseq), gsub, pattern = "^time_point$", replacement = "time_point_str") # time_point_str is the version that we'll use for DESeq2, so we can set a reference level (but we'll have to do that later)

# shortening experiment value
info_tagseq$experiment <- map_chr(info_tagseq$experiment, function(e) {
  if (e == "YPD to Low N") {
    return("LowN")
  }
  if (e == "CellCycle") {
    return("CC")
  }
  if (e == "HAP4andWTonYPD") {
    return("HAP4")
  }
  if (e == "SCtoLowPi") {
    return("LowPi")
  }
})
# preserving sample information that is non-unique for replicates
info_tagseq$condition <- paste(info_tagseq$genotype,
                               info_tagseq$experiment,
                               info_tagseq$time_point_num, sep="_")

# In case you're wondering, the collection date is when RNA was collected, 
# NOT when the living yeast sample was collected, 
# so we can't use it to filter for good wells 
# (only like 65 samples total have a collection date within 24 hours 
# for all 3 samples, and that's mostly because a lot of samples were collected on those dates)
info_tagseq <- select(info_tagseq, !collection_date)

### cleaning up some sample info column values in fay info prior to combining
setdiff(colnames(info_tagseq), colnames(info_fay))
info_fay$genotype <- "WT"
info_fay$time_point_str <- info_fay$time_point_num |> as.character()
info_fay$condition <- paste(info_fay$genotype,
                            info_fay$experiment,
                            info_fay$time_point_num, sep="_")
table(info_fay$condition, info_fay$organism)
info_fay <- arrange(info_fay, organism, experiment, time_point_num)
info_fay$well_flask_ID <- c("rep1", "rep2")
counts_fay <- counts_fay[,info_fay$sample_name] # also rearranging counts, because we index their columns by sample name
table(info_fay$condition, info_fay$well_flask_ID) # should have 4 entries for each rep: cer/par/hyc/hyp
setdiff(colnames(info_fay), colnames(info_tagseq))
info_fay <- select(info_fay, !parents_or_hybrid)

# last checks
setequal(colnames(info_tagseq), colnames(info_fay))
sum(colnames(counts_fay) == info_fay$sample_name)/ncol(counts_fay)
sum(colnames(counts_tagseq) == info_tagseq$sample_name)/ncol(counts_tagseq)
# combining
sample_info <- bind_rows(info_tagseq, info_fay)
counts <- cbind(counts_tagseq, counts_fay)
cat("percent sample info columns in same order as count data (before matching): ",
    sum(colnames(counts)==sample_info$sample_name)/ncol(counts))

#### Renaming replicates in LowN ####
sample_info |> group_by(organism, time_point_str, well_flask_ID,
                        experiment, genotype) |> 
  summarise(nreps = n()) |> filter((organism != "hyb" & nreps == 1) |
                                     organism == "hyb" & nreps == 2) |> nrow()
nrow(sample_info) - sum(sample_info$organism == "hyb")/2 # should be the same

# Every LowN should have at most 3 entries by the end:
sample_info |> filter(experiment == "LowN") |> group_by(organism,
                                                        well_flask_ID,
                                                        genotype) |> 
  summarise(nreps = n()) |> filter(nreps == 3)
sample_info |> filter(experiment == "LowN") |> group_by(organism,
                                                        well_flask_ID,
                                                        genotype) |> 
  summarise(nreps = n()) |> filter(nreps != 3)

# First of all, the LowN GS2018 samples have a different well_flask_ID# for the 1 hr sample than for 0 or 16 hr...
# For each row in the well plate (A-H), these are the numbers that are paired: 1-7, 2-8, 3-9, 4-10, 5-11, 6-12
col1 <- sapply(c("A", "B", "C", "D", "E", "F", "G", "H"), function(x) return(paste0(x, c(1:6)))) %>% as.vector()
col2 <- sapply(c("A", "B", "C", "D", "E", "F", "G", "H"), function(x) return(paste0(x, c(7:12))))  %>% as.vector()
gs2018_lookup <- tibble(TP1_TP3 = col1, TP2 = col2)

# arbitrarily assigning each ID its TP1/TP3 (column 1 in lookup table) value
standardizeGS2018ID <- function(id) {
  id_clipped <- gsub("_WT2_GS2018", "", id)
  id_clipped <- gsub("_G12_GS2018", "", id_clipped)
  id_clipped <- gsub("_G6_GS2018", "", id_clipped)
  rownum <- c(which(gs2018_lookup$TP1_TP3 == id_clipped), which(gs2018_lookup$TP2 == id_clipped))
  new_id <- gsub(gs2018_lookup$TP2[rownum], gs2018_lookup$TP1_TP3[rownum], id)
  return(new_id)
}
# tests for standardizeGS2018ID
standardizeGS2018ID("G8_G12_GS2018")

# before applying
sample_info |> filter(organism %in% c("cer", "par") & grepl("GS2018", well_flask_ID)) |> 
  select(well_flask_ID) |> table() 
sample_info |> filter(organism == "hyb" & grepl("GS2018", well_flask_ID)) |> 
  select(well_flask_ID) |> table() 

# applying to GS2018 samples
GS2018_idxs <- grepl("GS2018", sample_info$well_flask_ID)
sample_info$well_flask_ID[GS2018_idxs] <- sapply(sample_info$well_flask_ID[GS2018_idxs],
                                                 standardizeGS2018ID)

# after applying
sample_info |> filter(organism %in% c("cer", "par") & grepl("GS2018", well_flask_ID)) |> 
  select(well_flask_ID) |> table() 
sample_info |> filter(organism == "hyb" & grepl("GS2018", well_flask_ID)) |> 
  select(well_flask_ID) |> table() 
# hybrids should have 6, 2 alleles x 3 timepoints, Parents have just 3 timepoints. Except for H1, which is missing a timepoint
# giving the GS2018s their proper well IDs in sample_info

# getting rid of the "GS" part of the rep name, which can be different for the same
# sample at different timepoints and replacing it with the tag immediately before
# the well_flask_ID in the sample name for non-unique well_flask_IDs
sample_info$well_flask_ID <- gsub("_GS.*", "", sample_info$well_flask_ID)
sample_info$new_id <- map2(sample_info$well_flask_ID, 
             sample_info$sample_name, \(i, s) {
               ex <- sample_info |> filter(sample_name == s) |> 
                 select(experiment) |> pull()
               if (ex == "LowN") {
                 org <- sample_info |> filter(sample_name == s) |> 
                   select(organism) |> pull()
                 is_duplicate <- sample_info |> 
                   filter(experiment == "LowN" & 
                            organism == org &
                            well_flask_ID == i) |> 
                   group_by(time_point_str) |> 
                   summarise(n_per_tp = n())
                 if (any(is_duplicate$n_per_tp > 1)) {
                   s <- gsub("_GS.*", "", s)
                   new_tag <- gsub(i, "", s) |> strsplit(split = "_") |> 
                     unlist() |> tail(n = 1)
                   return(paste(new_tag, i, sep = "_"))
                 }
                 else {
                   return(i)
                 }
               }
               else {
                 return(i)
               }
             }) |> unlist()
# checking for non-unique IDs
sample_info |> 
  group_by(new_id, organism, experiment, time_point_str, genotype) |> 
  summarise(n_per_condition = n()) |> 
  filter((organism != "hyb" & n_per_condition > 1) |
           (organism == "hyb" & n_per_condition > 2)) # should be empty

# updating well_flask_ID 
sample_info$well_flask_ID <- sample_info$new_id
sample_info <- select(sample_info, -"new_id")

# checking example (two C5_A10s for the 960 timepoint)
sample_info |> filter(experiment == "LowN" & organism == "par" &
                        genotype == "GCN4delete")
# should have additional tag P1 or P2 on the C5_A10s

#### misc filtering of additional samples and genes ####

# Whittling LowPi to just WT --- not enough of them to warrant the complication of including a second genotype (plus they're only present in -5 and 180 min samples)
sample_info %>% filter(experiment == "LowPi") %>% select(genotype) %>% table()
# species
keep <- !(sample_info$experiment == "LowPi" & sample_info$genotype == "PHO4delete")
sample_info <- sample_info[keep,]
counts <- counts[,(colnames(counts) %in% sample_info$sample_name)]

# our 46 TF deletions
TFdel_lookup <- read_delim("data_files/downloaded_genomes_and_features/yeastract_46TFs.csv", col_names = FALSE, col_select = c(1,2), delim = ";") # gets some warnings, but so far has been fine
colnames(TFdel_lookup) <- c("common", "systematic")

# Check for missing (NA) values
geneHasNAs <- apply(counts, 1, function(x) {
  isNA <- sapply(x, is.na)
  return(any(isNA))
}) 
sum(geneHasNAs) # should have 0

#### Removing samples with small library sizes ####
# Exploring library sizes (un-normalized)
libsizes <- colSums(counts) # these counts include all samples, hybrid and parental
ggplot(tibble(libsize = libsizes), aes(x = libsize)) + geom_histogram()
min(libsizes)
max(libsizes)
median(libsizes)
sum(libsizes < 100000)
sort(libsizes)[c(1:20)] # a few of the CC have fairly tiny libraries, and will need to be removed before normalizing (or else you get ridiculously high outlier counts from small counts --- ex) (5/7000)*1e6 = 714 whereas (5/200000)*1e6 = 25)
# are the hybrid reps library sizes correlated between hyc and hyp alleles?
plotdf <- tibble(sample_name = names(libsizes[grepl("_hy[pc]", names(libsizes))]),
                 libsize = libsizes[grepl("_hy[pc]", names(libsizes))]) |> 
  left_join(y = select(filter(sample_info, organism == "hyb"), sample_name, experiment), by = "sample_name")
plotdf$allele <- if_else(grepl(pattern = "_hyc_", plotdf$sample_name),
                         true = "cer", false = "par")
plotdf$sample_name <- gsub("_hy[pc]_", "_hyb_", plotdf$sample_name)
plotdf <- plotdf |> pivot_wider(id_cols = c("sample_name", "experiment"), 
                                names_from = allele,
                                values_from = libsize)
ggplot(plotdf, aes(x = cer, y = par)) + 
  geom_point(aes(color = experiment), alpha = 0.5) +
  geom_text(data = filter(plotdf, cer > par*1.5 | par > cer*1.5),
            aes(label = sample_name), check_overlap = TRUE, color = "green") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_rect(xmin = 0, xmax = 100000, ymin = 0, ymax = 100000, color = "blue", alpha = 0) +
  geom_vline(xintercept = 100000, color = "blue", alpha = 0.5) +
  geom_hline(yintercept = 100000, color = "blue", alpha = 0.5)
# based on this plot, removing cell cycle 10 for having libraries too biased to one allele
# (the samples within the blue box will be removed by our library size threshold below)
# and no other samples should be below the blue lines but not within the blue box
hyblib_biased_samples <- c("WT_10_hyc_CellCycle_rep1",
                           "WT_10_hyp_CellCycle_rep1")
colSums(counts[,hyblib_biased_samples])
# removing
sample_info <- filter(sample_info, !(sample_name %in% hyblib_biased_samples))
counts <- counts[,colnames(counts) %in% sample_info$sample_name]

# generating log2-cpm counts
avgLibSizeInMillions <- mean(libsizes)/1e6
counts_lcpm <- log2((counts/libsizes)*1e6 + 2/avgLibSizeInMillions) # equivalent to edgeR::cpm with log=TRUE

# checking if paradoxus genes are *all* expressed less than cerevisiae pre-normalization
# (already checked this more elegantly in hybrid, where samples can be paired)
cer_par_mean_expr_diffs_parents <- rowMeans(counts[,sample_info$allele == "cer"]) - rowMeans(counts[,sample_info$allele == "par"])
hist(sign(cer_par_mean_expr_diffs_parents)*log(abs(cer_par_mean_expr_diffs_parents)), breaks = 50) 

# Does one species tend to have larger libsizes?
plotdf <- tibble(libsize_cer = colSums(counts[, grepl("_cer", colnames(counts))])[sample(c(1:400), 400)],
                 libsize_par = colSums(counts[, grepl("_par", colnames(counts))])[sample(c(1:400), 400)],
                 libsize_hyc = colSums(counts[, grepl("_hyc", colnames(counts))])[sample(c(1:400), 400)],
                 libsize_hyp = colSums(counts[, grepl("_hyp", colnames(counts))])[sample(c(1:400), 400)]) %>%
  pivot_longer(cols = c(libsize_cer, libsize_par, libsize_hyc, libsize_hyp))
ggplot(filter(plotdf, name %in% c("libsize_cer", "libsize_par")), aes(x = name, y = value)) + geom_boxplot(aes(fill = name)) 
t.test(value ~ name, filter(plotdf, name %in% c("libsize_cer", "libsize_par")))
ggplot(filter(plotdf, name %in% c("libsize_hyc", "libsize_hyp")), aes(x = name, y = value)) + geom_boxplot(aes(fill = name))
t.test(value ~ name, filter(plotdf, name %in% c("libsize_hyc", "libsize_hyp"))) 
# par libraries are smaller than cer, difference isn't significant between hybrid alleles,
# interestingly this is the opposite of what is seen in the Heat/Cold data alone:
plotdf <- tibble(libsize_cer = colSums(counts[, grepl("_cer", colnames(counts)) & sample_info$experiment %in% c("Heat", "Cold")]),
                 libsize_par = colSums(counts[, grepl("_par", colnames(counts)) & sample_info$experiment %in% c("Heat", "Cold")]),
                 libsize_hyc = colSums(counts[, grepl("_hyc", colnames(counts)) & sample_info$experiment %in% c("Heat", "Cold")]),
                 libsize_hyp = colSums(counts[, grepl("_hyp", colnames(counts)) & sample_info$experiment %in% c("Heat", "Cold")])) %>%
  pivot_longer(cols = c(libsize_cer, libsize_par, libsize_hyc, libsize_hyp))
ggplot(filter(plotdf, name %in% c("libsize_cer", "libsize_par")), aes(x = name, y = value)) + geom_boxplot(aes(fill = name)) 
t.test(value ~ name, filter(plotdf, name %in% c("libsize_cer", "libsize_par")))
ggplot(filter(plotdf, name %in% c("libsize_hyc", "libsize_hyp")), aes(x = name, y = value)) + geom_boxplot(aes(fill = name))
t.test(value ~ name, filter(plotdf, name %in% c("libsize_hyc", "libsize_hyp"))) 

#### normalizing ####
# normalizing counts to adjust for differences in library size
# @input: count matrix (genes are rows, columns are samples)
# @output: a count matrix normalzied for library size---integer counts in counts-per-million
countsPerMillion <- function(.cts) {
  librarySizes <- colSums(.cts, na.rm = TRUE)
  output <- apply(.cts, 1, function(x) {
    normalized <- (x/librarySizes)*1e6
    return(round(normalized))
  })
  return(t(output)) # For some unhinged reason, a vector output of apply across ROWS forms the COLUMNS of a new matrix
}
# tests for countsPerMillion
test_cts <- counts[,grepl("_par", colnames(counts))]
test_rowIdx <- sample(c(1:nrow(test_cts)), 1)
test_colIdx <- sample(c(1:ncol(test_cts)), 1)
test_count <- test_cts[test_rowIdx, test_colIdx]
test_cpm <- countsPerMillion(test_cts)
test_cpm_cpm <- countsPerMillion(test_cpm)
((test_count/colSums(test_cts, na.rm = TRUE)[test_colIdx])*1e6) %>% round() # what it should be
test_cts[test_rowIdx, test_colIdx] # what it is before using our function
test_cpm[test_rowIdx, test_colIdx] # what it is using our function
test_cpm_cpm[test_rowIdx, test_colIdx] # what it is if you run cpm too many times (should be the same as test_cpm)

### eliminating samples with too-small libraries
# Example of how small libraries skew data, YIL134W in CC before normalizing:
test_counts <- counts[, (grepl("_cer", colnames(counts)) | 
                           grepl("_par", colnames(counts))) &
                        grepl("_CellCycle_", colnames(counts))]
libSizes <- colSums(test_counts)
# before normalizing:
genedf <- tibble(expr = as.numeric(test_counts["YIL134W",]),
                 libsize = libSizes,
                 sample_name = colnames(test_counts))
ggplot(genedf, aes(x = libSizes, y = expr)) + geom_point(aes(color = sample_name == "WT_34_par_CellCycle_rep2")) + theme(legend.position = "none")
# after normalizing:
test_counts_cpm <- countsPerMillion(test_counts)
genedf <- tibble(expr = as.numeric(test_counts_cpm["YIL134W",]),
                 libsize = libSizes,
                 sample_name = colnames(test_counts))
ggplot(genedf, aes(x = libSizes, y = expr)) + geom_point(aes(color = sample_name == "WT_34_par_CellCycle_rep2")) + theme(legend.position = "none")
# At a certain lib size, there appears to begin to be a correlation between count and libsize, what size is that?
libSizes <- colSums(counts)
gene_idx <- sample(rownames(counts), 1) # rerun to make sure most genes have expr and libsize correlated past cutoff of 200k libsize
plotdf <- tibble(libsize = libSizes, 
                 expr = as.numeric(counts[gene_idx,]), 
                 gene_name = gene_idx,
                 sample_name = colnames(counts)) |> 
  left_join(select(sample_info, sample_name, experiment),
            by = "sample_name")
ggplot(plotdf, aes(x = libsize, y = expr)) + 
  geom_point(aes(color = experiment)) + geom_vline(xintercept = 100000, col = "red") + theme_classic() + ggtitle(gene_idx)

# Eliminating samples with libsize < 100,000
keep <- colSums(counts) > 100000
sum(keep)/length(keep)
sum(!keep)
keep_samples <- colnames(counts)[keep]
sample_info <- sample_info |> filter(sample_name %in% keep_samples)
counts <- counts[, sample_info$sample_name]

### Normalizing by length
# Normalizing Fay et al. 2023 data to cpm taking account of gene length
# following this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7373998/
# This is paried-end (one count = one pair), stranded, poly-A selected,
# so we should be safe to compare these counts to our tag-seq cpm counts
# (although if they weren't polyA selected or weren't stranded, direct comparison 
# isn't advisable b/c of biases in which subsets of the transcriptome are counted)

# normalizing with this strategy: 
# for each gene i, first obtain gene_i = 10^6 * (raw counts for gene_i/length(gene_i))
# then divide each gene by the new "library size" (the sum of these ^ across all genes in the sample)

# aka CPM = 10^6 * [(raw counts for gene_i/length(gene_i))/sum_i=1^i=nGenes(raw counts for gene_i/length(gene_i))]

# gene lengths from annotation file:
# Scer
gene_lens_scer <- read_tsv("data_files/downloaded_genomes_and_features/S288C.all_feature.gff",
                      col_names = FALSE, col_select = c(3,4,5,9)) |> 
  dplyr::rename("feature"="X3",
         "start"="X4",
         "end"="X5",
         "gene_name"="X9") |> 
  filter(feature == "gene") |> 
  as_tibble()
gene_lens_scer$gene_name <- sapply(gene_lens_scer$gene_name, \(g) {
  out_g <- gsub("Name=", "", strsplit(g, split = ";")[[1]][2])
  return(out_g)
}) |> as.character()
sum(gene_lens_scer$end > gene_lens_scer$start)
sum(gene_lens_scer$end < gene_lens_scer$start)
gene_lens_scer$length <- gene_lens_scer$end - gene_lens_scer$start
# Spar
gene_lens_spar <- read_tsv("data_files/downloaded_genomes_and_features/CBS432.all_feature.gff",
                           col_names = FALSE, col_select = c(3,4,5,9)) |> 
  dplyr::rename("feature"="X3",
         "start"="X4",
         "end"="X5",
         "gene_name"="X9") |> 
  filter(feature == "gene") |> 
  as_tibble()
gene_lens_spar$gene_name <- sapply(gene_lens_spar$gene_name, \(g) {
  out_g <- gsub("Name=", "", strsplit(g, split = ";")[[1]][2])
  return(out_g) |> as.character()
})
sum(gene_lens_spar$end > gene_lens_spar$start)
sum(gene_lens_spar$end < gene_lens_spar$start)
gene_lens_spar$length <- gene_lens_spar$end - gene_lens_spar$start
gene_lens_scer <- gene_lens_scer |> select(gene_name, length) |> 
  filter(gene_name %in% common_genes) |> unique()
gene_lens_spar <- gene_lens_spar |> select(gene_name, length) |> 
  filter(gene_name %in% common_genes) |> unique()
# check that duplicated genes have about the same length before removing duplicates
gene_lens_scer[gene_lens_scer$gene_name %in% gene_lens_scer$gene_name[duplicated(gene_lens_scer$gene_name)],] |> 
  arrange(gene_name)
gene_lens_scer <- gene_lens_scer[!duplicated(gene_lens_scer$gene_name),]
gene_lens_spar[gene_lens_spar$gene_name %in% gene_lens_spar$gene_name[duplicated(gene_lens_spar$gene_name)],] |> 
  arrange(gene_name)
# the ones in Scer are fine, but a few in Spar are concerning:
two_lengths_spar <- c("YAR028W", "YBR140C", 
                      "YDR025W/YBR048W", "YPL091W")
# what do they look like in Scer?
gene_lens_scer |> filter(gene_name %in% two_lengths_spar)
# looks like the larger length is most accurate
gene_lens_spar <- gene_lens_spar |> arrange(gene_name, desc(length))
gene_lens_spar <- gene_lens_spar[!duplicated(gene_lens_spar$gene_name),]
gene_lens_spar |> filter(gene_name %in% two_lengths_spar)
# how many genes have significantly different lengths in both species?
plotdf <- left_join(x = dplyr::rename(select(gene_lens_scer, gene_name, length), "Scer"="length"),
                    y = dplyr::rename(select(gene_lens_spar, gene_name, length), "Spar"="length"),
                    by = "gene_name")
ggplot(plotdf, aes(x = Scer, y = Spar)) + geom_point()
# most are the same exact length, a few have slightly different lengths
# arranging gene_lens into same gene order as count matrices
sum(gene_lens_scer$gene_name == common_genes) # pre-arranging
gene_lens_scer <- gene_lens_scer[sapply(common_genes, \(g) which(gene_lens_scer$gene_name == g)),]
sum(gene_lens_scer$gene_name == common_genes) # after
sum(gene_lens_spar$gene_name == common_genes) # pre-arranging
gene_lens_spar <- gene_lens_spar[sapply(common_genes, \(g) which(gene_lens_spar$gene_name == g)),]
sum(gene_lens_spar$gene_name == common_genes) # after

# @input: counts matrix, vector of lengths in same order as rows in .cts
countsPerMillionWithLength <- function(.cts, .lens) {
  rnames <- rownames(.cts)
  cnames <- colnames(.cts)
  if (length(.lens) != nrow(.cts)) {
    stop("gene lengths are not same length as nrow counts", length(.lens),
         "vs", nrow(.cts), "\n")
  }
  genes_over_length <- map(rownames(.cts), \(g) {
    gene_vec <- .cts[g,] |> as.numeric()
    gene_len <- .lens[which(rownames(.cts) == g)]
    return(gene_vec/gene_len)
  }) |> purrr::reduce(.f = rbind)
  lib_sums <- colSums(genes_over_length)
  output <- apply(genes_over_length, 1, \(x) {return(x/lib_sums)}) |> t()
  rownames(output) <- rnames
  colnames(output) <- cnames
  return(round(output*10^6))
}
# tests for countsPerMillionWithLength
test_cts <- counts_fay
test_rowIdx <- sample(c(1:nrow(test_cts)), 1)
test_colIdx <- sample(c(1:ncol(test_cts)), 1)
test_count <- test_cts[test_rowIdx, test_colIdx]
test_lens <- gene_lens_scer$length |> as.numeric()
sum(gene_lens_scer$gene_name == rownames(test_cts))/nrow(test_cts) # should be 100%
test_cpm <- countsPerMillionWithLength(test_cts, test_lens)
# test_cpm_cpm <- countsPerMillionWithLength(test_cpm, test_lens)
round(((test_count/test_lens[test_rowIdx])/sum(test_cts[,test_colIdx]/test_lens))*10^6) # what it should be
test_cts[test_rowIdx, test_colIdx] # what it is before using our function
test_cpm[test_rowIdx, test_colIdx] # what it is using our function
# test_cpm_cpm[test_rowIdx, test_colIdx] # what it is if you run cpm too many times (unlike the TagSeq normalization, this cannot be run over and over again)

counts_unnorm <- counts

# normalizing Fay et al. 2023 by length in each species
fay_cer <- counts[,sample_info$experiment %in% c("Heat", "Cold") &
                    sample_info$allele == "cer"]
sum(gene_lens_scer$gene_name == rownames(fay_cer))/nrow(fay_cer)
fay_cer_cpm <- countsPerMillionWithLength(.cts = fay_cer,
                                      .lens = gene_lens_scer$length)
fay_par <- counts[,sample_info$experiment %in% c("Heat", "Cold") &
                    sample_info$allele == "par"]
sum(gene_lens_spar$gene_name == rownames(fay_par))/nrow(fay_par)
fay_par_cpm <- countsPerMillionWithLength(.cts = fay_par,
                                      .lens = gene_lens_spar$length)
# normalizing tagseq not by length
tagseq_cpm <- countsPerMillion(counts[,!sample_info$experiment %in% c("Heat", "Cold")])
# re-combining
counts <- cbind(tagseq_cpm, fay_cer_cpm, fay_par_cpm)
counts <- counts[,sample_info$sample_name]

sum(rownames(counts_unnorm) == rownames(counts))/nrow(counts)
sum(colnames(counts_unnorm) == colnames(counts))/ncol(counts)

#### splitting off hybrid counts ####
counts_allele <- counts[, grepl("_hy[pc]", colnames(counts))]
counts <- counts[, !grepl("_hy[pc]", colnames(counts))]
sample_info_allele <- sample_info |> filter(organism == "hyb")
sample_info <- sample_info |> filter(organism != "hyb")
counts_unnorm_allele <- counts_unnorm[, grepl("_hy[pc]", colnames(counts_unnorm))]
counts_unnorm <- counts_unnorm[, !grepl("_hy[pc]", colnames(counts_unnorm))]

dim(counts)
dim(counts_unnorm)
dim(sample_info)
dim(counts_allele)
dim(counts_unnorm_allele)
dim(sample_info_allele)
sample_info |> select(sample_name, organism, allele) |> slice_sample(n = 10)
sample_info_allele |> select(sample_name, organism, allele) |> slice_sample(n = 10)
# checking that every hybrid sample name has exactly one cer and one par row
sample_info_allele |> select(sample_name, organism, allele) |> 
  group_by(sample_name) |> 
  summarise(ncernpar = unique(table(allele))) |> 
  select(ncernpar) |> table()
sample_info_allele |> select(allele) |> table()

#### visualizing expression ####

# one random gene
random_gene <- sample(rownames(counts)[which(rowSums(counts, na.rm = TRUE) > 100000)], 1)
oneGeneBoxplots <- function(.gene_idx) {
  expr_cer <- counts[.gene_idx, which(sample_info$organism == "cer")]
  expr_par <- counts[.gene_idx, which(sample_info$organism == "par")]
  expr_hyc <- counts_allele[.gene_idx, which(sample_info_allele$allele == "cer")]
  expr_hyp <- counts_allele[.gene_idx, which(sample_info_allele$allele == "par")]
  plotdf <- tibble(expr = c(expr_cer, expr_par, expr_hyc, expr_hyp),
                   allele = c(rep("cer", length(expr_cer)),
                              rep("par", length(expr_par)),
                              rep("hyc", length(expr_hyc)),
                              rep("hyp", length(expr_hyp))))
  p <- ggplot(plotdf, aes(x = allele, y = log2(expr + 1))) + geom_boxplot()
  return(p)
}
oneGeneBoxplots(random_gene)

# all genes
mean_expr_cer <- rowMeans(counts[, which(sample_info$organism == "cer")], na.rm = TRUE)
mean_expr_par <- rowMeans(counts[, which(sample_info$organism == "par")], na.rm = TRUE)
mean_expr_hyc <- rowMeans(counts_allele[, which(sample_info_allele$allele == "cer")], na.rm = TRUE)
mean_expr_hyp <- rowMeans(counts_allele[, which(sample_info_allele$allele == "par")], na.rm = TRUE)
plotdf <- tibble(cer_expr = c(mean_expr_cer, mean_expr_hyc),
                 par_expr = c(mean_expr_par, mean_expr_hyp),
                 type = c(rep("parent", length(rownames(counts))), rep("hybrid", length(rownames(counts)))))
ggplot(plotdf, aes(x = log(cer_expr), y = log(par_expr))) + geom_point(aes(color = type))
# Mean expression of most genes is highly correlated between species and between hybrid alleles

#### visualizing lowly expressed genes ####
plotdf <- tibble(meanExpr = apply(counts_lcpm, 1, mean),
                 sdExpr = apply(counts_lcpm, 1, sd))

# out of curiousity, here's the mean expression for TFs in the 
# genotypes where they are deleted (in cpm)
meansdel <- map(c(1:nrow(TFdel_lookup)), function(i) {
  genedel <- TFdel_lookup$common[i]
  gene_idx <- TFdel_lookup$systematic[i]
  delcounts <- counts_lcpm[rownames(counts_lcpm) == gene_idx, 
                           grepl(genedel, colnames(counts_lcpm))] |> as.numeric()
  return(mean(delcounts))
}) %>% unlist()
mean(meansdel, na.rm = TRUE)

# choosing cutoff expression (in not-log scale)
ggplot(data = plotdf, aes(x = meanExpr, y = sdExpr)) + geom_hex(bins = 70) + 
  geom_vline(xintercept = 5, color = "red")
# a good cutoff is where the lower bound of the sd stops being related to mean (i.e. becomes horizontal)
# lowly expressed genes will systematically have higher variance (or extremely low variance for that tiny tail of genes with high 0 counts)
2^5
cutoffExpr <- 30

# Note: we only want to filter genes that are lowly expressed in cer, par, hyc and hyp
# Example of why this matters: YPR199C
oneGeneBoxplots("YPR199C") # strongly expressed in cer and hyc
keep_overall <- apply(counts_lcpm, 1, mean) > 5
"YPR199C" %in% rownames(counts)[keep_overall] # wouldn't be kept by overall cutoff threshold

#### hybrid par allele deletion blocks ####

# these two blocks of contiguous genes are deleted in the F1 hybrid paradoxus allele in the CellCycle experiment,
# (will be removed from CellCycle specifically in data_for_figure_scripts.R)
omit_list <- c("YLR078C", "YLR077W", "YLR074C", "YLR072W", "YLR075W", "YLR073C", # large hybrid CC paradoxus haplotype deletion
               "YNL247W", "YNL244C")
# illustrating with one of these genes:
gene_idx <- sample(omit_list, 1)
# parents
genedf <- tibble(expr = as.numeric(counts[gene_idx,])) |> 
  bind_cols(sample_info) |> pivot_wider(id_cols = c("condition", "experiment"),
                                        names_from = "allele", values_from = "expr",
                                        values_fn = mean)
ggplot(genedf, aes(x = cer, y = par)) + geom_point(aes(color = experiment))
# hybrid
genedf <- tibble(expr = as.numeric(counts_allele[gene_idx,])) |> 
  bind_cols(sample_info_allele) |> pivot_wider(id_cols = c("condition", "experiment"),
                                               names_from = "allele", values_from = "expr",
                                               values_fn = mean)
ggplot(genedf, aes(x = cer, y = par)) + geom_point(aes(color = experiment))

#### Heat/Cold QC ####
# Comparing our in-house alignment to the Fay et al. 2023 alignment with a PCA
fay_inHouse <- cbind(counts[,sample_info$experiment %in% c("Heat", "Cold")],
                     counts_allele[,sample_info_allele$experiment %in% c("Heat", "Cold")])
load("data_files/Cleaned_Fay_Counts.RData")
load("data_files/Cleaned_Fay_Counts_Allele.RData")
fay_2023 <- cbind(fay, fay_allele)
colnames(fay_inHouse) <- colnames(fay_inHouse) |> 
  map(.f = \(nm) {
    nmnum <- parse_number(nm)
    if_else(grepl(pattern = "_cer", nm) | grepl(pattern = "_hyc", nm),
            true = paste0("Sc", nmnum),
            false = paste0("Sp", nmnum))
  }) |> unlist()
common_cols_fay_pca <- intersect(colnames(fay_2023), colnames(fay_inHouse))
fay_2023 <- fay_2023[, common_cols_fay_pca]
fay_inHouse <- fay_inHouse[, common_cols_fay_pca]
colnames(fay_2023) <- paste0(colnames(fay_2023), "_2023")
colnames(fay_inHouse) <- paste0(colnames(fay_inHouse), "_inHouse")
common_genes_fay_pca <- intersect(rownames(fay_2023),
                                  rownames(fay_inHouse))
# sample pca
pcamat <- cbind(fay_2023[common_genes_fay_pca,], 
                fay_inHouse[common_genes_fay_pca,])
# pcamat <- pcamat[rowMeans(pcamat) > 30 & rownames(pcamat) != "YFL014W",]
pcamat <- pcamat[rowMeans(pcamat) > 30,]
covmat <- cov(pcamat)
colnames(covmat) <- colnames(pcamat)
pca_res <- prcomp(covmat)
sample_info_pca <- sample_info |> filter(experiment %in% c("Heat", "Cold")) |> 
  mutate(sample_name = map(sample_name, .f = \(nm) {
    nmnum <- parse_number(nm)
    if_else((grepl(pattern = "_cer", nm) | grepl(pattern = "_hyc", nm)),
            true = paste0("Sc", nmnum),
            false = paste0("Sp", nmnum))
  }) |> unlist())
sample_info_pca <- sample_info_allele |> filter(experiment %in% c("Heat", "Cold")) |> 
  mutate(sample_name = map(sample_name, .f = \(nm) {
    nmnum <- parse_number(nm)
    if_else((grepl(pattern = "_cer", nm) | grepl(pattern = "_hyc", nm)),
            true = paste0("Sc", nmnum),
            false = paste0("Sp", nmnum))
  }) |> unlist()) |> bind_rows(sample_info_pca)
pcadf <- tibble(pc1 = pca_res$x[,1], 
                pc2 = pca_res$x[,2],
                sample_name = colnames(covmat)) |> 
  mutate(sample_name_nonunique = gsub(pattern = "_2023", 
                                      replacement = "", 
                                      gsub(pattern = "_inHouse", 
                                           replacement = "",
                                           sample_name))) |> 
  left_join(y = sample_info_pca, by = c("sample_name_nonunique"="sample_name"))
var_pct <- summary(pca_res)$importance[2, 1:2] # % variance explained
pcadf$sample_num <- parse_number(pcadf$sample_name_nonunique)
pcadf$alignment <- if_else(grepl(pattern = "_2023", pcadf$sample_name),
                           true = "2023", false = "inHouse")
pcadf$tag <- paste0(substring(pcadf$experiment, 1, 1),
                    pcadf$time_point_num)
# all 3
ggplot(pcadf,
       aes(x = pc1, y = pc2)) + 
  geom_line(aes(group = sample_name_nonunique,
                color = organism)) +
  geom_text(aes(label = tag,
                color = alignment)) +
  xlab(paste0("PC1, ", round(var_pct[1]*100, digits = 0), 
              "% of variance")) + 
  ylab(paste0("PC2, ", round(var_pct[2]*100, digits = 0), 
              "% of variance")) +
  theme(legend.title = element_blank())
# Those two Scer Heat 30min sample outliers are due to one gene, YFL014W (HSP12),
# (run the PCA without this gene and the outlier samples will no longer be outliers),
# which had ~100000 reads in our in-house alignment but an order of magnitude
# fewer reads in Fay et al. 2023. It is a highly expressed gene in heat shock,
# good to keep in mind this outlier

# Avg gene expression in YPD 0min samples of RNAseq vs tagseq
test_allele <- "cer"
plotdf <- tibble(mean0_LowN = rowMeans(counts[,sample_info$experiment == "LowN" &
                                                sample_info$time_point_num == 0 &
                                                sample_info$allele == test_allele]),
                 mean0_CC = rowMeans(counts[,sample_info$experiment == "CC" &
                                              sample_info$time_point_num == 0 &
                                              sample_info$allele == test_allele]),
                 mean0_Heat = rowMeans(counts[,sample_info$experiment == "Heat" &
                                                sample_info$time_point_num == 0 &
                                                sample_info$allele == test_allele]),
                 mean0_Cold = rowMeans(counts[,sample_info$experiment == "Cold" &
                                                sample_info$time_point_num == 0 &
                                                sample_info$allele == test_allele]),
                 var_LowN = rowVars(counts[,sample_info$experiment == "LowN" &
                                             sample_info$allele == test_allele]),
                 var_CC = rowVars(counts[,sample_info$experiment == "CC" &
                                           sample_info$allele == test_allele]),
                 var_Heat = rowVars(counts[,sample_info$experiment == "Heat" &
                                             sample_info$allele == test_allele]),
                 var_Cold = rowVars(counts[,sample_info$experiment == "Cold" &
                                             sample_info$allele == test_allele]))

# 1) each gene's variance across timepoints is less variable than in tag-seq
# mean at 0 min vs var across timepoints for Heat/Cold/LowN/CC
ggplot(plotdf, aes(x = log2(mean0_LowN), y = log2(var_LowN))) + geom_point() # LowN
ggplot(plotdf, aes(x = log2(mean0_CC), y = log2(var_CC))) + geom_point() # CC
ggplot(plotdf, aes(x = log2(mean0_Heat), y = log2(var_Heat))) + geom_point() # Heat
ggplot(plotdf, aes(x = log2(mean0_Cold), y = log2(var_Cold))) + geom_point() # Cold
# Conclusion: does not appear to be less variable

# 2) less variance in mean expression for all genes than in tag-seq
# each gene's mean expr between replicates at 0 min
# Heat/Cold vs other experiments
# LowN vs CC first as a control
ggplot(plotdf, aes(x = log2(mean0_LowN), y = log2(mean0_CC))) + geom_point() # LowN vs CC
ggplot(plotdf, aes(x = log2(mean0_Heat), y = log2(mean0_Cold))) + geom_point() # Heat vs Cold
ggplot(plotdf, aes(x = log2(mean0_LowN), y = log2(mean0_Heat))) + geom_point() # LowN vs Heat
ggplot(plotdf, aes(x = log2(mean0_LowN), y = log2(mean0_Cold))) + geom_point() # LowN vs Cold
# interestingly Spar has more similar expr between experiments
# (cerevisiae is the T73 wine strain that has the HGT from Torulaspora microellipsoides,
# paradoxus is the N17 that Fay et al. 2023 aligned to the CBS432 genome)

#### Splitting off TFdel datasets ####
# probably archive b/c we don't do DESeq2, but we might
# so keeping until we get to constructing that TFdel figure

# splitting off TFdel dataset b/c unnormalized counts with replicates are needed for DESeq2
infos_TFdel <- list(cer = sample_info[sample_info$experiment == "LowN" &
                                        sample_info$organism == "cer",],
                    par = sample_info[sample_info$experiment == "LowN" &
                                        sample_info$organism == "par",])
infos_TFdel_allele <- list(cer = sample_info_allele[sample_info_allele$experiment == "LowN" &
                                                      sample_info_allele$allele == "cer",],
                           par = sample_info_allele[sample_info_allele$experiment == "LowN" &
                                                      sample_info_allele$allele == "par",])
counts_TFdel <- list(cer = counts_unnorm[,infos_TFdel$cer$sample_name], 
                     par = counts_unnorm[,infos_TFdel$par$sample_name])
counts_TFdel_allele <- list(cer = counts_unnorm_allele[,infos_TFdel_allele$cer$sample_name], 
                            par = counts_unnorm_allele[,infos_TFdel_allele$par$sample_name])

# final number of genes and samples (should all be the same number of rows)
dim(counts)
dim(counts_allele)
dim(counts_TFdel$cer)
dim(counts_TFdel$par)
dim(counts_TFdel_allele$cer)
dim(counts_TFdel_allele$par)

#### saving to .RData file ####
save(counts, sample_info,
     cer_biased_genes,
     par_biased_genes,
     both_biased_genes,
     biased_samples, file = "data_files/Cleaned_Count_Data.RData")
save(counts_allele, sample_info_allele, file = "data_files/Cleaned_Count_Data_AlleleSpecific.RData")
save(counts_TFdel, counts_TFdel_allele, infos_TFdel, infos_TFdel_allele,
     file = "data_files/Cleaned_TFdel_Unnormalized_Counts.RData")

################################# Archive #######################################
#### TF filtering: TF must have at least 2 replicates at each timepoint ####
# Archived b/c we'll filter all these out later b/c a standard deviation of <2 points is NA
# # cer
# goodTable_cer <- infos_TFdel$cer |> filter(genotype != "WT") |> 
#   select(condition) |> table()
# cer_vec <- sapply(goodTable_cer, \(x) {return(x >= 2)})
# goodConditions_cer <- rownames(goodTable_cer)[cer_vec]
# # par
# goodTable_par <- infos_TFdel$par |> filter(genotype != "WT") |> 
#   select(condition) |> table()
# par_vec <- sapply(goodTable_par, \(x) {return(x >= 2)})
# goodConditions_par <- rownames(goodTable_par)[par_vec]
# # hyc
# goodTable_hyc <- infos_TFdel_allele$cer |> filter(genotype != "WT") |> 
#   select(condition) |> table()
# hyc_vec <- sapply(goodTable_hyc, \(x) {return(x >= 2)})
# goodConditions_hyc <- rownames(goodTable_hyc)[hyc_vec]
# # hyp
# goodTable_hyp <- infos_TFdel_allele$par |> filter(genotype != "WT") |> 
#   select(condition) |> table()
# hyp_vec <- sapply(goodTable_hyp, \(x) {return(x >= 2)})
# goodConditions_hyp <- rownames(goodTable_hyp)[hyp_vec]
# # all
# goodConditions <- purrr::reduce(list(goodConditions_cer, goodConditions_par, goodConditions_hyc, goodConditions_hyp), 
#                                 .f = intersect) # also ensuring the same set of TFs between species and hybrid
# # verifying these TFs do have all 6 timepoints x replicates in all 4 species/alleles
# goodTable_cer[goodConditions] |> min()
# goodTable_par[goodConditions] |> min()
# goodTable_hyc[goodConditions] |> min()
# goodTable_hyp[goodConditions] |> min() # should all have min 2
# 
# # filtering
# counts_TFdel$cer <- counts_TFdel$cer[,(infos_TFdel$cer$condition %in% goodConditions) | (infos_TFdel$cer$genotype == "WT")]
# counts_TFdel$par <- counts_TFdel$par[,(infos_TFdel$par$condition %in% goodConditions) | (infos_TFdel$par$genotype == "WT")]
# infos_TFdel$cer <- infos_TFdel$cer[(infos_TFdel$cer$condition %in% goodConditions) | (infos_TFdel$cer$genotype == "WT"),]
# infos_TFdel$par <- infos_TFdel$par[(infos_TFdel$par$condition %in% goodConditions) | (infos_TFdel$par$genotype == "WT"),]
# counts_TFdel_allele$cer <- counts_TFdel_allele$cer[,(infos_TFdel_allele$cer$condition %in% goodConditions) | (infos_TFdel_allele$cer$genotype == "WT")]
# counts_TFdel_allele$par <- counts_TFdel_allele$par[,(infos_TFdel_allele$par$condition %in% goodConditions) | (infos_TFdel_allele$par$genotype == "WT")]
# infos_TFdel_allele$cer <- infos_TFdel_allele$cer[(infos_TFdel_allele$cer$condition %in% goodConditions) | (infos_TFdel_allele$cer$genotype == "WT"),]
# infos_TFdel_allele$par <- infos_TFdel_allele$par[(infos_TFdel_allele$par$condition %in% goodConditions) | (infos_TFdel_allele$par$genotype == "WT"),]

#### Filtering LowN replicates ####
# archived b/c didn't seem necessary to know which samples came from the
# same well. We only use TFdel genotypes for the last paper figure, and
# the analysis is just happy to have however many replicates we can get
# # throwing out any well_flask_IDs that don't have all 3 timepoints for an organism
# LowNTimepoints <- sample_info %>% filter(experiment == "LowN") %>% select(time_point_str) %>% pull() %>% unique()
# wellHasAll3Timepoints <- function(wellID, org) {
#   if(wellID == "rep1" | wellID == "rep2") {
#     return(NA)
#   }
#   if (org != "hyb") {
#     timepoints <- filter(sample_info, well_flask_ID == wellID & organism == org) %>% select(time_point_str) %>% pull()
#   }
#   if (org == "hyb") {
#     timepoints_hyc <- filter(sample_info, well_flask_ID == wellID & organism == org & allele == "cer") %>% select(time_point_str) %>% pull()
#     timepoints_hyp <- filter(sample_info, well_flask_ID == wellID & organism == org & allele == "par") %>% select(time_point_str) %>% pull()
#     timepoints <- c(timepoints_hyc, timepoints_hyp)
#   }
#   return(setequal(unique(timepoints), LowNTimepoints))
# }
# # run it on all samples (takes a moment)
# well_has_all_3_timepoints <- map2_lgl(sample_info$well_flask_ID, sample_info$organism, .f = wellHasAll3Timepoints)
# 
# # tests for wellHasAll3Timepoints
# # positive test (should have all 3 timepoints)
# test <- sample_info[sapply(well_has_all_3_timepoints, isTRUE),] %>% select(well_flask_ID, organism) %>% slice_sample(n = 1) 
# sample_info %>% filter(well_flask_ID == test$well_flask_ID & organism == test$organism) %>% select(sample_name, time_point_str) 
# 
# # negative test (should not have all three timepoints)
# test <- sample_info[sapply(well_has_all_3_timepoints, isFALSE),] %>% select(well_flask_ID, organism) %>% slice_sample(n = 1) 
# sample_info %>% filter(well_flask_ID == test$well_flask_ID & organism == test$organism) %>% select(sample_name, time_point_str)
# 
# # or systematically, grouping by organism and replicate
# test <- sample_info[sapply(well_has_all_3_timepoints, isTRUE),] %>% 
#   select(time_point_str, well_flask_ID, organism) %>% group_by(organism, well_flask_ID) %>% 
#   summarise(nTimepointsRepresented = n_distinct(time_point_str))
# test %>% select(nTimepointsRepresented) %>% table()
# test <- sample_info[sapply(well_has_all_3_timepoints, isFALSE),] %>% 
#   select(time_point_str, well_flask_ID, organism) %>% group_by(organism, well_flask_ID) %>% 
#   summarise(nTimepointsRepresented = n_distinct(time_point_str))
# test %>% select(nTimepointsRepresented) %>% table()
# 
# # LowN should only have TRUE/FALSE
# well_has_all_3_timepoints[sample_info$experiment == "LowN"] %>% table()
# # and experiments besides LowN should only have NA
# well_has_all_3_timepoints[sample_info$experiment != "LowN"] %>% unique()
# 
# # removing samples without all 3 timepoints
# sum(!well_has_all_3_timepoints, na.rm = TRUE) # number of samples we're removing
# sample_info$sample_name[sapply(well_has_all_3_timepoints, isFALSE)] # names of samples we're removing
# 
# keep <- (well_has_all_3_timepoints | is.na(well_has_all_3_timepoints))
# sample_info <- sample_info[keep,]
# counts <- counts[,colnames(counts) %in% sample_info$sample_name]
# counts_unnorm <- counts_unnorm[,colnames(counts_unnorm) %in% sample_info$sample_name]
# 
# # replacing well_flask_IDs for the LowN samples from things like "C5_G12" to "rep1", "rep2", etc.
# well_flask_systematic_names <- sample_info$well_flask_ID # we need to create this vector separately because you can't overwrite the names until after every well has its new name or else unique includes each rep# in its set of unique values and you end up with things like rep313 as I found out
# for (i in 1:nrow(sample_info)) {
#   if (sample_info$experiment[i] == "LowN") {
#     org <- sample_info$organism[i]
#     gen <- sample_info$genotype[i]
#     well <- sample_info$well_flask_ID[i]
#     all_wells <- sample_info[sample_info$experiment == "LowN" &
#                                sample_info$organism == org &
#                                sample_info$genotype == gen,]$well_flask_ID
#     cat("index:", i, "Wells associated with", org, gen, ":", unique(all_wells), "\n")
#     welltab <- table(all_wells) %>% sort(decreasing = TRUE)
#     repNum <- which(names(welltab) == well) 
#     
#     well_flask_systematic_names[i] <- paste0("rep", repNum)
#   }
# }
# table(well_flask_systematic_names) %>% sort(decreasing = TRUE) # checking that they're all some sort of rep# now and that the lower reps have more samples (yes there literally are 34 hybrid WT wells)
# 
# sample_info$well_flask_ID <- well_flask_systematic_names
# 
# test <- sample_info %>% filter(experiment == "LowN") %>% select(well_flask_ID, organism, genotype) %>% mutate(DefOneWell = paste(well_flask_ID, organism, genotype))
# test_table <- table(test$DefOneWell) %>% as.list()
# # all three timepoints represented
# test_table[test_table == 3] %>% names() # all cer/par
# test_table[test_table == 6] %>% names() # these are all hybrids
# 
# # these are the remaining troublemakers that we need to get down to 3 samples (9 for hybrid b/c of 3 allele options)
# test_table[test_table != 3 & test_table != 6]
# 
# # whittle LowN replicate groups down so that every one has exactly 3 samples, 1 from each timepoint
# # removes samples where the name does not match at the middle value of the sample name (P1, P2, C1, C3, etc. what I'm calling the batch. It seems to match for three samples per replicate, one from each timepoint)
# # @input: a group of sample names from LowN that are supposed to be part of the same rep group
# # @output: the three sample names that are true replicates (their name matches besides the TP# and GS# parts)
# filterForTrueReps <- function(sample_names) {
#   batch <- sapply(sample_names, function(s) {
#     output <- strsplit(s, split = "_") %>% unlist()
#     return(output[4]) # The 4th section of the sample name has the number that doesn't vary for true reps
#   })
#   isTrueRep <- batch == names(which.max(table(batch)))
#   return(sample_names[isTrueRep])
# }
# 
# # tests for filterForTrueReps/equalizeNameForTrueReps
# test <- sample_info %>% filter(well_flask_ID == "rep2" & experiment == "LowN" & genotype == "MSN2delete" & organism == "par") %>% 
#   select(sample_name) %>% pull() 
# test # this par rep has 6 replicates and either three can be our true replicates
# filterForTrueReps(test) # this shoudl filter down to only 3 samples
# test <- sample_info %>% filter(well_flask_ID == "rep5" & experiment == "LowN" & genotype == "WT" & organism == "cer") %>% 
#   select(sample_name) %>% pull()
# test # this cer has 5 samples
# filterForTrueReps(test) # only 3 with matching names should remain
# test <- sample_info %>% filter(well_flask_ID == "rep3" & experiment == "LowN" & genotype == "WT" & organism == "hyb") %>% 
#   select(sample_name) %>% pull()
# test # shouldn't eliminate any samples from a good hybrid set of 6
# filterForTrueReps(test)
# 
# # now to apply it to all the troublemakers... with nested for loops!
# GoodLowNSampleNames <- vector(mode = "character", length = 0)
# for (org in c("cer", "par", "hyb")) {
#   for (gen in unique(sample_info$genotype)) {
#     for (repl in unique(sample_info$well_flask_ID)) {
#       cat("currently processing:", org, gen, repl, "\n")
#       sample_names <- sample_info %>% 
#         filter(well_flask_ID == repl & experiment == "LowN" & genotype == gen & organism == org) %>%
#         select(sample_name) %>% pull()
#       if (length(sample_names) == 0) {
#         next
#       }
#       if (org == "hyb") {
#         desired_length <- 6
#       }
#       if (org != "hyb") {
#         desired_length <- 3
#       }
#       if (length(sample_names) != desired_length) {
#         sample_names <- filterForTrueReps(sample_names)
#       }
#       if (length(sample_names) == desired_length) {
#         GoodLowNSampleNames <- c(GoodLowNSampleNames, sample_names)
#       }
#       else {
#         stop(paste("filtered samples are not desired length for samples:", unlist(sample_names), "in rep:", repl, "\n"))
#       }
#     }
#   }
# }
# 
# # filtering to true replicates (3 and only 3 timepoints per well ID)
# keep <- sample_info$sample_name %in% GoodLowNSampleNames | sample_info$experiment != "LowN"
# sum(!keep) # number of samples we're removing
# sample_info <- sample_info[keep,]
# keep_counts <- colnames(counts) %in% sample_info$sample_name
# counts <- counts[,keep_counts]
# counts_unnorm <- counts_unnorm[,keep_counts]
# 
# # checking all the reps are now only 3 samples (9 for hyb)
# test <- sample_info %>% filter(experiment == "LowN") %>% select(well_flask_ID, organism, genotype) %>% mutate(DefOneWell = paste(well_flask_ID, organism, genotype))
# test_table <- table(test$DefOneWell) %>% as.list()
# test_table[test_table == 3] %>% names() # all cer/par
# test_table[test_table == 6] %>% names() # all hyb
# test_table[test_table != 3 & test_table != 6] |> names() # shouldn't have any

#### more conservative TF filtering: TF must have all 3 timepoints in at least 2 replicates ####
# # cer
# goodTable_cer <- infos_TFdel$cer |> filter(genotype != "WT") |> 
#   select(genotype, time_point_str) |> table()
# cer_vec <- apply(goodTable_cer, 1, \(x) {all(x >= 2)})
# goodTFs_cer <- rownames(goodTable_cer)[cer_vec]
# # par
# goodTable_par <- infos_TFdel$par |> filter(genotype != "WT") |> 
#   select(genotype, time_point_str) |> table()
# par_vec <- apply(goodTable_par, 1, \(x) {all(x >= 2)})
# goodTFs_par <- rownames(goodTable_par)[par_vec]
# # hyc
# goodTable_hyc <- infos_TFdel_allele$cer |> filter(genotype != "WT") |> 
#   select(genotype, time_point_str) |> table()
# hyc_vec <- apply(goodTable_hyc, 1, \(x) {all(x >= 2)})
# goodTFs_hyc <- rownames(goodTable_hyc)[hyc_vec]
# # hyp
# goodTable_hyp <- infos_TFdel_allele$par |> filter(genotype != "WT") |> 
#   select(genotype, time_point_str) |> table()
# hyp_vec <- apply(goodTable_hyp, 1, \(x) {all(x >= 2)})
# goodTFs_hyp <- rownames(goodTable_hyp)[hyp_vec]
# # all
# goodTFs <- reduce(list(goodTFs_cer, goodTFs_par, goodTFs_hyc, goodTFs_hyp), 
#                   .f = intersect) # also ensuring the same set of TFs between species and hybrid
# # verifying these TFs do have all 6 timepoints x replicates in all 4 species/alleles
# goodTable_cer[goodTFs,] |> min()
# goodTable_par[goodTFs,] |> min()
# goodTable_hyc[goodTFs,] |> min()
# goodTable_hyp[goodTFs,] |> min() # should all have min 2
# 
# # filtering
# counts_TFdel$cer <- counts_TFdel$cer[(infos_TFdel$cer$genotype %in% goodTFs) | (infos_TFdel$cer$genotype == "WT"),]
# counts_TFdel$par <- counts_TFdel$par[(infos_TFdel$par$genotype %in% goodTFs) | (infos_TFdel$par$genotype == "WT"),]
# infos_TFdel$cer <- infos_TFdel$cer[(infos_TFdel$cer$genotype %in% goodTFs) | (infos_TFdel$cer$genotype == "WT"),]
# infos_TFdel$par <- infos_TFdel$par[(infos_TFdel$par$genotype %in% goodTFs) | (infos_TFdel$par$genotype == "WT"),]
# counts_TFdel_allele$cer <- counts_TFdel_allele$cer[(infos_TFdel_allele$cer$genotype %in% goodTFs) | (infos_TFdel_allele$cer$genotype == "WT"),]
# counts_TFdel_allele$par <- counts_TFdel_allele$par[(infos_TFdel_allele$par$genotype %in% goodTFs) | (infos_TFdel_allele$par$genotype == "WT"),]
# infos_TFdel_allele$cer <- infos_TFdel_allele$cer[(infos_TFdel_allele$cer$genotype %in% goodTFs) | (infos_TFdel_allele$cer$genotype == "WT"),]
# infos_TFdel_allele$par <- infos_TFdel_allele$par[(infos_TFdel_allele$par$genotype %in% goodTFs) | (infos_TFdel_allele$par$genotype == "WT"),]

#### Troubleshooting hybrid libsizes ####
# Archived b/c differences in libsize don't appear to be causing any biased counts
# #### QC: checking if hybrids have a bias to be higher raw counts or cpm than parents ####
# # motivation: hybrids tend to have higher LFC between WT and TFdel (and mostly positive). 
# # I'm wondering if this is a library size/normalization issue. 
# # normalizing shouldn't affect LFC, as long as there isn't a bias for WT or TFdel
# # to have larger libraries than the other one
# plotdf <- bind_cols(tibble(libsize = colSums(counts)), sample_info) |> 
#   filter(experiment == "LowN") |> 
#   mutate(isWT = genotype == "WT") |> 
#   mutate(type = paste(organism, allele, isWT))
# ggplot(plotdf, aes(x = type, y = libsize)) + 
#   geom_jitter(aes(color = time_point_str)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   geom_hline(yintercept = 100000)
# # there's a bias for the 16 h timepoint
# # to have larger libraries, but *only* in hyb WT, not even deletion
# # cer has a slight bias for 0 h having the largest library, but that's true of both
# # WT and tfdel, so it won't affect LFC between them
# 
# # This is most likely causing more genes to be ID'd as DE in hybrid b/c one TP is high
# # in WT and the other 2 are low versus TFdel, but that theoretically shouldn't happen
# # Because if each gene's raw expression level is directly proportional to its library size,
# # normalizing by lib size should make genes more comparable between samples
# # unless raw counts *aren't* directly proportional to libsize
# 
# # Do certain hybrid TF deletions have larger libraries than others?
# # We see the bias for higher LFC in hybrids is stronger in some TFs than others:
# # Worst offenders: CHA4 (TP2 and TP3), GCR2 (all timepoints), GLN3 (TP2 and TP3), GZF3 (TP1 and TP2),
# #                  HAP2 (TP1 and TP2), INO4 (TP3), TEC1 (TP1 and TP2)
# # least bad: LEU3 (all timepoints, maybe TP3 is a bit bad), MIG1 (TP2 and TP3), MSN2 (TP1 and TP2),
# #                 NRG1 (TP1 and TP2), ROX1, SOK2, AFT1
# # Is this pattern reflected in library size at all?
# # hyb, sum of both alleles, 16h
# plotdf_hyb <- plotdf |> filter(organism == "hyb" & time_point_num == 960) |> 
#   group_by(genotype, well_flask_ID, time_point_str) |>
#   summarise(hyb_libsize = sum(libsize))
# ggplot(plotdf_hyb, aes(x = genotype, y = hyb_libsize)) + 
#   geom_point(aes(color = genotype)) +
#   theme(legend.position = "none", axis.text.x = element_text(angle = 90))
# # cer 16h
# plotdf_cer <- plotdf |> filter(organism == "cer" & time_point_num == 960)
# ggplot(plotdf_cer, aes(x = genotype, y = libsize)) + 
#   geom_point(aes(color = genotype)) +
#   theme(legend.position = "none", axis.text.x = element_text(angle = 90))
# # par 16h
# plotdf_par <- plotdf |> filter(organism == "par" & time_point_num == 960)
# ggplot(plotdf_par, aes(x = genotype, y = libsize)) + 
#   geom_point(aes(color = genotype)) +
#   theme(legend.position = "none", axis.text.x = element_text(angle = 90))
# 
# 
# # let's look at an example of a gene we know has a problem --- YLR053C in NRG1 delete
# gene_idx <- "YLR053C"
# # raw counts
# genedf_wt <- bind_cols(tibble(expr = as.numeric(counts[gene_idx, sample_info$organism == "hyb" &
#                                                          sample_info$allele == "par" &
#                                                          sample_info$genotype == "WT" &
#                                                          sample_info$experiment == "LowN"])),
#                        sample_info[sample_info$organism == "hyb" &
#                                      sample_info$allele == "par" &
#                                      sample_info$genotype == "WT" &
#                                      sample_info$experiment == "LowN",])
# genedf_del <- bind_cols(tibble(expr = as.numeric(counts[gene_idx, sample_info$organism == "hyb" &
#                                                           sample_info$allele == "par" &
#                                                           sample_info$genotype == "NRG1delete" &
#                                                           sample_info$experiment == "LowN"])),
#                         sample_info[sample_info$organism == "hyb" &
#                                       sample_info$allele == "par" &
#                                       sample_info$genotype == "NRG1delete" &
#                                       sample_info$experiment == "LowN",])
# 
# plotdf <- bind_rows(genedf_wt, genedf_del)
# ggplot(plotdf, aes(x = time_point_num, y = expr)) + 
#   geom_point(aes(color = genotype)) +
#   geom_line(data = summarise(group_by(plotdf, genotype, time_point_num),
#                              mean_expr = mean(expr)), aes(x = time_point_num, y = mean_expr, group = genotype,
#                                                           color = genotype))
# # no difference btwn WT and hyc or hyp
# 
# # normalized counts
# countsPerMillion <- function(.cts) {
#   librarySizes <- colSums(.cts, na.rm = TRUE)
#   output <- apply(.cts, 1, function(x) {
#     normalized <- (x/librarySizes)*1e6
#     return(round(normalized))
#   })
#   return(t(output))
# }
# cpm <- countsPerMillion(counts)
# genedf_wt <- bind_cols(tibble(expr = as.numeric(cpm[gene_idx, sample_info$organism == "hyb" &
#                                                       sample_info$allele == "cer" &
#                                                       sample_info$genotype == "WT" &
#                                                       sample_info$experiment == "LowN"])),
#                        sample_info[sample_info$organism == "hyb" &
#                                      sample_info$allele == "cer" &
#                                      sample_info$genotype == "WT" &
#                                      sample_info$experiment == "LowN",])
# genedf_del <- bind_cols(tibble(expr = as.numeric(cpm[gene_idx, sample_info$organism == "hyb" &
#                                                        sample_info$allele == "cer" &
#                                                        sample_info$genotype == "NRG1delete" &
#                                                        sample_info$experiment == "LowN"])),
#                         sample_info[sample_info$organism == "hyb" &
#                                       sample_info$allele == "cer" &
#                                       sample_info$genotype == "NRG1delete" &
#                                       sample_info$experiment == "LowN",])
# 
# plotdf <- bind_rows(genedf_wt, genedf_del)
# ggplot(plotdf, aes(x = time_point_num, y = expr)) + 
#   geom_point(aes(color = genotype)) +
#   geom_line(data = summarise(group_by(plotdf, genotype, time_point_num),
#                              mean_expr = mean(expr)), aes(x = time_point_num, y = mean_expr, group = genotype,
#                                                           color = genotype))
# # well well well
# # so the issue appears to be the fact that the larger WT libsizes at 16 h 
# # tend to lower WT cpm counts at 16 h relative to the TFdel
# 
# # But if we assume WT and TFdel raw counts are comparable, where
# # is WT making up all that extra libsize from?
# # are certain genes increasing expr dramatically while most are staying the same?
# 
# # plot of unnormalized 0h/1h versus 16h WT hyb gene versus gene to see if certain
# # genes jump up a lot in the 16 h while most stay the same
# counts_0h1h <- counts[, sample_info$organism == "hyb" &
#                         sample_info$allele == "cer" &
#                         sample_info$genotype == "WT" &
#                         sample_info$experiment == "LowN" &
#                         sample_info$time_point_num %in% c(0, 60)]
# counts_16h <- counts[, sample_info$organism == "hyb" &
#                        sample_info$allele == "cer" &
#                        sample_info$genotype == "WT" &
#                        sample_info$experiment == "LowN" &
#                        sample_info$time_point_num == 960]
# plotdf <- tibble(expr_0h1h = rowMeans(counts_0h1h),
#                  expr_16h = rowMeans(counts_16h))
# plot(log2(plotdf$expr_0h1h), log2(plotdf$expr_16h))
# abline(a = 0, b = 1, col = "red")
# # there's no obvious genes that are causing counts to go up so much,
# # might be apparent for random pairs of samples instead of means across all samples
# # keep re-running to check different pairs of random samples
# plotdf <- tibble(expr_0h1h = counts_0h1h[, sample(c(1:ncol(counts_0h1h)), 1)],
#                  expr_16h = counts_16h[, sample(c(1:ncol(counts_16h)), 1)],
#                  gene_name = rownames(counts))
# plot(plotdf$expr_0h1h, plotdf$expr_16h)
# abline(a = 0, b = 1, col = "red")
# h_cutoff <- 1200
# v_cutoff <- 150
# abline(h = h_cutoff, col = "blue") # adjust these so they select for the outlier group
# abline(v = v_cutoff, col = "blue")
# problem_genes <- c(problem_genes, plotdf$gene_name[plotdf$expr_0h1h < v_cutoff & 
#                                                      plotdf$expr_16h > h_cutoff])
# names(table(problem_genes)[table(problem_genes) > 1])
# # looking these up results in mostly nitrogen metabolic genes, which
# # makes sense b/c those could really be activated in the 16h conditon versus the 0/1h
# 
# ### Comparing LFC in WT versus TFdel hybrids
# # calculating LFC in 0/1h versus 16h (1 gene LFC = log2(mean(16h)) - log2(mean(0/1h)))
# # and comparing these LFCs in hyb WT versus the deletion genotypes 
# # (all deletions treated the same as WT replicates b/c we're assuming)
# # that most genes actually responding to low nitrogen are going to
# # still respond in most TFdel mutants
# 
# # first collapsing gene counts into hyb --- hyc + hyp and taking means
# hyb_sample_info <- filter(sample_info, organism == "hyb") |> 
#   select(setdiff(colnames(sample_info), "allele"))
# hyb_sample_info$sample_name <- gsub(pattern = "_hy[pc]_", 
#                                     replacement = "_hyb_",
#                                     x = hyb_sample_info$sample_name)
# hyb_sample_info <- unique(hyb_sample_info)
# hyc_counts <- counts[,sample_info$organism == "hyb" & sample_info$allele == "cer"]
# colnames(hyc_counts) <- gsub(pattern = "_hyc_", 
#                              replacement = "_hyb_",
#                              x = colnames(hyc_counts))
# hyp_counts <- counts[,sample_info$organism == "hyb" & sample_info$allele == "par"]
# colnames(hyp_counts) <- gsub(pattern = "_hyp_", 
#                              replacement = "_hyb_",
#                              x = colnames(hyp_counts))
# sum(colnames(hyc_counts) == colnames(hyp_counts)) # should be 516
# hyb_counts <- hyc_counts + hyp_counts
# # now splitting into WT versus TFdel and 0h1h versus 16h, 
# # collapsing each gene into 1 mean in each group (4 means total per gene)
# plotdf <- tibble(LFC_WT = log2(rowMeans(hyb_counts[, hyb_sample_info$genotype == "WT" &
#                                                      hyb_sample_info$time_point_num == 960])) -
#                    log2(rowMeans(hyb_counts[, hyb_sample_info$genotype == "WT" &
#                                               hyb_sample_info$time_point_num != 960])),
#                  LFC_TFdel = log2(rowMeans(hyb_counts[, hyb_sample_info$genotype != "WT" &
#                                                         hyb_sample_info$time_point_num == 960])) -
#                    log2(rowMeans(hyb_counts[, hyb_sample_info$genotype != "WT" &
#                                               hyb_sample_info$time_point_num != 960])))
# plot(plotdf$LFC_WT, plotdf$LFC_TFdel, xlab = "LFC WT", ylab = "LFC TFdel")
# abline(a = 0, b = 1, col = "red")
# # no obvious outliers, but a clear bias for lower LFC in TFdel across the board
# # This makes sense given that counts are not yet normalized, so the fact that 
# # 16h libraries are systematically higher in WT means that LFC from 0/1h to 16h 
# # will be systematically higher
# # this also confirms that the increase in library size isn't the result of 
# # certain genes --- it's a systematic increase in every gene
# # so normalizing should work
# # normalizing
# hyb_counts <- countsPerMillion(hyb_counts)
# plotdf <- tibble(LFC_WT = log2(rowMeans(hyb_counts[, hyb_sample_info$genotype == "WT" &
#                                                      hyb_sample_info$time_point_num == 960]) + 1) -
#                    log2(rowMeans(hyb_counts[, hyb_sample_info$genotype == "WT" &
#                                               hyb_sample_info$time_point_num != 960]) + 1),
#                  LFC_TFdel = log2(rowMeans(hyb_counts[, hyb_sample_info$genotype != "WT" &
#                                                         hyb_sample_info$time_point_num == 960]) + 1) -
#                    log2(rowMeans(hyb_counts[, hyb_sample_info$genotype != "WT" &
#                                               hyb_sample_info$time_point_num != 960]) + 1))
# plot(plotdf$LFC_WT, plotdf$LFC_TFdel, xlab = "LFC WT", ylab = "LFC TFdel")
# abline(a = 0, b = 1, col = "red")
# # that mostly worked. The highest LFCs still appear to be higher
# # in WT however. Not the same sort of bias in TFdel
# 
# # here's the same data to make sure the density of genes isn't obsuring a bias:
# ggplot(plotdf, aes(x = LFC_WT, y = LFC_TFdel)) + 
#   geom_point(alpha = 0.02) + 
#   geom_abline(slope = 1, intercept = 0, color = "red")
# ggplot(filter(plotdf, abs(LFC_WT) > 0.5 | abs(LFC_TFdel) > 0.5), 
#        aes(x = LFC_WT, y = LFC_TFdel)) + 
#   geom_hex() + 
#   geom_abline(slope = 1, intercept = 0, color = "red")
# 
# sum(plotdf$LFC_WT > plotdf$LFC_TFdel, na.rm = TRUE)
# sum(plotdf$LFC_WT <= plotdf$LFC_TFdel, na.rm = TRUE)
# # about twice as many LFCs are higher in WT than TFdel
# # so *most* genes have higher counts in 16h WT than 16h TFdel
# # even after normalization.
# cor(x = plotdf$LFC_WT, y = plotdf$LFC_TFdel, use = "pairwise.complete.obs")
# 
# 
# # how does this compare to the cer and par parental plots? Do those follow y=x better?
# # do they have the same unevenness of + and - LFCs?
# # unnormalized, cer
# plotdf <- tibble(LFC_WT = log2(rowMeans(counts[, sample_info$genotype == "WT" &
#                                                  sample_info$organism == "cer" &
#                                                  sample_info$time_point_num == 960])) -
#                    log2(rowMeans(counts[, sample_info$genotype == "WT" &
#                                           sample_info$organism == "cer" &
#                                           sample_info$time_point_num != 960])),
#                  LFC_TFdel = log2(rowMeans(counts[, sample_info$genotype != "WT" &
#                                                     sample_info$organism == "cer" &
#                                                     sample_info$time_point_num == 960])) -
#                    log2(rowMeans(counts[, sample_info$genotype != "WT" &
#                                           sample_info$organism == "cer" &
#                                           sample_info$time_point_num != 960])))
# plot(plotdf$LFC_WT, plotdf$LFC_TFdel, xlab = "LFC WT", ylab = "LFC TFdel")
# abline(a = 0, b = 1, col = "red")
# # normalized, cer
# plotdf <- tibble(LFC_WT = log2(rowMeans(cpm[, sample_info$genotype == "WT" &
#                                               sample_info$organism == "cer" &
#                                               sample_info$time_point_num == 960]) + 1) -
#                    log2(rowMeans(cpm[, sample_info$genotype == "WT" &
#                                        sample_info$organism == "cer" &
#                                        sample_info$time_point_num != 960]) + 1),
#                  LFC_TFdel = log2(rowMeans(cpm[, sample_info$genotype != "WT" &
#                                                  sample_info$organism == "cer" &
#                                                  sample_info$time_point_num == 960]) + 1) -
#                    log2(rowMeans(cpm[, sample_info$genotype != "WT" &
#                                        sample_info$organism == "cer" &
#                                        sample_info$time_point_num != 960]) + 1))
# plot(plotdf$LFC_WT, plotdf$LFC_TFdel, xlab = "LFC WT", ylab = "LFC TFdel")
# abline(a = 0, b = 1, col = "red")
# cor(x = plotdf$LFC_WT, y = plotdf$LFC_TFdel, use = "pairwise.complete.obs")
# sum(plotdf$LFC_WT > plotdf$LFC_TFdel, na.rm = TRUE)
# sum(plotdf$LFC_WT <= plotdf$LFC_TFdel, na.rm = TRUE)
# 
# # unnormalized, par
# plotdf <- tibble(LFC_WT = log2(rowMeans(counts[, sample_info$genotype == "WT" &
#                                                  sample_info$organism == "par" &
#                                                  sample_info$time_point_num == 960])) -
#                    log2(rowMeans(counts[, sample_info$genotype == "WT" &
#                                           sample_info$organism == "par" &
#                                           sample_info$time_point_num != 960])),
#                  LFC_TFdel = log2(rowMeans(counts[, sample_info$genotype != "WT" &
#                                                     sample_info$organism == "par" &
#                                                     sample_info$time_point_num == 960])) -
#                    log2(rowMeans(counts[, sample_info$genotype != "WT" &
#                                           sample_info$organism == "par" &
#                                           sample_info$time_point_num != 960])))
# plot(plotdf$LFC_WT, plotdf$LFC_TFdel, xlab = "LFC WT", ylab = "LFC TFdel")
# abline(a = 0, b = 1, col = "red")
# 
# # normalized, par
# plotdf <- tibble(LFC_WT = log2(rowMeans(cpm[, sample_info$genotype == "WT" &
#                                               sample_info$organism == "par" &
#                                               sample_info$time_point_num == 960] + 1)) -
#                    log2(rowMeans(cpm[, sample_info$genotype == "WT" &
#                                        sample_info$organism == "par" &
#                                        sample_info$time_point_num != 960]) + 1),
#                  LFC_TFdel = log2(rowMeans(cpm[, sample_info$genotype != "WT" &
#                                                  sample_info$organism == "par" &
#                                                  sample_info$time_point_num == 960]) + 1) -
#                    log2(rowMeans(cpm[, sample_info$genotype != "WT" &
#                                        sample_info$organism == "par" &
#                                        sample_info$time_point_num != 960]) + 1))
# plot(plotdf$LFC_WT, plotdf$LFC_TFdel, xlab = "LFC WT", ylab = "LFC TFdel")
# abline(a = 0, b = 1, col = "red")
# # parent LFCs are perhaps better correlated between WT and TFdel, but it's hard to say
# cor(x = plotdf$LFC_WT, y = plotdf$LFC_TFdel, use = "pairwise.complete.obs")
# sum(plotdf$LFC_WT > plotdf$LFC_TFdel, na.rm = TRUE)
# sum(plotdf$LFC_WT <= plotdf$LFC_TFdel, na.rm = TRUE)
# 
# # conclusions: no obvious differences between parents and hybrids in terms of how most genes
# # change expression at 16h in WT versus mean across deletions
# # likely the 16 h WT hyb libsize bias doesn't have an effect
#### troubleshooting YBR085W in TFdel ####
# # YBR085W seems like a good poster child gene for how counts seem
# # to get more dramatically different in hybrid TFdel (any TFdel)
# # versus WT, and I want to make sure this isn't due to
# # a normalization issue with WT 16h having systematically higher libsizes
# # hyc, WT vs TFdel counts
# gene_idx <- "YBR085W"
# plotdf <- bind_rows(bind_cols(tibble(expr = as.numeric(counts[gene_idx, sample_info$organism == "hyb" &
#                                                                 sample_info$allele == "cer" &
#                                                                 sample_info$experiment == "LowN"]),
#                                      libsize = colSums(counts[, sample_info$organism == "hyb" &
#                                                                 sample_info$allele == "cer" &
#                                                                 sample_info$experiment == "LowN"])),
#                               sample_info[sample_info$organism == "hyb" &
#                                             sample_info$allele == "cer" &
#                                             sample_info$experiment == "LowN",]),
#                     bind_cols(tibble(expr = as.numeric(counts[gene_idx, sample_info$organism == "hyb" &
#                                                                 sample_info$allele == "par" &
#                                                                 sample_info$experiment == "LowN"]),
#                                      libsize = colSums(counts[, sample_info$organism == "hyb" &
#                                                                 sample_info$allele == "par" &
#                                                                 sample_info$experiment == "LowN"])),
#                               sample_info[sample_info$organism == "hyb" &
#                                             sample_info$allele == "par" &
#                                             sample_info$experiment == "LowN",])) |> 
#   group_by(allele, genotype, time_point_str) |> 
#   summarise(mean_expr = mean(expr),
#             mean_libsize = mean(libsize))
# plotdf$group <- if_else(plotdf$genotype == "WT",
#                         true = if_else(plotdf$allele == "cer",
#                                        true = "cer WT",
#                                        false = "par WT"),
#                         false = if_else(plotdf$allele == "cer",
#                                         true = "cer TFdel",
#                                         false = "par TFdel"))
# ggplot(plotdf, aes(x = mean_libsize, y = mean_expr)) +
#   geom_point(aes(color = group))
# ggplot(plotdf, aes(x = mean_libsize, y = mean_expr)) +
#   geom_line(aes(group = paste(allele, genotype), 
#                 color = group))
# # here's the weirdness - the cerevisiae TFdels are systematically
# # higher expressed than cer WT
# 
# # Is this a specifically cerevisiae problem? What about
# # one of the genes that's consistently more DE in paradoxus?
# # gene_idx <- "YAR050W"
# gene_idx <- sample(rownames(counts), 1)
# plotdf <- bind_rows(bind_cols(tibble(expr = as.numeric(counts[gene_idx, sample_info$organism == "hyb" &
#                                                                 sample_info$allele == "cer" &
#                                                                 sample_info$experiment == "LowN"]),
#                                      libsize = colSums(counts[, sample_info$organism == "hyb" &
#                                                                 sample_info$allele == "cer" &
#                                                                 sample_info$experiment == "LowN"])),
#                               sample_info[sample_info$organism == "hyb" &
#                                             sample_info$allele == "cer" &
#                                             sample_info$experiment == "LowN",]),
#                     bind_cols(tibble(expr = as.numeric(counts[gene_idx, sample_info$organism == "hyb" &
#                                                                 sample_info$allele == "par" &
#                                                                 sample_info$experiment == "LowN"]),
#                                      libsize = colSums(counts[, sample_info$organism == "hyb" &
#                                                                 sample_info$allele == "par" &
#                                                                 sample_info$experiment == "LowN"])),
#                               sample_info[sample_info$organism == "hyb" &
#                                             sample_info$allele == "par" &
#                                             sample_info$experiment == "LowN",])) |> 
#   group_by(allele, genotype, time_point_str) |> 
#   summarise(mean_expr = mean(expr),
#             mean_libsize = mean(libsize))
# plotdf$group <- if_else(plotdf$genotype == "WT",
#                         true = if_else(plotdf$allele == "cer",
#                                        true = "cer WT",
#                                        false = "par WT"),
#                         false = if_else(plotdf$allele == "cer",
#                                         true = "cer TFdel",
#                                         false = "par TFdel"))
# ggplot(plotdf, aes(x = mean_libsize, y = mean_expr)) +
#   geom_point(aes(color = group))
# ggplot(plotdf, aes(x = mean_libsize, y = mean_expr)) +
#   geom_line(aes(group = paste(allele, genotype), 
#                 color = group))
# # YPL225W identified randomly as chaperone that's expressed lower in cer WT
# # than any TF deletion (and it's highly expressed)
# # TDH3 WT systemtatically higher than all TFdeletions in par
# 
# # conclusions: this isn't a libsize issue, but some
# # genes have ALL TFdels expressed higher or lower than WT
