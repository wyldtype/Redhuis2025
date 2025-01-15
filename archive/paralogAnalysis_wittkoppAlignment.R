
###############################################################
################## Wittkopp Count Data Import + Normalization
###############################################################
# import gene count matrix
sampleNames <- list.files(pattern=".txt") %>% 
  lapply(gsub, pattern="_counts.txt", replacement="") %>% 
  flatten_chr()

geneNames <-
  read.table(file = "WT_08_cer_CellCycle_rep1_counts.txt", 
             colClasses=c(NA, "NULL")) %>% 
  as.matrix()

gcm <- list.files(pattern=".txt") %>% 
  purrr::map_dfc(read.table, sep="\t", colClasses=c("NULL", NA)) %>% 
  magrittr::set_rownames(geneNames) %>% 
  magrittr::set_colnames(sampleNames)

# create DGElist object
sample_info <- read.xlsx("../../bioSample1to999.xlsx", na.strings="not applicable") %>% 
  bind_rows(read.xlsx("../../bioSample1000toEnd.xlsx", na.strings="not applicable"))

wt_sample_info <- sample_info[sample_info$sample_name %in% colnames(gcm),]

dgelist0 <- DGEList(gcm, samples=wt_sample_info) %>% 
  calcNormFactors()

# filter out samples with less than a set number of total reads (need to do more research on what an appropriate cutoff is)
libCutoff <- 50000
dgelist1 <- dgelist0[,dgelist0$samples$lib.size > libCutoff]
hist(dgelist0$samples$lib.size)
hist(dgelist1$samples$lib.size)

# filter lowly expressed genes (< 1 cpm)
dgelist1$samples$group <- as.factor(dgelist1$samples$experiment)
keep <- filterByExpr(dgelist1, group = dgelist1$samples$group, lib.size = dgelist1$samples$lib.size)
dgelist <- dgelist1[keep,,]
cpm <- cpm(dgelist)
lcpm <- cpm(dgelist, log=TRUE) # used for visualization

##################### Wittkopp Count Data Analysis

# @input: gene names and correlation method ("pearson", "spearman", "kendall")
# as character vectors
# @output: expression correlation between those two genes across the environments in count matrix
correlate_paralogs <- function(gene1, gene2, corrMethod) {
  countMatrix <- cpm
  if (!is.na(gene2) & gene1 %in% rownames(cpm) & gene2 %in% rownames(cpm)) {
    exprVec1 <- countMatrix[rownames(countMatrix) == gene1,] %>% 
      as.numeric()
    exprVec2 <- countMatrix[rownames(countMatrix) == gene2,] %>% 
      as.numeric()
    return(cor(exprVec1, exprVec2, method=corrMethod))
  } else {
    return(NA)
  }
}

# example correlation. As expected, WGD paralogs TDH2 and TDH3 are much more 
# correlated than either is with SSD paralog TDH1
tdh3 <- "YGR192C"
tdh2 <-"YJR009C"
tdh1 <- "YJL052W"
correlate_paralogs(tdh3, tdh2, "pearson")
correlate_paralogs(tdh1, tdh2, "pearson")
correlate_paralogs(tdh3, tdh1, "pearson")
# Spearman ones are similar but not the same
correlate_paralogs(tdh3, tdh2, "spearman")
correlate_paralogs(tdh1, tdh2, "spearman")
correlate_paralogs(tdh3, tdh1, "spearman")

######## Data exploration
sample_info <- read.xlsx("../../bioSample1to999.xlsx", na.strings="not applicable") %>% 
  bind_rows(read.xlsx("../../bioSample1000toEnd.xlsx", na.strings="not applicable"))

# bringing data back to the tidyverse
lcpm_for_plotting <- as.data.frame(lcpm) %>% 
  rownames_to_column(var="gene_name") %>% 
  pivot_longer(-"gene_name") %>% 
  filter(!grepl("__", gene_name)) %>% 
  left_join(sample_info[,c("sample_name", "experiment", "time_point")],
            by=c("name"="sample_name"))

cpm_for_plotting <- as.data.frame(cpm) %>% 
  rownames_to_column(var="gene_name") %>% 
  pivot_longer(-"gene_name") %>% 
  filter(!grepl("__", gene_name)) %>% 
  left_join(sample_info[,c("sample_name", "experiment", "time_point")],
            by=c("name"="sample_name"))

# density plot of log-cpm values
p <- ggplot(data=lcpm_for_plotting, aes(x=value)) + geom_density(aes(color=experiment)) + ylim(0, 0.6)
p

# expression level of all genes across all environments (I think too many genes to color code)
paralogs <- read.table("../../marchant_paralogs.txt", header = TRUE)

p <- lcpm_for_plotting[lcpm_for_plotting$gene_name %in% paralogs$P1 | 
                         lcpm_for_plotting$gene_name %in% paralogs$P2, ] %>% 
  ggplot(aes(x=name, y=value, group=gene_name)) + 
  geom_line(aes(color= gene_name==tdh3 | gene_name==tdh2)) +
  theme(legend.position = "none")
p

p <- cpm_for_plotting[cpm_for_plotting$gene_name %in% paralogs$P1 | 
                        cpm_for_plotting$gene_name %in% paralogs$P2, ] %>% 
  ggplot(aes(x=name, y=value, group=gene_name)) + 
  geom_line(aes(color= gene_name==tdh3 | gene_name==tdh2)) +
  theme(legend.position = "none")
p

# spaghetti plot of pairs of paralogs (that only exist in single pairs)
single_pair_paralogs <- filter(paralogs, Dupli.SSD_WGD == 0 & MultiSSD == 0)
single_pair_paralogs$pair <- paste(single_pair_paralogs$P1, single_pair_paralogs$P2, sep="_")
single_pair_lcpm <- filter(lcpm_for_plotting, gene_name %in% single_pair_paralogs$P1 |
                             gene_name %in% single_pair_paralogs$P2) %>% 
  arrange(experiment, time_point)

# @input: gene systematic name
# @output: it and its paralog's systematic name concatenated with "_"
# Note: this only works with genes that only have one paralog
find_pair <- function(gene) {
  single_pair_paralogs[gene == single_pair_paralogs$P1 |
                         gene == single_pair_paralogs$P2, "pair"] %>%
    return()
}

single_pair_lcpm$pair <- map_chr(single_pair_lcpm$gene_name, function(x) {
  find_pair(x) %>% 
    return()
})

p <- ggplot(data=single_pair_lcpm, aes(x=name, y=value, group=gene_name)) + 
  geom_line(aes(color=pair)) + 
  theme(legend.position = "none")
p

# sample pair of paralogs' coexpression
sample_pair <- single_pair_lcpm[single_pair_lcpm$pair == "YDL081C_YDL130W",]
sample_pair$condition <- paste(sample_pair$experiment, sample_pair$time_point, sep = "_")
sample_pair <- sample_pair[order(sample_pair$condition),]
sample_pair[which(sample_pair$value == min(sample_pair$value)),"condition"]
sample_pair[which(sample_pair$value == max(sample_pair$value)),"condition"]

ggplot(data=sample_pair, aes(x=condition, y=value, group=gene_name)) + geom_line(aes(color=experiment))
ggplot(data=sample_pair[sample_pair$condition == "SCtoLowPi_135 min, No Pi" |
                          sample_pair$condition == "CellCycle_230 min, YPD",], 
       aes(x=condition, y=value, group=gene_name)) + geom_point(aes(color=gene_name))

# spaghetti plots are cool, but you can't really tell what's going on, so here are
# boxplots of each sample
boxplot(lcpm)
boxplot(log(dgelist$counts)) #un-normalized
ggplot(data=arrange(lcpm_for_plotting, experiment),
       aes(x=name, y=value)) + geom_boxplot() + ylab("log(cpm)") + xlab("sample")

# MDS plot
dgelist$samples$replicate <- 
  paste(dgelist$samples$experiment, 
        dgelist$samples$time_point, sep="_") %>% 
  as.factor()
library(questionr)
replicate_freqs <- questionr::freq(dgelist$samples$replicate) %>% 
  rownames_to_column(var="replicate_name")
table(dgelist$samples$replicate)
ggplot(data=dgelist$samples, aes(x=replicate)) + geom_histogram(stat="count") # best I could do
glMDSPlot(lcpm, labels=dgelist$samples$replicate, 
          groups=dgelist$samples[,c("group", "replicate")])

######## Kuzmin et al. 2020 trigenic interaciton fraction
tgi <- read.xlsx("../../Kuzmin_etAl_2020_trigenicFractions.xlsx", 
                 startRow=2, cols=c(4,5,11), na.strings="NaN")

# percentage of Kuzmin trigenic interaction genes present in my data
nrow(cpm[rownames(cpm) %in% tgi$ORF2,])/nrow(tgi)

# adding expression correlation values to each of the tgi pairs
tgi$expr_corr <- map2_dbl(tgi$ORF1, tgi$ORF2, ~ correlate_paralogs(.x, .y, "pearson"))

# using Fisher-z transformation to convert R values to normally distributed values
fisher_z_transform <- function(x) {
  return(0.5*(log(1 + x) - log(1 - x)))
}
tgi$z_score <- map_dbl(tgi$expr_corr, fisher_z_transform)

# assuming I'm doing this right, which is a big assumption, expression correlation
# is not significantly higher for more redundant (TGI>0.4) paralog pairs

#pvalue <- t.test(tgi$z_score ~ tgi$TGI.frac > 0.4)$p.value[1]
pvalue <- wilcox.test(tgi$expr_corr ~ tgi$TGI.frac > 0.4)$p.value
p <- ggplot(data=na.rm(tgi), aes(x=(TGI.frac > 0.4), y=expr_corr), na.rm=TRUE) + 
  geom_boxplot() +
  geom_dotplot(binaxis = 'y',
               stackdir = 'center', 
               alpha = .75,
               binpositions = 'all',
               binwidth = 0.03,
               fill = "lightseagreen", 
               colour = "lightseagreen") +
  labs(x = "trigenic interaction fraction > 0.4",
       y = "expression correlation between paralog pairs") +
  annotate(geom="text", label=paste0("wilcox.test, p-value: \n", formatC(pvalue, digits=2, format="e"))
           , x=1.5, y=0.75)
p

####### SSD and WGD comparison
paralogs <- read.table("../../marchant_paralogs.txt", header = TRUE)
nrow(cpm[rownames(cpm) %in% paralogs$P1,])/nrow(paralogs)
paralogs$expr_corr <- map2_dbl(paralogs$P1, paralogs$P2, ~ correlate_paralogs(.x, .y, "pearson"))
paralogs$z_score <- map_dbl(paralogs$expr_corr, fisher_z_transform)

#pvalue <- t.test(paralogs$z_score ~ paralogs$Duplication)$p.value[1]
pvalue <- wilcox.test(paralogs$expr_corr ~ paralogs$Duplication)$p.value
p <- ggplot(data=paralogs[paralogs$Duplication != "S",], 
            aes(x=Duplication, y=expr_corr)) + 
  geom_boxplot() +
  geom_dotplot(binaxis = 'y',
               stackdir = 'center', 
               alpha = .75,
               binpositions = 'all',
               binwidth = 0.01,
               fill = "lightseagreen", 
               colour = "lightseagreen") +
  labs(x = "duplication type",
       y = "expression correlation between paralog pairs") +
  annotate(geom="text", label=paste0("wilcox.test, p-value: \n", formatC(pvalue, digits = 2, format="e"))
           , x=1.5, y=0.75)
p