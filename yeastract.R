sapply(c("dplyr", "purrr", "tidyr", "ggpubr", "readr",
         "data.table", "ggplot2", "data.table"), require, character.only=TRUE)
setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2025/")
load("data_files/FinalDataframe3Disp.RData")
# Making idx lists
idxs <- finaldf |> filter(experiment == "LowN" & cer == 2 & par == 0) |> select(gene_name)
write_delim(idxs, delim = "\n", file = "gene_ontology/20_LowN.txt")

# conserved plastic
mat11 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_11.csv",
                    delim = ";")
mat22 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_22.csv",
                    delim = ";")
# conserved static
mat00 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_00.csv",
                    delim = ";")
# diverged plastic
mat12 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_12.csv",
                    delim = ";")
mat21 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_21.csv",
                    delim = ";")
# diverged, plastic in one species
mat10 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_10.csv",
                    delim = ";")
mat01 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_01.csv",
                    delim = ";")
mat20 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_20.csv",
                    delim = ";")
mat02 <- read_delim("../../aligning_the_molecular_phenotype/yeastract/RegulationMatrix_LowN_02.csv",
                    delim = ";")

### Connectivity
plotdf <- map2(list(mat11, mat22, mat12, mat21, mat10, mat01, mat20, mat02, mat00),
               list("11", "22", "12", "21", "10", "01", "20", "02", "00"), \(x, y) {
  output <- tibble(gene_name = x[1,-1],
                   nConnections = colSums(x[,-1]),
                   cluster = y)
  return(output)
}) |> purrr::reduce(.f = bind_rows)
ggplot(plotdf, 
       aes(x = cluster, y = nConnections)) +
  geom_violin(aes(color = cluster))


