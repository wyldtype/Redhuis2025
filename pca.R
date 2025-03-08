setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")
source("functions_for_figure_scripts.R")
load("data_files/Cleaned_Counts.RData")
load("data_files/Cleaned_Counts_Allele.RData")
load("data_files/FinalDataframe3Disp.RData")
# PCAs here are either sample (reasonable dimension matrices, run quickly) 
# or gene PCAs (5000x5000 matrices, take 30min to run)

#### Pan-experiment sample PCA ####
# principal component clustering of each experiment in cer, par, hyc, hyp
# with paths drawn to see progression over time
# princomps
covmat <- cov(counts[,sample_info$genotype == "WT"])
colnames(covmat) <- colnames(counts[,sample_info$genotype == "WT"])
pca_res <- prcomp(covmat)
pcadf <- tibble(pc1 = pca_res$x[,1], 
                pc2 = pca_res$x[,2],
                sample_name = colnames(covmat)) |> 
  left_join(sample_info, by = "sample_name") |> 
  mutate(condition = if_else(time_point_num == 0,
                             true = "baseline",
                             false = experiment))
var_pct <- summary(pca_res)$importance[2, 1:2] # % variance explained

# While there is variation between replicates, there are 
# common timepoint trajectories
pathdf <- pcadf |> group_by(time_point_num, experiment, organism, allele) |> 
  summarise(mean_pc1 = mean(pc1),
            mean_pc2 = mean(pc2)) |>
  ungroup() |> 
  arrange(time_point_num)
pathdf$timepoint_type <- apply(pathdf, 1, \(x) {
  e <- x["experiment"]
  tp <- x["time_point_num"] |> as.numeric()
  min_tp <- pathdf |> filter(experiment == e) |> 
    select(time_point_num) |> pull() |> min()
  max_tp <- pathdf |> filter(experiment == e) |> 
    select(time_point_num) |> pull() |> max()
  cat(e, tp, min_tp, max_tp, "\n")
  if (tp == min_tp) {
    return("start")
  }
  if (tp == max_tp) {
    return("end")
  }
  else {
    return("middle")
  }
})
# plotting by experiment
for (e in unique(sample_info$experiment)) {
  pdf(file = paste0("../../aligning_the_molecular_phenotype/paper_figures/Supplement/PCApaths_",
                    e, ".pdf"), height = 5, width = 7)
  p <- ggplot() + 
    geom_text(data = filter(pathdf, experiment == e),
              aes(x = mean_pc1, y = mean_pc2,
                  color = organism,
                  label = time_point_num)) +
    geom_line(data = filter(pcadf, experiment == e), 
              aes(x = pc1, y = pc2,
                  group = interaction(experiment, organism, allele, time_point_num)),
              color = "grey60") +
    xlab(paste0("PC1, ", round(var_pct[1]*100, digits = 0), 
                "% of variance")) + 
    ylab(paste0("PC2, ", round(var_pct[2]*100, digits = 0), 
                "% of variance")) +
    geom_path(data = filter(pathdf, experiment == e), 
              aes(x = mean_pc1, y = mean_pc2,
                  group = interaction(experiment, organism, allele),
                  color = organism)) +
    geom_point(data = filter(pcadf, experiment == e),
               aes(x = pc1, y = pc2),
               color = "grey60") +
    theme(legend.title = element_blank())
  print(p)
  dev.off()
}

# all experiments, colored by experiment
ggplot() + 
  geom_text(data = filter(pathdf, timepoint_type %in% c("start", "end")),
            aes(x = mean_pc1, y = mean_pc2,
                color = experiment,
                label = time_point_num)) +
  xlab(paste0("PC1, ", round(var_pct[1]*100, digits = 0), 
              "% of variance")) + 
  ylab(paste0("PC2, ", round(var_pct[2]*100, digits = 0), 
              "% of variance")) +
  geom_path(data = pathdf, 
            aes(x = mean_pc1, y = mean_pc2,
                group = interaction(experiment, organism, allele),
                color = experiment)) +
  theme(legend.title = element_blank())

# all experiments, colored by organism
ggplot() + 
  geom_text(data = filter(pathdf, timepoint_type %in% c("start", "end")),
            aes(x = mean_pc1, y = mean_pc2,
                color = organism,
                label = time_point_num)) +
  xlab(paste0("PC1, ", round(var_pct[1]*100, digits = 0), 
              "% of variance")) + 
  ylab(paste0("PC2, ", round(var_pct[2]*100, digits = 0), 
              "% of variance")) +
  geom_path(data = pathdf, 
            aes(x = mean_pc1, y = mean_pc2,
                group = interaction(experiment, organism, allele),
                color = organism)) +
  theme(legend.title = element_blank())

# plotting pairs of experiments to see how different 
# environments trigger different departures
# (no replicates indicated b/c it got chaotic)
# heat/cold
ggplot() + 
  geom_text(data = filter(pathdf, experiment %in% c("Heat", "Cold")),
            aes(x = mean_pc1, y = mean_pc2,
                color = interaction(organism, experiment),
                label = time_point_num)) +
  xlab(paste0("PC1, ", round(var_pct[1]*100, digits = 0), 
              "% of variance")) + 
  ylab(paste0("PC2, ", round(var_pct[2]*100, digits = 0), 
              "% of variance")) +
  geom_path(data = filter(pathdf, experiment %in% c("Heat", "Cold")), 
            aes(x = mean_pc1, y = mean_pc2,
                group = interaction(experiment, organism, allele),
                color = interaction(organism, experiment))) +
  theme(legend.title = element_blank())

# heat/LowPi
ggplot() + 
  geom_text(data = filter(pathdf, experiment %in% c("Heat", "LowPi")),
            aes(x = mean_pc1, y = mean_pc2,
                color = interaction(organism, experiment),
                label = time_point_num)) +
  xlab(paste0("PC1, ", round(var_pct[1]*100, digits = 0), 
              "% of variance")) + 
  ylab(paste0("PC2, ", round(var_pct[2]*100, digits = 0), 
              "% of variance")) +
  geom_path(data = filter(pathdf, experiment %in% c("Heat", "LowPi")), 
            aes(x = mean_pc1, y = mean_pc2,
                group = interaction(experiment, organism, allele),
                color = interaction(organism, experiment))) +
  theme(legend.title = element_blank())

# HAP4/LowPi
ggplot() + 
  geom_text(data = filter(pathdf, experiment %in% c("HAP4", "LowPi")),
            aes(x = mean_pc1, y = mean_pc2,
                color = interaction(organism, experiment),
                label = time_point_num)) +
  xlab(paste0("PC1, ", round(var_pct[1]*100, digits = 0), 
              "% of variance")) + 
  ylab(paste0("PC2, ", round(var_pct[2]*100, digits = 0), 
              "% of variance")) +
  geom_path(data = filter(pathdf, experiment %in% c("HAP4", "LowPi")), 
            aes(x = mean_pc1, y = mean_pc2,
                group = interaction(experiment, organism, allele),
                color = interaction(organism, experiment))) +
  theme(legend.title = element_blank())

# LowPi/LowN
ggplot() + 
  geom_text(data = filter(pathdf, experiment %in% c("LowN", "LowPi")),
            aes(x = mean_pc1, y = mean_pc2,
                color = interaction(organism, experiment),
                label = time_point_num)) +
  xlab(paste0("PC1, ", round(var_pct[1]*100, digits = 0), 
              "% of variance")) + 
  ylab(paste0("PC2, ", round(var_pct[2]*100, digits = 0), 
              "% of variance")) +
  geom_path(data = filter(pathdf, experiment %in% c("LowN", "LowPi")), 
            aes(x = mean_pc1, y = mean_pc2,
                group = interaction(experiment, organism, allele),
                color = interaction(organism, experiment))) +
  theme(legend.title = element_blank())

# HAP4/LowN
ggplot() + 
  geom_text(data = filter(pathdf, experiment %in% c("LowN", "HAP4")),
            aes(x = mean_pc1, y = mean_pc2,
                color = interaction(organism, experiment),
                label = time_point_num)) +
  xlab(paste0("PC1, ", round(var_pct[1]*100, digits = 0), 
              "% of variance")) + 
  ylab(paste0("PC2, ", round(var_pct[2]*100, digits = 0), 
              "% of variance")) +
  geom_path(data = filter(pathdf, experiment %in% c("LowN", "HAP4")), 
            aes(x = mean_pc1, y = mean_pc2,
                group = interaction(experiment, organism, allele),
                color = interaction(organism, experiment))) +
  theme(legend.title = element_blank())
#### Two experiment sample PCAs ####
### TODO: Pan-experiment shows us which environments are the most similar,
# but it can get messy, so the goal here is to separate PC1 into a time response
# and PC2 into an environment response for use in the leave out genes PCA below

#### leave-out genes sample PCAs ####
### TODO: Sample covariance matrix removes knowledge of which genes are 
# allowing PCAs to pick up on time and environmental challenge
# As both species tend to respond the same way, it's probably conserved genes,
# but are they conserved 1s or 2s or something else?
# To investigate, remove genes that are part of each cluster (1-1s then 2-2s)
# and see how they change the Two experiment sample PCA results



