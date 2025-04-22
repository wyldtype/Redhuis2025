sapply(c("dplyr", "purrr", "tidyr", "ggpubr", "readr",
         "data.table", "ggplot2", "data.table", "msir", 
         "WGCNA", "energy", "matrixStats"), require, character.only=TRUE)

# # for testing functions:
# setwd("/Users/annar/Documents/Wittkopp_Lab/networks/DDivergence/Redhuis2024/")
# load("data_files/Cleaned_Counts.RData")
# load("data_files/Cleaned_Counts_Allele.RData")
# load("data_files/CorrelationClustering.RData")
# load("data_files/TFdel.RData")

#### Taking means across replicates ####
# helper function that turns a counts and info pair into
# a tidy dataframe
makeDf <- function(.cts, .info, .join_by = "sample_name") {
  outdf <- as_tibble(.cts) |>
    bind_cols(tibble(gene_name = rownames(.cts))) |>
    pivot_longer(cols = colnames(.cts),
                 names_to = .join_by, values_to = "expr") |>
    left_join(.info,
              by = .join_by)
  return(outdf)
}

# takes mean expr in each condition (timepoint/experiment/allele)
collapseReplicates <- function(.info, .cts) {
  conditions_cer <- .info |> filter(allele == "cer") |> select(condition) |> pull() |> unique()
  conditions_par <- .info |> filter(allele == "par") |> select(condition) |> pull() |> unique()
  common_conditions <- intersect(conditions_cer, conditions_par)
  # cer
  collapsed_cts_cer <- map(common_conditions, \(cond) {
    samps <- .info |> filter(condition == cond & allele == "cer") |> 
      select(sample_name) |> pull()
    cts <- .cts[,samps, drop = FALSE] # in HAP4 and LowPi, which don't have replicates, this will take the mean of one sample
    return(rowMeans(cts))
  }) |> reduce(cbind)
  collapsed_cts_par <- map(common_conditions, \(cond) {
    samps <- .info |> filter(condition == cond & allele == "par") |> 
      select(sample_name) |> pull()
    cts <- .cts[,samps, drop = FALSE] 
    return(rowMeans(cts))
  }) |> reduce(cbind)
  colnames(collapsed_cts_cer) <- common_conditions
  colnames(collapsed_cts_par) <- common_conditions
  rownames(collapsed_cts_cer) <- rownames(.cts)
  rownames(collapsed_cts_par) <- rownames(.cts)
  collapsed <- list(cer = collapsed_cts_cer, par = collapsed_cts_par)
  return(collapsed)
}
# # tests for collapseReplicates
# # parents
# dim(counts)
# test_collapse <- collapseReplicates(sample_info, counts)
# dim(cbind(test_collapse$cer, test_collapse$par))
# # LowPi, shouldn't have changed much
# # (only a couple samples have replicates, and only in par)
# sample_info |> filter(experiment == "LowPi" & organism == "par") |>
#   select(condition) |> table() |> sort(decresing = TRUE)
# # cer first, should not have changed
# test1 <- counts[,sample_info$experiment == "LowPi" &
#                   sample_info$allele == "cer"]
# dim(test1)
# test2 <- test_collapse$cer[,grepl("LowPi", colnames(test_collapse$cer))]
# dim(test2)
# filter(sample_info, sample_name %in% colnames(test1)) |>
#   select(condition) |> pull() |> setdiff(y = colnames(test2))
# # now par, also should not have changed
# test1 <- counts[,sample_info$experiment == "LowPi" &
#                   sample_info$allele == "par"]
# dim(test1)
# test2 <- test_collapse$par[,grepl("LowPi", colnames(test_collapse$par))]
# dim(test2)
# filter(sample_info, sample_name %in% colnames(test1)) |>
#   select(condition) |> pull() |> setdiff(y = colnames(test2))
# # hybrid
# dim(counts_allele)
# test_collapse <- collapseReplicates(sample_info_allele, counts_allele)
# dim(cbind(test_collapse$cer, test_collapse$par))
# # HAP4 in hyc, should not have changed
# test1 <- counts_allele[,sample_info_allele$experiment == "HAP4" &
#                   sample_info_allele$allele == "cer"]
# dim(test1)
# test2 <- test_collapse$cer[,grepl("HAP4", colnames(test_collapse$cer))]
# dim(test2)
# filter(sample_info_allele, sample_name %in% colnames(test1)) |>
#   select(condition) |> pull() |> setdiff(y = colnames(test2))
#
# # checking for genes with significant reduction in variance post-collapse
# # additional tests for collapseReplicates
# # toy count matrix using random sample of genes (to test ability to separate low var genes out)
# toy_mat <- counts[, sample_info$organism == "par" & sample_info$experiment == "LowN"]
# toydf <- makeDf(toy_mat, sample_info)
# # assessing distribution of variance, before collapsing replicates
# plotdf0 <- toydf |> group_by(gene_name) |>
#   summarise(var_expr = var(expr),
#             mean_expr = mean(expr))
# p0 <- ggplot(plotdf0, aes(x = log2(var_expr/mean_expr))) + geom_density() +
#   geom_vline(xintercept = 0, color = "red") +
#   ggtitle("before collapsing replicates")
# # after collapsing replicates
# toy_mat_collapsed <- collapsed$par[, info$experiment == "LowN"]
# plotdf1 <- tibble(var_expr = apply(toy_mat_collapsed, 1, var),
#                   mean_expr = apply(toy_mat_collapsed, 1, mean),
#                   gene_name = rownames(toy_mat_collapsed))
# p1 <- ggplot(plotdf1, aes(x = log2(var_expr/mean_expr))) + geom_density() +
#   geom_vline(xintercept = 0, color = "red") +
#   ggtitle("after collapsing replicates")
# ggarrange(p0, p1, nrow = 1, ncol = 2)
# # mean expression pre and post (shouldn't have changed)
# plotdf <- left_join(plotdf0, plotdf1, by = "gene_name", suffix = c("_pre", "_post"))
# ggplot(plotdf, aes(x = log2(mean_expr_pre), y = log2(mean_expr_post))) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0, color = "gold")
# plotdf <- select(plotdf, var_expr_pre, var_expr_post, mean_expr_pre, gene_name) |>
#   rename("mean_expr"="mean_expr_pre")
# # difference in var/mean pre and post versus mean expr
# ggplot(plotdf, aes(x = log2(mean_expr + 1),
#                    y = log2(var_expr_pre) - log2(var_expr_post))) +
#   geom_hex() +
#   geom_hline(yintercept = 0, color = "gold")
# # as expected, there is lower variance post collapse (very few genes below y = 0)
# # most genes have fairly little difference in var when collapsing replicates,
# # but the difference can be extreme for a few genes, of middling expression level
# # For genes with a difference, do we trust the collapsed variance to be more accurate to expression shape?
# # first random gene that tend to have very little difference in variance
# gene_idx <- plotdf |> select(gene_name) |> pull() |> sample(1)
# ggplot(filter(toydf, gene_name == gene_idx), aes(x = time_point_str, y = log2(expr + 1))) +
#   geom_jitter() +
#   geom_line(data = summarise(group_by(filter(toydf, gene_name == gene_idx), time_point_str, gene_name), mean_expr = mean(expr)),
#             aes(x = time_point_str, y = log2(mean_expr + 1), group = gene_name),
#             color = "red")
# # second a gene with variance substantially reduced in collapsed form
# gene_idx <- plotdf |> filter(log2(var_expr_pre) - log2(var_expr_post) > 7) |> select(gene_name) |> pull() |> sample(1)
# ggplot(filter(toydf, gene_name == gene_idx), aes(x = time_point_str, y = log2(expr + 1))) +
#   geom_jitter() +
#   geom_line(data = summarise(group_by(filter(toydf, gene_name == gene_idx), time_point_str, gene_name), mean_expr = mean(expr)),
#             aes(x = time_point_str, y = log2(mean_expr + 1), group = gene_name),
#             color = "red") # these are genes with very little expression change attributable to timepoint,
# # the exact kind we want to have significantly reduced variance in collapsed counts

#### Taking moving averages of Low Phosphorus and Diauxic Shift ####

# rationale: Low Pi and Diauxic Shift experiments needs to have expression
# smoothed as moving average b/c there are not replicates

# taking moving average
# window size of 5 --- 2 on each side, fewer for edge cases
# note: this function expects counts to have columns in order
# from smallest timepoint (on the left) to largest timepoint (on the right)
getMovingAverage <- function(.cts) {
  cts_movavg <- map(.x = colnames(.cts), \(cond) {
    idx <- which(colnames(.cts) == cond)
    idxs <- c(idx - 2,
              idx - 1,
              idx,
              idx + 1,
              idx + 2)[c(idx - 2,
                         idx - 1,
                         idx,
                         idx + 1,
                         idx + 2) > 0 &
                         c(idx - 2,
                           idx - 1,
                           idx,
                           idx + 1,
                           idx + 2) < ncol(.cts)]
    if (ncol(.cts[, idxs, drop = FALSE]) < 2) {
      warning("only one timepoint for", cond, idx, "\n")
      return(.cts[, idxs])
    }
    return(rowMeans(.cts[, idxs, drop = FALSE]))
  }) |> reduce(.f = cbind)
  colnames(cts_movavg) <- colnames(.cts)
  rownames(cts_movavg) <- rownames(.cts)
  return(cts_movavg)
}

# # tests for getMovingAverage
# # random gene prior to moving average
# gene_idx <- sample(rownames(counts), 1)
# # change to try different experiments (you have to keep the same experiment/allele/species):
# test_collapsed <- collapsed$par
# test_info <- info
# test_experiment <- "LowPi"
# test_movavg <- getMovingAverage(test_collapsed[,test_info$experiment == test_experiment])
# 
# plotdf <- tibble(expr = test_collapsed[gene_idx, test_info$experiment == test_experiment],
#                  condition = colnames(test_collapsed[, test_info$experiment == test_experiment])) |>
#   left_join(test_info, by = "condition")
# ggplot(plotdf, aes(x = time_point_num, y = expr)) + geom_line()
# 
# # same random gene after moving average
# plotdf2 <- tibble(expr = test_movavg[gene_idx,],
#                   status = "after",
#                   condition = colnames(test_movavg)) |>
#   left_join(test_info, by = "condition")
# plotdf$status <- "before"
# plotdf <- bind_rows(plotdf, plotdf2)
# ggplot(plotdf, aes(x = time_point_num, y = expr)) +
#   geom_line(aes(group = status, color = status))

#### Plotting ####

# plotting function to visualize expression profiles of any 2 groups
# @input: counts, info, and names of two groups to compare
# @output: ggplot of both expression profiles with loess or line curves tracing the average expression
plotExpressionProfilePair <- function(.cts1, .cts2, 
                                      .info1, .info2, 
                                      .name1 = "S. cerevisiae", .name2 = "S. paradoxus",
                                      .color1 = "orange1", .color2 = "blue2",
                                      .method = "line", 
                                      .show_points = FALSE,
                                      .point_size = 0.1,
                                      .show_confidence_intervals = TRUE,
                                      .confidence_type = "mean",
                                      .legend = "right",
                                      .normalization = c("none", "log2", "scale", "center"),
                                      .plotlims = NULL,
                                      .plot_titles = "experiment") {
  if (.normalization == "none") {
    norm_func <- identity
    ylabel <- "Expression (counts per million)"
  }
  if (.normalization == "log2") {
    norm_func <- \(x) {log2(x + 1)}
    ylabel <- "Expression (log2)"
  }
  if (.normalization == "scale") {
    norm_func <- \(x) {t(scale(t(x)))}
    ylabel <- "Expression (centered and scaled)"
  }
  if (.normalization == "center") {
    norm_func <- \(x) {(x - rowMeans(x, na.rm = TRUE))}
    ylabel <- "Expression (centered counts per million)"
  }
  if (.normalization == "centered log2") {
    norm_func <- \(x) {(log2(x + 1) - rowMeans(log2(x + 1), na.rm = TRUE))}
    ylabel <- "Expression\n(centered log2)"
  }
  if (!setequal(unique(.info1$experiment), unique(.info2$experiment))) {
    stop("sample info dataframes do not contain same set of experiments\n")
  }
  ExperimentNames <- unique(.info1$experiment) # arbitrary to do info1 or info2
  nExperiments <- length(ExperimentNames)
  nGenes <- nrow(.cts1)
  info1 <- tibble(experiment = .info1$experiment,
                  time_point_num = .info1$time_point_num)
  info2 <- tibble(experiment = .info2$experiment,
                  time_point_num = .info2$time_point_num)
  expr1 <- norm_func(.cts1) |> t()
  colnames(expr1) <- rownames(.cts1)
  expr2 <- norm_func(.cts2) |> t()
  colnames(expr2) <- rownames(.cts2)
  gdf1 <- bind_cols(expr1, info1) |> 
    pivot_longer(cols = colnames(expr1), names_to = "gene_name", values_to = "expr")
  gdf1$group_id <- "1"
  gdf2 <- bind_cols(expr2, info2) |> 
    pivot_longer(cols = colnames(expr2), names_to = "gene_name", values_to = "expr")
  gdf2$group_id <- "2"
  # converting each gene's expression to its mean expression between replicates
  gdf <- bind_rows(gdf1, gdf2) |> 
    drop_na() |> # drops genes missing from an experiment (usually Heat/Cold)
    group_by(group_id, gene_name, experiment, time_point_num) |> 
    summarise(expr = mean(expr)) |> ungroup()
  plotdf <- gdf
  # creating consistent plotlims across all experiments
  max_expr <- max(gdf$expr, na.rm = TRUE)
  min_expr <- min(gdf$expr, na.rm = TRUE)
  plotlimdf <- gdf |> group_by(time_point_num, experiment, group_id) |>
    summarise(mean_expr = mean(expr),
              sd_expr = sd(expr)) 
  max_avg <- plotlimdf |> select(mean_expr) |>
    pull() |> max(na.rm = TRUE)
  min_avg <- plotlimdf |> select(mean_expr) |>
    pull() |> min(na.rm = TRUE)
  buffer <- plotlimdf |> select(sd_expr) |>
    pull() |> max(na.rm = TRUE)
  buffer <- 0.25
  max_avg <- max_avg + buffer
  min_avg <- min_avg - buffer
  if (is.null(.plotlims)) {
    if (!.show_points) {
      .plotlims <- c(min_avg, max_avg)
    }
    if (.show_points) {
      .plotlims <- c(min_expr, max_expr)
    }
  }
  # background color rectangles for differentiating the experiments
  rects <- data.frame(color = c("orchid", "lightgreen", "gold", "orange", "salmon", "lightblue"),
                  labels = c("Cell Cycle", "Diauxic Shift", "Low Nitrogen", "Low Phosphorus", "Heat Stress", "Cold Stress"),
                  experiment_names = c("CC", "HAP4", "LowN", "LowPi", "Heat", "Cold"))
  experiment_order <- c("HAP4", "CC", "LowN", "LowPi", "Heat", "Cold")
  # plotting
  plotlist <- vector(mode = "list", length = length(unique(plotdf$experiment)))
  names(plotlist) <- experiment_order[experiment_order %in% unique(plotdf$experiment)]
  for (e in unique(plotdf$experiment)) {
    plotdf_e <- filter(plotdf, experiment == e)
    p <- ggplot() + 
      theme_classic() +
      scale_color_discrete(type = c(.color1, .color2), labels = c(.name1, .name2)) +
      theme(legend.title = element_blank()) +
      # theme(panel.background = element_rect(fill = alpha(rects$color[rects$experiment_names == e], 0.3),
      #                                       color = alpha(rects$color[rects$experiment_names == e], 0.3),
      #                                       size = 0.5, linetype = "solid")) +
      ylab("") +
      xlab("") +
      ylim(.plotlims)
      # scale_y_continuous(breaks = seq(from = 0, to = ceiling(max_expr), by = 1),
      #                    limits = seq(from = 0, to = ceiling(max_expr), by = 1),
      #                    labels = seq(from = 0, to = ceiling(max_expr), by = 1))
    if (.plot_titles != "none" & .plot_titles == "experiment") {
      p <- p + ggtitle(rects$labels[rects$experiment_names == e])
    }
    if (.plot_titles != "none" & .plot_titles == "ngenes") {
      p <- p + ggtitle(paste(nGenes, "genes"))
    }
    if (.plot_titles != "none" & .plot_titles != "experiment" &
        .plot_titles != "ngenes") {
      p <- p + ggtitle(.plot_titles)
    }
    if (.show_points) {
      p <- p + geom_jitter(data = plotdf_e, aes(x = time_point_num, 
                                                y = expr, color = group_id), 
                           size = .point_size, alpha = 0.5)
    }
    if (.method == "loess") {
      loess1 <- loess.sd(x = plotdf_e %>% filter(group_id == 1) %>% select(time_point_num) %>% pull(),
                         y = plotdf_e %>% filter(group_id == 1) %>% select(expr) %>% pull(), nsigma = 1.96)
      loess2 <- loess.sd(x = plotdf_e %>% filter(group_id == 2) %>% select(time_point_num) %>% pull(),
                         y = plotdf_e %>% filter(group_id == 2) %>% select(expr) %>% pull(), nsigma = 1.96)
      # adding loess segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
      p <- p + 
        geom_smooth(aes(x =!!loess1$x, y =!!loess1$y), color = .color1, linewidth = 1) +
        geom_smooth(aes(x =!!loess2$x, y =!!loess2$y), color = .color2, linewidth = 1)
      if (.show_confidence_intervals) {
        p <- p + 
          geom_ribbon(aes(x =!!loess1$x, ymin =!!loess1$lower, ymax =!!loess1$upper), fill = .color1, alpha = 0.3) +
          geom_ribbon(aes(x =!!loess2$x, ymin =!!loess2$lower, ymax =!!loess2$upper), fill = .color2, alpha = 0.3)
      }
    }
    if (.method == "line") {
      # lines trace average expr at each timepoint/experiment for each group
      avgexpr1 <- plotdf_e %>% filter(group_id == 1) |> group_by(time_point_num) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                              sd_expr = sd(expr, na.rm = TRUE),
                                                                                              quant975 = quantile(expr, 0.975, na.rm = TRUE),
                                                                                              quant025 = quantile(expr, 0.025, na.rm = TRUE))
      avgexpr2 <- plotdf_e %>% filter(group_id == 2) |> group_by(time_point_num) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                              sd_expr = sd(expr, na.rm = TRUE),
                                                                                              quant975 = quantile(expr, 0.975, na.rm = TRUE),
                                                                                              quant025 = quantile(expr, 0.025, na.rm = TRUE))
      # adding line segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
      p <- p +
        geom_line(data = avgexpr1, aes(x = time_point_num, y = mean_expr), color = .color1, linewidth = 1) +
        geom_line(data = avgexpr2, aes(x = time_point_num, y = mean_expr), color = .color2, linewidth = 1) 
      if (.show_confidence_intervals) {
        if (.confidence_type == "mean") {
          # calculating 95% confidence in the mean interval
          avgexpr1$CI_upper <- 1.96*(avgexpr1$sd_expr/sqrt(nGenes))
          avgexpr2$CI_upper <- 1.96*(avgexpr2$sd_expr/sqrt(nGenes))
          avgexpr1$CI_lower <- 1.96*(avgexpr1$sd_expr/sqrt(nGenes))
          avgexpr2$CI_lower <- 1.96*(avgexpr2$sd_expr/sqrt(nGenes))
        }
        if (.confidence_type == "all") {
          # calculating bounds for 95% of genes
          avgexpr1$CI_upper <- abs(avgexpr1$quant975) - avgexpr1$mean_expr
          avgexpr2$CI_upper <- abs(avgexpr2$quant975) - avgexpr1$mean_expr
          avgexpr1$CI_lower <- abs(avgexpr1$quant025) + avgexpr1$mean_expr
          avgexpr2$CI_lower <- abs(avgexpr2$quant025) + avgexpr1$mean_expr
        }
        p <- p + 
          geom_ribbon(data = avgexpr1, aes(x = time_point_num, ymin = pmax(mean_expr - CI_lower, .plotlims[1]), ymax = pmin(mean_expr + CI_upper, .plotlims[2])),
                      fill = .color1, alpha = 0.3) +
          geom_ribbon(data = avgexpr2, aes(x = time_point_num, ymin = pmax(mean_expr - CI_lower, .plotlims[1]), ymax = pmin(mean_expr + CI_upper, .plotlims[2])),
                      fill = .color2, alpha = 0.3)
      }
    }
    plotlist[[e]] <- p
  }
  if (length(plotlist) == 1) {
    return(plotlist[[1]] + xlab("Timepoint (min)") + ylab(ylabel))
  }
  fullplot <- ggarrange(plotlist = plotlist, nrow = 1, ncol = length(unique(plotdf$experiment)),
                        common.legend = TRUE, legend = .legend)
  return(annotate_figure(fullplot, bottom = "Timepoint (min)", left = ylabel))
}
# # tests for plotExpressionProfilePair
# gene_idxs <- finaldf |> filter(experiment == "LowPi" & dynamics == "diverged" &
#                                  cer == 1 & par == 2) |>
#   select(gene_name) |> pull()
# plotExpressionProfilePair(collapsed$cer[gene_idxs,],
#                           collapsed$par[gene_idxs,],
#                           info,
#                           info,
#                           .method = "line", .show_points = FALSE,
#                           .normalization = "scale",
#                           .confidence_type = "mean")
# # 2-1 cluster dynamics-divergers in HAP4 and LowPi
# gene_idxs <- finaldf |> filter(experiment == "HAP4" & cer == 2 & par == 1) |>
#   select(gene_name) |> pull()
# plotExpressionProfilePair(collapsed$cer[gene_idxs, info$experiment == "HAP4", drop = FALSE],
#                           collapsed$par[gene_idxs, info$experiment == "HAP4", drop = FALSE],
#                           info[info$experiment == "HAP4",],
#                           info[info$experiment == "HAP4",],
#                           .method = "line", .show_points = FALSE,
#                           .normalization = "centered log2")
# plotExpressionProfilePair(collapsed$cer[gene_idxs, info$experiment == "LowPi"],
#                           collapsed$par[gene_idxs, info$experiment == "LowPi"],
#                           info[info$experiment == "LowPi",],
#                           info[info$experiment == "LowPi",],
#                           .method = "line", .show_points = FALSE,
#                           .normalization = "centered log2")
# # unclear why this would ever come up, but this makes sure that the
# # order in which the experiments appear in the dataset doesn't
# # affect what they're called in the plot
# plotExpressionProfilePair(collapsed$cer[gene_idxs, info$experiment %in% c("CC", "HAP4")],
#                           collapsed$par[gene_idxs, info$experiment  %in% c("CC", "HAP4")],
#                           info[info$experiment %in% c("CC", "HAP4"),],
#                           info[info$experiment %in% c("CC", "HAP4"),],
#                           .name1 = "S. cereviaise",
#                           .name2 = "S. paradoxus",
#                           .method = "line", .show_points = TRUE,
#                           .normalization = "log2")
# # what happens when we force HAP4 to come first in the dataset?
# plotExpressionProfilePair(cbind(collapsed$cer[gene_idxs, info$experiment == "HAP4"],
#                                 collapsed$cer[gene_idxs, info$experiment == "CC"]),
#                           cbind(collapsed$par[gene_idxs, info$experiment == "HAP4"],
#                                 collapsed$par[gene_idxs, info$experiment == "CC"]),
#                           bind_rows(info[info$experiment == "HAP4",],
#                                     info[info$experiment == "CC",]),
#                           bind_rows(info[info$experiment == "HAP4",],
#                                     info[info$experiment == "CC",]),
#                           .method = "line", .show_points = TRUE,
#                           .normalization = "log2")
# # it doesn't do anything. Set experiment order inside function. Sat Growth 
# # should still have the X for dynamics divergers
#
# # SHU1 and SHU2, just for curiosity
# # They're in a complex together
# plotExpressionProfilePair(collapsed$cer["YHL006C",, drop = FALSE],
#                           collapsed$cer["YDR078C",, drop = FALSE],
#                           info,
#                           info,
#                           .method = "line",
#                           .show_points = TRUE,
#                           .show_confidence_intervals = TRUE,
#                           .normalization = "log2",
#                           .name1 = "SHU1 cer",
#                           .name2 = "SHU2 cer",
#                           .color1 = "purple",
#                           .color2 = "red")
# plotExpressionProfilePair(collapsed$par["YHL006C",, drop = FALSE],
#                           collapsed$par["YDR078C",, drop = FALSE],
#                           info,
#                           info,
#                           .method = "line",
#                           .show_points = TRUE,
#                           .show_confidence_intervals = TRUE,
#                           .normalization = "log2",
#                           .name1 = "SHU1 par",
#                           .name2 = "SHU2 par",
#                           .color1 = "purple",
#                           .color2 = "red")

# plotting function to visualize expression profiles of any 2 groups
# @input: counts, info, and names of two groups to compare
# @output: ggplot of both expression profiles with loess or line curves tracing the average expression
plotExpressionRibbonsPair <- function(.cts1, .cts2, 
                                      .info1, .info2, 
                                      .name1 = "S. cerevisiae", .name2 = "S. paradoxus",
                                      .color1 = "orange1", .color2 = "blue2",
                                      .alpha = 0.2,
                                      .legend = "right",
                                      .normalization = c("none", "log2", "scale", "center"),
                                      .plotlims = NULL,
                                      .plot_titles = "experiment") {
  if (.normalization == "none") {
    norm_func <- identity
    ylabel <- "Expression (counts per million)"
  }
  if (.normalization == "log2") {
    norm_func <- \(x) {log2(x + 1)}
    ylabel <- "Expression (log2)"
  }
  if (.normalization == "scale") {
    norm_func <- \(x) {t(scale(t(x)))}
    ylabel <- "Expression (centered and scaled)"
  }
  if (.normalization == "center") {
    norm_func <- \(x) {(x - rowMeans(x, na.rm = TRUE))}
    ylabel <- "Expression (centered counts per million)"
  }
  if (.normalization == "centered log2") {
    norm_func <- \(x) {(log2(x + 1) - rowMeans(log2(x + 1), na.rm = TRUE))}
    ylabel <- "Expression\n(centered log2)"
  }
  if (!setequal(unique(.info1$experiment), unique(.info2$experiment))) {
    stop("sample info dataframes do not contain same set of experiments\n")
  }
  ExperimentNames <- unique(.info1$experiment) # arbitrary to do info1 or info2
  nExperiments <- length(ExperimentNames)
  nGenes <- nrow(.cts1)
  info1 <- tibble(experiment = .info1$experiment,
                  time_point_num = .info1$time_point_num)
  info2 <- tibble(experiment = .info2$experiment,
                  time_point_num = .info2$time_point_num)
  expr1 <- norm_func(.cts1) |> t()
  colnames(expr1) <- rownames(.cts1)
  expr2 <- norm_func(.cts2) |> t()
  colnames(expr2) <- rownames(.cts2)
  gdf1 <- bind_cols(expr1, info1) |> 
    pivot_longer(cols = colnames(expr1), names_to = "gene_name", values_to = "expr")
  gdf1$group_id <- "1"
  gdf2 <- bind_cols(expr2, info2) |> 
    pivot_longer(cols = colnames(expr2), names_to = "gene_name", values_to = "expr")
  gdf2$group_id <- "2"
  # converting each gene's expression to its mean expression between replicates
  gdf <- bind_rows(gdf1, gdf2) |> 
    drop_na() |> # drops genes missing from an experiment (usually Heat/Cold)
    group_by(group_id, gene_name, experiment, time_point_num) |> 
    summarise(expr = mean(expr)) |> ungroup()
  plotdf <- gdf
  # creating consistent plotlims across all experiments
  max_expr <- max(gdf$expr, na.rm = TRUE)
  min_expr <- min(gdf$expr, na.rm = TRUE)
  plotlimdf <- gdf |> group_by(time_point_num, experiment, group_id) |>
    summarise(mean_expr = mean(expr),
              sd_expr = sd(expr)) 
  max_avg <- plotlimdf |> select(mean_expr) |>
    pull() |> max(na.rm = TRUE)
  min_avg <- plotlimdf |> select(mean_expr) |>
    pull() |> min(na.rm = TRUE)
  buffer <- plotlimdf |> select(sd_expr) |>
    pull() |> max(na.rm = TRUE)
  buffer <- 0.25
  max_avg <- max_avg + buffer
  min_avg <- min_avg - buffer
  if (is.null(.plotlims)) {
    .plotlims <- c(min_avg, max_avg)
  }
  experiment_order <- c("HAP4", "CC", "LowN", "LowPi", "Heat", "Cold")
  # background color rectangles for differentiating the experiments
  rects <- data.frame(color = c("orchid", "lightgreen", "gold", "orange", "salmon", "lightblue"),
                      labels = c("Cell Cycle", "Diauxic Shift", "Low Nitrogen", "Low Phosphorus", "Heat Stress", "Cold Stress"),
                      experiment_names = c("CC", "HAP4", "LowN", "LowPi", "Heat", "Cold"))
  # plotting
  plotlist <- vector(mode = "list", length = length(unique(plotdf$experiment)))
  names(plotlist) <- experiment_order[experiment_order %in% unique(plotdf$experiment)]
  for (e in unique(plotdf$experiment)) {
    plotdf_e <- filter(plotdf, experiment == e)
    p <- ggplot() + 
      theme_classic() +
      scale_color_discrete(type = c(.color1, .color2), labels = c(.name1, .name2)) +
      theme(legend.title = element_blank()) +
      # theme(panel.background = element_rect(fill = alpha(rects$color[rects$experiment_names == e], 0.3),
      #                                       color = alpha(rects$color[rects$experiment_names == e], 0.3),
      #                                       size = 0.5, linetype = "solid")) +
      ylab("") +
      xlab("") +
      ylim(.plotlims)
    # scale_y_continuous(breaks = seq(from = 0, to = ceiling(max_expr), by = 1),
    #                    limits = seq(from = 0, to = ceiling(max_expr), by = 1),
    #                    labels = seq(from = 0, to = ceiling(max_expr), by = 1))
    if (.plot_titles != "none" & .plot_titles == "experiment") {
      p <- p + ggtitle(rects$labels[rects$experiment_names == e])
    }
    if (.plot_titles != "none" & .plot_titles == "ngenes") {
      p <- p + ggtitle(paste(nGenes, "genes"))
    }
    if (.plot_titles != "none" & .plot_titles != "experiment" &
        .plot_titles != "ngenes") {
      p <- p + ggtitle(.plot_titles[rects$experiment_names == e])
    }
    # collecting ribbons
    # lines trace average expr at each timepoint/experiment for each group
    quant_vec1 <- plotdf_e |> filter(group_id == 1) |> group_by(time_point_num) |> summarise(max_expr = max(expr, na.rm = TRUE),
                                                                                             quant90 = quantile(expr, 0.9, na.rm = TRUE),
                                                                                             quant80 = quantile(expr, 0.8, na.rm = TRUE),
                                                                                             quant60 = quantile(expr, 0.6, na.rm = TRUE),
                                                                                             quant55 = quantile(expr, 0.55, na.rm = TRUE),
                                                                                             quant45 = quantile(expr, 0.45, na.rm = TRUE),
                                                                                             quant40 = quantile(expr, 0.4, na.rm = TRUE),
                                                                                             quant20 = quantile(expr, 0.2, na.rm = TRUE),
                                                                                             quant10 = quantile(expr, 0.1, na.rm = TRUE),
                                                                                             min_expr = min(expr, na.rm = TRUE))
    quant_vec2 <- plotdf_e |> filter(group_id == 2) |> group_by(time_point_num) |> summarise(max_expr = max(expr, na.rm = TRUE),
                                                                                             quant90 = quantile(expr, 0.9, na.rm = TRUE),
                                                                                             quant80 = quantile(expr, 0.8, na.rm = TRUE),
                                                                                             quant60 = quantile(expr, 0.6, na.rm = TRUE),
                                                                                             quant55 = quantile(expr, 0.55, na.rm = TRUE),
                                                                                             quant45 = quantile(expr, 0.45, na.rm = TRUE),
                                                                                             quant40 = quantile(expr, 0.4, na.rm = TRUE),
                                                                                             quant20 = quantile(expr, 0.2, na.rm = TRUE),
                                                                                             quant10 = quantile(expr, 0.1, na.rm = TRUE),
                                                                                             min_expr = min(expr, na.rm = TRUE))
    # adding ribbons (the order we add them is the order they're arranged in the plot)
    p <- p +
      # # min and max, bounds of all genes
      #   geom_ribbon(data = quant_vec1, aes(x = time_point_num, 
      #                                      ymin = pmax(min_expr, .plotlims[1]), 
      #                                      ymax = pmin(max_expr, .plotlims[2])),
      #               fill = .color1, alpha = .alpha) +
      #   geom_ribbon(data = quant_vec2, aes(x = time_point_num, 
      #                                      ymin = pmax(min_expr, .plotlims[1]), 
      #                                      ymax = pmin(max_expr, .plotlims[2])),
      #               fill = .color2, alpha = .alpha) +
      # 90% and 10%, bounds of 80% of all genes
      geom_ribbon(data = quant_vec1, aes(x = time_point_num, 
                                         ymin = pmax(quant10, .plotlims[1]), 
                                         ymax = pmin(quant90, .plotlims[2])),
                  fill = .color1, alpha = .alpha) +
      geom_ribbon(data = quant_vec2, aes(x = time_point_num, 
                                         ymin = pmax(quant10, .plotlims[1]), 
                                         ymax = pmin(quant90, .plotlims[2])),
                  fill = .color2, alpha = .alpha) +
      # 80% and 20%, bounds of 60% of all genes
      geom_ribbon(data = quant_vec1, aes(x = time_point_num, 
                                         ymin = pmax(quant20, .plotlims[1]), 
                                         ymax = pmin(quant80, .plotlims[2])),
                  fill = .color1, alpha = .alpha) +
      geom_ribbon(data = quant_vec2, aes(x = time_point_num, 
                                         ymin = pmax(quant20, .plotlims[1]), 
                                         ymax = pmin(quant80, .plotlims[2])),
                  fill = .color2, alpha = .alpha) +
      # 60% and 40%, bounds of 20% of all genes
      geom_ribbon(data = quant_vec1, aes(x = time_point_num, 
                                         ymin = pmax(quant40, .plotlims[1]), 
                                         ymax = pmin(quant60, .plotlims[2])),
                  fill = .color1, alpha = .alpha) +
      geom_ribbon(data = quant_vec2, aes(x = time_point_num, 
                                         ymin = pmax(quant40, .plotlims[1]), 
                                         ymax = pmin(quant60, .plotlims[2])),
                  fill = .color2, alpha = .alpha) +
    # 55% and 45%, bounds of 10% of all genes
    geom_ribbon(data = quant_vec1, aes(x = time_point_num, 
                                       ymin = pmax(quant45, .plotlims[1]), 
                                       ymax = pmin(quant55, .plotlims[2])),
                fill = .color1, alpha = .alpha) +
      geom_ribbon(data = quant_vec2, aes(x = time_point_num, 
                                         ymin = pmax(quant45, .plotlims[1]), 
                                         ymax = pmin(quant55, .plotlims[2])),
                  fill = .color2, alpha = .alpha)
    plotlist[[e]] <- p
  }
  if (length(plotlist) == 1) {
    return(plotlist[[1]] + xlab("Timepoint (min)") + ylab(ylabel))
  }
  fullplot <- ggarrange(plotlist = plotlist, nrow = 1, ncol = length(unique(plotdf$experiment)),
                        common.legend = TRUE, legend = .legend)
  return(annotate_figure(fullplot, bottom = "Timepoint (min)", left = ylabel))
}
# # tests for plotExpressionProfilePair
# gene_idxs <- finaldf |> filter(experiment == "LowPi" & dynamics == "diverged" &
#                                  cer == 1 & par == 2) |>
#   select(gene_name) |> pull()
# plotExpressionRibbonsPair(collapsed$cer[gene_idxs,],
#                           collapsed$par[gene_idxs,],
#                           info,
#                           info,
#                           .color1 = "orange1",
#                           .color2 = "blue2",
#                           .normalization = "scale",
#                           .plotlims = c(-2.5, 2.5))

# Oh yes
# plots 4 groups of genes BUT you can't do this willy nilly
# for it to be interpretable, groups 1 and 2 are the main contrast
# and groups 3 and 4 are related to groups 1 and 2 respectively
# Typically 1 and 2 are the parental species, 3 is the hybrid allele of 1, 
# and 4 is the hybrid allele of 2
plotExpressionProfileQuartet <- function(.cts1, .cts2, .cts3, .cts4,
                                         .info1, .info2, .info3, .info4,
                                         .name1 = "S. cerevisiae",
                                         .name2 = "S. paradoxus",
                                         .name3 = "F1 hybrid, cerevisiae allele",
                                         .name4 = "F1 hybrid, paradoxus allele",
                                         .color1 = "orange1",
                                         .color2 = "blue2",
                                         .color3 = "orange4",
                                         .color4 = "blue4",
                                         .method = c("line", "loess"),
                                         .show_points = FALSE,
                                         .show_confidence_intervals = TRUE,
                                         .normalization = c("none", "log2", "scale", "center",
                                                            "centered log2"),
                                         .plotlims = NULL,
                                         .plot_titles = "experiment") {
  if (.normalization == "none") {
    norm_func <- identity
    ylabel <- "Expression (counts per million)"
  }
  if (.normalization == "log2") {
    norm_func <- \(x) {log2(x + 1)}
    ylabel <- "Expression (log2)"
  }
  if (.normalization == "scale") {
    norm_func <- \(x) {t(scale(t(x)))}
    ylabel <- "Expression (centered and scaled)"
  }
  if (.normalization == "center") {
    norm_func <- \(x) {
      return(x - rowMeans(x, na.rm = TRUE))}
    ylabel <- "Expression (centered counts per million)"
  }
  if (.normalization == "centered log2") {
    norm_func <- \(x) {(log2(x + 1) - rowMeans(log2(x + 1), na.rm = TRUE))}
    ylabel <- "Expression\n(centered log2)"
  }
  if (!setequal(unique(.info1$experiment), unique(.info2$experiment))) {
    stop("sample info dataframes do not contain same set of experiments\n")
  }
  ExperimentNames <- unique(.info1$experiment) # arbitrary which info to use
  nExperiments <- length(ExperimentNames)
  info1 <- tibble(experiment = .info1$experiment,
                  time_point_num = .info1$time_point_num)
  info2 <- tibble(experiment = .info2$experiment,
                  time_point_num = .info2$time_point_num)
  info3 <- tibble(experiment = .info3$experiment,
                  time_point_num = .info3$time_point_num)
  info4 <- tibble(experiment = .info4$experiment,
                  time_point_num = .info4$time_point_num)
  expr1 <- norm_func(.cts1) |> t()
  expr2 <- norm_func(.cts2) |> t()
  expr3 <- norm_func(.cts3) |> t()
  expr4 <- norm_func(.cts4) |> t()
  colnames(expr1) <- rownames(.cts1)
  colnames(expr2) <- rownames(.cts2)
  colnames(expr3) <- rownames(.cts3)
  colnames(expr4) <- rownames(.cts4)
  gdf1 <- bind_cols(expr1, info1) |> 
    pivot_longer(cols = colnames(expr1), names_to = "gene_name", values_to = "expr")
  gdf1$group_id <- "1"
  gdf2 <- bind_cols(expr2, info2) |> 
    pivot_longer(cols = colnames(expr2), names_to = "gene_name", values_to = "expr")
  gdf2$group_id <- "2"
  gdf3 <- bind_cols(expr3, info3) |> 
    pivot_longer(cols = colnames(expr3), names_to = "gene_name", values_to = "expr")
  gdf3$group_id <- "3"
  gdf4 <- bind_cols(expr4, info4) |> 
    pivot_longer(cols = colnames(expr4), names_to = "gene_name", values_to = "expr")
  gdf4$group_id <- "4"
  # converting each gene's expression to its mean expression between replicates
  gdf <- bind_rows(gdf1, gdf2, gdf3, gdf4) |> 
    drop_na() |> # drops genes missing from an experiment (usually Heat/Cold)
    group_by(group_id, gene_name, experiment, time_point_num) |>
    summarise(expr = mean(expr, na.rm = TRUE)) |> ungroup()
  plotlimdf <- gdf |> group_by(time_point_num, experiment, group_id) |>
    summarise(mean_expr = mean(expr, na.rm = TRUE),
              sd_expr = sd(expr, na.rm = TRUE)) 
  max_avg <- plotlimdf |> select(mean_expr) |>
    pull() |> max(na.rm = TRUE)
  min_avg <- plotlimdf |> select(mean_expr) |>
    pull() |> min(na.rm = TRUE)
  buffer <- 0.25
  max_avg <- max_avg + buffer
  min_avg <- min_avg - buffer
  # min/maxs when plotting points as well as averages:
  max_expr <- max(gdf$expr, na.rm = TRUE)
  min_expr <- min(gdf$expr, na.rm = TRUE)
  # setting plotlims if they haven't been manually set
  if (is.null(.plotlims)) {
    if (!.show_points) {
      .plotlims <- c(min_avg, max_avg)
    }
    if (.show_points) {
      .plotlims <- c(min_expr, max_expr)
    }
  }
  plotdf <- gdf
  # background color rectangles for differentiating the 4 experiments
  rects <- data.frame(color = c("lightgreen", "orchid", "gold", "orange", "salmon", "lightblue"),
                      labels = c("Diauxic Shift", "Cell Cycle", "Low Nitrogen", "Low Phosphorus", "Heat Stress", "Cold Stress"),
                      experiment_names = c("HAP4", "CC", "LowN", "LowPi", "Heat", "Cold"))
  # plotting
  plotlist <- vector(mode = "list", length = length(unique(plotdf$experiment)))
  names(plotlist) <- unique(plotdf$experiment)
  experiment_order <- intersect(rects$experiment_names, unique(plotdf$experiment))
  for (e in experiment_order) {
    plotdf_e <- filter(plotdf, experiment == e)
    p <- ggplot() + 
      theme_classic() +
      scale_color_discrete(type = c(.color1, .color2, .color3, .color4), 
                           labels = c(.name1, .name2, .name3, .name4)) +
      theme(legend.title = element_blank()) +
      # theme(panel.background = element_rect(fill = alpha(rects$color[rects$experiment_names == e], 0.3),
      #                                       color = alpha(rects$color[rects$experiment_names == e], 0.3),
      #                                       size = 0.5, linetype = "solid")) +
      ylab("") +
      xlab("") +
      ylim(.plotlims)
    if (.plot_titles != "none" & .plot_titles == "experiment") {
      p <- p + ggtitle(rects$labels[rects$experiment_names == e])
    }
    if (.plot_titles != "none" & .plot_titles == "ngenes") {
      p <- p + ggtitle(paste(nGenes, "genes"))
    }
    if (.plot_titles != "none" & .plot_titles != "experiment" &
        .plot_titles != "ngenes") {
      p <- p + ggtitle(.plot_titles)
    }
    if (.show_points) {
      p <- p + geom_jitter(data = plotdf_e, aes(x = time_point_num, y = expr, color = group_id), size = 0.1, alpha = 0.5)
    }
    if (.method == "loess") {
      loess1 <- loess.sd(x = plotdf_e %>% filter(group_id == 1) %>% select(time_point_num) %>% pull(),
                         y = plotdf_e %>% filter(group_id == 1) %>% select(expr) %>% pull(), nsigma = 1.96)
      loess2 <- loess.sd(x = plotdf_e %>% filter(group_id == 2) %>% select(time_point_num) %>% pull(),
                         y = plotdf_e %>% filter(group_id == 2) %>% select(expr) %>% pull(), nsigma = 1.96)
      loess3 <- loess.sd(x = plotdf_e %>% filter(group_id == 3) %>% select(time_point_num) %>% pull(),
                         y = plotdf_e %>% filter(group_id == 3) %>% select(expr) %>% pull(), nsigma = 1.96)
      loess4 <- loess.sd(x = plotdf_e %>% filter(group_id == 4) %>% select(time_point_num) %>% pull(),
                         y = plotdf_e %>% filter(group_id == 4) %>% select(expr) %>% pull(), nsigma = 1.96)
      # adding loess segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
      p <- p + 
        geom_smooth(aes(x =!!loess1$x, y =!!loess1$y), color = .color1, linewidth = 1) +
        geom_smooth(aes(x =!!loess2$x, y =!!loess2$y), color = .color2, linewidth = 1) +
        geom_smooth(aes(x =!!loess3$x, y =!!loess3$y), color = .color3, linewidth = 1) +
        geom_smooth(aes(x =!!loess4$x, y =!!loess4$y), color = .color4, linewidth = 1)
      if (.show_confidence_intervals) {
        p <- p + 
          geom_ribbon(aes(x =!!loess1$x, ymin =!!loess1$lower, ymax =!!loess1$upper), fill = .color1, alpha = 0.3) +
          geom_ribbon(aes(x =!!loess2$x, ymin =!!loess2$lower, ymax =!!loess2$upper), fill = .color2, alpha = 0.3)
      }
    }
    if (.method == "line") {
      # lines trace average expr at each timepoint/experiment for each group
      avgexpr1 <- plotdf_e %>% filter(group_id == 1) |> group_by(time_point_num) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                              sd_expr = sd(expr, na.rm = TRUE))
      avgexpr2 <- plotdf_e %>% filter(group_id == 2) |> group_by(time_point_num) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                              sd_expr = sd(expr, na.rm = TRUE))
      avgexpr3 <- plotdf_e %>% filter(group_id == 3) |> group_by(time_point_num) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                              sd_expr = sd(expr, na.rm = TRUE))
      avgexpr4 <- plotdf_e %>% filter(group_id == 4) |> group_by(time_point_num) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                              sd_expr = sd(expr, na.rm = TRUE))
      # adding line segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
      p <- p +
        geom_line(data = avgexpr1, aes(x = time_point_num, y = mean_expr), color = .color1, linewidth = 1, linetype = "solid") +
        geom_line(data = avgexpr2, aes(x = time_point_num, y = mean_expr), color = .color2, linewidth = 1, linetype = "solid") +
        geom_line(data = avgexpr3, aes(x = time_point_num, y = mean_expr), color = .color3, linewidth = 1, linetype = "dashed") +
        geom_line(data = avgexpr4, aes(x = time_point_num, y = mean_expr), color = .color4, linewidth = 1, linetype = "dashed")
      if (.show_confidence_intervals) {
        # calculating 95% confidence in the mean
        nGenes <- length(unique(plotdf$gene_name))
        avgexpr1$CI <- 1.96*(avgexpr1$sd_expr/sqrt(nGenes))
        avgexpr2$CI <- 1.96*(avgexpr2$sd_expr/sqrt(nGenes))
        avgexpr3$CI <- 1.96*(avgexpr3$sd_expr/sqrt(nGenes))
        avgexpr4$CI <- 1.96*(avgexpr4$sd_expr/sqrt(nGenes))
        
        p <- p + 
          geom_ribbon(data = avgexpr1, aes(x = time_point_num, ymin = pmax(mean_expr - CI, .plotlims[1]), ymax = pmin(mean_expr + CI, .plotlims[2])),
                      fill = .color1, alpha = 0.3) +
          geom_ribbon(data = avgexpr2, aes(x = time_point_num, ymin = pmax(mean_expr - CI, .plotlims[1]), ymax = pmin(mean_expr + CI, .plotlims[2])),
                      fill = .color2, alpha = 0.3) +
          geom_ribbon(data = avgexpr3, aes(x = time_point_num, ymin = pmax(mean_expr - CI, .plotlims[1]), ymax = pmin(mean_expr + CI, .plotlims[2])),
                      fill = .color1, alpha = 0.3) +
          geom_ribbon(data = avgexpr4, aes(x = time_point_num, ymin = pmax(mean_expr - CI, .plotlims[1]), ymax = pmin(mean_expr + CI, .plotlims[2])),
                      fill = .color2, alpha = 0.3)
      }
    }
    plotlist[[e]] <- p
  }
  if (length(plotlist) == 1) {
    return(plotlist[[1]] + xlab("Timepoint (min)") + ylab(ylabel))
  }
  fullplot <- ggarrange(plotlist = plotlist[experiment_order], nrow = 1, ncol = length(experiment_order),
                        common.legend = TRUE, legend = "right")
  return(annotate_figure(fullplot, bottom = "Timepoint (min)", left = ylabel))
}
# # tests for plotExpressionProfileQuartet
# # hardcoded for simplicity (subset of the unsigned module b (yellow) with positive and negatively correlated genes)
# conserved_idxs <- c("YKL013C", "YER009W", "YMR097C", "YJL189W", "YKL009W",
#                     "YEL054C", "YLR333C", "YBL050W", "YNL223W", "YNL162W")
# up_par_idxs <- c("YER102W", "YLR264W", "YMR304W", "YHR193C", "YEL034W",
#                  "YOR167C", "YBL072C", "YGL135W", "YDL191W", "YHR021C")
# up_cer_idxs <- c("YMR194C-B", "YHR161C", "YJL127C-B", "YDL027C", "YNL175C",
#                  "YHR104W", "YMR027W", "YDR479C", "YFR047C", "YJL055W")
# # first yellow cer vs par, with up_cer genes indicated
# test <- plotExpressionProfileQuartet(.cts1 = collapsed$cer[conserved_idxs,],
#                              .cts2 = collapsed$par[conserved_idxs,],
#                              .cts3 = collapsed$cer[up_cer_idxs,],
#                              .cts4 = collapsed$par[up_cer_idxs,],
#                              .info1 = info, .info2 = info, .info3 = info, .info4 = info,
#                              .method = "line",
#                              .show_points = TRUE,
#                              .show_confidence_intervals = FALSE,
#                              .normalization = "log2")
# # second yellow cer vs par, with up_par genes indicated
# plotExpressionProfileQuartet(.cts1 = collapsed$cer[conserved_idxs,],
#                              .cts2 = collapsed$par[conserved_idxs,],
#                              .cts3 = collapsed$cer[up_par_idxs,],
#                              .cts4 = collapsed$par[up_par_idxs,],
#                              .info1 = info, .info2 = info, .info3 = info, .info4 = info,
#                              .method = "line", .show_points = TRUE,
#                              .normalization = "log2")

# wrapper function for plotExpressionProfilePair/Quartet
# plot specific genes' expression in one single environment
plotGenes <- function(.gene_idxs,
                      .normalization = "log2",
                      .quartet = FALSE,
                      .plotlims = NULL,
                      .plot_titles = "none",
                      .collapsed = TRUE,
                      .experiment_name) {
  if (length(.experiment_name) != 1) {
    stop("plotGenes only accepts single environmental conditions\n")
  }
  if (.collapsed) {
    counts_cer <- collapsed$cer[.gene_idxs, info$experiment == .experiment_name, drop = FALSE]
    counts_par <- collapsed$par[.gene_idxs, info$experiment == .experiment_name, drop = FALSE]
    counts_hyc <- collapsed_allele$cer[.gene_idxs, info_allele$experiment == .experiment_name, drop = FALSE]
    counts_hyp <- collapsed_allele$par[.gene_idxs, info_allele$experiment == .experiment_name, drop = FALSE]
    info_cer <- info |> filter(experiment == .experiment_name)
    info_par <- info |> filter(experiment == .experiment_name)
    info_hyc <- info_allele |> filter(experiment == .experiment_name)
    info_hyp <- info_allele |> filter(experiment == .experiment_name)
  }
  if (!.collapsed) {
    counts_cer <- counts[.gene_idxs, sample_info$experiment == .experiment_name &
                           sample_info$allele == "cer", drop = FALSE]
    counts_par <- counts[.gene_idxs, sample_info$experiment == .experiment_name &
                           sample_info$allele == "par", drop = FALSE]
    counts_hyc <- counts_allele[.gene_idxs, sample_info_allele$experiment == .experiment_name &
                                  sample_info_allele$allele == "cer", drop = FALSE]
    counts_hyp <- counts_allele[.gene_idxs, sample_info_allele$experiment == .experiment_name &
                                  sample_info_allele$allele == "par", drop = FALSE]
    info_cer <- sample_info[sample_info$experiment == .experiment_name &
                              sample_info$allele == "cer",]
    info_par <- sample_info[sample_info$experiment == .experiment_name &
                              sample_info$allele == "par",]
    info_hyc <- sample_info_allele[sample_info_allele$experiment == .experiment_name &
                                     sample_info_allele$allele == "cer",]
    info_hyp <- sample_info_allele[sample_info_allele$experiment == .experiment_name &
                                     sample_info_allele$allele == "par",]
  }
  if (!.quartet) {
    p <- plotExpressionProfilePair(.cts1 = counts_cer,
                                   .cts2 = counts_par,
                                   .info1 = info_cer,
                                   .info2 = info_par,
                                   .name1 = "S. cer",
                                   .name2 = "S. par",
                                   .color1 = "orange1",
                                   .color2 = "blue2",
                                   .normalization = .normalization,
                                   .method = "line",
                                   .show_points = FALSE,
                                   .show_confidence_intervals = TRUE,
                                   .plotlims = .plotlims,
                                   .plot_titles = .plot_titles)
  }
  if (.quartet) {
    p <- plotExpressionProfileQuartet(.cts1 = counts_cer,
                                      .cts2 = counts_par,
                                      .cts3 = counts_hyc,
                                      .cts4 = counts_hyp,
                                      .info1 = info_cer,
                                      .info2 = info_par,
                                      .info3 = info_hyc,
                                      .info4 = info_hyp,
                                      .name1 = "S. cer",
                                      .name2 = "S. par",
                                      .name3 = "F1, cer allele",
                                      .name4 = "F1, par allele",
                                      .color1 = "orange1",
                                      .color2 = "blue2",
                                      .color3 = "orange4",
                                      .color4 = "blue4",
                                      .normalization = .normalization,
                                      .method = "line",
                                      .show_points = FALSE,
                                      .show_confidence_intervals = TRUE,
                                      .plotlims = .plotlims,
                                      .plot_titles = .plot_titles)
  }
  return(p)
}
# # tests for plotGenes
# # YIR041W (only remains in hybrid CC, filtered from all others by low expr)
# # but it does seem to be species-specific expression, and I want to see if that's
# # true across all experiments, just straddling the threshold of low expr
# test <- plotGenes("YIR041W", .experiment_name = "HAP4") # species specific
# # before adding drop_na to gdf, profiles didn't show up in Heat/Cold
# # b/c of few missing genes
# # Example: LowPi 1-2
# gene_idxs <- finaldf |> filter(experiment == "LowPi" & dynamics == "diverged" &
#                                  cer == 1 & par == 2) |>
#   select(gene_name) |> pull()
# plotGenes(gene_idxs, .experiment_name = "HAP4", .plot_titles = "HAP4")
# plotGenes(gene_idxs, .experiment_name = "LowPi", .plot_titles = "LowPi")
# plotGenes(gene_idxs, .experiment_name = "LowN", .plot_titles = "LowN")
# plotGenes(gene_idxs, .experiment_name = "Heat", .plot_titles = "Heat") # used to be missing
# plotGenes(gene_idxs, .experiment_name = "Cold", .plot_titles = "Cold") # used to be missing
# na_idxs <- which(is.na(collapsed$cer[gene_idxs, info$experiment == "Heat"]),
#                  arr.ind = TRUE)
# plotGenes(gene_idxs[-na_idxs[,1]], .experiment_name = "Heat")
# plotGenes(gene_idxs[-na_idxs[,1]], .experiment_name = "Cold")

# wrapper for plotExpressionProfile designed for all environments
plotEnvironments <- function(.gene_idxs,
                             .normalization = "log2",
                             .quartet = FALSE,
                             .plotlims = NULL,
                             .plot_titles = "none",
                             .collapsed = TRUE) {
  if (.collapsed) {
    counts_cer <- collapsed$cer[.gene_idxs,, drop = FALSE]
    counts_par <- collapsed$par[.gene_idxs,, drop = FALSE]
    counts_hyc <- collapsed_allele$cer[.gene_idxs,, drop = FALSE]
    counts_hyp <- collapsed_allele$par[.gene_idxs,, drop = FALSE]
    info_cer <- info
    info_par <- info
    info_hyc <- info_allele
    info_hyp <- info_allele
  }
  if (!.collapsed) {
    counts_cer <- counts[.gene_idxs, sample_info$allele == "cer", drop = FALSE]
    counts_par <- counts[.gene_idxs, sample_info$allele == "par", drop = FALSE]
    counts_hyc <- counts_allele[.gene_idxs, sample_info_allele$allele == "cer", drop = FALSE]
    counts_hyp <- counts_allele[.gene_idxs, sample_info_allele$allele == "par", drop = FALSE]
    info_cer <- sample_info[sample_info$allele == "cer",]
    info_par <- sample_info[sample_info$allele == "par",]
    info_hyc <- sample_info_allele[sample_info_allele$allele == "cer",]
    info_hyp <- sample_info_allele[sample_info_allele$allele == "par",]
  }
  if (!.quartet) {
    p <- plotExpressionProfilePair(.cts1 = counts_cer,
                                   .cts2 = counts_par,
                                   .info1 = info_cer,
                                   .info2 = info_par,
                                   .name1 = "S. cer",
                                   .name2 = "S. par",
                                   .color1 = "orange1",
                                   .color2 = "blue2",
                                   .normalization = .normalization,
                                   .method = "line",
                                   .show_points = FALSE,
                                   .show_confidence_intervals = TRUE,
                                   .plotlims = .plotlims,
                                   .plot_titles = .plot_titles)
  }
  if (.quartet) {
    p <- plotExpressionProfileQuartet(.cts1 = counts_cer,
                                      .cts2 = counts_par,
                                      .cts3 = counts_hyc,
                                      .cts4 = counts_hyp,
                                      .info1 = info_cer,
                                      .info2 = info_par,
                                      .info3 = info_hyc,
                                      .info4 = info_hyp,
                                      .name1 = "S. cer",
                                      .name2 = "S. par",
                                      .name3 = "F1, cer allele",
                                      .name4 = "F1, par allele",
                                      .color1 = "orange1",
                                      .color2 = "blue2",
                                      .color3 = "orange4",
                                      .color4 = "blue4",
                                      .normalization = .normalization,
                                      .method = "line",
                                      .show_points = FALSE,
                                      .show_confidence_intervals = TRUE,
                                      .plotlims = .plotlims,
                                      .plot_titles = .plot_titles)
  }
  return(p)
}
# # tests for plotEnvironments
# gene_idxs <- finaldf |> filter(experiment == "LowPi" & dynamics == "diverged" &
#                                  cer == 1 & par == 2) |>
#   select(gene_name) |> pull()
# plotEnvironments(gene_idxs)

# plots WT vs TFdel in both species
plotExpressionProfileTFdel <- function(.cts1, .cts2, .cts3, .cts4,
                                       .info1, .info2, .info3, .info4,
                                       .name1 = "S. cerevisiae WT",
                                       .name2 = "S. paradoxus WT",
                                       .name3 = "S. cerevisiae TFdel",
                                       .name4 = "S. paradoxus TFdel",
                                       .color1 = "orange1",
                                       .color2 = "blue2",
                                       .color3 = "orange4",
                                       .color4 = "blue4",
                                       .normalization = c("none", "log2", "scale", "centered log2"),
                                       .show_points_wt = TRUE,
                                       .show_points_tfdel = TRUE,
                                       .show_lines_wt = TRUE,
                                       .show_lines_tfdel = FALSE,
                                       .show_wt = TRUE,
                                       .show_tfdel = TRUE,
                                       .show_confidence_intervals = TRUE,
                                       .plotlims = NULL) {
  if (.normalization == "none") {
    norm_func <- identity
    ylabel <- "Expression"
  }
  if (.normalization == "log2") {
    norm_func <- \(x) {log2(x + 1)}
    ylabel <- "Expression (log2)"
  }
  if (.normalization == "scale") {
    norm_func <- \(x) {t(scale(t(x)))}
    ylabel <- "Expression (centered and scaled)"
  }
  if (.normalization == "centered log2") {
    norm_func <- \(x) {log2(x + 1) - rowMeans(log2(x + 1), na.rm = TRUE)}
    ylabel <- "Expression (centered log2)"
  }
  if (!setequal(unique(.info1$experiment), unique(.info2$experiment))) {
    stop("sample info dataframes do not contain same set of experiments\n")
  }
  info1 <- tibble(experiment = "LowN",
                  time_point_str = .info1$time_point_str,
                  genotype = .info1$genotype,
                  well_flask_ID = .info1$well_flask_ID)
  info2 <- tibble(experiment = "LowN",
                  time_point_str = .info2$time_point_str,
                  genotype = .info2$genotype,
                  well_flask_ID = .info2$well_flask_ID)
  info3 <- tibble(experiment = "LowN",
                  time_point_str = .info3$time_point_str,
                  genotype = .info3$genotype,
                  well_flask_ID = .info3$well_flask_ID)
  info4 <- tibble(experiment = "LowN",
                  time_point_str = .info4$time_point_str,
                  genotype = .info4$genotype,
                  well_flask_ID = .info4$well_flask_ID)
  expr1 <- norm_func(.cts1) |> t()
  expr2 <- norm_func(.cts2) |> t()
  expr3 <- norm_func(.cts3) |> t()
  expr4 <- norm_func(.cts4) |> t()
  colnames(expr1) <- rownames(.cts1)
  colnames(expr2) <- rownames(.cts2)
  colnames(expr3) <- rownames(.cts3)
  colnames(expr4) <- rownames(.cts4)
  gdf1 <- bind_cols(expr1, info1) |> 
    pivot_longer(cols = colnames(expr1), names_to = "gene_name", values_to = "expr")
  gdf1$group_id <- "1"
  gdf2 <- bind_cols(expr2, info2) |> 
    pivot_longer(cols = colnames(expr2), names_to = "gene_name", values_to = "expr")
  gdf2$group_id <- "2"
  gdf3 <- bind_cols(expr3, info3) |> 
    pivot_longer(cols = colnames(expr3), names_to = "gene_name", values_to = "expr")
  gdf3$group_id <- "3"
  gdf4 <- bind_cols(expr4, info4) |> 
    pivot_longer(cols = colnames(expr4), names_to = "gene_name", values_to = "expr")
  gdf4$group_id <- "4"
  # converting each gene's expression to its mean expression between replicates
  gdf <- bind_rows(gdf1, gdf2, gdf3, gdf4) |> 
    group_by(group_id, time_point_str, well_flask_ID, genotype) |> 
    summarise(expr = mean(expr, na.rm = TRUE)) |> ungroup()
  plotlimdf <- gdf |> group_by(time_point_str, genotype, group_id) |>
    summarise(mean_expr = mean(expr, na.rm = TRUE),
              sd_expr = sd(expr, na.rm = TRUE)) 
  max_avg <- plotlimdf |> select(mean_expr) |>
    pull() |> max(na.rm = TRUE)
  min_avg <- plotlimdf |> select(mean_expr) |>
    pull() |> min(na.rm = TRUE)
  buffer <- 0.25
  max_avg <- max_avg + buffer
  min_avg <- min_avg - buffer
  # min/maxs when plotting points as well as averages:
  max_expr <- max(gdf$expr, na.rm = TRUE)
  min_expr <- min(gdf$expr, na.rm = TRUE)
  # setting plotlims if they haven't been manually set
  if (is.null(.plotlims)) {
    if (!(.show_points_wt | .show_points_tfdel)) {
      .plotlims <- c(min_avg, max_avg)
    }
    if (.show_points_wt | .show_points_tfdel) {
      .plotlims <- c(min_expr - 0.1, max_expr + 0.1) # buffer b/c TFdel points are very large
    }
  }
  plotdf <- gdf
  # plotting
  p <- ggplot() + 
    theme_classic() +
    scale_color_discrete(type = c(.color1, .color2, .color3, .color4), 
                         labels = c(.name1, .name2, .name3, .name4),
                         limits = c("1", "2", "3", "4")) +
    theme(legend.title = element_blank(),
          legend.position = "none") +
    ylab("") +
    xlab("") +
    ylim(.plotlims)
  if (.show_points_wt & .show_wt) {
    p <- p + geom_jitter(data = filter(plotdf, genotype == "WT"),
                         aes(x = time_point_str, y = expr, color = group_id), 
                         size = 0.1, alpha = 0.5, shape = ".") 
  }
  if (.show_lines_wt & .show_wt) {
    # lines trace average expr at each timepoint/experiment for each group
    avgexpr1 <- plotdf %>% filter(group_id == 1) |> group_by(time_point_str, group_id) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                          sd_expr = sd(expr, na.rm = TRUE))
    avgexpr2 <- plotdf %>% filter(group_id == 2) |> group_by(time_point_str, group_id) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                          sd_expr = sd(expr, na.rm = TRUE))
    # adding line segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
    p <- p +
      geom_line(data = avgexpr1, aes(x = time_point_str, y = mean_expr, group = group_id), color = .color1, linewidth = 1, linetype = "solid") +
      geom_line(data = avgexpr2, aes(x = time_point_str, y = mean_expr, group = group_id), color = .color2, linewidth = 1, linetype = "solid")
    if (.show_confidence_intervals) {
      # calculating 95% confidence in the mean
      nGenes <- nrow(.cts1)
      avgexpr1$CI <- 1.96*(avgexpr1$sd_expr/sqrt(nGenes))
      avgexpr2$CI <- 1.96*(avgexpr2$sd_expr/sqrt(nGenes))
      p <- p + 
        geom_ribbon(data = avgexpr1, aes(x = time_point_str, ymin = pmax(mean_expr - CI, .plotlims[1]), ymax = pmin(mean_expr + CI, .plotlims[2])),
                    fill = .color1, alpha = 0.3) +
        geom_ribbon(data = avgexpr2, aes(x = time_point_str, ymin = pmax(mean_expr - CI, .plotlims[1]), ymax = pmin(mean_expr + CI, .plotlims[2])),
                    fill = .color2, alpha = 0.3)
    }
    if (.show_points_tfdel & .show_tfdel) {
      p <- p + geom_jitter(data = filter(plotdf, genotype != "WT"),
                           aes(x = time_point_str, y = expr, color = group_id),
                           shape = "+", size = 6, width = 0.1, height = 0) # height = 0 prevents vertical jitter
    }
    if (.show_lines_tfdel & .show_tfdel) {
      # lines trace average expr at each timepoint/experiment for each group
      avgexpr3 <- plotdf %>% filter(group_id == 3) |> group_by(time_point_str, group_id) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                                      sd_expr = sd(expr, na.rm = TRUE))
      avgexpr4 <- plotdf %>% filter(group_id == 4) |> group_by(time_point_str, group_id) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                                      sd_expr = sd(expr, na.rm = TRUE))
      # adding line segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
      p <- p +
        geom_line(data = avgexpr3, aes(x = time_point_str, y = mean_expr, group = group_id), color = .color3, linewidth = 1, linetype = "dashed") +
        geom_line(data = avgexpr4, aes(x = time_point_str, y = mean_expr, group = group_id), color = .color4, linewidth = 1, linetype = "dashed")
      if (.show_confidence_intervals) {
        # calculating 95% confidence in the mean
        nGenes <- nrow(.cts1)
        avgexpr3$CI <- 1.96*(avgexpr3$sd_expr/sqrt(nGenes))
        avgexpr4$CI <- 1.96*(avgexpr4$sd_expr/sqrt(nGenes))
        p <- p + 
          geom_ribbon(data = avgexpr3, aes(x = time_point_str, ymin = pmax(mean_expr - CI, .plotlims[1]), ymax = pmin(mean_expr + CI, .plotlims[2])),
                      fill = .color1, alpha = 0.3) +
          geom_ribbon(data = avgexpr4, aes(x = time_point_str, ymin = pmax(mean_expr - CI, .plotlims[1]), ymax = pmin(mean_expr + CI, .plotlims[2])),
                      fill = .color2, alpha = 0.3)
      }
    }
    return(p)
  }
}
# # tests for plotExpressionProfileTFdel
# TODO: ARG81 1-2 at TP1 has un-normalized TFdel replicates near the top of the WT pack,
# but log2 has them somewhere in the middle
# gene_idxs <- gene_idxs <- finaldf |> filter(experiment == "LowN" &
#                                               group == "dyn21") |> 
#   select(gene_name) |> pull()
# cer_wt_cols <- which(sample_info_tfdel_allele$genotype == "WT" & sample_info_tfdel_allele$allele == "cer")
# par_wt_cols <- which(sample_info_tfdel_allele$genotype == "WT" & sample_info_tfdel_allele$allele == "par")
# cer_tfdel_cols <- which(sample_info_tfdel_allele$genotype == "LEU3delete" & sample_info_tfdel_allele$allele == "cer")
# par_tfdel_cols <- which(sample_info_tfdel_allele$genotype == "LEU3delete" & sample_info_tfdel_allele$allele == "par")
# plotExpressionProfileTFdel(.cts1 = counts_tfdel_allele[gene_idxs, cer_wt_cols],
#                            .cts2 = counts_tfdel_allele[gene_idxs, par_wt_cols],
#                            .cts3 = counts_tfdel_allele[gene_idxs, cer_tfdel_cols],
#                            .cts4 = counts_tfdel_allele[gene_idxs, par_tfdel_cols],
#                            .info1 = sample_info_tfdel_allele[cer_wt_cols,],
#                            .info2 = sample_info_tfdel_allele[par_wt_cols,],
#                            .info3 = sample_info_tfdel_allele[cer_tfdel_cols,],
#                            .info4 = sample_info_tfdel_allele[par_tfdel_cols,],
#                            .normalization = "log2",
#                            .show_points_wt = TRUE,
#                            .show_points_tfdel = TRUE,
#                            .show_lines_tfdel = TRUE)

# Plots each genotype as a separate line
plotExpressionProfileTFdels <- function(.cts,
                                        .info,
                                        .name1 = "WT",
                                        .name2 = "TFdel",
                                        .color1 = "red",
                                        .color2 = "grey",
                                        .normalization = c("none", "log2", "scale", "centered log2"),
                                        .ylims = NULL,
                                        .ylab = TRUE) {
  if (.normalization == "none") {
    norm_func <- identity
    ylabel <- "Expression"
  }
  if (.normalization == "log2") {
    norm_func <- \(x) {log2(x + 1)}
    ylabel <- "Expression (log2)"
  }
  if (.normalization == "scale") {
    norm_func <- scale
    ylabel <- "Expression (centered and scaled)"
  }
  if (.normalization == "centered log2") {
    norm_func <- \(x) {log2(x + 1) - rowMeans(log2(x + 1), na.rm = TRUE)}
    ylabel <- "Expression (centered log2)"
  }
  if (!.ylab) {
    ylabel <- ""
  }
  gdf <- bind_cols(tibble(sample_name = colnames(.cts)),
                   t(.cts)) |> 
    left_join(.info, by = "sample_name") |> 
    pivot_longer(cols = rownames(.cts), 
                 names_to = "gene_name", values_to = "expr")
  plotdf <- gdf |> group_by(genotype, gene_name) |> 
    reframe(expr = norm_func(expr),
            time_point_str = time_point_str) |>
    ungroup() |>
    group_by(genotype, time_point_str) |> 
    summarise(mean_expr = mean(expr, na.rm = TRUE))
  # plotting
  p <- ggplot(plotdf, aes(x = time_point_str, y = mean_expr)) +
    geom_line(data = filter(plotdf, genotype != "WT"),
              aes(group = genotype), color = .color2, alpha = 0.5) +
    geom_line(data = filter(plotdf, genotype == "WT"),
              aes(group = genotype), color = .color1,
              linewidth = 1.5) +
    theme_classic() +
    ylab(ylabel) +
    xlab("")
  if (!is.null(.ylims)) {
    p <- p + ylim(.ylims)
  }
  return(p)
}
# Tests for plotExpressionProfileTFdels
# test_gene_idxs <- finaldf |> filter(experiment == "LowN" &
#                                       cer == 2 &
#                                       par == 2) |>
#   select(gene_name) |> pull()
# test_counts <- counts_tfdel[test_gene_idxs, sample_info_tfdel$organism == "cer"]
# test_info <- sample_info_tfdel[sample_info_tfdel$organism == "cer",]
# plotExpressionProfileTFdels(.cts = test_counts, .info = test_info,
#                             .normalization = "none")

plotSingleProfilesTFdel <- function(.cts1, .cts2, .cts3, .cts4,
                                    .info1, .info2, .info3, .info4,
                                    .name1 = "S. cerevisiae WT",
                                    .name2 = "S. paradoxus WT",
                                    .name3 = "S. cerevisiae TFdel",
                                    .name4 = "S. paradoxus TFdel",
                                    .color1 = "orange1",
                                    .color2 = "blue2",
                                    .color3 = "orange4",
                                    .color4 = "blue4",
                                    .normalization = c("none", "log2", "scale", "centered log2"),
                                    .show_points_wt = TRUE,
                                    .show_points_tfdel = TRUE,
                                    .show_lines_wt = TRUE,
                                    .show_lines_tfdel = FALSE,
                                    .show_wt = TRUE,
                                    .show_tfdel = TRUE,
                                    .show_confidence_intervals = TRUE,
                                    .plotlims = NULL) {
  if (.normalization == "none") {
    norm_func <- identity
    ylabel <- "Expression"
  }
  if (.normalization == "log2") {
    norm_func <- \(x) {log2(x + 1)}
    ylabel <- "Expression (log2)"
  }
  if (.normalization == "scale") {
    norm_func <- \(x) {t(scale(t(x)))}
    ylabel <- "Expression (centered and scaled)"
  }
  if (.normalization == "centered log2") {
    norm_func <- \(x) {log2(x + 1) - rowMeans(log2(x + 1), na.rm = TRUE)}
    ylabel <- "Expression (centered log2)"
  }
  if (!setequal(unique(.info1$experiment), unique(.info2$experiment))) {
    stop("sample info dataframes do not contain same set of experiments\n")
  }
  info1 <- tibble(experiment = "LowN",
                  time_point_str = .info1$time_point_str,
                  genotype = .info1$genotype,
                  well_flask_ID = .info1$well_flask_ID)
  info2 <- tibble(experiment = "LowN",
                  time_point_str = .info2$time_point_str,
                  genotype = .info2$genotype,
                  well_flask_ID = .info2$well_flask_ID)
  info3 <- tibble(experiment = "LowN",
                  time_point_str = .info3$time_point_str,
                  genotype = .info3$genotype,
                  well_flask_ID = .info3$well_flask_ID)
  info4 <- tibble(experiment = "LowN",
                  time_point_str = .info4$time_point_str,
                  genotype = .info4$genotype,
                  well_flask_ID = .info4$well_flask_ID)
  expr1 <- norm_func(.cts1) |> t()
  expr2 <- norm_func(.cts2) |> t()
  expr3 <- norm_func(.cts3) |> t()
  expr4 <- norm_func(.cts4) |> t()
  colnames(expr1) <- rownames(.cts1)
  colnames(expr2) <- rownames(.cts2)
  colnames(expr3) <- rownames(.cts3)
  colnames(expr4) <- rownames(.cts4)
  gdf1 <- bind_cols(expr1, info1) |> 
    pivot_longer(cols = colnames(expr1), names_to = "gene_name", values_to = "expr")
  gdf1$group_id <- "1"
  gdf2 <- bind_cols(expr2, info2) |> 
    pivot_longer(cols = colnames(expr2), names_to = "gene_name", values_to = "expr")
  gdf2$group_id <- "2"
  gdf3 <- bind_cols(expr3, info3) |> 
    pivot_longer(cols = colnames(expr3), names_to = "gene_name", values_to = "expr")
  gdf3$group_id <- "3"
  gdf4 <- bind_cols(expr4, info4) |> 
    pivot_longer(cols = colnames(expr4), names_to = "gene_name", values_to = "expr")
  gdf4$group_id <- "4"
  # enumerating replicates, so each one can stack right on top of each other
  plotdf <- bind_rows(gdf1, gdf2, gdf3, gdf4)
  plotdf <- plotdf |> group_by(group_id, genotype, time_point_str, gene_name) |> 
    reframe(rep = rank(gene_name, ties.method = "first"),
            gene_name = gene_name,
            time_point_str = time_point_str,
            genotype = genotype,
            expr = expr,
            group_id = group_id)
  plotdf$rep <- if_else(plotdf$genotype == "WT",
                        true = 1, false = plotdf$rep)
  # plotting
  p <- ggplot() + 
    theme_classic() +
    scale_color_discrete(type = c(.color1, .color2, .color3, .color4), 
                         labels = c(.name1, .name2, .name3, .name4),
                         limits = c("1", "2", "3", "4")) +
    theme(legend.title = element_blank(),
          legend.position = "none") +
    ylab("") +
    xlab("")
  # lines trace average expr at each timepoint/experiment for each group
  avgexpr1 <- plotdf %>% filter(group_id == 1) |> group_by(time_point_str, group_id, gene_name) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                                             sd_expr = sd(expr, na.rm = TRUE))
  avgexpr2 <- plotdf %>% filter(group_id == 2) |> group_by(time_point_str, group_id, gene_name) |> summarise(mean_expr = mean(expr, na.rm = TRUE),
                                                                                                             sd_expr = sd(expr, na.rm = TRUE))
  p <- p +
    geom_line(data = avgexpr1, aes(x = time_point_str, y = mean_expr, group = interaction(group_id, gene_name)), color = .color1) +
    geom_line(data = avgexpr2, aes(x = time_point_str, y = mean_expr, group = interaction(group_id, gene_name)), color = .color2) +
    geom_line(data =  filter(plotdf, group_id == 3),
              position = position_nudge(x = 0.05),
              aes(x = time_point_str, y = expr,
                  color = group_id, group = interaction(group_id, time_point_str, gene_name))) +
    geom_point(data =  filter(plotdf, group_id == 3),
               position = position_nudge(x = 0.05),
               size = 0.25,
               aes(x = time_point_str, y = expr,
                   color = group_id)) +
    geom_line(data =  filter(plotdf, group_id == 4),
              position = position_nudge(x = -0.05),
              aes(x = time_point_str, y = expr,
                  color = group_id, group = interaction(group_id, time_point_str, gene_name))) +
    geom_point(data =  filter(plotdf, group_id == 4),
               position = position_nudge(x = -0.05),
               size = 0.25,
               aes(x = time_point_str, y = expr,
                   color = group_id)) +
  facet_wrap(~gene_name)
  return(p)
}
### proportional area plots, set of 4
plotPropArea <- function(x1, x2, x3, x4, 
                         .colors = levdyn_colordf$type,
                         .buffer = 5, .size_bounds = 8,
                         .text_bounds = 10) {
  quad_center <- sqrt(max(x1, x2, x3, x4)) + .buffer
  # Quadrant 1, level divergers
  sqx <- sqrt(x1)
  boundx <- if_else(sqx > .size_bounds, true = sqx/2, false = .text_bounds)
  df1 <- tibble(box_x = c(sqx[1]/2, -sqx[2]/2, -sqx[3]/2, sqx[4]/2) + quad_center,
                box_y = c(sqx[1]/2, sqx[2]/2, -sqx[3]/2, -sqx[4]/2) + quad_center,
                text_x = c(boundx[1], -boundx[2], -boundx[3], boundx[4]) + quad_center,
                text_y = c(boundx[1], boundx[2], -boundx[3], -boundx[4]) + quad_center,
                size = sqx, color = .colors[1], label = x1)
  # Quadrant 2, conserved
  sqx <- sqrt(x2)
  boundx <- if_else(sqx > .size_bounds, true = sqx/2, false = .text_bounds)
  df2 <- tibble(box_x = c(sqx[1]/2, -sqx[2]/2, -sqx[3]/2, sqx[4]/2) - quad_center,
                box_y = c(sqx[1]/2, sqx[2]/2, -sqx[3]/2, -sqx[4]/2) + quad_center,
                text_x = c(boundx[1], -boundx[2], -boundx[3], boundx[4]) - quad_center,
                text_y = c(boundx[1], boundx[2], -boundx[3], -boundx[4]) + quad_center,
                size = sqx, color = .colors[2], label = x2)
  # Quadrant 3, dynamics divergers
  sqx <- sqrt(x3)
  boundx <- if_else(sqx > .size_bounds, true = sqx/2, false = .text_bounds)
  df3 <- tibble(box_x = c(sqx[1]/2, -sqx[2]/2, -sqx[3]/2, sqx[4]/2) - quad_center,
                box_y = c(sqx[1]/2, sqx[2]/2, -sqx[3]/2, -sqx[4]/2) - quad_center,
                text_x = c(boundx[1], -boundx[2], -boundx[3], boundx[4]) - quad_center,
                text_y = c(boundx[1], boundx[2], -boundx[3], -boundx[4]) - quad_center,
                size = sqx, color = .colors[3], label = x3)
  # Quadrant 4, level and dynamics divergers
  sqx <- sqrt(x4)
  boundx <- if_else(sqx > .size_bounds, true = sqx/2, false = .text_bounds)
  df4 <- tibble(box_x = c(sqx[1]/2, -sqx[2]/2, -sqx[3]/2, sqx[4]/2) + quad_center,
                box_y = c(sqx[1]/2, sqx[2]/2, -sqx[3]/2, -sqx[4]/2) - quad_center,
                text_x = c(boundx[1], -boundx[2], -boundx[3], boundx[4]) + quad_center,
                text_y = c(boundx[1], boundx[2], -boundx[3], -boundx[4]) - quad_center,
                size = sqx, color = .colors[4], label = x4)
  df <- bind_rows(df1, df2, df3, df4)
  mm <- max(df$size)*1.1
  ggplot(data=df, aes(x = box_x, y = box_y, width=size, height=size, 
                      group=factor(size))) +
    geom_tile(fill = df$color) +
    geom_text(data = filter(df, size > .size_bounds), 
              aes(label = label, x = text_x, y = text_y),
              col="white", size=5) +
    geom_text(data = filter(df, size <= .size_bounds), 
              aes(label=label, x = text_x, y = text_y),
              col="black", size=5) +
    geom_hline(aes(yintercept = quad_center), linewidth = 0.8) +
    geom_hline(aes(yintercept = -quad_center), linewidth = 0.8) +
    geom_vline(aes(xintercept = quad_center), linewidth = 0.8) +
    geom_vline(aes(xintercept = -quad_center), linewidth = 0.8) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none")
}
# # tests for plotPropArea
# plotPropArea(x1 = c(1870, 683, 23, 0), x2 = c(17, 32, 24, 1),
#              x3 = c(177, 436, 148, 89), x4 = c(17, 30, 16, 10),
#              .colors = c("red", "green", "blue", "yellow"))

### proportional area plots, just one set of 4
plotPropAreaSingle <- function(.counts,
                               .colors = levdyn_colordf$type,
                               .size_bounds = 14,
                               .text_bounds = 10,
                               .text_size = 8,
                               .xlims = NULL,
                               .ylims = NULL) {
  sqx <- sqrt(.counts)
  boundx <- if_else(sqx > .size_bounds, true = sqx/2, false = .text_bounds)
  df <- tibble(box_x = c(sqx[1]/2, -sqx[2]/2, -sqx[3]/2, sqx[4]/2),
               box_y = c(sqx[1]/2, sqx[2]/2, -sqx[3]/2, -sqx[4]/2),
               text_x = c(boundx[1], -boundx[2], -boundx[3], boundx[4]),
               text_y = c(boundx[1], boundx[2], -boundx[3], -boundx[4]),
               size = sqx, color = .colors, label = .counts)
  mm <- max(df$size)*1.1
  if (is.null(.xlims)) {
    .xlims <- c(-max(sqx), max(sqx))
  }
  if (is.null(.ylims)) {
    .ylims <- c(-max(sqx), max(sqx))
  }
  ggplot(data=df, aes(x = box_x, y = box_y, width=size, height=size,
                      group=factor(size))) +
    geom_tile(fill = df$color) +
    geom_text(data = filter(df, size > .size_bounds),
              aes(label = label, x = text_x, y = text_y),
              col = "white", size = .text_size) +
    geom_text(data = filter(df, size <= .size_bounds),
              aes(label=label, x = text_x, y = text_y),
              col = "black", size = .text_size) +
    geom_hline(aes(yintercept = 0), linewidth = 1) +
    geom_vline(aes(xintercept = 0), linewidth = 1) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none") +
    xlim(.xlims) +
    ylim(.ylims)
}
# # tests for plotPropAreaSingle
# plotPropAreaSingle(.counts = c(1870, 17, 177, 10),
#                    .colors = c("red", "green", "blue", "yellow"),
#                    .size_bounds = 12, .text_size = 5)

# Wrapper function for plotExpressionProfileTFdel
plotGenesTFdel <- function(.gene_idxs, .tf, .parents_or_hybrid = "parents",
                           .plotlims = NULL, .single_genes = FALSE,
                           .normalization = "log2") {
  if (.parents_or_hybrid == "parents" & !.single_genes) {
    p <- plotExpressionProfileTFdel(.cts1 = counts_tfdel[.gene_idxs,
                                                    sample_info_tfdel$organism == "cer" &
                                                      sample_info_tfdel$genotype == "WT",
                                                    drop = FALSE],
                               .cts2 = counts_tfdel[.gene_idxs,
                                                    sample_info_tfdel$organism == "par" &
                                                      sample_info_tfdel$genotype == "WT",
                                                    drop = FALSE],
                               .cts3 = counts_tfdel[.gene_idxs,
                                                    sample_info_tfdel$organism == "cer" &
                                                      sample_info_tfdel$genotype == paste0(.tf, "delete"),
                                                    drop = FALSE],
                               .cts4 = counts_tfdel[.gene_idxs,
                                                    sample_info_tfdel$organism == "par" &
                                                      sample_info_tfdel$genotype == paste0(.tf, "delete"),
                                                    drop = FALSE],
                               .info1 = sample_info_tfdel[sample_info_tfdel$organism == "cer" &
                                                            sample_info_tfdel$genotype == "WT",],
                               .info2 = sample_info_tfdel[sample_info_tfdel$organism == "par" &
                                                            sample_info_tfdel$genotype == "WT",],
                               .info3 = sample_info_tfdel[sample_info_tfdel$organism == "cer" &
                                                            sample_info_tfdel$genotype == paste0(.tf, "delete"),],
                               .info4 = sample_info_tfdel[sample_info_tfdel$organism == "par" &
                                                            sample_info_tfdel$genotype == paste0(.tf, "delete"),],
                               .normalization = .normalization, .plotlims = .plotlims)
  }
  if (.parents_or_hybrid == "hybrid" & !.single_genes) {
    p <- plotExpressionProfileTFdel(.cts1 = counts_tfdel_allele[.gene_idxs,
                                                           sample_info_tfdel_allele$allele == "cer" &
                                                             sample_info_tfdel_allele$genotype == "WT",
                                                           drop = FALSE],
                               .cts2 = counts_tfdel_allele[.gene_idxs,
                                                           sample_info_tfdel_allele$allele == "par" &
                                                             sample_info_tfdel_allele$genotype == "WT",
                                                           drop = FALSE],
                               .cts3 = counts_tfdel_allele[.gene_idxs,
                                                           sample_info_tfdel_allele$allele == "cer" &
                                                             sample_info_tfdel_allele$genotype == paste0(.tf, "delete"),
                                                           drop = FALSE],
                               .cts4 = counts_tfdel_allele[.gene_idxs,
                                                           sample_info_tfdel_allele$allele == "par" &
                                                             sample_info_tfdel_allele$genotype == paste0(.tf, "delete"),
                                                           drop = FALSE],
                               .info1 = sample_info_tfdel_allele[sample_info_tfdel_allele$allele == "cer" &
                                                                   sample_info_tfdel_allele$genotype == "WT",],
                               .info2 = sample_info_tfdel_allele[sample_info_tfdel_allele$allele == "par" &
                                                                   sample_info_tfdel_allele$genotype == "WT",],
                               .info3 = sample_info_tfdel_allele[sample_info_tfdel_allele$allele == "cer" &
                                                                   sample_info_tfdel_allele$genotype == paste0(.tf, "delete"),],
                               .info4 = sample_info_tfdel_allele[sample_info_tfdel_allele$allele == "par" &
                                                                   sample_info_tfdel_allele$genotype == paste0(.tf, "delete"),],
                               .normalization = .normalization, .plotlims = .plotlims)
  }
  if (.parents_or_hybrid == "parents" & .single_genes) {
    p <- plotSingleProfilesTFdel(.cts1 = counts_tfdel[.gene_idxs,
                                                      sample_info_tfdel$organism == "cer" &
                                                        sample_info_tfdel$genotype == "WT",
                                                      drop = FALSE],
                                 .cts2 = counts_tfdel[.gene_idxs,
                                                      sample_info_tfdel$organism == "par" &
                                                        sample_info_tfdel$genotype == "WT",
                                                      drop = FALSE],
                                 .cts3 = counts_tfdel[.gene_idxs,
                                                      sample_info_tfdel$organism == "cer" &
                                                        sample_info_tfdel$genotype == paste0(.tf, "delete"),
                                                      drop = FALSE],
                                 .cts4 = counts_tfdel[.gene_idxs,
                                                      sample_info_tfdel$organism == "par" &
                                                        sample_info_tfdel$genotype == paste0(.tf, "delete"),
                                                      drop = FALSE],
                                 .info1 = sample_info_tfdel[sample_info_tfdel$organism == "cer" &
                                                              sample_info_tfdel$genotype == "WT",],
                                 .info2 = sample_info_tfdel[sample_info_tfdel$organism == "par" &
                                                              sample_info_tfdel$genotype == "WT",],
                                 .info3 = sample_info_tfdel[sample_info_tfdel$organism == "cer" &
                                                              sample_info_tfdel$genotype == paste0(.tf, "delete"),],
                                 .info4 = sample_info_tfdel[sample_info_tfdel$organism == "par" &
                                                              sample_info_tfdel$genotype == paste0(.tf, "delete"),],
                                 .normalization = .normalization, .plotlims = .plotlims)
  }
  if (.parents_or_hybrid == "hybrid" & !.single_genes) {
    p <- plotSingleProfilesTFdel(.cts1 = counts_tfdel_allele[.gene_idxs,
                                                             sample_info_tfdel_allele$allele == "cer" &
                                                               sample_info_tfdel_allele$genotype == "WT",
                                                             drop = FALSE],
                                 .cts2 = counts_tfdel_allele[.gene_idxs,
                                                             sample_info_tfdel_allele$allele == "par" &
                                                               sample_info_tfdel_allele$genotype == "WT",
                                                             drop = FALSE],
                                 .cts3 = counts_tfdel_allele[.gene_idxs,
                                                             sample_info_tfdel_allele$allele == "cer" &
                                                               sample_info_tfdel_allele$genotype == paste0(.tf, "delete"),
                                                             drop = FALSE],
                                 .cts4 = counts_tfdel_allele[.gene_idxs,
                                                             sample_info_tfdel_allele$allele == "par" &
                                                               sample_info_tfdel_allele$genotype == paste0(.tf, "delete"),
                                                             drop = FALSE],
                                 .info1 = sample_info_tfdel_allele[sample_info_tfdel_allele$allele == "cer" &
                                                                     sample_info_tfdel_allele$genotype == "WT",],
                                 .info2 = sample_info_tfdel_allele[sample_info_tfdel_allele$allele == "par" &
                                                                     sample_info_tfdel_allele$genotype == "WT",],
                                 .info3 = sample_info_tfdel_allele[sample_info_tfdel_allele$allele == "cer" &
                                                                     sample_info_tfdel_allele$genotype == paste0(.tf, "delete"),],
                                 .info4 = sample_info_tfdel_allele[sample_info_tfdel_allele$allele == "par" &
                                                                     sample_info_tfdel_allele$genotype == paste0(.tf, "delete"),],
                                 .normalization = .normalization, .plotlims = .plotlims)
  }
  return(p)
}
# # test for plotGenesTFdel
# gene_idxs <- c("YBR163W", "YBR265W", "YGL048C", "YGR101W", "YJL085W", "YLR085C", "YLR118C",
#                "YLR274W", "YLR283W", "YOR047C", "YOR201C", "YOR315W", "YPR105C", "YPR129W",
#                "YPR198W")
# plotGenesTFdel(.tf = "GCR2", .gene_idxs = gene_idxs)
# plotGenesTFdel(.tf = "GCR2", .gene_idxs = gene_idxs, .single_genes = TRUE)

# Upset plot
# given a group name, which must be a column in .df,
# members of that group (possible values in that column, such as TFs in the deletion column),
# items in the group, other column names in .df (usually gene names, but can be gene-effect-TF combos if organism is group)
# creates an upset plot (3+ group Venn diagram) 
# which enumerates which combinations of groups share how many items
makeUpsetPlot <- function(.df, .group_name, .group_members, .item_names,
                          .min_comb_size = 5) {
  lt <- vector(mode = "list", length = 0)
  for (grpmem in .group_members) {
    lt[[grpmem]] <- filter(.df, .data[[.group_name]] == grpmem)
  }
  lt <- lt |> 
    map(.f = select, .item_names) |> 
    map(.f = \(x) {purrr::reduce(x, .f = paste0)})
  plotdf <- make_comb_mat(lt)
  plotdf <- plotdf[,comb_size(plotdf) >= .min_comb_size]
  p <- UpSet(plotdf, 
             set_order = .group_members,
             comb_order = order(comb_size(plotdf)),
             top_annotation = HeatmapAnnotation( 
               "number of genes" = anno_barplot(comb_size(plotdf), 
                                                ylim = c(0, max(comb_size(plotdf))*1.1),
                                                border = FALSE, 
                                                gp = gpar(fill = "black"), 
                                                height = unit(4, "cm")), 
               annotation_name_side = "left", 
               annotation_name_rot = 90))
  draw(p)
  suppressWarnings(decorate_annotation("number of genes", {
    grid.text(comb_size(plotdf)[column_order(p)], x = seq_along(comb_size(plotdf)), y = unit(comb_size(plotdf)[column_order(p)], "native") + unit(2, "pt"), 
              default.units = "native", just = c("left", "bottom"), 
              gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
  }))
}
# tests for makeUpsetPlot
# TF example within on organism
# makeUpsetPlot(.df = filter(effectdf, organism == "hyc" &
#                              effect == "dynamics"),
#               .group_name = "deletion",
#               .group_members = c("GAT1", "URE2", "GLN3"),
#               .item_names = "gene_name")
# # organism venn diagram of sharing gene-TF-effect items
# makeUpsetPlot(.df = effectdf,
#               .group_name = "organism",
#               .group_members = c("cer", "par", "hyc", "hyp"),
#               .item_names = c("gene_name", "deletion", "effect"))

# makes heatmap where rows are TF deletions, columns are groups,
# and counts are the number of significant effects
# @input: .df: dataframe with columns lfc, padj, deletion, gene_name, and grouping columns
# .groups: character vector of names of grouping columns
# @output: heatmap, through print function
# Note: groups should be in order from last-to-vary to first-to-vary.
# For example, if you have two groups: dynamics (cons, div) 
# and organism (cer, par), giving .groups = c("organism", "dynamics)
# would order columns left to right: cer_cons, cer_div, par_cons, par_div
makeGeneGroupHeatmap <- function(.df, .tf_order,
                                 .groups,
                                 .col_fun = colorRamp2(c(0, 10, 30, 100), c("blue", "yellow", "red", "magenta")),
                                 .legend = FALSE,
                                 .title = NULL) {
  griddf <- select(.df, c("deletion", .groups)) |> 
    filter(deletion %in% .tf_order) |>
    unique() |>
    select(-deletion) |> 
    expand_grid(deletion = .tf_order) |>
    unique() |> 
    arrange(across(.groups))
  effectsdf <- .df |> 
    filter(deletion %in% .tf_order) |> 
    group_by(across(all_of(c("deletion", .groups)))) |> 
    summarise(nGenes = sum(padj < p_thresh))
  plotdf <- left_join(griddf, effectsdf, by = colnames(griddf)) |> 
    mutate(nGenes = if_else(is.na(nGenes),
                            true = 0,
                            false = nGenes)) |> # group_by/summarise drops groups without any counts, and at the time of this code drop=FALSE wasn't implemented
    pivot_wider(id_cols = "deletion",
                names_from = .groups,
                values_from = "nGenes") |>
    ungroup()
  effects_mat <- plotdf |> select(-deletion) |> as.matrix()
  rownames(effects_mat) <- plotdf$deletion
  effects_mat <- effects_mat[rownames(effects_mat) %in% .tf_order,]
  Heatmap(effects_mat, col = .col_fun, na_col = "grey80",
          row_order = .tf_order,
          show_heatmap_legend = .legend,
          column_title = .title,
          column_order = colnames(effects_mat),
          cell_fun = function(j, i, x, y, width, height, fill) {
            output <- if_else(!(is.na(effects_mat[i, j])), 
                              true = as.character(effects_mat[i, j]), 
                              false = "-")
            grid.text(output, x, y, gp = gpar(fontsize = 10))
          })
}
# # tests for makeGeneGroupHeatmap
# # conserved TF effects TP1
# makeGeneGroupHeatmap(.df = filter(TFdeldf, 
#                                   tf_effect_conserved == TRUE &
#                                     timepoint == "TP1" &
#                                     organism %in% c("cer", "par")),
#                      .tf_order = tf_order,
#                      .groups = c("dynamics", "cer", "par", "lfc_sign", "organism"))
# # QC check that effects match what is seen in the heatmap:
# filter(TFdeldf, 
#        tf_effect_conserved == TRUE &
#          timepoint == "TP1" &
#          organism %in% c("cer", "par") &
#          cer == "0" & par == "0" & lfc_sign == -1 &
#          deletion == "GLN3") |> 
#   group_by(organism) |> summarise(nGenes = sum(padj < p_thresh))

# basic volcano plot to compare power
makeVolcanoPlot <- function(.tfdeldf, .tf, .org, .timepoint) {
  plotdf <- .tfdeldf |> filter(organism == .org &
                                 timepoint == .timepoint &
                                 deletion == .tf)
  ggplot(plotdf, aes(x = lfc, y = -log10(pval))) +
    geom_point(aes(color = padj < p_thresh)) +
    ylim(c(0, 15)) +
    xlim(c(-5, 5))
}
#### Gene Ontology Enrichment ####
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
    n_genes <- length(unique(finaldf$gene_name))
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
# test <- getGOSlimDf(.idxs = c("YGR192C", "YJR009C", "YJL052W"),
#                     .group_name = "tdhs", .min_hits = 1)
########################### Archive ###################################
################## WGCNA Module version ####################### 
# Archived when we switched to direct 
# correlation clustering of both species at once
# #### optional: including Fay et al. 2023 temp data ####
# # cer, par
# load("data_files/Cleaned_Fay_Counts.RData")
# fay <- fay[rownames(fay) %in% rownames(counts),] # ~500 mostly lowly expressed genes filtered out
# # hyc hyp
# load("data_files/Cleaned_Fay_Counts_allele.RData")
# fay_allele <- fay_allele[rownames(fay_allele) %in% rownames(counts),]
# # there's no easy way to left_join matrices unfortunately
# missing_fay_genes <- rownames(counts)[!(rownames(counts) %in% rownames(fay))]
# missing_fay <- matrix(NA, nrow = length(missing_fay_genes), ncol = ncol(fay))
# rownames(missing_fay) <- missing_fay_genes
# colnames(missing_fay) <- colnames(fay)
# missing_fay_allele <- matrix(NA, nrow = length(missing_fay_genes), ncol = ncol(fay_allele))
# rownames(missing_fay_allele) <- missing_fay_genes
# colnames(missing_fay_allele) <- colnames(fay_allele)
# fay <- rbind(fay, missing_fay)
# fay_allele <- rbind(fay_allele, missing_fay_allele)
# fay <- fay[rownames(counts),]
# fay_allele <- fay_allele[rownames(counts),]
# sum(rownames(fay) == rownames(counts))/nrow(counts)
# sum(rownames(fay_allele) == rownames(counts_allele))/nrow(counts_allele) # cbind doesn't check that rownames match, so make sure they match beforehand
# # joining
# counts <- cbind(counts, fay)
# counts_allele <- cbind(counts_allele, fay_allele)
# sample_info <- bind_rows(sample_info, sample_info_fay) # bind rows thankfully does check that the colnames are the same
# sample_info_allele <- bind_rows(sample_info_allele, sample_info_fay_allele)
# 
# # creating random sets of colors that have same number of each true color set
# random_colors00 <- list(cer = sample(colors00$cer, size = length(colors00$cer), replace = FALSE),
#                         par = sample(colors00$par, size = length(colors00$par), replace = FALSE))
# random_colors10 <- list(cer = sample(colors10$cer, size = length(colors10$cer), replace = FALSE),
#                         par = sample(colors10$par, size = length(colors10$par), replace = FALSE))
# random_colors25 <- list(cer = sample(colors25$cer, size = length(colors25$cer), replace = FALSE),
#                         par = sample(colors25$par, size = length(colors25$par), replace = FALSE))
# random_colors35 <- list(cer = sample(colors35$cer, size = length(colors35$cer), replace = FALSE),
#                         par = sample(colors35$par, size = length(colors35$par), replace = FALSE))
# 
# #### Condition-matched counts/info for parents and hybrids ####
# 
# ### condition-matching samples so that expression similarity can be calculated
# counts_cer <- counts[, sample_info$organism == "cer"] |> t()
# info_cer <- sample_info %>% filter(sample_name %in% rownames(counts_cer))
# counts_par <- counts[, sample_info$organism == "par"] |> t()
# info_par <- sample_info %>% filter(sample_name %in% rownames(counts_par))
# counts_hyc <- counts_allele[, sample_info_allele$organism == "hyb" & sample_info_allele$allele == "cer"] |> t()
# info_hyc <- sample_info_allele %>% filter(sample_name %in% rownames(counts_hyc))
# counts_hyp <- counts_allele[, sample_info_allele$organism == "hyb" & sample_info_allele$allele == "par"] |> t()
# info_hyp <- sample_info_allele %>% filter(sample_name %in% rownames(counts_hyp))
# 
# # Filtering out conditions (combinations of genotype, experiment, and timepoint) that are not represented in cer, par, and hyb
# common_conditions <- bind_rows(info_cer, info_par, info_hyc, info_hyp) %>% 
#   select(condition, allele, organism) |> unique() %>% 
#   pivot_wider(id_cols = c("condition"), names_from = c("allele", "organism"), values_from = "organism") %>% 
#   drop_na() %>% 
#   select(condition) |> pull()
# 
# # filtering down to common conditions
# info_cer <- filter(info_cer, condition %in% common_conditions) %>% arrange(experiment, condition, well_flask_ID)
# info_par <- filter(info_par, condition %in% common_conditions) %>% arrange(experiment, condition, well_flask_ID)
# info_hyc <- filter(info_hyc, condition %in% common_conditions) %>% arrange(experiment, condition, well_flask_ID)
# info_hyp <- filter(info_hyp, condition %in% common_conditions) %>% arrange(experiment, condition, well_flask_ID)
# 
# # checking conditions are the same (these should both be 1)
# sum(unique(info_cer$condition) == unique(info_par$condition))/length(unique(info_cer$condition))
# sum(unique(info_hyc$condition) == unique(info_hyp$condition))/length(unique(info_hyc$condition))
# 
# # paring down counts to match condition-matched info dfs
# counts_cer <- counts_cer[info_cer$sample_name,] # remember genes are COLUMNS in WGCNA
# counts_par <- counts_par[info_par$sample_name,]
# counts_hyc <- counts_hyc[info_hyc$sample_name,]
# counts_hyp <- counts_hyp[info_hyp$sample_name,]
# 
# # removing TFdel from WT datasets
# infos <- list(cer = info_cer[info_cer$genotype == "WT",],
#               par = info_par[info_par$genotype == "WT",])
# infos_allele <- list(cer = info_hyc[info_hyc$genotype == "WT",],
#                      par = info_hyp[info_hyp$genotype == "WT",])
# counts_cer <- counts_cer[info_cer$genotype == "WT",]
# counts_par <- counts_par[info_par$genotype == "WT",]
# counts_hyc <- counts_hyc[info_hyc$genotype == "WT",]
# counts_hyp <- counts_hyp[info_hyp$genotype == "WT",]
# 
# # taking mean count between replicates for WT datasets (for TFdel we definitely want replicates)
# WT_common_conditions <- infos$cer$condition |> unique()
# counts_all2 <- list(cer = matrix(0, nrow = length(WT_common_conditions), ncol = ncol(counts_cer),
#                                  dimnames = list(WT_common_conditions, colnames(counts_cer))),
#                     par = matrix(0, nrow = length(WT_common_conditions), ncol = ncol(counts_cer),
#                                  dimnames = list(WT_common_conditions, colnames(counts_cer))))
# 
# counts_all2_allele <- list(cer = matrix(0, nrow = length(WT_common_conditions), ncol = ncol(counts_cer),
#                                         dimnames = list(WT_common_conditions, colnames(counts_cer))),
#                            par = matrix(0, nrow = length(WT_common_conditions), ncol = ncol(counts_cer),
#                                         dimnames = list(WT_common_conditions, colnames(counts_cer))))
# # just to have a record of how many replicates each condition has
# replicatedf <- tibble(condition = WT_common_conditions,
#                       nreps_cer = 0,
#                       nreps_par = 0,
#                       nreps_hyc = 0,
#                       nreps_hyp = 0)
# for (cc in WT_common_conditions) {
#   cat(cc, "\n")
#   samps_cer <- infos$cer |> filter(condition == cc) |> select(sample_name) |> pull()
#   samps_par <- infos$par |> filter(condition == cc) |> select(sample_name) |> pull()
#   samps_hyc <- infos_allele$cer |> filter(condition == cc) |> select(sample_name) |> pull()
#   samps_hyp <- infos_allele$par |> filter(condition == cc) |> select(sample_name) |> pull()
#   # recording number of replicates per species per condition
#   replicatedf[replicatedf$condition == cc,]$nreps_cer <- length(samps_cer)
#   replicatedf[replicatedf$condition == cc,]$nreps_par <- length(samps_par)
#   replicatedf[replicatedf$condition == cc,]$nreps_hyc <- length(samps_hyc)
#   replicatedf[replicatedf$condition == cc,]$nreps_hyp <- length(samps_hyp)
#   # taking mean across replicates
#   gene_vec_cer <- counts_cer[samps_cer,,drop = FALSE] |> apply(MARGIN = 2, FUN = mean) |> as.numeric()
#   gene_vec_par <- counts_par[samps_par,,drop = FALSE] |> apply(MARGIN = 2, FUN = mean) |> as.numeric()
#   gene_vec_hyc <- counts_hyc[samps_hyc,,drop = FALSE] |> apply(MARGIN = 2, FUN = mean) |> as.numeric()
#   gene_vec_hyp <- counts_hyp[samps_hyp,,drop = FALSE] |> apply(MARGIN = 2, FUN = mean) |> as.numeric()
#   # assigning mean counts
#   counts_all2$cer[cc,] <- gene_vec_cer
#   counts_all2$par[cc,] <- gene_vec_par
#   counts_all2_allele$cer[cc,] <- gene_vec_hyc
#   counts_all2_allele$par[cc,] <- gene_vec_hyp
# }
# 
# # The only information we need from info now is condition, timepoint, experiment, and genotype,
# # which should be all the same info for each species
# # so joining all 4 is overkill, but it helps make sure they are all the same
# info <- inner_join(select(infos$cer, condition, time_point_num, time_point_str, experiment, genotype),
#                    select(infos$par, condition, time_point_num, time_point_str, experiment, genotype),
#                    by = c("condition", "time_point_num", "time_point_str", "experiment", "genotype"),
#                    relationship = "many-to-many") |> unique() |> 
#   left_join(select(infos_allele$cer, condition, time_point_num, time_point_str, experiment, genotype),
#             by = c("condition", "time_point_num", "time_point_str", "experiment", "genotype"),
#             relationship = "one-to-many") |> unique() |>  # don't ask my why I still need to filter out duplicates even though I did a left_join
#   left_join(select(infos_allele$par, condition, time_point_num, time_point_str, experiment, genotype),
#             by = c("condition", "time_point_num", "time_point_str", "experiment", "genotype"),
#             relationship = "one-to-many") |> unique()
# 
# # removing datastructures we won't use anymore for clarity
# rm(info_cer, info_par, info_hyc, info_hyp, sample_info, sample_info_allele, 
#    counts_cer, counts_par, counts_hyc, counts_hyp, counts, counts_allele)
# 
# #### Figure 1/2 bipartite construction and configuration model ####
# 
# # @input: 1) 2 vectors of same length containing WGCNA-derived color names for each gene in each species,
# #         2) name of current connection (color pair) to test
# # @output: TRUE/FALSE indicating whether test connection is stronger than would be
# #          expected by random chance, 95% confidence, based on permutation tests
# callConservedConnection <- function(.color_name1, .color_name2, .colors1, .colors2, .n_permutations_per_species = 50000) {
#   if (length(.colors1) != length(.colors2)) {
#     stop("color vectors are not the same length!\n")
#   }
#   n_connections <- sum(.colors1 == .color_name1 & .colors2 == .color_name2)
#   # helper function
#   scrambleConnections <- function(.steady_name, .scramble_name, .steady_colors, .scramble_colors) {
#     scrambled <- sample(.scramble_colors, size = length(.scramble_colors), replace = FALSE)
#     n_connections <- sum(.steady_colors == .steady_name & scrambled == .scramble_name)
#     return(n_connections)
#   }
#   # scrambling first set of colors
#   null_connections1 <- sapply(c(1:.n_permutations_per_species), \(dummy) {
#     output <- scrambleConnections(.steady_name = .color_name1, 
#                                   .scramble_name = .color_name2,
#                                   .steady_colors = .colors1,
#                                   .scramble_colors = .colors2)
#     return(output)
#   })
#   # scrambling second set of colors
#   null_connections2 <- sapply(c(1:.n_permutations_per_species), \(dummy) {
#     output <- scrambleConnections(.steady_name = .color_name2, 
#                                   .scramble_name = .color_name1,
#                                   .steady_colors = .colors2,
#                                   .scramble_colors = .colors1)
#     return(output)
#   })
#   null_connections <- c(null_connections1, null_connections2)
#   return(sum(null_connections >= n_connections)/length(null_connections))
# }
# # # tests for callConservedConnection
# # # positive control, known to be extremely conserved, should be 0
# # callConservedConnection(.color_name1 = "green",
# #                         .color_name2 = "green",
# #                         .colors1 = colors10$cer,
# #                         .colors2 = colors10$par)
# # 
# # # negative control, should be 1
# # callConservedConnection(.color_name1 = "red",
# #                         .color_name2 = "greenyellow",
# #                         .colors1 = rep("red", 30),
# #                         .colors2 = rep("greenyellow", 30))
# 
# # Function for constructing the matrix used by bipartite::plotweb
# # @input: two vectors (of same length, length nGenes) giving module 
# #         assignments (colors) of each gene in both comparison groups (usually cer vs par)
# # @output: matrix of nColors_top x nColors_bottom where entries are the
# #          number of connections between each color group
# makeBipartiteMatrix <- function(.top_colors, .bottom_colors, .simplify_thresh = 0) {
#   #module_names_top <- setdiff(unique(.top_colors), "grey")
#   #module_names_bottom <- setdiff(unique(.bottom_colors), "grey") # uncomment to exclude grey
#   module_names_top <- unique(.top_colors)
#   module_names_bottom <- unique(.bottom_colors)
#   color_counts_matrix <- matrix(nrow = length(module_names_bottom), ncol = length(module_names_top), dimnames = list(module_names_bottom, module_names_top))
#   for (col_top in module_names_top) {
#     for (col_bottom in module_names_bottom) {
#       x <- sum((.bottom_colors == col_bottom) & (.top_colors == col_top))
#       if (x < .simplify_thresh) {
#         x <- 0 # TODO: is this safe? Or does it rearrange the answer?
#       }
#       color_counts_matrix[col_bottom, col_top] <- x
#     }
#   }
#   return(list(matrix = color_counts_matrix, colors_high = module_names_top, colors_low = module_names_bottom))
# }
# # tests for makeBipartiteMatrix
# # makeBipartiteMatrix(c("red", "red", "green", "red", "grey"), c("green", "green", "green", "red", "grey"))
# 
# #### Figure 1/2 CCM dataframes ####
# 
# if (!file.exists("data_files/CCMs.RData")) {
#   # calling CCMs, which will be indicated on the bipartite graph
#   # @input: color vectors of length nGenes for pair of colors being compared
#   # @output: CCM dataframe giving the p-value of a configuration model for the 
#   #          given number of connections (block size) between each pair of colors
#   makeCCMDataFrame <- function(.clrs1, .clrs2, .clrs_name1 = "cer_color", .clrs_name2 = "par_color") {
#     each_color1 <- rep(unique(.clrs1), length(unique(.clrs2))) %>% as.character()
#     each_color2 <- sapply(unique(.clrs2), rep, times = length(unique(.clrs1))) %>% as.character()
#     ccmdf <- tibble(col1 = each_color1,
#                     col2 = each_color2)
#     names(ccmdf) <- c(.clrs_name1, .clrs_name2)
#     ccmdf$p_value <- map2(ccmdf[,.clrs_name1, drop = TRUE], ccmdf[,.clrs_name2, drop = TRUE], \(x, y) {
#       cat("starting on", x, "and", y, "\n")
#       return(callConservedConnection(.color_name1 = x,
#                                      .color_name2 = y,
#                                      .colors1 = .clrs1,
#                                      .colors2 = .clrs2))
#     }) %>% unlist()
#     ccmdf$block_size <- map2(ccmdf[,.clrs_name1, drop = TRUE], ccmdf[,.clrs_name2, drop = TRUE], \(x, y) {
#       return(sum(.clrs1 == x & .clrs2 == y))
#     }) %>% unlist()
#     return(ccmdf)
#   }
#   # # tests for make CCMdf
#   # test_colors_cer <- sample(c("brown","green","yellow"), size = length(colors10$cer), replace = TRUE)
#   # test_colors_par <- sample(c("brown","green","yellow"), size = length(colors10$par), replace = TRUE)
#   # test_CCM <- makeCCMDataFrame(.clrs1 = test_colors_cer,
#   #                              .clrs2 = test_colors_par)
#   
#   # no module merging (smallest modules)
#   CCMdf00 <- makeCCMDataFrame(colors00$cer, colors00$par)
#   
#   # colors 10 version (mid-size modules, where modules were merged if their eigengenes had over 0.9 correlation)
#   CCMdf10 <- makeCCMDataFrame(colors10$cer, colors10$par)
#   
#   # module merging at .75 (mid-size modules)
#   CCMdf25 <- makeCCMDataFrame(colors25$cer, colors25$par)
#   
#   # module merging at .65 (largest modules)
#   CCMdf35 <- makeCCMDataFrame(colors35$cer, colors35$par)
#   
#   # same but for random colors sets
#   random_CCMdf00 <- makeCCMDataFrame(random_colors00$cer, random_colors00$par)
#   random_CCMdf10 <- makeCCMDataFrame(random_colors10$cer, random_colors10$par)
#   random_CCMdf25 <- makeCCMDataFrame(random_colors25$cer, random_colors25$par)
#   random_CCMdf35 <- makeCCMDataFrame(random_colors35$cer, random_colors35$par)
#   
#   makeCCMGeneDataFrame <- function(.ccmdf, .cer_colors, .par_colors) {
#     colordf <- tibble(gene_name = colnames(counts_all2$cer), 
#                       cer_color = .cer_colors,
#                       par_color = .par_colors)
#     ccm_genedf <- left_join(genedf, colordf, by = "gene_name")
#     ccm_genedf <- left_join(ccm_genedf, .ccmdf, by = c("cer_color", "par_color"))
#     # because we're wondering if block size can predict divergence at all,
#     # we're looking for the maximum LFC across the 4 experiments
#     ccm_genedf$max_lfc <- map(c(1:nrow(ccm_genedf)), \(i) {
#       return(max(abs(c(ccm_genedf$effect_size_CC_species[i],
#                        ccm_genedf$effect_size_HAP4_species[i],
#                        ccm_genedf$effect_size_LowPi_species[i],
#                        ccm_genedf$effect_size_LowN_species[i]))))
#     }) %>% unlist()
#     return(ccm_genedf)
#   }
#   # tests for makeCCMGeneDataFrame
#   # test_CCM_genedf <- makeCCMGeneDataFrame(test_CCM, test_colors_cer, test_colors_par)
#   
#   CCM_genedf00 <- makeCCMGeneDataFrame(CCMdf00, colors00$cer, colors00$par)
#   CCM_genedf10 <- makeCCMGeneDataFrame(CCMdf10, colors10$cer, colors10$par)
#   CCM_genedf25 <- makeCCMGeneDataFrame(CCMdf25, colors25$cer, colors25$par)
#   CCM_genedf35 <- makeCCMGeneDataFrame(CCMdf35, colors35$cer, colors35$par)
#   
#   
#   random_CCM_genedf00 <- makeCCMGeneDataFrame(random_CCMdf00, random_colors00$cer, random_colors00$par)
#   random_CCM_genedf10 <- makeCCMGeneDataFrame(random_CCMdf10, random_colors10$cer, random_colors10$par)
#   random_CCM_genedf25 <- makeCCMGeneDataFrame(random_CCMdf25, random_colors25$cer, random_colors25$par)
#   random_CCM_genedf35 <- makeCCMGeneDataFrame(random_CCMdf35, random_colors35$cer, random_colors35$par)
#   
#   save(CCMdf00, CCMdf10, CCMdf25, CCMdf35, 
#        CCM_genedf00, CCM_genedf10, CCM_genedf25, CCM_genedf35, 
#        random_CCMdf00, random_CCMdf10, random_CCMdf25, random_CCMdf35,
#        random_CCM_genedf00, random_CCM_genedf10, random_CCM_genedf25, random_CCM_genedf35,
#        file = "data_files/CCMs.RData")
# }
# load(file = "data_files/CCMs.RData")
# 
# # Conservation tests for Leave One Out networks
# load(file = "data_files/LeaveOneOut.RData")
# if (!exists("leaveoneoutdf")) {
#   makeLOODataFrame <- function(.clrs1, .clrs2, .clrs_name1 = "loo_color", .clrs_name2 = "full_color") {
#     each_color1 <- rep(unique(.clrs1), length(unique(.clrs2))) %>% as.character()
#     each_color2 <- sapply(unique(.clrs2), rep, times = length(unique(.clrs1))) %>% as.character()
#     loodf <- tibble(col1 = each_color1,
#                     col2 = each_color2)
#     names(loodf) <- c(.clrs_name1, .clrs_name2)
#     loodf$p_value <- map2(loodf[,.clrs_name1, drop = TRUE], loodf[,.clrs_name2, drop = TRUE], \(x, y) {
#       cat("starting on", x, "and", y, "\n")
#       return(callConservedConnection(.color_name1 = x,
#                                      .color_name2 = y,
#                                      .colors1 = .clrs1,
#                                      .colors2 = .clrs2,
#                                      .n_permutations_per_species = 5000))
#     }) %>% unlist()
#     loodf$block_size <- map2(loodf[,.clrs_name1, drop = TRUE], loodf[,.clrs_name2, drop = TRUE], \(x, y) {
#       return(sum(.clrs1 == x & .clrs2 == y))
#     }) %>% unlist()
#     return(loodf)
#   }
#   leaveoneoutdf <- tibble("experiment" = character(0),
#                           "species" = character(0),
#                           "loo_color" = character(0),
#                           "full_color" = character(0),
#                           "p_value" = numeric(0),
#                           "block_size" = numeric(0))
#   netlist <- list(networks_noLowN, networks_noLowPi, networks_noHAP4, 
#                   networks_noCC, networks_noTemp)
#   names(netlist) <- c("LowN", "LowPi", "HAP4", "CC", "Temp")
#   for (sp in c("cer", "par")) {
#     for (e in c("LowN", "LowPi", "HAP4", "CC", "Temp")) {
#       net <- netlist[[e]]
#       ldf <- makeLOODataFrame(.clrs1 = net[["colors"]][[sp]],
#                               .clrs2 = colors00[[sp]])
#       ldf$experiment <- e
#       ldf$species <- sp
#       leaveoneoutdf <- bind_rows(leaveoneoutdf, ldf)
#     }
#   }
#   save(leaveoneoutdf, networks_noLowN, networks_noHAP4, 
#        networks_noCC, networks_noLowPi, networks_noTemp,
#        file = "data_files/LeaveOneOut.RData")
# }
# #### Figure 3 module dataframe ####
# 
# # for two matrices of counts, calculates
# # the expression similarity between them (as avg expression correlation and eigengene distance)
# # @input: counts1 and counts2 where genes are columns and samples are rows
# # @output: average expression and eigengene vectors for counts 1 and 2
# getExpressionSimilarity <- function(.cts1, .cts2) {
#   scaled_cts1 <- scale(.cts1)
#   scaled_cts2 <- scale(.cts2)
#   avg1 <- scaled_cts1 |> rowMeans(na.rm = TRUE)
#   avg2 <- scaled_cts2 |> rowMeans(na.rm = TRUE)
#   # eigen1 <- NA
#   # eigen2 <- NA
#   # if (.eigenToo) {
#   #   eigen1 <- .cts1[!is.na(.cts1)] |> svd(nu = 1, nv = 1)
#   #   eigen1 <- eigen1$u[,1] # module eigengene is first eigenvector of C*t(C) Scaled count matrix where genes are columns
#   #   eigen2 <- .cts2[!is.na(.cts2)] |> svd(nu= 1, nv = 1)
#   #   eigen2 <- eigen2$u[,1]
#   # }
#   return(list(avg1 = avg1,
#               avg2 = avg2,
#               # eigen1 = eigen1,
#               # eigen2 = eigen2,
#               expr1 = rowMeans(.cts1, na.rm = TRUE), # unscaled module expression, for comparing overall expression levels
#               expr2 = rowMeans(.cts2, na.rm = TRUE)))
# }
# # version using ModuleEigengenes for comparison
# # getExpressionSimilarityWGCNA <- function(.cts1, .cts2) {
# #   MEs1 <- moduleEigengenes(.cts1, colors = rep("blue", ncol(.cts1)), nPC = 1)
# #   MEs2 <- moduleEigengenes(.cts2, colors = rep("blue", ncol(.cts2)), nPC = 1)
# #   return(list(avg1 = MEs1$averageExpr$AEblue,
# #               avg2 = MEs2$averageExpr$AEblue,
# #               eigen1 = MEs1$eigengenes$MEblue,
# #               eigen2 = MEs2$eigengenes$MEblue))
# # }
# # # tests for getExpressionSimilarity
# # # portion of CCM
# # test_idxs <- c("YBR105C","YBR114W","YBR208C","YBR235W","YCL025C","YDL138W")
# # test_cts_cer <- counts_all2$cer[,test_idxs]
# # test_cts_par <- counts_all2$par[,test_idxs]
# # test_inHouse <- getExpressionSimilarity(test_cts_cer, test_cts_par)
# # #test_WGCNA <- getExpressionSimilarityWGCNA(test_cts_cer, test_cts_par)
# # cor(test_inHouse$avg1, test_inHouse$avg2)
# # #cor(test_WGCNA$avg1, test_WGCNA$avg2)
# # sum((test_inHouse$eigen1 - test_inHouse$eigen2)^2) %>% sqrt()
# # #sum((test_WGCNA$eigen1 - test_WGCNA$eigen2)^2) %>% sqrt()
# # # portion of cer-unique
# # test_idx <- c("YPL100W","YPL103C","YPL105C","YPL138C","YPL140C")
# # test_cts_cer <- counts_all2$cer[,test_idx]
# # test_cts_par <- counts_all2$par[,test_idx]
# # test_inHouse <- getExpressionSimilarity(test_cts_cer, test_cts_par)
# # #test_WGCNA <- getExpressionSimilarityWGCNA(test_cts_cer, test_cts_par)
# # cor(test_inHouse$avg1, test_inHouse$avg2)
# # #cor(test_WGCNA$avg1, test_WGCNA$avg2)
# # sum((test_inHouse$eigen1 - test_inHouse$eigen2)^2) %>% sqrt()
# # #sum((test_WGCNA$eigen1 - test_WGCNA$eigen2)^2) %>% sqrt()
# # # random module (single gene)
# # test_idx <- "YIL092W"
# # test_cts_cer <- counts_all2$cer[,test_idx, drop = FALSE]
# # test_cts_par <- counts_all2$par[,test_idx, drop = FALSE] # random modules can be single genes, which by default turns a count matrix into a 1D vector
# # test <- getExpressionSimilarity(test_cts_cer, test_cts_par)
# # cor(test$avg1, test$avg2)
# # sum((test$eigen1 - test$eigen2)^2) %>% sqrt()
# 
# ### moduledf
# if (!file.exists("data_files/modules.RData")) {
#   # matching genes to their modules in each species.
#   # CCMs get letters, divergent modules get numbers for each species (even for cer, odd for par)
#   CCMtestsdf00 <- CCMdf00 # these include pairs that didn't pass
#   CCMtestsdf10 <- CCMdf10 
#   CCMtestsdf25 <- CCMdf25
#   CCMtestsdf35 <- CCMdf35
#   
#   CCMdf00$is_CCM <- CCMdf00$p_value < 0.05/nrow(CCMdf00) & CCMdf00$cer_color != "grey" & CCMdf00$par_color != "grey"
#   CCMdf10$is_CCM <- CCMdf10$p_value < 0.05/nrow(CCMdf10) & CCMdf10$cer_color != "grey" & CCMdf10$par_color != "grey"
#   CCMdf25$is_CCM <- CCMdf25$p_value < 0.05/nrow(CCMdf25) & CCMdf25$cer_color != "grey" & CCMdf25$par_color != "grey"
#   CCMdf35$is_CCM <- CCMdf35$p_value < 0.05/nrow(CCMdf35) & CCMdf35$cer_color != "grey" & CCMdf35$par_color != "grey"
#   
#   random_CCMdf00$is_CCM <- random_CCMdf00$p_value < 0.05/nrow(random_CCMdf00) & CCMdf00$cer_color != "grey" & CCMdf00$par_color != "grey"
#   random_CCMdf10$is_CCM <- random_CCMdf10$p_value < 0.05/nrow(random_CCMdf10) & CCMdf10$cer_color != "grey" & CCMdf10$par_color != "grey"
#   random_CCMdf25$is_CCM <- random_CCMdf25$p_value < 0.05/nrow(random_CCMdf25) & CCMdf25$cer_color != "grey" & CCMdf25$par_color != "grey"
#   random_CCMdf35$is_CCM <- random_CCMdf35$p_value < 0.05/nrow(random_CCMdf35) & CCMdf35$cer_color != "grey" & CCMdf35$par_color != "grey"
#   
#   # making the module dataframes that record which cer/par color combination 
#   # corresponds to which module name (and random module name) and color (if CCM)
#   # this is hefty, so we need to create a few helper functions for the different steps
#   
#   # @input: a CCMdf enumerating the color in cer and par of each module (cer/par color combination)
#   # @output: moduledf giving module name and CCM color for each cer/par color combination
#   assignModuleNameAndColor <- function(.ccmdf, .ccm_genedf) {
#     # filtering for cer/par color combinations actually present in the data 
#     # (CCMdf has all combinations, CCM_genedf only has ones that actually exist)
#     mdf <- .ccmdf |> 
#       filter(paste(cer_color, par_color) %in% paste(.ccm_genedf$cer_color, .ccm_genedf$par_color))
#     
#     # assigning module coexpressed states: conserved co-expressed,
#     # diverged co-expressed, S. cerevisiae co-expressed, 
#     # S. paradoxus co-expressed, not co-expressed)
#     mdf$coexpressed <- apply(mdf, 1, \(x) {
#       cercol <- x["cer_color"] |> as.character()
#       parcol <- x["par_color"] |> as.character()
#       isccm <- x["is_CCM"] |> as.logical()
#       if (cercol == "grey" & parcol == "grey") {
#         return("never co-expressed")
#       }
#       if (cercol == "grey" & parcol != "grey") {
#         return("S. paradoxus co-expressed")
#       }
#       if (cercol != "grey" & parcol == "grey") {
#         return("S. cerevisiae co-expressed")
#       }
#       if (cercol != "grey" & parcol != "grey" & isccm) {
#         return("conserved co-expressed")
#       }
#       if (cercol != "grey" & parcol != "grey" & !isccm) {
#         return("diverged co-expressed")
#       }
#     }) |> unlist()
#     # assigning module names
#     module_name_vector <- vector(mode = "character", length = 0)
#     for (i in c(1:nrow(mdf))) {
#       coexpressed <- mdf[i, "coexpressed"] |> as.character()
#       if (coexpressed == "conserved co-expressed") {
#         module_name_vector <- c(module_name_vector, setdiff(letters, module_name_vector)[1])
#       }
#       if (coexpressed == "diverged co-expressed") {
#         module_name_vector <- c(module_name_vector, paste0("div", i))
#       }
#       if (coexpressed == "S. cerevisiae co-expressed") {
#         module_name_vector <- c(module_name_vector, paste0("cer", i))
#       }
#       if (coexpressed == "S. paradoxus co-expressed") {
#         module_name_vector <- c(module_name_vector, paste0("par", i))
#       }
#       if (coexpressed == "never co-expressed") {
#         module_name_vector <- c(module_name_vector, paste0("nev", i))
#       }
#     }
#     mdf$module_name <- module_name_vector
#     mdf$coexpressed <- factor(mdf$coexpressed, levels = c("conserved co-expressed",
#                                                           "diverged co-expressed",
#                                                           "S. cerevisiae co-expressed",
#                                                           "S. paradoxus co-expressed",
#                                                           "never co-expressed"))
#     # assigning colors to CCMs
#     # the rationale is that as much as possible, keep them consistent with the single-species colors
#     # first, straightforward scenarios:
#     # if a CCM is blue in both species, it should be blue
#     # if a CCM is grey in either species, it should not have a CCM color, because grey means it's not a CCM
#     mdf$CCM_color <- map(c(1:nrow(mdf)), \(i) {
#       x <- mdf$cer_color[i]
#       y <- mdf$par_color[i]
#       isccm <- mdf$is_CCM[i]
#       if (!isccm) {
#         return("none")
#       }
#       if (x == "grey" | y == "grey") {
#         return(NA)
#       }
#       if (x == y) {
#         return(x)
#       }
#       if (x != y) {
#         return("needs_color")
#       }
#     }) |> unlist()
#     # next we tackle the less-straightforward scenarios:
#     # if a CCM is purple in one species and yellow in another, but there's a
#     # different CCM that's purple in both species, then this CCM should be yellow
#     # (unless there's also a CCM that's yellow in both)
#     # you see why this is challenging
#     for (i in 1:nrow(mdf)) {
#       ccm_color <- mdf[i, "CCM_color"] |> as.character()
#       if (is.na(ccm_color)) {
#         next()
#       }
#       if (ccm_color != "needs_color") {
#         next()
#       }
#       if (ccm_color == "needs_color") {
#         cer_col <- mdf[i,"cer_color"] |> as.character()
#         par_col <- mdf[i,"par_color"] |> as.character()
#         if (cer_col %in% mdf$CCM_color) {
#           if (par_col %in% mdf$CCM_color) {
#             new_color <- setdiff(standardColors(sum(.ccmdf$is_CCM) + 2), 
#                                  c(mdf$CCM_color, 
#                                    "grey60", "cyan", "black", "brown", "tan")) |> # all the colors we want to exclude b/c they don't look nice (or are too close to other standard colors)
#               sample(1)
#           }
#           if (!(par_col %in% mdf$CCM_color)) {
#             new_color <- par_col
#           }
#         }
#         if (!(cer_col %in% mdf$CCM_color)) {
#           new_color <- cer_col
#         }
#         mdf[i, "CCM_color"] <- new_color
#       }
#     }
#     return(mdf)
#   }
#   # # tests for assignModuleNameAndColor
#   # assignModuleNameAndColor(CCMdf00[c(1,2,3),])
#   # assignModuleNameAndColor(CCMdf[c(1,2,3),])
#   # assignModuleNameAndColor(CCMdf25[c(1,2,3),])
#   # test_mdf <- assignModuleNameAndColor(CCMdf00, CCM_genedf00)
#   # test_m_genedf <-  left_join(CCM_genedf00, test_mdf, by = c("cer_color", "par_color"))
#   
#   # For each cer/par color combination, we collect the names of genes in that color combination
#   # then calculate expression similarity of those genes in each environment
#   addExpressionSimilarityToModuleDataFrame <- function(.mdf, .m_genedf) {
#     module_results <- map(.mdf$module_name, \(m) {
#       cat("working on module", m, "\n")
#       gene_idxs <- .m_genedf |> filter(module_name == m) |> select(gene_name) |> pull()
#       cts_cer <- counts_all2$cer[,gene_idxs, drop = FALSE]
#       cts_par <- counts_all2$par[,gene_idxs, drop = FALSE]
#       return(getExpressionSimilarity(cts_cer, cts_par))
#     })
#     names(module_results) <- .mdf$module_name
#     
#     module_results_by_experiment <- map2_dfr(module_results, names(module_results), \(x, y) {
#       avgCor <- cor(x$avg1, x$avg2, use = "pairwise.complete.obs")
#       avgCor_CC <- cor(x$avg1[info$experiment == "CC"], x$avg2[info$experiment == "CC"], use = "pairwise.complete.obs")
#       avgCor_LowN <- cor(x$avg1[info$experiment == "LowN"], x$avg2[info$experiment == "LowN"], use = "pairwise.complete.obs")
#       avgCor_LowPi <- cor(x$avg1[info$experiment == "LowPi"], x$avg2[info$experiment == "LowPi"], use = "pairwise.complete.obs")
#       avgCor_HAP4 <- cor(x$avg1[info$experiment == "HAP4"], x$avg2[info$experiment == "HAP4"], use = "pairwise.complete.obs")
#       avgCor_Heat <- cor(x$avg1[info$experiment == "Heat"], x$avg2[info$experiment == "Heat"], use = "pairwise.complete.obs")
#       avgCor_Cold <- cor(x$avg1[info$experiment == "Cold"], x$avg2[info$experiment == "Cold"], use = "pairwise.complete.obs")
#       # dcor_overall <- dcor(x$avg1, x$avg2)
#       # dcor_CC <- dcor(x$avg1[info$experiment == "CC"], x$avg2[info$experiment == "CC"])
#       # dcor_LowN <- dcor(x$avg1[info$experiment == "LowN"], x$avg2[info$experiment == "LowN"])
#       # dcor_LowPi <- dcor(x$avg1[info$experiment == "LowPi"], x$avg2[info$experiment == "LowPi"])
#       # dcor_HAP4 <- dcor(x$avg1[info$experiment == "HAP4"], x$avg2[info$experiment == "HAP4"])
#       # eigenDist <- sum((x$eigen1 - x$eigen2)^2) |> sqrt()
#       # eigenDist_CC <- sum((x$eigen1[info$experiment == "CC"] - x$eigen2[info$experiment == "CC"])^2) %>% sqrt()
#       # eigenDist_LowN <- sum((x$eigen1[info$experiment == "LowN"] - x$eigen2[info$experiment == "LowN"])^2) %>% sqrt()
#       # eigenDist_LowPi <- sum((x$eigen1[info$experiment == "LowPi"] - x$eigen2[info$experiment == "LowPi"])^2) %>% sqrt()
#       # eigenDist_HAP4 <- sum((x$eigen1[info$experiment == "HAP4"] - x$eigen2[info$experiment == "HAP4"])^2) %>% sqrt()
#       # abc <- sum(abs(x$avg1 - x$avg2))/nrow(info)
#       # abc_CC <- sum(abs(x$avg1[info$experiment == "CC"] - x$avg2[info$experiment == "CC"]))/sum(info$experiment == "CC") # dividing by number of observations (timepoints) per experiment otherwise the experiments with more timepoints will naturally have higher ABC
#       # abc_LowN <- sum(abs(x$avg1[info$experiment == "LowN"] - x$avg2[info$experiment == "LowN"]))/sum(info$experiment == "LowN")
#       # abc_LowPi <- sum(abs(x$avg1[info$experiment == "LowPi"] - x$avg2[info$experiment == "LowPi"]))/sum(info$experiment == "LowPi")
#       # abc_HAP4 <- sum(abs(x$avg1[info$experiment == "HAP4"] - x$avg2[info$experiment == "HAP4"]))/sum(info$experiment == "HAP4")
#       # mean_CC <- mean(c(x$expr1[info$experiment == "CC"], x$expr2[info$experiment == "CC"])) # mean including both comparison groups, not scaled
#       # mean_LowN <- mean(c(x$expr1[info$experiment == "LowN"], x$expr2[info$experiment == "LowN"]))
#       # mean_LowPi <- mean(c(x$expr1[info$experiment == "LowPi"], x$expr2[info$experiment == "LowPi"]))
#       # mean_HAP4 <- mean(c(x$expr1[info$experiment == "HAP4"], x$expr2[info$experiment == "HAP4"]))
#       return(list(module_name = y,
#                   avgCor = avgCor,
#                   avgCor_CC = avgCor_CC,
#                   avgCor_LowN = avgCor_LowN,
#                   avgCor_LowPi = avgCor_LowPi,
#                   avgCor_HAP4 = avgCor_HAP4,
#                   avgCor_Heat = avgCor_Heat,
#                   avgCor_Cold = avgCor_Cold
#                   # dcor = dcor_overall,
#                   # dcor_CC = dcor_CC,
#                   # dcor_LowN = dcor_LowN,
#                   # dcor_LowPi = dcor_LowPi,
#                   # dcor_HAP4 = dcor_HAP4,
#                   # eigenDist = eigenDist,
#                   # eigenDist_CC = eigenDist_CC,
#                   # eigenDist_LowN = eigenDist_LowN,
#                   # eigenDist_LowPi = eigenDist_LowPi,
#                   # eigenDist_HAP4 = eigenDist_HAP4,
#                   # abc = abc,
#                   # abc_CC = abc_CC,
#                   # abc_LowN = abc_LowN,
#                   # abc_LowPi = abc_LowPi,
#                   # abc_HAP4 = abc_HAP4,
#                   # mean_CC = mean_CC,
#                   # mean_LowN = mean_LowN,
#                   # mean_LowPi = mean_LowPi,
#                   # mean_HAP4 = mean_HAP4
#       ))
#     })
#     .mdf <- left_join(.mdf, module_results_by_experiment, by = "module_name")
#     # .mdf$mean <- pmean(.mdf$mean_CC,
#     #                    .mdf$mean_LowN,
#     #                    .mdf$mean_LowPi,
#     #                    .mdf$mean_HAP4)
#     
#     return(.mdf)
#   }
#   # tests for addExpressionSimilarityToModuleDataFrame
#   # test_mdf <- addExpressionSimilarityToModuleDataFrame(test_mdf, test_m_genedf)
#   
#   makeModuleDataFrame <- function(.ccmdf, .ccm_genedf, .colors_cer, .colors_par) {
#     mdf <- assignModuleNameAndColor(.ccmdf, .ccm_genedf)
#     # mdf_divergent <- assignDivergentModuleName(mdf, .ccmdf, .ccm_genedf, .colors_cer, .colors_par)
#     
#     # mdf <- bind_rows(mdf, mdf_divergent)
#     
#     # pairing each gene with its module name (necessary for calculating expression similarity below)
#     m_genedf <- left_join(.ccm_genedf, mdf, by = c("cer_color", "par_color"))
#     
#     mdf <- addExpressionSimilarityToModuleDataFrame(mdf, m_genedf)
#     
#     return(mdf)
#   }
#   # tests for makeModuleDataFrame
#   # test_mdf <- makeModuleDataFrame(CCMdf10, CCM_genedf10, colors10$cer, colors10$par)
#   
#   # making module data frames, each row is module
#   
#   moduledf00 <- makeModuleDataFrame(CCMdf00, CCM_genedf00, colors00$cer, colors00$par)
#   moduledf10 <- makeModuleDataFrame(CCMdf10, CCM_genedf10, colors10$cer, colors10$par)
#   moduledf25 <- makeModuleDataFrame(CCMdf25, CCM_genedf25, colors25$cer, colors25$par)
#   moduledf35 <- makeModuleDataFrame(CCMdf35, CCM_genedf35, colors35$cer, colors35$par)
#   random_moduledf00 <- makeModuleDataFrame(random_CCMdf00, random_CCM_genedf00, random_colors00$cer, random_colors00$par)
#   random_moduledf10 <- makeModuleDataFrame(random_CCMdf10, random_CCM_genedf10, random_colors10$cer, random_colors10$par)
#   random_moduledf25 <- makeModuleDataFrame(random_CCMdf25, random_CCM_genedf25, random_colors25$cer, random_colors25$par)
#   random_moduledf35 <- makeModuleDataFrame(random_CCMdf35, random_CCM_genedf35, random_colors35$cer, random_colors35$par)
#   # changing random module names and co-expression states, because they're meaningless and confusing to retain
#   random_moduledf00$module_name <- paste0("ran", c(1:nrow(random_moduledf00)))
#   random_moduledf10$module_name <- paste0("ran", c(1:nrow(random_moduledf10)))
#   random_moduledf25$module_name <- paste0("ran", c(1:nrow(random_moduledf25)))
#   random_moduledf35$module_name <- paste0("ran", c(1:nrow(random_moduledf35)))
#   random_moduledf00$coexpressed <- "random"
#   random_moduledf10$coexpressed <- "random"
#   random_moduledf25$coexpressed <- "random"
#   random_moduledf35$coexpressed <- "random"
#   # changing any additional colors from WGCNA that we don't like for publication modules
#   if (!("turquoise" %in% moduledf25$CCM_color)) {
#     moduledf25$CCM_color <- gsub("brown", "turquoise", moduledf25$CCM_color)
#   }
#   if (!("orange" %in% moduledf25$CCM_color)) {
#     moduledf25$CCM_color <- gsub("black", "orange", moduledf25$CCM_color)
#   }
#   
#   # making module gene-wise data frames, each row is gene (different from CCM_genedf because now there's module expression similarity information)
#   module_genedf00 <- left_join(CCM_genedf00, moduledf00, by = c("cer_color", "par_color", "block_size"))
#   module_genedf10 <- left_join(CCM_genedf10, moduledf10, by = c("cer_color", "par_color", "block_size"))
#   module_genedf25 <- left_join(CCM_genedf25, moduledf25, by = c("cer_color", "par_color", "block_size"))
#   module_genedf35 <- left_join(CCM_genedf35, moduledf35, by = c("cer_color", "par_color", "block_size"))
#   random_module_genedf00 <- left_join(random_CCM_genedf00, random_moduledf00, by = c("cer_color", "par_color", "block_size"))
#   random_module_genedf10 <- left_join(random_CCM_genedf10, random_moduledf10, by = c("cer_color", "par_color", "block_size"))
#   random_module_genedf25 <- left_join(random_CCM_genedf25, random_moduledf25, by = c("cer_color", "par_color", "block_size"))
#   random_module_genedf35 <- left_join(random_CCM_genedf35, random_moduledf35, by = c("cer_color", "par_color", "block_size"))
#   
#   
#   save(moduledf00, module_genedf00, moduledf10, module_genedf10, 
#        moduledf25, module_genedf25, moduledf35, module_genedf35,
#        random_moduledf00, random_module_genedf00, random_moduledf10, random_module_genedf10, 
#        random_moduledf25, random_module_genedf25, random_moduledf35, random_module_genedf35,
#        file = "data_files/modules.RData")
# }
# load("data_files/modules.RData")

# #### Figure 3: plotting expression profiles and adding allelic info to moduledf ####
# 
# # now adding the expression similarity between cer-par for each module (both CCM and diverged)
# if (!exists("moduledf_allele25")) {
#   addAllelicExpressionSimilarityToModuleDataFrame <- function(.mdf, .m_genedf) {
#     module_results_allele <- map(.mdf$module_name, \(m) {
#       gene_idxs <- .m_genedf |> filter(module_name == m) |> select(gene_name) |> pull()
#       cts_cer <- counts_all2_allele$cer[,gene_idxs, drop = FALSE]
#       cts_par <- counts_all2_allele$par[,gene_idxs, drop = FALSE]
#       return(getExpressionSimilarity(cts_cer, cts_par))
#     })
#     module_results_by_experiment_allele <- map_dfr(module_results_allele, \(x) {
#       cor_overall <- cor(x$avg1, x$avg2, use = "pairwise.complete.obs")
#       cor_CC <- cor(x$avg1[info$experiment == "CC"], x$avg2[info$experiment == "CC"], use = "pairwise.complete.obs")
#       cor_LowN <- cor(x$avg1[info$experiment == "LowN"], x$avg2[info$experiment == "LowN"], use = "pairwise.complete.obs")
#       cor_LowPi <- cor(x$avg1[info$experiment == "LowPi"], x$avg2[info$experiment == "LowPi"], use = "pairwise.complete.obs")
#       cor_HAP4 <- cor(x$avg1[info$experiment == "HAP4"], x$avg2[info$experiment == "HAP4"], use = "pairwise.complete.obs")
#       cor_Heat <- cor(x$avg1[info$experiment == "Heat"], x$avg2[info$experiment == "Heat"], use = "pairwise.complete.obs")
#       cor_Cold <- cor(x$avg1[info$experiment == "Cold"], x$avg2[info$experiment == "Cold"], use = "pairwise.complete.obs")
#       return(list(avgCor = cor_overall,
#                   avgCor_CC = cor_CC,
#                   avgCor_LowN = cor_LowN,
#                   avgCor_LowPi = cor_LowPi,
#                   avgCor_HAP4 = cor_HAP4,
#                   avgCor_Heat = cor_Heat,
#                   avgCor_Cold = cor_Cold))
#     })
#     moduledf_allele <- .mdf %>% select(module_name, is_CCM, CCM_color, coexpressed,
#                                        cer_color, par_color, block_size) %>% 
#       bind_cols(module_results_by_experiment_allele)
#     # moduledf_allele$cor <- pmean(moduledf_allele$cor_CC,
#     #                              moduledf_allele$cor_LowN,
#     #                              moduledf_allele$cor_LowPi,
#     #                              moduledf_allele$cor_HAP4,
#     #                              moduledf_allele$cor_Heat,
#     #                              moduledf_allele$cor_Cold)
#     return(moduledf_allele)
#   }
#   moduledf_allele00 <- addAllelicExpressionSimilarityToModuleDataFrame(moduledf00, module_genedf00)
#   moduledf_allele10 <- addAllelicExpressionSimilarityToModuleDataFrame(moduledf10, module_genedf10)
#   moduledf_allele25 <- addAllelicExpressionSimilarityToModuleDataFrame(moduledf25, module_genedf25)
#   moduledf_allele35 <- addAllelicExpressionSimilarityToModuleDataFrame(moduledf35, module_genedf35)
#   random_moduledf_allele00 <- addAllelicExpressionSimilarityToModuleDataFrame(random_moduledf00, random_module_genedf00)
#   random_moduledf_allele10 <- addAllelicExpressionSimilarityToModuleDataFrame(random_moduledf10, random_module_genedf10)
#   random_moduledf_allele25 <- addAllelicExpressionSimilarityToModuleDataFrame(random_moduledf25, random_module_genedf25)
#   random_moduledf_allele35 <- addAllelicExpressionSimilarityToModuleDataFrame(random_moduledf35, random_module_genedf35)
#   
#   save(moduledf00, module_genedf00, moduledf10, module_genedf10, 
#        moduledf25, module_genedf25, moduledf35, module_genedf35,
#        random_moduledf00, random_module_genedf00, random_moduledf10, random_module_genedf10, 
#        random_moduledf25, random_module_genedf25, random_moduledf35, random_module_genedf35, 
#        moduledf_allele00, moduledf_allele10, moduledf_allele25, moduledf_allele35,
#        random_moduledf_allele00, random_moduledf_allele10, random_moduledf_allele25, random_moduledf_allele35,
#        file = "data_files/modules.RData")
# }


############### Other archive ############### 
#### Figure 1/2 plotGenePair ####
# archived b/c plotExpressonProfilePair can do everything this can
# plotGenePair <- function(.gdf, .gene_pair_names = c("Gene 1", "Gene 2")) {
#   pearson_cor <- cor(.gdf$expr_1, .gdf$expr_2)
#   plotdf <- .gdf %>% pivot_longer(cols = c(expr_1, expr_2), names_to = "gene_id", values_to = "expr")
#   plotdf$sample_order <- frank(plotdf, experiment, time_point_num, ties.method = "first")
#   plotdf$experiment_plus <- map(c(1:nrow(plotdf)), \(i) {
#     if (plotdf$experiment[i] == "LowN")
#       return(paste(plotdf$experiment[i], plotdf$time_point_num[i]))
#     else
#       return(plotdf$experiment[i])
#   }) %>% unlist()
#   experiment_widths <- as.numeric(table(sort(plotdf$experiment_plus)))
#   experiment_names <- names(table(sort(plotdf$experiment_plus)))
#   experiment_breaks <- cumsum(experiment_widths)
#   rects <- tibble(xstart = c(0, experiment_breaks[-length(experiment_breaks)]),
#                   xend = experiment_breaks,
#                   experiment = experiment_names)
#   p <- ggplot() + 
#     geom_rect(data = rects, aes(ymin = 0, ymax = max(log2(plotdf$expr + 1)),
#                                 xmin = xstart, xmax = xend, fill = experiment_names), alpha = 0.5) +
#     geom_line(data = plotdf, aes(x = sample_order, y = log2(expr + 1), color = gene_id)) + 
#     theme_classic() +
#     scale_fill_discrete(type = c("orchid", "lightgreen","yellow2", "gold", "darkgoldenrod2", "lightblue"), 
#                         labels = c("Urea Shock", "No Transfers", "Low Nitrogen time 1", "Low Nitrogen time 2", "Low Nitrogen time 3", "Low Phosphorus")) +
#     scale_color_discrete(type = c(.color1, .color2), labels = .gene_pair_names) +
#     theme(legend.title = element_blank()) +
#     ylab("Expression (log2)") +
#     xlab("") +
#     scale_x_discrete(labels = element_blank()) +
#     ggtitle(paste("Pearson correlation = ", round(pearson_cor, 2)))
#   return(p)
# }
# # # tests for plotGenePair
# # test_gene_pair <- c("YJR010W", "YIL074C") # MET3/SER33
# # test_gdf <- tibble(expr_1 = counts_all2$cer[,test_gene_pair[1]] %>% as.numeric(),
# #                    expr_2 = counts_all2$cer[,test_gene_pair[2]] %>% as.numeric()) %>%
# #   bind_cols(info)
# # plotGenePair(test_gdf, test_gene_pair)
#### Stacked plotExpressionProfilePair/Quartet ####
# # plotting function to visualize expression profiles of any 2 groups
# # stacked over common x axis timepoint
# # @input: counts, info, and names of two groups to compare
# # @output: ggplot of both expression profiles with loess or line curves tracing the average expression
# plotExpressionProfilePairStack <- function(.cts1, .cts2, .info1, .info2, 
#                                            .name1 = "S. cerevisiae", .name2 = "S. paradoxus",
#                                            .color1 = "orange1", .color2 = "blue2",
#                                            .method = c("line", "loess"), 
#                                            .show_points = FALSE,
#                                            .show_confidence_intervals = TRUE,
#                                            .normalization = c("none", "log2", "scale")) {
#   if (.normalization == "none") {
#     norm_func <- identity
#     ylabel <- "Expression"
#   }
#   if (.normalization == "log2") {
#     norm_func <- \(x) {log2(x + 1)}
#     ylabel <- "Expression (log2)"
#   }
#   if (.normalization == "scale") {
#     norm_func <- scale
#     ylabel <- "Expression (centered and scaled)"
#   }
#   ExperimentNames <- unique(.info1$experiment) # arbitrary to do info1 or info2
#   nExperiments <- length(ExperimentNames)
#   info1 <- tibble(experiment = .info1$experiment,
#                   time_point_num = .info1$time_point_num)
#   info2 <- tibble(experiment = .info2$experiment,
#                   time_point_num = .info2$time_point_num)
#   # scaling timepoint so that each experiment starts at 0 and ends at 1
#   info1 <- info1 |> group_by(experiment) |> 
#     summarise(max_tp = max(time_point_num)) |> 
#     ungroup() |> right_join(y = info1, by = c("experiment"))
#   info1$scaled_tp <- info1$time_point_num/info1$max_tp
#   info2 <- info2 |> group_by(experiment) |> 
#     summarise(max_tp = max(time_point_num)) |> 
#     ungroup() |> right_join(y = info2, by = c("experiment"))
#   info2$scaled_tp <- info2$time_point_num/info2$max_tp
#   # normalizing expression
#   expr1 <- norm_func(.cts1)
#   expr2 <- norm_func(.cts2)
#   gdf1 <- bind_cols(expr1, info1) |> 
#     pivot_longer(cols = colnames(expr1), names_to = "gene_name", values_to = "expr")
#   gdf1$group_id <- "1"
#   gdf2 <- bind_cols(expr2, info2) |> 
#     pivot_longer(cols = colnames(expr2), names_to = "gene_name", values_to = "expr")
#   gdf2$group_id <- "2"
#   # converting each gene's expression to its mean expression between replicates
#   gdf <- bind_rows(gdf1, gdf2) |> 
#     group_by(group_id, gene_name, experiment, scaled_tp) |> 
#     summarise(expr = mean(expr)) |> ungroup()
#   # background color rectangles for differentiating the 4 experiments
#   rects <- data.frame(color = c("orchid", "lightgreen", "gold", "lightblue"),
#                       labels = c("Urea Shock", "Diauxic Shift", "Low Nitrogen", "Low Phosphorus"),
#                       experiment_names = c("CC", "HAP4", "LowN", "LowPi"))
#   rownames(rects) <- c("CC", "HAP4", "LowN", "LowPi")
#   rects <- rects[ExperimentNames,] # ensuring experiments are in same order as input dataset
#   min_expr <- min(gdf$expr)
#   max_expr <- max(gdf$expr)
#   # plotting
#   makeSingleExperimentPlot <- function(.plotdf, .experiment = c("CC", "HAP4", "LowN", "LowPi")) {
#     #background_color <- rects |> filter(experiment_names == .experiment) |> select(color) |> pull()
#     plot_title <- rects |> filter(experiment_names == .experiment) |> select(labels) |> pull()
#     rects <- filter(rects, experiment_names == .experiment)
#     .plotdf <- filter(.plotdf, experiment == .experiment)
#     p <- ggplot() + 
#       geom_rect(data = rects, aes(ymin = min_expr, ymax = max_expr, 
#                                   xmin = 0, xmax = 1, 
#                                   fill = experiment_names), alpha = 0.3) +
#       scale_fill_discrete(type = rects$color, labels = rects$labels, limits = rects$experiment_names) +
#       theme_classic() +
#       scale_color_discrete(type = c(.color1, .color2), labels = c(.name1, .name2)) +
#       #theme(legend.title = element_blank()) +
#       theme(legend.position = "none") +
#       ylab(ylabel) +
#       ymin(c(min_expr, max_expr)) +
#       xlab("Timepoint") +
#       theme(axis.ticks = element_blank(),
#             axis.text = element_blank()) +
#       ggtitle(plot_title)
#     if (.show_points) {
#       p <- p + geom_jitter(data = .plotdf, aes(x = scaled_tp, y = expr, color = group_id), size = 0.1, alpha = 0.5)
#     }
#     if (.method == "loess") {
#       loess1 <- loess.sd(x = .plotdf %>% filter(experiment == .experiment & group_id == 1) %>% select(scaled_tp) %>% pull(),
#                          y = .plotdf %>% filter(experiment == .experiment& group_id == 1) %>% select(expr) %>% pull(), nsigma = 1.96)
#       loess2 <- loess.sd(x = .plotdf %>% filter(experiment == .experiment& group_id == 2) %>% select(scaled_tp) %>% pull(),
#                          y = .plotdf %>% filter(experiment == .experiment& group_id == 2) %>% select(expr) %>% pull(), nsigma = 1.96)
#       # adding loess segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
#       p <- p + 
#         geom_smooth(aes(x =!!loess1$x, y =!!loess1$y), color = .color1, linewidth = 1) +
#         geom_smooth(aes(x =!!loess2$x, y =!!loess2$y), color = .color2, linewidth = 1)
#       if (.show_confidence_intervals) {
#         p <- p + 
#           geom_ribbon(aes(x =!!loess1$x, ymin =!!loess1$lower, ymax =!!loess1$upper), fill = .color1, alpha = 0.3) +
#           geom_ribbon(aes(x =!!loess2$x, ymin =!!loess2$lower, ymax =!!loess2$upper), fill = .color2, alpha = 0.3)
#       }
#     }
#     if (.method == "line") {
#       # lines trace average expr at each timepoint/experiment for each group
#       avgexpr1 <- .plotdf %>% filter(experiment == .experiment & group_id == 1) |> group_by(scaled_tp) |> summarise(mean_expr = mean(expr),
#                                                                                                                     sd_expr = sd(expr))
#       avgexpr2 <- .plotdf %>% filter(experiment == .experiment & group_id == 2) |> group_by(scaled_tp) |> summarise(mean_expr = mean(expr),
#                                                                                                                     sd_expr = sd(expr))
#       # adding line segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
#       p <- p + 
#         geom_line(data = avgexpr1, aes(x = scaled_tp, y = mean_expr), color = .color1, linewidth = 1) +
#         geom_line(data = avgexpr2, aes(x = scaled_tp, y = mean_expr), color = .color2, linewidth = 1)
#       if (.show_confidence_intervals) {
#         p <- p + 
#           geom_ribbon(data = avgexpr1, aes(x = scaled_tp, ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
#                       fill = .color1, alpha = 0.3) +
#           geom_ribbon(data = avgexpr2, aes(x = scaled_tp, ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
#                       fill = .color2, alpha = 0.3)
#       }
#     }
#     return(p)
#   }
#   plotList <- vector(mode = "list", length = length(unique(.info1$experiment)))
#   names(plotList) <- unique(.info1$experiment)
#   for (e in ExperimentNames) {
#     plotList[[e]] <- makeSingleExperimentPlot(gdf, .experiment = e)
#   }
#   return(ggarrange(plotlist = plotList, nrow = 4, ncol = 1, common.legend = FALSE))
# }
# 
# # tests for plotExpressionProfilePair
# # random_CCM_color <- moduledf25 |> filter(is_CCM & !is.na(CCM_color)) |> select(CCM_color) |> pull() |> sample(size = 1)
# # gene_idxs <-  module_genedf25  |> filter(CCM_color == random_CCM_color) |> select(gene_name) |> pull()
# # plotExpressionProfilePair(counts_all2$cer[info$experiment != "LowPi",gene_idxs],
# #                           counts_all2$par[info$experiment != "LowPi",gene_idxs],
# #                           info[info$experiment != "LowPi",],
# #                           info[info$experiment != "LowPi",],
# #                           .name1 = "S. cereviaise",
# #                           .name2 = "S. paradoxus",
# #                           .method = "line", .show_points = TRUE,
# #                           .normalization = "log2") + ggtitle(random_CCM_color)
# # p <- plotExpressionProfilePairStack(counts_all2$cer[info$experiment != "LowPi",gene_idxs],
# #                           counts_all2$par[info$experiment != "LowPi",gene_idxs],
# #                           info[info$experiment != "LowPi",],
# #                           info[info$experiment != "LowPi",],
# #                           .name1 = "S. cereviaise",
# #                           .name2 = "S. paradoxus",
# #                           .method = "line", 
# #                           .show_points = TRUE,
# #                           .show_confidence_intervals = FALSE,
# #                           .normalization = "log2")
# # annotate_figure(p, top = random_CCM_color)
# # 
# plotExpressionProfileQuartetStack <- function(.cts1, .cts2, .cts3, .cts4,
#                                               .info1, .info2, .info3, .info4,
#                                               .name1, .name2, .name3, .name4,
#                                               .color1, .color2, .color3, .color4,
#                                               .normalization = c("none", "log2", "scale")) {
#   
#   if (.normalization == "none") {
#     norm_func <- identity
#     ylabel <- "Expression"
#   }
#   if (.normalization == "log2") {
#     norm_func <- \(x) {log2(x + 1)}
#     ylabel <- "Expression (log2)"
#   }
#   if (.normalization == "scale") {
#     norm_func <- scale
#     ylabel <- "Expression (centered and scaled)"
#   }
#   ExperimentNames <- unique(.info3$experiment) # arbitrary to do info3 or info4
#   nExperiments <- length(ExperimentNames)
#   info1 <- tibble(experiment = .info1$experiment,
#                   time_point_num = .info1$time_point_num)
#   info2 <- tibble(experiment = .info2$experiment,
#                   time_point_num = .info2$time_point_num)
#   info3 <- tibble(experiment = .info3$experiment,
#                   time_point_num = .info3$time_point_num)
#   info4 <- tibble(experiment = .info4$experiment,
#                   time_point_num = .info4$time_point_num)
#   # scaling timepoint in [0,1] for each experiment
#   # scaling timepoint so that each experiment starts at 0 and ends at 1
#   info1 <- info1 |> group_by(experiment) |> 
#     summarise(max_tp = max(time_point_num)) |> 
#     ungroup() |> right_join(y = info1, by = c("experiment"))
#   info1$scaled_tp <- info1$time_point_num/info1$max_tp
#   info2 <- info2 |> group_by(experiment) |> 
#     summarise(max_tp = max(time_point_num)) |> 
#     ungroup() |> right_join(y = info2, by = c("experiment"))
#   info2$scaled_tp <- info2$time_point_num/info2$max_tp
#   info3 <- info3 |> group_by(experiment) |> 
#     summarise(max_tp = max(time_point_num)) |> 
#     ungroup() |> right_join(y = info3, by = c("experiment"))
#   info3$scaled_tp <- info3$time_point_num/info3$max_tp
#   info4 <- info4 |> group_by(experiment) |> 
#     summarise(max_tp = max(time_point_num)) |> 
#     ungroup() |> right_join(y = info4, by = c("experiment"))
#   info4$scaled_tp <- info4$time_point_num/info4$max_tp
#   # adding expression data
#   expr1 <- norm_func(.cts1)
#   expr2 <- norm_func(.cts2)
#   expr3 <- norm_func(.cts3)
#   expr4 <- norm_func(.cts4)
#   gdf1 <- bind_cols(expr1, info1) |> 
#     pivot_longer(cols = colnames(expr1), names_to = "gene_name", values_to = "expr")
#   gdf1$group_id <- "1"
#   gdf2 <- bind_cols(expr2, info2) |> 
#     pivot_longer(cols = colnames(expr2), names_to = "gene_name", values_to = "expr")
#   gdf2$group_id <- "2"
#   gdf3 <- bind_cols(expr3, info3) |> 
#     pivot_longer(cols = colnames(expr3), names_to = "gene_name", values_to = "expr")
#   gdf3$group_id <- "3"
#   gdf4 <- bind_cols(expr4, info4) |> 
#     pivot_longer(cols = colnames(expr4), names_to = "gene_name", values_to = "expr")
#   gdf4$group_id <- "4"
#   # converting each gene's expression to its mean expression between replicates
#   gdf <- bind_rows(gdf1, gdf2, gdf3, gdf4) |> 
#     group_by(group_id, gene_name, experiment, time_point_num, scaled_tp) |> 
#     summarise(expr = mean(expr)) |> ungroup()
#   # background color rectangles for differentiating the 4 experiments
#   rects <- data.frame(color = c("orchid", "lightgreen", "gold", "lightblue"),
#                       labels = c("Urea Shock", "Diauxic Shift", "Low Nitrogen", "Low Phosphorus"),
#                       experiment_names = c("CC", "HAP4", "LowN", "LowPi"))
#   rownames(rects) <- c("CC", "HAP4", "LowN", "LowPi")
#   rects <- rects[ExperimentNames,] # ensuring experiments are in same order as input dataset
#   min_expr <- min(gdf$expr)
#   max_expr <- max(gdf$expr)
#   # plotting
#   makeSingleExperimentPlot <- function(.plotdf, .experiment = c("CC", "HAP4", "LowN", "LowPi")) {
#     plot_title <- rects |> filter(experiment_names == .experiment) |> select(labels) |> pull()
#     rects <- filter(rects, experiment_names == .experiment)
#     .plotdf <- filter(.plotdf, experiment == .experiment)
#     p <- ggplot() + 
#       geom_rect(data = rects, aes(ymin = min_expr, ymax = max_expr, 
#                                   xmin = 0, xmax = 1, 
#                                   fill = experiment_names), alpha = 0.3) +
#       scale_fill_discrete(type = rects$color, labels = rects$labels, limits = rects$experiment_names) +
#       theme_classic() +
#       scale_color_discrete(type = c(.color1, .color2), labels = c(.name1, .name2)) +
#       #theme(legend.title = element_blank()) +
#       theme(legend.position = "none") +
#       ylab(ylabel) +
#       xlab("Timepoint") +
#       ylim(c(min_expr, max_expr)) +
#       theme(axis.ticks = element_blank(),
#             axis.text = element_blank()) +
#       ggtitle(plot_title)
#     
#     # lines trace average expr at each timepoint/experiment for each group
#     avgexpr1 <- .plotdf %>% filter(experiment == .experiment & group_id == 1) |> group_by(scaled_tp) |> summarise(mean_expr = mean(expr),
#                                                                                                                   sd_expr = sd(expr))
#     avgexpr2 <- .plotdf %>% filter(experiment == .experiment & group_id == 2) |> group_by(scaled_tp) |> summarise(mean_expr = mean(expr),
#                                                                                                                   sd_expr = sd(expr))
#     avgexpr3 <- .plotdf %>% filter(experiment == .experiment & group_id == 3) |> group_by(scaled_tp) |> summarise(mean_expr = mean(expr),
#                                                                                                                   sd_expr = sd(expr))
#     avgexpr4 <- .plotdf %>% filter(experiment == .experiment & group_id == 4) |> group_by(scaled_tp) |> summarise(mean_expr = mean(expr),
#                                                                                                                   sd_expr = sd(expr))
#     # adding line segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
#     p <- p + 
#       geom_ribbon(data = avgexpr1, aes(x = scaled_tp, ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
#                   fill = .color1, alpha = 0.3) +
#       geom_ribbon(data = avgexpr2, aes(x = scaled_tp, ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
#                   fill = .color2, alpha = 0.3) +
#       geom_line(data = avgexpr1, aes(x = scaled_tp, y = mean_expr), color = .color1, linewidth = 1) +
#       geom_line(data = avgexpr2, aes(x = scaled_tp, y = mean_expr), color = .color2, linewidth = 1) +
#       geom_line(data = avgexpr3, aes(x = scaled_tp, y = mean_expr), color = .color3, linewidth = 1) +
#       geom_line(data = avgexpr4, aes(x = scaled_tp, y = mean_expr), color = .color4, linewidth = 1)
#     return(p)
#   }
#   plotList <- vector(mode = "list", length = length(unique(.info1$experiment)))
#   names(plotList) <- unique(.info1$experiment)
#   for (e in ExperimentNames) {
#     plotList[[e]] <- makeSingleExperimentPlot(gdf, .experiment = e)
#   }
#   return(ggarrange(plotlist = plotList, nrow = 4, ncol = 1, common.legend = FALSE))
# }
# # tests for plotExpressionProfileQuartetStack
# # hardcoded for simplicity (subest of merge 25 yellow CCM)
# conserved_idxs <- c("YKL013C", "YER009W", "YMR097C", "YJL189W", "YKL009W",
#                     "YEL054C", "YLR333C", "YBL050W", "YNL223W", "YNL162W")
# up_par_idxs <- c("YER102W", "YLR264W", "YMR304W", "YHR193C", "YEL034W",
#                  "YOR167C", "YBL072C", "YGL135W", "YDL191W", "YHR021C")
# up_cer_idxs <- c("YMR194C-B", "YHR161C", "YJL127C-B", "YDL027C", "YNL175C",
#                  "YHR104W", "YMR027W", "YDR479C", "YFR047C", "YJL055W")
# # first yellow cer vs par, with up_cer genes indicated
# plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,conserved_idxs],
#                              .cts2 = counts_all2$par[,conserved_idxs],
#                              .cts3 = counts_all2$cer[,up_cer_idxs],
#                              .cts4 = counts_all2$par[,up_cer_idxs],
#                              .info1 = info, .info2 = info, .info3 = info, .info4 = info,
#                              .name1 = "S. cerevisiae, conserved",
#                              .name2 = "S. paradoxus, conserved",
#                              .name3 = "S. cerevisiae, up cer",
#                              .name4 = "S. paradoxus, up cer",
#                              .color1 = "orange1",
#                              .color2 = "blue2",
#                              .color3 = "orange4",
#                              .color4 = "blue4",
#                              .method = "line", .show_points = TRUE,
#                              .normalization = "log2")
# plotExpressionProfileQuartetStack(.cts1 = counts_all2$cer[,conserved_idxs],
#                                   .cts2 = counts_all2$par[,conserved_idxs],
#                                   .cts3 = counts_all2$cer[,up_cer_idxs],
#                                   .cts4 = counts_all2$par[,up_cer_idxs],
#                                   .info1 = info, .info2 = info, .info3 = info, .info4 = info,
#                                   .name1 = "S. cerevisiae, conserved",
#                                   .name2 = "S. paradoxus, conserved",
#                                   .name3 = "S. cerevisiae, up cer",
#                                   .name4 = "S. paradoxus, up cer",
#                                   .color1 = "orange1",
#                                   .color2 = "blue2",
#                                   .color3 = "orange4",
#                                   .color4 = "blue4",
#                                   .normalization = "log2")
# 
# # second yellow cer vs par, with up_par genes indicated
# plotExpressionProfileQuartetStack(.cts1 = counts_all2$cer[,conserved_idxs],
#                                   .cts2 = counts_all2$par[,conserved_idxs],
#                                   .cts3 = counts_all2$cer[,up_par_idxs],
#                                   .cts4 = counts_all2$par[,up_par_idxs],
#                                   .info1 = info, .info2 = info, .info3 = info, .info4 = info,
#                                   .name1 = "S. cerevisiae, conserved",
#                                   .name2 = "S. paradoxus, conserved",
#                                   .name3 = "S. cerevisiae, up par",
#                                   .name4 = "S. paradoxus, up par",
#                                   .color1 = "orange1",
#                                   .color2 = "blue2",
#                                   .color3 = "orange4",
#                                   .color4 = "blue4",
#                                   .normalization = "log2")


#### Versions of plotExpressionProfilePair/Quartet with continuous, unlabeled timepoint x axis ####
# plotting function to visualize expression profiles of any 2 groups
# @input: counts, info, and names of two groups to compare
# @output: ggplot of both expression profiles with loess or line curves tracing the average expression
# plotExpressionProfilePair <- function(.cts1, .cts2, .info1, .info2, 
#                                       .name1 = "S. cerevisiae", .name2 = "S. paradoxus",
#                                       .color1 = "orange1", .color2 = "blue2",
#                                       .method = c("line", "loess"), 
#                                       .show_points = FALSE,
#                                       .show_confidence_intervals = TRUE,
#                                       .normalization = c("none", "log2", "scale")) {
#   if (.normalization == "none") {
#     norm_func <- identity
#     ylabel <- "Expression"
#   }
#   if (.normalization == "log2") {
#     norm_func <- \(x) {log2(x + 1)}
#     ylabel <- "Expression (log2)"
#   }
#   if (.normalization == "scale") {
#     norm_func <- scale
#     ylabel <- "Expression (centered and scaled)"
#   }
#   ExperimentNames <- unique(.info1$experiment) # arbitrary to do info1 or info2
#   nExperiments <- length(ExperimentNames)
#   info1 <- tibble(experiment = .info1$experiment,
#                   time_point_num = .info1$time_point_num)
#   info2 <- tibble(experiment = .info2$experiment,
#                   time_point_num = .info2$time_point_num)
#   # creating sample order to organize observations on the x axis by both experiment and timepoint
#   # to make bins equidistant we need to rank samples by timepoint
#   # within each experiment then divide those ranks by the 
#   # number of timepoints in that experiment:
#   n_timepoint_table <- table(c(info1$experiment, info2$experiment),
#                              c(info1$time_point_num, info2$time_point_num))
#   info1$sample_order <- map2(info1$experiment, info1$time_point_num, \(e, tp) {
#     e_idx <- which(unique(info1$experiment) == e) - 1 # number from 0-3
#     n_timepoints <- sum(n_timepoint_table[e,] > 0)
#     local_timepoints <- info1 |> filter(experiment == e) |> 
#       select(time_point_num) |> unique() |> pull() |> sort()
#     local_rank <- which(local_timepoints == tp)
#     stopifnot(length(local_rank) == 1)
#     return(e_idx + (local_rank - 1)/(n_timepoints - 1))
#   }) |> unlist()
#   info2$sample_order <- map2(info2$experiment, info2$time_point_num, \(e, tp) {
#     e_idx <- which(unique(info2$experiment) == e) - 1 # number from 0-3
#     n_timepoints <- sum(n_timepoint_table[e,] > 0)
#     local_timepoints <- info2 |> filter(experiment == e) |> 
#       select(time_point_num) |> unique() |> pull() |> sort()
#     local_rank <- which(local_timepoints == tp)
#     stopifnot(length(local_rank) == 1)
#     return(e_idx + (local_rank - 1)/(n_timepoints - 1))
#   }) |> unlist()
#   expr1 <- norm_func(.cts1)
#   expr2 <- norm_func(.cts2)
#   gdf1 <- bind_cols(expr1, info1) |> 
#     pivot_longer(cols = colnames(expr1), names_to = "gene_name", values_to = "expr")
#   gdf1$group_id <- "1"
#   gdf2 <- bind_cols(expr2, info2) |> 
#     pivot_longer(cols = colnames(expr2), names_to = "gene_name", values_to = "expr")
#   gdf2$group_id <- "2"
#   # converting each gene's expression to its mean expression between replicates
#   gdf <- bind_rows(gdf1, gdf2) |> 
#     group_by(group_id, gene_name, experiment, time_point_num, sample_order) |> 
#     summarise(expr = mean(expr)) |> ungroup()
#   plotdf <- gdf
#   # background color rectangles for differentiating the 4 experiments
#   rects <- data.frame(color = c("orchid", "lightgreen", "gold", "lightblue"),
#                       labels = c("Urea Shock", "No Transfers", "Low Nitrogen", "Low Phosphorus"),
#                       experiment_names = c("CC", "HAP4", "LowN", "LowPi"))
#   rownames(rects) <- c("CC", "HAP4", "LowN", "LowPi")
#   rects <- rects[ExperimentNames,] # ensuring experiments are in same order as input dataset
#   rects$xstart <- c(0:(nExperiments - 1))
#   rects$xend <- c(1:nExperiments)
#   # plotting
#   p <- ggplot() + 
#     geom_rect(data = rects, aes(ymin = min(plotdf$expr), ymax = max(plotdf$expr), 
#                                 xmin = xstart, xmax = xend, 
#                                 fill = experiment_names), alpha = 0.3) +
#     theme_classic() +
#     scale_fill_discrete(type = rects$color, labels = rects$labels, limits = rects$experiment_names) +
#     scale_color_discrete(type = c(.color1, .color2), labels = c(.name1, .name2)) +
#     theme(legend.title = element_blank()) +
#     ylab(ylabel) +
#     xlab("Timepoint") +
#     theme(axis.ticks = element_blank(),
#           axis.text = element_blank())
#   if (.show_points) {
#     p <- p + geom_jitter(data = plotdf, aes(x = sample_order, y = expr, color = group_id), size = 0.1, alpha = 0.5)
#   }
#   if (.method == "loess") {
#     for (e in ExperimentNames) {
#       loess1 <- loess.sd(x = gdf %>% filter(experiment == e & group_id == 1) %>% select(sample_order) %>% pull(),
#                          y = gdf %>% filter(experiment == e & group_id == 1) %>% select(expr) %>% pull(), nsigma = 1.96)
#       loess2 <- loess.sd(x = gdf %>% filter(experiment == e & group_id == 2) %>% select(sample_order) %>% pull(),
#                          y = gdf %>% filter(experiment == e & group_id == 2) %>% select(expr) %>% pull(), nsigma = 1.96)
#       # adding loess segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
#       p <- p + 
#         geom_smooth(aes(x =!!loess1$x, y =!!loess1$y), color = .color1, linewidth = 1) +
#         geom_smooth(aes(x =!!loess2$x, y =!!loess2$y), color = .color2, linewidth = 1)
#       if (.show_confidence_intervals) {
#         p <- p + 
#           geom_ribbon(aes(x =!!loess1$x, ymin =!!loess1$lower, ymax =!!loess1$upper), fill = .color1, alpha = 0.3) +
#           geom_ribbon(aes(x =!!loess2$x, ymin =!!loess2$lower, ymax =!!loess2$upper), fill = .color2, alpha = 0.3)
#       }
#     }
#   }
#   if (.method == "line") {
#     # lines trace average expr at each timepoint/experiment for each group
#     for (e in ExperimentNames) {
#       avgexpr1 <- gdf %>% filter(experiment == e & group_id == 1) |> group_by(sample_order) |> summarise(mean_expr = mean(expr),
#                                                                                                          sd_expr = sd(expr))
#       avgexpr2 <- gdf %>% filter(experiment == e & group_id == 2) |> group_by(sample_order) |> summarise(mean_expr = mean(expr),
#                                                                                                          sd_expr = sd(expr))
#       # adding line segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
#       p <- p +
#         geom_line(data = avgexpr1, aes(x = sample_order, y = mean_expr), color = .color1, linewidth = 1) +
#         geom_line(data = avgexpr2, aes(x = sample_order, y = mean_expr), color = .color2, linewidth = 1)
#       if (.show_confidence_intervals) {
#         p <- p + 
#           geom_ribbon(data = avgexpr1, aes(x = sample_order, ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
#                       fill = .color1, alpha = 0.3) +
#           geom_ribbon(data = avgexpr2, aes(x = sample_order, ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
#                       fill = .color2, alpha = 0.3)
#       }
#     }
#   }
#   return(p)
# }
# # # tests for plotExpressionProfilePair
# # # random CCM in all but LowPi
# # random_CCM_color <- moduledf10 |> filter(is_CCM & !is.na(CCM_color)) |> select(CCM_color) |> pull() |> sample(size = 1)
# # gene_idxs <-  module_genedf10  |> filter(CCM_color == random_CCM_color) |> select(gene_name) |> pull()
# # plotExpressionProfilePair(counts_all2$cer[info$experiment != "LowPi",gene_idxs],
# #                           counts_all2$par[info$experiment != "LowPi",gene_idxs],
# #                           info[info$experiment != "LowPi",],
# #                           info[info$experiment != "LowPi",],
# #                           .name1 = "S. cereviaise",
# #                           .name2 = "S. paradoxus",
# #                           .method = "line",
# #                           .show_points = TRUE,
# #                           .show_confidence_intervals = FALSE,
# #                           .normalization = "log2") + ggtitle(random_CCM_color)
# #
# # # unclear why this would ever come up, but this makes sure that the
# # # order in which the experiments appear in the dataset doesn't
# # # affect what they're called in the plot
# # gene_idxs <-  module_genedf  |> filter(CCM_color == "cyan") |> select(gene_name) |> pull()
# # plotExpressionProfilePair(counts_all2$cer[info$experiment %in% c("CC", "HAP4"), gene_idxs],
# #                           counts_all2$par[info$experiment  %in% c("CC", "HAP4"), gene_idxs],
# #                           info[info$experiment %in% c("CC", "HAP4"),],
# #                           info[info$experiment %in% c("CC", "HAP4"),],
# #                           .name1 = "S. cereviaise",
# #                           .name2 = "S. paradoxus",
# #                           .method = "loess", .show_points = TRUE,
# #                           .normalization = "log2") + ggtitle(random_CCM_color) # cyan is steady across CC and drops dramatically in HAP4
# # # what happens when we force HAP4 to come first in the dataset?
# # plotExpressionProfilePair(rbind(counts_all2$cer[info$experiment == "HAP4",gene_idxs],
# #                                     counts_all2$cer[info$experiment == "CC",gene_idxs]),
# #                           rbind(counts_all2$par[info$experiment == "HAP4",gene_idxs],
# #                                     counts_all2$par[info$experiment == "CC",gene_idxs]),
# #                           bind_rows(info[info$experiment == "HAP4",],
# #                                     info[info$experiment == "CC",]),
# #                           bind_rows(info[info$experiment == "HAP4",],
# #                                     info[info$experiment == "CC",]),
# #                           .method = "loess", .show_points = TRUE,
# #                           .normalization = "log2") # the swoopy pattern should now come first and be called HAP4 (no transfers)
# # # # non-CCM
# # random_module <- moduledf |> filter(!is_CCM) |> select(module_name) |> pull() |> sample(1)
# # gene_idxs <- module_genedf |> filter(module_name == random_module) |> select(gene_name) |> pull()
# # length(gene_idxs)
# # plotExpressionProfilePair(counts_all2$cer[,gene_idxs, drop = FALSE],
# #                           counts_all2$par[,gene_idxs, drop = FALSE],
# #                           info,
# #                           info,
# #                           .method = "loess", .show_points = TRUE,
# #                           .normalization = "scale")
# # Oh yes
# # plots 4 groups of genes BUT you can't do this willy nilly
# # for it to be interpretable, groups 1 and 2 are the main contrast
# # and groups 3 and 4 are related to groups 1 and 2 respectively
# plotExpressionProfileQuartet <- function(.cts1, .cts2, .cts3, .cts4,
#                                          .info1, .info2, .info3, .info4,
#                                          .name1, .name2, .name3, .name4,
#                                          .color1, .color2, .color3, .color4,
#                                          .method = c("line", "loess"),
#                                          .show_points = FALSE,
#                                          .show_confidence_intervals = TRUE,
#                                          .normalization = c("none", "log2", "scale")) {
#   # first getting a plot of groups 1 and 2
#   p <- plotExpressionProfilePair(.cts1, .cts2,
#                                  .info1, .info2,
#                                  .name1, .name2,
#                                  .color1, .color2,
#                                  .method,
#                                  .show_points,
#                                  .show_confidence_intervals,
#                                  .normalization)
#   # then we'll add groups 3 and 4 to this plot
#   
#   # constructing plotdf for groups 3 and 4
#   if (.normalization == "none") {
#     norm_func <- identity
#     ylabel <- "Expression"
#   }
#   if (.normalization == "log2") {
#     norm_func <- \(x) {log2(x + 1)}
#     ylabel <- "Expression (log2)"
#   }
#   if (.normalization == "scale") {
#     norm_func <- scale
#     ylabel <- "Expression (centered and scaled)"
#   }
#   ExperimentNames <- unique(.info3$experiment) # arbitrary to do info3 or info4
#   nExperiments <- length(ExperimentNames)
#   info3 <- tibble(experiment = .info3$experiment,
#                   time_point_num = .info3$time_point_num)
#   info4 <- tibble(experiment = .info4$experiment,
#                   time_point_num = .info4$time_point_num)
#   # creating sample order to organize observations on the x axis by both experiment and timepoint
#   # to make bins equidistant we need to rank samples by timepoint
#   # within each experiment then divide those ranks by the 
#   # number of timepoints in that experiment:
#   n_timepoint_table <- table(c(info3$experiment, info4$experiment),
#                              c(info3$time_point_num, info4$time_point_num))
#   info3$sample_order <- map2(info3$experiment, info3$time_point_num, \(e, tp) {
#     e_idx <- which(unique(info3$experiment) == e) - 1 # number from 0-3
#     n_timepoints <- sum(n_timepoint_table[e,] > 0)
#     local_timepoints <- info3 |> filter(experiment == e) |> 
#       select(time_point_num) |> unique() |> pull() |> sort()
#     local_rank <- which(local_timepoints == tp)
#     stopifnot(length(local_rank) == 1)
#     return(e_idx + (local_rank - 1)/(n_timepoints - 1))
#   }) |> unlist()
#   info4$sample_order <- map2(info4$experiment, info4$time_point_num, \(e, tp) {
#     e_idx <- which(unique(info4$experiment) == e) - 1 # number from 0-3
#     n_timepoints <- sum(n_timepoint_table[e,] > 0)
#     local_timepoints <- info4 |> filter(experiment == e) |> 
#       select(time_point_num) |> unique() |> pull() |> sort()
#     local_rank <- which(local_timepoints == tp)
#     stopifnot(length(local_rank) == 1)
#     return(e_idx + (local_rank - 1)/(n_timepoints - 1))
#   }) |> unlist()
#   expr3 <- norm_func(.cts3)
#   expr4 <- norm_func(.cts4)
#   gdf3 <- bind_cols(expr3, info3) |> 
#     pivot_longer(cols = colnames(expr3), names_to = "gene_name", values_to = "expr")
#   gdf3$group_id <- "3"
#   gdf4 <- bind_cols(expr4, info4) |> 
#     pivot_longer(cols = colnames(expr4), names_to = "gene_name", values_to = "expr")
#   gdf4$group_id <- "4"
#   # converting each gene's expression to its mean expression between replicates
#   gdf <- bind_rows(gdf3, gdf4) |> 
#     group_by(group_id, gene_name, experiment, time_point_num, sample_order) |> 
#     summarise(expr = mean(expr)) |> ungroup()
#   plotdf <- gdf
#   # Adding groups 3 and 4 to the plot
#   p <- p + 
#     # this color scale will overwrite the one present in the pair plot:
#     scale_color_discrete(type = c(.color1, .color2, .color3, .color4), 
#                          labels = c(.name1, .name2, .name3, .name4),
#                          limits = c("1", "2", "3", "4"))
#   if (.show_points) {
#     p <- p + geom_jitter(data = plotdf, aes(x = sample_order, y = expr, color = group_id), size = 0.1, alpha = 0.5)
#   }
#   if (.method == "loess") {
#     for (e in ExperimentNames) {
#       loess3 <- loess.sd(x = gdf %>% filter(experiment == e & group_id == 3) %>% select(sample_order) %>% pull(),
#                          y = gdf %>% filter(experiment == e & group_id == 3) %>% select(expr) %>% pull(), nsigma = 1.96)
#       loess4 <- loess.sd(x = gdf %>% filter(experiment == e & group_id == 4) %>% select(sample_order) %>% pull(),
#                          y = gdf %>% filter(experiment == e & group_id == 4) %>% select(expr) %>% pull(), nsigma = 1.96)
#       # adding loess segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
#       p <- p + 
#         # geom_ribbon(aes(x =!!loess3$x, ymin =!!loess3$lower, ymax =!!loess3$upper), fill = .color3, alpha = 0.3) +
#         # geom_ribbon(aes(x =!!loess4$x, ymin =!!loess4$lower, ymax =!!loess4$upper), fill = .color4, alpha = 0.3) +
#         geom_smooth(aes(x =!!loess3$x, y =!!loess3$y), color = .color3, linewidth = 1) +
#         geom_smooth(aes(x =!!loess4$x, y =!!loess4$y), color = .color4, linewidth = 1)
#     }
#   }
#   if (.method == "line") {
#     # lines trace average expr at each timepoint/experiment for each group
#     for (e in ExperimentNames) {
#       avgexpr3 <- gdf %>% filter(experiment == e & group_id == 3) |> group_by(sample_order) |> summarise(mean_expr = mean(expr),
#                                                                                                          sd_expr = sd(expr))
#       avgexpr4 <- gdf %>% filter(experiment == e & group_id == 4) |> group_by(sample_order) |> summarise(mean_expr = mean(expr),
#                                                                                                          sd_expr = sd(expr))
#       # adding line segments (the =!! is from rlang and forces the mapping to not be lazily evaluated at the time of plotting):
#       p <- p + 
#         #geom_ribbon(data = avgexpr3, aes(x = sample_order, ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
#         #fill = .color3, alpha = 0.3) +
#         #geom_ribbon(data = avgexpr4, aes(x = sample_order, ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
#         #fill = .color4, alpha = 0.3) +
#         geom_line(data = avgexpr3, aes(x = sample_order, y = mean_expr), color = .color3, linewidth = 1) +
#         geom_line(data = avgexpr4, aes(x = sample_order, y = mean_expr), color = .color4, linewidth = 1)
#     }
#   }
#   return(p)
# }
# # tests for plotExpressionProfileQuartet
# # hardcoded for simplicity (subset of the unsigned module b (yellow) with positive and negatively correlated genes)
# # conserved_idxs <- c("YKL013C", "YER009W", "YMR097C", "YJL189W", "YKL009W",
# #                     "YEL054C", "YLR333C", "YBL050W", "YNL223W", "YNL162W")
# # up_par_idxs <- c("YER102W", "YLR264W", "YMR304W", "YHR193C", "YEL034W",
# #                  "YOR167C", "YBL072C", "YGL135W", "YDL191W", "YHR021C")
# # up_cer_idxs <- c("YMR194C-B", "YHR161C", "YJL127C-B", "YDL027C", "YNL175C",
# #                  "YHR104W", "YMR027W", "YDR479C", "YFR047C", "YJL055W")
# # # first yellow cer vs par, with up_cer genes indicated
# # plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,conserved_idxs],
# #                              .cts2 = counts_all2$par[,conserved_idxs],
# #                              .cts3 = counts_all2$cer[,up_cer_idxs],
# #                              .cts4 = counts_all2$par[,up_cer_idxs],
# #                              .info1 = info, .info2 = info, .info3 = info, .info4 = info,
# #                              .name1 = "S. cerevisiae, conserved",
# #                              .name2 = "S. paradoxus, conserved",
# #                              .name3 = "S. cerevisiae, up cer",
# #                              .name4 = "S. paradoxus, up cer",
# #                              .color1 = "orange1",
# #                              .color2 = "blue2",
# #                              .color3 = "orange4",
# #                              .color4 = "blue4",
# #                              .method = "line", 
# #                              .show_points = FALSE,
# #                              .show_confidence_intervals = FALSE,
# #                              .normalization = "log2")
# # # second yellow cer vs par, with up_par genes indicated
# # plotExpressionProfileQuartet(.cts1 = counts_all2$cer[,conserved_idxs],
# #                              .cts2 = counts_all2$par[,conserved_idxs],
# #                              .cts3 = counts_all2$cer[,up_par_idxs],
# #                              .cts4 = counts_all2$par[,up_par_idxs],
# #                              .info1 = info, .info2 = info, .info3 = info, .info4 = info,
# #                              .name1 = "S. cerevisiae, conserved",
# #                              .name2 = "S. paradoxus, conserved",
# #                              .name3 = "S. cerevisiae, up par",
# #                              .name4 = "S. paradoxus, up par",
# #                              .color1 = "orange1",
# #                              .color2 = "blue2",
# #                              .color3 = "orange4",
# #                              .color4 = "blue4",
# #                              .method = "line", .show_points = TRUE,
# #                              .normalization = "log2")
