#' @description Compare the enrichment of H3K9me3 and H3K27ac in transposable elements (TEs) across different tissues and ages.
#' @Author: Zhilin Long
#' @date 2026-10-27

# ~/longzhilin/Projects/AgingMouse/PairedTag/3.histoneLevel/H3K9me3/peak/H3K27ac_TE/TE8_all_valid_summary_20261027/category_K9-100_subfamily_all_valid

library(tidyverse)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(data.table)
set.seed(101)
source("~/longzhilin/ProjectCode/AgingMouse/Functions/colorAssign.R")
PT.info <- readRDS("~/longzhilin/Projects/AgingMouse/PairedTag/1.dataPreprocessing/PT.info.rds")

#### ---- category TE enrichment bubble plot ---- ####
check <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
## 100 cell types with H3K9me3 heatmap
setwd("~/longzhilin/Projects/AgingMouse/PairedTag/3.histoneLevel/H3K9me3/peak/H3K27ac_TE/TE8_all_valid_summary_20261027/category_K9-100_subfamily_all_valid")
out <- "~/longzhilin/Projects/AgingMouse/PairedTag/3.histoneLevel/H3K9me3/peak/H3K27ac_TE/TE8_all_valid_summary_20261027/category_K9-100_subfamily_all_valid"
category_order <- c("Glutamatergic", "GABAergic", "Glial", "Endothelial", "Stromal", "Epithelial", "Myeloid", "Lymphoid", "Adipose", "Muscle", "Germline")
family_path <- file.path(out, "tables/category_TE_family_enrichment.tsv")
subfamily_path <- file.path(out, "tables/category_TE_subfamily_enrichment.tsv.gz")

family <- fread(family_path)
subfamily <- fread(subfamily_path)
family_order <- c("LTR:ERVK", "LTR:ERV1", "LINE:L1", "SINE:Alu", "SINE:B2", "SINE:B4",
                  "LTR:ERVL-MaLR", "LTR:ERVL") # > 1% genome
te_classes <- c("LINE", "LINE?", "LTR", "LTR?", "SINE", "SINE?", "DNA", "DNA?", "RC", "RC?")
ranking <- subfamily[class %in% te_classes & category %in% category_order,
                     .(n_categories_fdr_lt_005 = sum(fdr_two_sided < .05),
                       n_significant_enriched = sum(fdr_two_sided < .05 & fold_enrichment > 1, na.rm = TRUE),
                       n_significant_depleted = sum(fdr_two_sided < .05 & fold_enrichment < 1, na.rm = TRUE),
                       total_observed = sum(observed), total_expected = sum(expected),
                       total_expected_exact = sum(expected_exact)),
                     by = .(feature, subfamily, family, class)]
setorderv(ranking, c("n_categories_fdr_lt_005", "total_observed", "feature"), c(-1L, -1L, 1L))
ranking[, rank := .I]
top <- ranking[1:30]
top[, display_label := paste(family, subfamily, sep = ":")]
data <- subfamily[feature %in% top$feature & category %in% category_order]
data[, `:=`(rank = match(feature, top$feature),
            display_label = paste(family, subfamily, sep = ":"),
            n_categories_fdr_lt_005 = top$n_categories_fdr_lt_005[match(feature, top$feature)])]
setorder(data, rank, category)

make_scale <- function(d, metric) {
  d <- copy(d)
  field <- if (metric == "pvalue") "p_two_sided" else "fdr_two_sided"
  prob <- d[[field]]
  check(all(is.finite(prob) & prob > 0 & prob <= 1), "Invalid P or FDR")
  cutoff <- -log10(.05)
  raw_color <- d$log2_fold_enrichment
  finite_color <- raw_color[is.finite(raw_color)]
  check(length(finite_color) > 0L, "No finite fold enrichments")
  # A symmetric integer color bound spans every finite observed log2 fold change.
  color_bound <- max(1, ceiling(max(abs(finite_color))))
  color_value <- pmax(-color_bound, pmin(color_bound, raw_color))
  size_value <- pmax(cutoff, -log10(prob))
  actual_size_max <- max(size_value)
  size_limit_max <- if (actual_size_max > cutoff) actual_size_max else cutoff + 1e-6
  half_step_breaks <- if (actual_size_max >= 1.5) seq(1.5, floor(actual_size_max * 2) / 2, by = .5) else numeric()
  # If no half-step tick fits, add one data-supported decimal tick when possible.
  decimal_break <- if (!length(half_step_breaks) && actual_size_max > cutoff) floor(actual_size_max * 10) / 10 else numeric()
  if (length(decimal_break) && decimal_break <= cutoff) decimal_break <- numeric()
  size_breaks <- c(cutoff, half_step_breaks, decimal_break)
  size_labels <- c("≤ 1.30", trimws(formatC(c(half_step_breaks, decimal_break), format = "fg", digits = 3)))
  check(all(size_breaks >= cutoff & size_breaks <= size_limit_max), "Size tick outside observed range")
  set(d, j = "color_value", value = color_value)
  set(d, j = "size_value", value = size_value)
  check(isTRUE(all.equal(d$size_value, size_value)) &&
          isTRUE(all.equal(d$color_value, color_value)), "Plot value assignment differs")
  list(data = d, color_bound = color_bound, size_max = actual_size_max,
       size_limit_max = size_limit_max, size_breaks = size_breaks, size_labels = size_labels,
       finite_color_min = min(finite_color), finite_color_max = max(finite_color),
       n_negative_infinite = sum(is.infinite(raw_color) & raw_color < 0, na.rm = TRUE),
       n_positive_infinite = sum(is.infinite(raw_color) & raw_color > 0, na.rm = TRUE),
       n_undefined = sum(is.na(raw_color)), n_minimum_size = sum(prob >= .05))
}

make_plot <- function(s, features, display, level, metric, subtitle) {
  d <- copy(s$data)
  d[, `:=`(feature = factor(feature, levels = features),
           category = factor(category, levels = rev(category_order)))]
  check(!anyNA(d$feature) && !anyNA(d$category) && !anyDuplicated(display),
        "Unexpected figure labels")
  color_labels <- c(if (s$n_negative_infinite > 0) paste0("≤ −", s$color_bound) else paste0("−", s$color_bound),
                    "0", if (s$n_positive_infinite > 0) paste0("≥ ", s$color_bound) else as.character(s$color_bound))
  ggplot(d, aes(feature, category)) +
    geom_point(data = d[!is.na(color_value)], aes(color = color_value, size = size_value)) +
    geom_point(data = d[is.na(color_value)], shape = 4, color = "#999999", size = 1.3) +
    scale_color_gradient2(low = "#3F78AC", mid = "#F1F1F1", high = "#B92635", midpoint = 0,
                          limits = c(-s$color_bound, s$color_bound),
                          breaks = c(-s$color_bound, 0, s$color_bound), labels = color_labels,
                          name = "log2(observed / expected)") +
    scale_size_continuous(range = c(1, 3), limits = c(-log10(.05), s$size_limit_max),
                          breaks = s$size_breaks, labels = s$size_labels,
                          name = if (metric == "pvalue") "−log10(P value)" else "−log10(FDR)") +
    scale_x_discrete(labels = setNames(display, features), drop = FALSE, expand = expansion(add = .55)) +
    labs(title = paste("TE", level, "enrichment across cell categories"),
         subtitle = subtitle, x = if (level == "family") "TE family" else "TE family:subfamily", y = NULL,
         caption = paste0("10,000 chromosome-stratified draws. ",
                          "Scale based on displayed data; zero observed maps to the negative color endpoint.\n",
                          "Size uses two-sided ", if (metric == "pvalue") "P" else "global FDR",
                          "; values ≥ 0.05 use minimum size. Cross: undefined ratio.")) +
    guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) +
    theme_minimal(base_size = 8, base_family = "Arial") +
    theme(panel.grid.major = element_line(color = "#EEEEEE", linewidth = .2),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = if (length(features) >= 14L) 90 else 55,
                                     hjust = 1, vjust = if (length(features) >= 14L) .5 else 1, size = 7),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 8),
          plot.caption = element_text(size = 6, hjust = 0),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.title = element_text(size = 7), legend.text = element_text(size = 7))
}

for (level in c("family", "subfamily")) {
  if(level == "family") {
    width = 90/25.4
  } else {
    width = 150/25.4
  }

  for (metric in c("pvalue", "fdr")) {
    input <- if (level == "family") family else subfamily
    features <- if (level == "family") family_order else top$feature
    display <- if (level == "family") sub("^[^:]+:", "", features) else top$display_label
    subtitle <- if (level == "family") "Eight TE families each covering >1% of the mm10 genome" else
      "Top 30 by number of categories with two-sided FDR < 0.05"
    if (anyDuplicated(display)) {
      display <- features
    }
    # Filter before calculating scales.
    input_plot <- input[feature %in% features & category %in% category_order]
    scale <- make_scale(input_plot, metric)
    plot <- make_plot(scale, features, display, level, metric, subtitle)
    grDevices::cairo_pdf(paste0("figures/final-", level, "-", metric, "-bubble.pdf"), width = width, height = 100/25.4, family = "Arial")
    print(plot)
    grDevices::dev.off()
  }
}
