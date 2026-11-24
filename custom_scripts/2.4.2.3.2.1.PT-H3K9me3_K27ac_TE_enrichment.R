#!/usr/bin/env Rscript

# Reproduce selected-100 category TE family and subfamily enrichment.
# 1. Read raw differential peaks, C1 domains, RepeatMasker and cell metadata.
# 2. Rebuild foreground peaks and the all-valid peak sampling universe.
# 3. Compile the embedded C++ kernel and perform 10,000 chromosome-stratified
#    randomizations for each selected cell type using its original global seed.
# 4. Sum the observed and matching iteration counts within each category.
# 5. Compute family and subfamily empirical tests and write result tables.
#

suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
options(stringsAsFactors = FALSE)

check <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

progress <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), ..., "\n")
  flush.console()
}

load_kernel <- function() {
  Rcpp::sourceCpp(
    code = paste(
      "#include <Rcpp.h>",
      "#include <algorithm>",
      "#include <numeric>",
      "#include <random>",
      "#include <vector>",
      "",
      "// [[Rcpp::plugins(cpp11)]]",
      "",
      "// Sample peak IDs without replacement and count each TE locus once per draw.",
      "// [[Rcpp::export]]",
      "Rcpp::IntegerMatrix simulate_te_counts_cpp(",
      "    Rcpp::List strata, Rcpp::IntegerVector take, Rcpp::List peak_loci,",
      "    Rcpp::IntegerVector locus_family, int n_families, int iterations, int seed) {",
      "  const int n_peaks = peak_loci.size();",
      "  std::vector<std::vector<int>> groups(strata.size());",
      "  std::vector<std::vector<int>> loci(n_peaks);",
      "  for (int g = 0; g < strata.size(); ++g) {",
      "    Rcpp::IntegerVector ids = strata[g];",
      "    groups[g].assign(ids.begin(), ids.end());",
      "  }",
      "  for (int p = 0; p < n_peaks; ++p) {",
      "    Rcpp::IntegerVector ids = peak_loci[p];",
      "    loci[p].assign(ids.begin(), ids.end());",
      "  }",
      "",
      "  Rcpp::IntegerMatrix output(n_families, iterations);",
      "  std::vector<int> seen(locus_family.size(), 0);",
      "  std::mt19937 generator(static_cast<unsigned int>(seed));",
      "  for (int b = 0; b < iterations; ++b) {",
      "    if (b % 1000 == 0) Rcpp::checkUserInterrupt();",
      "    const int stamp = b + 1;",
      "    for (int g = 0; g < static_cast<int>(groups.size()); ++g) {",
      "      auto& pool = groups[g];",
      "      for (int k = 0; k < take[g]; ++k) {",
      "        std::uniform_int_distribution<int> pick(k, pool.size() - 1);",
      "        std::swap(pool[k], pool[pick(generator)]);",
      "        for (const int locus_id : loci[pool[k] - 1]) {",
      "          const int t = locus_id - 1;",
      "          if (seen[t] != stamp) {",
      "            seen[t] = stamp;",
      "            ++output(locus_family[t] - 1, b);",
      "          }",
      "        }",
      "      }",
      "    }",
      "  }",
      "  return output;",
      "}",
      sep = "\n"
    ),
    showOutput = FALSE
  )
}

parse_regions <- function(regions) {
  check(!anyNA(regions) && !anyDuplicated(regions), "Missing or duplicate peak Regions")
  check(all(grepl("^[^:]+:[^-]+-[^-]+$", regions)), "Malformed peak Region")
  start <- as.numeric(sub("-.*$", "", sub("^[^:]+:", "", regions)))
  end <- as.numeric(sub("^.*-", "", regions))
  check(all(is.finite(start) & is.finite(end) & start >= 1 & end >= start &
              start == floor(start) & end == floor(end)),
        "Invalid 1-based inclusive peak coordinates")
  chr <- sub(":.*$", "", regions)
  data.table(Chr = chr, Start = start, End = end, Region = regions,
              canonical = paste0(chr, ":", sprintf("%.0f", start), "-", sprintf("%.0f", end)))
}

select_c1_candidates <- function(d, c1, valid_regions = NULL) {
  check(all(c("Region", "PValue", "logFC") %in% names(d)), "Missing differential columns")
  peaks <- parse_regions(as.character(d$Region))
  if (!is.null(valid_regions)) {
    check(all(peaks$canonical %in% valid_regions), "Peak absent from converted upstream BED")
  }
  valid <- is.finite(d$PValue) & d$PValue >= 0 & d$PValue <= 1 & is.finite(d$logFC)
  setkey(c1, Chr, Start, End)
  hits <- foverlaps(peaks[valid], c1, type = "any", nomatch = 0L)
  keep <- unique(hits[pmin(End, i.End) - pmax(Start, i.Start) + 1 >= 250, Region])
  peaks[, `:=`(PValue = d$PValue, logFC = d$logFC,
                logCPM = if ("logCPM" %in% names(d)) d$logCPM else NA_real_)]
  peaks[, Up := valid & PValue <= 0.05 & logFC > 0]
  peaks[Region %in% keep]
}

center_te_pairs <- function(peaks, te) {
  query <- peaks[, .(Chr, Start = floor((Start + End)/2),
                     End = floor((Start + End)/2), Region)]
  setkey(te, Chr, Start, End)
  hits <- foverlaps(query, te, type = "any", nomatch = 0L)
  check(all(hits$Start <= hits$i.Start & hits$i.Start <= hits$End), "Non-center TE hit")
  unique(hits[, .(Region, TE, subfamily)])
}

count_unique_te <- function(pairs, regions, subfamilies) {
  loci <- unique(pairs[Region %in% regions, .(TE, subfamily)])
  check(!anyDuplicated(loci$TE), "Ambiguous TE-to-subfamily mapping")
  family_ids <- match(loci$subfamily, subfamilies)
  check(!anyNA(family_ids), "Unknown TE subfamily")
  tabulate(family_ids, nbins = length(subfamilies))
}

build_simulation_input <- function(candidates, pairs) {
  check(!anyDuplicated(candidates$Region), "Duplicate sampling candidate")
  hits <- copy(pairs[Region %in% candidates$Region])
  loci <- unique(hits[, .(TE, subfamily)])
  check(!anyDuplicated(loci$TE), "Conflicting TE locus annotations")
  families <- sort(unique(loci$subfamily), method = "radix")
  peak_loci <- replicate(nrow(candidates), integer(), simplify = FALSE)
  if (nrow(hits)) {
    hits[, `:=`(peak_id = match(Region, candidates$Region), locus_id = match(TE, loci$TE))]
    mapping <- hits[, .(ids = list(unique(locus_id))), by = peak_id]
    peak_loci[mapping$peak_id] <- mapping$ids
  }
  strata <- split(seq_len(nrow(candidates)), candidates$Chr, drop = TRUE)
  take <- as.integer(vapply(strata, function(ids) sum(candidates$Up[ids]), integer(1)))
  list(strata = strata, take = take, peak_loci = peak_loci,
       locus_family = match(loci$subfamily, families), families = families,
       loci = loci)
}

simulate_null <- function(input, iterations, seed) {
  check(length(iterations) == 1L && iterations >= 1 && iterations == as.integer(iterations),
        "Invalid iteration count")
  check(length(input$strata) == length(input$take) &&
          all(input$take >= 0 & input$take <= lengths(input$strata)), "Invalid stratum sample sizes")
  peak_ids <- unlist(input$strata, use.names = FALSE)
  check(!anyDuplicated(peak_ids) && all(peak_ids >= 1 & peak_ids <= length(input$peak_loci)),
        "Invalid or repeated peak IDs across strata")
  ids <- unlist(input$peak_loci, use.names = FALSE)
  check(all(ids >= 1 & ids <= length(input$locus_family)), "Invalid TE locus IDs")
  check(all(input$locus_family >= 1 & input$locus_family <= length(input$families)),
        "Invalid subfamily IDs")
  simulate_te_counts_cpp(input$strata, input$take, input$peak_loci,
                         as.integer(input$locus_family), length(input$families),
                         as.integer(iterations), as.integer(seed))
}

exact_expectation <- function(input) {
  # A TE is missed only when all selected peaks avoid every candidate covering it.
  log_miss <- rep(0, length(input$locus_family))
  for (g in seq_along(input$strata)) {
    population <- length(input$strata[[g]])
    draws <- input$take[g]
    if (!draws) next
    k <- tabulate(unlist(input$peak_loci[input$strata[[g]]], use.names = FALSE),
                  nbins = length(input$locus_family))
    possible <- k <= population - draws
    contribution <- rep(-Inf, length(k))
    contribution[possible] <- lchoose(population - k[possible], draws) - lchoose(population, draws)
    log_miss <- log_miss + contribution
  }
  probability <- -expm1(log_miss)
  vapply(seq_along(input$families), function(s) sum(probability[input$locus_family == s]), numeric(1))
}

summarize_null <- function(observed, null) {
  check(nrow(null) == length(observed) && ncol(null) >= 2L, "Invalid null count matrix")
  iterations <- ncol(null)
  expected <- rowMeans(null)
  exceed <- vapply(seq_len(nrow(null)), function(s) sum(null[s, ] >= observed[s]), integer(1))
  deviation <- vapply(seq_len(nrow(null)), function(s) sd(null[s, ]), numeric(1))
  quantiles <- vapply(seq_len(nrow(null)), function(s) {
    quantile(null[s, ], c(.025, .975), type = 1, names = FALSE)
  }, numeric(2))
  fold <- observed / expected
  fold[observed == 0 & expected == 0] <- NA_real_
  data.table(expected = expected, null_sd = deviation,
             expected_mc_se = deviation/sqrt(iterations),
             null_q025 = quantiles[1L, ], null_q975 = quantiles[2L, ],
             n_exceed = exceed, p_enrichment = (exceed + 1)/(iterations + 1),
             fold_enrichment = fold, log2_fold_enrichment = log2(fold))
}

save_checkpoint <- function(object, path) {
  temporary <- paste0(path, ".incomplete")
  saveRDS(object, temporary, compress = "gzip")
  check(file.rename(temporary, path), paste("Cannot publish checkpoint", path))
}

validate_preparation <- function(prepared, input_hashes, source_hashes) {
  check(identical(prepared$input_hashes, input_hashes), "Preparation input hashes changed")
  check(identical(prepared$source_hashes, source_hashes), "Preparation source hashes changed")
  invisible(TRUE)
}

input_paths <- function() {
  root <- "~/suzhuojie/projects/pairedtag_27ac_4me1"
  list(
    differential = "~/longzhilin/Projects/AgingMouse/PairedTag/3.histoneLevel/H3K27ac/peak/mergePeak_zhuojie/differentialPeakAnalysis/diffPeak.list.rds",
    c1 = "~/longzhilin/public_html/AgingMouse/PairedTag_new/MergeIGV/H3K9me3/H3K9me3_tissue_peak_C1.bed",
    te = "~/suzhuojie/ref_data/TE_reference/mm10_TE_metadata_20250116.bed",
    gtf = "~/suzhuojie/ref_data/TE_reference/mm10_rmsk_TE.gtf",
    merged = file.path(root, "result/call_merge_peak_20260812/H3K27ac/peak_merge/H3K27ac_spm5_all_chr_1kb/all_celltype_summit_all_chr_merge_1kb.bed"),
    observed = file.path(root, "result/H3K27ac/C1_up_peak_TE_subfamily_counts/peak_center_TE_number/tables/counts/TE_subfamily_by_celltype_counts.tsv")
  )
}

prepare_inputs <- function(paths, run_dir, hashes, source_hashes) {
  reference <- fread(paths$observed, check.names = FALSE)
  labels <- names(reference)[-1L]
  check(nrow(reference) == 1146L && length(labels) == 195L, "Unexpected reference matrix dimensions")
  c1 <- fread(paths$c1, header = FALSE, select = 1:3, col.names = c("Chr", "Start", "End"))
  check(nrow(c1) == 172L && !anyNA(c1) && !anyDuplicated(c1), "Unexpected C1 intervals")
  check(all(c1$Start >= 0 & c1$End >= c1$Start), "Invalid C1 coordinates")
  merged <- fread(paths$merged, header = FALSE, select = 1:3)
  valid_regions <- paste0(merged[[1L]], ":", sprintf("%.0f", merged[[2L]] + 1),
                          "-", sprintf("%.0f", merged[[3L]]))
  check(!anyDuplicated(valid_regions), "Duplicate merged peak coordinates")
  rm(merged)

  progress("Reading differential peak list")
  differential <- readRDS(paths$differential)
  check(setequal(labels, names(differential)), "Reference and differential cell types differ")
  candidates <- setNames(vector("list", length(labels)), labels)
  summaries <- vector("list", length(labels))
  constraints <- vector("list", length(labels))
  for (i in seq_along(labels)) {
    d <- as.data.table(differential[[labels[i]]])
    candidates[[i]] <- select_c1_candidates(d, c1, valid_regions)
    u <- candidates[[i]]
    summaries[[i]] <- data.table(celltype = labels[i], n_tested = nrow(d),
                                 n_valid_tested = sum(is.finite(d$PValue) & d$PValue >= 0 &
                                                       d$PValue <= 1 & is.finite(d$logFC)),
                                 n_candidates = nrow(u), n_foreground = sum(u$Up),
                                 n_other = sum(!u$Up))
    constraints[[i]] <- u[, .(n_candidates = .N, n_foreground = sum(Up)), by = Chr][
      , `:=`(celltype = labels[i], fixed_stratum = n_candidates == n_foreground & n_foreground > 0)]
    if (i %% 25L == 0L || i == length(labels)) progress("C1 candidates", i, "/", length(labels))
  }
  summary <- rbindlist(summaries)
  check(sum(summary$n_foreground) == 8803L, "Original foreground peak count changed")
  rm(differential, d, valid_regions)
  invisible(gc())

  progress("Reading and checking TE reference BED/GTF")
  te <- fread(paths$te, header = FALSE, select = 1:4,
               col.names = c("Chr", "Start", "End", "TE"))
  gtf <- fread(paths$gtf, header = FALSE, sep = "\t", quote = "", select = c(1, 4, 5, 9),
                col.names = c("Chr", "Start", "End", "Attributes"))
  get_attr <- function(key) {
    value <- sub(paste0('.*', key, ' "([^"]+)".*'), '\\1', gtf$Attributes)
    check(all(nzchar(value) & value != gtf$Attributes), paste("Missing GTF attribute", key))
    value
  }
  ids <- get_attr("transcript_id")
  gene <- get_attr("gene_id")
  check(nrow(te) == 3725827L && nrow(gtf) == nrow(te), "Unexpected TE reference row count")
  check(!anyDuplicated(ids) && !anyDuplicated(te$TE), "Duplicate TE locus IDs")
  matching <- match(te$TE, ids)
  check(!anyNA(matching), "Unannotated TE locus")
  check(all(te$Chr == gtf$Chr[matching] & te$Start == gtf$Start[matching] &
              te$End == gtf$End[matching]), "TE BED/GTF coordinate mismatch")
  check(!anyNA(te) && all(te$Start >= 1 & te$End >= te$Start), "Invalid TE coordinates")
  te[, subfamily := gene[matching]]
  annotation <- unique(data.table(subfamily = gene, family = get_attr("family_id"), class = get_attr("class_id")))
  check(!anyDuplicated(annotation$subfamily), "Ambiguous subfamily classification")
  annotation[, annotated_subfamily := paste(class, family, subfamily, sep = ":")]
  annotation <- annotation[match(reference[[1L]], annotated_subfamily)]
  check(nrow(annotation) == 1146L && !anyNA(annotation), "Reference subfamily labels differ")
  subfamilies <- annotation$subfamily
  rm(gtf, ids, gene, matching)
  invisible(gc())

  all_peaks <- unique(rbindlist(lapply(candidates, function(u) u[, .(Chr, Start, End, Region)])))
  pairs <- center_te_pairs(all_peaks, te)
  foreground_regions <- unique(unlist(lapply(candidates, function(u) u[Up == TRUE, Region]), use.names = FALSE))
  check(length(foreground_regions) == 3522L, "Original distinct foreground peak set changed")
  observed <- vapply(candidates, function(u) count_unique_te(pairs, u[Up == TRUE, Region], subfamilies),
                      integer(length(subfamilies)))
  check(identical(unname(observed), unname(as.matrix(reference[, -1L]))),
        "Reconstructed foreground TE matrix differs from supplied matrix")
  check(sum(observed) == 2984L, "Original foreground TE count changed")

  fwrite(summary, file.path(run_dir, "celltype_summary.tsv"), sep = "\t")
  fwrite(rbindlist(constraints), file.path(run_dir, "chromosome_constraints.tsv"), sep = "\t")
  fwrite(rbindlist(candidates, idcol = "celltype"), file.path(run_dir, "c1_candidate_peaks.tsv.gz"), sep = "\t")
  fwrite(pairs, file.path(run_dir, "peak_center_TE_pairs.tsv.gz"), sep = "\t")
  fwrite(annotation, file.path(run_dir, "TE_subfamily_annotation.tsv"), sep = "\t")
  fwrite(reference, file.path(run_dir, "observed_TE_counts.tsv"), sep = "\t")
  progress("FOREGROUND VERIFIED: all 1146 x 195 counts match; 8803 peak calls; 3522 distinct peaks; 2984 TE locus-celltype calls")
  writeLines("All 223470 observed entries match the supplied matrix exactly.", file.path(run_dir, "FOREGROUND_VERIFIED"))
  list(candidates = candidates, pairs = pairs, annotation = annotation,
       observed = observed, labels = labels, summary = summary,
       input_hashes = hashes, source_hashes = source_hashes)
}

run_celltype <- function(i, prepared, config, run_dir) {
  label <- prepared$labels[i]
  path <- file.path(run_dir, "celltype_null", sprintf("%03d.rds", i))
  cell_seed <- as.integer(config$seed + i * 7919L)
  if (file.exists(path)) {
    saved <- readRDS(path)
    check(identical(saved$config, config) && identical(saved$celltype, label), "Incompatible checkpoint")
    return(saved$result)
  }
  universe <- prepared$candidates[[i]]
  observed <- prepared$observed[, i]
  n_foreground <- sum(universe$Up)
  result <- data.table(celltype = label, prepared$annotation, observed = observed,
                       n_candidates = nrow(universe), n_foreground = n_foreground,
                       n_iterations = if (n_foreground) config$iterations else 0L,
                       seed = cell_seed)
  numeric_columns <- c("expected", "expected_exact", "null_sd", "expected_mc_se", "null_q025", "null_q975",
                       "n_exceed", "p_enrichment", "fold_enrichment", "log2_fold_enrichment")
  for (column in numeric_columns) set(result, j = column, value = NA_real_)
  result[, status := "no_foreground"]
  null <- NULL
  supported_families <- character()
  if (n_foreground) {
    # Structurally unsupported families remain in the predefined testing family.
    result[, `:=`(expected = 0, expected_exact = 0, null_sd = 0, expected_mc_se = 0,
                   null_q025 = 0, null_q975 = 0, n_exceed = config$iterations,
                   p_enrichment = 1, status = "zero_null_support")]
    input <- build_simulation_input(universe, prepared$pairs)
    supported_families <- input$families
    if (length(input$families)) {
      family_index <- match(input$families, prepared$annotation$subfamily)
      forced <- input
      forced$strata <- lapply(input$strata, function(ids) ids[universe$Up[ids]])
      forced$take <- as.integer(lengths(forced$strata))
      check(identical(as.integer(simulate_null(forced, 1L, cell_seed)),
                      as.integer(observed[family_index])), paste("Real foreground kernel check failed", label))
      null <- simulate_null(input, config$iterations, cell_seed)
      stats <- summarize_null(observed[family_index], null)
      exact <- exact_expectation(input)
      # Rare hits can have zero sampled variance; retain their analytic Monte Carlo scale.
      check(all(abs(stats$expected - exact) <= 8 * pmax(stats$expected_mc_se,
                                                       sqrt(exact/config$iterations), 1/config$iterations) +
                  1/config$iterations), paste("Monte Carlo expectation calibration failed", label))
      for (column in names(stats)) set(result, i = family_index, j = column, value = as.numeric(stats[[column]]))
      set(result, i = family_index, j = "expected_exact", value = exact)
      result[expected_exact > 0, status := "tested"]
      result[observed > 0 & expected == 0, status := "null_zero_in_simulation"]
      result[expected_exact > 0 & null_sd == 0 & expected > 0, status := "fixed_in_simulation"]
      check(all(result$observed[result$expected_exact == 0] == 0), "Foreground outside null support")
      rownames(null) <- prepared$annotation$annotated_subfamily[family_index]
    }
  }
  result[, at_p_resolution := !is.na(p_enrichment) & n_exceed == 0]
  save_checkpoint(list(config = config, celltype = label, supported_subfamilies = supported_families,
                        result = result, null_counts = null), path)
  progress("Celltype finished", i, "/", length(prepared$labels), label,
             "foreground", n_foreground, "supported subfamilies", length(supported_families))
  result
}


# 1. Set input and output paths; preserve the original global cell-type order.
args <- commandArgs(trailingOnly = TRUE)
check(length(args) %in% 1:2,
      "Usage: Rscript SCRIPT.R NEW_OUTPUT_DIR [WORKERS]")
out <- args[1L]
workers <- if (length(args) == 2L) as.integer(args[2L]) else 4L
check(!is.na(workers) && workers >= 1L, "Invalid worker count")
check(!dir.exists(out) || length(list.files(out, all.files = TRUE,
                                           no.. = TRUE)) == 0L,
      "Output directory must be new or empty")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)
run <- file.path(out, "randomization")
dir.create(run)
root <- "~/longzhilin/Projects/AgingMouse/PairedTag"
heatmap_path <- file.path(
  root, "3.histoneLevel/H3K9me3/peak/heatmap.order.tissue.final.rds"
)
info_path <- file.path(root, "1.dataPreprocessing/PT.info.rds")
tables <- file.path(out, "tables")
null_dir <- file.path(out, "category_null")
dir.create(tables)
dir.create(null_dir)
paths <- input_paths()
all_inputs <- c(unlist(paths), heatmap = heatmap_path, info = info_path)
hashes <- tools::md5sum(all_inputs)
check(!anyNA(hashes), "A required raw input or metadata file is missing")
self_arg <- sub("^--file=", "", commandArgs(FALSE)[
  grepl("^--file=", commandArgs(FALSE))
])
check(length(self_arg) == 1L, "Cannot locate the running script")
source_hashes <- tools::md5sum(normalizePath(self_arg))
data.table::fwrite(
  data.table::data.table(path = names(hashes), md5 = unname(hashes)),
  file.path(out, "raw_input_manifest.tsv"), sep = "\t"
)
config <- list(
  iterations = 10000L, workers = workers, seed = 20261027L,
  null_model = "all_valid_tested_peaks_chromosome_stratified_without_replacement",
  input_hashes = hashes, source_hashes = source_hashes
)

# 2. Rebuild C1 foreground and validate all 195 original observed columns.
progress("Preparing foreground from raw peak and TE inputs")
old <- prepare_inputs(paths, run, hashes, source_hashes)
selected_cells <- names(readRDS(heatmap_path)$column_order)
check(length(selected_cells) == 100L && !anyNA(selected_cells) &&
        !anyDuplicated(selected_cells) && all(selected_cells %in% old$labels),
      "Invalid selected 100-cell-type list")
differential <- readRDS(paths$differential)
candidates <- lapply(seq_along(old$labels), function(i) {
  d <- as.data.table(differential[[old$labels[i]]])
  foreground <- old$candidates[[i]][Up == TRUE, Region]
  check(all(c("Region", "PValue", "logFC") %in% names(d)),
        "Missing differential peak columns")
  valid <- is.finite(d$PValue) & d$PValue >= 0 & d$PValue <= 1 &
    is.finite(d$logFC)
  regions <- as.character(d$Region[valid])
  check(!anyNA(regions) && !anyDuplicated(regions) &&
          all(foreground %in% regions), "Invalid all-valid sampling pool")
  data.table(Region = regions, Chr = sub(":.*$", "", regions),
             Up = regions %in% foreground)
})
names(candidates) <- old$labels
rm(differential)
summary <- rbindlist(lapply(seq_along(candidates), function(i) {
  pool <- candidates[[i]]
  data.table(celltype = old$labels[i], n_candidates = nrow(pool),
             n_foreground = sum(pool$Up), n_other = sum(!pool$Up))
}))
check(nrow(summary) == 195L && sum(summary$n_foreground) == 8803L &&
        sum(summary$n_candidates) == 20394367L,
      "All-valid candidates differ from the reference run")

# Recompute peak-center TE pairs from the input BED, including non-C1 peaks.
progress("Rebuilding all-valid peak-center TE overlaps")
regions <- unique(unlist(lapply(candidates, function(x) x$Region),
                         use.names = FALSE))
peaks <- parse_regions(regions)
merged <- fread(paths$merged, header = FALSE, select = 1:3)
merged_regions <- paste0(merged[[1L]], ":", sprintf("%.0f", merged[[2L]] + 1),
                         "-", sprintf("%.0f", merged[[3L]]))
check(!anyDuplicated(merged_regions) &&
        all(peaks$canonical %in% merged_regions),
      "A valid peak is absent from the upstream merged BED")
te <- fread(paths$te, header = FALSE, select = 1:4,
            col.names = c("Chr", "Start", "End", "TE"))
gtf <- fread(paths$gtf, header = FALSE, sep = "\t", quote = "",
             select = c(1, 4, 5, 9),
             col.names = c("Chr", "Start", "End", "Attributes"))
attr_value <- function(key) {
    value <- sub(paste0('.*', key, ' "([^"]+)".*'), '\\1', gtf$Attributes)
  check(all(nzchar(value) & value != gtf$Attributes),
        paste("Missing GTF attribute", key))
  value
}
ids <- attr_value("transcript_id")
gene <- attr_value("gene_id")
check(!anyDuplicated(ids) && !anyDuplicated(te$TE),
      "Duplicate TE locus identifiers")
index <- match(te$TE, ids)
check(!anyNA(index) && all(te$Chr == gtf$Chr[index] &
        te$Start == gtf$Start[index] & te$End == gtf$End[index]),
      "TE BED/GTF annotations or coordinates do not match")
te[, subfamily := gene[index]]
pairs <- center_te_pairs(peaks, te)
observed <- vapply(candidates, function(pool) {
  count_unique_te(pairs, pool[Up == TRUE, Region], old$annotation$subfamily)
}, integer(nrow(old$annotation)))
check(identical(unname(observed), unname(old$observed)),
      "Recomputed all-valid foreground counts differ from original matrix")
rm(regions, peaks, merged, merged_regions, te, gtf, ids, gene, index)
invisible(gc())
prepared <- list(
  candidates = candidates, pairs = pairs, annotation = old$annotation,
  observed = observed, labels = old$labels, summary = summary,
  input_hashes = hashes, source_hashes = source_hashes
)
fwrite(summary, file.path(run, "celltype_summary.tsv"), sep = "\t")
saveRDS(config, file.path(run, "run_config.rds"))
dir.create(file.path(run, "celltype_null"))

# 3. Keep global indices so each selected cell receives its historical seed.
load_kernel()
indices <- match(selected_cells, prepared$labels)
progress("Starting randomized sampling for", length(indices), "cell types")
results <- parallel::mclapply(
  indices, run_celltype, prepared = prepared, config = config, run_dir = run,
  mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE
)
failed <- vapply(results, function(x) {
  inherits(x, "try-error") || !is.data.table(x)
}, logical(1))
check(!any(failed), paste("Failed selected cell types:",
                          paste(selected_cells[failed], collapse = ", ")))
check(identical(hashes, tools::md5sum(all_inputs)),
      "Raw inputs changed during sampling")

# 4. Aggregate observed and corresponding random draws by category, and
# 5. Compute and export family and subfamily results.
check <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
order_categories <- c("Glutamatergic", "GABAergic", "Glial", "Endothelial", "Stromal",
                      "Epithelial", "Myeloid", "Lymphoid", "Adipose", "Muscle",
                      "Germline", "Others")
group_families <- function(a) {
  fam <- as.character(a$family)
  dna <- grepl("^DNA[?]?$", a$class)
  fam[dna & grepl("^hAT($|-)", fam)] <- "hAT"
  fam[dna & grepl("^TcMar($|-)", fam)] <- "TcMar"
  paste(a$class, fam, sep = ":")
}
group_counts <- function(mat, groups, levels) {
  check(is.matrix(mat) && nrow(mat) == length(groups), "Invalid matrix grouping")
  ans <- matrix(0, length(levels), ncol(mat), dimnames = list(levels, NULL))
  if (nrow(mat)) {
    pooled <- rowsum(mat, groups, reorder = FALSE)
    ans[match(rownames(pooled), levels), ] <- pooled
  }
  ans
}

calculate_stats <- function(observed, exact, null, features, category) {
  check(nrow(null) == length(features) && ncol(null) == 10000L,
        paste("Invalid randomization matrix:", category))
  n <- ncol(null)
  result <- data.table(category = category, feature = features,
                       observed = as.numeric(observed), expected_exact = as.numeric(exact),
                       expected = 0, fold_enrichment = NA_real_, log2_fold_enrichment = NA_real_,
                       null_sd = 0, null_q025 = 0, null_q975 = 0,
                       n_upper = n, n_lower = n, p_enrichment = 1, p_depletion = 1, p_two_sided = 1,
                       status = "zero_null_support")
  # Zero-only rows retain P = 1; exact support distinguishes finite-simulation zeros.
  keep <- which(exact > 0 | observed > 0 | rowSums(null) > 0)
  if (length(keep)) {
    x <- null[keep, , drop = FALSE]
    o <- observed[keep]
    mean_null <- rowMeans(x)
    upper <- rowSums(x >= o)
    lower <- rowSums(x <= o)
    p_upper <- (upper + 1) / (n + 1)
    p_lower <- (lower + 1) / (n + 1)
    fold <- o / mean_null
    fold[o == 0 & mean_null == 0] <- NA_real_
    intervals <- vapply(seq_len(nrow(x)), function(i)
      quantile(x[i, ], c(.025, .975), type = 1, names = FALSE), numeric(2))
    result[keep, `:=`(expected = mean_null, fold_enrichment = fold,
                      log2_fold_enrichment = log2(fold),
                      null_sd = vapply(seq_len(nrow(x)), function(i) sd(x[i, ]), numeric(1)),
                      null_q025 = intervals[1L, ], null_q975 = intervals[2L, ],
                      n_upper = upper, n_lower = lower,
                      p_enrichment = p_upper, p_depletion = p_lower,
                      p_two_sided = pmin(1, 2 * pmin(p_upper, p_lower)),
                      status = ifelse(exact[keep] > 0 & mean_null == 0, "zero_in_simulation", "tested"))]
  }
  result
}

cat("Reading selected cells and archived randomization metadata\n")
heatmap <- readRDS(heatmap_path)
selected <- names(heatmap$column_order)
check(length(selected) == 100L && !anyNA(selected) && !anyDuplicated(selected),
      "The heatmap selection must contain 100 uniquely named cell types")
annotation <- fread(file.path(run, "TE_subfamily_annotation.tsv"))
reference <- fread(file.path(run, "observed_TE_counts.tsv"), check.names = FALSE)
labels <- names(reference)[-1L]
check(nrow(annotation) == 1146L && length(labels) == 195L &&
        !anyDuplicated(annotation$annotated_subfamily) && !anyDuplicated(labels),
      "Unexpected annotation or reference dimensions")
check(all(selected %in% labels), paste("Missing selected cell types:",
                                       paste(setdiff(selected, labels), collapse = ", ")))
config <- readRDS(file.path(run, "run_config.rds"))
check(config$iterations == 10000L, "Expected 10,000 archived randomizations")
cell_summary <- fread(file.path(run, "celltype_summary.tsv"))
info <- as.data.table(readRDS(info_path))
check(all(c("tissue_cellType", "category") %in% names(info)), "Missing PT.info mapping fields")
map <- unique(info[, .(celltype = as.character(tissue_cellType), category = as.character(category))])
check(!anyNA(map) && !anyDuplicated(map$celltype), "Missing or conflicting category assignments")
check(all(selected %in% map$celltype) && all(selected %in% cell_summary$celltype),
      "Selected cell types missing from metadata")
map <- map[match(selected, celltype)]
check(all(map$category %in% order_categories), "Unexpected category")
map[, `:=`(heatmap_order = seq_along(selected), checkpoint_index = match(selected, labels))]
map <- merge(map, cell_summary[, .(celltype, n_foreground, n_candidates)],
             by = "celltype", sort = FALSE)
setorder(map, heatmap_order)
check(identical(map$celltype, selected), "Selection order changed")
fwrite(map, file.path(tables, "selected_celltype_category_mapping.tsv"), sep = "\t")

category_summary <- map[, .(n_celltypes = .N, n_active_celltypes = sum(n_foreground > 0L),
                             n_foreground = sum(n_foreground), n_candidate_calls = sum(n_candidates)),
                        by = category]
category_summary <- merge(data.table(category = order_categories), category_summary,
                          by = "category", all.x = TRUE, sort = FALSE)
for (col in c("n_celltypes", "n_active_celltypes", "n_foreground", "n_candidate_calls"))
  set(category_summary, which(is.na(category_summary[[col]])), col, 0)
category_summary[, category_index := match(category, order_categories)]
setorder(category_summary, category_index)
category_summary[, category_index := NULL]
fwrite(category_summary, file.path(tables, "category_summary.tsv"), sep = "\t")

annotation[, family_group := group_families(annotation)]
families <- sort(unique(annotation$family_group), method = "radix")
features <- annotation$annotated_subfamily
check(length(families) == 69L && setequal(unique(annotation$class),
  c("LINE", "SINE", "LTR", "DNA", "RC", "Other", "Unknown", "Satellite", "RNA",
    "LINE?", "SINE?", "LTR?", "DNA?", "RC?")), "Annotation groups changed")
fwrite(annotation, file.path(tables, "TE_subfamily_annotation.tsv"), sep = "\t")
checkpoint_paths <- file.path(run, "celltype_null", sprintf("%03d.rds", map$checkpoint_index))
paths <- c(heatmap_path, info_path, file.path(run, c("TE_subfamily_annotation.tsv",
             "observed_TE_counts.tsv", "run_config.rds", "celltype_summary.tsv")), checkpoint_paths)
hashes <- tools::md5sum(paths)
check(!anyNA(hashes), "One or more archived inputs are missing")
fwrite(data.table(path = names(hashes), md5 = unname(hashes)),
       file.path(out, "input_manifest.tsv"), sep = "\t")

sub_results <- vector("list", length(order_categories))
family_results <- vector("list", length(order_categories))
for (g in seq_along(order_categories)) {
  category <- order_categories[g]
  member_index <- which(map[["category"]] == category)
  members <- map[member_index]
  check(nrow(members) == category_summary$n_celltypes[match(category, category_summary$category)],
        paste("Incorrect selected cell count for", category))
  observed <- exact <- numeric(length(features))
  sub_null <- matrix(0L, length(features), config$iterations,
                     dimnames = list(features, NULL))
  for (j in seq_len(nrow(members))) {
    i <- members$checkpoint_index[j]
    saved <- readRDS(file.path(run, "celltype_null", sprintf("%03d.rds", i)))
    check(identical(saved$config, config) && identical(saved$celltype, labels[i]) &&
            identical(saved$result$annotated_subfamily, features), "Checkpoint provenance mismatch")
    check(identical(as.integer(saved$result$observed), as.integer(reference[[i + 1L]])),
          "Checkpoint observed counts do not match reference")
    observed <- observed + saved$result$observed
    ex <- saved$result$expected_exact
    ex[is.na(ex)] <- 0
    exact <- exact + ex
    if (!is.null(saved$null_counts)) {
      incoming <- saved$null_counts
      check(is.matrix(incoming) && ncol(incoming) == config$iterations &&
              !is.null(rownames(incoming)) && !anyDuplicated(rownames(incoming)),
            "Malformed checkpoint randomizations")
      idx <- match(rownames(incoming), features)
      check(!anyNA(idx), "Checkpoint contains unknown subfamily")
      sub_null[idx, ] <- sub_null[idx, , drop = FALSE] + incoming
    }
  }
  family_null <- group_counts(sub_null, annotation$family_group, families)
  family_obs <- group_counts(matrix(observed, ncol = 1L), annotation$family_group, families)[, 1L]
  family_exact <- group_counts(matrix(exact, ncol = 1L), annotation$family_group, families)[, 1L]
  check(sum(family_obs) == sum(observed) &&
          all(colSums(family_null) == colSums(sub_null)), "Family/subfamily totals disagree")
  sub_results[[g]] <- calculate_stats(observed, exact, sub_null, features, category)
  family_results[[g]] <- calculate_stats(family_obs, family_exact, family_null, families, category)
  supported <- which(exact > 0 | observed > 0 | rowSums(sub_null) > 0)
  saveRDS(list(category = category, selected_celltypes = members$celltype,
               subfamily_features = features, subfamily_observed = observed,
               subfamily_expected_exact = exact, supported_features = features[supported],
               subfamily_null_supported = sub_null[supported, , drop = FALSE],
               family_features = families, family_observed = family_obs, family_null = family_null,
               n_iterations = config$iterations),
          file.path(null_dir, paste0(sprintf("%02d", g), "_", category, ".rds")))
  cat("Completed category", category, "(", nrow(members), "cell types)\n")
  flush.console()
}

finish <- function(parts) {
  x <- rbindlist(parts)
  x[, `:=`(fdr_two_sided = p.adjust(p_two_sided, "BH"),
           fdr_enrichment = p.adjust(p_enrichment, "BH"),
           fdr_depletion = p.adjust(p_depletion, "BH"))]
  x[, fdr_two_sided_within_category := p.adjust(p_two_sided, "BH"), by = category]
  x <- merge(x, category_summary, by = "category", sort = FALSE)
  x[, n_iterations := config$iterations]
  x[, category_order := match(category, order_categories)]
  setorder(x, category_order, feature)
  x[, category_order := NULL]
  x
}
sub <- finish(sub_results)
sub <- merge(sub, annotation, by.x = "feature", by.y = "annotated_subfamily", sort = FALSE)
family <- finish(family_results)
check(nrow(sub) == 12L * 1146L && nrow(family) == 12L * 69L &&
        sum(sub$observed) == sum(family$observed), "Incomplete results")
check(identical(hashes, tools::md5sum(paths)), "Inputs changed during aggregation")
fwrite(sub, file.path(tables, "category_TE_subfamily_enrichment.tsv.gz"), sep = "\t", na = "NA")
fwrite(family, file.path(tables, "category_TE_family_enrichment.tsv"), sep = "\t", na = "NA")
