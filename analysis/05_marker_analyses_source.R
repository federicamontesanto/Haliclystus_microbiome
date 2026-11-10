# Haliclystus antarcticus 16S analysis.
# Reads host-only results created by 01_16S_DADA2_decontam_source.R.


############################
# 0. USER SETTINGS
############################

# Locate the repository automatically (Rscript, source(), or RStudio from root).
find_haliclystus_repo <- function() {
  candidates <- getwd()
  script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_file)) {
    candidates <- c(candidates, dirname(sub("^--file=", "", tail(script_file, 1L))))
  }
  sourced <- vapply(sys.frames(), function(f) {
    if (is.null(f$ofile)) "" else as.character(f$ofile)
  }, character(1))
  if (any(nzchar(sourced))) {
    candidates <- c(candidates, dirname(sourced[nzchar(sourced)]))
  }
  for (start in unique(candidates)) {
    p <- normalizePath(start, winslash = "/", mustWork = FALSE)
    for (i in 0:5) {
      if (file.exists(file.path(p, "metadata", "haliclystus_samples_65.csv"))) {
        return(normalizePath(p, winslash = "/", mustWork = TRUE))
      }
      parent <- dirname(p)
      if (parent == p) break
      p <- parent
    }
  }
  stop("Cannot find repository root. Open R/RStudio inside the cloned ",
       "Haliclystus_microbiome directory and run the script again.")
}
repo_dir <- find_haliclystus_repo()

base_dir <- repo_dir
output_dir <- file.path(repo_dir, "results", "16S")
rds_dir <- file.path(output_dir, "rds")
plots_dir <- file.path(output_dir, "plots")
tables_dir <- file.path(output_dir, "tables")

marker_dir <- file.path(output_dir, "updated_marker_analyses")
marker_plots_dir <- file.path(marker_dir, "plots")
marker_tables_dir <- file.path(marker_dir, "tables")
for (d in c(marker_dir, marker_plots_dir, marker_tables_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

poor_quality_samples <- c("CR2B", "CR3B", "CR4B", "CR5B")

# LEfSe-style settings.
lefse_firstalpha  <- 0.05
lefse_secondalpha <- 0.05
lefse_ldascore    <- 2
lefse_n_per_group <- 10
lefse_taxlevel_for_cladogram <- 6

# IndVal settings.
indval_permutations <- 999
indval_alpha <- 0.05
indval_n_per_group <- 10

# Pairwise compartment contrasts, with first group treated as positive.
compartment_pairs <- list(
  c("gut", "stalk"),
  c("gut", "tentacle"),
  c("stalk", "tentacle")
)

compartment_pal <- c(
  gut      = "#A51E2D",
  stalk    = "#F0B7A8",
  tentacle = "#177E89"
)
site_pal <- c(
  Baia_Byers = "#6E6E6E",
  Rizopatron = "#C46A1C"
)

############################
# 1. PACKAGES
############################

options(repos = c(CRAN = "https://cran.rstudio.com/"))

install_cran_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) install.packages(missing, dependencies = TRUE, type = "binary")
}

cran_pkgs <- c(
  "dplyr", "tidyr", "purrr", "stringr", "readr", "tibble", "ggplot2",
  "scales", "patchwork", "svglite", "forcats", "indicspecies", "permute",
  "phyloseq"
)
install_cran_if_missing(cran_pkgs)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", dependencies = TRUE, type = "binary")
if (!requireNamespace("MicrobiotaProcess", quietly = TRUE)) BiocManager::install("MicrobiotaProcess", update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(svglite)
  library(forcats)
  library(indicspecies)
  library(permute)
  library(phyloseq)
  library(MicrobiotaProcess)
})

############################
# 2. HELPERS
############################

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " ")))
  flush.console()
}

write_table <- function(obj, filename) {
  readr::write_csv(obj, file.path(marker_tables_dir, filename))
}

save_plot_multi <- function(plot_obj, filename_base, width = 8, height = 6) {
  ggplot2::ggsave(file.path(marker_plots_dir, paste0(filename_base, ".png")), plot_obj,
                  width = width, height = height, dpi = 450, bg = "white")
  ggplot2::ggsave(file.path(marker_plots_dir, paste0(filename_base, ".pdf")), plot_obj,
                  width = width, height = height, bg = "white", device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(marker_plots_dir, paste0(filename_base, ".svg")), plot_obj,
                  width = width, height = height, bg = "white", device = svglite::svglite)
}

clean_taxon <- function(x, fallback = "Unclassified") {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- fallback
  x
}

clean_taxon_label <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unclassified"
  x <- stringr::str_replace_all(x, "\\|", ";")
  x <- stringr::str_replace_all(x, "^.*;", "")
  x <- stringr::str_replace_all(x, "^[dkpcofgs]_{1,2}", "")
  x <- stringr::str_replace_all(x, "^[A-Za-z]:\\s*", "")
  x <- stringr::str_replace_all(x, "^g__", "")
  x <- stringr::str_replace_all(x, "^f__", "")
  x <- stringr::str_replace_all(x, "^o__", "")
  x <- stringr::str_replace_all(x, "^c__", "")
  x <- stringr::str_replace_all(x, "^p__", "")
  x <- stringr::str_squish(x)
  x[x == "" | is.na(x)] <- "Unclassified"
  x
}

# Pick the best label column from MicrobiotaProcess output. In some versions,
# Feature is an internal numeric ID and f contains the useful taxon name.
choose_label_column <- function(df) {
  # MicrobiotaProcess output differs by version. In some versions, `Feature` is
  # an internal numeric/node ID, whereas `f` or row names contain the readable
  # taxon label. Prefer columns containing taxonomic-looking text and penalise
  # purely numeric/internal IDs.
  candidates <- intersect(c("f", "taxon", "Taxon", "label", "node", "Feature", "feature", "rowname_from_object"), colnames(df))
  if (length(candidates) == 0) return(NA_character_)
  score_col <- function(col) {
    z <- as.character(df[[col]])
    z <- z[!is.na(z) & z != ""]
    if (length(z) == 0) return(-Inf)
    numeric_like <- stringr::str_detect(z, "^[-+]?[0-9.]+$")
    alpha_prop <- mean(stringr::str_detect(z, "[A-Za-z]"))
    non_numeric_prop <- mean(!numeric_like)
    tax_prefix_prop <- mean(stringr::str_detect(z, "(^|[;|])([dkpcofgs]__|[dkpcofgs]_)"))
    tax_word_prop <- mean(stringr::str_detect(stringr::str_to_lower(z), "bacteria|proteobacter|bacteroid|flavobacter|rhodobacter|polaribacter|sulfitobacter|stappiaceae|spongiibacter|genus|family|order|class|phylum"))
    mean_len <- mean(nchar(z))
    # Strongly penalise all-numeric columns.
    if (mean(numeric_like) > 0.8) return(-100)
    alpha_prop * 10 + non_numeric_prop * 5 + tax_prefix_prop * 15 + tax_word_prop * 8 + pmin(mean_len, 40) / 40
  }
  scores <- vapply(candidates, score_col, numeric(1))
  candidates[which.max(scores)]
}

theme_marker <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "grey25"),
      axis.ticks = ggplot2::element_line(linewidth = 0.30, colour = "grey25"),
      axis.text = ggplot2::element_text(colour = "grey15"),
      strip.background = ggplot2::element_rect(fill = "grey94", colour = "grey55", linewidth = 0.35),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
}

make_rank_count_matrix <- function(seqtab, tax, rank_col) {
  tax_vec <- as.character(tax[[rank_col]])
  tax_vec <- clean_taxon(tax_vec)
  tax_vec <- clean_taxon_label(tax_vec)
  names(tax_vec) <- rownames(tax)
  tax_vec <- tax_vec[colnames(seqtab)]
  out <- rowsum(t(seqtab), group = tax_vec, reorder = FALSE)
  out <- t(out)
  out
}

make_all_rank_relative_abundance_long <- function(seqtab, tax, meta) {
  rel <- sweep(seqtab, 1, rowSums(seqtab), "/")
  rel[is.na(rel)] <- 0
  tax_df <- as.data.frame(tax, stringsAsFactors = FALSE) |>
    tibble::rownames_to_column("ASV")

  rel_long <- as.data.frame(rel) |>
    tibble::rownames_to_column("SampleID") |>
    tidyr::pivot_longer(-SampleID, names_to = "ASV", values_to = "RelAbund") |>
    dplyr::filter(RelAbund > 0) |>
    dplyr::left_join(tax_df, by = "ASV")

  ranks <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(tax_df))

  out <- purrr::map_dfr(ranks, function(rank) {
    rel_long |>
      dplyr::mutate(
        Rank = rank,
        Taxon = clean_taxon_label(clean_taxon(.data[[rank]]))
      ) |>
      dplyr::filter(!Taxon %in% c("Unclassified", "NA")) |>
      dplyr::group_by(SampleID, Rank, Taxon) |>
      dplyr::summarise(RelAbund = sum(RelAbund), .groups = "drop")
  }) |>
    dplyr::left_join(
      tibble::tibble(
        SampleID = rownames(meta),
        Compartment = meta$Compartment,
        Site = meta$Site
      ),
      by = "SampleID"
    ) |>
    dplyr::mutate(TaxonClean = clean_taxon_label(Taxon))

  out
}

summarise_ra_for_features <- function(ra_long, feature_df, group_col, group_levels) {
  feature_map <- feature_df |>
    dplyr::distinct(TaxonLabel, TaxonClean)

  if (nrow(feature_map) == 0) return(tibble::tibble())

  out <- ra_long |>
    dplyr::mutate(Group = as.character(.data[[group_col]])) |>
    dplyr::filter(Group %in% group_levels) |>
    dplyr::inner_join(feature_map, by = "TaxonClean") |>
    # Keep TaxonClean through summarise(). The previous version grouped only by
    # SampleID, TaxonLabel and Group, which dropped TaxonClean. The later join
    # used TaxonClean and therefore failed with:
    # "Join columns in y must be present in the data: TaxonClean".
    dplyr::group_by(SampleID, TaxonLabel, TaxonClean, Group) |>
    dplyr::summarise(RelAbund = sum(RelAbund), .groups = "drop")

  sample_grid <- tibble::tibble(
    SampleID = rownames(meta_halo),
    Group = as.character(meta_halo[[group_col]])
  ) |>
    dplyr::filter(Group %in% group_levels) |>
    tidyr::crossing(feature_map)

  out <- sample_grid |>
    dplyr::left_join(out, by = c("SampleID", "Group", "TaxonLabel", "TaxonClean")) |>
    dplyr::mutate(RelAbund = tidyr::replace_na(RelAbund, 0)) |>
    dplyr::group_by(TaxonLabel, TaxonClean, Group) |>
    dplyr::summarise(
      mean_rel_abund = mean(RelAbund),
      sd_rel_abund = stats::sd(RelAbund),
      n = dplyr::n(),
      se_rel_abund = sd_rel_abund / sqrt(n),
      .groups = "drop"
    ) |>
    tidyr::complete(TaxonLabel, TaxonClean, Group = group_levels,
                    fill = list(mean_rel_abund = 0, sd_rel_abund = 0, n = 0, se_rel_abund = 0))

  out
}

extract_lefse_table <- function(diff_res, classgroup) {
  df <- tryCatch(as.data.frame(diff_res), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(tibble::tibble())
  df <- tibble::as_tibble(df, rownames = "rowname_from_object")

  label_col <- choose_label_column(df)
  if (is.na(label_col)) {
    df$FeatureLabel <- df$rowname_from_object
  } else {
    df$FeatureLabel <- as.character(df[[label_col]])
  }
  df$FeatureLabel <- clean_taxon_label(df$FeatureLabel)

  # If the chosen label column was still an internal numeric ID, try other
  # candidate columns and finally row names before accepting numeric labels.
  numeric_label_prop <- mean(stringr::str_detect(as.character(df$FeatureLabel), "^[-+]?[0-9.]+$"), na.rm = TRUE)
  if (is.finite(numeric_label_prop) && numeric_label_prop > 0.5) {
    fallback_candidates <- setdiff(intersect(c("f", "taxon", "Taxon", "label", "node", "Feature", "feature", "rowname_from_object"), colnames(df)), label_col)
    for (cand in fallback_candidates) {
      cand_label <- clean_taxon_label(as.character(df[[cand]]))
      cand_numeric_prop <- mean(stringr::str_detect(cand_label, "^[-+]?[0-9.]+$"), na.rm = TRUE)
      cand_alpha_prop <- mean(stringr::str_detect(cand_label, "[A-Za-z]"), na.rm = TRUE)
      if (is.finite(cand_numeric_prop) && cand_numeric_prop <= 0.5 && is.finite(cand_alpha_prop) && cand_alpha_prop > 0.5) {
        df$FeatureLabel <- cand_label
        break
      }
    }
  }

  lda_col <- intersect(c("LDAmean", "lda", "LDA", "ldascore", "effectsize"), colnames(df))[1]
  if (is.na(lda_col)) {
    numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    lda_col <- numeric_cols[1]
  }

  group_col <- if (classgroup %in% colnames(df)) classgroup else {
    hits <- intersect(c("Group", "group", "enrich_group", "enriched_group", "Sign_Group", "class", "Class"), colnames(df))
    if (length(hits) == 0) NA_character_ else hits[1]
  }

  p_col <- intersect(c("pvalue", "p_value", "Pvalue", "p"), colnames(df))[1]
  q_col <- intersect(c("fdr", "qvalue", "q_value", "padj"), colnames(df))[1]

  out <- df |>
    dplyr::mutate(
      LDA = as.numeric(.data[[lda_col]]),
      LDA_abs = abs(LDA),
      Group = if (!is.na(group_col)) as.character(.data[[group_col]]) else NA_character_,
      p_value = if (!is.na(p_col)) as.numeric(.data[[p_col]]) else NA_real_,
      q_value = if (!is.na(q_col)) as.numeric(.data[[q_col]]) else NA_real_,
      TaxonLabel = FeatureLabel,
      TaxonClean = clean_taxon_label(TaxonLabel)
    ) |>
    dplyr::filter(!is.na(LDA_abs), !is.na(TaxonLabel), TaxonLabel != "Unclassified") |>
    dplyr::arrange(dplyr::desc(LDA_abs))

  out
}

prepare_marker_top <- function(tbl, group_levels, n_per_group = 10, pair = NULL) {
  if (nrow(tbl) == 0) return(tibble::tibble())
  tbl <- tbl |>
    dplyr::mutate(
      Group = dplyr::recode(Group, gut_body = "gut", .default = Group),
      Group = factor(Group, levels = group_levels)
    ) |>
    dplyr::filter(!is.na(Group))

  out <- tbl |>
    dplyr::group_by(Group) |>
    dplyr::slice_max(LDA_abs, n = n_per_group, with_ties = FALSE) |>
    dplyr::ungroup()

  if (!is.null(pair)) {
    contrast <- paste(pair, collapse = "_vs_")
    out <- out |>
      dplyr::mutate(
        Contrast = contrast,
        signed_score = ifelse(as.character(Group) == pair[1], LDA_abs, -LDA_abs)
      )
  }

  out
}

plot_marker_ra_pair <- function(marker_df, ra_df, group_pal, score_col, score_label, title, subtitle = NULL, signed = FALSE) {
  if (nrow(marker_df) == 0) {
    return(ggplot2::ggplot() + theme_void() + ggplot2::labs(title = paste(title, "- no significant features")))
  }

  order_df <- marker_df |>
    dplyr::mutate(order_value = abs(.data[[score_col]])) |>
    dplyr::arrange(order_value) |>
    dplyr::pull(TaxonLabel)
  order_df <- unique(order_df)

  marker_df <- marker_df |>
    dplyr::mutate(TaxonLabel = factor(TaxonLabel, levels = order_df))

  ra_df <- ra_df |>
    dplyr::mutate(
      TaxonLabel = factor(TaxonLabel, levels = order_df),
      Group = factor(Group, levels = names(group_pal))
    )

  p_left <- ggplot2::ggplot(marker_df, ggplot2::aes(x = .data[[score_col]], y = TaxonLabel, colour = Group)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.35) +
    ggplot2::geom_point(size = 2.7, alpha = 0.95) +
    ggplot2::scale_colour_manual(values = group_pal, drop = FALSE) +
    theme_marker(9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "none"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = score_label, y = NULL)

  p_right <- ggplot2::ggplot(ra_df, ggplot2::aes(x = mean_rel_abund, y = TaxonLabel, fill = Group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.62, colour = "white", linewidth = 0.10) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = pmax(mean_rel_abund - se_rel_abund, 0), xmax = mean_rel_abund + se_rel_abund),
      position = ggplot2::position_dodge(width = 0.72), width = 0.18, linewidth = 0.28
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 0.1), expand = ggplot2::expansion(mult = c(0, 0.04))) +
    ggplot2::scale_fill_manual(values = group_pal, drop = FALSE) +
    theme_marker(9) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    ) +
    ggplot2::labs(title = "Mean relative abundance", x = "Relative abundance ± SE", y = NULL, fill = NULL)

  p_left + p_right + patchwork::plot_layout(widths = c(1.05, 1.2), guides = "collect")
}

plot_pairwise_summary <- function(df, score_col, score_label, title, group_pal) {
  if (nrow(df) == 0) return(ggplot2::ggplot() + theme_void() + ggplot2::labs(title = paste(title, "- no significant features")))
  df <- df |>
    dplyr::group_by(Contrast) |>
    dplyr::mutate(TaxonPanel = forcats::fct_reorder(TaxonLabel, .data[[score_col]])) |>
    dplyr::ungroup()

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[score_col]], y = TaxonPanel, colour = Group)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.35) +
    ggplot2::geom_point(size = 2.4, alpha = 0.95) +
    ggplot2::facet_wrap(~Contrast, scales = "free_y", nrow = 1) +
    ggplot2::scale_colour_manual(values = group_pal, drop = FALSE, name = "Enriched /\nindicator group") +
    theme_marker(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)) +
    ggplot2::labs(title = title, x = score_label, y = NULL)
}

############################
# 3. LOAD DATA AND ALIGN OBJECTS
############################

log_msg("Loading host data")
seqtab_halo_prev <- readRDS(file.path(rds_dir, "seqtab_haloclystus_prevfiltered.rds"))
taxa_halo_prev   <- readRDS(file.path(rds_dir, "taxa_haloclystus_prevfiltered.rds"))
meta_halo        <- readRDS(file.path(rds_dir, "metadata_haloclystus.rds")) |> as.data.frame()
ps_halo_tree     <- readRDS(file.path(rds_dir, "phyloseq_haloclystus_tree.rds"))

# The upstream input is restricted to the 65 host samples.
samples_to_exclude <- unique(poor_quality_samples)
found_to_exclude <- intersect(samples_to_exclude, meta_halo$Sample_ID)
log_msg("Excluding", length(found_to_exclude), "sample(s):", paste(found_to_exclude, collapse = ", "))

keep_mask <- !(meta_halo$Sample_ID %in% samples_to_exclude)
meta_halo <- meta_halo[keep_mask, , drop = FALSE]
if (!all(rownames(meta_halo) %in% rownames(seqtab_halo_prev)) && "Sample_ID" %in% colnames(meta_halo)) {
  rownames(meta_halo) <- meta_halo$Sample_ID
}
seqtab_halo_prev <- seqtab_halo_prev[rownames(meta_halo), , drop = FALSE]

taxa_halo_prev <- as.data.frame(taxa_halo_prev, stringsAsFactors = FALSE)
taxa_halo_prev <- taxa_halo_prev[colnames(seqtab_halo_prev), , drop = FALSE]

meta_halo$Compartment <- as.character(meta_halo$Compartment)
meta_halo$Compartment <- dplyr::recode(meta_halo$Compartment, gut_body = "gut", .default = meta_halo$Compartment)
meta_halo$Compartment <- factor(meta_halo$Compartment, levels = c("gut", "stalk", "tentacle"))
meta_halo$Site <- factor(meta_halo$Site, levels = c("Baia_Byers", "Rizopatron"))

ps_halo_tree <- phyloseq::prune_samples(rownames(meta_halo), ps_halo_tree)
phyloseq::sample_data(ps_halo_tree)$Compartment <- meta_halo[phyloseq::sample_names(ps_halo_tree), "Compartment"]
phyloseq::sample_data(ps_halo_tree)$Site <- meta_halo[phyloseq::sample_names(ps_halo_tree), "Site"]
phyloseq::sample_data(ps_halo_tree)$Compartment <- factor(phyloseq::sample_data(ps_halo_tree)$Compartment, levels = levels(meta_halo$Compartment))
phyloseq::sample_data(ps_halo_tree)$Site <- factor(phyloseq::sample_data(ps_halo_tree)$Site, levels = levels(meta_halo$Site))

all_rank_ra_long <- make_all_rank_relative_abundance_long(seqtab_halo_prev, taxa_halo_prev, meta_halo)
write_table(all_rank_ra_long, "06_all_rank_relative_abundance_long_for_marker_panels.csv")

############################
# 4. MULTI-CLASS LEfSe: COMPARTMENT AND SITE
############################

run_multiclass_lefse <- function(classgroup, group_pal, prefix, panel_title) {
  group_levels <- names(group_pal)
  log_msg("Running multi-class LEfSe for", classgroup)

  diff_res <- tryCatch({
    MicrobiotaProcess::diff_analysis(
      obj = ps_halo_tree,
      classgroup = classgroup,
      mlfun = "lda",
      filtermod = "fdr",
      firstcomfun = "kruskal.test",
      firstalpha = lefse_firstalpha,
      strictmod = TRUE,
      secondcomfun = "wilcox.test",
      subclmin = 3,
      subclwilcox = TRUE,
      secondalpha = lefse_secondalpha,
      ldascore = lefse_ldascore
    )
  }, error = function(e) {
    log_msg("LEfSe failed for", classgroup, ":", e$message)
    NULL
  })
  if (is.null(diff_res)) return(NULL)

  tbl <- extract_lefse_table(diff_res, classgroup = classgroup)
  write_table(tbl, paste0("06_lefse_results_multiclass_", prefix, ".csv"))

  group_summary <- tbl |>
    dplyr::mutate(Group = dplyr::recode(Group, gut_body = "gut", .default = Group)) |>
    dplyr::count(Group, name = "n_LEfSe_taxa") |>
    tidyr::complete(Group = group_levels, fill = list(n_LEfSe_taxa = 0))
  write_table(group_summary, paste0("06_lefse_group_summary_multiclass_", prefix, ".csv"))
  log_msg("LEfSe group summary for", classgroup, ":", paste(group_summary$Group, group_summary$n_LEfSe_taxa, sep = "=", collapse = "; "))

  top_df <- prepare_marker_top(tbl, group_levels = group_levels, n_per_group = lefse_n_per_group)
  write_table(top_df, paste0("06_lefse_top_features_multiclass_", prefix, ".csv"))

  if (nrow(top_df) > 0) {
    ra_df <- summarise_ra_for_features(all_rank_ra_long, top_df, group_col = classgroup, group_levels = group_levels)
    write_table(ra_df, paste0("06_lefse_relative_abundance_multiclass_", prefix, ".csv"))
    p <- plot_marker_ra_pair(
      marker_df = top_df,
      ra_df = ra_df,
      group_pal = group_pal,
      score_col = "LDA_abs",
      score_label = "LDA score",
      title = paste0(panel_title, " LEfSe taxa"),
      subtitle = paste0("Kruskal-Wallis FDR + Wilcoxon; LDA >= ", lefse_ldascore)
    )
    save_plot_multi(p, paste0("Fig_LEfSe_multiclass_", prefix, "_LDA_plus_relative_abundance"), width = 12.5, height = max(5.5, 0.32 * nrow(top_df) + 2.2))
  }

  p_clade <- tryCatch({
    MicrobiotaProcess::ggdiffclade(
      diff_res,
      alpha = 0.33,
      linewd = 0.20,
      skpointsize = 0.95,
      taxlevel = lefse_taxlevel_for_cladogram,
      cladetext = 2.1,
      setColors = TRUE
    ) +
      ggplot2::labs(title = paste0(panel_title, " LEfSe cladogram")) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }, error = function(e) {
    log_msg("LEfSe cladogram failed for", classgroup, ":", e$message)
    NULL
  })
  if (!is.null(p_clade)) save_plot_multi(p_clade, paste0("Fig_LEfSe_multiclass_", prefix, "_cladogram"), width = 9.5, height = 9.5)

  invisible(list(diff_res = diff_res, table = tbl, top = top_df))
}

lefse_compartment <- run_multiclass_lefse("Compartment", compartment_pal, "compartment", "Compartment")
lefse_site <- run_multiclass_lefse("Site", site_pal, "site", "Site")

############################
# 5. PAIRWISE COMPARTMENT LEfSe
############################

run_pairwise_lefse_one <- function(pair) {
  contrast <- paste(pair, collapse = "_vs_")
  log_msg("Running pairwise LEfSe:", contrast)
  sample_keep <- rownames(meta_halo)[as.character(meta_halo$Compartment) %in% pair]
  ps_pair <- phyloseq::prune_samples(sample_keep, ps_halo_tree)
  phyloseq::sample_data(ps_pair)$Compartment <- factor(as.character(phyloseq::sample_data(ps_pair)$Compartment), levels = pair)

  diff_res <- tryCatch({
    MicrobiotaProcess::diff_analysis(
      obj = ps_pair,
      classgroup = "Compartment",
      mlfun = "lda",
      filtermod = "fdr",
      firstcomfun = "kruskal.test",
      firstalpha = lefse_firstalpha,
      strictmod = TRUE,
      secondcomfun = "wilcox.test",
      subclmin = 3,
      subclwilcox = TRUE,
      secondalpha = lefse_secondalpha,
      ldascore = lefse_ldascore
    )
  }, error = function(e) {
    log_msg("Pairwise LEfSe failed for", contrast, ":", e$message)
    NULL
  })
  if (is.null(diff_res)) return(tibble::tibble())

  tbl <- extract_lefse_table(diff_res, classgroup = "Compartment") |>
    dplyr::mutate(Contrast = contrast)
  write_table(tbl, paste0("06_lefse_results_pairwise_", contrast, ".csv"))

  top_df <- prepare_marker_top(tbl, group_levels = pair, n_per_group = lefse_n_per_group, pair = pair)
  write_table(top_df, paste0("06_lefse_top_features_pairwise_", contrast, ".csv"))

  if (nrow(top_df) > 0) {
    ra_df <- summarise_ra_for_features(
      all_rank_ra_long |> dplyr::filter(as.character(Compartment) %in% pair),
      top_df,
      group_col = "Compartment",
      group_levels = pair
    )
    write_table(ra_df, paste0("06_lefse_relative_abundance_pairwise_", contrast, ".csv"))

    p <- plot_marker_ra_pair(
      marker_df = top_df,
      ra_df = ra_df,
      group_pal = compartment_pal[pair],
      score_col = "signed_score",
      score_label = paste0("Signed LDA score (+", pair[1], ", -", pair[2], ")"),
      title = paste0("LEfSe pairwise: ", contrast),
      subtitle = paste0("Kruskal-Wallis FDR + Wilcoxon; LDA >= ", lefse_ldascore)
    )
    save_plot_multi(p, paste0("Fig_LEfSe_pairwise_", contrast, "_LDA_plus_relative_abundance"), width = 12.5, height = max(5.0, 0.32 * nrow(top_df) + 2.2))
  }

  top_df
}

lefse_pairwise_top <- purrr::map_dfr(compartment_pairs, run_pairwise_lefse_one)
write_table(lefse_pairwise_top, "06_lefse_pairwise_top_features_all_contrasts.csv")

if (nrow(lefse_pairwise_top) > 0) {
  p_lefse_pair_summary <- plot_pairwise_summary(
    lefse_pairwise_top,
    score_col = "signed_score",
    score_label = "Signed LDA score",
    title = "Pairwise compartment LEfSe markers",
    group_pal = compartment_pal
  )
  save_plot_multi(p_lefse_pair_summary, "Fig_LEfSe_pairwise_compartment_signed_LDA_summary", width = 13.5, height = 6.8)
}

############################
# 6. INDVAL: SINGLE-COMPARTMENT AND PAIRWISE CONTRASTS
############################

log_msg("Building genus count matrix for IndVal")
genus_counts <- make_rank_count_matrix(seqtab_halo_prev, taxa_halo_prev, rank_col = "Genus")
genus_counts <- genus_counts[rownames(meta_halo), , drop = FALSE]
genus_rel <- sweep(genus_counts, 1, rowSums(genus_counts), "/")
genus_rel[is.na(genus_rel)] <- 0

genus_filter_tbl <- tibble::tibble(
  Genus = colnames(genus_counts),
  prevalence_prop = colMeans(genus_counts > 0),
  mean_rel_abund = colMeans(genus_rel)
) |>
  dplyr::arrange(dplyr::desc(mean_rel_abund))
write_table(genus_filter_tbl, "06_indval_genus_filter_summary.csv")
indval_keep <- genus_filter_tbl |>
  dplyr::filter(prevalence_prop > 0, mean_rel_abund > 0) |>
  dplyr::pull(Genus)
genus_counts_indval <- genus_counts[, indval_keep, drop = FALSE]

parse_multipatt_sign <- function(mp, group_levels, alpha = 0.05, contrast = NA_character_, pair = NULL) {
  tbl <- tibble::as_tibble(mp$sign, rownames = "Taxon")
  s_cols <- grep("^s\\.", colnames(tbl), value = TRUE)
  if (length(s_cols) == 0) stop("No s.* columns found in multipatt output.")
  s_group_names <- stringr::str_remove(s_cols, "^s\\.")
  s_mat <- tbl |> dplyr::select(dplyr::all_of(s_cols)) |> as.data.frame()
  single_group <- apply(s_mat, 1, function(z) {
    z <- as.numeric(z)
    if (all(is.na(z)) || max(z, na.rm = TRUE) <= 0) return(NA_character_)
    s_group_names[which.max(z)]
  })

  out <- tbl |>
    dplyr::mutate(
      Group = single_group,
      Group = dplyr::recode(Group, gut_body = "gut", .default = Group),
      TaxonLabel = clean_taxon_label(Taxon),
      TaxonClean = clean_taxon_label(Taxon),
      IndVal = stat,
      p_value = p.value,
      q_value = p.adjust(p.value, method = "BH"),
      Contrast = contrast
    ) |>
    dplyr::filter(!is.na(Group), Group %in% group_levels) |>
    dplyr::arrange(p_value, dplyr::desc(IndVal))

  if (!is.null(pair)) {
    out <- out |>
      dplyr::mutate(signed_score = ifelse(Group == pair[1], IndVal, -IndVal))
  }

  list(all = out, sig = out |> dplyr::filter(p_value <= alpha))
}

log_msg("Running multi-class single-compartment IndVal with duleg = TRUE")
indval_single <- indicspecies::multipatt(
  genus_counts_indval,
  meta_halo$Compartment,
  func = "IndVal.g",
  duleg = TRUE,
  control = permute::how(nperm = indval_permutations)
)
indval_parsed <- parse_multipatt_sign(indval_single, group_levels = levels(meta_halo$Compartment), alpha = indval_alpha)
indval_single_all <- indval_parsed$all
indval_single_sig <- indval_parsed$sig
write_table(indval_single_all, "06_indicator_species_genus_single_compartment_duleg_TRUE_all.csv")
write_table(indval_single_sig, "06_indicator_species_genus_single_compartment_duleg_TRUE_significant.csv")

indval_group_summary <- indval_single_sig |>
  dplyr::count(Group, name = "n_indicator_genera") |>
  tidyr::complete(Group = names(compartment_pal), fill = list(n_indicator_genera = 0))
write_table(indval_group_summary, "06_indicator_species_group_summary_duleg_TRUE.csv")
log_msg("IndVal single-compartment summary:", paste(indval_group_summary$Group, indval_group_summary$n_indicator_genera, sep = "=", collapse = "; "))

indval_top <- indval_single_sig |>
  dplyr::group_by(Group) |>
  dplyr::slice_max(IndVal, n = indval_n_per_group, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(Group = factor(Group, levels = names(compartment_pal)))
write_table(indval_top, "06_indicator_species_top_single_compartment_for_RA_panel.csv")

if (nrow(indval_top) > 0) {
  ra_indval <- summarise_ra_for_features(all_rank_ra_long, indval_top, group_col = "Compartment", group_levels = names(compartment_pal))
  write_table(ra_indval, "06_indicator_species_relative_abundance_single_compartment.csv")
  p_indval <- plot_marker_ra_pair(
    marker_df = indval_top,
    ra_df = ra_indval,
    group_pal = compartment_pal,
    score_col = "IndVal",
    score_label = "Indicator value",
    title = "Single-compartment indicator genera",
    subtitle = paste0("IndVal.g, duleg = TRUE, p <= ", indval_alpha, "; relative abundance shown as mean ± SE")
  )
  save_plot_multi(p_indval, "Fig_IndVal_genus_single_compartment_stat_plus_relative_abundance", width = 12.5, height = max(5.2, 0.32 * nrow(indval_top) + 2.2))
}

run_pairwise_indval_one <- function(pair) {
  contrast <- paste(pair, collapse = "_vs_")
  log_msg("Running pairwise IndVal:", contrast)
  sample_keep <- rownames(meta_halo)[as.character(meta_halo$Compartment) %in% pair]
  mat_pair <- genus_counts_indval[sample_keep, , drop = FALSE]
  grp_pair <- factor(as.character(meta_halo[sample_keep, "Compartment"]), levels = pair)

  # Remove columns that are all zero in this pair.
  mat_pair <- mat_pair[, colSums(mat_pair) > 0, drop = FALSE]
  if (ncol(mat_pair) < 2) {
    log_msg("Skipping pairwise IndVal", contrast, "because fewer than 2 genera remain.")
    return(tibble::tibble())
  }

  mp <- indicspecies::multipatt(
    mat_pair,
    grp_pair,
    func = "IndVal.g",
    duleg = TRUE,
    control = permute::how(nperm = indval_permutations)
  )
  parsed <- parse_multipatt_sign(mp, group_levels = pair, alpha = indval_alpha, contrast = contrast, pair = pair)
  all_tbl <- parsed$all
  sig_tbl <- parsed$sig
  write_table(all_tbl, paste0("06_indicator_species_pairwise_", contrast, "_all.csv"))
  write_table(sig_tbl, paste0("06_indicator_species_pairwise_", contrast, "_significant.csv"))

  top_df <- sig_tbl |>
    dplyr::group_by(Group) |>
    dplyr::slice_max(IndVal, n = indval_n_per_group, with_ties = FALSE) |>
    dplyr::ungroup()
  write_table(top_df, paste0("06_indicator_species_pairwise_", contrast, "_top_for_RA_panel.csv"))

  if (nrow(top_df) > 0) {
    ra_df <- summarise_ra_for_features(
      all_rank_ra_long |> dplyr::filter(as.character(Compartment) %in% pair),
      top_df,
      group_col = "Compartment",
      group_levels = pair
    )
    write_table(ra_df, paste0("06_indicator_species_pairwise_", contrast, "_relative_abundance.csv"))
    p <- plot_marker_ra_pair(
      marker_df = top_df,
      ra_df = ra_df,
      group_pal = compartment_pal[pair],
      score_col = "signed_score",
      score_label = paste0("Signed IndVal (+", pair[1], ", -", pair[2], ")"),
      title = paste0("Pairwise IndVal: ", contrast),
      subtitle = paste0("IndVal.g, duleg = TRUE, p <= ", indval_alpha)
    )
    save_plot_multi(p, paste0("Fig_IndVal_pairwise_", contrast, "_stat_plus_relative_abundance"), width = 12.5, height = max(5.0, 0.32 * nrow(top_df) + 2.2))
  }

  top_df
}

indval_pairwise_top <- purrr::map_dfr(compartment_pairs, run_pairwise_indval_one)
write_table(indval_pairwise_top, "06_indicator_species_pairwise_top_all_contrasts.csv")

if (nrow(indval_pairwise_top) > 0) {
  p_indval_pair_summary <- plot_pairwise_summary(
    indval_pairwise_top,
    score_col = "signed_score",
    score_label = "Signed IndVal",
    title = "Pairwise compartment indicator genera",
    group_pal = compartment_pal
  )
  save_plot_multi(p_indval_pair_summary, "Fig_IndVal_pairwise_compartment_signed_summary", width = 13.5, height = 6.8)
}

############################
# 7. FINISH
############################

notes_file <- file.path(marker_tables_dir, "06_marker_analysis_interpretation_notes_v3.txt")
cat(
  "Updated LEfSe + IndVal marker analysis v3\n",
  "========================================\n\n",
  "This script produces both multi-class and pairwise compartment marker analyses.\n",
  "Multi-class LEfSe/IndVal answer: which groups have discriminant/indicator taxa across all compartments?\n",
  "Pairwise LEfSe/IndVal answer: which taxa distinguish gut_vs_stalk, gut_vs_tentacle, and stalk_vs_tentacle?\n\n",
  "Pairwise marker analyses are useful for visual comparison with DESeq2 contrasts, but should be interpreted as exploratory marker support.\n",
  "Relative abundance panels show mean relative abundance ± SE for the same taxa shown in the marker-score panel.\n",
  file = notes_file
)

log_msg("Updated LEfSe + IndVal pairwise analysis complete.")
log_msg("Tables written to:", marker_tables_dir)
log_msg("Plots written to:", marker_plots_dir)
log_msg("Key pairwise plots:")
log_msg(" - Fig_LEfSe_pairwise_compartment_signed_LDA_summary.*")
log_msg(" - Fig_IndVal_pairwise_compartment_signed_summary.*")
log_msg(" - Fig_LEfSe_pairwise_<contrast>_LDA_plus_relative_abundance.*")
log_msg(" - Fig_IndVal_pairwise_<contrast>_stat_plus_relative_abundance.*")

############################
# COMBINED DESEQ2 PAIRWISE PLOT (3 PANELS IN ONE ROW)
############################

# keep only significant results already prepared in your object
# here I assume you already have a data frame like deseq_sig_df
# with columns:
# comparison, taxon_label, log2FoldChange, enriched_group

deseq_sig_df <- deseq_sig_df |>
  dplyr::mutate(
    comparison = factor(
      comparison,
      levels = c("gut_vs_stalk", "gut_vs_tentacle", "stalk_vs_tentacle")
    ),
    taxon_label = factor(taxon_label, levels = rev(unique(taxon_label)))
  )

# use a fixed x-range so the three panels align nicely
xmax_deseq <- max(abs(deseq_sig_df$log2FoldChange), na.rm = TRUE)
xmax_deseq <- ceiling(xmax_deseq / 5) * 5

p_deseq_pairwise_row <- ggplot2::ggplot(
  deseq_sig_df,
  ggplot2::aes(x = log2FoldChange, y = taxon_label, colour = enriched_group)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::facet_wrap(~ comparison, nrow = 1, scales = "fixed") +
  ggplot2::scale_x_continuous(
    limits = c(-xmax_deseq, xmax_deseq),
    breaks = seq(-xmax_deseq, xmax_deseq, by = 10),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::scale_colour_manual(values = compartment_pal, name = "Enriched in") +
  theme_fed(base_size = 10) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    axis.text.x = ggplot2::element_text(size = 8),
    axis.text.y = ggplot2::element_text(size = 8),
    axis.title.x = ggplot2::element_text(size = 10),
    axis.title.y = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(0.35, "lines"),
    legend.position = "right",
    plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
  ) +
  ggplot2::labs(
    title = "DESeq2 significant ASVs",
    subtitle = "padj <= 0.05",
    x = "log2 fold change"
  )

save_plot_multi(
  p_deseq_pairwise_row,
  "Fig_DESeq2_pairwise_ASVs_row",
  width = 14,
  height = 7
)
############################
# COMBINED PAIRWISE LEFSE PLOT (3 PANELS IN ONE ROW)
############################

# expected columns:
# comparison, TaxonClean, lda_score_signed, enriched_group

lefse_pairwise_plot_df <- lefse_pairwise_plot_df |>
  dplyr::mutate(
    comparison = factor(
      comparison,
      levels = c("gut_vs_stalk", "gut_vs_tentacle", "stalk_vs_tentacle")
    ),
    TaxonClean = factor(TaxonClean, levels = rev(unique(TaxonClean)))
  )

xmax_lefse <- max(abs(lefse_pairwise_plot_df$lda_score_signed), na.rm = TRUE)
xmax_lefse <- ceiling(xmax_lefse)

p_lefse_pairwise_row <- ggplot2::ggplot(
  lefse_pairwise_plot_df,
  ggplot2::aes(x = lda_score_signed, y = TaxonClean, colour = enriched_group)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::facet_wrap(~ comparison, nrow = 1, scales = "fixed") +
  ggplot2::scale_x_continuous(
    limits = c(-xmax_lefse, xmax_lefse),
    breaks = pretty(c(-xmax_lefse, xmax_lefse), n = 5),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::scale_colour_manual(values = compartment_pal, name = "Enriched in") +
  theme_fed(base_size = 10) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    axis.text.x = ggplot2::element_text(size = 8),
    axis.text.y = ggplot2::element_text(size = 8),
    axis.title.x = ggplot2::element_text(size = 10),
    axis.title.y = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(0.35, "lines"),
    legend.position = "right",
    plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
  ) +
  ggplot2::labs(
    title = "Pairwise LEfSe taxa",
    subtitle = "Kruskal-Wallis/Wilcoxon-derived pairwise discriminatory taxa",
    x = "Signed LDA score"
  )

save_plot_multi(
  p_lefse_pairwise_row,
  "Fig_LEfSe_pairwise_row",
  width = 14,
  height = 7
)
############################
# COMBINED PAIRWISE INDVAL PLOT (3 PANELS IN ONE ROW)
############################

# expected columns:
# comparison, TaxonClean, stat_signed, enriched_group

indval_pairwise_plot_df <- indval_pairwise_plot_df |>
  dplyr::mutate(
    comparison = factor(
      comparison,
      levels = c("gut_vs_stalk", "gut_vs_tentacle", "stalk_vs_tentacle")
    ),
    TaxonClean = factor(TaxonClean, levels = rev(unique(TaxonClean)))
  )

xmax_indval <- max(abs(indval_pairwise_plot_df$stat_signed), na.rm = TRUE)
xmax_indval <- ceiling(xmax_indval * 10) / 10

p_indval_pairwise_row <- ggplot2::ggplot(
  indval_pairwise_plot_df,
  ggplot2::aes(x = stat_signed, y = TaxonClean, colour = enriched_group)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::facet_wrap(~ comparison, nrow = 1, scales = "fixed") +
  ggplot2::scale_x_continuous(
    limits = c(-xmax_indval, xmax_indval),
    breaks = pretty(c(-xmax_indval, xmax_indval), n = 5),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::scale_colour_manual(values = compartment_pal, name = "Enriched in") +
  theme_fed(base_size = 10) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    axis.text.x = ggplot2::element_text(size = 8),
    axis.text.y = ggplot2::element_text(size = 8),
    axis.title.x = ggplot2::element_text(size = 10),
    axis.title.y = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(0.35, "lines"),
    legend.position = "right",
    plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
  ) +
  ggplot2::labs(
    title = "Pairwise indicator genera",
    subtitle = "IndVal pairwise summary",
    x = "Signed indicator value"
  )

save_plot_multi(
  p_indval_pairwise_row,
  "Fig_IndVal_pairwise_row",
  width = 14,
  height = 7
)
