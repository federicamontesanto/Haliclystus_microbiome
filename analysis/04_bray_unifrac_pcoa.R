# Haliclystus antarcticus 16S analysis.
################################################################################

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
rds_dir    <- file.path(output_dir, "rds")

ordination_dir <- file.path(output_dir, "ordination_bray_unifrac_pcoa")
plots_dir      <- file.path(ordination_dir, "plots")
tables_dir     <- file.path(ordination_dir, "tables")

for (d in c(ordination_dir, plots_dir, tables_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# Which ellipses to draw on ordination plots.
# Usually "Compartment" is best for your manuscript figure.
ellipse_group <- "Compartment"
ellipse_linetype <- "dotted"
ellipse_level <- 0.68       # 0.68 is compact; set 0.95 for wider conventional ellipses
min_ellipse_n <- 3          # no ellipse drawn for groups with fewer than this many samples

# UniFrac choices.
run_weighted_unifrac   <- TRUE
run_unweighted_unifrac <- TRUE

# PERMANOVA settings.
n_permutations <- 9999

# Plot styling.
point_size <- 3.2
point_alpha <- 0.92
save_width_single <- 7.2
save_height_single <- 5.8
save_width_combined <- 13.2
save_height_combined <- 5.8

# Colors used consistently with the rest of the project.
compartment_pal <- c(
  gut      = "#A51E2D",
  stalk    = "#F0B7A8",
  tentacle = "#177E89"
)
site_shape_values <- c(
  Baia_Byers = 21,
  Rizopatron = 24
)
site_pal <- c(
  Baia_Byers = "#5E5E5E",
  Rizopatron = "#1F1F1F"
)

############################
# 1. PACKAGES
############################

options(repos = c(CRAN = "https://cran.rstudio.com/"))

install_cran_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, dependencies = TRUE, type = "binary")
  }
}

cran_pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "readr", "tibble", "stringr",
  "vegan", "ape", "phyloseq", "patchwork", "scales", "svglite", "purrr"
)
install_cran_if_missing(cran_pkgs)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(stringr)
  library(vegan)
  library(ape)
  library(phyloseq)
  library(patchwork)
  library(scales)
  library(svglite)
  library(purrr)
})

############################
# 2. HELPERS
############################

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " ")))
  flush.console()
}

write_table <- function(obj, filename) {
  readr::write_csv(obj, file.path(tables_dir, filename))
}

save_plot_multi <- function(plot_obj, filename_base, width = 8, height = 6) {
  ggplot2::ggsave(
    filename = file.path(plots_dir, paste0(filename_base, ".png")),
    plot = plot_obj, width = width, height = height, dpi = 450, bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path(plots_dir, paste0(filename_base, ".pdf")),
    plot = plot_obj, width = width, height = height, bg = "white", device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    filename = file.path(plots_dir, paste0(filename_base, ".svg")),
    plot = plot_obj, width = width, height = height, bg = "white", device = svglite::svglite
  )
}

theme_ordination <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "grey20"),
      axis.ticks = ggplot2::element_line(linewidth = 0.30, colour = "grey20"),
      axis.text = ggplot2::element_text(colour = "grey15"),
      axis.title = ggplot2::element_text(colour = "grey10"),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_size - 1),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

standardise_compartment <- function(x) {
  x <- as.character(x)
  dplyr::recode(x, gut_body = "gut", .default = x)
}

make_individual_id <- function(meta) {
  if ("Individual_ID" %in% colnames(meta)) {
    return(factor(meta$Individual_ID))
  }

  digits <- stringr::str_extract(meta$Sample_ID, "\\d+")
  if (any(is.na(digits))) {
    bad <- meta$Sample_ID[is.na(digits)]
    stop(
      "Sample_ID(s) contain no digits, so an individual number could not be extracted: ",
      paste(bad, collapse = ", "),
      "\nCheck these IDs manually or add Individual_ID to metadata."
    )
  }
  factor(paste0(meta$Site, "_", digits))
}

as_dist_object <- function(x) {
  if (inherits(x, "dist")) return(x)
  stats::as.dist(as.matrix(x))
}

# Convert distance matrix to PCoA scores with metadata.
make_pcoa_scores <- function(dist_obj, meta, label) {
  dist_obj <- as_dist_object(dist_obj)
  ord <- ape::pcoa(as.matrix(dist_obj))

  vectors <- as.data.frame(ord$vectors[, 1:2, drop = FALSE]) |>
    tibble::rownames_to_column("SampleID") |>
    dplyr::rename(Axis1 = Axis.1, Axis2 = Axis.2) |>
    dplyr::left_join(
      meta |>
        tibble::rownames_to_column("SampleID"),
      by = "SampleID"
    ) |>
    dplyr::mutate(distance = label)

  eig_tbl <- ord$values |>
    as.data.frame() |>
    tibble::rownames_to_column("axis") |>
    dplyr::mutate(distance = label)

  # ape::pcoa() provides Relative_eig. Use it if available.
  rel_eig <- ord$values$Relative_eig
  if (is.null(rel_eig) || length(rel_eig) < 2 || any(!is.finite(rel_eig[1:2]))) {
    positive_eig <- ord$values$Eigenvalues
    positive_eig[positive_eig < 0] <- 0
    rel_eig <- positive_eig / sum(positive_eig)
  }

  list(
    ord = ord,
    scores = vectors,
    eig = eig_tbl,
    pct1 = round(100 * rel_eig[1], 1),
    pct2 = round(100 * rel_eig[2], 1)
  )
}

plot_pcoa <- function(scores, pct1, pct2, title, subtitle, filename_base = NULL) {
  scores <- scores |>
    dplyr::mutate(
      Compartment = factor(as.character(Compartment), levels = c("gut", "stalk", "tentacle")),
      Site = factor(as.character(Site), levels = c("Baia_Byers", "Rizopatron"))
    )

  ellipse_df <- scores |>
    dplyr::filter(!is.na(.data[[ellipse_group]])) |>
    dplyr::group_by(.data[[ellipse_group]]) |>
    dplyr::filter(dplyr::n() >= min_ellipse_n) |>
    dplyr::ungroup()

  p <- ggplot2::ggplot(
    scores,
    ggplot2::aes(x = Axis1, y = Axis2, colour = Compartment, shape = Site)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey88") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey88")

  if (nrow(ellipse_df) > 0) {
    p <- p +
      ggplot2::stat_ellipse(
        data = ellipse_df,
        ggplot2::aes(group = .data[[ellipse_group]], colour = Compartment),
        type = "t",
        level = ellipse_level,
        linewidth = 0.75,
        linetype = ellipse_linetype,
        alpha = 0.85,
        show.legend = FALSE
      )
  }

  p <- p +
    ggplot2::geom_point(size = point_size, alpha = point_alpha, stroke = 0.75) +
    ggplot2::scale_colour_manual(values = compartment_pal, drop = FALSE, name = "Compartment") +
    ggplot2::scale_shape_manual(values = site_shape_values, drop = FALSE, name = "Site") +
    theme_ordination(10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = paste0("PCoA1 (", pct1, "%)"),
      y = paste0("PCoA2 (", pct2, "%)")
    )

  if (!is.null(filename_base)) {
    save_plot_multi(p, filename_base, width = save_width_single, height = save_height_single)
  }

  p
}

run_permanova <- function(dist_obj, meta, label) {
  dist_obj <- as_dist_object(dist_obj)
  meta <- meta[labels(dist_obj), , drop = FALSE]

  # Unrestricted model: useful for overall Site + Compartment summary, but does
  # not respect repeated measures.
  perm_unrestricted <- vegan::adonis2(
    dist_obj ~ Site + Compartment,
    data = meta,
    permutations = n_permutations,
    by = "margin"
  ) |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    dplyr::mutate(distance = label, permutation_scheme = "unrestricted")

  # Repeated-measures-aware model for Compartment: permutations blocked within
  # Individual_ID. Site cannot be cleanly tested with Individual_ID blocks when
  # Individual_ID is nested in Site, so interpret the blocked output primarily
  # for Compartment.
  perm_ctrl <- permute::how(nperm = n_permutations, blocks = meta$Individual_ID)
  perm_blocked <- vegan::adonis2(
    dist_obj ~ Compartment,
    data = meta,
    permutations = perm_ctrl,
    by = "margin"
  ) |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    dplyr::mutate(distance = label, permutation_scheme = "blocked_by_Individual_ID")

  dplyr::bind_rows(perm_unrestricted, perm_blocked)
}

run_betadisper <- function(dist_obj, meta, label, group_col) {
  dist_obj <- as_dist_object(dist_obj)
  meta <- meta[labels(dist_obj), , drop = FALSE]
  group <- droplevels(factor(meta[[group_col]]))

  out <- tryCatch({
    bd <- vegan::betadisper(dist_obj, group = group)
    an <- anova(bd) |>
      as.data.frame() |>
      tibble::rownames_to_column("term") |>
      dplyr::mutate(distance = label, group = group_col)
    tuk <- tryCatch({
      tk <- TukeyHSD(bd)
      as.data.frame(tk$group) |>
        tibble::rownames_to_column("contrast") |>
        dplyr::mutate(distance = label, group = group_col)
    }, error = function(e) tibble::tibble())
    list(anova = an, tukey = tuk)
  }, error = function(e) {
    log_msg("betadisper failed for", label, group_col, ":", e$message)
    list(anova = tibble::tibble(), tukey = tibble::tibble())
  })

  out
}

############################
# 3. LOAD AND ALIGN HOST DATA
############################

log_msg("Loading RDS objects from", rds_dir)
seqtab_halo_prev <- readRDS(file.path(rds_dir, "seqtab_haloclystus_prevfiltered.rds"))
taxa_halo_prev   <- readRDS(file.path(rds_dir, "taxa_haloclystus_prevfiltered.rds"))
meta_halo        <- readRDS(file.path(rds_dir, "metadata_haloclystus.rds")) |> as.data.frame()
ps_halo_tree     <- readRDS(file.path(rds_dir, "phyloseq_haloclystus_tree.rds"))

# Ensure sample IDs are row names if needed.
if (!all(rownames(meta_halo) %in% rownames(seqtab_halo_prev)) && "Sample_ID" %in% colnames(meta_halo)) {
  rownames(meta_halo) <- meta_halo$Sample_ID
}

# The upstream input is restricted to the 65 host samples.
samples_to_exclude <- unique(poor_quality_samples)
found_to_exclude <- intersect(samples_to_exclude, meta_halo$Sample_ID)
log_msg("Excluding", length(found_to_exclude), "sample(s):", paste(found_to_exclude, collapse = ", "))

meta_halo <- meta_halo[!(meta_halo$Sample_ID %in% samples_to_exclude), , drop = FALSE]
seqtab_halo_prev <- seqtab_halo_prev[rownames(meta_halo), , drop = FALSE]

taxa_halo_prev <- as.data.frame(taxa_halo_prev, stringsAsFactors = FALSE)
taxa_halo_prev <- taxa_halo_prev[colnames(seqtab_halo_prev), , drop = FALSE]

meta_halo$Compartment <- standardise_compartment(meta_halo$Compartment)
meta_halo$Compartment <- factor(meta_halo$Compartment, levels = c("gut", "stalk", "tentacle"))
meta_halo$Site <- factor(meta_halo$Site, levels = c("Baia_Byers", "Rizopatron"))
meta_halo$Individual_ID <- make_individual_id(meta_halo)

# Report repeated samples per individual/compartment.
dup_check <- meta_halo |>
  dplyr::count(Individual_ID, Compartment) |>
  dplyr::filter(n > 1)
write_table(dup_check, "07_individual_compartment_replicates_for_ordination.csv")

individual_check <- meta_halo |>
  dplyr::count(Individual_ID, Compartment) |>
  tidyr::pivot_wider(names_from = Compartment, values_from = n, values_fill = 0)
write_table(individual_check, "07_CHECK_individual_id_assignment_for_ordination.csv")

# Align phyloseq object to retained samples and update sample_data.
common_ps_samples <- intersect(phyloseq::sample_names(ps_halo_tree), rownames(meta_halo))
if (length(common_ps_samples) == 0) {
  stop("No retained metadata samples matched sample_names(ps_halo_tree). Check sample IDs.")
}
if (length(common_ps_samples) < nrow(meta_halo)) {
  log_msg("WARNING:", nrow(meta_halo) - length(common_ps_samples), "metadata sample(s) absent from ps_halo_tree and will be omitted from UniFrac.")
}

ps_host <- phyloseq::prune_samples(common_ps_samples, ps_halo_tree)
ps_host <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_host) > 0, ps_host)

# Update sample_data in ps_host with cleaned metadata.
meta_for_ps <- meta_halo[phyloseq::sample_names(ps_host), , drop = FALSE]
phyloseq::sample_data(ps_host)$Compartment <- meta_for_ps$Compartment
phyloseq::sample_data(ps_host)$Site <- meta_for_ps$Site
phyloseq::sample_data(ps_host)$Individual_ID <- meta_for_ps$Individual_ID
phyloseq::sample_data(ps_host)$Sample_ID <- meta_for_ps$Sample_ID

# Bray-Curtis can use the cleaned count table directly.
meta_for_dist <- meta_halo[rownames(seqtab_halo_prev), , drop = FALSE]

write_table(
  tibble::tibble(
    SampleID = rownames(meta_for_dist),
    Sample_ID = meta_for_dist$Sample_ID,
    Site = as.character(meta_for_dist$Site),
    Compartment = as.character(meta_for_dist$Compartment),
    Individual_ID = as.character(meta_for_dist$Individual_ID),
    read_depth = rowSums(seqtab_halo_prev)
  ),
  "07_ordination_sample_metadata_used.csv"
)

############################
# 4. DISTANCE MATRICES
############################

log_msg("Computing Bray-Curtis distance")
dist_bray <- vegan::vegdist(seqtab_halo_prev, method = "bray")

# Ensure meta order follows distance labels.
meta_bray <- meta_for_dist[labels(dist_bray), , drop = FALSE]

distances <- list(BrayCurtis = dist_bray)
meta_by_distance <- list(BrayCurtis = meta_bray)

if (run_weighted_unifrac) {
  log_msg("Computing weighted UniFrac distance")
  dist_wu <- phyloseq::UniFrac(ps_host, weighted = TRUE, normalized = TRUE, fast = TRUE)
  dist_wu <- as_dist_object(dist_wu)
  distances$WeightedUniFrac <- dist_wu
  meta_by_distance$WeightedUniFrac <- meta_halo[labels(dist_wu), , drop = FALSE]
}

if (run_unweighted_unifrac) {
  log_msg("Computing unweighted UniFrac distance")
  dist_uwu <- phyloseq::UniFrac(ps_host, weighted = FALSE, normalized = TRUE, fast = TRUE)
  dist_uwu <- as_dist_object(dist_uwu)
  distances$UnweightedUniFrac <- dist_uwu
  meta_by_distance$UnweightedUniFrac <- meta_halo[labels(dist_uwu), , drop = FALSE]
}

############################
# 5. PCoA ORDINATION PLOTS
############################

ordination_results <- list()
plots <- list()

for (nm in names(distances)) {
  pretty_label <- dplyr::case_when(
    nm == "BrayCurtis" ~ "Bray-Curtis dissimilarity",
    nm == "WeightedUniFrac" ~ "Weighted UniFrac distance",
    nm == "UnweightedUniFrac" ~ "Unweighted UniFrac distance",
    TRUE ~ nm
  )

  log_msg("Running PCoA for", pretty_label)
  res <- make_pcoa_scores(distances[[nm]], meta_by_distance[[nm]], pretty_label)
  ordination_results[[nm]] <- res

  write_table(res$scores, paste0("07_pcoa_scores_", nm, ".csv"))
  write_table(res$eig, paste0("07_pcoa_eigenvalues_", nm, ".csv"))

  p <- plot_pcoa(
    scores = res$scores,
    pct1 = res$pct1,
    pct2 = res$pct2,
    title = paste0("PCoA - ", pretty_label),
    subtitle = paste0("Dotted ellipses: ", ellipse_group, " groups; descriptive ordination"),
    filename_base = paste0("Fig_ordination_", nm, "_pcoa_dotted_ellipses")
  )
  plots[[nm]] <- p
}

# Combined Bray + weighted UniFrac 2-panel figure.
if (all(c("BrayCurtis", "WeightedUniFrac") %in% names(plots))) {
  p_combined <- patchwork::wrap_plots(
    plots$BrayCurtis + ggplot2::theme(legend.position = "bottom"),
    plots$WeightedUniFrac + ggplot2::theme(legend.position = "bottom"),
    nrow = 1,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = "Host-associated microbiome ordination",
      subtitle = "PCoA based on Bray-Curtis and weighted UniFrac distances"
    ) &
    ggplot2::theme(legend.position = "bottom")

  save_plot_multi(
    p_combined,
    "Fig_ordination_BrayCurtis_WeightedUniFrac_2panel_dotted_ellipses",
    width = save_width_combined,
    height = save_height_combined
  )
}

# Optional 3-panel combined figure if unweighted UniFrac was also computed.
if (all(c("BrayCurtis", "WeightedUniFrac", "UnweightedUniFrac") %in% names(plots))) {
  p_combined3 <- patchwork::wrap_plots(
    plots$BrayCurtis + ggplot2::theme(legend.position = "bottom"),
    plots$WeightedUniFrac + ggplot2::theme(legend.position = "bottom"),
    plots$UnweightedUniFrac + ggplot2::theme(legend.position = "bottom"),
    nrow = 1,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = "Host-associated microbiome ordination",
      subtitle = "PCoA based on Bray-Curtis, weighted UniFrac, and unweighted UniFrac distances"
    ) &
    ggplot2::theme(legend.position = "bottom")

  save_plot_multi(
    p_combined3,
    "Fig_ordination_BrayCurtis_Weighted_Unweighted_UniFrac_3panel_dotted_ellipses",
    width = 18.0,
    height = save_height_combined
  )
}

############################
# 6. PERMANOVA AND DISPERSION DIAGNOSTICS
############################

log_msg("Running PERMANOVA summaries")
permanova_tbl <- purrr::imap_dfr(distances, function(d, nm) {
  run_permanova(d, meta_by_distance[[nm]], nm)
})
write_table(permanova_tbl, "07_permanova_BrayCurtis_UniFrac.csv")

log_msg("Running betadisper dispersion diagnostics")
bd_anova <- list()
bd_tukey <- list()
for (nm in names(distances)) {
  for (grp in c("Compartment", "Site")) {
    bd <- run_betadisper(distances[[nm]], meta_by_distance[[nm]], nm, grp)
    bd_anova[[paste(nm, grp, sep = "_")]] <- bd$anova
    bd_tukey[[paste(nm, grp, sep = "_")]] <- bd$tukey
  }
}
write_table(dplyr::bind_rows(bd_anova), "07_betadisper_anova_BrayCurtis_UniFrac.csv")
write_table(dplyr::bind_rows(bd_tukey), "07_betadisper_tukey_BrayCurtis_UniFrac.csv")

############################
# 7. CAPTION / INTERPRETATION NOTES
############################

notes_file <- file.path(tables_dir, "07_ordination_interpretation_notes.txt")
cat(
  "Bray-Curtis and UniFrac PCoA ordination notes\n",
  "==============================================\n\n",
  "This script generates distance-based PCoA ordinations for host-associated Haloclystus microbiomes.\n",
  "PCoA is used because Bray-Curtis and UniFrac are dissimilarity/distance matrices.\n\n",
  "Recommended figure language:\n",
  "Principal coordinates analysis based on Bray-Curtis dissimilarity and weighted UniFrac distance showed host-associated microbiome structure by body compartment and site. Dotted ellipses indicate compartment-level grouping and are descriptive only.\n\n",
  "Recommended methods language:\n",
  "Bray-Curtis dissimilarity was calculated from ASV count tables. Weighted and unweighted UniFrac distances were calculated using the host phylogeny in the phyloseq object. PCoA was performed using ape::pcoa. PERMANOVA was run using vegan::adonis2, with an unrestricted model for Site + Compartment and a repeated-measures-aware blocked model for Compartment using Individual_ID as permutation blocks. Homogeneity of multivariate dispersion was assessed using vegan::betadisper.\n\n",
  "Interpretation caution:\n",
  "Blocked PERMANOVA is most appropriate for testing Compartment because compartments can be repeated within individual animals. Site is not cleanly tested under Individual_ID blocks when Individual_ID is nested within Site; use the unrestricted Site term descriptively and report dispersion diagnostics.\n\n",
  "Key outputs:\n",
  "07_permanova_BrayCurtis_UniFrac.csv\n",
  "07_betadisper_anova_BrayCurtis_UniFrac.csv\n",
  "Fig_ordination_BrayCurtis_WeightedUniFrac_2panel_dotted_ellipses.*\n",
  "Fig_ordination_BrayCurtis_pcoa_dotted_ellipses.*\n",
  "Fig_ordination_WeightedUniFrac_pcoa_dotted_ellipses.*\n",
  file = notes_file
)

############################
# 8. SESSION INFO
############################

sink(file.path(ordination_dir, "session_info_ordination_bray_unifrac_pcoa.txt"))
print(sessionInfo())
sink()

log_msg("Standalone Bray-Curtis / UniFrac PCoA ordination complete.")
log_msg("Tables written to:", tables_dir)
log_msg("Plots written to:", plots_dir)
log_msg("Main figure candidate: Fig_ordination_BrayCurtis_WeightedUniFrac_2panel_dotted_ellipses.*")
