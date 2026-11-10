# Haliclystus antarcticus 16S analysis.
# Reads host-only results created by 01_16S_DADA2_decontam_source.R.


############################
# 0. PACKAGES
############################

options(repos = c(CRAN = "https://cran.rstudio.com/"))

install_cran_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, dependencies = TRUE, type = "binary")
  }
}

cran_pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "readr", "tibble", "purrr", "stringr",
  "vegan", "ape", "phyloseq", "patchwork", "RColorBrewer", "scales",
  "lme4", "lmerTest", "emmeans", "UpSetR", "indicspecies",
  "ggrepel"
)

install_cran_if_missing(cran_pkgs)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(vegan)
  library(ape)
  library(phyloseq)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(UpSetR)
  library(indicspecies)
  library(ggrepel)
})

############################
# 1. REPOSITORY-RELATIVE PATHS
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
plots_dir  <- file.path(output_dir, "plots")
tables_dir <- file.path(output_dir, "tables")

dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

############################
# 2. HELPERS (same style as your existing scripts)
############################

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " ")))
  flush.console()
}

write_table <- function(obj, filename) {
  readr::write_csv(obj, file.path(tables_dir, filename))
}

save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  ggplot2::ggsave(
    filename = file.path(plots_dir, filename),
    plot = plot_obj, width = width, height = height, dpi = 300, bg = "white"
  )
}

theme_fed <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.key = ggplot2::element_blank()
    )
}

clr_transform <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat); storage.mode(mat) <- "numeric"
  t(apply(mat, 1, function(x) { z <- x + pseudocount; logz <- log(z); logz - mean(logz) }))
}

############################
# 3. LOAD YOUR EXISTING OBJECTS
############################

seqtab_halo_prev <- readRDS(file.path(rds_dir, "seqtab_haloclystus_prevfiltered.rds"))
taxa_halo_prev   <- readRDS(file.path(rds_dir, "taxa_haloclystus_prevfiltered.rds"))
meta_halo        <- readRDS(file.path(rds_dir, "metadata_haloclystus.rds")) |> as.data.frame()
ps_halo_tree     <- readRDS(file.path(rds_dir, "phyloseq_haloclystus_tree.rds"))  # from beta-div script

meta_halo$Compartment <- factor(meta_halo$Compartment)
meta_halo$Site        <- factor(meta_halo$Site)
seqtab_halo_prev <- seqtab_halo_prev[rownames(meta_halo), , drop = FALSE]

############################
# 4. DERIVE / VERIFY INDIVIDUAL (ANIMAL) ID  <-- CHECK THIS CAREFULLY
############################
# Repeated-measures structure: gut_body, stalk, tentacle almost certainly come
# from the SAME individual animal. If you have an explicit Individual_ID /
# Animal_ID column in your original metadata, use that instead of this guess.

if ("Individual_ID" %in% colnames(meta_halo)) {
  meta_halo$Individual_ID <- factor(meta_halo$Individual_ID)
} else {
  # IMPORTANT: your two sites do NOT use the same naming convention.
  # Baia_Byers is clean: <prefix><B|S|T><number>, e.g. CBB1 / CBS1 / CBT1.
  # Rizopatron is inconsistent: CR1S, CR1T, CR2, CR2B, CRB1, CRB1B, CRT2 ...
  # A single "letter-then-number" regex breaks on Rizopatron IDs like CR2B
  # (number-then-letter) or CR2 (no compartment letter at all).
  #
  # Robust fix: extract ONLY the digits from Sample_ID and combine with the
  # existing (already-curated) Site column. This correctly reconstructs the
  # individual regardless of where/whether a compartment letter appears.
  
  digits <- stringr::str_extract(meta_halo$Sample_ID, "\\d+")
  
  if (any(is.na(digits))) {
    bad <- meta_halo$Sample_ID[is.na(digits)]
    stop(
      "Sample_ID(s) contain no digits, so an individual number could not be extracted: ",
      paste(bad, collapse = ", "),
      "\nCheck these IDs manually."
    )
  }
  
  meta_halo$Individual_ID <- factor(paste0(meta_halo$Site, "_", digits))
}

# Rizopatron gut_body samples were each extracted twice from the same
# individual (confirmed genuine replicates, not data-entry duplicates). Both
# are kept as extra within-individual replication - this is handled correctly
# by the stratified PERMANOVA (strata = Individual_ID permutes within blocks
# regardless of block size) and by the LMM (Individual_ID as random effect
# absorbs the extra replication). No merging/averaging needed.
dup_check <- meta_halo |>
  dplyr::count(Individual_ID, Compartment) |>
  dplyr::filter(n > 1)

if (nrow(dup_check) > 0) {
  write_table(dup_check, "individual_compartment_replicates.csv")
  log_msg(
    nrow(dup_check), " Individual_ID x Compartment combination(s) have >1 sample ",
    "(confirmed Rizopatron gut_body re-extractions) - kept as extra replication. ",
    "See individual_compartment_replicates.csv for the list."
  )
}

log_msg("Individuals detected:", nlevels(meta_halo$Individual_ID))
individual_check <- meta_halo |>
  dplyr::count(Individual_ID, Compartment) |>
  tidyr::pivot_wider(names_from = Compartment, values_from = n, values_fill = 0)
write_table(individual_check, "CHECK_individual_id_assignment.csv")
# >>> Open this CSV: each row (individual) should have 0/1 in each compartment
# >>> column (0 = that compartment wasn't sampled for this animal - expected;
# >>> anything >1 means a duplicate ID or a parsing error worth investigating).

############################
# 5. DISTANCE MATRICES (recomputed to match your beta-div script exactly)
############################

clr_mat  <- clr_transform(seqtab_halo_prev, pseudocount = 1)
dist_ait <- dist(clr_mat, method = "euclidean")
dist_bray <- vegan::vegdist(seqtab_halo_prev, method = "bray")
dist_wu   <- phyloseq::UniFrac(ps_halo_tree, weighted = TRUE, normalized = TRUE, fast = TRUE)

dist_list <- list(Aitchison = dist_ait, BrayCurtis = dist_bray, WeightedUniFrac = dist_wu)

############################
# 6. STRATIFIED PERMANOVA (blocks permutations within Individual_ID)
############################
# strata = Individual_ID restricts permutations to occur only WITHIN each
# animal, which is the correct way to test Compartment when compartments are
# repeated-measures on the same individual. Site is still a free fixed effect.

run_stratified_permanova <- function(dist_obj, meta, label) {
  perm_ctrl <- how(nperm = 9999, blocks = meta$Individual_ID)
  
  formula_txt <- if (nlevels(meta$Site) > 1) {
    "dist_obj ~ Site + Compartment"
  } else {
    "dist_obj ~ Compartment"
  }
  
  res <- vegan::adonis2(as.formula(formula_txt), data = meta,
                        permutations = perm_ctrl, by = "margin")
  
  out <- as.data.frame(res) |> tibble::rownames_to_column("term") |>
    dplyr::mutate(distance = label)
  write_table(out, paste0("permanova_stratified_", label, ".csv"))
  out
}

stratified_permanova_results <- purrr::imap_dfr(dist_list, run_stratified_permanova, meta = meta_halo)
log_msg("Stratified PERMANOVA complete for all three distance metrics")

############################
# 7. ALPHA DIVERSITY - MIXED MODELS (Individual_ID as random effect)
############################

alpha_df <- tibble::tibble(
  SampleID = rownames(seqtab_halo_prev),
  Compartment = meta_halo$Compartment,
  Site = meta_halo$Site,
  Individual_ID = meta_halo$Individual_ID,
  Shannon = vegan::diversity(seqtab_halo_prev, index = "shannon"),
  InvSimpson = vegan::diversity(seqtab_halo_prev, index = "invsimpson"),
  Observed = vegan::specnumber(seqtab_halo_prev),
  Fisher = vegan::fisher.alpha(seqtab_halo_prev)
)
write_table(alpha_df, "alpha_diversity_haloclystus.csv")

fit_alpha_lmm <- function(metric_name, df) {
  d <- df |> dplyr::select(Compartment, Site, Individual_ID, value = all_of(metric_name)) |>
    dplyr::filter(is.finite(value))
  
  # if every individual has exactly one sample per compartment, the random
  # effect will be estimable; if not, this will warn (singular fit) - check.
  model <- tryCatch(
    lmerTest::lmer(value ~ Compartment + Site + (1 | Individual_ID), data = d),
    error = function(e) NULL, warning = function(w) suppressWarnings(
      lmerTest::lmer(value ~ Compartment + Site + (1 | Individual_ID), data = d))
  )
  if (is.null(model)) return(list(anova = tibble(), posthoc = tibble()))
  
  an <- as.data.frame(anova(model)) |> tibble::rownames_to_column("term") |>
    dplyr::mutate(metric = metric_name)
  
  emm <- emmeans::emmeans(model, pairwise ~ Compartment, adjust = "BH")
  posthoc <- as.data.frame(emm$contrasts) |> dplyr::mutate(metric = metric_name)
  
  list(anova = an, posthoc = posthoc)
}

alpha_metrics <- c("Shannon", "InvSimpson", "Observed", "Fisher")
alpha_lmm_results <- purrr::map(alpha_metrics, fit_alpha_lmm, df = alpha_df)
names(alpha_lmm_results) <- alpha_metrics

alpha_anova_tbl <- purrr::map_dfr(alpha_lmm_results, "anova")
alpha_posthoc_tbl <- purrr::map_dfr(alpha_lmm_results, "posthoc")

write_table(alpha_anova_tbl, "alpha_diversity_LMM_anova.csv")
write_table(alpha_posthoc_tbl, "alpha_diversity_LMM_posthoc.csv")

# figure with significance stars from the LMM posthoc, same layout as before
alpha_long <- alpha_df |>
  tidyr::pivot_longer(cols = all_of(alpha_metrics), names_to = "Metric", values_to = "Value")

p_alpha <- ggplot2::ggplot(alpha_long, ggplot2::aes(x = Compartment, y = Value, fill = Compartment)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
  ggplot2::facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  theme_fed() +
  ggplot2::labs(title = "Haloclystus alpha diversity by compartment (LMM-tested)",
                subtitle = "See alpha_diversity_LMM_posthoc.csv for pairwise significance",
                x = "Compartment", y = NULL)

save_plot(p_alpha, "alpha_diversity_haloclystus_LMM.png", width = 10, height = 8)

############################
# 8. RAREFACTION / EXTRAPOLATION CURVES (iNEXT, Hill number q=0)
############################

############################
# 8. RAREFACTION CURVES (vegan::rarecurve - fast, interpolation only)
############################
# NOTE: switched from iNEXT to vegan::rarecurve. iNEXT's per-sample bootstrap
# extrapolation can become pathologically slow (or appear to hang) when a
# handful of samples have unusually high read depth or a very skewed
# abundance distribution. rarecurve() gives the same core information
# (species accumulation vs. sequencing depth) without bootstrapping, and
# runs in seconds even for hundreds of samples - a standard, reviewer-
# accepted alternative for the "was sequencing depth sufficient" figure.

log_msg("Computing rarefaction curves (vegan::rarecurve)...")

raremax <- min(rowSums(seqtab_halo_prev))
log_msg("Rarefying to a step size based on the smallest library (", raremax, "reads )")

rare_list <- vegan::rarecurve(
  as.matrix(seqtab_halo_prev),
  step = max(1, round(raremax / 50)),  # ~50 points per curve
  tidy = TRUE
)

rare_df <- rare_list |>
  dplyr::rename(SampleID = Site, m = Sample, qD = Species) |>
  dplyr::left_join(
    tibble::tibble(SampleID = rownames(meta_halo), Compartment = meta_halo$Compartment,
                   Site = meta_halo$Site),
    by = "SampleID"
  )

p_rare <- ggplot2::ggplot(rare_df, ggplot2::aes(x = m, y = qD, group = SampleID, color = Compartment)) +
  ggplot2::geom_line(linewidth = 0.4, alpha = 0.7) +
  ggplot2::geom_vline(xintercept = raremax, linetype = 2, color = "grey40") +
  ggplot2::facet_wrap(~Site) +
  theme_fed() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)) +
  ggplot2::labs(title = "Rarefaction curves (interpolation only)",
                subtitle = paste0("Dashed line = smallest library size (", raremax, " reads)"),
                x = "Number of sequences", y = "Species richness")

save_plot(p_rare, "rarefaction_extrapolation_haloclystus.png", width = 10, height = 6)
write_table(rare_df, "rarefaction_extrapolation_data.csv")
log_msg("Rarefaction curves complete")

############################
# 9. CORE MICROBIOME / SHARED ASVs (UpSet plot)
############################
# ASV considered "present" in a compartment if found in >= presence_threshold
# fraction of that compartment's samples (default 33%: adjust as needed)

presence_threshold <- 1/3

core_sets <- purrr::map(levels(meta_halo$Compartment), function(cmp) {
  samp <- rownames(meta_halo)[meta_halo$Compartment == cmp]
  sub <- seqtab_halo_prev[samp, , drop = FALSE]
  prev <- colSums(sub > 0) / nrow(sub)
  names(prev)[prev >= presence_threshold]
})
names(core_sets) <- levels(meta_halo$Compartment)

upset_input <- UpSetR::fromList(core_sets)

png(file.path(plots_dir, "core_microbiome_upset.png"), width = 8, height = 5.5, units = "in", res = 300)
print(UpSetR::upset(upset_input, order.by = "freq", main.bar.color = "grey30",
                    sets.bar.color = c("#F8766D", "#00BA38", "#619CFF")[seq_along(core_sets)],
                    text.scale = 1.3))
dev.off()

core_summary <- tibble::tibble(
  Compartment = names(core_sets),
  n_ASVs_present = lengths(core_sets)
)
write_table(core_summary, "core_microbiome_summary.csv")

shared_all_three <- Reduce(intersect, core_sets)
write_table(tibble::tibble(ASV = shared_all_three), "core_ASVs_shared_all_compartments.csv")
log_msg("ASVs shared across all three compartments (>=", presence_threshold, "prevalence):",
        length(shared_all_three))

############################
# 10. INDICATOR SPECIES ANALYSIS (genus level, indicspecies::multipatt)
############################

genus_vec <- taxa_halo_prev[, "Genus"]
genus_vec[is.na(genus_vec) | genus_vec == ""] <- "Unclassified"

genus_counts <- rowsum(t(seqtab_halo_prev), group = genus_vec, reorder = FALSE)
genus_counts <- t(genus_counts)
genus_counts <- genus_counts[rownames(meta_halo), , drop = FALSE]

log_msg("Running indicator species analysis (multipatt, 999 permutations)...")

indval <- indicspecies::multipatt(
  genus_counts, meta_halo$Compartment,
  func = "IndVal.g", control = how(nperm = 999)
)

log_msg("Indicator species analysis complete")

indval_summary <- summary(indval, indvalcomp = TRUE)  # printed summary, also build a table below

indval_tbl <- tibble::as_tibble(indval$sign, rownames = "Genus") |>
  dplyr::filter(p.value <= 0.05) |>
  dplyr::arrange(p.value)

write_table(indval_tbl, "indicator_species_genus_multipatt.csv")

# bubble plot of top indicators
indval_plot_df <- indval_tbl |>
  dplyr::mutate(
    Compartment = levels(meta_halo$Compartment)[
      apply(dplyr::select(indval_tbl, dplyr::starts_with("s.")), 1, which.max)
    ]
  ) |>
  dplyr::group_by(Compartment) |>
  dplyr::slice_min(p.value, n = 10) |>
  dplyr::ungroup()

p_indval <- ggplot2::ggplot(indval_plot_df, ggplot2::aes(x = stat, y = reorder(Genus, stat), color = Compartment)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(title = "Indicator genera by compartment (IndVal, p <= 0.05)",
                x = "Indicator value (stat)", y = NULL) +
  theme_fed() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))

save_plot(p_indval, "indicator_species_genus.png", width = 8, height = 7)




############################
# 11. COMPOSITION BARPLOTS (Family and Genus level, mean relative abundance)
############################
# Same logic as your original composition_family/genus figures, rebuilt here
# so this script is self-contained. Adjust n_top_* to show more/fewer taxa
# before collapsing the rest into "Other".

plot_composition <- function(rank_col, n_top = 13, facet_site = FALSE) {
  tax_vec <- taxa_halo_prev[, rank_col]
  tax_vec[is.na(tax_vec) | tax_vec == ""] <- "Unclassified"
  
  rel <- sweep(seqtab_halo_prev, 1, rowSums(seqtab_halo_prev), "/")
  
  long_df <- as.data.frame(rel) |>
    tibble::rownames_to_column("SampleID") |>
    tidyr::pivot_longer(-SampleID, names_to = "ASV", values_to = "RelAbund") |>
    dplyr::mutate(Taxon = tax_vec[ASV]) |>
    dplyr::left_join(
      tibble::tibble(SampleID = rownames(meta_halo), Compartment = meta_halo$Compartment,
                     Site = meta_halo$Site),
      by = "SampleID"
    )
  
  mean_by_taxon <- long_df |>
    dplyr::group_by(Taxon) |>
    dplyr::summarise(MeanAbund = mean(RelAbund), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(MeanAbund))
  
  top_taxa <- mean_by_taxon |> dplyr::slice_head(n = n_top) |> dplyr::pull(Taxon)
  
  group_vars <- if (facet_site) c("Compartment", "Site", "Taxon") else c("Compartment", "Taxon")
  
  plot_df <- long_df |>
    dplyr::mutate(Taxon = ifelse(Taxon %in% top_taxa, Taxon, "Other")) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_vars, "SampleID")))) |>
    dplyr::summarise(RelAbund = sum(RelAbund), .groups = "drop") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(MeanRelAbund = mean(RelAbund), .groups = "drop")
  
  taxon_levels <- c(setdiff(top_taxa, "Other"), "Other")
  plot_df$Taxon <- factor(plot_df$Taxon, levels = taxon_levels)
  
  n_colors <- length(taxon_levels)
  palette_taxa <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(12, n_colors), "Paired"))(n_colors)
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Compartment, y = MeanRelAbund, fill = Taxon)) +
    ggplot2::geom_col() +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_fill_manual(values = palette_taxa) +
    theme_fed() +
    ggplot2::labs(title = paste0("Haloclystus composition - ", rank_col),
                  x = "Compartment", y = "Mean relative abundance", fill = "Taxon")
  
  if (facet_site) p <- p + ggplot2::facet_wrap(~Site)
  
  list(plot = p, data = plot_df)
}

log_msg("Building composition barplots (Family, Genus)...")

comp_family <- plot_composition("Family", n_top = 13, facet_site = FALSE)
save_plot(comp_family$plot, "composition_family_haloclystus_v2.png", width = 8, height = 6.5)
write_table(comp_family$data, "composition_family_data.csv")

comp_genus <- plot_composition("Genus", n_top = 15, facet_site = FALSE)
save_plot(comp_genus$plot, "composition_genus_haloclystus_v2.png", width = 8, height = 6.5)
write_table(comp_genus$data, "composition_genus_data.csv")

# same two, but split by Site as well - useful supplementary figure to check
# whether compartment composition is consistent across your two locations
comp_family_site <- plot_composition("Family", n_top = 13, facet_site = TRUE)
save_plot(comp_family_site$plot, "composition_family_by_site.png", width = 11, height = 6.5)

comp_genus_site <- plot_composition("Genus", n_top = 15, facet_site = TRUE)
save_plot(comp_genus_site$plot, "composition_genus_by_site.png", width = 11, height = 6.5)

log_msg("Composition barplots complete")

############################
# 12. PCA BIPLOT (CLR-transformed genus table, with taxon loading arrows)
############################


log_msg("Running PCA on CLR-transformed genus table...")

genus_clr <- clr_transform(genus_counts, pseudocount = 1)  # genus_counts built in Section 10
pca_res <- prcomp(genus_clr, center = TRUE, scale. = FALSE)

pct_var <- round(100 * (pca_res$sdev^2) / sum(pca_res$sdev^2), 1)

pca_scores <- as.data.frame(pca_res$x[, 1:2]) |>
  tibble::rownames_to_column("SampleID") |>
  dplyr::left_join(
    tibble::tibble(SampleID = rownames(meta_halo), Compartment = meta_halo$Compartment,
                   Site = meta_halo$Site),
    by = "SampleID"
  )

# top contributing genera by loading vector length (PC1 & PC2 combined)
loadings <- as.data.frame(pca_res$rotation[, 1:2]) |>
  tibble::rownames_to_column("Genus") |>
  dplyr::mutate(len = sqrt(PC1^2 + PC2^2)) |>
  dplyr::arrange(dplyr::desc(len)) |>
  dplyr::slice_head(n = 12)

# scale arrows to roughly fill the sample-score plotting space
arrow_scale <- 0.8 * max(abs(pca_scores$PC1), abs(pca_scores$PC2)) / max(loadings$len)
loadings <- loadings |> dplyr::mutate(PC1_scaled = PC1 * arrow_scale, PC2_scaled = PC2 * arrow_scale)

p_pca <- ggplot2::ggplot(pca_scores, ggplot2::aes(x = PC1, y = PC2)) +
  ggplot2::geom_point(ggplot2::aes(color = Compartment, shape = Site), size = 3, alpha = 0.85) +
  ggplot2::geom_segment(
    data = loadings,
    ggplot2::aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
    arrow = grid::arrow(length = grid::unit(0.2, "cm")), color = "grey30", inherit.aes = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = loadings, ggplot2::aes(x = PC1_scaled, y = PC2_scaled, label = Genus),
    size = 3, color = "grey20", inherit.aes = FALSE, max.overlaps = 20
  ) +
  theme_fed() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)) +
  ggplot2::labs(title = "Haloclystus PCA (CLR-transformed genus abundances)",
                x = paste0("PC1 (", pct_var[1], "%)"),
                y = paste0("PC2 (", pct_var[2], "%)"))

save_plot(p_pca, "pca_biplot_genus_clr.png", width = 9, height = 7.5)
write_table(pca_scores, "pca_scores_genus_clr.csv")
write_table(loadings |> dplyr::select(Genus, PC1, PC2, len), "pca_loadings_genus_clr.csv")

log_msg("PCA biplot complete")

############################
# 13. HEATMAP (top genera x samples, annotated by Compartment/Site)
############################

install_cran_if_missing("pheatmap")
suppressPackageStartupMessages(library(pheatmap))

log_msg("Building heatmap of top genera...")

n_top_heatmap <- 30
genus_rel <- sweep(genus_counts, 1, rowSums(genus_counts), "/")
top_genera_heatmap <- names(sort(colMeans(genus_rel), decreasing = TRUE))[1:min(n_top_heatmap, ncol(genus_rel))]

heatmap_mat <- t(genus_rel[, top_genera_heatmap, drop = FALSE])
heatmap_mat_log <- log10(heatmap_mat + 1e-4)  # log scale, small pseudocount for zeros

annotation_col <- data.frame(
  Compartment = meta_halo$Compartment,
  Site = meta_halo$Site,
  row.names = rownames(meta_halo)
)

ann_colors <- list(
  Compartment = c(gut_body = "#F8766D", stalk = "#00BA38", tentacle = "#619CFF"),
  Site = c(Baia_Byers = "#8C8C8C", Rizopatron = "#2C2C2C")
)

pheatmap::pheatmap(
  heatmap_mat_log,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  clustering_method = "ward.D2",
  show_colnames = FALSE,
  fontsize_row = 8,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))(100),
  main = paste0("Top ", length(top_genera_heatmap), " genera (log10 relative abundance)"),
  filename = file.path(plots_dir, "heatmap_top_genera.png"),
  width = 11, height = 8
)

write_table(
  as.data.frame(heatmap_mat) |> tibble::rownames_to_column("Genus"),
  "heatmap_top_genera_data.csv"
)

log_msg("Heatmap complete")

############################
# 14. DESeq2 DIFFERENTIAL ABUNDANCE 
############################


install_cran_if_missing("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages(library(DESeq2))

log_msg("Running DESeq2 differential abundance (ASV level)...")

dds <- phyloseq::phyloseq_to_deseq2(ps_halo_tree, ~ Site + Compartment)

# geometric-mean size factors robust to zeros (standard fix for microbiome
# count tables, since DESeq2's default median-ratio method fails when most
# taxa are zero in most samples)
geo_means <- apply(DESeq2::counts(dds), 1, function(x) {
  if (all(x == 0)) return(0)
  exp(mean(log(x[x > 0])))
})
dds <- DESeq2::estimateSizeFactors(dds, geoMeans = geo_means)
dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric")

compartment_pairs <- utils::combn(levels(meta_halo$Compartment), 2, simplify = FALSE)

deseq2_results <- purrr::map_dfr(compartment_pairs, function(pair) {
  res <- DESeq2::results(dds, contrast = c("Compartment", pair[1], pair[2]), alpha = 0.05)
  as.data.frame(res) |>
    tibble::rownames_to_column("ASV") |>
    dplyr::mutate(
      contrast = paste0(pair[1], "_vs_", pair[2]),
      Genus = taxa_halo_prev[ASV, "Genus"],
      Family = taxa_halo_prev[ASV, "Family"]
    )
})

write_table(deseq2_results, "deseq2_results_all_contrasts.csv")

deseq2_sig <- deseq2_results |>
  dplyr::filter(!is.na(padj), padj <= 0.05) |>
  dplyr::mutate(TaxonLabel = ifelse(is.na(Genus) | Genus == "", paste0("ASV_", substr(ASV, 1, 8)), Genus))

write_table(deseq2_sig, "deseq2_results_significant.csv")
log_msg("DESeq2 complete -", nrow(deseq2_sig), "significant ASVs (padj <= 0.05) across all contrasts")

# volcano plot, one panel per pairwise contrast
p_volcano <- ggplot2::ggplot(deseq2_results, ggplot2::aes(x = log2FoldChange, y = -log10(pvalue))) +
  ggplot2::geom_point(ggplot2::aes(color = !is.na(padj) & padj <= 0.05), alpha = 0.6, size = 1.2) +
  ggplot2::scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey70"),
                              labels = c("TRUE" = "padj <= 0.05", "FALSE" = "n.s."), name = NULL) +
  ggplot2::facet_wrap(~contrast) +
  theme_fed() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)) +
  ggplot2::labs(title = "DESeq2 differential abundance (ASV level)",
                x = "log2 fold-change", y = "-log10(p-value)")

save_plot(p_volcano, "deseq2_volcano.png", width = 11, height = 5)

# effect-size dot plot of top significant genera per contrast (same visual
# language as the ANCOM-BC2 plot in Section 12, for easy side-by-side comparison)
deseq2_plot_df <- deseq2_sig |>
  dplyr::group_by(contrast) |>
  dplyr::slice_max(abs(log2FoldChange), n = 15) |>
  dplyr::ungroup()

if (nrow(deseq2_plot_df) > 0) {
  p_deseq2 <- ggplot2::ggplot(deseq2_plot_df, ggplot2::aes(x = log2FoldChange, y = reorder(TaxonLabel, log2FoldChange), color = log2FoldChange > 0)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::facet_wrap(~contrast, scales = "free_y") +
    ggplot2::scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"), guide = "none") +
    theme_fed() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)) +
    ggplot2::labs(title = "DESeq2 significant ASVs (padj <= 0.05)", x = "log2 fold-change", y = NULL)
  
  save_plot(p_deseq2, "deseq2_effect_size.png", width = 10, height = 7)
}

# cross-check: how many taxa are called significant by BOTH ANCOM-BC2 and DESeq2
if (exists("ancom_plot_df") && nrow(ancom_plot_df) > 0) {
  log_msg("Compare deseq2_results_significant.csv against ancombc2 output to report ",
          "cross-method agreement in your methods/results section")
}

log_msg("DESeq2 differential abundance complete")

