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

install_bioc_if_missing <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", type = "binary")
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    BiocManager::install(missing, ask = FALSE, update = FALSE, type = "binary")
  }
}

cran_pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "readr", "tibble", "purrr",
  "vegan", "ape", "phangorn", "patchwork", "RColorBrewer", "scales"
)

bioc_pkgs <- c(
  "phyloseq", "DECIPHER", "Biostrings"
)

optional_bioc_pkgs <- c("ANCOMBC", "microbiome")

install_cran_if_missing(cran_pkgs)
install_bioc_if_missing(bioc_pkgs)
install_bioc_if_missing(optional_bioc_pkgs)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(purrr)
  library(vegan)
  library(ape)
  library(phangorn)
  library(phyloseq)
  library(DECIPHER)
  library(Biostrings)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
})

run_ancombc <- requireNamespace("ANCOMBC", quietly = TRUE) &&
  requireNamespace("microbiome", quietly = TRUE)
############################
# 1. PATHS
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

if (!dir.exists(base_dir)) stop("Base directory does not exist: ", base_dir)
if (!dir.exists(output_dir)) stop("Output directory does not exist: ", output_dir)
if (!dir.exists(rds_dir)) stop("RDS directory does not exist: ", rds_dir)

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

cat("Using base_dir:\n", base_dir, "\n")
cat("Using output_dir:\n", output_dir, "\n")
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

save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  ggplot2::ggsave(
    filename = file.path(plots_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
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
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  t(apply(mat, 1, function(x) {
    z <- x + pseudocount
    logz <- log(z)
    logz - mean(logz)
  }))
}

safe_ellipse <- function(...) {
  tryCatch(
    ggplot2::stat_ellipse(...),
    error = function(e) NULL,
    warning = function(w) suppressWarnings(ggplot2::stat_ellipse(...))
  )
}

pairwise_permanova <- function(dist_obj, meta, group_var, adjust_by = NULL, permutations = 9999) {
  groups <- levels(factor(meta[[group_var]]))
  if (length(groups) < 2) return(tibble())
  
  combs <- combn(groups, 2, simplify = FALSE)
  
  purrr::map_dfr(combs, function(cc) {
    keep <- meta[[group_var]] %in% cc
    meta_sub <- droplevels(meta[keep, , drop = FALSE])
    
    if (nrow(meta_sub) < 3) {
      return(tibble(
        group1 = cc[1],
        group2 = cc[2],
        F = NA_real_,
        R2 = NA_real_,
        p = NA_real_
      ))
    }
    
    dist_mat <- as.matrix(dist_obj)
    dist_sub <- as.dist(dist_mat[keep, keep])
    
    formula_txt <- if (is.null(adjust_by) || !(adjust_by %in% colnames(meta_sub))) {
      paste0("dist_sub ~ ", group_var)
    } else {
      paste0("dist_sub ~ ", adjust_by, " + ", group_var)
    }
    
    ad <- tryCatch(
      vegan::adonis2(as.formula(formula_txt), data = meta_sub, permutations = permutations, by = "margin"),
      error = function(e) NULL
    )
    
    if (is.null(ad)) {
      return(tibble(
        group1 = cc[1],
        group2 = cc[2],
        F = NA_real_,
        R2 = NA_real_,
        p = NA_real_
      ))
    }
    
    term_row <- rownames(ad) == group_var
    
    tibble(
      group1 = cc[1],
      group2 = cc[2],
      F = unname(ad$F[term_row]),
      R2 = unname(ad$R2[term_row]),
      p = unname(ad$`Pr(>F)`[term_row])
    )
  }) |>
    mutate(p_adj = p.adjust(p, method = "BH"))
}

get_top_taxa_palette <- function(n) {
  if (n <= 12) {
    RColorBrewer::brewer.pal(max(3, n), "Paired")[seq_len(n)]
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n)
  }
}

build_tree_from_asvs <- function(seqtab_matrix, tree_rds_path) {
  if (file.exists(tree_rds_path)) return(readRDS(tree_rds_path))
  
  seqs <- Biostrings::DNAStringSet(colnames(seqtab_matrix))
  names(seqs) <- colnames(seqtab_matrix)
  
  alignment <- DECIPHER::AlignSeqs(seqs, anchor = NA, verbose = FALSE)
  phang_align <- phangorn::phyDat(as.matrix(alignment), type = "DNA")
  dm <- phangorn::dist.ml(phang_align)
  treeNJ <- phangorn::NJ(dm)
  fit <- phangorn::pml(treeNJ, data = phang_align)
  
  fitGTR <- phangorn::optim.pml(
    fit,
    model = "GTR",
    optInv = TRUE,
    optGamma = TRUE,
    rearrangement = "stochastic",
    control = phangorn::pml.control(trace = 0)
  )
  
  tree_rooted <- phangorn::midpoint(fitGTR$tree)
  saveRDS(tree_rooted, tree_rds_path)
  tree_rooted
}

extract_ancom_tables <- function(ancom_obj) {
  out <- list()
  
  for (nm in names(ancom_obj)) {
    if (is.data.frame(ancom_obj[[nm]]) || tibble::is_tibble(ancom_obj[[nm]])) {
      out[[nm]] <- as.data.frame(ancom_obj[[nm]])
    }
  }
  
  if ("res" %in% names(ancom_obj) && is.list(ancom_obj$res)) {
    for (nm in names(ancom_obj$res)) {
      if (is.data.frame(ancom_obj$res[[nm]]) || tibble::is_tibble(ancom_obj$res[[nm]])) {
        out[[paste0("res_", nm)]] <- as.data.frame(ancom_obj$res[[nm]])
      }
    }
  }
  
  out
}

############################
# 3. LOAD INPUT OBJECTS
############################

seqtab_file <- file.path(rds_dir, "seqtab_haloclystus_prevfiltered.rds")
taxa_file   <- file.path(rds_dir, "taxa_haloclystus_prevfiltered.rds")
meta_file   <- file.path(rds_dir, "metadata_haloclystus.rds")

if (!file.exists(seqtab_file)) stop("Missing file: ", seqtab_file)
if (!file.exists(taxa_file)) stop("Missing file: ", taxa_file)
if (!file.exists(meta_file)) stop("Missing file: ", meta_file)

seqtab_halo_prev <- readRDS(seqtab_file)
taxa_halo_prev <- readRDS(taxa_file)
meta_halo <- readRDS(meta_file)

meta_halo <- as.data.frame(meta_halo)

if (!("Compartment" %in% colnames(meta_halo))) stop("Compartment column not found in metadata")
if (!("Site" %in% colnames(meta_halo))) stop("Site column not found in metadata")

meta_halo$Compartment <- factor(meta_halo$Compartment)
meta_halo$Site <- factor(meta_halo$Site)

if (!all(rownames(meta_halo) %in% rownames(seqtab_halo_prev))) {
  stop("Metadata rownames do not match sequence table rownames")
}

seqtab_halo_prev <- seqtab_halo_prev[rownames(meta_halo), , drop = FALSE]

log_msg("Loaded seqtab:", paste(dim(seqtab_halo_prev), collapse = " x "))
log_msg("Loaded taxa:", paste(dim(taxa_halo_prev), collapse = " x "))
log_msg("Loaded metadata:", nrow(meta_halo), "samples")

############################
# 4. PHYLOSEQ OBJECT
############################

otu <- phyloseq::otu_table(seqtab_halo_prev, taxa_are_rows = FALSE)
tax <- phyloseq::tax_table(as.matrix(taxa_halo_prev))
sam <- phyloseq::sample_data(meta_halo)

ps_halo <- phyloseq::phyloseq(otu, tax, sam)
saveRDS(ps_halo, file.path(rds_dir, "phyloseq_haloclystus.rds"))

############################
# 5. ALPHA DIVERSITY
############################

alpha_df <- tibble::tibble(
  SampleID = rownames(seqtab_halo_prev),
  Compartment = meta_halo$Compartment,
  Site = meta_halo$Site,
  Shannon = vegan::diversity(seqtab_halo_prev, index = "shannon"),
  InvSimpson = vegan::diversity(seqtab_halo_prev, index = "invsimpson"),
  Observed = vegan::specnumber(seqtab_halo_prev),
  Fisher = vegan::fisher.alpha(seqtab_halo_prev)
)

write_table(alpha_df, "alpha_diversity_haloclystus.csv")

alpha_long <- alpha_df |>
  tidyr::pivot_longer(
    cols = c(Shannon, InvSimpson, Observed, Fisher),
    names_to = "Metric",
    values_to = "Value"
  )

p_alpha <- ggplot2::ggplot(alpha_long, ggplot2::aes(x = Compartment, y = Value, fill = Compartment)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
  ggplot2::facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  theme_fed() +
  ggplot2::labs(title = "Haloclystus alpha diversity by compartment", x = "Compartment", y = NULL)

save_plot(p_alpha, "alpha_diversity_haloclystus.png", width = 10, height = 8)

alpha_tests <- purrr::map_dfr(c("Shannon", "InvSimpson", "Observed", "Fisher"), function(metric) {
  df <- alpha_df |>
    dplyr::select(Compartment, Site, value = all_of(metric)) |>
    dplyr::filter(is.finite(value))
  
  if (nlevels(factor(df$Compartment)) < 2) return(tibble::tibble())
  
  kw <- kruskal.test(value ~ Compartment, data = df)
  
  tibble::tibble(
    metric = metric,
    test = "Kruskal-Wallis",
    statistic = unname(kw$statistic),
    p = kw$p.value
  )
})

write_table(alpha_tests, "alpha_diversity_tests.csv")

############################
# 6–8. BETA DIVERSITY
# Replace old Sections 6, 7, and 8 with this whole block
############################

# Helper used only for beta-diversity plots
add_compartment_ellipses <- function(df) {
  valid_groups <- df |>
    dplyr::count(Compartment) |>
    dplyr::filter(n >= 3) |>
    dplyr::pull(Compartment)
  
  if (length(valid_groups) == 0) return(NULL)
  
  ggplot2::stat_ellipse(
    data = df |> dplyr::filter(Compartment %in% valid_groups),
    mapping = ggplot2::aes(
      x = Axis1,
      y = Axis2,
      color = Compartment,
      group = Compartment
    ),
    inherit.aes = FALSE,
    type = "norm",
    linewidth = 0.5,
    linetype = 2,
    show.legend = FALSE
  )
}

############################
# 6. BETA DIVERSITY - AITCHISON
############################

clr_mat <- clr_transform(seqtab_halo_prev, pseudocount = 1)
dist_ait <- dist(clr_mat, method = "euclidean")
pcoa_ait <- ape::pcoa(dist_ait)

pcoa_ait_df <- tibble::tibble(
  Axis1 = pcoa_ait$vectors[, 1],
  Axis2 = pcoa_ait$vectors[, 2],
  Compartment = meta_halo$Compartment,
  Site = meta_halo$Site,
  SampleID = rownames(meta_halo)
)

p_ait <- ggplot2::ggplot(pcoa_ait_df, ggplot2::aes(x = Axis1, y = Axis2)) +
  ggplot2::geom_point(
    ggplot2::aes(color = Compartment, shape = Site),
    size = 3
  ) +
  add_compartment_ellipses(pcoa_ait_df) +
  theme_fed() +
  ggplot2::labs(
    title = "Haloclystus Aitchison PCoA",
    x = "Axis 1",
    y = "Axis 2"
  )

save_plot(p_ait, "haloclystus_pcoa_aitchison.png", width = 8, height = 6)

ait_formula <- if (nlevels(meta_halo$Site) > 1) {
  dist_ait ~ Site + Compartment
} else {
  dist_ait ~ Compartment
}

ait_adonis <- vegan::adonis2(
  ait_formula,
  data = meta_halo,
  permutations = 9999,
  by = "margin"
)

write_table(
  as.data.frame(ait_adonis) |> tibble::rownames_to_column("term"),
  "permanova_aitchison.csv"
)

ait_bd <- vegan::betadisper(dist_ait, meta_halo$Compartment)
ait_bd_perm <- permutest(ait_bd)

betadisper_ait_tbl <- tibble::tibble(
  test = "betadisper_aitchison",
  F = unname(ait_bd_perm$tab[1, "F"]),
  p = unname(ait_bd_perm$tab[1, "Pr(>F)"])
)

write_table(betadisper_ait_tbl, "betadisper_aitchison.csv")

pair_ait <- pairwise_permanova(
  dist_ait,
  meta_halo,
  group_var = "Compartment",
  adjust_by = if (nlevels(meta_halo$Site) > 1) "Site" else NULL
)

write_table(pair_ait, "pairwise_permanova_aitchison.csv")

############################
# 7. BETA DIVERSITY - BRAY-CURTIS
############################

bray <- vegan::vegdist(seqtab_halo_prev, method = "bray")
pcoa_bray <- ape::pcoa(bray)

pcoa_bray_df <- tibble::tibble(
  Axis1 = pcoa_bray$vectors[, 1],
  Axis2 = pcoa_bray$vectors[, 2],
  Compartment = meta_halo$Compartment,
  Site = meta_halo$Site,
  SampleID = rownames(meta_halo)
)

p_bray <- ggplot2::ggplot(pcoa_bray_df, ggplot2::aes(x = Axis1, y = Axis2)) +
  ggplot2::geom_point(
    ggplot2::aes(color = Compartment, shape = Site),
    size = 3
  ) +
  add_compartment_ellipses(pcoa_bray_df) +
  theme_fed() +
  ggplot2::labs(
    title = "Haloclystus Bray-Curtis PCoA",
    x = "Axis 1",
    y = "Axis 2"
  )

save_plot(p_bray, "haloclystus_pcoa_bray.png", width = 8, height = 6)

bray_formula <- if (nlevels(meta_halo$Site) > 1) {
  bray ~ Site + Compartment
} else {
  bray ~ Compartment
}

bray_adonis <- vegan::adonis2(
  bray_formula,
  data = meta_halo,
  permutations = 9999,
  by = "margin"
)

write_table(
  as.data.frame(bray_adonis) |> tibble::rownames_to_column("term"),
  "permanova_bray.csv"
)

bray_bd <- vegan::betadisper(bray, meta_halo$Compartment)
bray_bd_perm <- permutest(bray_bd)

betadisper_bray_tbl <- tibble::tibble(
  test = "betadisper_bray",
  F = unname(bray_bd_perm$tab[1, "F"]),
  p = unname(bray_bd_perm$tab[1, "Pr(>F)"])
)

write_table(betadisper_bray_tbl, "betadisper_bray.csv")

pair_bray <- pairwise_permanova(
  bray,
  meta_halo,
  group_var = "Compartment",
  adjust_by = if (nlevels(meta_halo$Site) > 1) "Site" else NULL
)

write_table(pair_bray, "pairwise_permanova_bray.csv")

############################
# 8. BETA DIVERSITY - WEIGHTED UNIFRAC
############################

tree_file <- file.path(rds_dir, "haloclystus_tree.rds")
phy_tree_obj <- build_tree_from_asvs(seqtab_halo_prev, tree_file)

ps_halo_tree <- phyloseq::phyloseq(
  phyloseq::otu_table(seqtab_halo_prev, taxa_are_rows = FALSE),
  phyloseq::tax_table(as.matrix(taxa_halo_prev)),
  phyloseq::sample_data(meta_halo),
  phyloseq::phy_tree(phy_tree_obj)
)

saveRDS(ps_halo_tree, file.path(rds_dir, "phyloseq_haloclystus_tree.rds"))

wu <- phyloseq::UniFrac(
  ps_halo_tree,
  weighted = TRUE,
  normalized = TRUE,
  parallel = FALSE,
  fast = TRUE
)

pcoa_wu <- ape::pcoa(wu)

pcoa_wu_df <- tibble::tibble(
  Axis1 = pcoa_wu$vectors[, 1],
  Axis2 = pcoa_wu$vectors[, 2],
  Compartment = meta_halo$Compartment,
  Site = meta_halo$Site,
  SampleID = rownames(meta_halo)
)

p_wu <- ggplot2::ggplot(pcoa_wu_df, ggplot2::aes(x = Axis1, y = Axis2)) +
  ggplot2::geom_point(
    ggplot2::aes(color = Compartment, shape = Site),
    size = 3
  ) +
  add_compartment_ellipses(pcoa_wu_df) +
  theme_fed() +
  ggplot2::labs(
    title = "Haloclystus weighted UniFrac PCoA",
    x = "Axis 1",
    y = "Axis 2"
  )

save_plot(p_wu, "haloclystus_pcoa_weighted_unifrac.png", width = 8, height = 6)

wu_formula <- if (nlevels(meta_halo$Site) > 1) {
  wu ~ Site + Compartment
} else {
  wu ~ Compartment
}

wu_adonis <- vegan::adonis2(
  wu_formula,
  data = meta_halo,
  permutations = 9999,
  by = "margin"
)

write_table(
  as.data.frame(wu_adonis) |> tibble::rownames_to_column("term"),
  "permanova_weighted_unifrac.csv"
)

wu_bd <- vegan::betadisper(wu, meta_halo$Compartment)
wu_bd_perm <- permutest(wu_bd)

betadisper_wu_tbl <- tibble::tibble(
  test = "betadisper_weighted_unifrac",
  F = unname(wu_bd_perm$tab[1, "F"]),
  p = unname(wu_bd_perm$tab[1, "Pr(>F)"])
)

write_table(betadisper_wu_tbl, "betadisper_weighted_unifrac.csv")

pair_wu <- pairwise_permanova(
  wu,
  meta_halo,
  group_var = "Compartment",
  adjust_by = if (nlevels(meta_halo$Site) > 1) "Site" else NULL
)

write_table(pair_wu, "pairwise_permanova_weighted_unifrac.csv")

############################
# 9. TAXONOMIC COMPOSITION
############################

rel_abund <- sweep(seqtab_halo_prev, 1, rowSums(seqtab_halo_prev), "/")

plot_tax_level <- function(rel_abund, taxa_mat, meta_df, level = "Family", top_n = 15, filename = NULL) {
  tax_df <- as.data.frame(taxa_mat, stringsAsFactors = FALSE) |>
    tibble::rownames_to_column("ASV")
  
  long_df <- as.data.frame(rel_abund) |>
    tibble::rownames_to_column("SampleID") |>
    tidyr::pivot_longer(-SampleID, names_to = "ASV", values_to = "RelAbund") |>
    dplyr::left_join(
      meta_df |> tibble::rownames_to_column("SampleID") |> dplyr::select(SampleID, Compartment),
      by = "SampleID"
    ) |>
    dplyr::left_join(tax_df, by = "ASV") |>
    dplyr::mutate(Taxon = .data[[level]], Taxon = ifelse(is.na(Taxon) | Taxon == "", "Unclassified", Taxon))
  
  mean_df <- long_df |>
    dplyr::group_by(Compartment, Taxon) |>
    dplyr::summarise(RelAbund = mean(RelAbund), .groups = "drop")
  
  top_taxa <- mean_df |>
    dplyr::group_by(Taxon) |>
    dplyr::summarise(total = sum(RelAbund), .groups = "drop") |>
    dplyr::arrange(desc(total)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(Taxon)
  
  plot_df <- mean_df |>
    dplyr::mutate(Taxon = ifelse(Taxon %in% top_taxa, Taxon, "Other")) |>
    dplyr::group_by(Compartment, Taxon) |>
    dplyr::summarise(RelAbund = sum(RelAbund), .groups = "drop")
  
  pal <- get_top_taxa_palette(length(unique(plot_df$Taxon)))
  names(pal) <- unique(plot_df$Taxon)
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Compartment, y = RelAbund, fill = Taxon)) +
    ggplot2::geom_col() +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_fill_manual(values = pal) +
    theme_fed() +
    ggplot2::labs(title = paste("Haloclystus composition -", level), x = "Compartment", y = "Mean relative abundance")
  
  if (!is.null(filename)) save_plot(p, filename, width = 10, height = 7)
  p
}

plot_tax_level(rel_abund, taxa_halo_prev, meta_halo, level = "Family", top_n = 15, filename = "composition_family_haloclystus.png")
plot_tax_level(rel_abund, taxa_halo_prev, meta_halo, level = "Genus", top_n = 15, filename = "composition_genus_haloclystus.png")

