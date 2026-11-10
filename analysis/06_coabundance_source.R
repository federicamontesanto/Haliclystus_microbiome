# Haliclystus antarcticus 16S analysis.

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

network_dir <- file.path(output_dir, "host_coabundance_network")
network_tables_dir <- file.path(network_dir, "tables")
network_plots_dir <- file.path(network_dir, "plots")
for (d in c(network_dir, network_tables_dir, network_plots_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Analysis level. Genus is recommended. Missing genus labels are represented as
# Family_unclassified so informative family-level unresolved taxa are not all merged.
taxon_rank <- "Genus"

# Sample exclusions matched to the rest of your host pipeline.
poor_quality_samples <- c("CR2B", "CR3B", "CR4B", "CR5B")

# Genus filtering before correlations. These conservative filters reduce sparse,
# unstable correlations.
min_prevalence_prop <- 0.20        # taxon must occur in at least 20% of host samples
min_mean_rel_abund <- 0.001        # taxon mean relative abundance must be >= 0.1%
max_taxa_for_network <- 50         # keep top taxa by mean relative abundance after filters

# CLR and correlation settings.
clr_pseudocount <- 1e-6            # applied to relative abundances before CLR
cor_method <- "spearman"           # rank-based, robust to non-normality
edge_abs_rho_threshold <- 0.55      # retain strong associations
edge_q_threshold <- 0.05            # FDR threshold across all taxon pairs

# If strict thresholds produce very few edges, this optional relaxed table/plot is useful
# for exploration. It is clearly labelled as relaxed.
make_relaxed_network_if_needed <- TRUE
min_edges_for_strict_network <- 8
relaxed_abs_rho_threshold <- 0.45
relaxed_q_threshold <- 0.10

# Compartment association of nodes.
# A genus is labelled gut/stalk/tentacle-associated if its highest mean abundance
# is at least this fold above the second-highest compartment mean.
compartment_ratio_threshold <- 1.5

# Plotting controls.
label_top_n_network_nodes <- 30
network_layout <- "fr"             # ggraph layout; "fr" is Fruchterman-Reingold

############################
# 1. PACKAGES
############################

options(repos = c(CRAN = "https://cran.rstudio.com/"))
options(timeout = max(3600, getOption("timeout")))

install_cran_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) install.packages(missing, dependencies = TRUE, type = "binary")
}

cran_pkgs <- c(
  "dplyr", "tidyr", "purrr", "stringr", "readr", "tibble", "ggplot2",
  "scales", "svglite", "ggrepel", "igraph", "tidygraph", "ggraph",
  "RColorBrewer"
)
install_cran_if_missing(cran_pkgs)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(svglite)
  library(ggrepel)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(RColorBrewer)
})

############################
# 2. HELPERS
############################

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " ")))
  flush.console()
}

write_table <- function(obj, filename) {
  readr::write_csv(obj, file.path(network_tables_dir, filename))
}

save_plot_multi <- function(plot_obj, filename_base, width = 8, height = 6) {
  ggplot2::ggsave(file.path(network_plots_dir, paste0(filename_base, ".png")), plot_obj,
                  width = width, height = height, dpi = 450, bg = "white")
  ggplot2::ggsave(file.path(network_plots_dir, paste0(filename_base, ".pdf")), plot_obj,
                  width = width, height = height, bg = "white", device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(network_plots_dir, paste0(filename_base, ".svg")), plot_obj,
                  width = width, height = height, bg = "white", device = svglite::svglite)
}

clean_taxon <- function(x, fallback = "Unclassified") {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- fallback
  x
}

theme_net <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "grey20"),
      axis.ticks = ggplot2::element_line(linewidth = 0.30, colour = "grey20"),
      axis.text = ggplot2::element_text(colour = "grey15"),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey94", colour = "grey55", linewidth = 0.35),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

clr_transform <- function(mat, pseudocount = 1e-6) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat <- mat + pseudocount
  log_mat <- log(mat)
  sweep(log_mat, 1, rowMeans(log_mat), "-")
}

safe_cor_test <- function(x, y, method = "spearman") {
  out <- tryCatch(
    suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(out)) return(c(rho = NA_real_, p = NA_real_))
  c(rho = unname(out$estimate), p = out$p.value)
}

make_taxon_label <- function(tax_df, rank_col = "Genus") {
  genus <- clean_taxon(tax_df$Genus)
  family <- clean_taxon(tax_df$Family)
  order <- clean_taxon(tax_df$Order)
  phylum <- clean_taxon(tax_df$Phylum)

  label <- dplyr::case_when(
    rank_col == "Genus" & !genus %in% c("Unclassified", "uncultured", "unknown") ~ genus,
    rank_col == "Genus" & !family %in% c("Unclassified", "uncultured", "unknown") ~ paste0(family, "_unclassified"),
    rank_col == "Genus" & !order  %in% c("Unclassified", "uncultured", "unknown") ~ paste0(order, "_unclassified"),
    rank_col == "Genus" & !phylum %in% c("Unclassified", "uncultured", "unknown") ~ paste0(phylum, "_unclassified"),
    rank_col == "Family" & !family %in% c("Unclassified", "uncultured", "unknown") ~ family,
    rank_col == "Family" & !order  %in% c("Unclassified", "uncultured", "unknown") ~ paste0(order, "_unclassified"),
    TRUE ~ "Unclassified"
  )

  label <- stringr::str_replace_all(label, "[^A-Za-z0-9_. -]+", "_")
  label <- stringr::str_squish(label)
  label
}

compartment_pal <- c(
  gut = "#A51E2D",
  stalk = "#F0B7A8",
  tentacle = "#177E89",
  shared = "#7A7A7A",
  unclear = "#C9C9C9"
)

edge_pal <- c(positive = "#7A1F2B", negative = "#2E6F73")

############################
# 3. LOAD HOST DATA
############################

log_msg("Loading host RDS objects from", rds_dir)
seqtab_halo_prev <- readRDS(file.path(rds_dir, "seqtab_haloclystus_prevfiltered.rds"))
taxa_halo_prev   <- readRDS(file.path(rds_dir, "taxa_haloclystus_prevfiltered.rds"))
meta_halo        <- readRDS(file.path(rds_dir, "metadata_haloclystus.rds")) |> as.data.frame()

# Defensive sample exclusions, matching the rest of your downstream workflow.
# The upstream input is restricted to the 65 host samples.
samples_to_exclude <- unique(poor_quality_samples)
found_to_exclude <- intersect(samples_to_exclude, meta_halo$Sample_ID)
log_msg("Excluding", length(found_to_exclude), "sample(s):", paste(found_to_exclude, collapse = ", "))

keep_mask <- !(meta_halo$Sample_ID %in% samples_to_exclude)
meta_halo <- meta_halo[keep_mask, , drop = FALSE]
seqtab_halo_prev <- seqtab_halo_prev[rownames(meta_halo), , drop = FALSE]

meta_halo$Compartment <- as.character(meta_halo$Compartment)
meta_halo$Compartment <- dplyr::recode(meta_halo$Compartment, gut_body = "gut", .default = meta_halo$Compartment)
meta_halo$Compartment <- factor(meta_halo$Compartment, levels = c("gut", "stalk", "tentacle"))
meta_halo$Site <- factor(meta_halo$Site)

# Ensure taxa rows follow ASV columns.
taxa_halo_prev <- as.data.frame(taxa_halo_prev, stringsAsFactors = FALSE)
taxa_halo_prev <- taxa_halo_prev[colnames(seqtab_halo_prev), , drop = FALSE]

############################
# 4. BUILD GENUS-LEVEL HOST TABLE
############################

host_tax <- taxa_halo_prev |>
  tibble::rownames_to_column("ASV") |>
  dplyr::mutate(
    Kingdom = clean_taxon(Kingdom),
    Phylum = clean_taxon(Phylum),
    Class = clean_taxon(Class),
    Order = clean_taxon(Order),
    Family = clean_taxon(Family),
    Genus = clean_taxon(Genus),
    Taxon = make_taxon_label(dplyr::cur_data_all(), rank_col = taxon_rank)
  ) |>
  dplyr::filter(
    !stringr::str_detect(stringr::str_to_lower(Taxon), "chloroplast|mitochondria"),
    !stringr::str_detect(stringr::str_to_lower(Family), "chloroplast|mitochondria"),
    !stringr::str_detect(stringr::str_to_lower(Order), "chloroplast|mitochondria")
  )

# Keep only ASVs with usable taxonomic labels.
seqtab_use <- seqtab_halo_prev[, host_tax$ASV, drop = FALSE]

asv_to_taxon <- host_tax |>
  dplyr::select(ASV, Taxon, Family, Genus, Phylum, Class, Order)

genus_counts_long <- as.data.frame(seqtab_use) |>
  tibble::rownames_to_column("SampleID") |>
  tidyr::pivot_longer(-SampleID, names_to = "ASV", values_to = "Count") |>
  dplyr::filter(Count > 0) |>
  dplyr::left_join(asv_to_taxon, by = "ASV") |>
  dplyr::group_by(SampleID, Taxon) |>
  dplyr::summarise(Count = sum(Count), .groups = "drop")

genus_counts_wide <- genus_counts_long |>
  tidyr::pivot_wider(names_from = Taxon, values_from = Count, values_fill = 0) |>
  as.data.frame()

rownames(genus_counts_wide) <- genus_counts_wide$SampleID
genus_counts_wide$SampleID <- NULL
genus_counts_mat <- as.matrix(genus_counts_wide)
storage.mode(genus_counts_mat) <- "numeric"

genus_rel_mat <- sweep(genus_counts_mat, 1, rowSums(genus_counts_mat), "/")
genus_rel_mat[is.na(genus_rel_mat)] <- 0

sample_info <- tibble::tibble(
  SampleID = rownames(meta_halo),
  Site = meta_halo$Site,
  Compartment = meta_halo$Compartment
)

# Taxon summary before filtering.
taxon_summary <- tibble::tibble(
  Taxon = colnames(genus_counts_mat),
  prevalence_n = colSums(genus_counts_mat > 0),
  prevalence_prop = colMeans(genus_counts_mat > 0),
  mean_rel_abund = colMeans(genus_rel_mat),
  max_rel_abund = apply(genus_rel_mat, 2, max)
) |>
  dplyr::left_join(
    asv_to_taxon |>
      dplyr::group_by(Taxon) |>
      dplyr::summarise(
        Family = dplyr::first(Family[!is.na(Family) & Family != "Unclassified"]),
        Genus = dplyr::first(Genus[!is.na(Genus) & Genus != "Unclassified"]),
        Phylum = dplyr::first(Phylum[!is.na(Phylum) & Phylum != "Unclassified"]),
        n_asvs_collapsed = dplyr::n_distinct(ASV),
        .groups = "drop"
      ),
    by = "Taxon"
  ) |>
  dplyr::mutate(
    Family = clean_taxon(Family),
    Genus = clean_taxon(Genus),
    Phylum = clean_taxon(Phylum)
  ) |>
  dplyr::arrange(dplyr::desc(mean_rel_abund))

write_table(taxon_summary, "05_coabundance_taxon_summary_before_filtering.csv")

selected_taxa <- taxon_summary |>
  dplyr::filter(prevalence_prop >= min_prevalence_prop, mean_rel_abund >= min_mean_rel_abund) |>
  dplyr::slice_head(n = max_taxa_for_network) |>
  dplyr::pull(Taxon)

if (length(selected_taxa) < 4) {
  stop("Fewer than 4 taxa passed the filters. Lower min_prevalence_prop or min_mean_rel_abund.")
}

log_msg("Taxa retained for co-abundance network:", length(selected_taxa))
write_table(
  taxon_summary |> dplyr::filter(Taxon %in% selected_taxa),
  "05_coabundance_taxa_used_for_network.csv"
)

mat_rel <- genus_rel_mat[, selected_taxa, drop = FALSE]
clr_mat <- clr_transform(mat_rel, pseudocount = clr_pseudocount)

write_table(
  as.data.frame(clr_mat) |> tibble::rownames_to_column("SampleID"),
  "05_coabundance_clr_matrix_used_for_correlations.csv"
)

############################
# 5. CORRELATION MATRIX AND EDGE TABLE
############################

log_msg("Computing", cor_method, "correlations on CLR-transformed", taxon_rank, "profiles")
cor_mat <- stats::cor(clr_mat, method = cor_method, use = "pairwise.complete.obs")

taxon_pairs <- utils::combn(colnames(clr_mat), 2, simplify = FALSE)
cor_pairs <- purrr::map_dfr(taxon_pairs, function(pair) {
  res <- safe_cor_test(clr_mat[, pair[1]], clr_mat[, pair[2]], method = cor_method)
  tibble::tibble(
    taxon1 = pair[1],
    taxon2 = pair[2],
    rho = as.numeric(res["rho"]),
    p_value = as.numeric(res["p"])
  )
}) |>
  dplyr::mutate(
    q_value = p.adjust(p_value, method = "BH"),
    association = dplyr::case_when(
      rho > 0 ~ "positive",
      rho < 0 ~ "negative",
      TRUE ~ "zero"
    ),
    abs_rho = abs(rho)
  ) |>
  dplyr::arrange(q_value, dplyr::desc(abs_rho))

write_table(cor_pairs, "05_coabundance_correlations_all_pairs.csv")

strict_edges <- cor_pairs |>
  dplyr::filter(!is.na(q_value), q_value <= edge_q_threshold, abs_rho >= edge_abs_rho_threshold) |>
  dplyr::arrange(dplyr::desc(abs_rho))

write_table(strict_edges, "05_coabundance_edges_strict.csv")

relaxed_edges <- cor_pairs |>
  dplyr::filter(!is.na(q_value), q_value <= relaxed_q_threshold, abs_rho >= relaxed_abs_rho_threshold) |>
  dplyr::arrange(dplyr::desc(abs_rho))
write_table(relaxed_edges, "05_coabundance_edges_relaxed_exploratory.csv")

edges_to_use <- strict_edges
threshold_label <- paste0("strict: |rho| >= ", edge_abs_rho_threshold, ", FDR q <= ", edge_q_threshold)
network_suffix <- "strict"

if (make_relaxed_network_if_needed && nrow(strict_edges) < min_edges_for_strict_network && nrow(relaxed_edges) >= min_edges_for_strict_network) {
  log_msg("Strict network has", nrow(strict_edges), "edges; using relaxed exploratory network for plotting.")
  edges_to_use <- relaxed_edges
  threshold_label <- paste0("relaxed exploratory: |rho| >= ", relaxed_abs_rho_threshold, ", FDR q <= ", relaxed_q_threshold)
  network_suffix <- "relaxed_exploratory"
} else {
  log_msg("Strict network edges:", nrow(strict_edges))
}

############################
# 6. NODE ANNOTATION
############################

# Mean abundance by compartment for node annotation.
rel_long <- as.data.frame(mat_rel) |>
  tibble::rownames_to_column("SampleID") |>
  tidyr::pivot_longer(-SampleID, names_to = "Taxon", values_to = "RelAbund") |>
  dplyr::left_join(sample_info, by = "SampleID")

compartment_means <- rel_long |>
  dplyr::group_by(Taxon, Compartment) |>
  dplyr::summarise(mean_rel_abund = mean(RelAbund), .groups = "drop") |>
  tidyr::complete(Taxon, Compartment = levels(meta_halo$Compartment), fill = list(mean_rel_abund = 0))

compartment_assoc <- compartment_means |>
  dplyr::group_by(Taxon) |>
  dplyr::arrange(dplyr::desc(mean_rel_abund), .by_group = TRUE) |>
  dplyr::summarise(
    dominant_compartment = as.character(Compartment[1]),
    dominant_mean_rel_abund = mean_rel_abund[1],
    second_mean_rel_abund = dplyr::if_else(dplyr::n() >= 2, mean_rel_abund[2], 0),
    dominance_ratio = dominant_mean_rel_abund / (second_mean_rel_abund + 1e-12),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    compartment_class = dplyr::case_when(
      dominant_mean_rel_abund <= 0 ~ "unclear",
      dominance_ratio >= compartment_ratio_threshold ~ dominant_compartment,
      TRUE ~ "shared"
    )
  )

# Kruskal-Wallis compartment association as a descriptive node-level test.
kw_tbl <- rel_long |>
  dplyr::group_by(Taxon) |>
  dplyr::summarise(
    kw_p = tryCatch(stats::kruskal.test(RelAbund ~ Compartment)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) |>
  dplyr::mutate(kw_q = p.adjust(kw_p, method = "BH"))

node_base <- taxon_summary |>
  dplyr::filter(Taxon %in% selected_taxa) |>
  dplyr::left_join(compartment_assoc, by = "Taxon") |>
  dplyr::left_join(kw_tbl, by = "Taxon") |>
  dplyr::mutate(
    compartment_class = factor(compartment_class, levels = c("gut", "stalk", "tentacle", "shared", "unclear")),
    label_priority = dplyr::dense_rank(dplyr::desc(mean_rel_abund))
  )

############################
# 7. BUILD NETWORK AND MODULE STATISTICS - igraph VERSION
############################

if (nrow(edges_to_use) == 0) {
  log_msg("No edges passed selected thresholds. Network plot will be skipped; heatmap still produced.")
  nodes_final <- node_base |>
    dplyr::mutate(
      module = NA_integer_,
      degree = 0L,
      strength = 0,
      positive_degree = 0L,
      negative_degree = 0L,
      Zi = 0,
      Pi = 0,
      role = "isolated/no selected edge"
    )
  write_table(nodes_final, paste0("05_coabundance_nodes_", network_suffix, ".csv"))
  module_summary <- tibble::tibble()
} else {
  edge_df <- edges_to_use |>
    dplyr::transmute(
      from = as.character(taxon1),
      to = as.character(taxon2),
      rho = as.numeric(rho),
      q_value = as.numeric(q_value),
      association = as.character(association),
      weight = as.numeric(abs_rho)
    )

  vertex_df <- node_base |>
    dplyr::mutate(Taxon = as.character(Taxon)) |>
    dplyr::rename(name = Taxon)

  g <- igraph::graph_from_data_frame(
    d = edge_df,
    directed = FALSE,
    vertices = vertex_df
  )

  node_names <- igraph::V(g)$name

  # Module detection uses only non-isolated nodes, but isolates are retained in
  # the final node table and in the graph object.
  deg_all <- igraph::degree(g, v = igraph::V(g), mode = "all")
  names(deg_all) <- node_names

  non_isolate_vertices <- igraph::V(g)[deg_all > 0]
  module_vec <- rep(NA_integer_, length(node_names))
  names(module_vec) <- node_names

  g_for_modules <- igraph::induced_subgraph(g, vids = non_isolate_vertices)

  if (igraph::vcount(g_for_modules) >= 3 && igraph::ecount(g_for_modules) >= 2) {
    comm <- tryCatch(
      igraph::cluster_louvain(g_for_modules, weights = igraph::E(g_for_modules)$weight),
      error = function(e) {
        log_msg("cluster_louvain failed; assigning all connected taxa to module 1:", e$message)
        NULL
      }
    )

    if (is.null(comm)) {
      module_vec[igraph::V(g_for_modules)$name] <- 1L
    } else {
      module_vec[igraph::V(g_for_modules)$name] <- as.integer(igraph::membership(comm))
    }
  } else if (igraph::vcount(g_for_modules) > 0) {
    module_vec[igraph::V(g_for_modules)$name] <- 1L
  }

  igraph::V(g)$module <- module_vec[node_names]

  str_all <- igraph::strength(g, v = igraph::V(g), weights = igraph::E(g)$weight)
  names(str_all) <- node_names

  # Robust positive/negative edge counts by node.
  edge_assoc <- igraph::edge_attr(g, "association")
  positive_degree <- setNames(integer(length(node_names)), node_names)
  negative_degree <- setNames(integer(length(node_names)), node_names)

  for (v_name in node_names) {
    v_id <- match(v_name, node_names)
    eids <- as.integer(igraph::incident(g, v = v_id, mode = "all"))
    if (length(eids) > 0) {
      assoc_v <- edge_assoc[eids]
      positive_degree[v_name] <- sum(assoc_v == "positive", na.rm = TRUE)
      negative_degree[v_name] <- sum(assoc_v == "negative", na.rm = TRUE)
    }
  }

  # Zi/Pi calculation using absolute rho as edge strength.
  modules <- module_vec[node_names]
  names(modules) <- node_names

  adj_abs <- as.matrix(igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  rownames(adj_abs) <- node_names
  colnames(adj_abs) <- node_names

  within_strength <- setNames(numeric(length(node_names)), node_names)
  for (v_name in node_names) {
    if (is.na(modules[v_name])) {
      within_strength[v_name] <- 0
    } else {
      same_mod <- names(modules)[!is.na(modules) & modules == modules[v_name]]
      within_strength[v_name] <- sum(adj_abs[v_name, same_mod, drop = TRUE], na.rm = TRUE)
    }
  }

  total_strength <- rowSums(adj_abs, na.rm = TRUE)
  total_strength <- total_strength[node_names]

  Zi <- setNames(numeric(length(node_names)), node_names)
  for (m in unique(stats::na.omit(modules))) {
    idx <- names(modules)[!is.na(modules) & modules == m]
    vals <- within_strength[idx]
    if (length(vals) > 1 && is.finite(stats::sd(vals)) && stats::sd(vals) > 0) {
      Zi[idx] <- (vals - mean(vals)) / stats::sd(vals)
    } else {
      Zi[idx] <- 0
    }
  }

  Pi <- setNames(numeric(length(node_names)), node_names)
  module_levels <- unique(stats::na.omit(modules))
  for (v_name in node_names) {
    k_total <- total_strength[v_name]
    if (!is.finite(k_total) || k_total <= 0 || length(module_levels) == 0) {
      Pi[v_name] <- 0
    } else {
      frac_sq <- vapply(module_levels, function(m) {
        idx <- names(modules)[!is.na(modules) & modules == m]
        (sum(adj_abs[v_name, idx], na.rm = TRUE) / k_total)^2
      }, numeric(1))
      Pi[v_name] <- 1 - sum(frac_sq)
    }
  }

  role <- dplyr::case_when(
    Zi >= 2.5 & Pi > 0.62 ~ "network hub",
    Zi >= 2.5 & Pi <= 0.62 ~ "module hub",
    Zi < 2.5 & Pi > 0.62 ~ "connector",
    TRUE ~ "peripheral"
  )

  # Isolated nodes should not be classified as peripheral hubs/connectors.
  role[deg_all[node_names] == 0] <- "isolated/no selected edge"

  node_metrics <- tibble::tibble(
    Taxon = node_names,
    module = as.integer(modules[node_names]),
    degree = as.integer(deg_all[node_names]),
    strength = as.numeric(str_all[node_names]),
    positive_degree = as.integer(positive_degree[node_names]),
    negative_degree = as.integer(negative_degree[node_names]),
    Zi = as.numeric(Zi[node_names]),
    Pi = as.numeric(Pi[node_names]),
    role = as.character(role)
  )

  nodes_final <- node_base |>
    dplyr::left_join(node_metrics, by = "Taxon") |>
    dplyr::mutate(
      module = dplyr::if_else(is.na(module), NA_integer_, module),
      degree = tidyr::replace_na(degree, 0L),
      strength = tidyr::replace_na(strength, 0),
      positive_degree = tidyr::replace_na(positive_degree, 0L),
      negative_degree = tidyr::replace_na(negative_degree, 0L),
      Zi = tidyr::replace_na(Zi, 0),
      Pi = tidyr::replace_na(Pi, 0),
      role = tidyr::replace_na(role, "isolated/no selected edge")
    )

  # Reattach the cleaned node attributes to g for plotting.
  for (col in setdiff(colnames(nodes_final), "Taxon")) {
    vals <- nodes_final[[col]][match(node_names, nodes_final$Taxon)]
    igraph::vertex_attr(g, col, index = igraph::V(g)) <- vals
  }

  write_table(nodes_final, paste0("05_coabundance_nodes_", network_suffix, ".csv"))
  write_table(edges_to_use, paste0("05_coabundance_edges_used_for_network_", network_suffix, ".csv"))

  module_summary <- nodes_final |>
    dplyr::filter(!is.na(module), degree > 0) |>
    dplyr::group_by(module) |>
    dplyr::summarise(
      n_taxa = dplyr::n(),
      mean_degree = mean(degree),
      mean_strength = mean(strength),
      dominant_compartment_classes = paste(names(sort(table(compartment_class), decreasing = TRUE)), collapse = "; "),
      top_taxa_by_strength = paste(Taxon[order(strength, decreasing = TRUE)][seq_len(min(8, dplyr::n()))], collapse = "; "),
      .groups = "drop"
    ) |>
    dplyr::arrange(module)

  write_table(module_summary, paste0("05_coabundance_module_summary_", network_suffix, ".csv"))
}

############################
# 8. PLOTS
############################

# 8A. Correlation heatmap.
log_msg("Building correlation heatmap")
heat_order <- colnames(cor_mat)
if (ncol(cor_mat) >= 3) {
  dist_mat <- as.dist(1 - cor_mat)
  hc <- hclust(dist_mat, method = "average")
  heat_order <- hc$labels[hc$order]
}

cor_long <- as.data.frame(cor_mat[heat_order, heat_order, drop = FALSE]) |>
  tibble::rownames_to_column("Taxon_y") |>
  tidyr::pivot_longer(-Taxon_y, names_to = "Taxon_x", values_to = "rho") |>
  dplyr::mutate(
    Taxon_x = factor(Taxon_x, levels = heat_order),
    Taxon_y = factor(Taxon_y, levels = rev(heat_order))
  )

p_heat <- ggplot2::ggplot(cor_long, ggplot2::aes(x = Taxon_x, y = Taxon_y, fill = rho)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.05) +
  ggplot2::scale_fill_gradient2(
    low = "#2E6F73", mid = "white", high = "#7A1F2B", midpoint = 0,
    limits = c(-1, 1), name = "Spearman\nrho"
  ) +
  ggplot2::coord_fixed() +
  theme_net(8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = ggplot2::element_text(size = 6),
    axis.title = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    title = "Host genus-level CLR-Spearman co-abundance matrix",
    subtitle = paste0(length(selected_taxa), " prevalent taxa; correlations are exploratory co-abundance, not direct interactions")
  )

save_plot_multi(p_heat, "FigS_host_genus_coabundance_correlation_heatmap", width = 10, height = 9)

# 8B. Network plot.
if (exists("g") && igraph::ecount(g) > 0) {
  log_msg("Building network plot with", igraph::vcount(g), "nodes and", igraph::ecount(g), "edges")

  # Attach updated node attributes to graph for ggraph.
  # Use vertex_attr() rather than V(g)[name][[col]] to avoid igraph
  # subsetting errors in some Windows/igraph versions.
  node_attr <- nodes_final |> dplyr::rename(name = Taxon)
  graph_node_names <- igraph::V(g)$name
  for (col in setdiff(colnames(node_attr), "name")) {
    vals <- node_attr[[col]][match(graph_node_names, node_attr$name)]
    igraph::vertex_attr(g, col, index = igraph::V(g)) <- vals
  }

  # Label strongest / most abundant nodes.
  label_nodes <- nodes_final |>
    dplyr::filter(degree > 0) |>
    dplyr::arrange(dplyr::desc(strength), dplyr::desc(mean_rel_abund)) |>
    dplyr::slice_head(n = label_top_n_network_nodes) |>
    dplyr::pull(Taxon)

  set.seed(1001)
  p_net <- ggraph::ggraph(g, layout = network_layout, weights = igraph::E(g)$weight) +
    ggraph::geom_edge_link(
      ggplot2::aes(edge_alpha = weight, edge_width = weight, edge_colour = association),
      show.legend = TRUE
    ) +
    ggraph::scale_edge_colour_manual(values = edge_pal, name = "Association") +
    ggraph::scale_edge_alpha(range = c(0.25, 0.85), guide = "none") +
    ggraph::scale_edge_width(range = c(0.25, 1.6), guide = "none") +
    ggraph::geom_node_point(
      ggplot2::aes(size = prevalence_prop, fill = compartment_class),
      shape = 21, colour = "grey20", stroke = 0.25, alpha = 0.95
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = ifelse(name %in% label_nodes, name, "")),
      repel = TRUE, size = 2.6, colour = "grey10", max.overlaps = 100
    ) +
    ggplot2::scale_fill_manual(values = compartment_pal, drop = FALSE, name = "Dominant\ncompartment") +
    ggplot2::scale_size_continuous(range = c(2.5, 8), labels = scales::percent_format(accuracy = 1), name = "Prevalence") +
    ggraph::theme_graph(base_size = 10) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9)
    ) +
    ggplot2::labs(
      title = "Exploratory host genus co-abundance network",
      subtitle = paste0(threshold_label, "; node colour indicates dominant host compartment")
    )

  save_plot_multi(p_net, paste0("FigS_host_genus_coabundance_network_", network_suffix), width = 10.5, height = 8)

  # 8C. Zi-Pi plot.
  p_zipi <- nodes_final |>
    dplyr::filter(degree > 0) |>
    ggplot2::ggplot(ggplot2::aes(x = Pi, y = Zi, fill = compartment_class, size = strength, label = Taxon)) +
    ggplot2::geom_hline(yintercept = 2.5, linetype = 2, colour = "grey55", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0.62, linetype = 2, colour = "grey55", linewidth = 0.35) +
    ggplot2::geom_point(shape = 21, colour = "grey20", stroke = 0.25, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = function(d) d |> dplyr::arrange(dplyr::desc(strength)) |> dplyr::slice_head(n = 20),
      size = 2.6, max.overlaps = 50
    ) +
    ggplot2::scale_fill_manual(values = compartment_pal, drop = FALSE, name = "Dominant\ncompartment") +
    ggplot2::scale_size_continuous(range = c(2.0, 7.0), name = "Node strength") +
    theme_net(10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)) +
    ggplot2::labs(
      title = "Network role classification of prevalent host genera",
      subtitle = "Zi = within-module connectivity; Pi = among-module participation; exploratory only",
      x = "Participation coefficient (Pi)",
      y = "Within-module connectivity (Zi)"
    )

  save_plot_multi(p_zipi, paste0("FigS_host_genus_coabundance_ZiPi_", network_suffix), width = 8, height = 6.2)
} else {
  log_msg("Skipping network and Zi-Pi plots because no selected edges were available.")
}

############################
# 9. SUMMARY NOTES
############################

notes_file <- file.path(network_tables_dir, "05_host_genus_coabundance_interpretation_notes.txt")
cat(
  "Host genus co-abundance network - interpretation notes\n",
  "====================================================\n\n",
  "This is an exploratory co-abundance analysis, not direct evidence of ecological interaction.\n",
  "The analysis uses prevalent host-associated taxa, relative abundances transformed with CLR, and Spearman correlations across host samples.\n",
  "Edges are retained after FDR correction and absolute rho filtering.\n\n",
  "Recommended wording:\n",
  "Exploratory co-abundance analysis of prevalent host-associated genera revealed recurrent correlation structure among dominant taxa. Because these networks are based on compositional 16S profiles, they were interpreted as co-abundance patterns rather than direct microbial interactions.\n\n",
  "Main settings:\n",
  paste0("Taxonomic rank: ", taxon_rank, "\n"),
  paste0("Minimum prevalence proportion: ", min_prevalence_prop, "\n"),
  paste0("Minimum mean relative abundance: ", min_mean_rel_abund, "\n"),
  paste0("Maximum taxa retained: ", max_taxa_for_network, "\n"),
  paste0("Strict edge threshold: |rho| >= ", edge_abs_rho_threshold, ", q <= ", edge_q_threshold, "\n"),
  paste0("Network plotted: ", network_suffix, " (", threshold_label, ")\n\n"),
  "Key files:\n",
  "05_coabundance_taxa_used_for_network.csv\n",
  "05_coabundance_correlations_all_pairs.csv\n",
  paste0("05_coabundance_edges_used_for_network_", network_suffix, ".csv\n"),
  paste0("05_coabundance_nodes_", network_suffix, ".csv\n"),
  paste0("05_coabundance_module_summary_", network_suffix, ".csv\n"),
  "FigS_host_genus_coabundance_correlation_heatmap.*\n",
  paste0("FigS_host_genus_coabundance_network_", network_suffix, ".*\n"),
  paste0("FigS_host_genus_coabundance_ZiPi_", network_suffix, ".*\n"),
  file = notes_file
)

log_msg("Host genus co-abundance network analysis complete.")
log_msg("Tables written to:", network_tables_dir)
log_msg("Plots written to:", network_plots_dir)
log_msg("Most important plots: correlation heatmap, network plot, Zi-Pi plot.")
