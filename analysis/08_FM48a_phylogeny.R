################################################################################
# FM48a PHYLOGENOMIC TREE
# Windows RStudio + Ubuntu-24.04 / WSL
#
# INPUTS FROM THE PREVIOUS ANI ANALYSIS
#   results/genomics/FM48a_genome_taxonomy_ANI/
#       FM48a_Oceanisphaera_ANI_ranked.csv
#       Oceanisphaera_reference_genome_files.csv
#
#   results/genomics_work/FM48a_ANI_WSL_work/
#       FM48a.fasta
#       refs/*.fna
#
# SOFTWARE IN WSL
#   GToTree
#   IQ-TREE
#
# WORKFLOW
#   1) Deduplicate paired GCA/GCF records so each biological strain is used once.
#   2) Prefer RefSeq (GCF) when the same assembly is present as GCA + GCF.
#   3) For repeated strain labels, retain the highest-completeness genome.
#   4) Add FM48a.
#   5) Download type-material genomes from the sibling genera
#      Oceanimonas and Dongshaea for rooting.
#   6) Run GToTree with the bacterial single-copy marker set.
#   7) Infer a maximum-likelihood tree with IQ-TREE using ModelFinder,
#      1000 ultrafast bootstrap replicates and 1000 SH-aLRT replicates.
#   8) Root the final tree on the non-Oceanisphaera outgroups if possible.
#
################################################################################


############################
# 0. SETTINGS
############################

wsl_distro <- Sys.getenv("HALICLYSTUS_WSL_DISTRO", unset = "Ubuntu-24.04")
conda_env <- "fm48_phylo"

cpus_detected <- parallel::detectCores(logical = TRUE)
if (is.na(cpus_detected)) cpus_detected <- 2L
cpus <- max(1L, min(8L, cpus_detected - 1L))

# If TRUE, fresh copies of the outgroup genome packages are downloaded.
force_redownload_outgroups <- FALSE

# Closely related sibling genera used to root the Oceanisphaera tree.
outgroup_genera <- c("Oceanimonas", "Dongshaea")


############################
# 1. FIND PROJECT + PREVIOUS ANI OUTPUTS
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
downloads_dir <- file.path(repo_dir, "data", "genomics")
ani_dir <- file.path(repo_dir, "results", "genomics", "FM48a_genome_taxonomy_ANI")

ani_csv <- file.path(
  ani_dir,
  "FM48a_Oceanisphaera_ANI_ranked.csv"
)

ref_map_csv <- file.path(
  ani_dir,
  "Oceanisphaera_reference_genome_files.csv"
)

if (!file.exists(ani_csv)) {
  stop("Missing previous ANI result:\n  ", ani_csv)
}

if (!file.exists(ref_map_csv)) {
  stop("Missing previous reference map:\n  ", ref_map_csv)
}

ani <- read.csv(
  ani_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ref_map <- read.csv(
  ref_map_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_ani <- c(
  "accession",
  "organism_name",
  "strain",
  "checkm_completeness",
  "checkm_contamination",
  "is_type_material"
)

missing_ani <- setdiff(
  required_ani,
  colnames(ani)
)

if (length(missing_ani) > 0) {
  stop(
    "ANI table is missing columns: ",
    paste(missing_ani, collapse = ", ")
  )
}

if (!all(c("accession", "staged_path") %in% colnames(ref_map))) {
  stop(
    "Reference map must contain accession and staged_path."
  )
}


############################
# 2. OUTPUT / WORK FOLDERS
############################

final_dir <- file.path(
  base_dir,
  "FM48a_phylogenomics"
)

dir.create(
  final_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Keep WSL work files in an output folder under this checkout.
work_dir <- file.path(repo_dir, "results", "genomics_work", "FM48a_phylogenomics_WSL_work")

stage_dir <- file.path(
  work_dir,
  "genomes"
)

outgroup_pkg_dir <- file.path(
  work_dir,
  "outgroup_packages"
)

gtotree_dir <- file.path(
  work_dir,
  "GToTree_output"
)

iqtree_dir <- file.path(
  work_dir,
  "IQTREE_output"
)

for (d in c(
  work_dir,
  stage_dir,
  outgroup_pkg_dir,
  iqtree_dir
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

cat("\n============================================================\n")
cat("FM48a PHYLOGENOMICS\n")
cat("============================================================\n")
cat("Final project folder:\n  ",
    normalizePath(final_dir, winslash = "/"), "\n", sep = "")
cat("WSL work folder:\n  ",
    normalizePath(work_dir, winslash = "/"), "\n\n", sep = "")


############################
# 3. FIND WINDOWS + WSL TOOLS
############################

datasets_candidates <- Sys.which(c("datasets.exe", "datasets"))
datasets_bin <- unname(datasets_candidates[nzchar(datasets_candidates)])[1]
if (is.na(datasets_bin) || !nzchar(datasets_bin)) {
  stop("NCBI datasets CLI is required on PATH (datasets or datasets.exe).")
}

wsl_bin <- Sys.which("wsl.exe")
if (!nzchar(wsl_bin)) {
  wsl_bin <- Sys.which("wsl")
}
if (!nzchar(wsl_bin)) {
  stop("Could not find wsl.exe.")
}


############################
# 4. WINDOWS -> WSL PATH CONVERTER
############################

windows_to_wsl <- function(path) {

  p <- normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )

  if (!grepl("^[A-Za-z]:/", p)) {
    stop("Cannot convert Windows path to WSL:\n", p)
  }

  drive <- tolower(substr(p, 1, 1))
  rest <- substr(p, 3, nchar(p))

  paste0(
    "/mnt/",
    drive,
    rest
  )
}


############################
# 5. FIND CONDA INSIDE WSL
############################

# Discover conda in the current WSL user's environment, without assuming a login name.
conda_command <- paste(
  "command -v conda 2>/dev/null ||",
  "for d in \"$HOME/miniforge3\" \"$HOME/mambaforge\" \"$HOME/miniconda3\"; do",
  "if test -x \"$d/bin/conda\"; then echo \"$d/bin/conda\"; break; fi;",
  "done"
)
conda_check <- suppressWarnings(system2(
  wsl_bin, args = c("-d", wsl_distro, "--", "bash", "-lc", shQuote(conda_command)),
  stdout = TRUE, stderr = TRUE
))
conda_paths <- trimws(conda_check[grepl("^/", trimws(conda_check))])
if (!length(conda_paths)) {
  stop("Cannot locate conda inside the selected WSL distribution. ",
       "Install conda and the fm48_phylo environment before running this script.")
}
conda_wsl <- conda_paths[1]

cat("Using WSL conda:\n  ", conda_wsl, "\n\n", sep = "")


############################
# 6. CHECK PHYLOGENOMICS ENVIRONMENT
############################

wsl_capture <- function(command) {

  system2(
    wsl_bin,
    args = c(
      "-d",
      wsl_distro,
      "--",
      "bash",
      "-lc",
      shQuote(command)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
}

gtotree_check <- wsl_capture(
  paste(
    shQuote(conda_wsl),
    "run -n",
    shQuote(conda_env),
    "bash -lc",
    shQuote("command -v GToTree || true")
  )
)

if (!any(grepl("GToTree", gtotree_check))) {
  stop(
    "GToTree was not found in conda environment '",
    conda_env,
    "'.\n",
    "Run the Ubuntu setup script first."
  )
}

iqtree_name <- NULL

for (candidate in c("iqtree3", "iqtree2", "iqtree")) {

  chk <- wsl_capture(
    paste(
      shQuote(conda_wsl),
      "run -n",
      shQuote(conda_env),
      "bash -lc",
      shQuote(
        paste(
          "command -v",
          candidate,
          "|| true"
        )
      )
    )
  )

  if (any(nzchar(chk))) {
    iqtree_name <- candidate
    break
  }
}

if (is.null(iqtree_name)) {
  stop(
    "IQ-TREE was not found in conda environment '",
    conda_env,
    "'."
  )
}

cat("GToTree: OK\n")
cat("IQ-TREE executable:", iqtree_name, "\n\n")


############################
# 7. DEDUPLICATE OCEANISPHAERA GENOMES
############################

dat <- merge(
  ani,
  ref_map[, c("accession", "staged_path")],
  by = "accession",
  all.x = TRUE,
  sort = FALSE
)

dat <- dat[
  !is.na(dat$staged_path) &
    file.exists(dat$staged_path),
  ,
  drop = FALSE
]

# Biological grouping:
# When strain exists, species + strain defines the biological isolate.
# When strain is absent, retain accessions separately.
strain_clean <- trimws(
  ifelse(
    is.na(dat$strain),
    "",
    dat$strain
  )
)

dat$biological_key <- ifelse(
  nzchar(strain_clean),
  paste(
    dat$organism_name,
    strain_clean,
    sep = " || "
  ),
  paste(
    dat$organism_name,
    dat$accession,
    sep = " || "
  )
)

# Prefer:
#   1) GCF over GCA
#   2) higher CheckM completeness
#   3) lower CheckM contamination
dat$is_gcf <- grepl(
  "^GCF_",
  dat$accession
)

dat$completeness_for_sort <- ifelse(
  is.na(dat$checkm_completeness),
  -Inf,
  dat$checkm_completeness
)

dat$contamination_for_sort <- ifelse(
  is.na(dat$checkm_contamination),
  Inf,
  dat$checkm_contamination
)

dat <- dat[
  order(
    dat$biological_key,
    -dat$is_gcf,
    -dat$completeness_for_sort,
    dat$contamination_for_sort
  ),
  ,
  drop = FALSE
]

ocean_selected <- dat[
  !duplicated(dat$biological_key),
  ,
  drop = FALSE
]

cat(
  "Oceanisphaera ANI rows:",
  nrow(ani),
  "\n"
)

cat(
  "Biological Oceanisphaera genomes after deduplication:",
  nrow(ocean_selected),
  "\n\n"
)

write.csv(
  ocean_selected,
  file.path(
    final_dir,
    "Oceanisphaera_phylogenomics_genomes_selected.csv"
  ),
  row.names = FALSE
)


############################
# 8. SAFE TREE LABELS
############################

safe_label <- function(x) {

  x <- gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    x
  )

  x <- gsub(
    "_+",
    "_",
    x
  )

  gsub(
    "^_|_$",
    "",
    x
  )
}

ocean_selected$tree_label <- vapply(
  seq_len(nrow(ocean_selected)),
  function(i) {

    org <- ocean_selected$organism_name[i]
    strain <- ocean_selected$strain[i]
    acc <- ocean_selected$accession[i]

    if (is.na(strain) || !nzchar(trimws(strain))) {
      strain <- acc
    }

    safe_label(
      paste0(
        org,
        "__",
        strain,
        "__",
        acc
      )
    )
  },
  FUN.VALUE = character(1)
)


############################
# 9. CLEAN + STAGE OCEANISPHAERA GENOMES
############################

# Remove old staged genome files.
old_stage <- list.files(
  stage_dir,
  full.names = TRUE
)

if (length(old_stage) > 0) {
  unlink(
    old_stage,
    recursive = TRUE,
    force = TRUE
  )
}

ocean_selected$phylo_path <- file.path(
  stage_dir,
  paste0(
    ocean_selected$tree_label,
    ".fna"
  )
)

copy_ok <- mapply(
  file.copy,
  from = ocean_selected$staged_path,
  to = ocean_selected$phylo_path,
  MoreArgs = list(overwrite = TRUE)
)

if (!all(copy_ok)) {
  stop(
    "At least one Oceanisphaera reference genome could not be staged."
  )
}


############################
# 10. ADD FM48a
############################

fm48_candidates <- unique(c(
  file.path(repo_dir, "results", "genomics_work", "FM48a_ANI_WSL_work", "FM48a.fasta"),
  file.path(repo_dir, "data", "genomics", "ZTRMXK_1_FM_48a_1_reference.fasta")
))

fm48_candidates <- fm48_candidates[
  file.exists(fm48_candidates)
]

if (length(fm48_candidates) == 0) {
  stop("Could not locate the FM48a genome FASTA.")
}

fm48_stage <- file.path(
  stage_dir,
  "FM48a__Oceanisphaera_ostreae_candidate.fna"
)

file.copy(
  fm48_candidates[1],
  fm48_stage,
  overwrite = TRUE
)


############################
# 11. DOWNLOAD TYPE-MATERIAL OUTGROUPS
############################

download_outgroup <- function(genus) {

  genus_dir <- file.path(
    outgroup_pkg_dir,
    genus
  )

  dir.create(
    genus_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  zip_path <- file.path(
    genus_dir,
    paste0(
      genus,
      "_type_material.zip"
    )
  )

  extract_dir <- file.path(
    genus_dir,
    "package"
  )

  if (
    force_redownload_outgroups ||
    !file.exists(zip_path)
  ) {

    if (file.exists(zip_path)) {
      unlink(zip_path)
    }

    cat(
      "\nDownloading type-material genomes for ",
      genus,
      "...\n",
      sep = ""
    )

    status <- system2(
      datasets_bin,
      args = c(
        "download",
        "genome",
        "taxon",
        genus,
        "--from-type",
        "--include",
        "genome",
        "--mag",
        "exclude",
        "--exclude-atypical",
        "--exclude-multi-isolate",
        "--filename",
        zip_path,
        "--no-progressbar"
      )
    )

    if (!identical(status, 0L)) {
      warning(
        "NCBI download failed for outgroup genus ",
        genus
      )
      return(character())
    }
  }

  if (
    force_redownload_outgroups &&
    dir.exists(extract_dir)
  ) {
    unlink(
      extract_dir,
      recursive = TRUE,
      force = TRUE
    )
  }

  if (!dir.exists(extract_dir)) {
    dir.create(
      extract_dir,
      recursive = TRUE
    )
    unzip(
      zip_path,
      exdir = extract_dir
    )
  }

  fastas <- list.files(
    extract_dir,
    pattern = "_genomic\\.fna$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(fastas) == 0) {
    return(character())
  }

  # Keep one genome per accession directory.
  get_acc <- function(p) {

    parts <- strsplit(
      normalizePath(
        p,
        winslash = "/"
      ),
      "/",
      fixed = TRUE
    )[[1]]

    h <- grep(
      "^GC[AF]_[0-9]+\\.[0-9]+$",
      parts,
      value = TRUE
    )

    if (length(h) > 0) {
      h[length(h)]
    } else {
      tools::file_path_sans_ext(
        basename(p)
      )
    }
  }

  acc <- vapply(
    fastas,
    get_acc,
    FUN.VALUE = character(1)
  )

  keep <- !duplicated(acc)

  fastas <- fastas[keep]
  acc <- acc[keep]

  staged <- character(
    length(fastas)
  )

  for (i in seq_along(fastas)) {

    out <- file.path(
      stage_dir,
      paste0(
        "OUTGROUP__",
        safe_label(genus),
        "__",
        safe_label(acc[i]),
        ".fna"
      )
    )

    file.copy(
      fastas[i],
      out,
      overwrite = TRUE
    )

    staged[i] <- out
  }

  staged
}

outgroup_files <- unlist(
  lapply(
    outgroup_genera,
    download_outgroup
  ),
  use.names = FALSE
)

if (length(outgroup_files) == 0) {
  warning(
    "No outgroup genome could be downloaded. ",
    "The tree can still be inferred but will not be outgroup-rooted."
  )
}

cat(
  "\nOutgroup genomes staged:",
  length(outgroup_files),
  "\n"
)


############################
# 12. BUILD GENOME MANIFEST
############################

all_genomes <- list.files(
  stage_dir,
  pattern = "\\.fna$",
  full.names = TRUE
)

if (length(all_genomes) < 5) {
  stop(
    "Too few genomes were staged for phylogenomics."
  )
}

manifest <- data.frame(
  tree_label = tools::file_path_sans_ext(
    basename(all_genomes)
  ),
  windows_path = normalizePath(
    all_genomes,
    winslash = "/"
  ),
  stringsAsFactors = FALSE
)

manifest$group <- ifelse(
  grepl(
    "^FM48a",
    manifest$tree_label
  ),
  "FM48a",
  ifelse(
    grepl(
      "^OUTGROUP",
      manifest$tree_label
    ),
    "Outgroup",
    "Oceanisphaera_reference"
  )
)

manifest$wsl_path <- vapply(
  manifest$windows_path,
  windows_to_wsl,
  FUN.VALUE = character(1)
)

write.csv(
  manifest,
  file.path(
    final_dir,
    "FM48a_phylogenomics_manifest.csv"
  ),
  row.names = FALSE
)

genome_list_win <- file.path(
  work_dir,
  "genome_files.txt"
)

writeLines(
  manifest$wsl_path,
  genome_list_win
)

genome_list_wsl <- windows_to_wsl(
  genome_list_win
)

cat(
  "Total genomes entering phylogenomics:",
  nrow(manifest),
  "\n"
)


############################
# 13. RUN GToTree
############################

if (dir.exists(gtotree_dir)) {
  unlink(
    gtotree_dir,
    recursive = TRUE,
    force = TRUE
  )
}

gtotree_wsl <- windows_to_wsl(
  gtotree_dir
)

cat("\n============================================================\n")
cat("RUNNING GToTree\n")
cat("============================================================\n")

gtotree_command <- paste(
  shQuote(conda_wsl),
  "run -n",
  shQuote(conda_env),
  "GToTree",
  "-f",
  shQuote(genome_list_wsl),
  "-H Bacteria",
  "-j",
  cpus,
  "-o",
  shQuote(gtotree_wsl)
)

cat(
  gtotree_command,
  "\n\n"
)

gtotree_log <- wsl_capture(
  gtotree_command
)

writeLines(
  gtotree_log,
  file.path(
    final_dir,
    "GToTree_console_log.txt"
  )
)

cat(
  paste(
    tail(
      gtotree_log,
      50
    ),
    collapse = "\n"
  ),
  "\n"
)

alignment_win <- file.path(
  gtotree_dir,
  "Aligned_SCGs.faa"
)

if (!file.exists(alignment_win)) {
  stop(
    "GToTree did not create Aligned_SCGs.faa.\n",
    "Inspect GToTree_console_log.txt."
  )
}


############################
# 14. GToTree QC FILES
############################

gtt_summary <- file.path(
  gtotree_dir,
  "Genomes_summary_info.tsv"
)

gtt_hits <- file.path(
  gtotree_dir,
  "SCG_hit_counts.tsv"
)

if (file.exists(gtt_summary)) {
  file.copy(
    gtt_summary,
    file.path(
      final_dir,
      "Genomes_summary_info.tsv"
    ),
    overwrite = TRUE
  )
}

if (file.exists(gtt_hits)) {
  file.copy(
    gtt_hits,
    file.path(
      final_dir,
      "SCG_hit_counts.tsv"
    ),
    overwrite = TRUE
  )
}

citations_file <- file.path(
  gtotree_dir,
  "citations.txt"
)

if (file.exists(citations_file)) {
  file.copy(
    citations_file,
    file.path(
      final_dir,
      "GToTree_citations.txt"
    ),
    overwrite = TRUE
  )
}


############################
# 15. RUN IQ-TREE
############################

alignment_wsl <- windows_to_wsl(
  alignment_win
)

iqtree_prefix_win <- file.path(
  iqtree_dir,
  "FM48a_Oceanisphaera_phylogenomics"
)

iqtree_prefix_wsl <- windows_to_wsl(
  iqtree_prefix_win
)

cat("\n============================================================\n")
cat("RUNNING IQ-TREE\n")
cat("============================================================\n")

iqtree_command <- paste(
  shQuote(conda_wsl),
  "run -n",
  shQuote(conda_env),
  iqtree_name,
  "-s",
  shQuote(alignment_wsl),
  "-m MFP",
  "-B 1000",
  "--alrt 1000",
  "-T AUTO",
  "--prefix",
  shQuote(iqtree_prefix_wsl)
)

cat(
  iqtree_command,
  "\n\n"
)

iqtree_log <- wsl_capture(
  iqtree_command
)

writeLines(
  iqtree_log,
  file.path(
    final_dir,
    "IQTREE_console_log.txt"
  )
)

cat(
  paste(
    tail(
      iqtree_log,
      50
    ),
    collapse = "\n"
  ),
  "\n"
)

treefile_win <- paste0(
  iqtree_prefix_win,
  ".treefile"
)

if (!file.exists(treefile_win)) {
  stop(
    "IQ-TREE did not create the final .treefile.\n",
    "Inspect IQTREE_console_log.txt."
  )
}

file.copy(
  treefile_win,
  file.path(
    final_dir,
    "FM48a_Oceanisphaera_phylogenomics_unrooted.treefile"
  ),
  overwrite = TRUE
)

iqtree_report <- paste0(
  iqtree_prefix_win,
  ".iqtree"
)

if (file.exists(iqtree_report)) {
  file.copy(
    iqtree_report,
    file.path(
      final_dir,
      "FM48a_Oceanisphaera_phylogenomics.iqtree"
    ),
    overwrite = TRUE
  )
}


############################
# 16. ROOT TREE ON OUTGROUPS
############################

if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}

tree <- ape::read.tree(
  treefile_win
)

outgroup_tips <- tree$tip.label[
  grepl(
    "^OUTGROUP",
    tree$tip.label
  )
]

rooted_tree <- tree

if (length(outgroup_tips) > 0) {

  rooted_tree <- ape::root(
    tree,
    outgroup = outgroup_tips,
    resolve.root = TRUE
  )

  cat(
    "\nTree rooted using",
    length(outgroup_tips),
    "outgroup tip(s).\n"
  )

} else {

  warning(
    "No OUTGROUP labels were found in the final tree. ",
    "Leaving the IQ-TREE topology unrooted."
  )
}

rooted_tree_path <- file.path(
  final_dir,
  "FM48a_Oceanisphaera_phylogenomics_ROOTED.tree"
)

ape::write.tree(
  rooted_tree,
  file = rooted_tree_path
)


############################
# 17. SIMPLE TREE PDF + PNG
############################

# The final styling can later be improved for the manuscript.
# This preview intentionally preserves all tree labels and node support values.

plot_tree <- function(device_fun, filename, width, height) {

  device_fun(
    filename,
    width = width,
    height = height
  )

  par(
    mar = c(1, 1, 3, 1)
  )

  ape::plot.phylo(
    rooted_tree,
    type = "phylogram",
    show.tip.label = TRUE,
    cex = 0.65,
    no.margin = FALSE,
    main = "Phylogenomic placement of FM48a within Oceanisphaera"
  )

  if (!is.null(rooted_tree$node.label)) {
    ape::nodelabels(
      rooted_tree$node.label,
      frame = "none",
      cex = 0.45,
      adj = c(1.1, -0.3)
    )
  }

  dev.off()
}

pdf(
  file.path(
    final_dir,
    "FM48a_Oceanisphaera_phylogenomics_preview.pdf"
  ),
  width = 12,
  height = 9
)

par(
  mar = c(1, 1, 3, 1)
)

ape::plot.phylo(
  rooted_tree,
  type = "phylogram",
  show.tip.label = TRUE,
  cex = 0.65,
  no.margin = FALSE,
  main = "Phylogenomic placement of FM48a within Oceanisphaera"
)

if (!is.null(rooted_tree$node.label)) {
  ape::nodelabels(
    rooted_tree$node.label,
    frame = "none",
    cex = 0.45,
    adj = c(1.1, -0.3)
  )
}

dev.off()

png(
  file.path(
    final_dir,
    "FM48a_Oceanisphaera_phylogenomics_preview.png"
  ),
  width = 2400,
  height = 1800,
  res = 200
)

par(
  mar = c(1, 1, 3, 1)
)

ape::plot.phylo(
  rooted_tree,
  type = "phylogram",
  show.tip.label = TRUE,
  cex = 0.65,
  no.margin = FALSE,
  main = "Phylogenomic placement of FM48a within Oceanisphaera"
)

if (!is.null(rooted_tree$node.label)) {
  ape::nodelabels(
    rooted_tree$node.label,
    frame = "none",
    cex = 0.45,
    adj = c(1.1, -0.3)
  )
}

dev.off()


############################
# 18. TIP ANNOTATION TABLE
############################

tip_annotation <- data.frame(
  tree_tip = rooted_tree$tip.label,
  group = ifelse(
    grepl(
      "^FM48a",
      rooted_tree$tip.label
    ),
    "FM48a",
    ifelse(
      grepl(
        "^OUTGROUP",
        rooted_tree$tip.label
      ),
      "Outgroup",
      "Oceanisphaera"
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  tip_annotation,
  file.path(
    final_dir,
    "FM48a_phylogenomics_tree_tip_annotations.csv"
  ),
  row.names = FALSE
)


############################
# 19. FINAL SUMMARY
############################

summary_lines <- c(
  "FM48a PHYLOGENOMIC TREE",
  "=======================",
  "",
  paste0(
    "Oceanisphaera biological genomes after deduplication: ",
    nrow(ocean_selected)
  ),
  paste0(
    "Outgroup genomes included: ",
    length(outgroup_files)
  ),
  paste0(
    "Total genomes submitted to GToTree: ",
    nrow(manifest)
  ),
  "",
  "Marker workflow:",
  "  GToTree bacterial single-copy marker set",
  "  Individual marker alignment + trimming + concatenation",
  "",
  "Maximum-likelihood inference:",
  paste0(
    "  IQ-TREE executable: ",
    iqtree_name
  ),
  "  ModelFinder: MFP",
  "  Ultrafast bootstrap: 1000",
  "  SH-aLRT: 1000",
  "",
  paste0(
    "Outgroup-rooted: ",
    length(outgroup_tips) > 0
  ),
  "",
  "Main outputs:",
  "  FM48a_Oceanisphaera_phylogenomics_ROOTED.tree",
  "  FM48a_Oceanisphaera_phylogenomics_unrooted.treefile",
  "  FM48a_Oceanisphaera_phylogenomics.iqtree",
  "  FM48a_Oceanisphaera_phylogenomics_preview.pdf",
  "  FM48a_Oceanisphaera_phylogenomics_preview.png",
  "  FM48a_phylogenomics_manifest.csv",
  "  Genomes_summary_info.tsv",
  "  SCG_hit_counts.tsv",
  "  GToTree_citations.txt"
)

writeLines(
  summary_lines,
  file.path(
    final_dir,
    "FM48a_PHYLOGENOMICS_SUMMARY.txt"
  )
)

cat("\n============================================================\n")
cat("PHYLOGENOMICS COMPLETE\n")
cat("============================================================\n")
cat(
  paste(
    summary_lines,
    collapse = "\n"
  ),
  "\n"
)

cat(
  "\nResults folder:\n  ",
  normalizePath(
    final_dir,
    winslash = "/"
  ),
  "\n",
  sep = ""
)
