# FM48a analysis script — repository-relative input/output directories.
# External executables and genome data must be installed/provided. 


################################################################################
# FM48a GENOME TAXONOMY — WINDOWS R + UBUNTU/WSL FASTANI
#
# Requires Windows with WSL and standard command-line tools:
#   WSL distro:   Ubuntu-24.04 (default; configurable)
#   NCBI tools:   datasets and dataformat available on PATH
#   FastANI:      required inside Ubuntu/WSL
#
# WHAT THIS SCRIPT DOES
#   1) Finds the FM48a assembly.
#   2) Downloads current Oceanisphaera genomes from NCBI.
#   3) Extracts NCBI metadata, including type-material information.
#   4) Copies all reference genomes into a simple Windows/WSL staging folder.
#   5) Runs FastANI through Ubuntu/WSL in BOTH directions.
#   6) Calculates reciprocal ANI and reciprocal alignment fractions.
#   7) Ranks all Oceanisphaera references.
#   8) Separately reports the best type-material genome.
#
################################################################################


############################
# 0. SETTINGS
############################

taxon_name <- "Oceanisphaera"
wsl_distro <- Sys.getenv("HALICLYSTUS_WSL_DISTRO", unset = "Ubuntu-24.04")

ani_screen_threshold <- 95
af_screen_threshold <- 0.50

cpus_detected <- parallel::detectCores(logical = TRUE)
if (is.na(cpus_detected)) cpus_detected <- 2L
cpus <- max(1L, min(8L, cpus_detected - 1L))

force_redownload <- FALSE


############################
# 1. WINDOWS PATHS
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

# The genome FASTA goes in data/genomics/; generated ANI results are under results/.
base_dir <- repo_dir
downloads_dir <- file.path(repo_dir, "data", "genomics")
analysis_dir <- file.path(repo_dir, "results", "genomics", "FM48a_genome_taxonomy_ANI")

# NCBI Datasets CLI programs must be installed and available on PATH.
find_executable <- function(candidates) {
  paths <- Sys.which(candidates)
  found <- unname(paths[nzchar(paths)])
  if (!length(found)) stop("Required executable not on PATH: ",
                            paste(candidates, collapse = " or "))
  found[1]
}
datasets_bin <- find_executable(c("datasets.exe", "datasets"))
dataformat_bin <- find_executable(c("dataformat.exe", "dataformat"))
dir.create(
  analysis_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Persistent temporary files for the ANI workflow; never use user-specific paths.
work_dir <- file.path(repo_dir, "results", "genomics_work", "FM48a_ANI_WSL_work")

dir.create(
  work_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

refs_stage_dir <- file.path(
  work_dir,
  "refs"
)

dir.create(
  refs_stage_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n============================================================\n")
cat("FM48a GENOME TAXONOMY — WINDOWS + WSL FASTANI\n")
cat("============================================================\n")
cat("Project:\n  ", normalizePath(base_dir, winslash = "/"), "\n", sep = "")
cat("Final results:\n  ", normalizePath(analysis_dir, winslash = "/"), "\n", sep = "")
cat("WSL work folder:\n  ", normalizePath(work_dir, winslash = "/"), "\n\n", sep = "")


############################
# 2. FIND FM48a GENOME
############################

genome_filename <- "ZTRMXK_1_FM_48a_1_reference.fasta"

genome_candidates <- unique(c(
  file.path(base_dir, genome_filename),
  file.path(downloads_dir, genome_filename),
  file.path(getwd(), genome_filename)
))

genome_candidates <- genome_candidates[
  file.exists(genome_candidates)
]

if (length(genome_candidates) == 0) {
  stop(
    "Could not find the FM48a genome FASTA:\n  ",
    genome_filename
  )
}

query_source <- normalizePath(
  genome_candidates[1],
  winslash = "/",
  mustWork = TRUE
)

# Copy to a simple no-space path for WSL.
query_stage <- file.path(
  work_dir,
  "FM48a.fasta"
)

file.copy(
  query_source,
  query_stage,
  overwrite = TRUE
)

query_stage <- normalizePath(
  query_stage,
  winslash = "/",
  mustWork = TRUE
)

cat("FM48a genome source:\n  ", query_source, "\n", sep = "")
cat("FM48a staged genome:\n  ", query_stage, "\n\n", sep = "")


############################
# 3. BASIC FM48a ASSEMBLY QC
############################

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop(
    "Biostrings is required. Install it before continuing."
  )
}

genome <- Biostrings::readDNAStringSet(query_source)

contig_lengths <- nchar(as.character(genome))
assembly_size <- sum(contig_lengths)
n_contigs <- length(contig_lengths)

sorted_lengths <- sort(
  contig_lengths,
  decreasing = TRUE
)

n50 <- sorted_lengths[
  which(
    cumsum(sorted_lengths) >= assembly_size / 2
  )[1]
]

seq_chars <- strsplit(
  toupper(
    paste(
      as.character(genome),
      collapse = ""
    )
  ),
  "",
  fixed = TRUE
)[[1]]

valid_acgt <- seq_chars %in% c("A", "C", "G", "T")

gc_percent <- 100 * mean(
  seq_chars[valid_acgt] %in% c("G", "C")
)

assembly_stats <- data.frame(
  genome = "FM48a",
  assembly_size_bp = assembly_size,
  n_contigs = n_contigs,
  N50_bp = n50,
  GC_percent = gc_percent,
  stringsAsFactors = FALSE
)

write.csv(
  assembly_stats,
  file.path(
    analysis_dir,
    "FM48a_assembly_stats.csv"
  ),
  row.names = FALSE
)

cat("FM48a assembly size:", format(assembly_size, big.mark = ","), "bp\n")
cat("FM48a contigs:", n_contigs, "\n")
cat("FM48a N50:", format(n50, big.mark = ","), "bp\n")
cat("FM48a GC:", sprintf("%.3f%%", gc_percent), "\n\n")


############################
# 4. CHECK WINDOWS + WSL TOOLS
############################

wsl_bin <- Sys.which("wsl.exe")
if (!nzchar(wsl_bin)) {
  wsl_bin <- Sys.which("wsl")
}

if (!nzchar(wsl_bin)) {
  stop("Windows could not find wsl.exe.")
}

run_capture <- function(command, args) {
  tryCatch(
    system2(
      command,
      args = args,
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) {
      paste0("ERROR: ", conditionMessage(e))
    }
  )
}

datasets_test <- run_capture(
  datasets_bin,
  "--version"
)

dataformat_test <- run_capture(
  dataformat_bin,
  c("tsv", "genome", "--help")
)

fastani_test <- run_capture(
  wsl_bin,
  c(
    "-d",
    wsl_distro,
    "--",
    "fastANI",
    "--version"
  )
)

cat("NCBI datasets:\n")
cat(paste(datasets_test, collapse = "\n"), "\n\n")

cat("NCBI dataformat functional test:\n")
if (any(grepl("Convert", dataformat_test, ignore.case = TRUE)) ||
    any(grepl("Usage", dataformat_test, ignore.case = TRUE)) ||
    any(grepl("genome", dataformat_test, ignore.case = TRUE))) {
  cat("  OK\n\n")
} else {
  cat(paste(head(dataformat_test, 5), collapse = "\n"), "\n\n")
}

cat("FastANI in WSL:\n")
cat(paste(fastani_test, collapse = "\n"), "\n\n")

if (!any(grepl("1\\.33|version", fastani_test, ignore.case = TRUE))) {
  stop(
    "FastANI could not be called through WSL.\n",
    "Expected Ubuntu distro: ",
    wsl_distro
  )
}


############################
# 5. WINDOWS -> WSL PATH CONVERSION
############################

windows_to_wsl <- function(path) {

  p <- normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )

  # Expected form: a Windows drive-letter path such as X:/folder
  if (!grepl("^[A-Za-z]:/", p)) {
    stop(
      "Cannot convert this Windows path to WSL:\n",
      p
    )
  }

  drive <- tolower(substr(p, 1, 1))
  rest <- substr(p, 3, nchar(p))

  paste0(
    "/mnt/",
    drive,
    rest
  )
}

wsl_query <- windows_to_wsl(query_stage)
wsl_work_dir <- windows_to_wsl(work_dir)


############################
# 6. DOWNLOAD OCEANISPHAERA GENOMES FROM NCBI
############################

ncbi_zip <- file.path(
  work_dir,
  "Oceanisphaera_NCBI_genomes.zip"
)

ncbi_extract_dir <- file.path(
  work_dir,
  "Oceanisphaera_NCBI_package"
)

if (
  force_redownload ||
  !file.exists(ncbi_zip)
) {

  if (file.exists(ncbi_zip)) {
    unlink(ncbi_zip)
  }

  cat("\n============================================================\n")
  cat("DOWNLOADING OCEANISPHAERA GENOMES FROM NCBI\n")
  cat("============================================================\n")

  datasets_args <- c(
    "download",
    "genome",
    "taxon",
    taxon_name,
    "--include",
    "genome",
    "--mag",
    "exclude",
    "--exclude-atypical",
    "--exclude-multi-isolate",
    "--filename",
    ncbi_zip,
    "--no-progressbar"
  )

  status <- system2(
    datasets_bin,
    args = datasets_args
  )

  if (!identical(status, 0L)) {
    stop(
      "NCBI genome download failed with exit status ",
      status,
      "."
    )
  }
}

if (!file.exists(ncbi_zip)) {
  stop("NCBI genome package was not created.")
}

cat(
  "NCBI package size:",
  round(file.info(ncbi_zip)$size / 1024^2, 2),
  "MB\n"
)

# Fresh extraction if required.
if (
  force_redownload &&
  dir.exists(ncbi_extract_dir)
) {
  unlink(
    ncbi_extract_dir,
    recursive = TRUE,
    force = TRUE
  )
}

if (!dir.exists(ncbi_extract_dir)) {

  dir.create(
    ncbi_extract_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  unzip(
    ncbi_zip,
    exdir = ncbi_extract_dir
  )
}


############################
# 7. EXTRACT NCBI METADATA
############################

metadata_tsv <- file.path(
  work_dir,
  "Oceanisphaera_NCBI_metadata.tsv"
)

metadata_fields <- paste(
  c(
    "accession",
    "organism-name",
    "organism-infraspecific-strain",
    "type_material-label",
    "type_material-display_text",
    "checkm-completeness",
    "checkm-contamination"
  ),
  collapse = ","
)

metadata_stderr <- file.path(
  work_dir,
  "dataformat_metadata_stderr.txt"
)

metadata_status <- system2(
  dataformat_bin,
  args = c(
    "tsv",
    "genome",
    "--package",
    ncbi_zip,
    "--fields",
    metadata_fields
  ),
  stdout = metadata_tsv,
  stderr = metadata_stderr
)

# system2 may return either status or output depending on platform.
if (!file.exists(metadata_tsv) ||
    file.info(metadata_tsv)$size == 0) {
  warning(
    "NCBI metadata table could not be created. ",
    "ANI will still run, but organism/type labels may be unavailable."
  )
  ncbi_metadata <- data.frame()
} else {

  ncbi_metadata <- read.delim(
    metadata_tsv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write.csv(
    ncbi_metadata,
    file.path(
      analysis_dir,
      "Oceanisphaera_NCBI_metadata.csv"
    ),
    row.names = FALSE
  )

  cat(
    "\nNCBI metadata rows:",
    nrow(ncbi_metadata),
    "\n"
  )

  cat("NCBI metadata columns:\n")
  cat(
    paste0(
      "  - ",
      colnames(ncbi_metadata)
    ),
    sep = "\n"
  )
  cat("\n")
}


############################
# 8. FIND AND STAGE REFERENCE GENOMES
############################

reference_fastas <- list.files(
  ncbi_extract_dir,
  pattern = "\\.(fna|fa|fasta)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

# Retain genomic FASTA files; exclude README-like or non-genome files.
reference_fastas <- reference_fastas[
  grepl(
    "genomic\\.fna$|\\.fna$|\\.fasta$|\\.fa$",
    basename(reference_fastas),
    ignore.case = TRUE
  )
]

if (length(reference_fastas) == 0) {
  stop(
    "No genome FASTA files were found after extracting the NCBI package."
  )
}

get_accession_from_path <- function(path) {

  p <- normalizePath(
    path,
    winslash = "/",
    mustWork = TRUE
  )

  pieces <- strsplit(
    p,
    "/",
    fixed = TRUE
  )[[1]]

  hit <- grep(
    "^GC[AF]_[0-9]+\\.[0-9]+$",
    pieces,
    value = TRUE
  )

  if (length(hit) > 0) {
    return(hit[length(hit)])
  }

  b <- basename(path)

  m <- regexpr(
    "GC[AF]_[0-9]+\\.[0-9]+",
    b
  )

  if (m[1] > 0) {
    return(regmatches(b, m))
  }

  paste0(
    "REF_",
    sprintf(
      "%04d",
      match(
        path,
        reference_fastas
      )
    )
  )
}

reference_map <- data.frame(
  source_path = reference_fastas,
  accession = vapply(
    reference_fastas,
    get_accession_from_path,
    FUN.VALUE = character(1)
  ),
  stringsAsFactors = FALSE
)

# Keep a single genome file per accession.
reference_map <- reference_map[
  !duplicated(reference_map$accession),
  ,
  drop = FALSE
]

# Clear old staged references so stale files do not enter the analysis.
old_ref_files <- list.files(
  refs_stage_dir,
  full.names = TRUE
)

if (length(old_ref_files) > 0) {
  unlink(old_ref_files, force = TRUE)
}

reference_map$staged_path <- file.path(
  refs_stage_dir,
  paste0(
    reference_map$accession,
    ".fna"
  )
)

copy_ok <- mapply(
  FUN = function(from, to) {
    file.copy(
      from,
      to,
      overwrite = TRUE
    )
  },
  from = reference_map$source_path,
  to = reference_map$staged_path
)

if (!all(copy_ok)) {
  stop("At least one reference genome could not be staged.")
}

reference_map$staged_path <- normalizePath(
  reference_map$staged_path,
  winslash = "/",
  mustWork = TRUE
)

reference_map$wsl_path <- vapply(
  reference_map$staged_path,
  windows_to_wsl,
  FUN.VALUE = character(1)
)

write.csv(
  reference_map,
  file.path(
    analysis_dir,
    "Oceanisphaera_reference_genome_files.csv"
  ),
  row.names = FALSE
)

cat(
  "\nOceanisphaera reference genomes staged:",
  nrow(reference_map),
  "\n"
)


############################
# 9. CREATE FASTANI REFERENCE LIST
############################

ref_list_windows <- file.path(
  work_dir,
  "Oceanisphaera_reference_list.txt"
)

writeLines(
  reference_map$wsl_path,
  ref_list_windows
)

wsl_ref_list <- windows_to_wsl(
  ref_list_windows
)

forward_out_windows <- file.path(
  work_dir,
  "FM48a_vs_Oceanisphaera_fastANI.tsv"
)

reverse_out_windows <- file.path(
  work_dir,
  "Oceanisphaera_vs_FM48a_fastANI.tsv"
)

wsl_forward_out <- windows_to_wsl(
  forward_out_windows
)

wsl_reverse_out <- windows_to_wsl(
  reverse_out_windows
)


############################
# 10. RUN FASTANI: FM48a -> REFERENCES
############################

cat("\n============================================================\n")
cat("FASTANI: FM48a -> OCEANISPHAERA REFERENCES\n")
cat("============================================================\n")

if (file.exists(forward_out_windows)) {
  unlink(forward_out_windows)
}

forward_args <- c(
  "-d",
  wsl_distro,
  "--",
  "fastANI",
  "-q",
  wsl_query,
  "--rl",
  wsl_ref_list,
  "-o",
  wsl_forward_out,
  "-t",
  as.character(cpus)
)

forward_log <- system2(
  wsl_bin,
  args = forward_args,
  stdout = TRUE,
  stderr = TRUE
)

cat(
  paste(
    forward_log,
    collapse = "\n"
  ),
  "\n"
)

if (!file.exists(forward_out_windows)) {
  stop(
    "FastANI forward output was not created.\n",
    "See the console output above."
  )
}


############################
# 11. RUN FASTANI: REFERENCES -> FM48a
############################

cat("\n============================================================\n")
cat("FASTANI: OCEANISPHAERA REFERENCES -> FM48a\n")
cat("============================================================\n")

if (file.exists(reverse_out_windows)) {
  unlink(reverse_out_windows)
}

reverse_args <- c(
  "-d",
  wsl_distro,
  "--",
  "fastANI",
  "--ql",
  wsl_ref_list,
  "-r",
  wsl_query,
  "-o",
  wsl_reverse_out,
  "-t",
  as.character(cpus)
)

reverse_log <- system2(
  wsl_bin,
  args = reverse_args,
  stdout = TRUE,
  stderr = TRUE
)

cat(
  paste(
    reverse_log,
    collapse = "\n"
  ),
  "\n"
)

if (!file.exists(reverse_out_windows)) {
  stop(
    "FastANI reverse output was not created.\n",
    "See the console output above."
  )
}


############################
# 12. READ FASTANI OUTPUT
############################

read_fastani <- function(path) {

  if (
    !file.exists(path) ||
    file.info(path)$size == 0
  ) {
    return(data.frame())
  }

  x <- read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(x) < 5) {
    stop(
      "Unexpected FastANI output format in:\n",
      path
    )
  }

  x <- x[, 1:5, drop = FALSE]

  colnames(x) <- c(
    "query_path",
    "reference_path",
    "ANI",
    "mapped_fragments",
    "total_fragments"
  )

  x$alignment_fraction <-
    x$mapped_fragments /
    x$total_fragments

  x
}

forward <- read_fastani(
  forward_out_windows
)

reverse <- read_fastani(
  reverse_out_windows
)

if (
  nrow(forward) == 0 &&
  nrow(reverse) == 0
) {
  stop(
    "FastANI returned no reportable Oceanisphaera matches.\n",
    "This can occur if all genomes are too divergent or mapping is insufficient."
  )
}


############################
# 13. ACCESSION FROM FASTANI PATH
############################

accession_from_fastani_path <- function(path) {

  b <- basename(
    gsub(
      "\\\\",
      "/",
      path
    )
  )

  m <- regexpr(
    "GC[AF]_[0-9]+\\.[0-9]+",
    b
  )

  if (m[1] > 0) {
    return(regmatches(b, m))
  }

  tools::file_path_sans_ext(b)
}


############################
# 14. MERGE RECIPROCAL FASTANI RESULTS
############################

if (nrow(forward) > 0) {

  forward$accession <- vapply(
    forward$reference_path,
    accession_from_fastani_path,
    FUN.VALUE = character(1)
  )

  fwd <- forward[
    ,
    c(
      "accession",
      "ANI",
      "alignment_fraction",
      "mapped_fragments",
      "total_fragments"
    ),
    drop = FALSE
  ]

  names(fwd)[-1] <- paste0(
    names(fwd)[-1],
    "_FM48a_query"
  )

} else {

  fwd <- data.frame(
    accession = character(),
    stringsAsFactors = FALSE
  )
}

if (nrow(reverse) > 0) {

  reverse$accession <- vapply(
    reverse$query_path,
    accession_from_fastani_path,
    FUN.VALUE = character(1)
  )

  rev <- reverse[
    ,
    c(
      "accession",
      "ANI",
      "alignment_fraction",
      "mapped_fragments",
      "total_fragments"
    ),
    drop = FALSE
  ]

  names(rev)[-1] <- paste0(
    names(rev)[-1],
    "_reference_query"
  )

} else {

  rev <- data.frame(
    accession = character(),
    stringsAsFactors = FALSE
  )
}

ani <- merge(
  fwd,
  rev,
  by = "accession",
  all = TRUE
)

ani_cols <- grep(
  "^ANI_",
  names(ani),
  value = TRUE
)

af_cols <- grep(
  "^alignment_fraction_",
  names(ani),
  value = TRUE
)

ani$ANI_mean_reciprocal <- apply(
  ani[, ani_cols, drop = FALSE],
  1,
  function(z) {
    z <- as.numeric(z)
    if (all(is.na(z))) return(NA_real_)
    mean(z, na.rm = TRUE)
  }
)

ani$ANI_min_reciprocal <- apply(
  ani[, ani_cols, drop = FALSE],
  1,
  function(z) {
    z <- as.numeric(z)
    if (all(is.na(z))) return(NA_real_)
    min(z, na.rm = TRUE)
  }
)

ani$alignment_fraction_min <- apply(
  ani[, af_cols, drop = FALSE],
  1,
  function(z) {
    z <- as.numeric(z)
    if (all(is.na(z))) return(NA_real_)
    min(z, na.rm = TRUE)
  }
)

ani$passes_95ANI_50AF_screen <-
  ani$ANI_mean_reciprocal >= ani_screen_threshold &
  ani$alignment_fraction_min >= af_screen_threshold


############################
# 15. CLEAN NCBI METADATA COLUMN NAMES
############################

if (nrow(ncbi_metadata) > 0) {

  find_metadata_col <- function(patterns) {

    for (pat in patterns) {

      z <- grep(
        pat,
        colnames(ncbi_metadata),
        ignore.case = TRUE,
        value = TRUE
      )

      if (length(z) > 0) {
        return(z[1])
      }
    }

    NA_character_
  }

  accession_col <- find_metadata_col(
    c(
      "^Assembly Accession$",
      "^accession$"
    )
  )

  organism_col <- find_metadata_col(
    c(
      "^Organism Name$",
      "organism.*name"
    )
  )

  strain_col <- find_metadata_col(
    c(
      "Infraspecific.*Strain",
      "strain"
    )
  )

  type_label_col <- find_metadata_col(
    c(
      "^Type Material Label$",
      "type.*label"
    )
  )

  type_text_col <- find_metadata_col(
    c(
      "^Type Material Display Text$",
      "type.*display"
    )
  )

  completeness_col <- find_metadata_col(
    c(
      "CheckM.*Completeness",
      "completeness"
    )
  )

  contamination_col <- find_metadata_col(
    c(
      "CheckM.*Contamination",
      "contamination"
    )
  )

  meta_clean <- data.frame(
    accession = if (!is.na(accession_col)) {
      ncbi_metadata[[accession_col]]
    } else {
      NA_character_
    },
    organism_name = if (!is.na(organism_col)) {
      ncbi_metadata[[organism_col]]
    } else {
      NA_character_
    },
    strain = if (!is.na(strain_col)) {
      ncbi_metadata[[strain_col]]
    } else {
      NA_character_
    },
    type_material_label = if (!is.na(type_label_col)) {
      ncbi_metadata[[type_label_col]]
    } else {
      NA_character_
    },
    type_material_display_text = if (!is.na(type_text_col)) {
      ncbi_metadata[[type_text_col]]
    } else {
      NA_character_
    },
    checkm_completeness = if (!is.na(completeness_col)) {
      suppressWarnings(
        as.numeric(
          ncbi_metadata[[completeness_col]]
        )
      )
    } else {
      NA_real_
    },
    checkm_contamination = if (!is.na(contamination_col)) {
      suppressWarnings(
        as.numeric(
          ncbi_metadata[[contamination_col]]
        )
      )
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )

  meta_clean <- meta_clean[
    !duplicated(meta_clean$accession),
    ,
    drop = FALSE
  ]

  ani <- merge(
    ani,
    meta_clean,
    by = "accession",
    all.x = TRUE,
    sort = FALSE
  )

} else {

  ani$organism_name <- NA_character_
  ani$strain <- NA_character_
  ani$type_material_label <- NA_character_
  ani$type_material_display_text <- NA_character_
  ani$checkm_completeness <- NA_real_
  ani$checkm_contamination <- NA_real_
}


############################
# 16. IDENTIFY TYPE-MATERIAL ROWS
############################

type_text_combined <- paste(
  ani$type_material_label,
  ani$type_material_display_text
)

ani$is_type_material <-
  !is.na(type_text_combined) &
  nzchar(
    trimws(
      gsub(
        "NA",
        "",
        type_text_combined,
        fixed = TRUE
      )
    )
  )

# More conservative fallback: detect common type-related wording.
ani$is_type_material <-
  ani$is_type_material |
  grepl(
    "type|neotype|reference strain",
    type_text_combined,
    ignore.case = TRUE
  )

ani <- ani[
  order(
    -ani$ANI_mean_reciprocal,
    -ani$alignment_fraction_min
  ),
  ,
  drop = FALSE
]

rownames(ani) <- NULL


############################
# 17. WRITE MAIN ANI TABLE
############################

write.csv(
  ani,
  file.path(
    analysis_dir,
    "FM48a_Oceanisphaera_ANI_ranked.csv"
  ),
  row.names = FALSE
)

type_hits <- ani[
  ani$is_type_material %in% TRUE,
  ,
  drop = FALSE
]

write.csv(
  type_hits,
  file.path(
    analysis_dir,
    "FM48a_Oceanisphaera_TYPE_MATERIAL_ANI.csv"
  ),
  row.names = FALSE
)


############################
# 18. CONSOLE RESULTS
############################

cat("\n============================================================\n")
cat("TOP OCEANISPHAERA ANI MATCHES\n")
cat("============================================================\n")

show_cols <- c(
  "organism_name",
  "strain",
  "accession",
  "ANI_mean_reciprocal",
  "alignment_fraction_min",
  "is_type_material",
  "type_material_display_text",
  "passes_95ANI_50AF_screen"
)

show_cols <- intersect(
  show_cols,
  names(ani)
)

print(
  head(
    ani[, show_cols, drop = FALSE],
    20
  ),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("BEST TYPE-MATERIAL MATCHES\n")
cat("============================================================\n")

if (nrow(type_hits) > 0) {

  print(
    head(
      type_hits[, show_cols, drop = FALSE],
      20
    ),
    row.names = FALSE
  )

} else {

  cat(
    "No type-material labels were recovered in the NCBI metadata table.\n"
  )
}


############################
# 19. ANI FIGURE
############################

top_n <- min(
  20L,
  nrow(ani)
)

if (top_n > 0) {

  plot_dat <- ani[
    seq_len(top_n),
    ,
    drop = FALSE
  ]

  labels <- ifelse(
    !is.na(plot_dat$organism_name) &
      nzchar(plot_dat$organism_name),
    paste0(
      plot_dat$organism_name,
      ifelse(
        !is.na(plot_dat$strain) &
          nzchar(plot_dat$strain),
        paste0(" [", plot_dat$strain, "]"),
        ""
      ),
      "\n",
      plot_dat$accession
    ),
    plot_dat$accession
  )

  png(
    file.path(
      analysis_dir,
      "FM48a_Oceanisphaera_ANI_top20.png"
    ),
    width = 2200,
    height = 1500,
    res = 180
  )

  par(
    mar = c(5, 20, 4, 2)
  )

  barplot(
    rev(plot_dat$ANI_mean_reciprocal),
    names.arg = rev(labels),
    horiz = TRUE,
    las = 1,
    xlab = "Mean reciprocal ANI (%)",
    main = "FM48a vs Oceanisphaera genomes",
    xlim = c(
      max(
        0,
        floor(
          min(
            plot_dat$ANI_mean_reciprocal,
            na.rm = TRUE
          )
        ) - 2
      ),
      100
    ),
    cex.names = 0.70
  )

  abline(
    v = ani_screen_threshold,
    lty = 2
  )

  dev.off()
}


############################
# 20. FINAL SUMMARY
############################

best <- ani[1, , drop = FALSE]

best_type <- if (nrow(type_hits) > 0) {
  type_hits[1, , drop = FALSE]
} else {
  NULL
}

summary_lines <- c(
  "FM48a GENOME TAXONOMY — FASTANI SUMMARY",
  "=======================================",
  "",
  paste0(
    "FM48a assembly size: ",
    format(assembly_size, big.mark = ","),
    " bp"
  ),
  paste0(
    "FM48a contigs: ",
    n_contigs
  ),
  paste0(
    "FM48a N50: ",
    format(n50, big.mark = ","),
    " bp"
  ),
  paste0(
    "FM48a GC: ",
    sprintf("%.3f%%", gc_percent)
  ),
  "",
  paste0(
    "Oceanisphaera genomes downloaded/staged: ",
    nrow(reference_map)
  ),
  paste0(
    "FastANI reportable comparisons: ",
    nrow(ani)
  ),
  "",
  "Best overall match:",
  paste0(
    "  Organism: ",
    ifelse(
      is.na(best$organism_name),
      "<metadata unavailable>",
      best$organism_name
    )
  ),
  paste0(
    "  Strain: ",
    ifelse(
      is.na(best$strain),
      "<not reported>",
      best$strain
    )
  ),
  paste0(
    "  Accession: ",
    best$accession
  ),
  paste0(
    "  Mean reciprocal ANI: ",
    sprintf(
      "%.4f%%",
      best$ANI_mean_reciprocal
    )
  ),
  paste0(
    "  Minimum reciprocal alignment fraction: ",
    sprintf(
      "%.4f",
      best$alignment_fraction_min
    )
  ),
  paste0(
    "  Type material: ",
    best$is_type_material
  ),
  paste0(
    "  Passes 95% ANI + 0.50 AF screen: ",
    best$passes_95ANI_50AF_screen
  ),
  ""
)

if (!is.null(best_type)) {

  summary_lines <- c(
    summary_lines,
    "Best type-material match:",
    paste0(
      "  Organism: ",
      ifelse(
        is.na(best_type$organism_name),
        "<metadata unavailable>",
        best_type$organism_name
      )
    ),
    paste0(
      "  Strain: ",
      ifelse(
        is.na(best_type$strain),
        "<not reported>",
        best_type$strain
      )
    ),
    paste0(
      "  Accession: ",
      best_type$accession
    ),
    paste0(
      "  Mean reciprocal ANI: ",
      sprintf(
        "%.4f%%",
        best_type$ANI_mean_reciprocal
      )
    ),
    paste0(
      "  Minimum reciprocal alignment fraction: ",
      sprintf(
        "%.4f",
        best_type$alignment_fraction_min
      )
    ),
    paste0(
      "  Type-material annotation: ",
      best_type$type_material_display_text
    ),
    paste0(
      "  Passes 95% ANI + 0.50 AF screen: ",
      best_type$passes_95ANI_50AF_screen
    ),
    ""
  )
}

summary_lines <- c(
  summary_lines,
  "Interpretation:",
  paste0(
    "  ANI >= ",
    ani_screen_threshold,
    "% and reciprocal alignment fraction >= ",
    af_screen_threshold,
    " are used here as a species-level screening criterion."
  ),
  "  Species naming should prioritize comparison to type-material genomes.",
  "  A later GTDB-Tk/phylogenomic analysis can provide an independent genome-tree placement.",
  "",
  "Main output files:",
  "  FM48a_Oceanisphaera_ANI_ranked.csv",
  "  FM48a_Oceanisphaera_TYPE_MATERIAL_ANI.csv",
  "  FM48a_Oceanisphaera_ANI_top20.png",
  "  Oceanisphaera_NCBI_metadata.csv",
  "  FM48a_FASTANI_SUMMARY.txt"
)

writeLines(
  summary_lines,
  file.path(
    analysis_dir,
    "FM48a_FASTANI_SUMMARY.txt"
  )
)

cat("\n============================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================\n")
cat(
  paste(
    summary_lines,
    collapse = "\n"
  ),
  "\n"
)

cat("\nResults are here:\n  ")
cat(
  normalizePath(
    analysis_dir,
    winslash = "/",
    mustWork = TRUE
  )
)
cat("\n\nDONE.\n")
