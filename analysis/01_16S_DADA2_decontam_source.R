# Haliclystus antarcticus — host-only 16S V4 DADA2 workflow
#

## ---------------------------- 1. Configuration ----------------------------
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

data_dir <- file.path(repo_dir, "data")
metadata_path <- file.path(repo_dir, "metadata", "haliclystus_samples_65.csv")
fastq_dir <- file.path(data_dir, "fastq")
silva_dir <- file.path(data_dir, "silva")
silva_train <- file.path(silva_dir, "silva_nr99_v138.1_train_set.fa.gz")
silva_species <- file.path(silva_dir, "silva_species_assignment_v138.1.fa.gz")

# The same relative output folder is used by scripts 01–06.
output_dir <- file.path(repo_dir, "results", "16S")
rds_dir <- file.path(output_dir, "rds")
tables_dir <- file.path(output_dir, "tables")
plots_dir <- file.path(output_dir, "plots")
filtered_dir <- file.path(output_dir, "filtered")

# Original V4 primers and retained preprocessing settings.
primer_fwd <- "GTGYCAGCMGCCGCGGTAA"
primer_rev <- "GGACTACNVGGGTWTCTAAT"
trim_left <- c(nchar(primer_fwd), nchar(primer_rev))
trunc_len <- c(220L, 180L)
max_ee <- c(1, 1)
trunc_q <- 5L
min_overlap <- 12L

# ASV and sample filtering.
length_window <- 8L
pre_chimera_min_total <- 10L
pre_chimera_min_prevalence <- 3L
max_pre_chimera_asvs <- 120000L
min_library_size <- 1000L
min_prevalence_fraction <- 0.05
contam_threshold <- 0.1
use_multithread <- FALSE
set.seed(123)

## ------------------------------ 2. Packages ------------------------------
needed <- c("dada2", "decontam", "phyloseq")
not_installed <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(not_installed) > 0) {
  stop("Install the required R/Bioconductor packages first: ",
       paste(not_installed, collapse = ", "))
}

for (p in c(metadata_path, silva_train, silva_species)) {
  if (!file.exists(p)) stop("Missing required input: ", p)
}
if (!dir.exists(fastq_dir)) stop("FASTQ directory not found: ", fastq_dir,
                                    "\nCreate data/fastq/ in this repository and place the paired FASTQ files there.")
for (p in c(rds_dir, tables_dir, plots_dir, filtered_dir)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}
write_csv <- function(x, filename) {
  utils::write.csv(x, file.path(tables_dir, filename), row.names = FALSE, na = "")
}
report <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""))

## ------------------------ 3. Strict sample selection ------------------------
metadata <- utils::read.csv(
  metadata_path, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
required <- c("Sample_ID", "16S_ID", "Species", "Site", "Compartment",
              "DNA_concentration_ng_uL_original",
              "FASTQ_R1_expected", "FASTQ_R2_expected")
missing_cols <- setdiff(required, names(metadata))
if (length(missing_cols)) stop("Missing metadata fields: ", paste(missing_cols, collapse = ", "))

if (nrow(metadata) != 65L) {
  stop("The host-only metadata must contain exactly 65 rows (found ", nrow(metadata), ").")
}
for (col in c("Sample_ID", "16S_ID", "FASTQ_R1_expected", "FASTQ_R2_expected")) {
  if (anyNA(metadata[[col]]) || any(!nzchar(trimws(metadata[[col]]))) || anyDuplicated(metadata[[col]])) {
    stop("Invalid, missing or repeated entries in metadata field: ", col)
  }
}
if (!all(metadata$Species == "Haliclystus antarcticus")) {
  stop("The metadata must contain only Haliclystus antarcticus tissue samples.")
}
if (!all(metadata$Site %in% c("Byers Bay", "Risopatron")) ||
    !all(metadata$Compartment %in% c("gut", "peduncle", "tentacle"))) {
  stop("Unexpected collection site or anatomical compartment in host metadata.")
}

# Use the exact paired-end FASTQ basenames supplied in the metadata. This
# avoids scanning and processing unrelated FASTQ files within the directory.
expected_F <- paste0(metadata$`16S_ID`, "_R1_001.fastq.gz")
expected_R <- paste0(metadata$`16S_ID`, "_R2_001.fastq.gz")
if (!identical(metadata$FASTQ_R1_expected, expected_F) ||
    !identical(metadata$FASTQ_R2_expected, expected_R)) {
  stop("FASTQ filenames in metadata do not match the 16S_ID naming convention.")
}
files_F <- file.path(fastq_dir, metadata$FASTQ_R1_expected)
files_R <- file.path(fastq_dir, metadata$FASTQ_R2_expected)
missing_files <- c(files_F[!file.exists(files_F)], files_R[!file.exists(files_R)])
if (length(missing_files)) {
  stop("Missing ", length(missing_files), " expected FASTQ(s). First missing files: ",
       paste(basename(head(missing_files, 8L)), collapse = ", "))
}

metadata$Site_label <- metadata$Site
metadata$Compartment_label <- metadata$Compartment
# Retain category names expected by the archived figure-analysis scripts.
metadata$Site <- factor(ifelse(metadata$Site == "Byers Bay", "Baia_Byers", "Rizopatron"),
                        levels = c("Baia_Byers", "Rizopatron"))
metadata$Compartment <- factor(ifelse(metadata$Compartment == "peduncle", "stalk",
                                      metadata$Compartment),
                               levels = c("gut", "stalk", "tentacle"))
metadata$Sample_type <- "Haliclystus"
metadata$dna_conc_original <- suppressWarnings(
  as.numeric(metadata$DNA_concentration_ng_uL_original)
)

# Preserve the concentration substitution used in the original frequency
# analysis, recording each affected sample without editing the source CSV.
valid_conc <- is.finite(metadata$dna_conc_original) & metadata$dna_conc_original > 0
if (sum(valid_conc) < 3L) stop("Too few positive DNA concentrations for decontam.")
min_valid <- min(metadata$dna_conc_original[valid_conc])
metadata$dna_conc <- metadata$dna_conc_original
metadata$dna_conc[!valid_conc] <- min_valid / 2
write_csv(
  data.frame(
    Sample_ID = metadata$Sample_ID,
    `16S_ID` = metadata$`16S_ID`,
    original = metadata$DNA_concentration_ng_uL_original,
    used_ng_uL = metadata$dna_conc,
    substituted = !valid_conc,
    check.names = FALSE
  ),
  "dna_concentrations_used.csv"
)

sample_ids <- metadata$`16S_ID`
names(files_F) <- sample_ids
names(files_R) <- sample_ids
write_csv(data.frame(sample_id = sample_ids,
                     file_F = basename(files_F), file_R = basename(files_R)),
          "fastq_pairing_table.csv")
report("Selected ", length(sample_ids), " Haliclystus samples from metadata.")

## ------------------------- 4. Filter and denoise -------------------------
filtFs <- file.path(filtered_dir, paste0(sample_ids, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample_ids, "_R_filt.fastq.gz"))
names(filtFs) <- sample_ids
names(filtRs) <- sample_ids

filt_out <- dada2::filterAndTrim(
  fwd = files_F, filt = filtFs, rev = files_R, filt.rev = filtRs,
  truncLen = trunc_len, trimLeft = trim_left, maxN = 0,
  maxEE = max_ee, truncQ = trunc_q, rm.phix = TRUE,
  compress = TRUE, multithread = use_multithread,
  matchIDs = TRUE, verbose = TRUE
)
if (nrow(filt_out) != length(sample_ids)) stop("Unexpected number of filterAndTrim outputs.")

qc <- data.frame(sample_id = sample_ids,
                 input = as.numeric(filt_out[, 1]),
                 filtered = as.numeric(filt_out[, 2]))
write_csv(qc, "filter_and_trim_summary.csv")
saveRDS(filt_out, file.path(rds_dir, "filter_and_trim_summary.rds"))

# Samples with no usable reads cannot be denoised. Record exactly which drop.
retained <- qc$filtered > 0
if (any(!retained)) {
  report("No filtered reads for: ", paste(qc$sample_id[!retained], collapse = ", "))
}
if (sum(retained) < 3L) stop("Fewer than three host samples remain after filtering.")
sample_ids <- sample_ids[retained]
filtFs <- filtFs[retained]
filtRs <- filtRs[retained]

report("Learning error profiles for ", length(sample_ids), " samples")
errF <- dada2::learnErrors(filtFs, multithread = use_multithread, verbose = TRUE)
errR <- dada2::learnErrors(filtRs, multithread = use_multithread, verbose = TRUE)
saveRDS(errF, file.path(rds_dir, "errF.rds"))
saveRDS(errR, file.path(rds_dir, "errR.rds"))

# Inference is per sample, matching the archived processing strategy.
dadaFs <- setNames(lapply(filtFs, function(f) dada2::dada(f, err = errF,
                                                        multithread = use_multithread)), sample_ids)
dadaRs <- setNames(lapply(filtRs, function(f) dada2::dada(f, err = errR,
                                                        multithread = use_multithread)), sample_ids)
saveRDS(dadaFs, file.path(rds_dir, "dadaFs.rds"))
saveRDS(dadaRs, file.path(rds_dir, "dadaRs.rds"))

mergers <- dada2::mergePairs(dadaFs, unname(filtFs), dadaRs, unname(filtRs),
                             trimOverhang = TRUE, minOverlap = min_overlap, verbose = TRUE)
saveRDS(mergers, file.path(rds_dir, "mergers.rds"))
seqtab_raw <- dada2::makeSequenceTable(mergers)
if (!all(sample_ids %in% rownames(seqtab_raw))) {
  stop("Sequence table row names do not match expected sample identifiers.")
}
seqtab_raw <- seqtab_raw[sample_ids, , drop = FALSE]
saveRDS(seqtab_raw, file.path(rds_dir, "seqtab_raw.rds"))
if (ncol(seqtab_raw) == 0L) stop("No ASVs survived read-pair merging.")

# Same amplicon-length and pre-chimera filters as the archived script.
asv_lengths <- nchar(colnames(seqtab_raw))
main_length <- as.integer(names(sort(table(asv_lengths), decreasing = TRUE)[1]))
keep_length <- abs(asv_lengths - main_length) <= length_window
seqtab_len <- seqtab_raw[, keep_length, drop = FALSE]
saveRDS(seqtab_len, file.path(rds_dir, "seqtab_length_filtered.rds"))

asv_total <- colSums(seqtab_len)
asv_prevalence <- colSums(seqtab_len > 0)
keep_asv <- asv_total >= pre_chimera_min_total &
  asv_prevalence >= pre_chimera_min_prevalence
seqtab <- seqtab_len[, keep_asv, drop = FALSE]
if (ncol(seqtab) == 0L) stop("No ASVs survived the pre-chimera filter.")
if (ncol(seqtab) > max_pre_chimera_asvs) {
  stop("Unexpectedly large pre-chimera ASV table: ", ncol(seqtab))
}
saveRDS(seqtab, file.path(rds_dir, "seqtab.rds"))

report("Removing chimeras")
seqtab_nochim <- dada2::removeBimeraDenovo(
  seqtab, method = "consensus", multithread = use_multithread, verbose = TRUE
)
if (ncol(seqtab_nochim) == 0L) stop("No ASVs remain after chimera removal.")
saveRDS(seqtab_nochim, file.path(rds_dir, "seqtab_nochim.rds"))

getN <- function(x) sum(dada2::getUniques(x))
track <- data.frame(
  sample_id = sample_ids,
  input = qc$input[match(sample_ids, qc$sample_id)],
  filtered = qc$filtered[match(sample_ids, qc$sample_id)],
  denoisedF = vapply(dadaFs, getN, numeric(1)),
  denoisedR = vapply(dadaRs, getN, numeric(1)),
  merged = vapply(mergers, getN, numeric(1)),
  nonchim = rowSums(seqtab_nochim)
)
write_csv(track, "read_tracking.csv")
saveRDS(track, file.path(rds_dir, "read_tracking.rds"))

## ------------------- 5. Taxonomy and bacterial ASVs -------------------
report("Assigning SILVA taxonomy")
taxa <- dada2::assignTaxonomy(seqtab_nochim, silva_train,
                             multithread = use_multithread)
taxa <- dada2::addSpecies(taxa, silva_species)
saveRDS(taxa, file.path(rds_dir, "taxa_raw.rds"))
write_csv(data.frame(ASV = rownames(taxa), taxa, check.names = FALSE),
          "taxonomy_raw.csv")

keep_tax <- rep(TRUE, nrow(taxa))
if ("Kingdom" %in% colnames(taxa)) {
  keep_tax <- keep_tax & !is.na(taxa[, "Kingdom"]) & taxa[, "Kingdom"] == "Bacteria"
}
if ("Order" %in% colnames(taxa)) {
  keep_tax <- keep_tax & (is.na(taxa[, "Order"]) |
                            !grepl("Chloroplast", taxa[, "Order"], ignore.case = TRUE))
}
if ("Family" %in% colnames(taxa)) {
  keep_tax <- keep_tax & (is.na(taxa[, "Family"]) |
                            !grepl("Chloroplast|Mitochondria", taxa[, "Family"], ignore.case = TRUE))
}
seqtab_taxclean <- seqtab_nochim[, keep_tax, drop = FALSE]
taxa_taxclean <- taxa[colnames(seqtab_taxclean), , drop = FALSE]
if (ncol(seqtab_taxclean) == 0L) stop("No bacterial ASVs remain after taxonomy filtering.")
saveRDS(seqtab_taxclean, file.path(rds_dir, "seqtab_taxclean.rds"))
saveRDS(taxa_taxclean, file.path(rds_dir, "taxa_taxclean.rds"))

## ----------------- 6. Frequency-based decontamination -----------------

meta_matched <- metadata[match(rownames(seqtab_taxclean), metadata$`16S_ID`), , drop = FALSE]
if (anyNA(meta_matched$`16S_ID`)) stop("Metadata matching failed.")
rownames(meta_matched) <- meta_matched$`16S_ID`

contam_frequency <- decontam::isContaminant(
  seqtab_taxclean,
  method = "frequency",
  conc = meta_matched$dna_conc,
  threshold = contam_threshold
)
contam_tbl <- data.frame(
  ASV = colnames(seqtab_taxclean),
  contaminant_frequency = as.logical(contam_frequency$contaminant)
)
if (anyNA(contam_tbl$contaminant_frequency)) {
  stop("decontam returned NA classifications. Inspect DNA concentration inputs.")
}
# Keep the original reporting field for downstream compatibility.
contam_tbl$contaminant_any <- contam_tbl$contaminant_frequency
write_csv(contam_tbl, "decontam_results.csv")

seqtab_decontam <- seqtab_taxclean[,
  !contam_tbl$contaminant_frequency, drop = FALSE]
if (ncol(seqtab_decontam) == 0L) stop("No ASVs remain after decontamination.")
taxa_decontam <- taxa_taxclean[colnames(seqtab_decontam), , drop = FALSE]
saveRDS(seqtab_decontam, file.path(rds_dir, "seqtab_decontam.rds"))
saveRDS(taxa_decontam, file.path(rds_dir, "taxa_decontam.rds"))

## ----------------- 7. Final host sample/ASV filtering -----------------
meta_halo <- meta_matched
seqtab_halo <- seqtab_decontam[meta_halo$`16S_ID`, , drop = FALSE]
keep_library <- rowSums(seqtab_halo) >= min_library_size
write_csv(data.frame(sample_id = rownames(seqtab_halo),
                     reads_after_decontam = rowSums(seqtab_halo),
                     passed_library_filter = keep_library),
          "host_library_size_qc.csv")
seqtab_halo <- seqtab_halo[keep_library, , drop = FALSE]
meta_halo <- meta_halo[keep_library, , drop = FALSE]
if (nrow(seqtab_halo) < 3L) stop("Too few samples after library-size filtering.")
rownames(meta_halo) <- meta_halo$`16S_ID`

# Keep the original prevalence cutoff formula used in downstream analyses.
prevalence_cutoff <- max(1L, floor(nrow(seqtab_halo) * min_prevalence_fraction))
keep_prev <- colSums(seqtab_halo > 0) >= prevalence_cutoff
seqtab_halo_prev <- seqtab_halo[, keep_prev, drop = FALSE]
if (ncol(seqtab_halo_prev) == 0L) stop("No ASVs survive prevalence filtering.")
taxa_halo_prev <- taxa_decontam[colnames(seqtab_halo_prev), , drop = FALSE]

saveRDS(seqtab_halo, file.path(rds_dir, "seqtab_haloclystus_raw.rds"))
saveRDS(seqtab_halo_prev, file.path(rds_dir, "seqtab_haloclystus_prevfiltered.rds"))
saveRDS(meta_halo, file.path(rds_dir, "metadata_haloclystus.rds"))
saveRDS(taxa_halo_prev, file.path(rds_dir, "taxa_haloclystus_prevfiltered.rds"))

ps_halo <- phyloseq::phyloseq(
  phyloseq::otu_table(seqtab_halo_prev, taxa_are_rows = FALSE),
  phyloseq::tax_table(as.matrix(taxa_halo_prev)),
  phyloseq::sample_data(meta_halo)
)
saveRDS(ps_halo, file.path(rds_dir, "phyloseq_haloclystus.rds"))

write_csv(data.frame(
  metric = c("samples_in_metadata", "samples_processed_after_read_filter",
             "samples_after_library_filter", "asvs_after_chimera_filter",
             "frequency_contaminants", "asvs_after_prevalence_filter"),
  value = c(nrow(metadata), length(sample_ids), nrow(seqtab_halo),
            ncol(seqtab_nochim), sum(contam_tbl$contaminant_frequency),
            ncol(seqtab_halo_prev))
), "host_pipeline_summary.csv")
report("Complete: ", nrow(seqtab_halo_prev), " host samples; ",
       ncol(seqtab_halo_prev), " retained ASVs.")
report("Run 02_phyloseq_beta_source.R next to produce the phylogenetic tree object.")
