# First R script
library(tidyverse)
library(reticulate)
library(magrittr)
library(glue)
library(seqinr)
library(future)
library(future.batchtools)

tmpdir <- "/central/scratch/jbok/tmp"
homedir <- "/central/groups/MazmanianLab/joeB"
wkdir <- glue(
  "{homedir}/",
  "Microbiota-Immunomodulation/Microbiota-Immunomodulation"
)

src_dir <- glue("{wkdir}/notebooks")
source(glue("{src_dir}/R-scripts/helpers_general.R"))
source(glue("{src_dir}/R-scripts/helpers_sequence-screens.R"))

reticulate::use_condaenv(condaenv = "pdmbsR", required = TRUE)


#______________________________________________________________________________
# GGET Blast ----
library(batchtools)
blastr <- function(sequences_aa,
                   ensembl_id,
                   db = "refseq_protein",
                   wkdir = "/central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation",
                   output_limit = 10000) {
  require(reticulate)
  require(glue)
  reticulate::use_condaenv(condaenv = "pdmbsR", required = TRUE)
  gget <- reticulate::import("gget")
  blast_results <- gget$blast(sequences_aa, database = db, limit = output_limit)
  saveRDS(
    blast_results,
    glue("{wkdir}/data/interim/blast_results/{ensembl_id}.rds")
  )
  return(blast_results)
}

# configure registry ----
breg <- makeRegistry(
  file.dir = glue(
    "{wkdir}/.cluster_runs/",
    "{Sys.Date()}_BLAST_{sample(10000:99999, 1)}/"
  ),
  seed = 42
)
breg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
  template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
  scheduler.latency = 2,
  fs.latency = 65
)
# Submit Jobs ----
jobs <- batchMap(
  fun = blastr,
  args = dplyr::select(keyname_hits, sequences_aa, ensembl_id),
  reg = breg
)
setJobNames(jobs, paste0("BLAST_", keyname_hits$ensembl_id), reg = breg)
submitJobs(jobs,
  resources = list(walltime = 10800,
    memory = 1024,
    ncpus = 8,
    max.concurrent.jobs = 9999)
)

# getJobNames()$job.id
# waitForJobs(sleep = 5)
# tst <- loadResult(1)
# tst %>% glimpse




#_______________________________________________________________________________
# Chunking protein catalog fasta files ----

# INITATE FUTURE-SLURM BACKEND FOR PARALLEL PROCESSING
future::plan(
  future.batchtools::batchtools_slurm,
  template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
  resources = list(
    name = "chunk-fasta-catalog",
    memory = "50G",
    ncpus = 2,
    walltime = 36000
  )
)

# divide original fasta into 10 smaller fasta files, roughly even in size
chunk_fasta_file <- function(fasta_dir,
                             fasta_path,
                             working_dir = wkdir) {
  require(seqinr)
  require(fs)
  require(glue)
  require(purrr)
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  fasta_name <- fs::path_ext_remove(basename(fasta_path))
  f <- read.fasta(fasta_path)
  fc <- chunk_func(f, 10)
  for (i in 1:length(fc)) {
    write.fasta(
      purrr::map(getSequence(fc[[i]]), toupper),
      getAnnot(fc[[i]]),
      glue("{fasta_dir}/{fasta_name}_BIN_{i}.fasta")
    )
  }
}

protein_cat_dir <- glue("{homedir}/Downloads/protein_catalogs/tmp_mim_catalog")
fout %<-% {
  chunk_fasta_file(
    protein_cat_dir,
    glue("{protein_cat_dir}/uhgp-95_v1.fasta")
  )
}


#______________________________________________________________________________
# Alignment Algorithums ----

keyname_hits <- readRDS(glue("{wkdir}/data/interim/tmp/2023-01-27_gget_proteins-of-interest.rds"))
protein_catalog_path <- "/central/groups/MazmanianLab/joeB/Downloads/protein_catalogs/tmp_mim_catalog"
chunked_catalog_paths <- list.files(protein_catalog_path, full.names = TRUE) %>% keep(grepl("BIN", .))

# create a datatable mapping each ensembl_id to each chunked catalog
binned_eids <- data.frame(
  "ensembl_id" = rep(keyname_hits$ensembl_id, each = 10),
  "refdb_path" = rep(chunked_catalog_paths, times = length(keyname_hits$ensembl_id))
)

alignment_df_list <- bind_rows(
  keyname_hits %>% 
    dplyr::select(ensembl_id, sequences_aa) %>%
    mutate(method = "needle"),
  keyname_hits %>% 
    dplyr::select(ensembl_id, sequences_aa) %>%
    mutate(method = "water")
    ) %>% 
    full_join(binned_eids, by = "ensembl_id")

alignment_df_list %>% glimpse()
View(alignment_df_list)


#______________________________________________________________________________
# EXECUTE ALIGNMENT JOBS ----

# configure registry ----
cluster_run <- glue("{Sys.Date()}_EMBOSS-Alignment_ID-{rand_string()}/")
message("\n\nRUNNING:  ", cluster_run, "\n")
breg <- makeRegistry(
  file.dir = glue(
    "{wkdir}/.cluster_runs/",
    cluster_run
  ),
  seed = 42
)
breg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
  template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
  scheduler.latency = 0.5,
  fs.latency = 65
)
# Submit Jobs ----
jobs <- batchMap(
  fun = align_fasta_sequences,
  args = alignment_df_list,
  reg = breg
)
setJobNames(jobs,
  paste("EMBOSS",
  alignment_df_list$ensembl_id, 
  alignment_df_list$method, 
  sep = "_"),
  reg = breg
)
getJobNames(jobs, reg = breg)
submitJobs(jobs,
  resources = list(
    walltime = 60000,
    memory = 4096,
    ncpus = 8,
    max.concurrent.jobs = 9999
  )
)

waitForJobs(sleep = 3)


#______________________________________________________________________________
# Aggregate Alignment Files ----


file_dir <- glue("{wkdir}/data/interim/emboss_alignments")
file_paths <- list.files(file_dir, pattern = "rds$", full.names = TRUE)
needle_files <- file_paths %>%
  keep(grepl("ENSG00000169194_uhgp-95_v1_BIN_.*needle.rds", .))
water_files <- file_paths %>%
  keep(grepl("ENSG00000169194_uhgp-95_v1_BIN_.*water.rds", .))

# INITATE FUTURE-SLURM BACKEND FOR PARALLEL PROCESSING
library(future.batchtools)
future::plan(
  future.batchtools::batchtools_slurm,
  template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
  resources = list(
    name = "aggregate-rds",
    memory = "10G",
    ncpus = 1,
    walltime = 16384
  )
)
aggregate_alignment_files <- function(input_files,
                                      outputdir,
                                      id,
                                      working_dir = wkdir) {
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  message(glue("{Sys.time()}: Compiling {length(input_files)} files"))
  compiled_data <- load_rds_files(input_files)
  message(glue("{Sys.time()}: Saving compiled data to {outputdir}/{id}"))
  saveRDS(compiled_data, glue("{outputdir}/{id}"))
}

keyname_hits <- readRDS(
  glue("{wkdir}/data/interim/tmp/2023-01-27_gget_proteins-of-interest.rds")
)
eids <- keyname_hits$ensembl_id
file_dir <- glue("{wkdir}/data/interim/emboss_alignments")
file_paths <- list.files(file_dir, pattern = "rds$", full.names = TRUE)
output_dir <- glue("{wkdir}/data/processed/emboss_alignments")

for (eid in eids) {
  print(eid)
  for (method in c("needle", "water")) {
    name_out <- glue("{Sys.Date()}_{eid}_uhgp-95_v1.{method}.rds")
    files <- file_paths %>%
      keep(grepl(glue("{eid}_uhgp-95_v1_BIN_.*{method}.rds"), .))
    agg_files_needle %<-% {
      aggregate_alignment_files(
        files,
        output_dir,
        name_out)
    }
  }
}


# IL-10
df_aligned <-
  readRDS(glue("{output_dir}/2023-02-05_ENSG00000136634_uhgp-95_v1.water.rds"))
df_aligned %>% glimpse()
# calculate summary stats 

df_aligned %>%
  dplyr::summarize(
    alignment_pairs = n(),
    mean_score = mean(Score),
    min_score = min(Score),
    max_score = max(Score),
    mean_percent_identity = mean(percent_identity),
    min_percent_identity = min(percent_identity),
    max_percent_identity = max(percent_identity),
    mean_percent_similarity = mean(percent_similarity),
    min_percent_similarity = min(percent_similarity),
    max_percent_similarity = max(percent_similarity)
  )


df_aligned_l50 <- df_aligned %>%
  filter(Length > 50)
df_aligned_l50 %>% glimpse()

df_aligned_l10 %>%
  ggplot(aes(Score)) +
  geom_histogram(bins = 100)

df_aligned_s150 %>%
  ggplot(aes(percent_similarity, percent_identity)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm") +
  theme_bw()







#_______________________________________________________________________________
# Downloading protein catalogs ----





# Hadza ----
# Downloading Hadza Mag metadata
shell_do(glue(
  "mkdir -p {wkdir}/data/input/protein_catalog_metadata"
))
shell_do(glue(
  "wget -P {wkdir}/data/input/protein_catalog_metadata/",
  " https://www.biorxiv.org/content/biorxiv/early/2022/10/07/2022.03.30.486478/DC2/embed/media-2.xlsx",
))
shell_do(glue(
  "mv {wkdir}/data/input/protein_catalog_metadata/media-2.xlsx",
  " {wkdir}/data/input/protein_catalog_metadata/Hadza-MAG-metadata.xlsx"
))

# Read in metadata, select accession IDs for MAGs which are contamination and completeion criteria
#' Consistent with the definition of Medium Quality MAGs in the paper https://www.nature.com/articles/nbt.3893/tables/1
#' We will filter MAGs with a completion <= 50% and contamination < 10%

# hadza_mag_metadata <-

# hadza_meta <- list()
hadza_meta <- c("prokaryotes", "viruses", "eukaryotes") %>%
  purrr::set_names() %>%
  purrr::map(
  ~ readxl::read_excel(
    glue("{wkdir}/data/input/protein_catalog_metadata/Hadza-MAG-metadata.xlsx"),
    sheet = .
  )
)

#Explore quality metrics for each MAG type
hadza_meta$viruses$checkv_quality %>% table()
hadza_meta$prokaryotes$MAG_quality %>% table()
hadza_meta$eukaryotes$MAG_quality %>% table()

# Remove low-quality viruse MAGs
hadza_meta$viruses <- hadza_meta$viruses %>%
  filter(checkv_quality != "Low-quality")
hadza_meta$prokaryotes %>% glimpse()

col <- "mag_bin_acc"
hadza_meta$prokaryotes[[col]] %>% length()
hadza_meta$prokaryotes[[col]] %>% unique() %>% length()



# slurm_out <-  glue("{wkdir}/.cluster_runs/{Sys.Date()}_Kraken2_paired_RNASEQ")
# data_out <- glue("{pdmbs_dir}/workflow/RNASEQ/results/kraken2")
# dir.create(slurm_out, recursive = TRUE)
# dir.create(data_out, recursive = TRUE)


# for (sampleID in sampleIDs_rna) {
#   if (!file.exists(
#     file.path(data_out, glue("UHGG_reports/{sampleID}_UHGG_report.tsv"))
#   )) {
#     jobname <- glue("kraken2_{sampleID}")
#     # pause if there are more than 9,999 slurm jobs submitted
#     check_slurm_overload()
#     shell_do(
#       glue(
#         "sbatch",
#         " --job-name={jobname}",
#         " --output={slurm_out}/{jobname}.out",
#         " --error={slurm_out}/{jobname}.err",
#         " notebooks/shell_scripts/kraken2-paired.sh",
#         " -s {sampleID}",
#         " -r '/central/groups/MazmanianLab/joeB/PDMBS/workflow/RNASEQ/results/kraken2'",
#         " -f '/central/groups/MazmanianLab/joeB/PDMBS/workflow/RNASEQ/clean_reads'"
#       )
#     )
#     # to prevent memory alloc issues
#     Sys.sleep(0.02)
#   }
# }



# fastq-dl ERZ14760497 ENA



# enaDataGet -f fastq -d {outdir} {acc}

# enaGroupGet -f submitted SAMEA2591084
# conda install -c bioconda enabrowsertools
# conda install -c bioconda fastq-dl

# enaDataGet -f submitted ERR164409
# enaDataGet -f fasta -d /central/groups/MazmanianLab/joeB/Downloads/scrap ERZ4561419
# enaDataGet ERZ4561419
# /central/groups/MazmanianLab/joeB/Downloads/sratoolkit.3.0.1-centos_linux64/bin/fastq-dump ERZ4561419



# Hadza MAG directory
# /central/groups/MazmanianLab/joeB/Downloads/mag_library/PRJEB49206_HADZA