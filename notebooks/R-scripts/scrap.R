# First R script
library(tidyverse)
library(reticulate)
library(magrittr)
library(glue)
library(seqinr)
library(future)
library(future.batchtools)
library(tictoc)
library(listenv)

tmpdir <- "/central/scratch/jbok/tmp"
homedir <- "/central/groups/MazmanianLab/joeB"
wkdir <- glue(
  "{homedir}/",
  "Microbiota-Immunomodulation/Microbiota-Immunomodulation"
)

src_dir <- glue("{wkdir}/notebooks")
source(glue("{src_dir}/R-scripts/helpers_general.R"))
source(glue("{src_dir}/R-scripts/helpers_sequence-screens.R"))

# reticulate::use_condaenv(condaenv = "pdmbsR", required = TRUE)


#______________________________________________________________________________
# GGET Blast ----
library(batchtools)
blastr <- function(sequences_aa,
                   ensembl_id,
                   db = "refseq_protein",
                   wkdir = wkdir,
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


# According 

# Collecting Human Immune proteome data
IDB_IRIS <-
  read.delim(
    glue("{wkdir}/data/input/InnateDB_metadata/",
    "InnateDB_genes_Immunogenetic-Related-Information-Source.txt"
    ),
    stringsAsFactors = FALSE, header = TRUE
  )

IDB_IRIS %>% colnames()
IDB_IRIS$ensembl %>% length()


gget <- import("gget")
# proteins_of_interest <- c("interleukin", "interferons", "TGFB", "TNFA")
keyname_search <- gget$search(IDB_IRIS$ensembl, "homo_sapiens")

keyname_hits <- keyname_search %>% filter(biotype == "protein_coding") %>% 
  filter(gene_name %in% c("IL1RN", "IL4", "IL6", "IL10", "IL11", "IL13"))

# collect amino acid sequences for each protein using the ensembl_id
gget_seq <- function(ensembl_id,  amino_acid=TRUE){
    seq_results <- gget$seq(ensembl_id, translate = amino_acid)
    return(seq_results[2])
}
keyname_hits %<>% mutate(sequences_aa =  map_chr(ensembl_id, gget_seq))
keyname_hits

saveRDS(keyname_hits, glue("{wkdir}/data/interim/tmp/{Sys.Date()}_gget_proteins-of-interest.rds"))





#______________________________________________________________________________
# Alignment Algorithms ----
# Obtaining 

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



# Selecting 

/central/groups/MazmanianLab/joeB/Downloads/protein_catalogs/clustered_catalogs/merged/2023-02-13_catalog_MMSeq2-95_rep_seq.fasta





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





# Read in metadata, select accession IDs for MAGs which are contamination and completeion criteria
#' Consistent with the definition of Medium Quality MAGs in the paper https://www.nature.com/articles/nbt.3893/tables/1
#' We will filter MAGs with a completion <= 50% and contamination < 10%

# hadza_mag_metadata <-



hadza_metagenomes_dir <- glue("/central/scratch/jbok/PRJEB49206_HADZA")
hadza_meta <- 
  c("prokaryotes", "viruses", "eukaryotes") %>%
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


hadza_meta$prokaryotes$assembly_acc %>% unique()  %>% length()
hadza_meta$prokaryotes$mag_sample_acc %>% unique()  %>% length()
hadza_meta$prokaryotes$mag_bin_acc %>% unique()  %>% length()

col <- "mag_bin_acc"
hadza_meta$prokaryotes[[col]] %>% length()
hadza_meta$prokaryotes[[col]] %>% unique() %>% length()



# catalog_paths <- list(
#   "UHGP" = glue("{protein_catalogs}/UHGP_v2.0.1/uhgp-95/uhgp-95.faa"),
#   "ELGG" = glue("{protein_catalogs}/ELGP_95.faa"),
#   "KIJ" = glue("{protein_catalogs}/KIJ_CD-HIT-100_Proteins.faa"),
#   "Hadza" = glue("{protein_catalogs}/hadza.fasta"),
#   "RUMMETA" = glue("{protein_catalogs}/RGMGC.geneSet.faa"),
#   "RUGS" = glue("{protein_catalogs}/cow-rumen-v1.0/protein_catalogue-95/protein_catalogue-95.faa"),
#   "MGV" = glue("{protein_catalogs}/mgv_proteins.faa"),
#   "VEuPathDB" = glue("{protein_catalogs}/VEuPathDB_v61.fasta"),
#   "Wormbase" = glue("{protein_catalogs}/wormbase-v17.0.fasta")
# )

# slurm_shell_do_mmseq2_clust(
#   input_fasta = glue("{protein_catalogs}/TEST.fasta"),
#   output = glue("{protein_catalogs}/clustered_catalogs/TEST/TEST"),
# )

# mmseqs easy-linclust --cov-mode 1 -c 0.8 \
# --kmer-per-seq 200 \
# --min-seq-id 0.95 \
# --threads 4 TEST.fasta \
# /central/groups/MazmanianLab/joeB/Downloads/protein_catalogs/clustered_catalogs/TEST/TEST tmp



# glue(
#   "~/bbmap/bbrename.sh",
#   " in={}",
#   " out={}_renamed_TEST.fasta",
#   " prefix=CAT_{catalog}",
#   " addprefix=t",
#   " fixjunk=t"
# )


# ~/bbmap/bbrename.sh in=TEST.fasta out=renamed_TEST.fasta prefix=TEST addprefix=t fixjunk=t






