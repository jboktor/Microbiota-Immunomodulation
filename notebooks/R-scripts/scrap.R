
# First R script
library(reticulate)
library(tidyverse)
library(magrittr)
# renv::use_python(python= '/home/jboktor/miniconda3/envs/pdmbsR' )

reticulate::use_condaenv(condaenv = 'pdmbsR',  required = TRUE)
# reticulate::conda_install("pdmbsR", "gget", channel = "bioconda")


# What are the immune cell proteins of interest to screen? 

#' Maybe for expedience we can start off with a set list of cytokines, ie IL2, IL...
#' But for a more rigorous approach we can mine existing datasets to more appropriately answer this question
#' 
#' What type of dataset will provide us with the info we want? Ideally proteins that are implicated in active immune 
#' cell regulation are of interest, these proteins will be upgregulated and heavily involived in immune cell
#' perturbation assays.
#' An ideal dataset for this may include single cell transcriptional profiles of immune cells (PBMCs and gut tissues would be great)
#' with a microbial perturbation. This way we can explore the upgregulation/downregulation of genes and select
#' protein targets of interest from there. 
#' 






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






# needle -asequence test.fasta \
# -bsequence human-gut/v1.0/uhgp_catalogue/uhgp-95/uhgp-95.faa \
# -gapopen 10.0 \
# -gapextend 0.5 \
# -outfile 'test-ughp95v1.needle'

# needle -asequence test.fasta \
# -bsequence test-db.fasta \
# -gapopen 5.0 \
# -gapextend 0.5 \
# -datafile EBLOSUM62 \
# -aformat score \
# -outfile 'test-ughp95v1.needle'

# water -asequence test.fasta \
# -bsequence test-db.fasta \
# -gapopen 5.0 \
# -gapextend 0.5 \
# -datafile EBLOSUM62 \
# -outfile 'test-ughp95v1_srspair.water'




 parse_emboss_alignment_file <- function(file_path) {
  # Script to read a text file, select only lines that start with a hash, and parse the data into a data frame
  # Step 1: Read in the text file
  text_file <- readLines(con <- file(file_path))
  hash_lines <- grep("^#", text_file, value = TRUE)
  df_results <- tibble()
  col_vars <- c(
    "1", "2", "Matrix", "Gap_penalty", "Extend_penalty", 
    "Length", "Identity", "Similarity", "Gaps", "Score"
    )
  for (line in hash_lines) {
    # Split the line by ":" and extract the second element (the value)
    line_split <- strsplit(line, ":")[[1]]
    key <-  gsub("# ", "", line_split[1])
    val <- gsub(" ", "", line_split[2])

    if (key %in% col_vars){
      if (key == "1") {
        # initiate list
        new_row <- list()
        new_row[[key]] = val
      } else if (key == "Score") {
        # terminate list and bind to datafarme
        new_row[[key]] = val
        df_results %<>% bind_rows(new_row)
      } else {
        new_row[[key]] = val
      }
    }
  }
  df_results %<>% 
    mutate(
      percent_identity =  map_chr(Identity, 
          ~ gsub("[\\(\\%)]", "", regmatches(., gregexpr("\\(.*?\\)", .))[[1]])),
      percent_similarity = map_chr(Similarity, 
        ~ gsub("[\\(\\%)]", "", regmatches(., gregexpr("\\(.*?\\)", .))[[1]])),
      percent_gaps = map_chr(Gaps, 
        ~ gsub("[\\(\\%)]", "", regmatches(., gregexpr("\\(.*?\\)", .))[[1]]))
    ) %>% 
    mutate_at(vars(Gap_penalty, Extend_penalty, Length, Score, contains("percent_")),  as.numeric) %>% 
    rename(protein_1 = `1`, protein_2 = `2` )
  return(df_results)
 }


parse_emboss_alignment_file("/central/groups/MazmanianLab/joeB/Downloads/blastdb/test-ughp95v1_srspair.water")
parse_emboss_alignment_file("/central/groups/MazmanianLab/joeB/Downloads/blastdb/test-ughp95v1_srspair.needle")








# blastr <- function(sequences_aa,
#                    ensembl_id,
#                    db = "refseq_protein",
#                    wkdir = "/central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation",
#                    output_limit = 10000) {
#   require(reticulate)
#   require(glue)
#   reticulate::use_condaenv(condaenv = "pdmbsR", required = TRUE)
#   gget <- reticulate::import("gget")
#   blast_results <- gget$blast(sequences_aa, database = db, limit = output_limit)
#   saveRDS(
#     blast_results,
#     glue("{wkdir}/data/interim/blast_results/{ensembl_id}.rds")
#   )
#   return(blast_results)
# }


# first save fasta file in temp 
library(batchtools)
library(seqinr)
library(glue)
library(dplyr)

keyname_hits <- readRDS(glue("{wkdir}/data/interim/tmp/2023-01-27_gget_proteins-of-interest.rds"))
keyname_hits %>% glimpse()

wkdir <- "/central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation"
tmpdir <- "/central/scratch/jbok/tmp"

eid <- keyname_hits$ensembl_id[1]
seq <- keyname_hits$sequences_aa[1]
keyname_hits %>% glimpse()

alignment_df_list <- bind_rows(
  keyname_hits %>% 
    dplyr::select(ensembl_id, sequences_aa) %>% 
    mutate(method = "needle"),
  keyname_hits %>% 
    dplyr::select(ensembl_id, sequences_aa) %>% 
    mutate(method = "water")
)
alignment_df_list %>% glimpse()


#' write an R to write a fasta file in a temp directory and then run needleman wunsch and smith waterman alignment 
#' algorithm on the fasta file using the emboss package and UHGP 95v1 database

#' @param sequences_aa a character vector of amino acid sequences
#' @param ensembl_id a character vector of ensembl ids
#' @param db a character vector of the database to use for the blast search 
#' @param wkdir a character vector of the working directory
#' @param output_limit a numeric vector of the number of hits to return
#' @param method a character vector of the alignment algorithm to use
#' @return a data frame of the blast results
#' @export 

align_fasta_sequences <- function(ensembl_id, sequences_aa, method) {
  require(seqinr)
  require(glue)
  require(dplyr)
  require(stringr)
  tmpdir <- "/central/scratch/jbok/tmp"
  wkdir <- "/central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation"
  source(glue("{wkdir}/notebooks/R-scripts/helpers_general.R"))
  
  fasta_tmp <- glue("{tmpdir}/{ensembl_id}.fasta")
  write.fasta(
    sequences = sequences_aa, 
    names = ensembl_id, 
    file.out = fasta_tmp,
    open = "w", 
    nbchar = 60, 
    as.string = FALSE
    ) 

  refdb_path <- "/central/groups/MazmanianLab/joeB/Downloads/human-gut/v1.0/uhgp_catalogue/uhgp-95/uhgp-95.faa"
  output_dir <- glue("{wkdir}/data/interim/emboss_alignments")
  output_file <- glue("{output_dir}/{ensembl_id}_UHGP95v1_alignment.{method}")

  if (!file.exists(output_file)){
    shell_do(
      glue(
        "{method} -asequence {fasta_tmp}",
        " -bsequence {refdb_path}",
        " -gapopen 5.0",
        " -gapextend 0.5",
        " -datafile EBLOSUM62",
        " -outfile '{output_file}'"
      )
    )
  }
  # delete temp fasta file
  unlink(fasta_tmp)
}


# configure registry ----
breg <- makeRegistry(
  file.dir = glue(
    "{wkdir}/.cluster_runs/",
    "{Sys.Date()}_EMBOSS-Alignment_{sample(10000:99999, 1)}/"
  ),
  seed = 42
)
breg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
  template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
  scheduler.latency = 1,
  fs.latency = 65
)
# Submit Jobs ----
jobs <- batchMap(
  fun = align_fasta_sequences,
  args = alignment_df_list,
  reg = breg
)
setJobNames(jobs, paste("EMBOSS", alignment_df_list$ensembl_id, alignment_df_list$method, sep = "_"), reg = breg)
submitJobs(jobs,
  resources = list(walltime = 10800,
    memory = 1024,
    ncpus = 8,
    max.concurrent.jobs = 9999)
)
waitForJobs(sleep = 3)






#_______________________________________________________________________________\
# Magnify dataset 
























