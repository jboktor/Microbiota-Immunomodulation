
library(tidyverse)
library(reticulate)
library(magrittr)
library(glue)
library(seqinr)
library(future)
library(future.batchtools)
library(fs)
library(tictoc)
library(listenv)
library(progress)

tmpdir <- "/central/scratch/jbok/tmp"
homedir <- "/central/groups/MazmanianLab/joeB"
wkdir <- glue(
  "{homedir}/",
  "Microbiota-Immunomodulation/Microbiota-Immunomodulation"
)
src_dir <- glue("{wkdir}/notebooks")
source(glue("{src_dir}/R-scripts/helpers_general.R"))
source(glue("{src_dir}/R-scripts/helpers_sequence-screens.R"))
protein_catalogs <- glue("{homedir}/", "Downloads/protein_catalogs")
shell_do(glue("mkdir -p {wkdir}/data/interim/fastas/raw/monomers"))
shell_do(glue("mkdir -p {wkdir}/data/interim/fastas/processed/monomers"))
shell_do(glue("mkdir -p {wkdir}/data/interim/fastas/processed/complexes"))

# ______________________________________________________________________________
#                            TGF-B & HP-TGM Analysis
# ______________________________________________________________________________
# Sequence-level similarity

gget_seq <- function(ensembl_id, amino_acid = TRUE) {
  seq_results <- gget$seq(ensembl_id, translate = amino_acid)
  return(seq_results[2])
}
tgf_keys <- c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3")
tgfb_search <- gget$search(tgf_keys, "homo_sapiens")
tgfb_search %<>% filter(biotype == "protein_coding", gene_name %in% tgf_keys)
# collect amino acid sequences for each protein using the ensembl_id
tgfb_search %<>% mutate(sequences_aa =  map_chr(ensembl_id, gget_seq))

HP_TMG_SEQ <- "DDSGCMPFSDEAATYKYVAKGPKNIEIPAQIDNSGMYPDYTHVKRFCKGLHGEDTTGWFVGICLASQWYYYEGVQECDDRRCSPLPTNDTVSFEYLKATVNPGIIFNITVHPDASGKYPELTYIKRICKNFPTDSNVQGHIIGMCYNAEWQFSSTPTCPASGCPPLPDDGIVFYEYYGYAGDRHTVGPVVTKDSSGNYPSPTHARRRCRALSQEADPGEFVAICYKSGTTGESHWEYYKNIGKCPDPRCKPLEANESVHYEYFTMTNETDKKKGPPAKVGKSGKYPEHTCVKKVCSKWPYTCSTGGPIFGECIGATWNFTALMECINARGCSSDDLFDKLGFEKVIVRKGEGSDSYKDDFARFYATGSKVIAECGGKTVRLECSNGEWHEPGTKTVHRCTKDGIRTL"
TGFB1 <- "ALDTNYCFSSTEKNCCVRQLYIDFRKDLGWKWIHEPKGYHANFCLGPCPYIWSLDTQYSKVLALYNQHNPGASAAPCCVPQALEPLPIVYYVGRKPKVEQLSNMIVRSCKCS"
TGFB2 <- "ALDAAYCFRNVQDNCCLRPLYIDFKRDLGWKWIHEPKGYNANFCAGACPYLWSSDTQHSRVLSLYNTINPEASASPCCVSQDLEPLTILYYIGKTPKIEQLSNMIVKSCKCS"
TGFB3 <- "ALDTNYCFRNLEENCCVRPLYIDFRQDLGWKWVHEPKGYYANFCSGPCPYLRSADTTHSTVLGLYNTLNPEASASPCCVPQDLEPLTILYYVGRTPKVEQLSNMVVKSCKCS"

hp_mimic_data <-
  tribble(
    ~gene_name, ~sequences_aa,
    "HP-TGM",   HP_TMG_SEQ
  )
tgfb_search %<>% bind_rows(hp_mimic_data)

saveRDS(tgfb_search, glue(
  "{wkdir}/data/interim/tmp/",
  "{Sys.Date()}_TGFB_gget.rds"
))



#______________________________________________________________________________

# Save individual fasta files for each protein -----
tgfb_search <- readRDS(glue(
  "{wkdir}/data/interim/tmp/",
  "2023-02-27_TGFB_gget.rds"
))

# Manually editing the fasta sequences to remove the signal peptide and other non-active domains
tgfb_search %<>%
  mutate(sequences_aa = case_when(
    gene_name == "TGFB1" ~ TGFB1,
    gene_name == "TGFB2" ~ TGFB2,
    gene_name == "TGFB3" ~ TGFB3,
    TRUE ~ sequences_aa
  ))

fastas_raw_dir <- glue("{wkdir}/data/interim/fastas/raw/monomers")
for (tf in tgfb_search$gene_name) {
  tf_aa <- tgfb_search %>%
    filter(gene_name == tf) %>%
    pull(sequences_aa)
  fasta_path <- glue("{fastas_raw_dir}/{tf}.fasta")
  write.fasta(
    sequences = tf_aa,
    names = tf,
    file.out = fasta_path,
    open = "w",
    nbchar = 60,
    as.string = FALSE
  )
}

# ______________________________________________________________________________
# Using SignalP 6.0 identify and remove signal peptides from fastas
fastas_raw_paths <- paste0(
  fastas_raw_dir, "/",
  list.files(fastas_raw_dir, pattern = "fasta")
)
names(fastas_raw_paths) <- fastas_raw_paths %>%
  purrr::map(~ fs::path_ext_remove(basename(.)))

fastas_processed_tgfb <- list()
for (id in names(fastas_raw_paths)) {
  fasta_path <- fastas_raw_paths[[id]]
  fastas_processed_tgfb[[id]] <- glue(
    "{wkdir}/data/interim/fastas/processed/monomers/{id}"
  )
  if (dir.exists(fastas_processed_tgfb[[id]])) {
    next
  }
  message("Processing ", id, " fasta file...")
  slurm_shell_do(
    glue(
      "conda run -n signalp6",
      " signalp6",
      " --organism eukarya",
      " --fastafile {fasta_path}",
      " --output_dir {fastas_processed_tgfb[[id]]}",
      " --write_procs 8",
      " --format all",
      " --mode fast"
    ),
    ncpus = 8
  )
}

# IF no signal peptide is found, the fasta file is copied to the processed dir
for (id in names(fastas_raw_paths)) {
  raw_fasta <- fastas_raw_paths[[id]]
  processed_fasta <- glue(
    "{fastas_processed_tgfb[[id]]}/",
    "processed_entries.fasta"
    )
  if (file.info(processed_fasta)$size == 0){
    message("Overwriting ", id, " processed fasta file...")
    shell_do(glue("cp {raw_fasta} {processed_fasta}"))
  }
}



# ' Signal peptide prediction model based on a Bert protein language 
# ' model encoder and a conditional random field (CRF) decoder.
# ' Methods 
# ' organism as eukarya to trigger post-processing 
# ' of the SP predictions to prevent spurious results 
# ' (only predicts type Sec/SPI).
# ' using fast mode 

#______________________________________________________________________________
# Saving protein complexes (post-processing) ----
tgfb_search <- readRDS(glue(
  "{wkdir}/data/interim/tmp/",
  "2023-02-27_TGFB_gget.rds"
))
# Split into TFs and Receptors and combine all possible pairs
tgf_receptors <- tgfb_search %>% 
  filter(grepl("receptor", ensembl_description)) %>% 
  pull(gene_name)
tgf_tf <- tgfb_search %>% 
  filter(!grepl("receptor", ensembl_description)) %>% 
  pull(gene_name)

format_seq <- function(fasta) {
  fasta %>% 
    seqinr::getSequence() %>% 
    unlist() %>%
    paste(collapse = "")
}

# Saving paired fasta files
fastas_merged_path <- glue("{wkdir}/data/interim/fastas/processed/complexes")
for (tf in tgf_tf) {
  tf_aa <- seqinr::read.fasta(
      glue("{fastas_processed_tgfb[[tf]]}/processed_entries.fasta"),
      seqtype = "AA"
      )
  for (receptor in tgf_receptors) {
    rec_aa <- seqinr::read.fasta(
      glue("{fastas_processed_tgfb[[receptor]]}/processed_entries.fasta"),
      seqtype = "AA"
      )
    merged_name <- paste(
      seqinr::getName(tf_aa), 
      seqinr::getName(rec_aa), 
      sep = "_"
      )
    merged_fasta_path <- glue("{fastas_merged_path}/{merged_name}.fasta")
    if (file.exists(merged_fasta_path)) {
      next
    }
    message("Saving: ", merged_name)
    write.fasta(
      sequences = list(format_seq(tf_aa), format_seq(rec_aa)),
      names =list(seqinr::getName(tf_aa), seqinr::getName(rec_aa)),
      file.out = merged_fasta_path,
      open = "w",
      nbchar = 60,
      as.string = FALSE
    )
  }
}

# Saving trioed fasta files
rec1 <- seqinr::read.fasta(
  glue("{fastas_processed_tgfb[['TGFBR1']]}/processed_entries.fasta"),
  seqtype = "AA"
  )
rec2 <- seqinr::read.fasta(
  glue("{fastas_processed_tgfb[['TGFBR2']]}/processed_entries.fasta"),
  seqtype = "AA"
  )
for (tf in tgf_tf) {
  tf_aa <- seqinr::read.fasta(
      glue("{fastas_processed_tgfb[[tf]]}/processed_entries.fasta"),
      seqtype = "AA"
      )
  merged_name <- paste(seqinr::getName(tf_aa), "TGFBR1_TGFBR2", sep = "_")
  merged_fasta_path <- glue("{fastas_merged_path}/{merged_name}.fasta")
  if (file.exists(merged_fasta_path)) {
    next
  }
  message("Saving: ", merged_name)
  write.fasta(
    sequences = list(
      format_seq(tf_aa), 
      format_seq(rec1),
      format_seq(rec2)
      ),
    names =list(
      seqinr::getName(tf_aa),
      seqinr::getName(rec1), 
      seqinr::getName(rec2)
      ),
    file.out = merged_fasta_path,
    open = "w",
    nbchar = 60,
    as.string = FALSE
  )
}



#______________________________________________________________________________

# Sequence alignment

#' For each pair of samples calculate the needleman-wunch
#'  and smith waterman alignment scores

# Creating a slurm batchtools dataframe of job metadata
sequence_analysis_dir <- glue(
  "{wkdir}/data/interim/",
  "tgfb_analyses/sequence-alignment"
)
tgfb_ligand_dir <- glue("{wkdir}/data/interim/alphafold-multimer/TGFB_ligands")
shell_do(glue("mkdir -p {sequence_analysis_dir}"))
shell_do(glue("mkdir -p {tgfb_ligand_dir}"))

seqinr::write.fasta(
  sequences = as.list(tgf_tf$sequences_aa),
  names = as.list(tgf_tf$gene_name),
  file.out = glue("{tgfb_ligand_dir}/tgfb_ligands.fasta"),
  open = "w",
  nbchar = 60,
  as.string = FALSE
)

alignment_df_list <- bind_rows(
  tgf_tf %>%
    mutate(ensembl_id = gene_name) %>%
    dplyr::select(ensembl_id, sequences_aa) %>%
    mutate(method = "needle"),
  tgf_tf %>%
    mutate(ensembl_id = gene_name) %>%
    dplyr::select(ensembl_id, sequences_aa) %>%
    mutate(method = "water")
) %>%
  mutate(
    refdb_path = glue("{tgfb_ligand_dir}/tgfb_ligands.fasta"),
    output_dir = sequence_analysis_dir
  )

alignment_df_list %>% glimpse()


#______________________________________________________________________________
# EXECUTE ALIGNMENT JOBS ----
library(batchtools)

# configure registry ----
cluster_run <- glue("{Sys.Date()}_TGFB-EMBOSS-Alignment_ID-{rand_string()}/")
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
  paste("TGFB-EMBOSS",
  alignment_df_list$ensembl_id,
  alignment_df_list$method,
  sep = "_"),
  reg = breg
)
getJobNames(jobs, reg = breg)
submitJobs(jobs,
  resources = list(
    walltime = 3600,
    memory = 1024,
    ncpus = 1,
    max.concurrent.jobs = 9999
  )
)



#______________________________________________________________________________
# Aggregate Alignment Files ----

file_paths <- list.files(sequence_analysis_dir, pattern = "rds$", full.names = TRUE)
needle_files <- file_paths %>%
  keep(grepl("HP-TGM_tgfb_ligands.needle.rds", .))
water_files <- file_paths %>%
  keep(grepl("HP-TGM_tgfb_ligands.water.rds", .))


tgfb_alignments <-  
  bind_rows(
    file_paths %>% 
      keep(grepl("needle", .)) %>% 
      load_rds_files() %>% 
      mutate(method = "Needleman-Wunsch (global) "),
    file_paths %>% 
      keep(grepl("water", .)) %>% 
      load_rds_files() %>% 
      mutate(method = "Smith-Waterman (local)")
  )


# View(tgfb_alignments)
tgfb_alignments %>% glimpse()

library(viridis) 

alignment_heatmap <- tgfb_alignments %>% 
  ggplot(aes(x = protein_1, y = protein_2, fill = Score)) +
  geom_tile() +
  geom_text(aes(label = paste0(percent_identity, "%\n", percent_similarity, "%")), 
    color = "white", size = 3) +
  theme_bw() +
  facet_wrap( ~ method) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  scale_fill_viridis_c(option = "G") +
  labs(x = "", y = "", fill = "Alignment\nScore") +
  theme(axis.text.x = element_text(angle=45, hjust=1) )


ggsave(alignment_heatmap, filename = glue("{wkdir}/figures/tgfb-analysis/sequence-alignment-heatmap.png"),
dpi = 300, width = 7, height = 3
)



# HMMER Approach _____________________________________________________

# # Save individual fasta files for each protein
# tgf_tf <- readRDS(glue(
#   "{wkdir}/data/interim/tmp/",
#   "2023-02-25_TGFBR-ligands_gget.rds"
# ))

# tgfb_fasta_dir <- glue("{wkdir}/data/interim/alphafold-multimer/TGFB_fasta")
# for (tf in tgf_tf$gene_name) {
#   tf_aa <- tgf_tf %>%
#     filter(gene_name == tf) %>%
#     pull(sequences_aa)
#   fasta_path <- glue("{tgfb_fasta_dir}/{tf}.fasta")
#   write.fasta(
#     sequences = tf_aa,
#     names = tf,
#     file.out = fasta_path,
#     open = "w",
#     nbchar = 60,
#     as.string = FALSE
#   )
# }

# tgfb_fasta_paths <- paste0(
#   tgfb_fasta_dir, "/", 
#   list.files(tgfb_fasta_dir, pattern = "fasta") %>% keep(!grepl("_", .))
# )



# Create HMM MSA profiles for each protein using jackhmmer against UniProt
for (fasta_path in tgfb_fasta_paths) {
  id <- fs::path_ext_remove(basename(fasta_path))
  output_dir <- glue(
    "{wkdir}/data/processed/",
    "jackhmmer_results/{id}"
  )
  shell_do(glue("mkdir -p {output_dir}"))
  if (!file.exists(glue("{output_dir}/{id}_msa.fasta"))) {
    slurm_shell_do(
      cmd = glue(
        "jackhmmer ",
        " --cpu 8",
        " -N 5",
        " -E 0.05",
        " -A {output_dir}/{id}_msa.fasta",
        " -o {output_dir}/{id}_jackhmmer.out",
        " --tblout {output_dir}/{id}_seqhits.txt",
        " --domtblout {output_dir}/{id}_domainhits.txt",
        " --noali",
        " {fasta_path}",
        " /central/software/alphafold/data/uniprot/uniprot.fasta"
      ),
      jobname = glue("jkhmr_{fasta_path}_{rand_string()}"),
      working_dir = wkdir,
      ncpus = 8,
      memory = "4G",
      walltime = 36000
    )
  }
}


# Hmmsearch run on MSA profile against custom database
protein_catalog_path <- glue(
  "{homedir}/Downloads/protein_catalogs/clustered_catalogs",
  "/merged/2023-02-13_catalog_MMSeq2-95_rep_seq.fasta"
)
ids <- tgfb_fasta_paths %>%
  purrr::map(~ fs::path_ext_remove(basename(.)))

for (id in ids){
  input_dir <- glue(
    "{wkdir}/data/processed/",
    "jackhmmer_results/{id}"
  )
  output_dir <- glue(
    "{wkdir}/data/processed/",
    "hmmersearch_results/{id}"
  )
  shell_do(glue("mkdir -p {output_dir}"))
  # make hmm profile from MSA
  fasta_path <- glue("{input_dir}/{id}_msa.fasta")
  shell_do(glue("hmmbuild {output_dir}/{id}.hmm {fasta_path}"))
  # hmmsearch on protein catalog using hmm profile
  slurm_shell_do(
    cmd = glue(
      "hmmsearch ",
      " --cpu 8",
      " -A {output_dir}/{id}_msa.hmm",
      " -o {output_dir}/{id}_hmmsearch.out",
      " --tblout {output_dir}/{id}_tblout.txt",
      " --domtblout {output_dir}/{id}_domtblout.txt",
      " --noali",
      " {output_dir}/{id}.hmm",
      " {protein_catalog_path}"
    ),
    jobname = glue("hmrsearch_{fasta_path}_{rand_string()}"),
    working_dir = wkdir,
    ncpus = 8,
    memory = "4G",
    walltime = 36000
  )
}




# _____________________________________________________________________________

# Read in hmmsearch results
#' Source code for rhmmer package manually copied from: 
#' https://github.com/arendsee/rhmmer/blob/master/R/parse.R
library(readr)

source(glue("{src_dir}/R-scripts/rhmmer-package.R"))


domain_hits <- read_domtblout(
    glue("{wkdir}/data/processed/hmmersearch_results/test/test_domainhits.txt")
)
sequence_hits <- read_tblout(
  glue("{wkdir}/data/processed/hmmersearch_results/test/test_seqhits.txt")
)
domain_hits %<>%
  dplyr::rename_at(vars(-"domain_name"), ~ paste0(., "_domain_hits"))
sequence_hits %<>%
  dplyr::rename_at(vars(-"domain_name"), ~ paste0(., "_sequence_hits"))

hmmer_hits <- domain_hits %>%
  full_join(sequence_hits)


catalogs <- c(
  "ELGG",
  "Hadza",
  "KIJ",
  "MGV",
  "RUGS",
  "RUMMETA",
  "UHGP",
  "VEuPathDB",
  "Wormbase"
)

catalogs_cols <- pal_npg()(length(catalogs))
names(catalogs_cols) <- catalogs

hmmer_results_dir <- glue("{wkdir}/data/interim/tgfb_analyses/hmmersearch_results/")
hmmer_results_paths <- list.dirs(
  glue("{wkdir}/data/interim/tgfb_analyses/hmmersearch_results/"), 
  recursive = FALSE)


hmmer_hit_plots <- list()
for (id in basename(hmmer_results_paths)) {
  sequence_hits <- read_tblout(
  glue("{hmmer_results_dir}/{id}/{id}_tblout.txt")
  )
  hmmer_hit_plots[[id]] <- sequence_hits %>% 
    mutate(catalog = map_chr(
        str_before_nth(domain_name, "_", 2), 
        ~ str_remove(., "CATID_"))
        ) %>%
      ggplot(aes(
        x = -log10(sequence_evalue), 
        y = -log10(best_domain_evalue),
        group = description
        )) +
      scale_color_manual(values = catalogs_cols) +
      geom_point(aes(color = catalog)) +
      theme_bw()
}

ggplotly(hmmer_hit_plots$`HP-TGM` )

ggplotly(hmmer_hit_plots$TGFB1)
ggplotly(hmmer_hit_plots$TGFB2)
ggplotly(hmmer_hit_plots$TGFB3)


# Color 


# library(bio3d)
# tst_jack <- bio3d::hmmer(
#   seq = bio3d::as.fasta(matrix(tgf_tf$sequences_aa[4])),
#   type = "jackhmmer",
#   db = "uniprotkb",
#   verbose = TRUE,
#   timeout = 120
# )

# tst_jack2 <- bio3d::hmmer(
#   seq = bio3d::read.fasta("seqs.fasta"),
#   type = "jackhmmer",
#   db = "uniprotkb",
#   verbose = TRUE,
#   timeout = 120
# )

# jackhmmer_wrapper <- function(sequence) {
#   jkhmodut <- bio3d::hmmer(
#     seq = sequence,
#     type = "jackhmmer",
#     db = "uniprotkb",
#     verbose = TRUE,
#     timeout = 120
#   )
#   aa <- bio3d::get.seq(
#     ids = jkhmout$hit.tbl$acc2,
#     outfile = glue(
#       "{wkdir}/data/interim/",
#       "tmp/seqs_{rand_string()}.fasta"
#     ),
#     db = "uniprot",
#     verbose = FALSE
#   )
#   return(aa)
# }





#______________________________________________________________________________

# library(UniprotR) 
# #Read Accessions from csv file , Note : Accessions must be in the first column. 
# Accessions <-GetAccessionList("https://s3.amazonaws.com/csvpastebin/uploads/9571fa356c67a0c7c95e8431799a051a/Accessions.csv") 
# #Get Taxanomy Information 
# TaxaObj <- GetNamesTaxa(Accessions) 
# #Visualize Chromosomes localization
# PlotChromosomeInfo(TaxaObj)
# #Visualize protein's gene name as Network 
# PlotGenesNetwork(TaxaObj)
