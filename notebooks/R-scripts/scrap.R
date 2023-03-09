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

reticulate::use_condaenv(condaenv = "pdmbsR", required = TRUE)
gget <- import("gget")


gene_df <- read.delim(
  glue("{wkdir}/data/input/Immport_gene_lists/ImmportGeneList.txt"),
      stringsAsFactors = FALSE, header = TRUE
      )
go_gene_df <- read.delim(
  glue("{wkdir}/data/input/Immport_gene_lists/ImmportGeneListGOAnnotation.txt"),
      stringsAsFactors = FALSE, header = TRUE
      )

glimpse(gene_df)
glimpse(go_gene_df)

View(gene_df)

gene_df %>% glimpse()


immune_goi <- gene_df %>%  
  select(-c(Chromosome, Category)) %>% 
  filter(
    !grepl("HLA", Symbol),
    !grepl("immunoglobulin", Name),
    !grepl("IGK|IGL", Symbol),
    !grepl("T cell receptor", Name)
    ) %>% 
  distinct()
immune_goi %>% pull(Symbol)
immune_goi %>% glimpse()


collect_gget_search <- 


# Download Ensembl IDs for all genes of interest using gget
pb <- progress_bar$new(total = nrow(immune_goi))
gget_results <- tibble()
for (gene in immune_goi$Symbol) {
  pb$tick()
  ggout <- try({ 
    gget$search(gene, "homo_sapiens", 'gene') %>% 
      filter(gene_name == gene, biotype == "protein_coding") %>% 
      mutate_all(as.character)
    })
  if (class(ggout) == "data.frame") {
    gget_results %<>% bind_rows(ggout)
  }
}
saveRDS(
  gget_results,
  glue("{wkdir}/data/interim/tmp/2023-03-07_gget_proteins-of-interest-noseq.rds")
)

tst <- readRDS(
  glue("{wkdir}/data/interim/tmp/2023-03-07_gget_proteins-of-interest-noseq.rds"
  ))

immune_goi %>% glimpse
gget_results %>% glimpse
View(gget_results)
View(immune_goi)

gget_results <- readRDS(
  glue("{wkdir}/data/interim/tmp/2023-03-07_gget_proteins-of-interest-noseq.rds")
)

# double check that all genes of interest are in the gget_results
gget_results$gene_name %>% unique() %>% length()
immune_goi$Symbol %>% unique() %>% length()
setdiff(immune_goi$Symbol %>% unique(),
gget_results$gene_name %>% unique()
)


# Download the sequences for all genes of interest using gget$seq
  gget_seq <- function(ensembl_id, amino_acid = TRUE) {
    seq_results <- gget$seq(ensembl_id, translate = amino_acid)
    return(seq_results[2])
  }

gget_eid_translated <- gget_results$ensembl_id %>% 
  purrr::set_names() %>% 
  purrr::map(., possibly(gget_seq, otherwise = NA))

saveRDS(
  gget_eid_translated,
  glue("{wkdir}/data/interim/tmp/2023-03-07_gget_eid-translated.rds")
)


# Append the sequences to the gget_results data.frame
# TODO


# Save the gget_results data.frame
saveRDS(keyname_hits, glue(
  "{wkdir}/data/interim/tmp/",
  "{Sys.Date()}_gget_proteins-of-interest.rds"
))



# processed_eids <- list()
# # Collecting Human Immune proteome data
# IDB_IRIS <-
#   read.delim(
#     glue("{wkdir}/data/input/InnateDB_metadata/",
#     "InnateDB_genes_Immunogenetic-Related-Information-Source.txt"
#     ),
#     stringsAsFactors = FALSE, header = TRUE
#   )

# query <- c("IL2")
# tst <- gget$search(query, "homo_sapiens", 'gene')

# tst %>% 
#   filter(gene_name == query, biotype == "protein_coding")

# View(tst)


# IDB_IRIS %>% colnames()
# IDB_IRIS$ensembl %>% length()

# slurm_gget_info <- function(eids) {
#   require(data.table)
#   require(dplyr)
#   require(magrittr)
#   require(purrr)

#   reticulate::use_condaenv(condaenv = "pdmbsR", required = TRUE)
#   gget <- reticulate::import("gget")
#   #' function to collect amino acid sequences for each protein
#   #' using ensembl_id
  # gget_seq <- function(ensembl_id, amino_acid = TRUE) {
  #   seq_results <- gget$seq(ensembl_id, translate = amino_acid)
  #   return(seq_results[2])
  # }
#   eid_results <- list() #data.table()
#   for (eid in eids) {
#     result <- try(
#       {
#         gget$info(eid)
#       },
#       silent = TRUE
#     )
#     if (class(result) == "try-error") {
#       message("ERROR with: ", eid)
#     } else {
#       result_aa <- result %>%
#         mutate(
#           sequences_aa =
#             purrr::map_chr(ensembl_id, gget_seq)
#         )
#       eid_results[[eid]] <- result_aa
#     }
#   }
#   return(eid_results)
# }

# # Initiate future.batchtools backend for parallel processing
# future::plan(
#   future.batchtools::batchtools_slurm,
#   template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
#   resources = list(
#     name = "gget-info-IRIS",
#     memory = "50",
#     ncpus = 2,
#     walltime = 36000
#   )
# )

# eids_list <- IDB_IRIS$ensembl
# # eids_list %<>% keep(. %nin% names(processed_eids))

# # # Chunk files (10 per job) and download
# # n_jobs <- ceiling(length(eids_list) / 5)
# # download_runs <- listenv()
# # for (job in 1:n_jobs) {
# #   eids_chunk <- chunk_func(eids_list, n_jobs)[[job]]
# #   download_runs[[job]] %<-% slurm_gget_info(eids_chunk)
# # }

# # idb_df <- download_runs %>% as.list()
# # processed_eids %<>% append(unlist(idb_df, recursive = FALSE, use.names = TRUE))
# # saveRDS(
# #   processed_eids,
# #   glue(
# #     "{wkdir}/data/interim/tmp/",
# #     "{gsub(' ', '_', Sys.time())}",
# #     "_gget_proteins-of-interest.rds"
# #   )
# # )


# eids_translated <- eids_list %>% 
#   purrr::set_names() %>% 
#   purrr::map( ~ gget$seq(., translate = TRUE))

# saveRDS(eids_translated, glue(
#   "{wkdir}/data/interim/tmp/",
#   "{Sys.Date()}_IDB-IRIS_ensembl-ID_translated.rds"
# ))


# # tstme <- gget$info("ENSG00000158270") %>%
# #   mutate(
# #     sequences_aa =
# #       purrr::map_chr(ensembl_id, gget_seq)
# #   )
# # tstme
# # View(tstme)


# idb_df %>%
#   mutate(TGFB_hit = if_else(grepl("TGFB", protein_names), 1, 0)) %>%
#   mutate(TNFA_hit = if_else(grepl("TNFA", protein_names), 1, 0)) %>%
#   mutate(IL_hit = if_else(grepl("interleukin", protein_names), 1, 0)) %>%
#   mutate(IF_hit = if_else(grepl("interferons", protein_names), 1, 0)) %>%
#   arrange(TGFB_hit, TNFA_hit, IL_hit, IF_hit)

# saveRDS(keyname_hits, glue(
#   "{wkdir}/data/interim/tmp/",
#   "{Sys.Date()}_gget_proteins-of-interest-preAA.rds"
# ))


# library(progress)
# pb <- progress_bar$new(
#   total = length(IDB_IRIS$ensembl),
#   format = "[:bar] :current/:total (:percent)"
# )
# idb_df <- tibble()
# for (eid in IDB_IRIS$ensembl) {
#   message("Processing: ", eid)
#   pb$tick()
#   result <- try(
#     {
#       gget$info(eid)
#     },
#     silent = TRUE
#   )
#   # Process any error messages
#   if (class(result) == "try-error") {
#     # Ignore warnings while processing errors
#     message("ERROR with: ", eid)
#   } else {
#     idb_df %<>% bind_rows()
#   }
# }


# keyname_hits <- keyname_search %>%
  # filter(biotype == "protein_coding") %>%
  # filter(gene_name %in%
  #   c("IL1RN", "IL4", "IL6", "IL10", "IL11", "IL13"))

# # collect amino acid sequences for each protein using the ensembl_id
# gget_seq <- function(ensembl_id,  amino_acid=TRUE){
#     seq_results <- gget$seq(ensembl_id, translate = amino_acid)
#     return(seq_results[2])
# }
# keyname_hits %<>% mutate(sequences_aa =  map_chr(ensembl_id, gget_seq))

# saveRDS(keyname_hits, glue(
#   "{wkdir}/data/interim/tmp/",
#   "{Sys.Date()}_gget_proteins-of-interest.rds"
# ))





#______________________________________________________________________________
# Alignment Algorithms ----


proteins_of_interest <- c("interleukin", "interferons", "TGFB", "TNFA")
keyname_search <- gget$search(proteins_of_interest, "homo_sapiens")
keyname_hits <- keyname_search %>% filter(
  biotype == "protein_coding",
  !grepl("receptor", ensembl_description)
)
# collect amino acid sequences for each protein using the ensembl_id
gget_seq <- function(ensembl_id, amino_acid = TRUE) {
  seq_results <- gget$seq(ensembl_id, translate = amino_acid)
  return(seq_results[2])
}
keyname_hits %<>% mutate(sequences_aa =  map_chr(ensembl_id, gget_seq))
keyname_hits %>% glimpse
saveRDS(keyname_hits, glue("{wkdir}/data/interim/tmp/",
"{Sys.Date()}_gget_proteins-of-interest.rds")
)

#______________________________________________________________________________

# Save individual fasta files for each protein -----
keyname_hits <- readRDS(
  glue("{wkdir}/data/interim/tmp/",
  "2023-02-24_gget_proteins-of-interest.rds")
  )

fastas_raw_dir <- glue("{wkdir}/data/interim/fastas/raw/monomers")
for (id in keyname_hits$ensembl_id) {
  id_aa <- keyname_hits %>%
    filter(ensembl_id == id) %>%
    pull(sequences_aa)
  fasta_path <- glue("{fastas_raw_dir}/{id}.fasta")
  seqinr::write.fasta(
    sequences = id_aa,
    names = id,
    file.out = fasta_path,
    open = "w",
    nbchar = 60,
    as.string = FALSE
  )
}

# ______________________________________________________________________________
# Using SignalP 6.0 identify and remove signal peptides from fastas
fastas_raw_paths <-
  list.files(fastas_raw_dir, pattern = "fasta", full.names = TRUE) %>% 
  keep(grepl("ENSG", .,))

names(fastas_raw_paths) <- fastas_raw_paths %>%
  purrr::map(~ fs::path_ext_remove(basename(.)))

fastas_processed_keyname_hits <- list()
for (id in names(fastas_raw_paths)) {
  fasta_path <- fastas_raw_paths[[id]]
  fastas_processed_keyname_hits[[id]] <- glue(
    "{wkdir}/data/interim/fastas/processed/monomers/{id}"
  )
  if (dir.exists(fastas_processed_keyname_hits[[id]])) {
    next
  }
  message("Processing ", id, " fasta file...")
  slurm_shell_do(
    glue(
      "conda run -n signalp6",
      " signalp6",
      " --organism eukarya",
      " --fastafile {fasta_path}",
      " --output_dir {fastas_processed_keyname_hits[[id]]}",
      " --write_procs 4",
      " --format all",
      " --mode fast"
    ),
    ncpus = 4
  )
}

# IF no signal peptide is found, the fasta file is copied to the processed dir
for (id in names(fastas_raw_paths)) {
  raw_fasta <- fastas_raw_paths[[id]]
  processed_fasta <- glue(
    "{fastas_processed_keyname_hits[[id]]}/",
    "processed_entries.fasta"
    )
  if (file.info(processed_fasta)$size == 0){
    message("Overwriting ", id, " processed fasta file...")
    shell_do(glue("cp {raw_fasta} {processed_fasta}"))
  }
}

saveRDS(fastas_processed_keyname_hits, glue(
  "{wkdir}/data/interim/tmp/",
  "{Sys.Date()}_fasta-paths_processed_keyname_hits.rds"
))





#______________________________________________________________________________

# load in data 
keyname_hits <- readRDS(
  glue("{wkdir}/data/interim/tmp/",
  "2023-02-24_gget_proteins-of-interest.rds")
  )
protein_catalog_path <- "/central/groups/MazmanianLab/joeB/Downloads/protein_catalogs"
chunked_catalog_paths <- list.files(
  glue("{protein_catalog_path}/clustered_catalogs/merged/chunked_fasta"),
  full.names = TRUE) %>% 
  keep(grepl("chunk", .))

# create a datatable mapping each ensembl_id to each chunked catalog
binned_eids <- data.frame(
  "ensembl_id" = rep(
    keyname_hits$ensembl_id,
    each = length(chunked_catalog_paths)
  ),
  "refdb_path" = rep(
    chunked_catalog_paths,
    times = length(keyname_hits$ensembl_id)
  )
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


#______________________________________________________________________________
# EXECUTE ALIGNMENT JOBS ----
library(batchtools)

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
  scheduler.latency = 0.3,
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
    walltime = 172800,
    memory = 50000,
    ncpus = 1,
    max.concurrent.jobs = 9999
  )
)

# waitForJobs(sleep = 3)


#______________________________________________________________________________
# Aggregate Alignment Files ----


# # INITATE FUTURE-SLURM BACKEND FOR PARALLEL PROCESSING
# future::plan(
#   future.batchtools::batchtools_slurm,
#   template = glue(
#     "{wkdir}/batchtools_templates/",
#     "batchtools.slurm.tmpl"
#     ),
#   resources = list(
#     name = glue("aggregate-rds_{rand_string()}"),
#     memory = "128G",
#     ncpus = 1,
#     walltime = 10800
#   )
# )


collect_summary_metrics <- function(input_db_path, outputdir) {
  require(DBI)
  require(RSQLite)
  require(glue)
  require(dplyr)
  require(future)
  require(batchtools)
  require(listenv)
  require(stringr)

  seq_con <- DBI::dbConnect(RSQLite::SQLite(), input_db_path)
  message(glue(
    "{Sys.time()}: ",
    "Connecting to {DBI::dbListTables(seq_con)} .."
  ))

  # subfunctions for collecting data of interest
  get_subsample <- function(n_rand = 100000) {
    DBI::dbConnect(RSQLite::SQLite(), input_db_path) %>%
      dplyr::tbl(., DBI::dbListTables(.)[1]) %>%
      dplyr::slice_sample(n = n_rand) %>%
      dplyr::collect()
  }
  get_top_features <- function(col, top_n = 10000) {
    DBI::dbConnect(RSQLite::SQLite(), input_db_path) %>%
      dplyr::tbl(., DBI::dbListTables(.)[1]) %>% 
      dplyr::slice_max({{ col }}, n = top_n) %>%
      dplyr::collect()
  }
  get_summary_stats <- function() {
    DBI::dbConnect(RSQLite::SQLite(), input_db_path) %>%
    dplyr::tbl(., DBI::dbListTables(.)[1]) %>% 
      dplyr::summarise(
        across(
          c(Score, percent_identity, percent_similarity),
          list(mean = mean, median = median, sd = sd)
        )
      ) %>%
      dplyr::collect()
  }

  results <- listenv()
  message(glue("{Sys.time()}: Collecting DB subsample"))
  results[["subsample"]] %<-% {
    get_subsample()
  }

  message(glue(
    "{Sys.time()}: ",
    "Collecting top 10,0000 percent identity alignments"
  ))
  results[["top_percent_identity"]] %<-% {
    get_top_features(col = "percent_identity")
  }

  message(glue(
    "{Sys.time()}: ",
    "Collecting top 10,0000 percent similarity alignments"
  ))
  results[["top_percent_similarity"]] %<-% {
    get_top_features(col = "percent_similarity")
  }

  message(glue(
    "{Sys.time()}: ",
    "Collecting summary stats"
  ))
  results[["summary_stats"]] %<-% {
    get_summary_stats()
  }

  message(glue("{Sys.time()}: Disconnecting from DB"))
  DBI::dbDisconnect(seq_con)
  message(glue("{Sys.time()}: Collecting all data ..."))
  results_list <- as.list(results)

  output_path <- glue(
      "{outputdir}/",
      "{fs::path_ext_remove(basename(input_db_path))}.rds"
      )
  message(glue("{Sys.time()}: Saving data to \n{output_path}\n"))
  saveRDS(results_list, output_path)
}

aggregate_alignment_files <- function(input_files,
                                      output_dir_db,
                                      output_path_results
                                      id,
                                      working_dir = wkdir) {
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  message(glue("{Sys.time()}: Compiling {length(input_files)} files"))
  compiled_data <- load_rds_files(input_files)
  message(glue("{Sys.time()}: Saving compiled data to {output_dir_db}/{id}"))

  # create SQLite DB
  create_sql_db(
    table = compiled_data,
    table_name = id,
    db_loc = glue("{output_dir_db}/{id}.db"),
    overwrite_lg = TRUE,
    append_lg = FALSE
  )
  # collect summary metrics from sql db
  collect_summary_metrics(
    input_db_path = glue("{output_dir_db}/{id}.db"),
    outputdir = output_path_results
  )
  # remove SQLite DB
  shell_do("rm {output_dir_db}/{id}.db")
}

future::plan(
  list(
    tweak(
      future.batchtools::batchtools_slurm,
      template = glue(
        "{wkdir}/batchtools_templates/",
        "batchtools.slurm.tmpl"
      ),
      resources = list(
        name = glue("aggregate-rds_{rand_string()}"),
        memory = "64G",
        ncpus = 4,
        walltime = 108000
      )
    ),
    multisession
  )
)



# load in data 
keyname_hits <- readRDS(
  glue("{wkdir}/data/interim/tmp/",
  "2023-02-24_gget_proteins-of-interest.rds")
  )

eids <- keyname_hits$ensembl_id
file_dir <- glue("/central/scratch/jbok/emboss_alignments")
file_paths <- list.files(file_dir, pattern = "rds$", full.names = TRUE)
sql_dir <- "/central/scratch/jbok/emboss_alignments-SQL-DB"
results_dir <- glue({"{wkdir}/data/processed/emboss_alignments/2023-03-07"})
shell_do(glue("mkdir -p {sql_dir}"))

for (eid in eids) {
  message("Processing: ", eid)
  for (method in c("needle", "water")) {
    name_out <-
      glue(
        "{eid}_2023-02-13_catalog_",
        "MMSeq2-95_rep_seq.{method}"
      )
    # check if sql db exists
    output_check <- file.size(glue("{sql_dir}/{name_out}.db"))
    if (!is.na(output_check)) {
      if (output_check > 0) {
        message("Skipping, file exists: ", name_out, "\n")
        next
      }
    }
    message("Aggregating: " ,method, " ", eid)
    files <- file_paths %>%
      keep(grepl(
        glue(
          "{eid}_2023-02-13_catalog_",
          "MMSeq2-95_rep_seq_.*{method}.rds"
        ), .
      ))
    agg_files_needle %<-% {
      aggregate_alignment_files(
        input_files = files,
        output_dir_db = sql_dir,
        output_path_results = results_dir,
        id = name_out
      )
    }
  }
}










# # make a test DB
# fname <- "TEST_2"
# create_sql_db(
#   table = readRDS(glue(
#     "{wkdir}/data/interim/emboss_alignments-TEST/",
#     "ENSG00000095752_tst2-SPLIT2.water.rds"
#   )),
#   table_name = fname,
#   db_loc = glue("{wkdir}/{fname}.db")
# )

# # TESTRUN ON TEST DB
# input_db <- glue("{wkdir}/{fname}.db")
# testrun %<-% {
#   collect_summary_metrics(
#     input_db_path = input_db,
#     outputdir = wkdir
#   )
# }

# asum$subsample %>% glimpse()

# library(ggpointdensity)
# library(ggside)
# # read TEST OUTPUT
# asum$subsample %>%
#   ggplot(aes(Length, percent_similarity + 1)) +
#   geom_pointdensity(adjust = 0.7) +
#   scale_color_viridis(option = "G") +
#   geom_ysidedensity(aes(x = stat(density))) +
#   geom_xsidedensity(aes(y = stat(density))) +
#   geom_vline(aes(xintercept = 25)) +
#   labs(
#     y = "Percent Similarity (Needleman-Wunsch)", 
#     x = "Peptide Length",
#     title = "") +
#   theme_bw() +
#   scale_x_log10() +
#   scale_ysidex_continuous(guide = guide_axis(angle = 90), minor_breaks = NULL) +
#   scale_xsidey_continuous(minor_breaks = NULL) +
#     theme(
#       plot.title = element_text(hjust = 0.5),
#       ggside.panel.scale.x = .17,
#       ggside.panel.scale.y = .17,
#       legend.position = "bottom"
#     )



# asum <- readRDS(glue("{wkdir}/TEST.rds"))
# names(asum)
# asum %>% glimpse
# asum$top_percent_identity %>% glimpse
# asum$subsample %>% glimpse







#______________________________________________________________________________
# Visualizing Alignment Data ----


# # Function for visualization of alignment data ----
# summarize_alignment <- function(id, 
#                                 needle_rds_path, 
#                                 water_rds_path,
#                                 data_out_dir,
#                                 fig_out_dir
#                                 ) {
#   require(tidyverse)
#   require(glue)
#   require(ggpointdensity)
#   require(ggside)
#   require(viridis)
#   results <- list()
#   # Summary functions ----
#   get_summary_stats <- function(df) {
#     df %>%
#       dplyr::summarize(
#         alignment_pairs = n(),
#         mean_score = mean(Score),
#         min_score = min(Score),
#         max_score = max(Score),
#         mean_percent_identity = mean(percent_identity),
#         min_percent_identity = min(percent_identity),
#         max_percent_identity = max(percent_identity),
#         mean_percent_similarity = mean(percent_similarity),
#         min_percent_similarity = min(percent_similarity),
#         max_percent_similarity = max(percent_similarity)
#       )
#   }
#   plot_table_top_hits <- function(df, metric, nhits = 10000) {
#     df %>%
#       dplyr::select(contains(c("protein", "Length", "Score", "percent"))) %>%
#       dplyr::top_n(nhits, percent_similarity) %>%
#       dplyr::arrange(-percent_similarity) %>%
#       dplyr::mutate(algorithm = metric)
#   }
#   plot_alignment_distribution <- function(df, title_id = id) {
#   p <- df %>%
#     ggplot(aes(percent_similarity_needle, percent_similarity_water)) +
#     geom_pointdensity() +
#     scale_color_viridis(option = "G") +
#     geom_ysidedensity(aes(x = stat(density))) +
#     geom_xsidedensity(aes(y = stat(density))) +
#     labs(
#       x = "Percent Similarity (Needleman-Wunsch)", 
#       y = "Percent Similarity (Smith-Waterman)",
#       title = title_id) +
#     theme_bw() +
#     scale_ysidex_continuous(guide = guide_axis(angle = 90), minor_breaks = NULL) +
#     scale_xsidey_continuous(minor_breaks = NULL) +
#       theme(
#         plot.title = element_text(hjust = 0.5),
#         ggside.panel.scale.x = .17,
#         ggside.panel.scale.y = .17,
#         legend.position = "bottom"
#       )
#     return(p)
#   }

#   message(Sys.time(), " Reading in data...")
#   df_aligned_needle <-
#     readRDS(glue("{needle_rds_path}")) %>% filter(Length >= 25)
#   df_aligned_water <-
#     readRDS(glue("{water_rds_path}")) %>% filter(Length >= 25)

#   message(Sys.time(), " Summarizing in data...")
#   # 1) generate summary statistics
#   results[["stat_summary_needle"]] <- df_aligned_needle %>% get_summary_stats()
#   results[["stat_summary_water"]] <- df_aligned_water %>% get_summary_stats()
#   # 2) extract top hits 
#   results[["top_hits_needle"]] <- plot_table_top_hits(df_aligned_needle, "Needleman-Wunch")
#   results[["top_hits_water"]] <- plot_table_top_hits(df_aligned_water, "Smith-Waterman")
  
#   shell_do(glue("mkdir -p {fig_out_dir}"))
#   shell_do(glue("mkdir -p {data_out_dir}"))
#   message(Sys.time(), " Saving data summary...")
#   saveRDS(results, 
#     glue("{data_out_dir}/{id}_alignment_summary.rds")
#     )
  
#   # 3) plot alignment distribution
#   message(Sys.time(), " Formating data for plotting...")
#   needle_hits <- df_aligned_needle %>%
#     select(contains(c("protein", "Length", "Score", "percent"))) %>%
#     dplyr::rename_at(vars(!contains("protein")), ~ paste0(., "_needle"))
#   water_hits <- df_aligned_water %>%
#     select(contains(c("protein", "Length", "Score", "percent"))) %>%
#     dplyr::rename_at(vars(!contains("protein")), ~ paste0(., "_water"))
  
#   message(Sys.time(), " Plotting alignment distribution data...")
#   ggsave(
#     filename = glue("{fig_out_dir}/{id}_alignment_distribution.png"),
#     plot = 
#       plot_alignment_distribution(
#         full_join(needle_hits, water_hits), 
#         title_id = id
#         ),
#     width = 6, height = 6, dpi = 600
#   )
# }


alignment_results_dir <-
  glue("{wkdir}/data/processed/emboss_alignments/2023-03-02")
alignment_results_paths <- list.files(alignment_results_dir, full.names = TRUE)
# read in keyname_hits metadata
keyname_hits <- readRDS(
  glue(
    "{wkdir}/data/interim/tmp/",
    "2023-02-24_gget_proteins-of-interest.rds"
  )
)

gene_annots <- keyname_hits %>%
  select(ensembl_id, gene_name)
  
batch_info_df <- 
  tibble(paths = alignment_results_paths) %>%
  mutate(ensembl_id = paths %>% 
    purrr::map_chr( ~ basename(.) %>% 
      strex::str_before_nth(., "_", 2) %>% 
      strex::str_after_nth(., "_", 1)
    )) %>% 
  mutate(name_col = case_when(
    grepl("needle", paths) ~ "needle_rds_path",
    grepl("water", paths) ~ "water_rds_path",
    TRUE ~ "ERROR"
  )) %>% 
  pivot_wider(names_from = name_col, values_from = "paths") %>%
  left_join(gene_annots, by = "ensembl_id") %>%
  mutate(id = paste0(ensembl_id, "_", gene_name)) %>%
  mutate(data_out_dir = 
    glue(
      "{wkdir}/data/processed/emboss_alignments/",
      "2023-03-02_alignment_summaries/{id}"
      ),
      fig_out_dir = glue("{wkdir}/figures/sequence-alignment/{id}")) %>%
  select(-c(ensembl_id, gene_name))

batch_info_df %>% glimpse()
# View(batch_info_df)

# Launch batch jobs ----

# configure registry ----
cluster_run <- glue("{Sys.Date()}_Alignment-Summary-{rand_string()}/")
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
  fun = summarize_alignment,
  args = batch_info_df,
  reg = breg
)
# setJobNames(jobs,
#   paste("TGFB-EMBOSS",
#   alignment_df_list$ensembl_id,
#   alignment_df_list$method,
#   sep = "_"),
#   reg = breg
# )
# getJobNames(jobs, reg = breg)
submitJobs(jobs,
  resources = list(
    walltime = 36000,
    memory = 32000,
    ncpus = 1,
    max.concurrent.jobs = 9999
  )
)






df <- 
readRDS(glue(
  "{wkdir}/data/processed/emboss_alignments/2023-03-03/",
  "2023-03-03_ENSG00000112115_2023-02-13_catalog_MMSeq2-95_rep_seq.needle.rds"
))


sql_dir <- glue("{wkdir}/data/processed/emboss_alignments/SQL-DB")
shell_do(glue("mkdir -p {sql_dir}"))

fname <- "2023-03-03_ENSG00000112115_2023-02-13_catalog_MMSeq2-95_rep_seq.needle"
input_dir <- glue("{wkdir}/data/processed/emboss_alignments/2023-03-03")

create_sql_db(
  table = readRDS(glue(
    "{input_dir}/",
    "2023-03-03_ENSG00000112115_2023-02-13_catalog_MMSeq2-95_rep_seq.needle.rds"
  )),
  table_name = fname,
  db_loc = glue("{sql_dir}/{fname}.db")
)





tstdf <- batch_info_df %>% 
  filter(id == "ENSG00000136634_IL10") %>% 
  glimpse()

tic()
  df_aligned_needle <-
    readRDS(tstdf$needle_rds_path) %>% filter(Length >= 25)
toc()

tic() 
df_aligned_water <-
    readRDS(tstdf$water_rds_path) %>% filter(Length >= 25)
toc()

# Sample Alignment Summary ----
  needle_hits <- df_aligned_needle %>%
    select(contains(c("protein", "Length", "Score", "percent"))) %>%
    dplyr::rename_at(vars(!contains("protein")), ~ paste0(., "_needle"))
  water_hits <- df_aligned_water %>%
    select(contains(c("protein", "Length", "Score", "percent"))) %>%
    dplyr::rename_at(vars(!contains("protein")), ~ paste0(., "_water"))


#_______________________________________________________________________________
# HMMER Analysis ----

fastas_processed_keyname_hits <- readRDS(glue(
  "{wkdir}/data/interim/tmp/",
  "2023-03-02_fasta-paths_processed_keyname_hits.rds"
))

# Create HMM MSA profiles for each protein using jackhmmer against UniProt
for (fasta_path in fastas_processed_keyname_hits) {
  id <- fs::path_ext_remove(basename(fasta_path))
  output_dir <- glue(
    "{wkdir}/data/processed/",
    "jackhmmer_results/{id}"
  )
  shell_do(glue("mkdir -p {output_dir}"))
  if (!file.exists(glue("{output_dir}/{id}_msa.fasta")) |
    file.size(glue("{output_dir}/{id}_msa.fasta")) == 0) {
    message("Running Jackhmmer on: ", id)
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
        " {fasta_path}/processed_entries.fasta",
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

# _____________________________________________________________________________
# HMMER Search ----
# Hmmsearch run on MSA profile against custom database
protein_catalog_path <- glue(
  "{homedir}/Downloads/protein_catalogs/clustered_catalogs",
  "/merged/2023-02-13_catalog_MMSeq2-95_rep_seq.fasta"
)

# TODO: temporary fix for missing files ENSG00000134352 ENSG00000113302
ids <- fastas_processed_keyname_hits %>%
  purrr::map(~ fs::path_ext_remove(basename(.))) %>% 
  keep(. %nin% c("ENSG00000134352", "ENSG00000113302"))


for (id in ids) {
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






#_______________________________________________________________________________
# Downloading protein catalogs ----
# Read in metadata, select accession IDs for MAGs which are contamination and completeion criteria
#' Consistent with the definition of Medium Quality MAGs in the paper https://www.nature.com/articles/nbt.3893/tables/1
#' We will filter MAGs with a completion <= 50% and contamination < 10%

# hadza_metagenomes_dir <- glue("/central/scratch/jbok/PRJEB49206_HADZA")
# hadza_meta <- 
#   c("prokaryotes", "viruses", "eukaryotes") %>%
#   purrr::set_names() %>%
#   purrr::map(
#   ~ readxl::read_excel(
#     glue("{wkdir}/data/input/protein_catalog_metadata/Hadza-MAG-metadata.xlsx"),
#     sheet = .
#   )
# )

# #Explore quality metrics for each MAG type
# hadza_meta$viruses$checkv_quality %>% table()
# hadza_meta$prokaryotes$MAG_quality %>% table()
# hadza_meta$eukaryotes$MAG_quality %>% table()

# # Remove low-quality viruse MAGs
# hadza_meta$viruses <- hadza_meta$viruses %>%
#   filter(checkv_quality != "Low-quality")
# hadza_meta$prokaryotes %>% glimpse()


# hadza_meta$prokaryotes$assembly_acc %>% unique()  %>% length()
# hadza_meta$prokaryotes$mag_sample_acc %>% unique()  %>% length()
# hadza_meta$prokaryotes$mag_bin_acc %>% unique()  %>% length()

# col <- "mag_bin_acc"
# hadza_meta$prokaryotes[[col]] %>% length()
# hadza_meta$prokaryotes[[col]] %>% unique() %>% length()




#__________________________________________________________________________

# gget <- import("gget")
# proteins_of_interest <- c("interleukin", "interferons", "TGFB", "TNFA")
# keyname_search <- gget$search(proteins_of_interest, "homo_sapiens")
# keyname_hits <- keyname_search %>% filter(biotype == "protein_coding") %>% 
#   filter(gene_name %in% c("IL1RN", "IL4", "IL6", "IL10", "IL11", "IL13"))

# # collect amino acid sequences for each protein using the ensembl_id
# gget_seq <- function(ensembl_id,  amino_acid=TRUE){
#     seq_results <- gget$seq(ensembl_id, translate = amino_acid)
#     return(seq_results[2])
# }
# keyname_hits %<>% mutate(sequences_aa =  map_chr(ensembl_id, gget_seq))
# keyname_hits

# saveRDS(keyname_hits, glue("{wkdir}/data/interim/tmp/{Sys.Date()}_gget_proteins-of-interest.rds"))








library(bio3d)
library(r3dmol)

pdb_data <- bio3d::read.pdb(
  glue(
  "{wkdir}/data/interim/alphafold-multimer/",
  "HP-TGM_TGFBR1_TGFBR2/HP-TGM_TGFBR1_TGFBR2/ranked_0.pdb"
  ),
  multi = TRUE
)

r3dmol(
    viewer_spec = m_viewer_spec(
      cartoonQuality = 2,
      lowerZoomLimit = 50,
      upperZoomLimit = 2000
    )
  ) %>%
    m_add_model(data = m_bio3d(pdb_data)) %>%
    m_add_surface(style = m_style_surface(opacity = 0.4)) %>%
    m_set_style(style = m_style_cartoon()) %>%
    m_set_style(sel = m_sel(chain = "B"),
                style = m_style_cartoon(color = "#FFD105")) %>%
    m_set_style(sel = m_sel(chain = "C"),
                style = m_style_cartoon(color = "#FC027E")) %>% 
    m_zoom_to() %>% 
    m_rotate(angle = 90, axis = "y") %>%
    m_spin()


library(reticulate)
reticulate::use_condaenv(condaenv = "pdmbsR", required = TRUE)


source_python(glue("{src_dir}/python-scripts/pickle_reader.py"))
pickle_data <- read_pickle_file(glue(
  "{wkdir}/data/interim/alphafold-multimer/",
  "HP-TGM_TGFBR2/HP-TGM_TGFBR2/features.pkl"
))

pickle_data <- read_pickle_file(glue(
  "{wkdir}/data/interim/alphafold-multimer/",
  "2023-02-2"
  "HP-TGM_TGFBR2/HP-TGM_TGFBR2/"
))


str(pickle_data)

"data/interim/alphafold-multimer/2023-02-28_TGFB_complexes/HP-TGM_TGFBR3/HP-TGM_TGFBR3/"


# #______________________________________________________________________________
# # GGET Blast ----
# library(batchtools)
# blastr <- function(sequences_aa,
#                    ensembl_id,
#                    db = "refseq_protein",
#                    wkdir = wkdir,
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

# # configure registry ----
# breg <- makeRegistry(
#   file.dir = glue(
#     "{wkdir}/.cluster_runs/",
#     "{Sys.Date()}_BLAST_{sample(10000:99999, 1)}/"
#   ),
#   seed = 42
# )
# breg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
#   template = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
#   scheduler.latency = 2,
#   fs.latency = 65
# )
# # Submit Jobs ----
# jobs <- batchMap(
#   fun = blastr,
#   args = dplyr::select(keyname_hits, sequences_aa, ensembl_id),
#   reg = breg
# )
# setJobNames(jobs, paste0("BLAST_", keyname_hits$ensembl_id), reg = breg)
# submitJobs(jobs,
#   resources = list(walltime = 10800,
#     memory = 1024,
#     ncpus = 8,
#     max.concurrent.jobs = 9999)
# )

# # getJobNames()$job.id
# # waitForJobs(sleep = 5)
# # tst <- loadResult(1)
# # tst %>% glimpse
