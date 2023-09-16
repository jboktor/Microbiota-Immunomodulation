
#'Script to read a text file, select only lines that start with a hash, 
#' and parse the data into a data frame
#' @param file_path a character vector of the path to the text file
#' @return a data frame of the parsed data
parse_emboss_alignment_file_vect <- function(file_path) {
  suppressPackageStartupMessages(require(magrittr))
  suppressPackageStartupMessages(require(tibble))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(purrr))

  text_file <- readLines(file_path, n = -1)
  hash_lines <- grep("^#", text_file, value = TRUE) %>%
    grep(":", ., value = TRUE)

  df_results <- data.table(matrix(
    ncol = 10,
    nrow = length(grep("Similarity", hash_lines, value = TRUE)
    )))
  colnames(df_results) <-
    c("1", "2", "Matrix", "Gap_penalty", "Extend_penalty",
      "Length", "Identity", "Similarity", "Gaps", "Score"
    )

  line_split <- purrr::map(hash_lines, ~ strsplit(., ":")[[1]]) %>%
    purrr::map(~ gsub("# ", "", .)) %>%
    purrr::map(~ gsub(" ", "", .))

  # select list elements of interest
  line_split %<>%
    keep(unlist(purrr::map(
      line_split, ~ .[c(1)] %in% colnames(df_results)
      )))

  get_matching_nlist_elements <- function(l, i, match) {
    loi <- l %>% keep(unlist(purrr::map(
    l, ~ .[c(1)] == match
    )))
    loi %>% map_chr(c(i))
  }
  for (col in colnames(df_results)) {
    df_results[, col] <- get_matching_nlist_elements(line_split, 2, col)
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
    mutate_at(
      vars(Gap_penalty, Extend_penalty, Length, Score, contains("percent_")),
      as.numeric
      ) %>%
      rename(protein_1 = `1`, protein_2 = `2`)

  return(df_results)
  }


#' Funtion to write a fasta file in a temp directory and then run
#' needleman wunsch and smith waterman alignment
#' algorithms on the fasta file using the emboss package
#' and a binned reference fasta catalog
#' @param sequences_aa a character vector of amino acid sequences
#' @param ensembl_id a character vector of ensembl ids
#' @param refdb_path a character vector of the database
#'  to use for the blast search
#' @param wkdir a character vector of the working directory
#' @param output_limit a numeric vector of the number of hits to return
#' @param method a character vector of the alignment algorithm to use
#' @return a data frame of the blast results
align_fasta_sequences <- function(ensembl_id,
                                  fasta_path,
                                  method,
                                  refdb_path,
                                  output_dir =
                                    "/central/scratch/jbok/emboss_alignments") {
  suppressPackageStartupMessages(require(seqinr))
  suppressPackageStartupMessages(require(glue))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(stringr))
  tmpdir <- "/central/scratch/jbok/tmp"
  group_dir <- "/central/groups/MazmanianLab/joeB"
  wkdir <-
    glue("{group_dir}/Microbiota-Immunomodulation/Microbiota-Immunomodulation")
  source(glue("{wkdir}/notebooks/R-scripts/helpers_general.R"))
  source(glue("{wkdir}/notebooks/R-scripts/helpers_sequence-screens.R"))

  refdb_name <- fs::path_ext_remove(basename(refdb_path))
  alignment_file <- glue("{output_dir}/{ensembl_id}_{refdb_name}.{method}")
  rds_out <- glue("{output_dir}/{ensembl_id}_{refdb_name}.{method}.rds")

  if (file.exists(rds_out)) {
    message(
      Sys.time(),
      ": Skipping alignment \n", rds_out, "\n output file exists"
    )
    if (file.exists(alignment_file)) {
      shell_do(glue("rm {alignment_file}"))
    }
    next
  }

  if (!file.exists(alignment_file)) {
    # fasta_tmp <- glue("{tmpdir}/{ensembl_id}.fasta")
    # message(Sys.time(), ": Writing fasta file to ", fasta_tmp)
    # write.fasta(
    #   sequences = sequences_aa,
    #   names = ensembl_id,
    #   file.out = fasta_tmp,
    #   open = "w",
    #   nbchar = 60,
    #   as.string = FALSE
    # )
    message(Sys.time(), ": Executing ", method, " on ", fasta_path)
    cmd <- glue(
      "{method} -asequence {fasta_path}",
      " -bsequence {refdb_path}",
      " -gapopen 5.0",
      " -gapextend 0.5",
      " -datafile EBLOSUM62",
      " -outfile '{alignment_file}'"
    )
    message(Sys.time(), ": Executing command:  ", cmd)
    shell_do(cmd)
    message(Sys.time(), ": Removing tmp fasta file: ", fasta_path)
    # shell_do(glue("rm {fasta_path}"))
  } else {
    message(
      Sys.time(),
      ": Skipping ", method, " : ", alignment_file, " already file exists"
    )
  }
  if (!file.exists(rds_out)) {
    message(Sys.time(), ": Parsing ", alignment_file)
    parsed_df <- parse_emboss_alignment_file_vect(alignment_file)
    saveRDS(parsed_df, rds_out)
    # delete aligment file
    shell_do(glue("rm {alignment_file}"))
  }
}

# Download the sequences for all genes of interest using gget$seq
gget_seq <- function(ensembl_id, amino_acid = TRUE) {
  seq_results <- gget$seq(ensembl_id, translate = amino_acid)
  return(seq_results[2])
}

#' Function to translate ensembl IDs into AA sequences using gget
slurm_gget_seq <- function(eids) {
  require(purrr)
  require(reticulate)
  reticulate::use_condaenv(condaenv = "pdmbsR", required = TRUE)
  gget <- reticulate::import("gget")
  gget_seq <- function(ensembl_id, amino_acid = TRUE) {
    seq_results <- gget$seq(ensembl_id, translate = amino_acid)
    return(seq_results[2])
  }
  eids %>%
    purrr::set_names() %>%
    purrr::map(., purrr::possibly(gget_seq, otherwise = NA))
}



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
  get_subsample <- function() {
    # read in random protein catalog selection
    rand_p <- readRDS(
      glue(
        "/central/groups/MazmanianLab/joeB/",
        "Microbiota-Immunomodulation/Microbiota-Immunomodulation/",
        "data/interim/emboss_alignments/random-protein-2_N99968.rds"
      )
    )
    query <- sprintf(
      "SELECT * FROM results WHERE protein_2 IN ('%s')",
      paste0(rand_p, collapse = "','")
    )
    DBI::dbGetQuery(
      DBI::dbConnect(RSQLite::SQLite(), input_db_path),
      query
    )
  }
  get_top_features <- function(col, top_n = 1000) {
    DBI::dbConnect(RSQLite::SQLite(), input_db_path) %>%
      dplyr::tbl(., DBI::dbListTables(.)[1]) %>%
      dplyr::slice_max(order_by = {{ col }}, n = top_n, with_ties = FALSE) %>%
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
  results[["subsample"]] %<-% {
    get_subsample()
  }
  results[["top_percent_identity"]] %<-% {
    get_top_features(col = "percent_identity")
  }
  results[["top_percent_similarity"]] %<-% {
    get_top_features(col = "percent_similarity")
  }
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
                                      output_path_results,
                                      id,
                                      working_dir = wkdir) {
  require(glue)
  require(purrr)
  require(data.table)
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  message(glue(
    "{Sys.time()}: Saving {length(input_files)}",
    " files to {output_dir_db}/{id}"
  ))
  input_files %>%
    purrr::map(
      ~ append_sql_db(
        table = as.data.table(readRDS(.)),
        db_loc = glue("{output_dir_db}/{id}.db"),
        table_name = "results"
      )
    )
  # collect summary metrics from sql db
  collect_summary_metrics(
    input_db_path = glue("{output_dir_db}/{id}.db"),
    outputdir = output_path_results
  )
  # remove SQLite DB
  shell_do("rm {output_dir_db}/{id}.db")
}

#_____________________________________________________________________________
# HHsuite utility functions

slurm_do_msa_trim <- function(
    fname_list,
    out_dir,
    working_dir) {
  require(glue)
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  message(get_time(), " Data dir: ", out_dir)
  for (path in fname_list) {
    fname <- fs::path_ext_remove(basename(path))
    if (file.exists(glue("{out_dir}/{fname}.trimD30.fasta"))) {
      message(glue("{get_time()} {fname} already exists, skipping..."))
      next
    }
    message(glue("{get_time()} {fname} processing..."))
    shell_do(
      glue(
        "hhfilter -i {path}",
        " -o {out_dir}/{fname}.trimD30.a3m -diff 30 && ",
        "reformat.pl -r {out_dir}/{fname}.trimD30.a3m",
        " {out_dir}/temp_{fname}.trimD30.fasta && ",
        "mamba run -n pdmbsR seqkit rename",
        " {out_dir}/temp_{fname}.trimD30.fasta >",
        " {out_dir}/{fname}.trimD30.fasta"
      )
    )
    shell_do(glue("rm {out_dir}/temp_{fname}.trimD30.fasta"))
  }
}

shell_do_convert_a3m_msa_hmm <- function(
    fname_list,
    out_dir,
    working_dir) {
  require(glue)
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  for (path in fname_list) {
    fname <- fs::path_ext_remove(basename(path))
    if (file.exists(glue("{out_dir}/{fname}.hmm"))) {
      message(glue("{get_time()} {fname} already exists, skipping..."))
      next
    }
    message(glue("{get_time()} {fname} processing..."))
    shell_do(
      glue(
        "reformat.pl -r {out_dir}/{fname}.a3m ",
        "{out_dir}/temp_{fname}.fasta &&",
        " mamba run -n pdmbsR seqkit rename ",
        "{out_dir}/temp_{fname}.fasta > ",
        "{out_dir}/{fname}.fasta &&",
        " hmmbuild {out_dir}/{fname}.hmm {out_dir}/{fname}.fasta"
      )
    )
    shell_do(glue("rm {out_dir}/temp_{fname}.fasta"))
  }
}

shell_do_hmmsearch <- function(
    input_paths,
    output_dir,
    db_path,
    working_dir = wkdir) {
  require(glue)
  require(fs)
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  for (hmm_path in input_paths) {
    id <- fs::path_ext_remove(basename(hmm_path))
    if (file.exists(glue("{output_dir}/{id}_hmmsearch.out"))) {
      message(glue("{id} already exists, skipping..."))
      next
    }
    hmmsearch_cmd <- glue(
      "hmmsearch ",
      " --cpu 24",
      " -o {output_dir}/{id}_hmmsearch.out",
      " --tblout {output_dir}/{id}_tblout.txt",
      " --domtblout {output_dir}/{id}_domtblout.txt",
      " --noali",
      " {hmm_path}",
      " {db_path}"
    )
    shell_do(hmmsearch_cmd)
  }
}
