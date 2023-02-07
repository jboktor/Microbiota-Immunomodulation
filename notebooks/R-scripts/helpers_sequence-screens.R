
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
#' @param db a character vector of the database to use for the blast search
#' @param wkdir a character vector of the working directory
#' @param output_limit a numeric vector of the number of hits to return
#' @param method a character vector of the alignment algorithm to use
#' @return a data frame of the blast results
align_fasta_sequences <- function(ensembl_id, sequences_aa, method, refdb_path) {
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

  fasta_tmp <- glue("{tmpdir}/{ensembl_id}.fasta")
  message(Sys.time(), ": Writing fasta file to ", fasta_tmp)
  write.fasta(
    sequences = sequences_aa,
    names = ensembl_id,
    file.out = fasta_tmp,
    open = "w",
    nbchar = 60,
    as.string = FALSE
  )

  refdb_name <- fs::path_ext_remove(basename(refdb_path))
  output_dir <- glue("{wkdir}/data/interim/emboss_alignments")
  alignment_file <- glue("{output_dir}/{ensembl_id}_{refdb_name}.{method}")

  if (!file.exists(alignment_file)) {
    message(Sys.time(), ": Executing ", method, " on ", fasta_tmp)
    cmd <- glue(
      "{method} -asequence {fasta_tmp}",
      " -bsequence {refdb_path}",
      " -gapopen 5.0",
      " -gapextend 0.5",
      " -datafile EBLOSUM62",
      " -outfile '{alignment_file}'"
    )
    message(Sys.time(), ": Executing command:  ", cmd)
    shell_do(cmd)
  } else {
    message(
      Sys.time(),
      ": Skipping ", method, " on ", fasta_tmp, " because output file exists"
    )
  }
  rds_out <- glue("{output_dir}/{ensembl_id}_{refdb_name}.{method}.rds")
  if (!file.exists(rds_out)) {
    message(Sys.time(), ": Parsing ", alignment_file)
    parsed_df <- parse_emboss_alignment_file_vect(alignment_file)
    saveRDS(parsed_df, rds_out)
  } else {
    message(
      Sys.time(),
      ": Skipping parsing \n", alignment_file, "\n output file exists"
    )
  }
  # delete temp fasta file
  unlink(fasta_tmp)
}
