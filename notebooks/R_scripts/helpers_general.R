# ______________________________________________________________________________
# Utility Functions ----

shell_do <- function(command_string, stdout_path = "", stderr_path = "") {
  inputs <- unlist(stringr::str_split(command_string, " "))
  system2(
    command = inputs[1],
    args = inputs[-1],
    stdout = stdout_path,
    stderr = stderr_path
  )
}

# Function to run any shell command using slurm and future.batchtools
slurm_shell_do <- function(cmd,
                           jobname = glue("slurm-shell-{rand_string()}"),
                           working_dir = wkdir,
                           template_path = glue("{wkdir}/batchtools_templates/batchtools.slurm.tmpl"),
                           memory = "1G",
                           ncpus = 1,
                           walltime = 3600) {
  require(magrittr)
  require(future)
  require(future.batchtools)
  # Initiate future.batchtools backend for parallel processing
  future::plan(
    future.batchtools::batchtools_slurm,
    template = template_path,
    resources = list(
      name = jobname,
      memory = memory,
      ncpus = ncpus,
      walltime = walltime
    )
  )
  job %<-% shell_do(cmd)
}

chunk_func <- function(x, n) {
  split(x, cut(seq_along(x), n, labels = FALSE))
}

# from https://stackoverflow.com/questions/42734547/generating-random-strings
rand_string <- function(characters = 5,
                        numbers = 5,
                        symbols = 0,
                        lowerCase = 0,
                        upperCase = 0) {
  ASCII <- NULL
  if (symbols > 0) {
    ASCII <- c(ASCII, sample(c(33:47, 58:34, 91:96, 123:126), symbols))
  }
  if (numbers > 0) {
    ASCII <- c(ASCII, sample(48:57, numbers))
  }
  if (upperCase > 0) {
    ASCII <- c(ASCII, sample(65:90, upperCase))
  }
  if (lowerCase > 0) {
    ASCII <- c(ASCII, sample(97:122, lowerCase))
  }
  if (characters > 0) {
    ASCII <- c(ASCII, sample(c(65:90, 97:122), characters))
  }
  return(rawToChar(as.raw(sample(ASCII, length(ASCII)))))
}

load_rds_files <- function(input_files) {
  require(purrr)
  require(data.table)
  require(dplyr)
  files_out <- purrr::map(
    input_files, ~ as.data.table(readRDS(.))
  ) %>%
    bind_rows()
  return(files_out)
}


# Downloading UHGP datasets via wget
# Function to download files via wget and submit to slurm
wget_download_slurm <- function(jobname,
                                download_link,
                                slurm_out,
                                output_dir,
                                threads = 1,
                                walltime = "5:00:00",
                                mem_per_cpu = "1G") {
  shell_do(
    glue(
      "sbatch",
      " --job-name={jobname}",
      " --ntasks={threads}",
      " --output={slurm_out}/{jobname}.out",
      " --error={slurm_out}/{jobname}.err",
      " --time={walltime}",
      " --mem-per-cpu={mem_per_cpu}",
      " /central/home/jboktor/slurm_wget.sh",
      " -u {download_link}",
      " -o {output_dir}"
    )
  )
}


# Wget a list of download paths
wget_path_list <- function(download_list,
                           download_dir,
                           working_dir = wkdir) {
  suppressPackageStartupMessages(require(glue))
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(purrr))
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  purrr::map(
    download_list,
    ~ shell_do(glue("wget -P {download_dir} {.}"))
  )
}

shell_do_prodigal <- function(input_genome,
                              output_fasta,
                              addtn_flags = "", # add space before flag
                              working_dir = wkdir) {
  suppressPackageStartupMessages(require(glue))
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  shell_do(
    glue(
      "prodigal{addtn_flags}",
      " -i {input_genome}",
      " -a {output_fasta}"
    )
  )
}

shell_do_prodigal_list <- function(input_genome_list,
                                   output_fasta_list,
                                   addtn_flags = "",
                                   working_dir = wkdir) {
  suppressPackageStartupMessages(require(glue))
  suppressPackageStartupMessages(require(purrr))
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  purrr::map(input_genome_list, ~ shell_do_unzip(.)) %>%
    unlist() %>%
    purrr::map2(
      .x = .,
      .y = output_fasta_list,
      ~ shell_do_prodigal(
        input_genome = .x,
        output_fasta = .y,
        addtn_flags = addtn_flags,
        working_dir = working_dir
      )
    )
  purrr::map(input_genome_list, ~ shell_do_zip(.))
}

# ' Run brute force six frame translation
shell_do_transeq <- function(input_genome_list,
                             output_fasta_list,
                             working_dir = wkdir) {
  require(glue)
  require(purrr)
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  transeq_cmd <- function(input_genome,
                          output_fasta) {
    shell_do(
      glue(
        "conda run -n pdmbsR",
        " transeq {input_genome} {output_fasta}",
        " -frame 6 -trim -clean"
      )
    )
  }
  run <- purrr::map(input_genome_list, ~ shell_do_unzip(.)) %>%
    unlist() %>%
    purrr::map2(
      .x = .,
      .y = output_fasta_list,
      ~ transeq_cmd(
        input_genome = .x,
        output_fasta = .y
      )
    )
  rezip <- purrr::map(input_genome_list, ~shell_do_zip(.))
}


# function to decompress gz files
shell_do_unzip <- function(path) {
  if (grepl(".gz", path)) {
    try(
      {
        shell_do(glue("gunzip {path}"))
      },
      silent = TRUE
    )
    return(gsub(".gz", "", path))
  }
  return(path)
}
# function to compress files
shell_do_zip <- function(path, working_dir = wkdir) {
  if (!grepl("\\.gz$", path)) {
    shell_do(glue("gzip {path}"))
  }
}

shell_do_concat_list <- function(path_list, output_file, working_dir = wkdir) {
  require(glue)
  require(purrr)
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  for (f in path_list) {
    shell_do(glue(
      "cat {f} >> {output_file}"
    ))
  }
}

shell_do_diamond <- function(
    input_faa, output_m8, dmnd_db, wkdir, sensitivity) {
  require(glue)
  source(glue("{wkdir}/notebooks/R-scripts/helpers_general.R"))
  dmnd_cmd <- glue(
    "diamond blastp",
    " -d {dmnd_db}",
    " -q {input_faa}",
    " -o {output_m8}",
    " --max-target-seqs 0",
    "{sensitivity}"
  )
  shell_do(dmnd_cmd)
}

# inverse of logical statement
`%nin%` <- Negate(`%in%`)

# function to count ongoing slurm jobs
count_slurm_jobs <- function(params = c("-u", "jboktor")) {
  queue <- system2(
    command = "squeue",
    args = params,
    stdout = TRUE
  )
  length(queue) - 1
}

#' sleep until a condition is true
wait_until <- function(conditional, interval = 2) {
  # Keep looping until the condition is met
  while (!conditional()) {
    Sys.sleep(interval)
  }
}

check_slurm_overload <- function(njobs = 9999, interval = 2) {
  wait_until(function() {
    count_slurm_jobs() < njobs
  }, interval)
}

slurm_run_alphafold <- function(jobname,
                                slurm_out,
                                input_fasta,
                                output_dir,
                                mode = "multimer",
                                walltime = "1-00:00",
                                use_gpu = "True",
                                mem = "64G",
                                gpus = 1,
                                cpus_per_task = 8) {
  shell_do(
    glue(
      "sbatch",
      " --job-name={jobname}",
      " --output={slurm_out}/{jobname}.out",
      " --error={slurm_out}/{jobname}.err",
      " --time={walltime}",
      " --gres=gpu:{gpus}",
      " --mem={mem}",
      " --cpus-per-task={cpus_per_task}",
      " {src_dir}/shell-scripts/alphafold.sub",
      " -g {use_gpu}",
      " -i {input_fasta}",
      " -o {output_dir}",
      " -m {mode}",
      " -n 2" # number of models to generate
    )
  )
}


#' Function to create SQL databases from data-tables
#' @param table a data-table to save
#' @param table_name name of the table
#' @param db_loc full file path of db to save
create_sql_db <- function(db_loc,
                          table_name,
                          table,
                          overwrite_lg = TRUE,
                          append_lg = FALSE) {
  require(DBI)
  require(RSQLite)
  require(dplyr)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)
  DBI::dbWriteTable(con,
    name = table_name,
    value = table,
    overwrite = overwrite_lg,
    append = append_lg
  )
  DBI::dbListTables(con) %>% print()
  DBI::dbDisconnect(con)
}


#' Function to iteratively append data-tables into a SQL database
#' @param table a data-table to save
#' @param table_name name of the table
#' @param db_loc full file path of db to save
append_sql_db <- function(db_loc,
                          table_name,
                          table) {
  require(DBI)
  require(RSQLite)
  require(dplyr)
  require(purrr)
  require(data.table)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)

  DBI::dbWriteTable(
    con,
    value = table,
    name = table_name,
    overwrite = FALSE,
    append = TRUE
  )
  DBI::dbListTables(con) %>% print()
  DBI::dbDisconnect(con)
}


list_items_to_df <- function(list_items, name_col, value_col) {
  list_items %>%
    bind_rows(.id = "id") %>%
    tidyr::pivot_longer(!id,
      names_to = name_col, values_to = value_col
    ) %>%
    dplyr::select(-id)
}

qc_check_files <- function(df, col) {
  col <- sym(col)
  df %>%
    dplyr::summarize(
      n_eids = n_distinct(ensembl_id),
      n_files = n_distinct(!!col),
      n_files_exist = sum(fs::file_exists(!!col)),
      n_files_not_empty = sum(fs::file_size(!!col) > 0, na.rm = TRUE),
      min_file_size = min(fs::file_size(!!col), na.rm = TRUE),
      max_file_size = max(fs::file_size(!!col), na.rm = TRUE)
    )
}


# # Plotting functions ----
# require(ggplot2)
# require(ggpackets)
# require(ggside)
# require(viridis)
# require(ggpointdensity)

# ggpk_dist <- ggpackets::ggpacket() +
#   geom_histogram(bins = 150, fill = "lightblue", alpha = 0.7) +
#   facet_wrap(~method) +
#   theme_light() +
#   theme(
#     strip.background = element_rect(fill="white", size=1, color="grey"),
#     strip.text.x = element_text(size = 13, color = "black")
#   )

# ggpk_length_vs_sim <- ggpackets::ggpacket() +
#   geom_pointdensity(size = 0.25) +
#   scale_color_viridis(option = "G") +
#   geom_ysidedensity(aes(x = stat(density))) +
#   geom_xsidedensity(aes(y = stat(density))) +
#   labs(
#     x = "Percent Similarity",
#     y = "Peptide Length") +
#   facet_wrap(~ method) +
#   theme_light() +
#   scale_y_log10() +
#   scale_ysidex_continuous(guide = "none") +
#   scale_xsidey_continuous(guide = "none") +
#   theme(
#     ggside.panel.scale.x = .17,
#     ggside.panel.scale.y = .17,
#     strip.background = element_rect(fill="white", size=1, color="grey"),
#     strip.text.x = element_text(size = 13, color = "black"),
#     legend.position = "bottom",
#     legend.key.width = unit(1, 'cm')
#   )


  extract_seq_from_catalog <- function(fasta_header_list,
                                       output_fasta_path, catalog_path) {
    #' Function to extract sequences from a fasta file into a new fasta file
    fasta_headers_tmp <-
      tempfile(fileext = ".txt", tmpdir = "/central/scratch/jbok/tmp")
    write.table(
      as.data.frame(fasta_header_list),
      file = fasta_headers_tmp,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    shell_do(
      glue(
        "conda run -n pdmbsR seqtk subseq {catalog_path}",
        " {fasta_headers_tmp} > {output_fasta_path}"
      )
    )
    unlink(fasta_headers_tmp)
  }


get_time <- function(){
  print(format(Sys.time(), "%Y-%m-%d_%H:%M:%S"))
}

# Function to retrieve FASTA sequence from UniProt
get_fasta_from_UniProt <- function(accession) {
  require(httr)
  url <- paste0("https://www.uniprot.org/uniprot/", accession, ".fasta")
  response <- GET(url)
  fasta <- content(response, "text")
  return(fasta)
}


# defining blastp function
run_blastp_custom <- function(
    input_fasta, output_path, working_dir,
    refdb = glue(
      "{wkdir}/data/interim/refseq_genomes/",
      "genus_sampling/blastdb/2023-06-04_RefSeq-genus"
    ),
    output_cols = glue(
      "qaccver saccver staxids",
      " pident length mismatch gapopen qstart",
      " qend sstart send evalue bitscore"
    )) {
  require(glue)
  require(strex)
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  blastp_cmd <- glue(
    "blastp -query {input_fasta}",
    " -db {refdb}",
    " -outfmt '6 {output_cols}'",
    " -evalue 1",
    " -num_threads 16",
    " -max_target_seqs 100000",
    " -out {output_path}"
  )
  print(blastp_cmd)
  shell_do(blastp_cmd)
}

# GGET Blast ----
blastr <- function(sequences_aa,
                  #  ensembl_id,
                   type = "blastp",
                   db = "refseq_protein",
                   wkdir = wkdir,
                   evalue = 10,
                   output_limit = 1000) {
  require(reticulate)
  require(glue)
  use_condaenv("/home/jboktor/miniconda3/envs/pdmbsR/bin/python",
    required = TRUE
  )
  # reticulate::use_condaenv(condaenv = "pdmbsR", required = TRUE)
  gget <- reticulate::import("gget")
  blast_results <- gget$blast(
    sequences_aa,
    program = type,
    database = db,
    limit = output_limit,
    expect = evalue
  )
  # saveRDS(
  #   blast_results,
  #   glue("{wkdir}/data/interim/blast_results/{ensembl_id}.rds")
  # )
  return(blast_results)
}

plot_msa <- function(fa_path, outfile, wk = wkdir) {
  require(ggmsa)
  require(glue)
  p_msa <- ggmsa(fa_path, char_width = 0.5, seq_name = FALSE) +
    geom_seqlogo() +
    geom_msaBar()
  ggsave(outfile, p_msa,
    width = 15, height = 10, dpi = 300
  )
}

plot_msa_wrapper <- function(path_list, output_list) {
  for (i in 1:length(path_list)) {
    plot_msa(path_list[i], output_list[i])
  }
}

seriate_matrix_rows <- function(
    mat, seriate_method = "OLO", dist_method = "euclidean") {
  order <- mat %>%
    stats::dist(method = dist_method) %>%
    seriation::seriate(method = seriate_method) %>%
    seriation::get_order()
  ranked_order <- rownames(mat)[order]
  return(ranked_order)
}

#' This function accepts a pairwise alignment object from the protr package
#' and returns metrics for percent identity and percent similarity
extract_sw_fident  <- function(aln){
  100 * nchar(aln@pattern) / aln@pattern@unaligned@ranges@width
}
