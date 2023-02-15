#______________________________________________________________________________
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

chunk_func <- function(x, n){
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

# Function to run any shell command using slurm and future.batchtools
slurm_shell_do <- function(cmd,
                           jobname = glue("slurm-shell-{rand_string()}"),
                           working_dir = wkdir,
                           memory = "1G",
                           ncpus = 1,
                           walltime = 3600) {
  require(magrittr)
  require(future)
  require(future.batchtools)
  # Initiate future.batchtools backend for parallel processing
  future::plan(
    future.batchtools::batchtools_slurm,
    template = glue("{working_dir}/batchtools_templates/batchtools.slurm.tmpl"),
    resources = list(
      name = jobname,
      memory = memory,
      ncpus = ncpus,
      walltime = walltime
    )
  )
  job %<-% shell_do(cmd)
}

# Downloading UHGP datasets via wget
# Function to download files via wget and submit to slurm
wget_download_slurm <- function(jobname,
                                download_link,
                                slurm_out,
                                output_dir,
                                walltime = "5:00:00",
                                mem_per_cpu = "5G") {
  shell_do(
    glue(
      "sbatch",
      " --job-name={jobname}",
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
                              working_dir = wkdir) {
  suppressPackageStartupMessages(require(glue))
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  shell_do(
    glue(
      "prodigal -i {input_genome}",
      " -a {output_fasta}"
    )
  )
}

shell_do_prodigal_list <- function(input_genome_list,
                                   output_fasta_list,
                                   working_dir = wkdir) {
  suppressPackageStartupMessages(require(glue))
  suppressPackageStartupMessages(require(purrr))
  source(glue("{working_dir}/notebooks/R-scripts/helpers_general.R"))
  purrr::map(
    input_genome_list,
    ~ shell_do_unzip(.)
  ) %>%
  unlist() %>% 
  purrr::map2(
    .x = .,
    .y = output_fasta_list,
    ~ shell_do_prodigal(
      input_genome = .x,
      output_fasta = .y,
      working_dir = working_dir
    )
  )
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

