

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

# Downloading UHGP datasets via wget
# Function to download files via wget and submit to slurm
wget_download_slurm <- function(jobname,
                                download_link,
                                slurm_out,
                                output_dir) {
    shell_do(
      glue(
        "sbatch",
        " --job-name={jobname}",
        " --output={slurm_out}/{jobname}.out",
        " --error={slurm_out}/{jobname}.err",
        " --time=4-00:00:00",
        " --mem-per-cpu=100G",
        " /central/home/jboktor/slurm_wget.sh",
        " -u {download_link}",
        " -o {output_dir}"
      )
  )
}
