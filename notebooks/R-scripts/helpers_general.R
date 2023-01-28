

shell_do <- function(command_string, stdout_path = "", stderr_path = "") {
  inputs <- unlist(stringr::str_split(command_string, " "))
  system2(
    command = inputs[1],
    args = inputs[-1],
    stdout = stdout_path,
    stderr = stderr_path
  )
}
