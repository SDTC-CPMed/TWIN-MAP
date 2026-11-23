# Shared path settings for the R scripts.
#
# Every script sources this file through the short block near its top, which finds this folder and sets
# EARLYT_CODE. Folder locations are read from environment variables of the same name or, when a variable is
# not set, from paths.env in this folder (copy paths.env.example to paths.env and fill it in).
# A script asks only for the folders it uses, through earlyt_path("<NAME>").

local({
  f <- file.path(Sys.getenv("EARLYT_CODE", "."), "paths.env")
  if (!file.exists(f)) return(invisible())
  for (line in trimws(readLines(f, warn = FALSE, encoding = "UTF-8"))) {
    if (!nzchar(line) || startsWith(line, "#") || !grepl("=", line, fixed = TRUE)) next
    key <- trimws(sub("=.*$", "", line))
    val <- gsub('^["\']|["\']$', "", trimws(sub("^[^=]*=", "", line)))
    if (!nzchar(Sys.getenv(key))) do.call(Sys.setenv, setNames(list(val), key))   # environment variables win
  }
})

earlyt_path <- function(name) {
  v <- Sys.getenv(name)
  if (!nzchar(v)) stop(name, " is not set. Copy paths.env.example to paths.env and fill it in, or set the environment variable.", call. = FALSE)
  path.expand(v)
}
