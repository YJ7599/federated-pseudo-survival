## Install the R packages used by the simulation and plotting scripts.
## Missing packages are installed into the first library in .libPaths(), or into the
## folder given by the environment variable FEDSURV_RLIB if it is set.
## Run from the project root:  Rscript install_packages.R
lib <- Sys.getenv("FEDSURV_RLIB")
if (nzchar(lib)) { dir.create(lib, recursive = TRUE, showWarnings = FALSE); .libPaths(c(lib, .libPaths())) }
pkgs <- c("survival", "splines2", "intccr", "pseudo", "rmutil", "dplyr", "tidyr", "ggplot2",
          "scales", "pda", "jsonlite", "glmtrans", "glmnet", "readr", "data.table", "Rcpp")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing, lib = .libPaths()[1], repos = "https://cloud.r-project.org")
ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
print(data.frame(package = pkgs, installed = ok,
                 version = vapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_, "")),
      row.names = FALSE)
if (!all(ok)) stop("Packages not available: ", paste(pkgs[!ok], collapse = ", "))
