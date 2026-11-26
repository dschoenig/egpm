# module load StdEnv/2026 gcc/15.3.0 gdal/3.13.3 geos/3.14.1 proj/9.9.0 udunits/2.2.28 python/3.15.2 r/4.6.1

# Required for all scripts up to 3_*.R
pkgs.cran <- c(
  # utilities.R (sourced by all scripts)
  "data.table", "GA", "ggplot2", "kohonen", "lwgeom", "Matrix", "memoise",
  "mgcv", "mvnfast", "sf", "spdep", "stars", "stringi", "tweedie",
  # matching and model checks
  "DHARMa", "marginaleffects", "MatchIt"
)

pkgs.installed <- rownames(installed.packages())

pkgs.missing <- setdiff(pkgs.cran, pkgs.installed)
if (length(pkgs.missing) > 0) install.packages(pkgs.missing)

# Optional: only needed for the later scripts (4_*.R onwards)
install.optional <- TRUE

pkgs.optional <- c(
  # parallel processing in 5_comparison_*.R and 5_performance_*.R
  "doFuture", "foreach", "future", "progressr",
  # 6_plot_*.R and 6_table_*.R
  "ggdist", "gt", "scico"
)

if (install.optional) {
  pkgs.missing <- setdiff(pkgs.optional, pkgs.installed)
  if (length(pkgs.missing) > 0) install.packages(pkgs.missing)
}

# Installs the development version of data.table only if it is newer.
data.table::update_dev_pkg()

# RandomFields and RandomFieldsUtils are archived on CRAN. On R 4.6 their
# sources fail to compile, and RFsimulate() segfaults without a patch, so they
# are built by a script. Needs internet: on the cluster, run on a login node.
# Skipped if RFsimulate() already works; tested in a separate R process, so an
# unpatched install that segfaults does not take down this session.
rf.ok <- system2("Rscript",
                 c("-e", shQuote(paste("suppressMessages(library(RandomFields));",
                                       "invisible(RFsimulate(RMexp(), x = 1:10))"))),
                 stdout = FALSE, stderr = FALSE) == 0
if (!rf.ok) {
  rf.install <- c("randomfields_install/install.sh",
                  "src/randomfields_install/install.sh")
  rf.install <- rf.install[file.exists(rf.install)][1]
  stopifnot("Run from the project root or src/" = !is.na(rf.install))
  stopifnot("RandomFields installation failed" =
              system2("bash", rf.install) == 0)
}
