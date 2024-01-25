# module load StdEnv/2023 gcc/12.3 gdal/3.7.2 geos/3.12.0 python/3.11.5 udunits/2.2.28 r/4.2.2
# module load StdEnv/2023 gcc/12.3 gdal/3.7.2 geos/3.12.0 python/3.11.5 udunits/2.2.28 r/4.3.1

c(
  "mgcv",
  # "brms",
  # "rstan",
  # "RandomFields",
  # "RandomFieldsUtils",
  "sf",
  "stars",
  "lwgeom",
  "spdep",
  "igraph",
  "data.table",
  "mvnfast",
  "colorspace",
  "stringi",
  "ggplot2",
  "patchwork",
  "kohonen",
  "GA",
  "nloptr",
  "memoise",
  "Matrix",
  "MatchIt",
  "marginaleffects",
  "sandwich",
  "tweedie") |>
install.packages()

install.packages(c("rstan", "brms"),
                 repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

data.table::update_dev_pkg()

# install.packages("~/projects/pkgarchive/RandomFieldsUtils_1.2.5.tar.gz",
#                  repos = NULL, type = "source")
# install.packages("~/projects/pkgarchive/RandomFields_3.3.14.tar.gz",
#                  repos = NULL, type = "source")

devtools::install_github("cran/RandomFieldsUtils")
devtools::install_github("cran/RandomFields")
