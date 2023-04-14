# module load StdEnv/2020 gcc/9.3.0 gdal/3.5.1 geos/3.10.2 python/3.10 udunits/2.2.28 r/4.2.2

c(
  "mgcv",
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
  "memoise",
  "Matrix") |>
install.packages()

data.table::update_dev_pkg()

# install.packages("~/projects/pkgarchive/RandomFieldsUtils_1.2.5.tar.gz",
#                  repos = NULL, type = "source")
# install.packages("~/projects/pkgarchive/RandomFields_3.3.14.tar.gz",
#                  repos = NULL, type = "source")

devtools::install_github("cran/RandomFields")
devtools::install_github("cran/RandomFieldsUtils")
