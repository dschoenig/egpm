# module load StdEnv/2020 gcc/9.3.0 gdal/3.4.3 geos/3.10.2 python/3.10 udunits/2.2.28 r/4.2.1

install.packages(c(
                   "bayesplot",
                   "colorspace",
                   "data.table",
                   "devtools",
                   "doParallel",
                   "dplyr",
                   "ggdist",
                   "ggpattern",
                   "ggplot2",
                   "igraph",
                   "KernSmooth",
                   "kohonen",
                   "MatchIt",
                   "mgcv",
                   "mvnfast",
                   "patchwork",
                   "pdp",
                   "posterior",
                   "ranger",
                   "RandomFields",
                   "RandomFieldsUtils",
                   "raster",
                   "sf",
                   "spdep",
                   "stars",
                   "stringi"
                   ))
                   
 devtools::install_github("cran/RandomFields")
