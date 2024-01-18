
library(data.table)
library(mgcv)
library(stringi)

source("utilities.R")

args <- commandArgs(trailingOnly = TRUE)
ls.type <- args[1]
mod.type <- args[2]

ls.type <- "binary_imbalance_high"
mod.type <- "egp_som25"

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
path.mod <- paste0(path.base, "models/", ls.type, "/")
path.results <- paste0(path.base, "results/", ls.type, "/")
if(!dir.exists(path.results)) dir.create(path.results, recursive = TRUE)

file.results <- paste0(path.results, mod.type, "_imbalance.rds")

files.mod <-
  paste0(path.mod, list.files(path.mod)) |>
  stri_subset(regex = paste0(mod.type, "_\\d+.rds$"))

# ids <- as.integer(stri_match_last_regex(files.mod, "\\d+")[,1])

imb.i <- list()

for(i in 1:length(files.mod)) {

  ta <- Sys.time()

  message(paste0("Summarizing results for landscape ", i, " …"))

  mod.res <- readRDS(files.mod[i])

  imb.j <- list()

  mod.para <- mod.res$parameters

  for(j in 1:nrow(mod.para)) {

    message(paste0("Parameter combination ", j, "/", nrow(mod.para)), " …")

    cf.def <- mod.res$models[[j]]$egp.def
    data <- mod.res$sample

    imb.mod <-
      egp_imbalance(data,
                    group = NULL,
                    variables = c("z1", "z2", "z3", "z4"),
                    cf.def = cf.def,
                    js.type = "discrete")

    imb.j[[j]] <- 
      merge(cf.def$groups[, c(cf.def$group.var, cf.def$group.by.c), with = FALSE],
            imb.mod)

  }

  imb.i[[i]] <- rbindlist(imb.j, idcol = "mod.id")

  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

imb <- rbindlist(imb.i, idcol = "landscape")

saveRDS(imb, file.results)

