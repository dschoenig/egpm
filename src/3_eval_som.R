library(data.table)
library(mgcv)
library(stringi)

source("utilities.R")

args <- commandArgs(trailingOnly = TRUE)
ls.type <- args[1]
mod.type <- args[2]
som.suffix <- args[3]

# ls.type <- "imbalance_high"
# mod.type <- "egp_som25"

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
path.mod <- paste0(path.base, "models/", ls.type, "/")
path.results <- paste0(path.base, "results/", ls.type, "/")
if(!dir.exists(path.results)) dir.create(path.results, recursive = TRUE)

file.results <- paste0(path.results, mod.type, "_", som.suffix, ".rds")

files.mod <-
  paste0(path.mod, list.files(path.mod)) |>
  stri_subset(regex = paste0(mod.type, "_\\d+.rds$"))

ids <- as.integer(stri_match_last_regex(files.mod, "\\d+")[,1])

eval.i <- list()

for(i in ids) {

  ta <- Sys.time()

  message(paste0("Summarizing results for landscape ", i, " …"))

  mod.res <- readRDS(files.mod[i])

  eval.j <- list()

  for(j in seq_along(mod.res$models)) {

    som <- mod.res$models[[j]]$som

    eval.j[[j]] <-
      data.table(quant.error = quantization_error(som),
                 topo.error = topological_error(som),
                 var.expl = variance_explained(som))

  }

  eval.i[[i]] <- rbindlist(eval.j, idcol = "mod.id")

  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

eval <- rbindlist(eval.i, idcol = "landscape")

saveRDS(eval, file.results)
