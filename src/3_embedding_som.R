library(data.table)
library(mgcv)
library(stringi)

source("utilities.R")

args <- commandArgs(trailingOnly = TRUE)
ls.type <- args[1]
mod.type <- args[2]
som.suffix <- args[3]

# ls.type <- "noeff_imbalance_high"
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

# ids <- as.integer(stri_match_last_regex(files.mod, "\\d+")[,1])

eval.i <- list()

for(i in 1:length(files.mod)) {

  ta <- Sys.time()

  message(paste0("Summarizing results for landscape ", i, " …"))

  mod.res <- readRDS(files.mod[i])

  eval.j <- list()

  mod.para <- mod.res$parameters

  for(j in 1:nrow(mod.para)) {

    message(paste0("Parameter combination ", j, "/", nrow(mod.para)), " …")

    if(j > 1) {
      som.same <-
        mod.para[c(j-1, j),
                 all(unlist(lapply(.SD, \(x) x[1] == x[2]))),
                 .SDcols = c("som.topology", "som.dim", "som.epochs"),
                 nomatch = NULL]
    } else {
      som.same <- FALSE
    }

    if(som.same) {
      message("Same as previous SOM …")
      eval.j[[j]] <- eval.j[[j-1]]
    } else {
      message("Evaluating SOM …")
      som <- mod.res$models[[j]]$som
      mod <- mod.res$models[[j]]$model$model
      units.diff <- unit_diff(unit.var = som$unit.classif,
                              group.var = mod$type)
      eval.j[[j]] <-
        data.table(quant.error = quantization_error(som),
                   topo.error = topological_error(som),
                   var.expl = variance_explained(som),
                   units.ov_frac = units.diff$ov_frac,
                   units.js_div = units.diff$js_div)
    }

  }

  eval.i[[i]] <- rbindlist(eval.j, idcol = "mod.id")

  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

eval <- rbindlist(eval.i, idcol = "landscape")

saveRDS(eval, file.results)

