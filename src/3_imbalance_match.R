library(data.table)
library(MatchIt)
library(stringi)

source("utilities.R")

args <- commandArgs(trailingOnly = TRUE)
ls.type <- args[1]
mod.type <- args[2]

# ls.type <- "tweedie_imbalance_high"
# mod.type <- "match"

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


    matched <- mod.res$models[[j]]$matched

    if(is.null(matched)) {

      data <- mod.res$sample

      id.trt = which(data$type == "treatment")
      id.ref = which(data$type == "reference")

      imb.j[[j]] <-
        imbalance(data,
                  id.trt = id.trt,
                  id.ref = id.ref,
                  type = "raw",
                  variables = c("z1", "z2", "z3", "z4"),
                  js.type = "discrete")

    } else {

      data <-
        match.data(matched,
                   data = mod.res$sample,
                   drop.unmatched = FALSE,
                   weights = ".w")
      data[, .w := .w / sum(.w), by = "type"]

      id.trt = which(data$type == "treatment")
      id.ref = which(data$type == "reference")
      w.trt = data$.w[id.trt]
      w.ref = data$.w[id.ref]

      imb.j[[j]] <-
        imbalance(data,
                  id.trt = id.trt,
                  id.ref = id.ref,
                  w.trt = w.trt,
                  w.ref = w.ref,
                  type = "both",
                  variables = c("z1", "z2", "z3", "z4"),
                  js.type = "discrete")

    }

  }

  imb.i[[i]] <- rbindlist(imb.j, idcol = "mod.id", fill = TRUE)

  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

imb <- rbindlist(imb.i, idcol = "landscape")

saveRDS(imb, file.results)

