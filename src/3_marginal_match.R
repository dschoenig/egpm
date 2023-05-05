library(data.table)
library(MatchIt)
library(marginaleffects)

source("utilities.R")

args <- commandArgs(trailingOnly = TRUE)
ls.type <- args[1]
mod.type <- args[2]

# ls.type <- "imbalance_high"
# mod.type <- "match"

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
path.mod <- paste0(path.base, "models/", ls.type, "/")
path.results <- paste0(path.base, "results/", ls.type, "/")
if(!dir.exists(path.results)) dir.create(path.results, recursive = TRUE)

file.results <- paste0(path.results, mod.type, ".rds")

files.mod <-
  paste0(path.mod, list.files(path.mod)) |>
  stri_subset(regex = paste0(mod.type, "_\\d+.rds$"))

ids <- as.integer(stri_match_last_regex(files.mod, "\\d+")[,1])

params.i <- list()
marginals.i <- list()
dev.expl.i <- list()

# SUBSETTING IDS REMOVE LATER
task_id <- as.integer(args[3])
task_count <- as.integer(args[4])
row.chunks <- chunk_seq(1, length(ids), ceiling(length(ids)/ task_count))
chunk <- row.chunks$from[task_id]:row.chunks$to[task_id]
ids <- chunk

for(i in ids) {

  ta <- Sys.time()

  message(paste0("Summarizing results for model ", i, " …"))

  mod.res <- readRDS(files.mod[i])

  marginals.j <- list()
  dev.expl.j <- list()

  for(j in seq_along(mod.res$models)) {

    if(mod.res$parameters$type[j] == "lm") {
      ls.fit <-
        as.data.table(mod.res$models[[j]]$model$model)
      mod.res$models[[j]]$marginal <-
        avg_comparisons(mod.res$models[[j]]$model,
                        variables = "type",
                        vcov = "HC3",
                        newdata = ls.fit[type == "treatment"])
    }

    marginals.j[[j]] <- 
      as.data.table(mod.res$models[[j]]$marginal) |>
      DT(, .(mean = estimate,
             q2.5 = conf.low,
             q97.5 = conf.high))

    # UPDATE WITH GLM VERSION
    dev.expl.j[[j]] <-
      data.table(dev.expl = unname(summary(mod.res$models[[j]]$model)$r.squared))

  }

  params.i[[i]] <- mod.res$parameters
  params.i[[i]][, mod.id := 1:.N]
  marginals.i[[i]] <- rbindlist(marginals.j, idcol = "mod.id")
  dev.expl.i[[i]] <- rbindlist(dev.expl.j, idcol = "mod.id")

  saveRDS(mod.res, files.mod[i])

  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

params <- rbindlist(params.i)
setcolorder(params.i[[i]], c("landscape", "mod.id"))
marginals <- rbindlist(marginals.i, idcol = "landscape")
setcolorder(marginals, c("landscape", "mod.id"))
dev.expl <- rbindlist(dev.expl.i, idcol = "landscape")
setcolorder(dev.expl, c("landscape", "mod.id"))

results <-
  list(parameters = params,
       marginal = marginals,
       deviance = dev.expl)

saveRDS(results, file.results)


