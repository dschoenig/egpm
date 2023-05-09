library(data.table)

source("utilities.R")

args <- commandArgs(trailingOnly = TRUE)
ls.type <- args[1]

# ls.type <- "imbalance_high"

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")

file.results <- paste0(path.ls, "summary.rds")

file.par <- paste0(path.ls, "parameters.rds")

parameters <- readRDS(file.par)

ids <- parameters$id

results.l <- list()

for(i in ids) {

  if(file.exists(parameters[id == i, file.path])) {

    ta <- Sys.time()

    message(paste0("Summarizing results for landscape ", i, " …"))
    
    ls <- readRDS(parameters[id == i, file.path])
   
    imb <- ls$optim[names(ls$optim) %like% "imbalance"]

    results.l[[i]] <- 
      data.table(imbalance = list(imb),
                 imbalance.mean = mean(imb),
                 area.prop = ls$optim["area.prop"])

    tb <- Sys.time()
    te <- tb-ta
    print(te)

  }
}

results <- rbindlist(results.l, idcol = "id")

saveRDS(results, file.results)


