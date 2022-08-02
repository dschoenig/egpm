library(data.table)
library(mgcv)
library(kohonen)
library(MatchIt)

source("utilities_fcne.R")
source("utilities.R")


args <- commandArgs(trailingOnly = TRUE)
ls.type <- args[1]
mod.type <- args[2]

# ls.type <- "1000_4cov_nl"
# mod.type <- "rf_sam0.01_t1000"

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
path.mod <- paste0(path.base, "models/", ls.type, "/")
path.results <- paste0(path.base, "results/", ls.type, "/")
if(!dir.exists(path.results)) dir.create(path.results, recursive = TRUE)

file.results <- paste0(path.results, mod.type, ".rds")

files.mod <- paste0(path.mod,
                    list.files(path.mod, pattern = mod.type))

ids <- as.integer(stri_match_last_regex(files.mod, "\\d{4}"))

effects.noint <- list()
r.sq.noint <- list()

effects.int <- list()
r.sq.int <- list()

for(i in ids) {

  ta <- Sys.time()

  message(paste0("Summarizing results for random forest model ", i, " …"))

  mod.res <- readRDS(files.mod[i])
  
  mod.eff.noint <-
    as.data.table(mod.res$estimates.noint$effects)

  effects.noint[[i]] <-
    data.table(id = i,
               mean = mod.eff.noint[type == "treatment", yhat] -
                      mod.eff.noint[type == "control", yhat])

  r.sq.noint[[i]] <- data.table(id = i,
                           r.sq = mod.res$estimates.noint$rf$r.squared)

  mod.eff.int <-
    as.data.table(mod.res$estimates.int$effects)

  effects.int[[i]] <-
    data.table(id = i,
               mean = mod.eff.int[type == "treatment", yhat] -
                      mod.eff.int[type == "control", yhat])

  r.sq.int[[i]] <- data.table(id = i,
                           r.sq = mod.res$estimates.int$rf$r.squared)

  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

results <-
  list(noint = list(effects = rbindlist(effects.noint),
                    r.sq = rbindlist(r.sq.noint)),
       int = list(effects = rbindlist(effects.int),
                  r.sq = rbindlist(r.sq.int)))


saveRDS(results, file.results)

