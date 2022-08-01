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
# mod.type <- "egp_sam0.01_som50"

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
terms.noint <- list()
dev.expl.noint <- list()

effects.int <- list()
terms.int <- list()
dev.expl.int <- list()

for(i in ids) {

  ta <- Sys.time()

  message(paste0("Summarizing results for EGP model ", i, " …"))

  mod.res <- readRDS(files.mod[i])
  
  mod.egp <- mod.res$estimates.noint$gam
 
  mod.eff.noint <-
    summary(mod.res$estimates.noint$effects,
            mean, \(x) quantile2(x, probs = c(0.025, 0.975))) |>
    as.data.table() |>
    subset(variable == "treatment")
  mod.eff.noint[, id := i]
  mod.eff.noint <- mod.eff.noint[, .(id, mean, q2.5, q97.5)]

  edf.noint <- NA
  k.noint <- NA
  labels <- NA
  n.smooth <- length(mod.egp$smooth)
  for (j in 1:n.smooth) {
    labels[j] <- mod.egp$smooth[[j]]$label
    edf.noint[j] <-
      sum(mod.egp$edf[mod.egp$smooth[[j]]$first.para:mod.egp$smooth[[j]]$last.para])
    k.noint[j] <- 
      mod.egp$smooth[[j]]$bs.dim
  }

  effects.noint[[i]] <- mod.eff.noint

  terms.noint[[i]] <-
    data.table(id = i,
               term = labels,
               edf = edf.noint,
               k = k.noint)

  dev.expl.noint[[i]] <-
    data.table(id = i,
               dev.expl = (mod.egp$null.deviance - mod.egp$deviance) /
                          mod.egp$null.deviance * 100)


  mod.egp <- mod.res$estimates.int$gam
 
  mod.eff.int <-
    summary(mod.res$estimates.int$effects,
            mean, \(x) quantile2(x, probs = c(0.025, 0.975))) |>
    as.data.table() |>
    subset(variable == "treatment")
  mod.eff.int[, id := i]
  mod.eff.int <- mod.eff.int[, .(id, mean, q2.5, q97.5)]

  edf.int <- NA
  k.int <- NA
  labels <- NA
  n.smooth <- length(mod.egp$smooth)
  for (j in 1:n.smooth) {
    labels[j] <- mod.egp$smooth[[j]]$label
    edf.int[j] <-
      sum(mod.egp$edf[mod.egp$smooth[[j]]$first.para:mod.egp$smooth[[j]]$last.para])
    k.int[j] <- 
      mod.egp$smooth[[j]]$bs.dim
  }

  effects.int[[i]] <- mod.eff.int

  terms.int[[i]] <-
    data.table(id = i,
               term = labels,
               edf = edf.int,
               k = k.int)

  dev.expl.int[[i]] <-
    data.table(id = i,
               dev.expl = (mod.egp$null.deviance - mod.egp$deviance) /
                           mod.egp$null.deviance * 100)

  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

results <-
  list(noint = list(effects = rbindlist(effects.noint),
                    terms = rbindlist(terms.noint),
                    dev.expl = rbindlist(dev.expl.noint)),
       int = list(effects = rbindlist(effects.int),
                  terms = rbindlist(terms.int),
                  dev.expl = rbindlist(dev.expl.int)))


saveRDS(results, file.results)

