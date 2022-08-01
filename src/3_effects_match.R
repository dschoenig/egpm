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
# mod.type <- "match_sam0.01"

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
dev.expl.noint <- list()

effects.int <- list()
dev.expl.int <- list()

extract_lm <- function(x) {
  x.sum <- summary(x)
  x.coef <- x.sum$coefficients["typetreatment",1:2]
  x.ci <- unname(x.coef[1] + x.coef[2] * qnorm(c(0.025, 0.975)))
  dev.expl <- unname(x.sum$r.squared)
  return(list(mean = unname(x.coef[1]), ci = x.ci, dev.expl = dev.expl))
}

for(i in ids) {

  ta <- Sys.time()

  message(paste0("Summarizing results for matching models ", i, " …"))

  mod.res <- readRDS(files.mod[i])
 
  mod.est.noint <- lapply(mod.res$estimates.noint, extract_lm)
  mod.eff.noint <- list()
  for(j in seq_along(mod.est.noint)) {
    est.noint <- mod.est.noint[[j]]
    mod.eff.noint[[j]] <-
      data.table(id = i,
                 method = names(mod.est.noint)[[j]],
                 mean = est.noint$mean,
                 q2.5 = est.noint$ci[1],
                 q97.5 = est.noint$ci[2],
                 dev.expl = est.noint$dev.expl 
                 )
  }
  mod.eff.noint <- rbindlist(mod.eff.noint)

  mod.est.int <- lapply(mod.res$estimates.int, extract_lm)
  mod.eff.int <- list()
  for(j in seq_along(mod.est.int)) {
    est.int <- mod.est.int[[j]]
    mod.eff.int[[j]] <-
      data.table(id = i,
                 method = names(mod.est.int)[[j]],
                 mean = est.int$mean,
                 q2.5 = est.int$ci[1],
                 q97.5 = est.int$ci[2],
                 dev.expl = est.int$dev.expl 
                 )
  }
  mod.eff.int <- rbindlist(mod.eff.int)

  effects.noint[[i]] <- mod.eff.noint[, .(id, method, mean, q2.5, q97.5)]
  dev.expl.noint[[i]] <- mod.eff.noint[, .(id, method, dev.expl)]
  effects.int[[i]] <- mod.eff.int[, .(id, method, mean, q2.5, q97.5)]
  dev.expl.int[[i]] <- mod.eff.int[, .(id, method, dev.expl)]

  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

results <-
  list(noint = list(effects = rbindlist(effects.noint),
                    dev.expl = rbindlist(dev.expl.noint)),
       int = list(effects = rbindlist(effects.int),
                  dev.expl = rbindlist(dev.expl.int)))


saveRDS(results, file.results)

