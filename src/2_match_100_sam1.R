args <- commandArgs(trailingOnly = TRUE)

library(mgcv)
library(kohonen)
library(MatchIt)

source("utilities_fcne.R")
source("utilities.R")

n.threads <- as.integer(args[1])

# n.threads <- 4

path.base <- "../"
ls.type <- "100_4cov_nl"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
mod.type <- "match_sam1"
path.mod <- paste0(path.base, "models/", ls.type, "/")
if(!dir.exists(path.mod)) dir.create(path.mod, recursive = TRUE)

file.par <- paste0(path.ls, "parameters.rds")

sam.frac <- 1

parameters <- readRDS(file.par)
# parameters <- parameters[1:2]


# Helper function

extract_lm <- function(x) {
  x.sum <- summary(x)
  x.coef <- x.sum$coefficients["typetreatment",1:2]
  x.ci <- unname(x.coef[1] + x.coef[2] * qnorm(c(0.025, 0.975)))
  dev.expl <- unname(x.sum$r.squared)
  return(list(mean = unname(x.coef[1]), ci = x.ci, dev.expl = dev.expl))
}


for(i in 1:nrow(parameters)) {

  ta <- Sys.time()
  
  message(paste0("Fitting matching models for landscape ", i, " / ", nrow(parameters), " …"))

  results.mod <- list()

  ls.par <- 
    as.list(parameters[i,]) |>
    lapply(unlist)
  file.ls <- paste0(path.ls.data,
                    stri_pad_left(ls.par$id, 4, 0), ".rds")
  file.mod <- paste0(path.mod, mod.type, "_", stri_pad_left(ls.par$id, 4, 0), ".rds")
  
  ls <- readRDS(file.ls)$landscape

  set.seed(ls.par$seed) 
  sam <- sample(1:nrow(ls), round(sam.frac * nrow(ls)))
  ls.sam <- na.omit(ls[sam,])

  mod.lm <- lm(response ~ type, data = ls.sam)
  mod.lmcov <- lm(response ~ type + z1 + z2 + z3 + z4, data = ls.sam)

  matched.cem <- matchit(type ~ z1 + z2 + z3 + z4, method = "cem", data = ls.sam)
  md.cem <- match.data(matched.cem)
  mod.cem <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.cem)

  matched.nn.ps <-
    matchit(type ~ z1 + z2 + z3 + z4, method = "nearest", distance = "glm", data = ls.sam)
  md.nn.ps <- match.data(matched.nn.ps)
  mod.nn.ps <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.ps)

  matching <- 
    list(lm = extract_lm(mod.lm),
         lmcov = extract_lm(mod.lmcov),
         cem = extract_lm(mod.cem),
         nn.ps = extract_lm(mod.nn.ps))

  results.mod <- list()
  for(i in seq_along(matching)) {
    est <- matching[[i]]
    results.mod[[i]] <-
      data.table(id = ls.par$id,
                 method = names(matching)[i],
                 mean = est$mean,
                 q2.5 = est$ci[1],
                 q97.5 = est$ci[2],
                 dev.expl = est$dev.expl 
                 )
  }
  
  results.mod <- rbindlist(results.mod)

  # Export results
  saveRDS(results.mod, file.mod)

  tb <- Sys.time()
  te <- tb-ta
  print(te)

  }

