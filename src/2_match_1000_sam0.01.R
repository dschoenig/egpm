args <- commandArgs(trailingOnly = TRUE)

library(mgcv)
library(kohonen)
library(MatchIt)

source("utilities_fcne.R")
source("utilities.R")

n.threads <- as.integer(args[1])

# n.threads <- 4

path.base <- "../"
ls.type <- "1000_4cov_nl"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
mod.type <- "match_sam0.01"
path.mod <- paste0(path.base, "models/", ls.type, "/")
if(!dir.exists(path.mod)) dir.create(path.mod, recursive = TRUE)

file.par <- paste0(path.ls, "parameters.rds")

sam.frac <- 0.01

parameters <- readRDS(file.par)
# parameters <- parameters[1:2]


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
 

  if(!file.exists(file.mod)) {

    ls <- readRDS(file.ls)$landscape

    set.seed(ls.par$seed) 
    sam <- sample(1:nrow(ls), round(sam.frac * nrow(ls)))
    ls.sam <- na.omit(ls[sam,])
    ls.sam[, id := cell]
    ls.sam <- ls.sam[, !"cell"]
    setcolorder(ls.sam, "id")

    results.mod[["sample"]] <- ls.sam


    ## MODELS (LANDSCAPE WITHOUT INTERACTIONS) ###################################


    mod.lm <- lm(response ~ type, data = ls.sam)
    mod.lmcov <- lm(response ~ type + z1 + z2 + z3 + z4, data = ls.sam)

    matched.cem <- matchit(type ~ z1 + z2 + z3 + z4, method = "cem", data = ls.sam)
    md.cem <- match.data(matched.cem)
    mod.cem <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.cem)

    matched.nn.ps <-
      matchit(type ~ z1 + z2 + z3 + z4, method = "nearest", distance = "glm", data = ls.sam)
    md.nn.ps <- match.data(matched.nn.ps)
    mod.nn.ps <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.ps)

    results.mod[["estimates.noint"]] <- 
      list(lm = mod.lm,
           lmcov = mod.lmcov,
           cem = mod.cem,
           nn.ps = mod.nn.ps)

    rm(mod.lm, mod.lmcov,
       matched.cem, md.cem, mod.cem,
       matched.nn.ps, md.nn.ps, mod.nn.ps)


    ## MODELS (LANDSCAPE WITH INTERACTIONS) ######################################


    mod.lm <- lm(response.int ~ type, data = ls.sam)
    mod.lmcov <- lm(response.int ~ type + z1 + z2 + z3 + z4, data = ls.sam)

    matched.cem <- matchit(type ~ z1 + z2 + z3 + z4, method = "cem", data = ls.sam)
    md.cem <- match.data(matched.cem)
    mod.cem <- lm(response.int ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.cem)

    matched.nn.ps <-
      matchit(type ~ z1 + z2 + z3 + z4, method = "nearest", distance = "glm", data = ls.sam)
    md.nn.ps <- match.data(matched.nn.ps)
    mod.nn.ps <- lm(response.int ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.ps)

    results.mod[["estimates.int"]] <- 
      list(lm = mod.lm,
           lmcov = mod.lmcov,
           cem = mod.cem,
           nn.ps = mod.nn.ps)

    rm(mod.lm, mod.lmcov,
       matched.cem, md.cem, mod.cem,
       matched.nn.ps, md.nn.ps, mod.nn.ps)

    # Export results
    saveRDS(results.mod, file.mod)

    rm(ls, ls.sam)
    rm(results.mod)

  }

  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

