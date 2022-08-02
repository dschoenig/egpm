args <- commandArgs(trailingOnly = TRUE)

library(mgcv)
library(kohonen)
library(MatchIt)
library(ranger)
library(pdp)

source("utilities_fcne.R")
source("utilities.R")

n.threads <- as.integer(args[1])
task_id <- as.integer(args[2])
task_count <- as.integer(args[3])

n.threads <- 4
task_id <- 1
task_count <- 100

path.base <- "../"
ls.type <- "1000_4cov_nl"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
mod.type <- "rf_sam0.01_t1000"
path.mod <- paste0(path.base, "models/", ls.type, "/")
if(!dir.exists(path.mod)) dir.create(path.mod, recursive = TRUE)

file.par <- paste0(path.ls, "parameters.rds")

sam.frac <- 0.01
rf.num.trees <- 1000
rf.min.node.size <- 10
rf.splitrule <- "variance"

parameters <- readRDS(file.par)
# parameters <- parameters[1:2]


row.chunks <- chunk_seq(1, nrow(parameters), ceiling(nrow(parameters) / task_count))
chunk <- row.chunks$from[task_id]:row.chunks$to[task_id]


# chunk <- 1
for(i in chunk) {

  ta <- Sys.time()
  
  message(paste0("Fitting random forest models for landscape ", i, " / ", nrow(parameters), " …"))

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
  ls.sam[, id := cell]
  ls.sam <- ls.sam[, !"cell"]
  setcolorder(ls.sam, "id")

  results.mod[["sample"]] <- ls.sam

  ## RANDOM FOREST #############################################################


  mod.rf <- ranger(response ~ type + z1 + z2 + z3 + z4,
                   data = ls.sam, splitrule = rf.splitrule, num.trees = rf.num.trees,
                   sample.fraction = 1, num.threads = n.threads)
  rf.pd <- partial(mod.rf, pred.var = "type")


  results.mod[["estimates.noint"]] <- list(rf = mod.rf,
                                           effects = rf.pd)

  rm(mod.rf, rf.pd)


  mod.rf <- ranger(response.int ~ type + z1 + z2 + z3 + z4,
                   data = ls.sam, splitrule = rf.splitrule, num.trees = rf.num.trees,
                   sample.fraction = 1)
  rf.pd <- partial(mod.rf, pred.var = "type")


  results.mod[["estimates.int"]] <- list(rf = mod.rf,
                                         effects = rf.pd)

  rm(mod.rf, rf.pd)

  saveRDS(results.mod, file.mod)

  rm(results.mod) 
  
  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

