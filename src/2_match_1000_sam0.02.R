args <- commandArgs(trailingOnly = TRUE)

library(mgcv)
library(kohonen)
library(MatchIt)

source("utilities_fcne.R")
source("utilities.R")

n.threads <- as.integer(args[1])
task_id <- as.integer(args[2])
task_count <- as.integer(args[3])
if(length(args) > 3) {
  chunks.bypass <- as.integer(args[4:length(args)])
} else {
  chunks.bypass <- NULL
}

path.base <- "../"
ls.type <- "1000_4cov_nl"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
mod.type <- "match_sam0.02"
path.mod <- paste0(path.base, "models/", ls.type, "/")
if(!dir.exists(path.mod)) dir.create(path.mod, recursive = TRUE)

file.par <- paste0(path.ls, "parameters.rds")

sam.frac <- 0.02
overwrite <- FALSE

parameters <- readRDS(file.par)
# parameters <- parameters[1:2]
ls.total <- nrow(parameters)

row.chunks <- chunk_seq(1, nrow(parameters), ceiling(nrow(parameters) / task_count))
if(!is.null(chunks.bypass)) {
  row.chunks <- lapply(row.chunks, \(x) x[chunks.bypass])
}
chunk <- row.chunks$from[task_id]:row.chunks$to[task_id]


files.tmp <- paste0(paste0(tempdir(), "/", mod.type, "_",
                           stri_pad_left(parameters[chunk, id], 4, 0)),
                    ".rds")
files.res <- paste0(paste0(path.mod, mod.type, "_",
                           stri_pad_left(parameters[chunk, id], 4, 0)),
                    ".rds")

for(i in chunk) {

  ta <- Sys.time()
  
  i.step <- which(chunk == i)
  
  message(paste0("Fitting matching models for landscape ", parameters[i, id],
                 "/", ls.total, " (",
                 i.step, "/", length(chunk),
                 " in chunk)", " …"))

  results.mod <- list()

  ls.par <- 
    as.list(parameters[i,]) |>
    lapply(unlist)
  file.ls <- paste0(path.ls.data,
                    stri_pad_left(ls.par$id, 4, 0), ".rds")

  if(overwrite | !file.exists(files.res[i.step])) {

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

    matched.cem.st <- matchit(type ~ z1 + z2 + z3 + z4,
                           method = "cem", cutpoints = "st",
                           data = ls.sam)
    md.cem.st <- match.data(matched.cem.st)
    mod.cem.st <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights,
                          data = md.cem.st)
    rm(matched.cem.st, md.cem.st)

    matched.cem.sc <- matchit(type ~ z1 + z2 + z3 + z4,
                           method = "cem", cutpoints = "sc",
                           data = ls.sam)
    md.cem.sc <- match.data(matched.cem.sc)
    mod.cem.sc <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights,
                  data = md.cem.sc)
    rm(matched.cem.sc, md.cem.sc)

    matched.cem.fd <- matchit(type ~ z1 + z2 + z3 + z4,
                           method = "cem", cutpoints = "fd",
                           data = ls.sam)
    md.cem.fd <- match.data(matched.cem.fd)
    mod.cem.fd <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights,
                     data = md.cem.fd)
    rm(matched.cem.fd, md.cem.fd)

    matched.nn.ps.nr <-
      matchit(type ~ z1 + z2 + z3 + z4,
              method = "nearest", distance = "glm", replace = FALSE,
              data = ls.sam)
    md.nn.ps.nr <- match.data(matched.nn.ps.nr)
    mod.nn.ps.nr <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.ps.nr)
    rm(matched.nn.ps.nr, md.nn.ps.nr)

    matched.nn.ps.re <-
      matchit(type ~ z1 + z2 + z3 + z4,
              method = "nearest", distance = "glm", replace = TRUE,
              data = ls.sam)
    md.nn.ps.re <- match.data(matched.nn.ps.re)
    mod.nn.ps.re <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.ps.re)
    rm(matched.nn.ps.re, md.nn.ps.re)

    matched.nn.mh.nr <-
      matchit(type ~ z1 + z2 + z3 + z4, replace = FALSE,
              method = "nearest", distance = "mahalanobis",
              data = ls.sam)
    md.nn.mh.nr <- match.data(matched.nn.mh.nr)
    mod.nn.mh.nr <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.mh.nr)
    rm(matched.nn.mh.nr, md.nn.mh.nr)

    matched.nn.mh.re <-
      matchit(type ~ z1 + z2 + z3 + z4, replace = TRUE,
              method = "nearest", distance = "mahalanobis",
              data = ls.sam)
    md.nn.mh.re <- match.data(matched.nn.mh.re)
    mod.nn.mh.re <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.mh.re)
    rm(matched.nn.mh.re, md.nn.mh.re)

    results.mod[["estimates.noint"]] <- 
      list(lm = mod.lm,
           lmcov = mod.lmcov,
           cem.st = mod.cem.st,
           cem.sc = mod.cem.sc,
           cem.fd = mod.cem.fd,
           nn.ps.nr = mod.nn.ps.nr,
           nn.ps.re = mod.nn.ps.re,
           nn.mh.nr = mod.nn.mh.nr,
           nn.mh.re = mod.nn.mh.re)

    rm(mod.lm, mod.lmcov,
       mod.cem.st, mod.cem.sc, mod.cem.fd,
       mod.nn.ps.nr, mod.nn.ps.re,
       mod.nn.mh.nr, mod.nn.mh.re)


    ## MODELS (LANDSCAPE WITH INTERACTIONS) ######################################


    mod.lm <- lm(response.int ~ type, data = ls.sam)
    mod.lmcov <- lm(response.int ~ type + z1 + z2 + z3 + z4, data = ls.sam)

    matched.cem.st <- matchit(type ~ z1 + z2 + z3 + z4,
                           method = "cem", cutpoints = "st",
                           data = ls.sam)
    md.cem.st <- match.data(matched.cem.st)
    mod.cem.st <- lm(response.int ~ type + z1 + z2 + z3 + z4, weights = weights,
                          data = md.cem.st)
    rm(matched.cem.st, md.cem.st)

    matched.cem.sc <- matchit(type ~ z1 + z2 + z3 + z4,
                           method = "cem", cutpoints = "sc",
                           data = ls.sam)
    md.cem.sc <- match.data(matched.cem.sc)
    mod.cem.sc <- lm(response.int ~ type + z1 + z2 + z3 + z4, weights = weights,
                  data = md.cem.sc)
    rm(matched.cem.sc, md.cem.sc)

    matched.cem.fd <- matchit(type ~ z1 + z2 + z3 + z4,
                           method = "cem", cutpoints = "fd",
                           data = ls.sam)
    md.cem.fd <- match.data(matched.cem.fd)
    mod.cem.fd <- lm(response.int ~ type + z1 + z2 + z3 + z4, weights = weights,
                     data = md.cem.fd)
    rm(matched.cem.fd, md.cem.fd)

    matched.nn.ps.nr <-
      matchit(type ~ z1 + z2 + z3 + z4,
              method = "nearest", distance = "glm", replace = FALSE,
              data = ls.sam)
    md.nn.ps.nr <- match.data(matched.nn.ps.nr)
    mod.nn.ps.nr <- lm(response.int ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.ps.nr)
    rm(matched.nn.ps.nr, md.nn.ps.nr)

    matched.nn.ps.re <-
      matchit(type ~ z1 + z2 + z3 + z4,
              method = "nearest", distance = "glm", replace = TRUE,
              data = ls.sam)
    md.nn.ps.re <- match.data(matched.nn.ps.re)
    mod.nn.ps.re <- lm(response.int ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.ps.re)
    rm(matched.nn.ps.re, md.nn.ps.re)

    matched.nn.mh.nr <-
      matchit(type ~ z1 + z2 + z3 + z4, replace = FALSE,
              method = "nearest", distance = "mahalanobis",
              data = ls.sam)
    md.nn.mh.nr <- match.data(matched.nn.mh.nr)
    mod.nn.mh.nr <- lm(response.int ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.mh.nr)
    rm(matched.nn.mh.nr, md.nn.mh.nr)

    matched.nn.mh.re <-
      matchit(type ~ z1 + z2 + z3 + z4, replace = TRUE,
              method = "nearest", distance = "mahalanobis",
              data = ls.sam)
    md.nn.mh.re <- match.data(matched.nn.mh.re)
    mod.nn.mh.re <- lm(response.int ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.mh.re)
    rm(matched.nn.mh.re, md.nn.mh.re)

    results.mod[["estimates.int"]] <- 
      list(lm = mod.lm,
           lmcov = mod.lmcov,
           cem.st = mod.cem.st,
           cem.sc = mod.cem.sc,
           cem.fd = mod.cem.fd,
           nn.ps.nr = mod.nn.ps.nr,
           nn.ps.re = mod.nn.ps.re,
           nn.mh.nr = mod.nn.mh.nr,
           nn.mh.re = mod.nn.mh.re)

    rm(mod.lm, mod.lmcov,
       mod.cem.st, mod.cem.sc, mod.cem.fd,
       mod.nn.ps.nr, mod.nn.ps.re,
       mod.nn.mh.nr, mod.nn.mh.re)

    # Export results

    message("Copying results to final destination …")
    saveRDS(results.mod, files.tmp[i.step])
    file.copy(files.tmp[i.step], files.res[i.step], overwrite = TRUE)
    file.remove(files.tmp[i.step])

    rm(ls, ls.sam)
    rm(results.mod)

  }

  tb <- Sys.time()
  te <- tb-ta
  print(te)

}
