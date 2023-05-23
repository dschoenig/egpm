args <- commandArgs(trailingOnly = TRUE)

library(MatchIt)
library(marginaleffects)
library(glmmTMB)

source("utilities.R")

ls.type <- args[1]
mod.type <- args[2]
resp.type <- args[3]
sam.frac <- as.numeric(args[4])
task_id <- as.integer(args[5])
task_count <- as.integer(args[6])
overwrite <- as.logical(args[7])
if(is.na(overwrite)) overwrite <- TRUE

# ls.type <- "binary_imbalance_high"
# mod.type <- "match"
# resp.type <- "binary"
# sam.frac <- 0.01
# task_id <- 60
# task_count <- 1000
# overwrite <- TRUE

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
path.mod <- paste0(path.base, "models/", ls.type, "/")
if(!dir.exists(path.mod)) dir.create(path.mod, recursive = TRUE)

file.par <- paste0(path.ls, "parameters.rds")
file.log <- paste0(path.mod, mod.type, ".log")

parameters <- readRDS(file.par)
ls.total <- nrow(parameters)

row.chunks <- chunk_seq(1, nrow(parameters), ceiling(nrow(parameters) / task_count))
chunk <- row.chunks$from[task_id]:row.chunks$to[task_id]

if(!overwrite & file.exists(file.log)) {
  simulated <- as.integer(read.table(file.log)[,1])
} else {
  simulated <- integer(0)
}

chunk <- setdiff(chunk, simulated)

files.tmp <- paste0(tempdir(), "/", mod.type, "_",
                    parameters[id %in% chunk, file.name],
                    ".rds")

files.res <- paste0(path.mod, mod.type, "_",
                    parameters[id %in% chunk, file.name],
                    ".rds")

# chunk <- 5
for(i in chunk) {

  ta <- Sys.time()

  i.step <- which(chunk == i)
  
  message(paste0("Fitting models for landscape ", i,
                 "/", ls.total, " (",
                 i.step, "/", length(chunk),
                 " in chunk)", " …"))

  results.ls <- list()

  ls.par <- 
    as.list(parameters[id == i,]) |>
    lapply(unlist)

  file.ls <- ls.par$file.path

  ls <- readRDS(file.ls)$landscape

  set.seed(ls.par$seed) 
  sam <- sample(1:nrow(ls), round(sam.frac * nrow(ls)))
  ls.sam <- na.omit(ls[sam,])

  mod.para <-
    rbind(CJ(landscape = ls.par$id,
             sam.frac,
             type = "glm",
             mod.cov = factor(c("exclude", "include", "interact")),
             trt.int = c(FALSE, TRUE),
             sorted = FALSE),
          CJ(landscape = ls.par$id,
             sam.frac,
             type = "cem",
             cutpoints = factor(c("st", "sc", "fd")),
             mod.cov = factor(c("exclude", "include", "interact")),
             trt.int = c(FALSE, TRUE),
             sorted = FALSE),
          CJ(landscape = ls.par$id,
             sam.frac,
             type = "nearest",
             distance = factor(c("glm", "mahalanobis")),
             replace = c(FALSE, TRUE),
             mod.cov = factor(c("exclude", "include", "interact")),
             trt.int = c(FALSE, TRUE),
             sorted = FALSE),
          fill = TRUE)

  mod.para[, type := factor(type)]

  mod.fam <-
    switch(resp.type,
           normal = gaussian(link = "identity"),
           binary = binomial(link = "logit"),
           tweedie = tweedie(link = log))

  results.mod <- list()

  for(j in 1:nrow(mod.para)) {

    message(paste0("Parameter combination ", j, "/", nrow(mod.para)), " …")

    ls.fit <- copy(ls.sam)
  
    if(mod.para$trt.int[j] == FALSE) {

      mod.form <-
        switch(as.character(mod.para$mod.cov[j]),
               exclude = response ~ type,
               include = response ~ type + z1 + z2 + z3 + z4,
               interact = response ~ type + 
                            z1 + z2 + z3 + z4 +
                            z1:z2 + z1:z3 + z1:z4 + z2:z3 + z2:z4 + z3:z4)
    } else {
      mod.form <-
        switch(as.character(mod.para$mod.cov[j]),
               exclude = response ~ type,
               include = response ~ type * (z1 + z2 + z3 + z4),
               interact = response ~ type * ( 
                            z1 + z2 + z3 + z4 +
                            z1:z2 + z1:z3 + z1:z4 + z2:z3 + z2:z4 + z3:z4))
    }

    match.form <-
      switch(as.character(mod.para$mod.cov[j]),
             exclude = type ~ z1 + z2 + z3 + z4,
             include = type ~ z1 + z2 + z3 + z4,
             interact = type ~
                          z1 + z2 + z3 + z4 +
                          z1:z2 + z1:z3 + z1:z4 + z2:z3 + z2:z4 + z3:z4)

    if(mod.para$type[j] == "glm") {

      mod.glm <- NULL
      while(is.null(mod.glm)) {
        try({
          mod.glm <-
            glmmTMB(mod.form,
                    family = mod.fam,
                    data = ls.fit)
        })
        if(is.null(mod.glm)) message("Model failed. Trying again …")
      }

      marginal <-
        avg_comparisons(mod.glm,
                        variables = "type",
                        vcov = vcov(mod.glm, full = TRUE),
                        newdata = ls.fit[type == "treatment"]) |>
        as.data.table()
      marginal[, poly := NA]
   
      if(ls.fit[type == "treatment", length(unique(poly))] > 1) {
        marginal <-
          avg_comparisons(mod.glm,
                          variables = "type",
                          by = "poly",
                          vcov = vcov(mod.glm, full = TRUE),
                          newdata = ls.fit[type == "treatment"][order(-poly)]) |>
          as.data.table() |>
          rbind(marginal, fill = TRUE)
      }

      marginal <-
        marginal[order(poly, na.last = FALSE),
                 .(poly, estimate, std.error, conf.low, conf.high)]

      results.mod[[j]] <-
        list(marginal = marginal)

      rm(mod.glm, marginal)

    } else {

      if(mod.para$type[j] == "cem") {
        matched <-
          matchit(match.form,
                  method = as.character(mod.para$type[j]),
                  cutpoints = as.character(mod.para$cutpoints[j]),
                  data = ls.fit)
      }
      if(mod.para$type[j] == "nearest") {
        matched <-
          matchit(match.form,
                  method = as.character(mod.para$type[j]),
                  distance = as.character(mod.para$distance[j]),
                  replace = mod.para$replace[j],
                  data = ls.fit)
      }

      ls.match <- match.data(matched)

      mod.match <- NULL
      while(is.null(mod.match)) {
        try({
          mod.match <-
            glmmTMB(mod.form,
                    family = mod.fam,
                    data = ls.match,
                    weights = weights)
        })
        if(is.null(mod.match)) message("Model failed. Trying again …")
      }


      marginal <-
        avg_comparisons(mod.match,
                        variables = "type",
                        vcov = vcov(mod.match, full = TRUE),
                        wts = "weights",
                        newdata = ls.match[type == "treatment"]) |>
        as.data.table()
      marginal[, poly := NA]
   
      if(ls.fit[type == "treatment", length(unique(poly))] > 1) {
        marginal <-
          avg_comparisons(mod.match,
                          variables = "type",
                          by = "poly",
                          vcov = vcov(mod.match, full = TRUE),
                          wts = "weights",
                          newdata = ls.match[type == "treatment"][order(-poly)]) |>
          as.data.table() |>
          rbind(marginal, fill = TRUE)
      }

      marginal <-
        marginal[order(poly, na.last = FALSE),
                 .(poly, estimate, std.error, conf.low, conf.high)]

      results.mod[[j]] <-
        list(matched = matched,
             marginal = marginal)

      rm(matched,
         ls.match,
         mod.match,
         marginal)
      gc()
    }
  }

  # mars <- numeric(0)
  # for(k in seq_along(results.mod)) {
  #   mars[k] <-
  #     results.mod[[k]]$marginal
  # }
  
  results.ls <-
    list(parameters = mod.para,
         models = results.mod)


  saveRDS(object = results.ls, file = files.tmp[i.step])

  rm(results.mod) 
  gc()
  
  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

message("Copying results to final destination …")

for(i in seq_along(files.tmp)) {
  file.copy(files.tmp[i], files.res[i], overwrite = TRUE)
}

paste0('echo "',
       paste0(chunk, collapse = "\n"),
       '" >> ',
       file.log) |>
system()
