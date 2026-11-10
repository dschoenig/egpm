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

# n.threads <- 4
# task_id <- 1
# task_count <- 200
# chunks.bypass <- c(12, 13)

sam.frac <- 0.005
som.dim <- 50
som.rlen <- 1000
egp.k.som <- 500
egp.k.geo <- 250
egp.max.knots.som <- som.dim^2
egp.max.knots.geo <- egp.k.geo*10
egp.approx <- TRUE
egp.basis <- "gp"
egp.select <- TRUE
overwrite <- TRUE

path.base <- "../"
ls.type <- "1000_4cov_nl"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
mod.type <- "egp_sam0.005_som50"
path.mod <- paste0(path.base, "models/", ls.type, "/")
if(!dir.exists(path.mod)) dir.create(path.mod, recursive = TRUE)

file.par <- paste0(path.ls, "parameters.rds")

parameters <- readRDS(file.par)
# parameters <- parameters[1:2]
ls.total <- nrow(parameters)

if(!overwrite) {
  files <- list.files(path.mod, pattern = mod.type)
  if(length(files > 0)) {
    files.mod <- paste0(path.mod,
                        list.files(path.mod, pattern = mod.type))
    ids <- as.integer(stri_match_last_regex(files.mod, "\\d{4}"))
    parameters <- parameters[!id %in% ids]
  }
}

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

# chunk <- 1
# results <- list()
for(i in chunk) {

  ta <- Sys.time()

  i.step <- which(chunk == i)

  message(paste0("Fitting EGP models for landscape ", parameters[i, id],
                 "/", ls.total, " (",
                 i.step, "/", length(chunk),
                 " in chunk)", " …"))

  results.mod <- readRDS(file = files.res[i.step])
  ls.sam <- results.mod$sample
  som.fit <- results.mod$som


  message("Estimating marginal effect …")

  mod.egp <- results.mod$estimates.noint$gam
  post <- results.mod$estimates.noint$post

  data.bl <- assign_bl_som(ls.sam,
                           som = som.fit,
                           cov.col = c("z1", "z2", "z3", "z4"),
                           id.col = "id")
  lp <-
    evaluate_posterior(model = mod.egp,
                       posterior = post,
                       newdata = ls.sam,
                       id.col = "id",
                       type = "response",
                       progress = FALSE)

  groups.type <-
      ls.sam |>
      ids_by_group(id.col = "id", group.vars = "type")

  id.list.type <- groups.type$ids
  names(id.list.type) <- groups.type$type
  yhat.type <- aggregate_variables(lp,
                      agg.fun = mean,
                      ids = id.list.type,
                      progress = FALSE
                      )

  groups.som <-
      ls.sam[type == "control"] |>
      ids_by_group(id.col = "id", group.vars = "som_bmu")
  id.list.som <- groups.som$ids
  names(id.list.som) <- groups.som$som_bmu
  yhat.som <- aggregate_variables(lp,
                      agg.fun = mean,
                      ids = id.list.som,
                      agg.size = 1e5,
                      progress = FALSE
                      )

  ids.units <- data.bl[,
                         .(id,
                           som_bmu.bl,
                           som_bmu.bl.w)
                         ][,
                           lapply(.SD, unlist), c("id")
                           ][,
                             .(id,
                               som_bmu.bl,
                               som_bmu.bl.w)]
  setkey(ids.units, id)

  # Reweigh baseline SOM units for each group, based on what points they where
  # assigned to
  w.points <-
    lapply(id.list.type,
           \(x) {
                 extract_weights(ids.units[.(x)],
                                 w.col = "som_bmu.bl.w",
                                 by.col = "som_bmu.bl",
                                 standardize = TRUE)
                })
  yhat.bl.type <- reweigh_posterior(yhat.som, w = w.points)

  eff.mar <- arc(yhat.type, yhat.bl.type)

  # Export results

  results.mod$estimates.noint$effects <- eff.mar

  rm(mod.egp, post, eff.mar)


  message("Estimating marginal effect …")

  mod.egp <- results.mod$estimates.int$gam
  post <- results.mod$estimates.int$post

  data.bl <- assign_bl_som(ls.sam,
                           som = som.fit,
                           cov.col = c("z1", "z2", "z3", "z4"),
                           id.col = "id")
  lp <-
    evaluate_posterior(model = mod.egp,
                       posterior = post,
                       newdata = ls.sam,
                       id.col = "id",
                       type = "response",
                       progress = FALSE)

  groups.type <-
      ls.sam |>
      ids_by_group(id.col = "id", group.vars = "type")

  id.list.type <- groups.type$ids
  names(id.list.type) <- groups.type$type
  yhat.type <- aggregate_variables(lp,
                      agg.fun = mean,
                      ids = id.list.type,
                      progress = FALSE
                      )

  groups.som <-
      ls.sam[type == "control"] |>
      ids_by_group(id.col = "id", group.vars = "som_bmu")
  id.list.som <- groups.som$ids
  names(id.list.som) <- groups.som$som_bmu
  yhat.som <- aggregate_variables(lp,
                      agg.fun = mean,
                      ids = id.list.som,
                      agg.size = 1e5,
                      progress = FALSE
                      )

  ids.units <- data.bl[,
                         .(id,
                           som_bmu.bl,
                           som_bmu.bl.w)
                         ][,
                           lapply(.SD, unlist), c("id")
                           ][,
                             .(id,
                               som_bmu.bl,
                               som_bmu.bl.w)]
  setkey(ids.units, id)

  # Reweigh baseline SOM units for each group, based on what points they where
  # assigned to
  w.points <-
    lapply(id.list.type,
           \(x) {
                 extract_weights(ids.units[.(x)],
                                 w.col = "som_bmu.bl.w",
                                 by.col = "som_bmu.bl",
                                 standardize = TRUE)
                })
  yhat.bl.type <- reweigh_posterior(yhat.som, w = w.points)

  eff.mar <- arc(yhat.type, yhat.bl.type)

  # Export results

  results.mod$estimates.int$effects <- eff.mar

  rm(mod.egp, post, eff.mar)

  # results[[i]] <- results.mod
  saveRDS(object = results.mod, file = files.tmp[i.step])

  rm(results.mod) 
  
  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

message("Copying results to final destination …")

for(i in seq_along(files.tmp)) {
  file.copy(files.tmp[i], files.res[i], overwrite = TRUE)
}
