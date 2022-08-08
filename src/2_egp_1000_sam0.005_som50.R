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

files.mod <- paste0(path.mod, mod.type, "_",
                    stri_pad_left(parameters[chunk, id], 4, 0),
                    ".rds")

# chunk <- 1
results <- list()
for(i in chunk) {

  ta <- Sys.time()
  
  message(paste0("Fitting EGP models for landscape ", parameters[i, id],
                 "/", ls.total, " (",
                 which(chunk == i), "/", length(chunk),
                 " in chunk)", " …"))

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

  ## SOM #######################################################################

  message(paste0("Fitting SOM …"))

  grid <- somgrid(xdim = som.dim, ydim = som.dim, 
                  topo = "rectangular", 
                  neighbourhood.fct = "gaussian")
  cov.z <- scale(ls.sam[, .(z1, z2, z3, z4)],
                 center = TRUE, scale = TRUE)

  som.fit <- som(cov.z,
                 grid = grid, 
                 rlen = som.rlen,
                 radius = som.dim,
                 init = init_som(na.omit(cov.z), som.dim, som.dim),
                 mode = "pbatch", 
                 cores = n.threads,
                 normalizeDataLayers = FALSE)
  som.fit$scale <- list(mean = attr(cov.z, "scaled:center"),
                        sd = attr(cov.z, "scaled:scale"))
  # file.som <- paste0(path.mod, "egp_sam0.01_som50_", stri_pad_left(ls.par$id, 4, 0), ".rds")
  # som.fit <- readRDS(file.som)$som

  mapped <- 
      ls.sam[, .(z1, z2, z3, z4)] |>
      scale_data_som(som = som.fit) |>
      embed_som(som = som.fit,
                grid.coord = TRUE)
  ls.sam[,
           `:=`(som_bmu = mapped$bmu[,1],
                som_x = mapped$grid.coordinates$bmu.1[,"x"],
                som_y = mapped$grid.coordinates$bmu.1[,"y"])
           ]

  ls.sam[, id := cell]
  ls.sam <- ls.sam[, !"cell"]
  setcolorder(ls.sam, "id")
  
  results.mod[["sample"]] <- ls.sam
  results.mod[["som"]] <- som.fit

  egp.k.som <- min(nrow(unique(mapped[[1]])), egp.k.som)

  # if(som.eval) {
  #   quality <-
  #     evaluate_embedding(ls.sam[, .(z1, z2, z3, z4)],
  #                        mapped = as.matrix(ls.sam[, .(som_x, som_y)]),
  #                        k.max = egp.k,
  #                        combined = FALSE)
  #   som.quality[[i]] <-
  #     data.table(id = parameters[i, id],
  #                rec = mean(quality$dev.expl^-1)^-1,
  #                ve = variance_explained(som.fit),
  #                te = topological_error(som.fit))
  # }


  ## EGP (LANDSCAPE WITHOUT INTERACTIONS) ######################################

  # if(is.null(egp.k)) {
  #   egp.k <- min(nrow(unique(ls.sam[, .(som_x, som_y)]) - 1),
  #                floor((nrow(ls.sam) * (0.25))))
  # }


  message("Fitting GAM (response without interactions) …")

 if(egp.approx) {
    mod.egp <- bam(response ~
                   s(x, y, by = type, bs = egp.basis, k = egp.k.geo,
                     xt = list(max.knots = egp.max.knots.geo)) +
                   s(som_x, som_y, bs = egp.basis, k = egp.k.som,
                     xt = list(max.knots = egp.max.knots.som)),
                   data = ls.sam,
                   select = TRUE,
                   discrete = TRUE,
                   nthreads = n.threads
                   )
  } else {
    mod.egp <- gam(response ~
                   s(x, y, by = type, bs = egp.basis, k = egp.k.geo,
                     xt = list(max.knots = egp.max.knots.geo)) +
                   s(som_x, som_y, bs = egp.basis, k = egp.k.som,
                     xt = list(max.knots = egp.max.knots.som)),
                   data = ls.sam,
                   select = egp.select,
                   method= "REML",
                   optimizer = "efs"
                  )
  }
  # summary(mod.egp)
  # AIC(mod.egp)
  # summary(mod.egp)
  # AIC(mod.egp)

  
  # Estimate marginal effect

  message("Estimating marginal effect …")

  data.bl <- assign_bl_som(ls.sam,
                           som = som.fit,
                           cov.col = c("z1", "z2", "z3", "z4"),
                           id.col = "id")

  post <- rmvn(1000, coef(mod.egp), vcov(mod.egp, unconditional = TRUE))
  colnames(post) <- names(coef(mod.egp))
  post <- as_draws_matrix(post)

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
      ls.sam |>
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

  results.mod[["estimates.noint"]] <- list(gam = mod.egp,
                                           posterior = post,
                                           effects = eff.mar)

  rm(mod.egp, post, eff.mar)


  ## EGP (LANDSCAPE WITH INTERACTIONS) #########################################

  # if(is.null(egp.k)) {
  #   egp.k <- min(nrow(unique(ls.sam[, .(som_x, som_y)]) - 1),
  #                floor((nrow(ls.sam) * (0.25))))
  # }

  message("Fitting GAM (response with interactions) …")

 if(egp.approx) {
    mod.egp <- bam(response.int ~
                   s(x, y, by = type, bs = egp.basis, k = egp.k.geo,
                     xt = list(max.knots = egp.max.knots.geo)) +
                   s(som_x, som_y, bs = egp.basis, k = egp.k.som,
                     xt = list(max.knots = egp.max.knots.som)),
                   data = ls.sam,
                   select = TRUE,
                   discrete = TRUE,
                   nthreads = n.threads
                   )
  } else {
    mod.egp <- gam(response.int ~
                   s(x, y, by = type, bs = egp.basis, k = egp.k.geo,
                     xt = list(max.knots = egp.max.knots.geo)) +
                   s(som_x, som_y, bs = egp.basis, k = egp.k.som,
                     xt = list(max.knots = egp.max.knots.som)),
                   data = ls.sam,
                   select = egp.select,
                   method= "REML",
                   optimizer = "efs"
                  )
  }
  # summary(mod.egp)
  # AIC(mod.egp)

  
  # Estimate marginal effect

  message("Estimating marginal effect …")

  data.bl <- assign_bl_som(ls.sam,
                           som = som.fit,
                           cov.col = c("z1", "z2", "z3", "z4"),
                           id.col = "id")

  post <- rmvn(1000, coef(mod.egp), vcov(mod.egp, unconditional = TRUE))
  colnames(post) <- names(coef(mod.egp))
  post <- as_draws_matrix(post)

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
      ls.sam |>
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

  results.mod[["estimates.int"]] <- list(gam = mod.egp,
                                         posterior = post,
                                         effects = eff.mar)

  rm(mod.egp, post, eff.mar)

  results[[i]] <- results.mod

  rm(results.mod) 
  
  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

message("Copying results to final destination …")

for(i in seq_along(results)) {
  saveRDS(results[[i]], files.mod[i])
}
