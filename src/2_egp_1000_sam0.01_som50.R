args <- commandArgs(trailingOnly = TRUE)

library(mgcv)
library(kohonen)
library(MatchIt)

source("utilities_fcne.R")
source("utilities.R")

n.threads <- as.integer(args[1])
task_id <- as.integer(args[2])
task_count <- as.integer(args[3])

# n.threads <- 4
# task_id <- 1
# task_count <- 100

path.base <- "../"
ls.type <- "1000_4cov_nl"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
mod.type <- "egp_sam0.01_som25"
path.mod <- paste0(path.base, "models/", ls.type, "/")
if(!dir.exists(path.mod)) dir.create(path.mod, recursive = TRUE)

file.par <- paste0(path.ls, "parameters.rds")

sam.frac <- 0.01
som.dim <- 50
som.rlen <- 1000
# som.eval <- TRUE
egp.k <- 200
egp.max.knots <- som.dim^2
egp.approx <- TRUE
egp.basis <- "gp"
egp.select <- TRUE


parameters <- readRDS(file.par)
# parameters <- parameters[1:2]


row.chunks <- chunk_seq(1, nrow(parameters), ceiling(nrow(parameters) / task_count))
chunk <- row.chunks$from[task_id]:row.chunks$to[task_id]


# som.quality <- list()

for(i in chunk) {

  ta <- Sys.time()
  
  message(paste0("Fitting EGP model for landscape ", i, " / ", nrow(parameters), " …"))

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

  results.mod[["som"]] <- som.fit

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


  ## EGP (NO INTERACTIONS) ####################################################

  # if(is.null(egp.k)) {
  #   egp.k <- min(nrow(unique(ls.sam[, .(som_x, som_y)]) - 1),
  #                floor((nrow(ls.sam) * (0.25))))
  # }

  if(egp.approx) {
    mod.egp <- bam(response ~
                   s(x, y, by = type, bs = egp.basis, k = egp.k) +
                   s(som_x, som_y, bs = egp.basis, k = egp.k,
                     xt = list(max.knots = egp.max.knots)),
                     # xt = list(max.knots = 2000)),
                   data = ls.sam,
                   select = TRUE,
                   discrete = TRUE,
                   nthreads = n.threads
                   )
  } else {
    mod.egp <- gam(response ~
                   s(x, y, by = type, bs = egp.basis, k = egp.k) +
                   s(som_x, som_y, bs = egp.basis, k = egp.k,
                     xt = list(max.knots = egp.max.knots)),
                   data = ls.sam,
                   select = egp.select,
                   method= "REML",
                   optimizer = "efs"
                  )
  }
  # summary(mod.egp)
  # AIC(mod.egp)

  ## Extract model info

  dev.expl <- (mod.egp$null.deviance - mod.egp$deviance)/mod.egp$null.deviance * 100

  egp.edf <- NA
  labels <- NA
  n.smooth <- length(mod.egp$smooth)
  for (j in 1:n.smooth) {
    egp.edf[j] <-
      sum(mod.egp$edf[mod.egp$smooth[[j]]$first.para:mod.egp$smooth[[j]]$last.para])
  }
  for (j in 1:n.smooth) {
    labels[j] <- mod.egp$smooth[[j]]$label
  }
  names(egp.edf) <- labels

  
  # Estimate marginal effect

  ls.sam[, id := cell]
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
                       progress = TRUE)

  groups.type <-
      ls.sam |>
      ids_by_group(id.col = "id", group.vars = "type")

  id.list.type <- groups.type$ids
  names(id.list.type) <- groups.type$type
  yhat.type <- aggregate_variables(lp,
                      agg.fun = mean,
                      ids = id.list.type,
                      )

  groups.som <-
      ls.sam |>
      ids_by_group(id.col = "id", group.vars = "som_bmu")
  id.list.som <- groups.som$ids
  names(id.list.som) <- groups.som$som_bmu
  yhat.som <- aggregate_variables(lp,
                      agg.fun = mean,
                      ids = id.list.som,
                      agg.size = 1e5
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

  eff.mar <- extract_variable(arc(yhat.type, yhat.bl.type), "treatment")

  eff.mar <-
    summary(arc(yhat.type, yhat.bl.type), mean, \(x) quantile2(x, c(0.025, 0.975))) |>
    as.data.table() |>
    subset(variable == "treatment")

  # Export results

  results.mod[["estimates"]] <-
    data.table(id = ls.par$id,
               mean = eff.mar$mean,
               q2.5 = eff.mar$q2.5,
               q97.5 = eff.mar$q97.5,
               dev.expl = dev.expl 
               )
  results.mod[["edf"]] <- as.data.table(as.list(c(egp.edf, "k" = egp.k, "id" = ls.par$id)))

  rm(mod.egp, eff.mar, dev.expl, egp.edf)


  ## EGP (INCLUDING INTERACTIONS) ##############################################

  # if(is.null(egp.k)) {
  #   egp.k <- min(nrow(unique(ls.sam[, .(som_x, som_y)]) - 1),
  #                floor((nrow(ls.sam) * (0.25))))
  # }

  if(egp.approx) {
    mod.egp <- bam(response.int ~
                   s(x, y, by = type, bs = egp.basis, k = egp.k) +
                   s(som_x, som_y, bs = egp.basis, k = egp.k,
                     xt = list(max.knots = egp.max.knots)),
                     # xt = list(max.knots = 2000)),
                   data = ls.sam,
                   select = TRUE,
                   discrete = TRUE,
                   nthreads = n.threads
                   )
  } else {
    mod.egp <- gam(response.int ~
                   s(x, y, by = type, bs = egp.basis, k = egp.k) +
                   s(som_x, som_y, bs = egp.basis, k = egp.k,
                     xt = list(max.knots = egp.max.knots)),
                   data = ls.sam,
                   select = egp.select,
                   method= "REML",
                   optimizer = "efs"
                  )
  }
  # summary(mod.egp)
  # AIC(mod.egp)

  ## Extract model info

  dev.expl <- (mod.egp$null.deviance - mod.egp$deviance)/mod.egp$null.deviance * 100

  egp.edf <- NA
  labels <- NA
  n.smooth <- length(mod.egp$smooth)
  for (j in 1:n.smooth) {
    egp.edf[j] <-
      sum(mod.egp$edf[mod.egp$smooth[[j]]$first.para:mod.egp$smooth[[j]]$last.para])
  }
  for (j in 1:n.smooth) {
    labels[j] <- mod.egp$smooth[[j]]$label
  }
  names(egp.edf) <- labels

  
  # Estimate marginal effect

  ls.sam[, id := cell]
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
                       progress = TRUE)

  groups.type <-
      ls.sam |>
      ids_by_group(id.col = "id", group.vars = "type")

  id.list.type <- groups.type$ids
  names(id.list.type) <- groups.type$type
  yhat.type <- aggregate_variables(lp,
                      agg.fun = mean,
                      ids = id.list.type,
                      )

  groups.som <-
      ls.sam |>
      ids_by_group(id.col = "id", group.vars = "som_bmu")
  id.list.som <- groups.som$ids
  names(id.list.som) <- groups.som$som_bmu
  yhat.som <- aggregate_variables(lp,
                      agg.fun = mean,
                      ids = id.list.som,
                      agg.size = 1e5
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

  eff.mar <- extract_variable(arc(yhat.type, yhat.bl.type), "treatment")

  eff.mar <-
    summary(arc(yhat.type, yhat.bl.type), mean, \(x) quantile2(x, c(0.025, 0.975))) |>
    as.data.table() |>
    subset(variable == "treatment")

  # Export results

  results.mod[["estimates.int"]] <-
    data.table(id = ls.par$id,
               mean = eff.mar$mean,
               q2.5 = eff.mar$q2.5,
               q97.5 = eff.mar$q97.5,
               dev.expl = dev.expl 
               )
  results.mod[["edf.int"]] <- as.data.table(as.list(c(egp.edf, "k" = egp.k, "id" = ls.par$id)))

  saveRDS(results.mod, file.mod)

  rm(results.mod)
  
  tb <- Sys.time()
  te <- tb-ta
  print(te)

  }

