args <- commandArgs(trailingOnly = TRUE)

library(RandomFields)
library(RandomFieldsUtils)
library(raster)
library(sf)
library(stars)
library(spdep)
library(igraph)
library(data.table)
library(mvnfast)
library(colorspace)
library(stringi)
library(ggplot2)
library(patchwork)
library(kohonen)
library(mgcv)
library(MatchIt)
library(optmatch)

RFoptions(install="no")

source("~/projects/fcne_analysis/src/utilities.R")
source("utilities.R")

path.base <- "../"
ls.type <- "100_4cov_lin"
mod.type <- "egp_25"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
path.models <- paste0(path.base, "models/", ls.type, "/", mod.type, "/")

task_id <- as.integer(args[1])
task_count <- as.integer(args[2])

task_id <- 1
task_count <- 1

# if(length(args) < 3) {
#   n.threads <- c(2,1)
# } else {
#   n.threads <- as.integer(args[3])
# }

file.parameters <- paste0(path.ls, "parameters.rds")


parameters <- readRDS(file.parameters)
parameters <- parameters[1:100,]

# Construct chunk overview
row.chunks <- chunk_seq(1, nrow(parameters), ceiling(nrow(parameters) / task_count))

# Subset parameters
parameters <- parameters[row.chunks$from[task_id]:row.chunks$to[task_id],]
silence <- gc()

estimates <- list()
estimates.regions <- list()
som.quality <- list()

system.time({
for(i in 1:nrow(parameters)) {
  ls.id <- parameters[i, id]

  message(paste0("Evaluating landscape ", ls.id, " …"))

  file.ls <- paste0(path.ls.data, stri_pad_left(ls.id, 4, 0), ".rds")

  ls <- readRDS(file.ls)


  # sam <- sample(1:nrow(ls), 10000)
  sam <- 1:10000
  # som.dim <- 25
  som.dim <- 25
  grid <- somgrid(xdim = som.dim, ydim = som.dim, 
                  topo = "rectangular", 
                  neighbourhood.fct = "gaussian")
  cov.z <- scale(ls[sam, .(z1, z2, z3, z4)],
                 center = TRUE, scale = TRUE)
  som.fit <- som(cov.z,
                 grid = grid, 
                 rlen = 1000,
                 radius = som.dim,
                 init = init_som(cov.z, som.dim, som.dim),
                 mode = "pbatch", 
                 cores = 4,
                 normalizeDataLayers = FALSE)
  som.fit$scale <- list(mean = attr(cov.z, "scaled:center"),
                        sd = attr(cov.z, "scaled:scale"))

  mapped <- 
      ls[, .(z1, z2, z3, z4)] |>
      # land[sam, .(z1, z2, z3, z4)] |>
      scale_data_som(som = som.fit) |>
      embed_som(som = som.fit,
                grid.coord = TRUE)
  ls[,
           `:=`(som_bmu = mapped$bmu[,1],
                som_x = mapped$grid.coordinates$bmu.1[,"x"],
                som_y = mapped$grid.coordinates$bmu.1[,"y"])
           ]
  
  quality <-
    evaluate_embedding(ls[, .(z1, z2, z3, z4)],
                       mapped = as.matrix(ls[, .(som_x, som_y)]),
                       k.max = 200,
                       combined = FALSE)
 
  som.quality[[i]] <-
    data.table(id = ls.id,
               rec = mean(quality$dev.expl^-1)^-1,
               ve = variance_explained(som.fit),
               te = topological_error(som.fit))


  ls[, id := cell]


  mod.lm <- lm(response ~ type + z1 + z2 + z3 + z4,
               data = ls)
    
  mod.lmcov <- lm(response ~ type,
                   data = ls)
    
  mod.gam <- gam(response ~
             type + s(som_x, som_y, bs = "gp", k = 200),
             # type + s(som_x, som_y, bs = "gp", k = 25),
             # type + z1 + z2 + z3 + z4,
             # type,
             data = ls,
             select = TRUE,
             method= "REML",
             optimizer = "efs"
            )
  summary(mod.gam)
  AIC(mod.gam)

  mod.egp <- gam(response ~
             s(x, y, by = type, bs = "gp", k = 200) + s(som_x, som_y, bs = "gp", k = 200),
             # type + s(som_x, som_y, bs = "gp", k = 25),
             # type + z1 + z2 + z3 + z4,
             # type,
             data = ls,
             select = TRUE,
             method= "REML",
             optimizer = "efs"
            )
  summary(mod.egp)
  AIC(mod.egp)

  X <- predict(mod.egp, ls[,], type = "lpmatrix")
  post <- rmvn(1000, coef(mod.egp), vcov(mod.egp, unconditional = TRUE))
  sel <- grep("type", names(coef(mod.egp)))
  post[, -sel] <- 0
  lp <- X %*% t(post)
  id.treat <- which(ls[,]$type == "treatment")
  mean(apply(lp[id.treat,], 2, mean))
  Mode(apply(lp[id.treat,], 2, mean))
  quantile2(apply(lp[id.treat,], 2, mean))

  data.bl <- assign_bl_som(ls,
                           som = som.fit,
                           cov.col = c("z1", "z2", "z3", "z4"),
                           id.col = "id")

  post <- rmvn(1000, coef(mod.egp), vcov(mod.egp, unconditional = TRUE))
  colnames(post) <- names(coef(mod.egp))
  post <- as_draws_matrix(post)

  lp <-
    evaluate_posterior(model = mod.egp,
                       posterior = post,
                       newdata = ls,
                       id.col = "id",
                       type = "link",
                       progress = TRUE)

  groups.type <-
      ls |>
      ids_by_group(id.col = "id", group.vars = "type")

  id.list.type <- groups.type$ids
  names(id.list.type) <- groups.type$type
  yhat.type <- aggregate_variables(lp,
                      agg.fun = E,
                      ids = id.list.type,
                      )

  groups.som <-
      ls |>
      ids_by_group(id.col = "id", group.vars = "som_bmu")
  id.list.som <- groups.som$ids
  names(id.list.som) <- groups.som$som_bmu
  yhat.som <- aggregate_variables(lp,
                      agg.fun = E,
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
  eff.bl.type <- reweigh_posterior(yhat.som, w = w.points)

  eff.arc <- arc(yhat.type, eff.bl.type)
  summary(eff.arc)

  matched.cem <- matchit(type ~ z1 + z2 + z3 + z4, method = "cem", data = ls)
  md.cem <- match.data(matched.cem)
  mod.cem <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.cem)
  # lm(response ~ type, weights = weights, data = md.cem) |>
  # summary(mod.cem)


  matched.nn.ps <-
    matchit(type ~ z1 + z2 + z3 + z4, method = "nearest", distance = "glm", data = ls)
  md.nn.ps <- match.data(matched.nn.ps)
  mod.nn.ps <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.ps)
  # lm(response ~ type, weights = weights, data = md.nn.ps) |>
  # glmmTMB(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.ps) |>
  # summary(mod.nn.ps)


  # matched.opt.ps <-
  #   matchit(type ~ z1 + z2 + z3 + z4, method = "full", distance = "glm", data = ls)
  # md.opt.ps <- match.data(matched.opt.ps)

  # lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.opt.ps) |>
  # # lm(response ~ type, weights = weights, data = md.ps) |>
  # summary()


  results[[i]] <-
    data.table(
               id = ls.id,
               lm = coef(mod.lm)["typetreatment"],
               lmcov = coef(mod.lmcov)["typetreatment"],
               gam = coef(mod.gam)["typetreatment"],
               egp = mean(extract_variable(eff.arc, "treatment")),
               cem = coef(mod.cem)["typetreatment"],
               nn = coef(mod.nn.ps)["typetreatment"]
               )


  regions <-
    bin_cols(ls, c("x", "y"), bin.min = c(0, 0), bin.res = c(50,50), append = TRUE) |>
    ids_by_group(id.col = "id", group.vars = c("x.bin", "y.bin")) 
  regions[, group.label := factor(toupper(letters[1:nrow(regions)]))]
  setnames(regions, "group.label", "region")

  ls <-
  regions[,
          .(region, ids)
          ][, lapply(.SD, \(x) as.numeric(unlist(x))), by = region] |>
  merge(x = ls, y = _, by.x = "id", by.y = "ids")

  groups.region <-
      ls |>
      ids_by_group(id.col = "id", group.vars = "region")

  id.list.region <- groups.region$ids
  names(id.list.region) <- groups.region$region
  yhat.region <- aggregate_variables(lp,
                      agg.fun = E,
                      ids = id.list.region,
                      )

  w.points <-
    lapply(id.list.region,
           \(x) {
                 extract_weights(ids.units[.(x)],
                                 w.col = "som_bmu.bl.w",
                                 by.col = "som_bmu.bl",
                                 standardize = TRUE)
                })
  eff.bl.region <- reweigh_posterior(yhat.som, w = w.points)

  eff.arc <- arc(yhat.region, eff.bl.region)
  
  eff.regions <-
  merge(
        ls[, .(treatment = mean(treatment)), by = region],
        as.data.table(summary(eff.arc, egp = mean)),
        by.x = "region", by.y = "variable")
  eff.regions[, id := ls.id]

  results.regions[[i]] <- eff.regions

  if(i == 10) {
    saveRDS(list(results, results.regions, som.quality), "first10.rds")
  }

}
})


saveRDS(list(results, results.regions, som.quality), "test.egp.rds")
