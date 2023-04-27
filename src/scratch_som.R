# library(RandomFields)
# library(RandomFieldsUtils)
# library(raster)
# library(sf)
# library(stars)
# library(spdep)
# library(igraph)
# library(data.table)
# library(mvnfast)
# library(colorspace)
# library(stringi)
# library(ggplot2)
# library(patchwork)
# library(kohonen)
# library(MatchIt)


# source("utilities_fcne.R")
source("utilities.R")

ls.dim <- 100
ls.imbalance <- 0.1
parallel <- FALSE
n <- 1000

set.seed(19010511) # Rose Ausländer
ls.seeds <- round(runif(n, 0, ls.imbalance) * 1e8)
cov.effect.mu <- c(-2, 2)
cov.int.effect.mu <- c(-2, 2)
cov.effect.range <- c(0.5, 2)
cov.int.effect.range <- c(0.5, 2)
z1.mix.w <- runif(n, 0.2, 0.5)
z2.mix.w <- runif(n, 0.2, 0.5)
z3.mix.prop <- runif(n, 0.2, 0.5)
z4.mix.w <- runif(n, 0.2, 0.5)
opt.iter <- ifelse(ls.imbalance <= 0.3, 250, 500)
opt.run <- ifelse(ls.imbalance <= 0.3, 100, 200)
parameters <-
  data.table(
             id = 1:n,
             seed = sample(1:1e8, n),
             x.dim = ls.dim,
             y.dim = ls.dim,
             # Imbalance
             areas.imbalance = ls.imbalance,
             # Effects
             treatment.eff.mean = 1,
             z1.effect.type = sample(c("sigmoid", "minimum", "unimodal", "bimodal"), n,
                                     replace = TRUE),
             z1.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]),
             z1.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             z2.effect.type = sample(c("sigmoid", "minimum", "unimodal", "bimodal"), n,
                                     replace = TRUE),
             z2.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]),
             z2.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             z3.effect.type = sample(c("sigmoid", "minimum", "unimodal", "bimodal"), n,
                                     replace = TRUE),
             z3.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]),
             z3.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             z4.effect.type = sample(c("sigmoid", "minimum", "unimodal", "bimodal"), n,
                                     replace = TRUE),
             z4.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]),
             z4.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             int12.effect.range = runif(n, cov.int.effect.range[1], cov.int.effect.range[2]) ,
             int12.effect.mu = runif(n, cov.int.effect.mu[1], cov.int.effect.mu[2]),
             int12.effect.nuclei = sample(20:50, n, replace = TRUE),
             int13.effect.range = runif(n, cov.int.effect.range[1], cov.int.effect.range[2]) ,
             int13.effect.mu = runif(n, cov.int.effect.mu[1], cov.int.effect.mu[2]),
             int13.effect.nuclei = sample(20:50, n, replace = TRUE),
             int14.effect.range = runif(n, cov.int.effect.range[1], cov.int.effect.range[2]) ,
             int14.effect.mu = runif(n, cov.int.effect.mu[1], cov.int.effect.mu[2]),
             int14.effect.nuclei = sample(20:50, n, replace = TRUE),
             int23.effect.range = runif(n, cov.int.effect.range[1], cov.int.effect.range[2]) ,
             int23.effect.mu = runif(n, cov.int.effect.mu[1], cov.int.effect.mu[2]),
             int23.effect.nuclei = sample(20:50, n, replace = TRUE),
             int24.effect.range = runif(n, cov.int.effect.range[1], cov.int.effect.range[2]) ,
             int24.effect.mu = runif(n, cov.int.effect.mu[1], cov.int.effect.mu[2]),
             int24.effect.nuclei = sample(20:50, n, replace = TRUE),
             int34.effect.range = runif(n, cov.int.effect.range[1], cov.int.effect.range[2]) ,
             int34.effect.mu = runif(n, cov.int.effect.mu[1], cov.int.effect.mu[2]),
             int34.effect.nuclei = sample(20:50, n, replace = TRUE),
             # Residual variation
             e.exp.var = 0.5,
             e.exp.scale = runif(n, 0.01*ls.dim, ls.dim),
             e.nug.var = 0.5,
             # Parameters for generating functions
             treatment.mat.nu = 1,
             treatment.mat.scale = runif(n, 0.01*ls.dim, 0.5*ls.dim),
             treatment.mat.var = 1,
             treatment.damp.type = "asymmetric",
             treatment.damp.infl = 0.5,
             treatment.damp.scale = runif(n, 1/(2*pi), 1/(0.5*pi)),
             mix.alpha = runif(n, 0.5, 1.5),
             z1.fbm.alpha = runif(n, 0.5, 1.5),
             z1.mix.w = z1.mix.w,
             z2.grad.phi = runif(n, 0, 2*pi),
             z2.mix.w = z2.mix.w,
             z3.dist.n = sample(5:10, n, replace = TRUE),
             z3.dist.acc = ls.dim/100,
             z3.mix.prop = z3.mix.prop,
             z4.seg.n = sample(10:20, n, replace = TRUE),
             z4.mat.nu = runif(n, 1, 1.5),
             z4.mat.scale = runif(n, 0.1, 10),
             z4.mix.w = z4.mix.w,
             # Parameters for defining treatment and reference areas
             score.type = "mahalanobis",
             areas.seg.seed = 3,
             areas.score.sam = 1e5,
             areas.imbalance.tol = 0,
             areas.area.prop = runif(n, 0.35, 0.65),
             areas.area.tol = 0.15,
             areas.area.exact = FALSE,
             areas.seg.res = 0.05 * ls.dim,
             areas.seg.min.dist = 0.025 * ls.dim,
             areas.seg.min.area = (ls.dim/10)^2,
             areas.seg.even = runif(n, 1, 2),
             areas.seg.prec = 5e-4 * ls.dim,
             areas.min.bound.dist = 0,
             areas.opt.imp.imb = 5,
             areas.opt.imp.even = 1,
             areas.opt.imp.area = 1,
             areas.opt.imb.agg = mean,
             areas.opt.pop = 250,
             areas.opt.prec = 1e-3,
             areas.opt.pcrossover = 0.9,
             areas.opt.pmutation = 0.5,
             areas.opt.max.iter = opt.iter,
             areas.opt.run = opt.run,
             areas.opt.parallel = parallel,
             areas.opt.fine = TRUE,
             areas.opt.fine.max.iter = 100,
             areas.opt.fine.constr = TRUE,
             areas.opt.fine.tol = 1e-6,
             areas.opt.cache = TRUE
             )


# attach(as.list(parameters[id == 1]))
# detach(as.list(parameters[id == 1]))

# ls.par <-
# parameters[z1.effect.type != z2.effect.type &
#            z1.effect.type != z3.effect.type &
#            z1.effect.type != z4.effect.type &
#            z2.effect.type != z3.effect.type &
#            z2.effect.type != z4.effect.type &
#            z3.effect.type != z4.effect.type][1] |>
#          as.list()

ls.par <- as.list(parameters[id == 1])
ls.par <- lapply(ls.par, unlist)
ls.par$areas.imbalance <- 0.3
ls.par$areas.opt.imb.agg <- ls.par$areas.opt.imb.agg[[1]]
ls.par$areas.opt.parallel <- 4

attach(ls.par)
detach(ls.par)

system.time({
  ls <- do.call(generate_landscape_4cov_nl, ls.par)
})
# ls <- readRDS("landscape.rds")


plot_landscape_4cov_nl(ls$landscape, select = "all")
plot_functions_4cov_nl(ls$fun)

ls$landscape[, mean(response), by = type]


som.egp2 <-
      egp_som(ls$landscape,
              topo = "hexagonal",
              x.dim = 10,
              y.dim = 10,
              epochs = 1000,
              vars = c("z1", "z2", "z3", "z4"),
              parallel = 4)

som.egp <-
      egp_som(ls$landscape,
              topo = "hexagonal",
              x.dim = 10,
              y.dim = 10,
              epochs = 1000,
              vars = c("z1", "z2", "z3", "z4"),
              parallel = 4)

som.egp$grid.nb == som.egp2$grid.nb


som.egp$grid.nb[1,]
som.egp2$grid.nb[1,]

d_cohen(ls$landscape[type == "reference", z1], ls$landscape[type == "treatment", z1])
d_cohen(ls$landscape[type == "reference", z2], ls$landscape[type == "treatment", z2])
d_cohen(ls$landscape[type == "reference", z3], ls$landscape[type == "treatment", z3])
d_cohen(ls$landscape[type == "reference", z4], ls$landscape[type == "treatment", z4])


ks_stat(ls$landscape[type == "reference", z1], ls$landscape[type == "treatment", z1])
ks_stat(ls$landscape[type == "reference", z2], ls$landscape[type == "treatment", z2])
ks_stat(ls$landscape[type == "reference", z3], ls$landscape[type == "treatment", z3])
ks_stat(ls$landscape[type == "reference", z4], ls$landscape[type == "treatment", z4])



data <- ls$landscape
data[, type.o := factor(type, levels = levels(type), ordered = TRUE)]
# data[, zone := fifelse(x < 50, 1, 2)]

ls.sam <- sample(1:nrow(ls$landscape), 1e4)
data <- data[ls.sam]


# EGP NEW -----------------------------------------------------------

som.fit <-
  egp_som(data,
          x.dim = 25,
          y.dim = 25,
          epochs = 1000,
          vars = c("z1", "z2", "z3", "z4"),
          parallel = 4)


data[, c("som.unit", "som.x", "som.y") := get_bmu(som.fit, coord = TRUE, list = TRUE)]

mod.egp <- bam(response ~
               s(x, y, bs = "gp", k = 250,
                 xt = list(max.knots = 2500)) +
               s(x, y, by = type.o, bs = "gp", k = 250,
                 xt = list(max.knots = 2500)) +
               # s(x, y, by = type, bs = "gp", k = 250,
               #   xt = list(max.knots = 2500)) +
               s(som.x, som.y, bs = "gp", k = 250,
                 xt = list(max.knots = 625)),
               data = data,
               select = TRUE,
               discrete = TRUE,
               nthreads = 4
               )

summary(mod.egp)


# ls$landscape[, type := factor(type, levels = levels(type), ordered = FALSE)]
# mod.gam <- 
#   gam(response.int ~
#         s(x, y, bs = "sz", k = 100, xt = list(bs = "gp", m = 3)) +
#         s(type, x, y, bs = "sz", k = 100, xt = list(bs = "gp", m = 3)) +
#         s(som.x, som.y, bs = "gp", m = 3, k = 100),
#       method = "REML",
#       optimizer = "efs",
#       data = ls$landscape)

# summary(mod.gam)

mod.post <-
  egp_posterior_draw(mod.egp, 1000, unconditional = TRUE, package = "mgcv")


pred <-
  egp_posterior_predict(
                        model = mod.egp,
                        data = data,
                        id.var = "cell",
                        # pred.name = "response",
                        # ids = sam2,
                        posterior = mod.post)

# data[, cell2 := cell]

cf <-
  egp_define_counterfactual(
                            data = data,
                            som = som.fit,
                            fac.ids = data[type == "treatment", cell],
                            cf.ids = data[type == "reference", cell],
                            # compare.by = "zone",
                            # nb.strategy = "sequential",
                            nb.strategy = "expand",
                            group.by = list(NULL, "poly"),
                            som.var = "som.unit",
                            id.var = "cell",
                            deg.max = NULL,
                            n.min = 1)

units <- egp_summarize_units(pred, cf)
fac <- egp_evaluate_factual(pred, cf, agg.size = 1e2)
count <- egp_evaluate_counterfactual(pred, cf, units)
mar <- egp_marginal(fac, count)

mar[, .(mean(marginal), sd(marginal)), by =  .(group.id = group.id)]


cf2$assignment
cf$assignment

data[, mean(treatment), by = poly]

cf2 <-
  egp_define_counterfactual(
                            data = data,
                            som = som.fit,
                            fac.ids = data[type == "treatment", cell],
                            cf.ids = data[type == "reference", cell],
                            # compare.by = "zone",
                            nb.strategy = "expand",
                            group.by = list(NULL, "poly"),
                            som.var = "som.unit",
                            id.var = "cell",
                            deg.max = NULL,
                            n.min = 1)


units <- egp_summarize_units(pred, cf2)
fac <- egp_evaluate_factual(pred, cf2, agg.size = 1e2)
count <- egp_evaluate_counterfactual(pred, cf2, units)
mar <- egp_marginal(fac, count)

mar[, .(mean(ATE), sd(ATE)), by =  .(group.id = group.id)]


# EGP NEW -----------------------------------------------------------

som2 <-
  egp_som(data,
          x.dim = 10,
          y.dim = 10,
          epochs = 1000,
          vars = c("z1", "z2", "z3", "z4"),
          parallel = 4)


data[, c("som.unit", "som.x", "som.y") := get_bmu(som2, coord = TRUE, list = TRUE)]

mod.egp <- bam(response ~
               s(x, y, bs = "gp", k = 250,
                 xt = list(max.knots = 2500)) +
               s(x, y, by = type.o, bs = "gp", k = 250,
                 xt = list(max.knots = 2500)) +
               s(som.x, som.y, bs = "gp", k = 250,
                 xt = list(max.knots = 625)),
               data = data,
               select = TRUE,
               discrete = TRUE,
               nthreads = 4
               )

summary(mod.egp)


# ls$landscape[, type := factor(type, levels = levels(type), ordered = FALSE)]
# mod.gam <- 
#   gam(response.int ~
#         s(x, y, bs = "sz", k = 100, xt = list(bs = "gp", m = 3)) +
#         s(type, x, y, bs = "sz", k = 100, xt = list(bs = "gp", m = 3)) +
#         s(som.x, som.y, bs = "gp", m = 3, k = 100),
#       method = "REML",
#       optimizer = "efs",
#       data = ls$landscape)

# summary(mod.gam)

mod.post <-
  egp_posterior_draw(mod.egp, 1000, unconditional = TRUE, package = "mgcv")


pred <-
  egp_posterior_predict(
                        model = mod.egp,
                        data = data,
                        id.var = "cell",
                        # pred.name = "response",
                        # ids = sam2,
                        posterior = mod.post)

# data[, cell2 := cell]

cf <-
  egp_define_counterfactual(
                            data = data,
                            som = som.fit,
                            fac.ids = data[,cell],
                            cf.ids = data[type == "reference", cell],
                            # compare.by = "zone",
                            nb.strategy = "expand",
                            group.by = list("type", "poly"),
                            som.var = "som.unit",
                            id.var = "cell",
                            deg.max = NULL,
                            n.min = 1)

units <- egp_summarize_units(pred, cf)
fac <- egp_evaluate_factual(pred, cf, agg.size = 1e2)
count <- egp_evaluate_counterfactual(pred, cf, units)
mar <- egp_marginal(fac, count)

mar[, .(mean(marginal), sd(marginal)), by =  .(group.id = group.id)]

# EGP old -----------------------------------------------------------


som.dim <- 25
som.rlen <- 1000
egp.k.som <- 250
egp.k.geo <- 250
egp.max.knots.som <- som.dim^2
egp.max.knots.geo <- egp.k.geo*10
egp.approx <- TRUE
egp.basis <- "gp"
egp.select <- TRUE
n.threads = 4




  ## SOM #######################################################################

  message(paste0("Fitting SOM …"))

  ls.sam <- data

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
  
  # # results.mod[["sample"]] <- ls.sam
  # results.mod[["som"]] <- som.fit

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

    mod.egp <- bam(response ~
                   s(x, y, by = type.o, bs = egp.basis, k = egp.k.geo,
                     xt = list(max.knots = egp.max.knots.geo)) +
                   s(som_x, som_y, bs = egp.basis, k = egp.k.som,
                     xt = list(max.knots = egp.max.knots.som)),
                   data = ls.sam,
                   select = TRUE,
                   discrete = TRUE,
                   nthreads = n.threads
                   )
  # summary(mod.egp)
  # AIC(mod.egp)
  # summary(mod.egp)
  # AIC(mod.egp)

  
  # Estimate marginal effect

  message("Estimating marginal effect …")

  data.bl <- assign_bl_som(ls.sam,
                           som = som.fit,
                           bl.level = "reference",
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
      ls.sam[type == "reference"] |>
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

ls.sam[, type := factor(type, ordered = FALSE)]

# Match --------------------------------------------------------------


  ls.sam <- data

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











dev.new()

merge(data[, .(x, y, cell, treatment)]) |>
ggplot() +
  geom_raster(aes(x = x, y = y, fill = treatment)) +
  scale_fill_viridis_c()

mar

matched.nn.mh.re <- matchit(type ~ z1 + z2 + z3 + z4,
                            # data = ls$landscape[sample(1:1e4, 500)], method = "cem")
                            data = ls$landscape, method = "cem")
    md.nn.mh.re <- match.data(matched.nn.mh.re)
    mod.nn.mh.re <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.mh.re)


data[, mean(treatment), by = poly]

cf$groups

fac[, mean(factual), by = group.id]


cf$units[, n := unlist(lapply(cell, length))]

predictions = pred
cf.def <- cf
name = "counterfactual"


factual <- fac
counterfactual <- count


get_influence(cf)[.influence > 0]


data = data.dt,
id.ol


predictions <- as.data.table(predictions)
counterfactual <- cf




fac[, .(mean(factual), sd(factual)), by = poly]




# ARGS: predictions, data, group.vars,











                      model = mod.gam
                      data = data[sample(1:1e4, 10)]
                      id.var = "cell"
                      pred.name = "response"
                      ids = 2:5
                      posterior = mod.post


system.time({


sam <- sample(1:1e4, 10)
sam2 <- sam[1:2]
sam


pred

pred.m <- pred[, .(response.pred = mean(response)), by = cell]

cbind(data, pred.m[, -"cell"]) |>
ggplot() +
  geom_raster(aes(x = x, y = y, fill = response.int)) +
  scale_fill_continuous_divergingx("Roma", rev = TRUE)


# Prepare


set.seed(19010511+1) # Rose Ausländer
n <- 1000
ls.dim <- 100
cov.effect.mu <- c(-2, 2)
cov.effect.range <- c(0.2, 5)
z1.mix.w <- runif(n, 0.2, 0.5)
z2.mix.w <- runif(n, 0.2, 0.5)
z3.mix.prop <- runif(n, 0.2, 0.5)
z4.mix.w <- runif(n, 0.2, 0.5)
parameters <-
  data.table(
             id = 1:n,
             seed = sample(1:1e8, n),
             x.dim = ls.dim,
             y.dim = ls.dim,
             # Imbalance
             areas.imbalance = 0.1,
             # Effects
             treatment.eff.mean = 1,
             z1.effect.type = sample(c("sigmoid", "minimum", "unimodal", "bimodal"), n,
                                     replace = TRUE),
             z1.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]),
             z1.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             z2.effect.type = sample(c("sigmoid", "minimum", "unimodal", "bimodal"), n,
                                     replace = TRUE),
             z2.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]),
             z2.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             z3.effect.type = sample(c("sigmoid", "minimum", "unimodal", "bimodal"), n,
                                     replace = TRUE),
             z3.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]),
             z3.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             z4.effect.type = sample(c("sigmoid", "minimum", "unimodal", "bimodal"), n,
                                     replace = TRUE),
             z4.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]),
             z4.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             int12.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]) ,
             int12.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             int12.effect.nuclei = sample(20:50, n, replace = TRUE),
             int13.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]) ,
             int13.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             int13.effect.nuclei = sample(20:50, n, replace = TRUE),
             int14.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]) ,
             int14.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             int14.effect.nuclei = sample(20:50, n, replace = TRUE),
             int23.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]) ,
             int23.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             int23.effect.nuclei = sample(20:50, n, replace = TRUE),
             int24.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]) ,
             int24.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             int24.effect.nuclei = sample(20:50, n, replace = TRUE),
             int34.effect.range = runif(n, cov.effect.range[1], cov.effect.range[2]) ,
             int34.effect.mu = runif(n, cov.effect.mu[1], cov.effect.mu[2]),
             int34.effect.nuclei = sample(20:50, n, replace = TRUE),
             # Residual variation
             e.exp.var = 0.2,
             e.exp.scale = runif(n, 0.1*ls.dim, ls.dim),
             e.nug.var = 0.2,
             # Parameters for generating functions
             treatment.mat.nu = runif(n, 1, 1.5),
             treatment.mat.scale = runif(n, 0.05*ls.dim, 0.25*ls.dim),
             treatment.mat.var = runif(n, 0.5, 2),
             treatment.damp.infl = 0.2,
             treatment.damp.scale = runif(n, 1/(5*pi), 1/(2*pi)),
             mix.alpha = runif(n, 0.5, 1.5),
             z1.fbm.alpha = runif(n, 0.5, 1.5),
             z1.mix.w = z1.mix.w,
             z2.grad.phi = runif(n, 0, 2*pi),
             z2.mix.w = z2.mix.w,
             z3.dist.n = sample(10:25, n, replace = TRUE),
             z3.dist.acc = ls.dim/100,
             z3.mix.prop = z3.mix.prop,
             z4.seg.n = sample(10:25, n, replace = TRUE),
             z4.mat.nu = runif(n, 1, 2),
             z4.mat.scale = runif(n, 0.1, 1),
             z4.mix.w = z4.mix.w,
             # Parameters for defining treatment and reference areas
             score.type = "mahalanobis",
             areas.seg.seed = 3,
             areas.score.sam = 1e4,
             areas.area.prop = runif(n, 0.35, 0.65),
             areas.area.exact = FALSE,
             areas.seg.res = 0.05 * ls.dim,
             areas.seg.min.dist = 0.03*ls.dim,
             areas.seg.min.area = (ls.dim/10)^2,
             areas.seg.even = 2,
             areas.seg.prec = 5e-4 * ls.dim,
             areas.min.bound.dist = 0,
             areas.opt.imp.imb = 5,
             areas.opt.imp.area = 1,
             areas.opt.pop = 100,
             areas.opt.prec = 1e-4,
             areas.opt.pcrossover = 0.9,
             areas.opt.pmutation = 0.5,
             areas.opt.max.iter = 500,
             areas.opt.run = 100,
             areas.opt.parallel = 4,
             areas.opt.fine = TRUE,
             areas.opt.fine.max.iter = 100,
             areas.opt.fine.constr = TRUE,
             areas.opt.fine.tol = 1e-6,
             verbose = 2
             )
rm(z1.mix.w, z2.mix.w, z3.mix.prop, z4.mix.w)

# attach(as.list(parameters[id == 1]))
# detach(as.list(parameters[id == 1]))

ls.par <-
parameters[z1.effect.type != z2.effect.type &
           z1.effect.type != z3.effect.type &
           z1.effect.type != z4.effect.type &
           z2.effect.type != z3.effect.type &
           z2.effect.type != z4.effect.type &
           z3.effect.type != z4.effect.type][1] |>
         as.list()
# ls.par <- as.list(parameters[id == 1])
ls.par$areas.imbalance <- 0.1
ls <- do.call(generate_landscape_4cov_nl, ls.par)



plot_landscape_4cov_nl(ls$landscape, select = "all", interactions = TRUE)
plot_functions_4cov_nl(ls$fun)

som.fit <-
  egp_som(ls$landscape,
          x.dim = 10,
          y.dim = 10,
          vars = c("z1", "z2", "z3", "z4"),
          parallel = 4)

get_bmu(som.fit)
get_codes(som.fit)
get_nb(som.fit)
get_poly(som.fit)

som.fit$scale$mean

x <- ls$landscape
som <- som.fit




data <- ls$landscape
data[, c("som.bmu", "som.x", "som.y") := get_bmu(som.fit, coord = TRUE, list = TRUE)]


data[order(som.bmu),
     .(som.mean = mean(response.int), n = .N),
     by = .(type, som.bmu)] |>
dcast(som.bmu ~ type, value.var = "som.mean")  |>
DT(,.(som.bmu, diff = treatment - reference)) |>
merge(get_poly(som.fit), by.x = "som.bmu", by.y = "id") |>
st_as_sf() |>
ggplot() +
  geom_sf(aes(fill = diff)) +
  # geom_sf(aes(fill = som.mean)) +
  # geom_sf(aes(fill = n)) +
  scale_fill_continuous_divergingx("Roma", rev = TRUE)
  # scale_fill_viridis_c()



get_bmu(som.fit, coord = FALSE, list = TRUE)

data
cf.ids = data[type == "reference", cell]
id.var = "cell"
# by = "z4"
# by.vars = NULL
by.vars = "subgroup"
som.var = "som.bmu"
som = som.fit
n.min = 1
deg.max = NULL
cf.name = "cf.bmu"
list = FALSE


cf <-
egp_define_counterfactual(data = data,
                          cf.ids = data[type == "reference", cell],
                          # by.vars = "subgroup",
                          som.var = "som.bmu",
                          id.var = "cell",
                          deg.max = NULL,
                          n.min = 1)

get_influence(cf)


mod.post <-
egp_posterior_draw(mod.gam, 1000, unconditional = TRUE, package = "mgcv")

system.time({
pred <-
egp_posterior_predict(mod.gam,
                      data = data,
                      id.var = "cell",
                      posterior = mod.post)
})


som.sum <-
lapply(cf$ids$cell, \(x) apply(pred[,paste0("cell:", x[[1]])], 1, mean))

x <- cf$ids$cell[1]

som.pred <- list()

as.data.table(as_draws_df(pred))

data.table(draw = 1:ndraws,
           post.pred = apply(pred[,paste0("cell:", x[[1]])], 1, mean))


lapply(som.sum, length)

pred[, "cell:8449"]

f.mod[8449]

cf2 <- cf$ids

cf2[, post.pred := som.sum]
ndraws <- nrow(pred)

cf

cf.long <-
cf2[, .(cf.bmu, post.pred = unlist(post.pred)), by = .I][, .(m.pred = mean(post.pred)), by =cf.bmu]

# unlist and include draw column



data.dt2 <- data
data.dt2$fitted <- f.mod

cf$ids


merge(
data.dt2[type == "reference", .(m = mean(fitted)), som.bmu],
cf.long[, .(som.bmu = cf.bmu, m.pred)],
by = "som.bmu")


cf.long[,
        .(pred.mean = mean(post.pred)),
        by = .(som.bmu = cf.bmu)])

data.dt3
|>



str(som.sum)

cf

pred.m <- apply(pred, 2, mean)
f.mod <- fitted(mod.gam)

summary(f.mod - pred.m)

ggplot() +
  geom_point(aes(x = f.mod, y = pred.m))


model = mod.gam
# data = NULL
data = ls$landscape
posterior = mod.post
id.var = "cell"
type = "response"
ids = c(1, 10)
coef = NULL
predict.chunk = NULL
post.chunk = NULL
progress = TRUE
model = model
posterior = posterior
newdata = data
id.col = id.var
type = type
obs = ids
predict.chunk = predict.chunk
post.chunk = post.chunk
progress = progress
discrete = TRUE
                             


# OKAY!
egp_posterior_draw
# ARGS: model, n, unconditional

# OKAY!
egp_posterior_predict
# ARGS: model, data, coef = NULL

egp_summarize_som
# ARGS: predictions, counterfactual.def, pred.name = "predicted"

egp_evaluate_factual
# ARGS: predictions, data, group.vars,

egp_evaluate_counterfactual
# ARGS: predictions, counterfactual.def, data, group.vars, som.summary

egp_evaluate_marginal
# ARGS: factual, counterfactual, type = "absolute" "relative"





x <- cf



x$ids[, .(cell = unlist(cell)), by = .(subgroup, cf.bmu)][cell == 730]

x$ids[subgroup == 1 & cf.bmu == 87]

counts





    ids.cf



    data.dt
    
    data.cf[som.bmu == 96]

    groups.cf[, c(".grkey", by.vars, som.var, cf.name), with = FALSE] |>
    merge(data.dt, by = by.vars)

    data.grp <- merge(data.dt, by = by)



     nbh.l[[data.grp$.grkey]]

.grkey <- 1

    cf <- 

      data.grp[,
              .(id.col,
                # cf.col = nbh.l[[.grkey]][som.col]),
                cf.col = .grkey),
              env = list(id.col = id.var,
                         cf.col = cf.name,
                         som.col = som.var)]

  }

  if(list == TRUE) {
    return(as.list(cf))
  } else {
    return(cf)
  }
}




  n <- n.bmu
  n.min <- 1
  deg.max <- 100

  data[, .(som.bmu, nbh$neighbourhood[som.bmu])][som.bmu == 97]

  nbh <- find_neighbourhood(A, n)


find_neighbourhood <- function(A,
                               n,
                               n.min = 1,
                               deg.max = NULL) {

  if(is.null(deg.max)) {
    deg.max <- length(n)
  }

  nbh <- list() 

  n.nbh <- n
  deg.nbh <- integer(length(n))

  # deg 0 nbh

  bmu.sel <- which(n.nbh >= n.min)
  nbh[bmu.sel] <- bmu.sel

  deg <- 0
  while(any(n.nbh < n.min)) {
    deg <- deg + 1
    # deg 1 nbh
    bmu.sel <- which(n.nbh < n.min)
    if(deg < 2) {
      A.nbh <- A
    } else {
      A.nbh <- A.nbh %*% A
    }
    nbh.sel <- apply(A.nbh[bmu.sel,,drop = FALSE], 1,
                    \(x) which(x > 0),
                    simplify = FALSE)
    n.nbh[bmu.sel] <- unlist(lapply(nbh.sel, \(x) sum(n[x])))
    deg.nbh[bmu.sel] <- deg
    nbh[bmu.sel] <- nbh.sel
    if(deg == deg.max) {
      nbh[which(n.nbh < n.min)] <- NA
      break
    }
  }
  
  neighbourhood <-
    list(neighbourhood = nbh,
         n = n.nbh,
         degree = deg.nbh)

  return(neighbourhood)
}


}



# if(som.var length of rows: use as bmus; if length 1, use as var
# selector

# if id.col not present set to NULL; if id.col NULL, create "id" and
# rename later to row.id

# Function to get neighbours for degree n, then number of ref units in
                              # these
# Find units: while condition, compute size in neighbours, if >= n.min,
                              # store, if not, leave NA, continoure only
                              # with the NA units; abort when max no of
# neighbours reached



x <- ls$landscape



grid$pts

  hex_poly(grid$pts)

  pts.l <-
    apply(grid$pts, 1, \(x) st_point(x, dim = "XY"),
          simplify = FALSE)


grid_rectangular_poly(x.dim = 10, y.dim = 10, d = 1) |>
st_centroid()

|>
plot()

  st_centroid(grid_hex_poly(x.dim = 10, y.dim = 10, d = 1))[16:20,]



st_sf(pts)[16:25,]
  
  pts <-
    pts.l |>
    st_sfc()

    |>
    st_union() |>
    st_voronoi() |>
    st_cast() |>
    # st_intersection(poly.lim) |>
    st_as_sf() |>
    st_set_geometry("geometry") |>
    plot()

  pts$id <- 1:nrow(pts)

plot(pts)

st_sample(envelope, 10)


grid.x.min <- min(grid$pts[,1])
grid.x.max <- max(grid$pts[,1])
grid.y.min <- min(grid$pts[,2])
grid.y.max <- max(grid$pts[,2])

grid.x.min <- 1
grid.x.max <- 10
grid.y.min <- 0
grid.y.max <- 10

poly.lim <-
  c(grid.x.min, grid.y.min,
    grid.x.max, grid.y.min,
    grid.x.max, grid.y.max,
    grid.x.min, grid.y.max,
    grid.x.min, grid.y.min) |>
  matrix(ncol = 2, byrow = TRUE) |>
  list() |>
  st_polygon()


  poly.lim <-
    matrix(c(0, 0,
             min(grid$pts[,1], 0,
             min(grid$pts[,1], max(grid.,
             0, grid$ydim,
             0, 0),
           ncol = 2, byrow = 2) |>
    list() |>
    st_polygon()

    plot(poly.lim)
    plot(pts, add = TRUE)

  st_make_grid(poly.lim, n = c(10, 10)) |>
  st_centroid()

  |>
  st_intersects(poly.lim, sparse = FALSE) |>
  sum()

  |>


  plot()

  st_centroid()


  plot()


  st_voronoi(pts, envelope) |>
  st_intersection(envelope) |>
  plot()

    |>
    st_voronoi() |>
    plot()

    st_sf() |>

    do.call(c, pts.l)

  st_point($pts[1,])



  init.pca

  init = som_init_pc
  init <- init(grid, x.scaled)
  apply(init, 2, range)
  apply(x.scaled, 2, range)
