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
library(MatchIt)

RFoptions(install="no")

source("utilities_fcne.R")
source("utilities.R")


# Prepare


# new treatment algorithm

n <- 1000
ls.dim <- 100


set.seed(19010511+1) # Rose Ausländer
parameters <-
  data.table(
             id = 1:n,
             seed = sample(1:1e8, n),
             x.dim = ls.dim,
             y.dim = ls.dim,
             # Effect sizes
             treatment.effect.size.sp = 1,
             treatment.effect.size.bd = 0,
             z1.effect.type = sample(c("sigmoid", "minimum", "unimodal", "bimodal"), n,
                                     replace = TRUE),
             z1.effect.range = runif(n, 0.5, 4),
             z1.effect.mu = runif(n, -0.5, 0.5),
             z2.effect.type = sample(c("sigmoid", "minimum", "unimodal", "bimodal"), n,
                                     replace = TRUE),
             z2.effect.range = runif(n, 0.5, 4),
             z2.effect.mu = runif(n, -0.5, 0.5),
             z3.effect.type = sample(c("sigmoid", "minimum", "unimodal", "bimodal"), n,
                                     replace = TRUE),
             z3.effect.range = runif(n, 0.5, 4),
             z3.effect.mu = runif(n, -0.5, 0.5),
             z4.effect.type = sample(c("sigmoid", "minimum", "unimodal", "bimodal"), n,
                                     replace = TRUE),
             z4.effect.range = runif(n, 0.5, 4),
             z4.effect.mu = runif(n, -0.5, 0.5),
             int12.effect.range = runif(n, 0.25, 2) ,
             int12.effect.mu = runif(n, -0.25, 0.25),
             int12.effect.nuclei = sample(2:5, n, replace = TRUE),
             int13.effect.range = runif(n, 0.25, 2) ,
             int13.effect.mu = runif(n, -0.25, 0.25),
             int13.effect.nuclei = sample(2:5, n, replace = TRUE),
             int14.effect.range = runif(n, 0.25, 2) ,
             int14.effect.mu = runif(n, -0.25, 0.25),
             int14.effect.nuclei = sample(2:5, n, replace = TRUE),
             int23.effect.range = runif(n, 0.25, 2) ,
             int23.effect.mu = runif(n, -0.25, 0.25),
             int23.effect.nuclei = sample(2:5, n, replace = TRUE),
             int24.effect.range = runif(n, 0.25, 2) ,
             int24.effect.mu = runif(n, -0.25, 0.25),
             int24.effect.nuclei = sample(2:5, n, replace = TRUE),
             int34.effect.range = runif(n, 0.25, 2) ,
             int34.effect.mu = runif(n, -0.25, 0.25),
             int34.effect.nuclei = sample(2:5, n, replace = TRUE),
             # Parameters for generating functions
             # treatment.nuclei = sample(3:7, n, replace = TRUE),
             # treatment.phi.range = list(c(-pi/4, pi/4)),
             # treatment.shift.range = list(c(-0.25, 0.25)),
             # treatment.x.scale.range = list(c(0.5, 1)),
             # treatment.y.scale.range = list(c(0.25, 0.5)),
             # treatment.nuc.eff.range = list(c(0.5,2)),
             treatment.nuclei = 10,
             treatment.phi.range = list(c(-pi/4, pi/4)),
             treatment.shift.range = list(c(0, 0.5)),
             treatment.x.scale.range = list(c(0.1, 0.5)),
             treatment.y.scale.range = list(c(0.1, 0.5)),
             treatment.nuc.eff.range = list(c(0.5,2)),
             z1.fbm.alpha = runif(n, 0.5, 1.5),
             z1.fbm.var = runif(n, 0.1, 1),
             z1.fbm.scale = runif(n, 0.1, 1),
             # z1.fbm.w = 0.8,
             z1.fbm.w = 1,
             # z1.fbm2.alpha = runif(n, 0.5, 1.5),
             # z1.fbm2.var = runif(n, 0.1, 1),
             # z1.fbm2.scale = runif(n, 0.1, 1),
             # z1.fbm.ratio = runif(n, 1, 2),
             # z1.grad.phi = runif(n, -1, 1) * 2 * pi,
             z1.grad.phi = sample(c(-pi/8, pi/8), n, replace = TRUE) +
                           sample(c(0, pi), n, replace = TRUE),
             z1.grad.w = 0.2,
             z2.fbm.alpha = runif(n, 0.5, 1.5),
             z2.fbm.var = runif(n, 0.1, 1),
             z2.fbm.scale = runif(n, 0.1, 1),
             z2.fbm.w = 0.2,
             # z2.grad.phi = runif(n, -1, 1) * 2 * pi,
             # z2.grad.phi = runif(n, -pi/4, pi/4) + sample(c(0, pi), n),
             z2.grad.phi = sample(c(-pi/8, pi/8), n, replace = TRUE) +
                           sample(c(0, pi), n, replace = TRUE),
             z2.grad.shift = 0,
             z2.grad.w = 0.8,
             z3.dist.n = sample(10:25, n, replace = TRUE),
             # z3.grad.phi = runif(n, -1, 1) * 2 * pi,
             z3.grad.phi = runif(n, -pi/8, pi/8) +
                           sample(c(0, pi), n, replace = TRUE),
             z3.grad.prop = 1,
             z3.acc = 1,
             z4.seg.n = sample(10:25, n, replace = TRUE),
             z4.mat.nu = runif(n, 1, 2),
             z4.mat.var = runif(n, 0.1, 1),
             z4.mat.scale = runif(n, 0.1, 1),
             z4.mat.w = 1,
             # z4.grad.phi = runif(n, -1, 1) * 2 * pi,
             z4.grad.phi = runif(n, -pi/8, pi/8) +
                           sample(c(0, pi), n, replace = TRUE),
             z4.grad.w = 0,
             split.n = 200,
             split.prop = runif(n, 0.4, 0.6),
             e.exp.var = 0.1,
             e.exp.scale = 100,
             e.nug.var = 0.1,
             e.rand.var = 0.1
             )

attach(as.list(parameters[id == 1]))

  z1 <- generate_z1(x.dim = x.dim,
                    y.dim = y.dim,
                    fbm.alpha = z1.fbm.alpha,
                    fbm.var = z1.fbm.var,
                    fbm.scale = z1.fbm.scale,
                    fbm.w = z1.fbm.w,
                    grad.phi = z1.grad.phi,
                    grad.w = z1.grad.w,
                    rescale = c(0,1),
                    name = "z1")

  z2 <- generate_z2(
                    x.dim = x.dim,
                    y.dim = y.dim,
                    fbm.alpha = z2.fbm.alpha,
                    fbm.var = z2.fbm.var,
                    fbm.scale = z2.fbm.scale,
                    fbm.w = z2.fbm.w,
                    grad.phi = z2.grad.phi,
                    grad.shift = z2.grad.shift,
                    grad.w = z2.grad.w,
                    rescale = c(0,1),
                    name = "z2")

  z3 <- generate_z3(x.dim = x.dim,
                    y.dim = y.dim,
                    dist.n = z3.dist.n,
                    grad.phi = z3.grad.phi,
                    grad.prop = z3.grad.prop,
                    acc = z3.acc,
                    rescale = c(0,1),
                    name = "z3")

  z4 <- generate_z4(x.dim = x.dim,
                    y.dim = y.dim,
                    seg.n = z4.seg.n,
                    mat.nu = z4.mat.nu,
                    mat.var = z4.mat.var,
                    mat.scale = z4.mat.scale,
                    mat.w = z4.mat.w,
                    grad.phi = z4.grad.phi,
                    grad.w = z4.grad.w,
                    rescale = c(0,1),
                    name = "z4")

covariates <- c(z1, z2, z3, z4)

cov.inv <-
  inverse(covariates,
          attributes = sample(c(TRUE, FALSE), size = length(covariates), replace = TRUE)) |>
  scale_int()

scores <- list()
scores$prod <-
  score(cov.inv, type = "prod") |>
  scale_int()
scores$euclidean <-
  score(cov.inv, type = "euclidean") |>
  scale_int()
scores$sumofsquares <-
  score(cov.inv, type = "sumofsquares") |>
  scale_int()
scores$mahalanobis <-
  score(cov.inv, type = "mahalanobis") |>
  scale_int()



score.cov <- list()
for(i in seq_along(covariates)) {
  sign <- sample(c(-1, 1), 1)
  score.cov[[i]] <-
    (covariates[i] * sign) |>
    scale_int(c(0, 1))
}

prob.cov
  do.call(c, score.cov) |>
  merge() |>
  st_apply(1:2, prod, .fname = "score") |>
  scale_int(c(0, 1))

score.cov <- 
  do.call(c, score.cov) |>
  merge() |>
  st_apply(1:2, \(x) sum(x^2), .fname = "score") |>
  scale_int(c(0, 1))

zplot(prob.cov)

# END preparation 

load("poly.RData")


x.dim = 100
y.dim = 100
imbalance = 0.3
# imbalance.tol = NULL
imbalance.tol = 0.01
area.prop = 0.5
# area.tol = NULL
area.tol = c(0.4, 0.6)
area.exact = FALSE
score = scores$mahalanobis
score.sam = 1e4
seg.seed = 3
seg.res = 5
seg.min.dist = 5
seg.min.area = 100
seg.even = 2
seg.prec = 0.1
min.bound.dist = 0.1
verbose = TRUE
opt.imp.imb = 5
opt.imp.area = 1
opt.pop = 100
opt.prec = 1e-4
opt.pcrossover = 0.9
opt.pmutation = 0.3
opt.max.iter = 100
opt.run = 25
opt.parallel = 4
opt.fine = TRUE
opt.fine.max.iter = 100
opt.fine.constr = TRUE
opt.fine.tol = 1e-6


trt.ar <-
generate_areas_poly2(
x.dim = 100,
y.dim = 100,
imbalance = 0.1,
# imbalance.tol = NULL,
imbalance.tol = 0.01,
area.prop = 0.5,
area.tol = NULL,
# area.tol = c(0.4, 0.6),
area.exact = FALSE,
score = scores$mahalanobis,
score.sam = 1e4,
seg.seed = 3,
seg.res = 5,
seg.min.dist = 5,
seg.min.area = 100,
seg.even = 1,
seg.prec = 0.1,
min.bound.dist = 0,
verbose = TRUE,
opt.imp.imb = 5,
opt.imp.area = 1,
opt.pop = 100,
opt.prec = 1e-4,
opt.pcrossover = 0.9,
opt.pmutation = 0.3,
opt.max.iter = 250,
opt.run = 25,
opt.parallel = 4,
opt.fine = TRUE,
opt.fine.max.iter = 100,
opt.fine.constr = TRUE,
opt.fine.tol = 1e-6)

plot(trt.ar$shape)

save.image("poly.RData")


stopImplicitCluster()
registerDoParallel(4)

n <- 50

balance.test <-
  expand.grid(
              imbalance = c(0.1, 0.2, 0.3),
              # distance = c("prod", "sumofsquares", "euclidean", "mahalanobis"))
              distance = c("mahalanobis"))
setDT(balance.test)
balance.test[, par.id := rleidv(balance.test)]

cov.inv <-
  inverse(covariates,
          attributes = sample(c(TRUE, FALSE), size = length(covariates), replace = TRUE)) |>
  scale_int()

scores <- list()
scores$prod <-
  score(cov.inv, type = "prod") |>
  scale_int()
scores$euclidean <-
  score(cov.inv, type = "euclidean") |>
  scale_int()
scores$sumofsquares <-
  score(cov.inv, type = "sumofsquares") |>
  scale_int()
scores$mahalanobis <-
  score(cov.inv, type = "mahalanobis") |>
  scale_int()
scores$qprod <-
  score(cov.inv, type = "qprod") |>
  scale_int()
scores$qss <-
  score(cov.inv, type = "qss") |>
  scale_int()

zplot(scores$qprod)

imb.test <- c(0.1, 0.2, 0.3)

maha.l <- list()
for(i in seq_along(imb.test)) {
ar.trt <-
  generate_areas_poly2(
    x.dim = 100,
    y.dim = 100,
    imbalance = imb.test[i],
    # imbalance.tol = NULL,
    imbalance.tol = 0.01,
    area.prop = 0.5,
    # area.tol = NULL,
    area.tol = c(0.4, 0.6),
    area.exact = FALSE,
    score = scores$mahalanobis,
    score.sam = 1e4,
    seg.seed = 3,
    seg.res = 5,
    seg.min.dist = 5,
    seg.min.area = 100,
    seg.even = 1,
    seg.prec = 0.1,
    min.bound.dist = 0,
    verbose = TRUE,
    opt.imp.imb = 5,
    opt.imp.area = 1,
    opt.pop = 100,
    opt.prec = 1e-4,
    opt.pcrossover = 0.9,
    opt.pmutation = 0.3,
    opt.max.iter = 250,
    opt.run = 25,
    opt.parallel = 4,
    opt.fine = TRUE,
    opt.fine.max.iter = 100,
    opt.fine.constr = TRUE,
    opt.fine.tol = 1e-6)
# plot(ar.trt$shape)
ar.shape <- ar.trt$shape
res.l <- list()
for(j in seq_along(covariates)) {
  ref.obs <- na.omit(as.vector(covariates[j][ar.shape[ar.shape$type == "reference",]][[1]]))
  trt.obs <- na.omit(as.vector(covariates[j][ar.shape[ar.shape$type == "treatment",]][[1]]))
  res.l[[j]] <-
    data.table(variable = names(covariates)[j],
               smd = d_cohen(trt.obs, ref.obs),
               vr = var_ratio(trt.obs, ref.obs),
               ks = ks_stat(trt.obs, ref.obs),
               js = js_div(trt.obs, ref.obs, type = "discrete")
    )
}
maha.l[[i]] <- rbindlist(res.l)
}

qss.l <- list()
for(i in seq_along(imb.test)) {
ar.trt <-
  generate_areas_poly2(
    x.dim = 100,
    y.dim = 100,
    imbalance = imb.test[i],
    # imbalance.tol = NULL,
    imbalance.tol = 0.01,
    area.prop = 0.5,
    # area.tol = NULL,
    area.tol = c(0.4, 0.6),
    area.exact = FALSE,
    score = scores$qss,
    score.sam = 1e4,
    seg.seed = 3,
    seg.res = 5,
    seg.min.dist = 5,
    seg.min.area = 100,
    seg.even = 1,
    seg.prec = 0.1,
    min.bound.dist = 0,
    verbose = TRUE,
    opt.imp.imb = 5,
    opt.imp.area = 1,
    opt.pop = 100,
    opt.prec = 1e-4,
    opt.pcrossover = 0.9,
    opt.pmutation = 0.3,
    opt.max.iter = 250,
    opt.run = 25,
    opt.parallel = 4,
    opt.fine = TRUE,
    opt.fine.max.iter = 100,
    opt.fine.constr = TRUE,
    opt.fine.tol = 1e-6)
# plot(ar.trt$shape)
ar.shape <- ar.trt$shape
res.l <- list()
for(j in seq_along(covariates)) {
  ref.obs <- na.omit(as.vector(covariates[j][ar.shape[ar.shape$type == "reference",]][[1]]))
  trt.obs <- na.omit(as.vector(covariates[j][ar.shape[ar.shape$type == "treatment",]][[1]]))
  res.l[[j]] <-
    data.table(variable = names(covariates)[j],
               smd = d_cohen(trt.obs, ref.obs),
               vr = var_ratio(trt.obs, ref.obs),
               ks = ks_stat(trt.obs, ref.obs),
               js = js_div(trt.obs, ref.obs, type = "discrete")
    )
}
qss.l[[i]] <- rbindlist(res.l)
}

qprod.l <- list()
for(i in seq_along(imb.test)) {
ar.trt <-
  generate_areas_poly2(
    x.dim = 100,
    y.dim = 100,
    imbalance = imb.test[i],
    # imbalance.tol = NULL,
    imbalance.tol = 0.01,
    area.prop = 0.5,
    # area.tol = NULL,
    area.tol = c(0.4, 0.6),
    area.exact = FALSE,
    score = scores$qprod,
    score.sam = 1e4,
    seg.seed = 3,
    seg.res = 5,
    seg.min.dist = 5,
    seg.min.area = 100,
    seg.even = 1,
    seg.prec = 0.1,
    min.bound.dist = 0,
    verbose = TRUE,
    opt.imp.imb = 5,
    opt.imp.area = 1,
    opt.pop = 100,
    opt.prec = 1e-4,
    opt.pcrossover = 0.9,
    opt.pmutation = 0.3,
    opt.max.iter = 250,
    opt.run = 25,
    opt.parallel = 4,
    opt.fine = TRUE,
    opt.fine.max.iter = 100,
    opt.fine.constr = TRUE,
    opt.fine.tol = 1e-6)
# plot(ar.trt$shape)
ar.shape <- ar.trt$shape
res.l <- list()
for(j in seq_along(covariates)) {
  ref.obs <- na.omit(as.vector(covariates[j][ar.shape[ar.shape$type == "reference",]][[1]]))
  trt.obs <- na.omit(as.vector(covariates[j][ar.shape[ar.shape$type == "treatment",]][[1]]))
  res.l[[j]] <-
    data.table(variable = names(covariates)[j],
               smd = d_cohen(trt.obs, ref.obs),
               vr = var_ratio(trt.obs, ref.obs),
               ks = ks_stat(trt.obs, ref.obs),
               js = js_div(trt.obs, ref.obs, type = "discrete")
    )
}
qprod.l[[i]] <- rbindlist(res.l)
}


ar.trt <-
  generate_areas_poly2(
    x.dim = 100,
    y.dim = 100,
    imbalance = 0.5,
    # imbalance.tol = NULL,
    imbalance.tol = 0.01,
    area.prop = 0.5,
    # area.tol = NULL,
    area.tol = c(0.4, 0.6),
    area.exact = FALSE,
    score = scores$qprod,
    score.sam = 1e4,
    seg.seed = 3,
    seg.res = 5,
    seg.min.dist = 5,
    seg.min.area = 100,
    seg.even = 1,
    seg.prec = 0.1,
    min.bound.dist = 0,
    verbose = TRUE,
    opt.imp.imb = 5,
    opt.imp.area = 1,
    opt.pop = 100,
    opt.prec = 1e-4,
    opt.pcrossover = 0.9,
    opt.pmutation = 0.3,
    opt.max.iter = 250,
    opt.run = 25,
    opt.parallel = 4,
    opt.fine = TRUE,
    opt.fine.max.iter = 100,
    opt.fine.constr = TRUE,
    opt.fine.tol = 1e-6)

plot(ar.trt$shape)

ar.shape <- ar.trt$shape
res.l <- list()
for(j in seq_along(covariates)) {
  ref.obs <- na.omit(as.vector(covariates[j][ar.shape[ar.shape$type == "reference",]][[1]]))
  trt.obs <- na.omit(as.vector(covariates[j][ar.shape[ar.shape$type == "treatment",]][[1]]))
  res.l[[j]] <-
    data.table(variable = names(covariates)[j],
               smd = d_cohen(trt.obs, ref.obs),
               vr = var_ratio(trt.obs, ref.obs),
               ks = ks_stat(trt.obs, ref.obs),
               js = js_div(trt.obs, ref.obs, type = "discrete")
    )
}
rbindlist(res.l)

