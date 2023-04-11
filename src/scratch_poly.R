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

ls.par <- as.list(parameters[id == 1])
ls.par$areas.imbalance <- 0.1
ls <- do.call(generate_landscape_4cov_nl, ls.par)


plot_landscape_4cov_nl(ls$landscape, select = "all", interactions = TRUE)
plot_functions_4cov_nl(ls$fun)

ls$landscape[, .(mean(response.int)), type]

range(ls$landscape$f.int23)

ls_theme <-
  theme_minimal(base_family = "IBMPlexSansCondensed", base_size = 11) +
  theme(
        plot.title = element_text(hjust = 0,
                                  face = "bold",
                                  margin = margin(l = 0, b = 5, t = 11)),
        plot.margin = margin(3, 3, 3, 3),
        strip.text = element_text(size = rel(0.8),
                                  hjust = 0.5),
        strip.background = element_rect(fill = "grey90", color = NA),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        legend.justification = "top"
  )



  guide_fill <-
    guides(fill = guide_colorbar(
                                 ticks.colour = "grey5",
                                 ticks.linewidth = 0.2,
                                 frame.colour = "grey5",
                                 frame.linewidth = 0.2,
                                 barwidth = 1,
                                 barheight = 5,
                                 label.position = "right",
                                 label.hjust = 1,
                                 draw.ulim = TRUE,
                                 draw.llim = TRUE
                                 ))




plot_landscape_4cov_nl(ls$landscape, select = "all", interactions = TRUE)

plot_functions_4cov_nl(ls$fun)

ls$landscape[, type := factor(type, levels = levels(type), ordered = TRUE)]

ls$landscape[, .(mean(response.int)), type]
ls$landscape[, .(mean(response.int)), poly]

mod.gam <- 
  bam(response.int ~
        s(x, y, bs = "gp", m = 3, k = 100) +
        s(x, y, by = type, bs = "gp", m = 3, k = 100),
      discrete = TRUE,
      data = ls$landscape)

summary(mod.gam)
AIC(mod.gam)

mod.pred <-
  predict(mod.gam, type = "terms", terms = "s(x,y):typetreatment") |>
  as.data.table() |>
  setnames(c("pred"))


field1 = st_as_stars(ls$landscape[, .(x, y, z1)])
                            field2 = st_as_stars(ls$landscape[, .(x, y, z2)])
                            range = 2
                            mu = 0.25
                            nuclei = sample(20:50, 1)
                            f.acc = 0.01


generate_nonlinear_effect(field = field1,
                          c("sigmoid", "minimum", "unimodal", "bimodal")[4],
                          range = 1,
                          mu = 0)$fun |>
ggplot() +
  geom_line(aes(x = z1, y = f.z1))

int.eff <-
generate_interaction_effect(field1 = st_as_stars(ls$landscape[, .(x, y, z1)]),
                            field2 = st_as_stars(ls$landscape[, .(x, y, z2)]),
                            range = 2,
                            mu = 0.25,
                            nuclei = sample(20:50, 1),
                            f.acc = 0.01)


int.eff$fun[exists == TRUE, .(val1 = z1, val2 = z2, f = f.z1z2)]  |>
      ggplot() +
        aes(x = val1, y = val2) +
        geom_raster(aes(fill = f)) +
        geom_contour(aes(z = f),
                     colour = "gray35",
                     linewidth = 0.2) +
        scale_fill_continuous_diverging(palette = "Tropic") +
        scale_x_continuous(n.breaks = 3) +
        scale_y_continuous(n.breaks = 3) +
        fun_theme

mod.pred <-
  cbind(ls$landscape, mod.pred)

mod.pred[, .(mean(pred)), by = type]

ls$landscape[, type := factor(type, levels = levels(type), ordered = FALSE)]
ls$landscape[, type := factor(type, levels = levels(type), ordered = FALSE)]

matched.nn.mh.re <- matchit(type ~ z1 + z2 + z3 + z4,
                            # data = ls$landscape[sample(1:1e4, 500)], method = "cem")
                            data = ls$landscape, method = "cem")
    md.nn.mh.re <- match.data(matched.nn.mh.re)
    mod.nn.mh.re <- lm(response.int ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.mh.re)

lm(response ~ type + z1 + z2 + z3 + z4, data = ls$landscape)


ls$landscape[, .(type = as.numeric(type == "treatment"), z1, z2, z3, z4)] |>
cor(method = "kendall")


ggplot(mod.pred) +
  geom_raster(aes(x = x, y=y, fill = treatment)) +
  scale_fill_viridis_c()

plot(mod.gam)

  theme(legend.justification = "top")

    layout.overview <-
      "ACC
       BCC"
    p.select <-
      plots$type +
      plots$response +
      plots$cov +
      plot_layout(design = layout.overview, widths = 1, heights = 1, guides = "collect") &
      ls_theme

    layout.eff <-
      "ABB
       DBB
       CCC
       CCC"
    p.select <-
      plots$treatment +
      plots$eff.main +
      plots$eff.int +
      plots$error +
      plot_layout(design = layout.eff, widths = 1, guides = "collect") &
      ls_theme

    layout.eff.noint <-
      "ABB
       CBB"
    p.select <-
      plots$treatment +
      plots$eff.main +
      plots$error +
      plot_layout(design = layout.eff.main, widths = 1, heights = 1, guides = "collect") &
     ls_theme

    layout.all <-
      "ACC
       BCC
       DEE
       GEE
       FFF
       FFF"
    p.select <-
      plots$type +
      plots$response +
      plots$cov +
      plots$treatment +
      plots$eff.main +
      plots$eff.int +
      plots$error +
      plot_layout(design = layout.all, widths = 1, guides = "collect") &
      ls_theme

    layout.all.noint <-
      "ACC
       BCC
       DEE
       FEE"
    p.select <-
      plots$type +
      plots$response +
      plots$cov +
      plots$treatment +
      plots$eff.main +
      plots$error +
      plot_layout(design = layout.all.noint, widths = 1, guides = "collect") &
      ls_theme


    plots$eff.int

    p.select <- 


      plots$type +
      plots$z1 + plots$z2 + plots$z3 + plots$z4 +
      plots$response +
      plot_layout(design = layout.overview, guides = "collect") &
      ls_theme

ls$landscape[, .(mean(response)), by = poly]
ls$landscape[, .(mean(treatment)), by = poly]



x<-ls$fun
main <- x$main
interactions <- x$interactions

plot_functions_4cov_nl(x = ls$fun)

x<-ls$fun
main = NULL
select = "all"
interactions = NULL
title.main = "Covariate effects (main)"
title.int = "Covariate effects (interactions)"
contour.bins = 11
fun_theme = NULL



  
options(warn =0)







  p.select <-
  plots$f.int12 + plots$f.int13 + plots$f.int14 +
  plots$f.int23 + plots$f.int24 + plots$f.int34 +
  plot_layout(design = layout.effects, guides = "collect") &
  ls_theme


  p.select <-
  plots$treatment + plots$error +
  plots$f.z1 + plots$f.z2 + plots$f.z3 + plots$f.z4 +
  plots$f.int12 + plots$f.int13 + plots$f.int14 +
  plots$f.int23 + plots$f.int24 + plots$f.int34 +
  plot_layout(design = layout.effects.all, guides = "collect") &
  ls_theme


  plots$treatment + plots$error + plots$f.z1 + plots$f.z2 + plots$f.z3 + plots$f.z4 +
  plot_layout(design = layout2, guides = "collect") &
  ls_theme

  p2

  p3 <-
  plots$response + plot_layout(design = "A##", guides = "keep") +
  ls_theme

  pw <- wrap_plots(p1, p2, p3, ncol = 1, heights = c(2.06, 2.06, 1))

  return(pw)
}
}


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


# x.dim = 100
# y.dim = 100
# imbalance = 0.3
# # imbalance.tol = NULL
# imbalance.tol = 0.01
# area.prop = 0.5
# # area.tol = NULL
# area.tol = c(0.4, 0.6)
# area.exact = FALSE
# score = scores$mahalanobis
# score.sam = 1e4
# seg.seed = 3
# seg.res = 5
# seg.min.dist = 5
# seg.min.area = 100
# seg.even = 2
# seg.prec = 0.1
# min.bound.dist = 0.1
# verbose = TRUE
# opt.imp.imb = 5
# opt.imp.area = 1
# opt.pop = 100
# opt.prec = 1e-4
# opt.pcrossover = 0.9
# opt.pmutation = 0.3
# opt.max.iter = 100
# opt.run = 25
# opt.parallel = 4
# opt.fine = TRUE
# opt.fine.max.iter = 100
# opt.fine.constr = TRUE
# opt.fine.tol = 1e-6


x.dim = 100
y.dim = 100
imbalance = 0.5
# imbalance.tol = NULL
imbalance.tol = 0.01
area.prop = 0.5
# area.tol = NULL
area.tol = c(0.4, 0.6)
area.exact = FALSE
score = scores$qprod
score.sam = 1e4
seg.seed = 3
seg.res = 5
elitism = 20
seg.min.dist = 3
seg.min.area = 100
seg.even = 2
seg.prec = 0.1
min.bound.dist = 0
verbose = TRUE
opt.imp.imb = 5
opt.imp.area = 1
opt.pop = 100
opt.prec = 1e-4
opt.pcrossover = 0.9
opt.pmutation = 0.5
opt.max.iter = 500
opt.run = 100
opt.parallel = 4
opt.fine = TRUE
opt.fine.max.iter = 100
opt.fine.constr = TRUE
opt.fine.tol = 1e-6


trt.ar <-
generate_areas_poly(
x.dim = 100,
y.dim = 100,
imbalance = 0.1,
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
elitism = 20,
seg.min.dist = 3,
seg.min.area = 100,
seg.even = 2,
seg.prec = 0.1,
min.bound.dist = 0,
verbose = TRUE,
opt.imp.imb = 5,
opt.imp.area = 1,
opt.pop = 100,
opt.prec = 1e-4,
opt.pcrossover = 0.9,
opt.pmutation = 0.5,
opt.max.iter = 500,
opt.run = 100,
opt.parallel = 4,
opt.fine = TRUE,
opt.fine.max.iter = 100,
opt.fine.constr = TRUE,
opt.fine.tol = 1e-6)

trt.ar$table$poly
plot(trt.ar$shape)

ar.shape <- trt.ar$shape
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



ref.poly <-
  trt.ar$shape[trt.ar$shape$type == "reference", "geometry"] |>
  st_cast("POLYGON")

trt.poly <-
  trt.ar$shape[trt.ar$shape$type == "treatment", "geometry"] |>
  st_cast("POLYGON")


sam.dens <- 0.5


treat.l <- list()
for(i in 1:nrow(trt.poly)) {
  sin.cen <-
    with(trt.poly[i,],
         st_sample(geometry, max(20, round(sam.dens * st_area(geometry)))))
  sin.cen.dt <-
    cbind(st_coordinates(sin.cen),
          st_distance(sin.cen, st_cast(trt.poly[i,], "MULTILINESTRING"))) |>
    as.data.table() |>
    setNames(c("x", "y", "dist"))
  trt.poly.l <- list()
  if(all(sin.cen.dt$dist < 0.5)) next
  for(j in 1:nrow(sin.cen.dt)) {
    scale.max <- sin.cen.dt[j, dist]
    if(scale.max < 0.5) next
    scale1 <- runif(1, 0.5, scale.max)
    scale2 <- runif(1, max(0.5 * scale1, 0.5), min(2 * scale1, scale.max))
    scales <- sample(c(scale1, scale2))
    trt.poly.l[[j]] <-
      generate_sin_gradient(x.dim = x.dim,
                            y.dim = y.dim,
                            centre = sin.cen.dt[j, c(x, -y+y.dim)],
                            x.scale = scales[1],
                            y.scale = scales[2],
                            rescale = c(0.25, runif(1, 0.5, 2)),
                            name = paste0("trt", j))
  }
  treat.l[[i]] <-
    do.call(c, trt.poly.l[!unlist(lapply(trt.poly.l, is.null))]) |>
    merge(name = "treatment") |>
    st_apply(1:2, FUN = mean, na.rm = TRUE) |>
    scale_int(c(0, runif(1, 0.5, 2)))
}




zplot(dist.rast)


plot(trt.ar$shape)

trt <- generate_treatment(x.dim = x.dim, y.dim = y.dim,
                          st_cast(st_geometry(trt.ar$shape[2,]), "POLYGON"),
                          mat.scale = 20)
# trt <- generate_treatment(x.dim = x.dim, y.dim = y.dim, nu =1.5, mat.scale = 0.05,
#                           trt.ar$shape[2,])
ggplot() +
  geom_stars(data = trt) +
  geom_sf(data = trt.ar$shape, fill = NA, colour = "white") +
  scale_fill_viridis_c()



trt.ar$table

x <- seq(0, 1, 0.01)
plot(x, plogis(x, 0.5, 1/ (2*pi)))

?plogis

[trt.poly[i,]] <- 0.25

?plogis

|>
zplot(n = 1)

summary(dist.dt$damp)

treat <-
  do.call(c, treat.l[!unlist(lapply(treat.l, is.null))]) |>
  merge(name = "treatment") |>
  st_apply(1:2, FUN = sum, na.rm = TRUE) |>
  scale_int()


zplot(treat)




generate_treatment(x.dim, y.dim, trt.poly, mat.nu = 1, mat.var = 0.1, mat.scale = 25, name = "stuff") |>
zplot()

  with(dist.mat.dt[poly != 0, dist.lim, .(dist.lim, treatment)],
       plot(dist.lim, treatment))

  dist.mat.dt[poly != 0, .(dist.lim, treatment)] |>
  ggplot(aes(x = dist.lim, y = treatment)) +
    geom_smooth()


  dist.mat.dt[, .(x, y, treatment)] |>
  st_as_stars()|>
  zplot()

  dist.dt[poly == 0, damp := 0]

  x <- seq(0, 1, 0.01)
  y <- plogis(x, location = 0.2, scale = 1/(4*pi))
  plot(x,y, ylim = c(0, 1))

  plogis


  dist.dt[, damp := mat.m * (damp/(damp.trt.m - damp.ref.m))]
  dist.dt[, damp := damp/mean(damp)]

  dist.dt[poly != 0, mean(damp)]
  dist.dt[poly != 0, .N]

  dist.dt[, damp := (effect.mean / mat.m)*(damp/(damp.trt.m - damp.ref.m))]

  dist.rast <- st_as_stars(dist.dt[, .(x, y, damp)])
  

  treat <- matern * dist.rast
  names(treat) <- "treatment"
 

  treat <- treat - min(treat$treatment, na.rm = TRUE)

  effect.mean / mat.m

  zplot(treat)

  treat.dt <-
    as.data.table(treat) |>
    merge(trt.dt, by = c("x", "y"))

  trt.mean <- treat.dt[poly != 0, var(treatment, na.rm = TRUE)]
  trt.mean <- treat.dt[poly != 0, mean(treatment, na.rm = TRUE)]


  ref.mean <- treat.dt[poly == 0, mean(treatment, na.rm = TRUE)]

  eff.scale <-
    effect.mean /
    (trt.mean - ref.mean)
  
  treat <- treat * eff.scale

  mean(treat[-trt.poly][["treatment"]], na.rm = TRUE)
  
  

  zplot(treat)


  treat <-
    do.call(c, treat.l[!unlist(lapply(treat.l, is.null))]) |>
    merge(name = "treatment") |>
    st_apply(1:2, FUN = sum, na.rm = TRUE) |>
    scale_int()


  }










(generate_matern(x.dim = 100, y.dim = 100, scale = min(x.dim, y.dim)/5) * dist.rast) |>
zplot(n = 1024)

plot(trt.poly)


generate_matern(x.dim = 100, y.dim = 100, scale = min(x.dim, y.dim)/5) |>
zplot()

# save.image("poly.RData")

load("poly.RData")

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

