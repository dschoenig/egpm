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



data <- ls$landscape
ls$landscape[, type := factor(type, levels = levels(type), ordered = TRUE)]
data[, zone := fifelse(x < 50, 1, 2)]

som.fit <-
  egp_som(ls$landscape,
          x.dim = 20,
          y.dim = 20,
          epochs = 1000,
          vars = c("z1", "z2", "z3", "z4"),
          parallel = 4)



data[, c("som.unit", "som.x", "som.y") := get_bmu(som.fit, coord = TRUE, list = TRUE)]


ls$landscape[, type := factor(type, levels = levels(type), ordered = TRUE)]

mod.gam <- 
  bam(response.int ~
        s(x, y, bs = "gp", m = 3, k = 100) +
        s(x, y, by = type, bs = "gp", m = 3, k = 100) + 
        s(som.x, som.y, bs = "gp", m = 3, k = 100),
      discrete = TRUE,
      data = ls$landscape)

summary(mod.gam)


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
  egp_posterior_draw(mod.gam, 1000, unconditional = TRUE, package = "mgcv")


pred <-
  egp_posterior_predict(
                        model = mod.gam,
                        data = data,
                        id.var = "cell",
                        # pred.name = "response",
                        # ids = sam2,
                        posterior = mod.post)

data[, cell2 := cell]

cf <-
  egp_define_counterfactual(
                            data = data,
                            fac.ids = data[type == "treatment", cell],
                            cf.ids = data[type == "reference", cell],
                            # compare.by = "zone",
                            group.by = list(NULL, "poly"),
                            som.var = "som.unit",
                            id.var = "cell",
                            deg.max = NULL,
                            n.min = 1)


units <- egp_summarize_units(pred, cf)
fac <- egp_evaluate_factual(pred, cf, agg.size = 1e2)
count <- egp_evaluate_counterfactual(pred, cf, units)
mar <- egp_marginal(fac, count)

mar[, .(mean(ATE), sd(ATE)), by =  .(group.id = group.id)]

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
    mod.nn.mh.re <- lm(response.int ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.mh.re)


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
