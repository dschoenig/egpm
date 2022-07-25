library(NLMR)
library(landscapetools)
library(RandomFields)
library(RandomFieldsUtils)
library(MatchIt)
library(raster)
library(sf)
library(stars)
library(spdep)
library(igraph)
library(zoo)
library(data.table)
library(mvnfast)

RFoptions(install="no")

source("~/projects/fcne_analysis/src/utilities.R")
source("utilities.R")

x.dim <- 100
y.dim <- 100

fbm.alpha = runif(1, 0.5, 1.5)
fbm.var = runif(1, 0.1, 10)
fbm.scale = runif(1, 0.1, 10)
grad.phi = runif(1, 0, 2*pi)
fbm.alpha
fbm.var
fbm.scale
z1 <- generate_z1(x.dim = 100,
                  y.dim = 100,
                  fbm.alpha,
                  fbm.var,
                  fbm.scale,
                  # fbm.scale = 0.2 ,
                  # fbm.alpha = runif(1, 0.5, 1.5),
                  # fbm.var = runif(1, 0.1, 10),
                  # fbm.scale = runif(1, 0.1, 10),
                  fbm.w = 0.9,
                  grad.phi = grad.phi,
                  grad.w = 0.1,
                  # rescale = c(0,1),
                  name = "z1")
zplot(z1)


 mod <- RMfbm(alpha = 1.5, var = 10, scale = 10)
  fbm <- 
    RFsimulate(model = mod, x = 1:200, y = 1:200) |>
    as.matrix() |>
    t() |>
    matrix2stars() |>
    zplot()


data.table(a = 1:5,
           b = a-5)


z1 <- generate_z1(x.dim = x.dim,
                  y.dim = y.dim,
                  fbm.alpha = 1.5,
                  fbm.var = 0.1,
                  fbm.scale = 0.1,
                  fbm.w = 0.9,
                  grad.phi = pi/3,
                  grad.w = 0.1,
                  rescale = c(0,1),
                  name = "z1")
zplot(z1)

z2 <- generate_z2(
                  x.dim = x.dim,
                  y.dim = y.dim,
                  fbm.alpha = 1,
                  fbm.var = 1,
                  fbm.scale = 1,
                  fbm.w = 0.2,
                  grad.phi = runif(1, 0, 2*pi),
                  grad.shift = 1,
                  grad.w = 0.8,
                  rescale = c(0,1),
                  name = "z2")
zplot(z2)

z3 <- generate_z3(x.dim = x.dim,
                  y.dim = y.dim,
                  dist.n = round(runif(1, 5, 25)),
                  dist.w = 0.9,
                  grad.phi = pi/3,
                  grad.w = 0.1,
                  acc = 0.1,
                  rescale = c(0,1),
                  name = "z3")
zplot(z3)

z4 <- generate_z4(x.dim = x.dim,
                  y.dim = y.dim,
                  seg.n = round(runif(1, 10, 50)),
                  seg.nu = runif(1, 1, 2),
                  seg.var = runif(1, 0.1, 1),
                  seg.scale = runif(1, 0.1, 1),
                  seg.w = 0.9,
                  grad.phi = pi/3,
                  grad.w = 0.1,
                  rescale = c(0,1),
                  name = "z4")
zplot(z4)

split <- generate_split(x.dim = x.dim, 
                        y.dim = y.dim,
                        n = 200,
                        prop = 0.5,
)
type <- st_rasterize(split, template = generate_empty(x.dim = x.dim, y.dim = y.dim))
type <- generate_empty(x.dim, y.dim) + type
type <- setNames(type, "type")
# zplot(type, 2)

treatment <-
  generate_treatment(x.dim = x.dim, y.dim = y.dim,
                     effect.size.sp = 0.9,
                     effect.size.bd = 0.0,
                     nuclei = 10,
                     poly = split,
                     phi.range = c(-pi/4, pi/4),
                     shift.range = c(-0.25, 0.25),
                     x.scale.range = c(0.5, 1.5),
                     y.scale.range = c(0.25, 0.5),
                     nuc.eff.range = c(1, 2),
                     name = "treatment")
zplot(treatment)

  effect.size.sp = 0.9
                     effect.size.bd = 0.1
                     nuclei = 5
                     poly = split
                     phi.range = c(-pi/4, pi/4)
                     shift.range = c(-0.25, 0.25)
                     x.scale.range = c(1, 1.5)
                     y.scale.range = c(0.25, 0.5)
                     nuc.eff.range = c(1, 2)
                     name = "treatment"



# zplot(treatment)


covariates <- c(z1, z2, z3, z4)

effect.size <- round(runif(4, 0.5, 2) * sample(c(1, -1), 4, replace = TRUE), 1)
effects <- list()
for(i in seq_along(effect.size)){
  effects[[i]] <-
    matrix2stars(covariates[[i]] * effect.size[i])
}

effects <-
  do.call(c, effects) |>
  setNames(paste0("f", 1:4))

mod.error <-
  # RMmatern(nu = 1, var = 0.1, scale = 100) +
  RMexp(var = 0.1, scale = 100) +
  RMnugget(var = 0.1)
error.sp <- RFsimulate(mod.error, x = 1:x.dim, y = 1:y.dim)
error.sp <-
  generate_empty(x.dim, y.dim) + st_as_stars(error.sp)
error.sp <- setNames(error.sp, "error")

error.rand <-
  matrix(rnorm(prod(x.dim, y.dim), 0, 0.1)) |>
  matrix2stars()

error <- error.sp + error.rand

# error <- generate_fbm(200, 200, 0.5, rescale = c(-0.5, 0.5))

# error <-
#   matrix(rnorm(200, 0, 0.1), 200, 200) |>
#   matrix2stars()

response <-
  c(effects, treatment, error) |>
  merge() |>
  st_apply(1:2, sum) |>
  setNames("response")

# zplot(response)
# zplot(treatment)

landscape <- c(response, type, covariates, treatment, effects, error)


land <- as.data.frame(landscape)
setDT(land)

land[, type := factor(ifelse(type == 0, "control", "treatment"), levels = c("control", "treatment"))]
land$id <- 1:nrow(land)

library(mgcv)

sam <- sample(1:nrow(land), 400)


library(kohonen)
som.dim <- 20
cov.z <- scale(land[sam, .(z1, z2, z3, z4)],
               center = TRUE, scale = TRUE)
grid <- somgrid(xdim = som.dim, ydim = som.dim, 
                topo = "rectangular", 
                neighbourhood.fct = "gaussian")
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
    land[, .(z1, z2, z3, z4)] |>
    # land[sam, .(z1, z2, z3, z4)] |>
    scale_data_som(som = som.fit) |>
    embed_som(som = som.fit,
              grid.coord = TRUE)
land[,
         `:=`(som_bmu = mapped$bmu[,1],
              som_x = mapped$grid.coordinates$bmu.1[,"x"],
              som_y = mapped$grid.coordinates$bmu.1[,"y"])
         ]


land[sam, .(mean(treatment), mean(response)), type]

merge(landscape) |> zplot()


mod <- gam(response ~ type,
           data = land[sam,],
           select = TRUE, method = "REML")
summary(mod)

mod <- gam(response ~
           s(x, y, k = 30) + type,
           data = land[sam,],
           select = TRUE, method = "REML")
summary(mod)

mod <- gam(response ~
           s(x, y, by = type, bs = "tp", k = 50) + type,
           data = land[sam,],
           select = TRUE, method = "REML")
summary(mod)


mod <- gam(response ~ 
           s(x, y, by = type, bs = "gp", k = 30) + s(z1) + s(z2) +s(z3) + s(z4),
           data = land[sam,], select = TRUE, method = "REML", optimizer = "efs")
summary(mod)
AIC(mod)

mod <- gam(response ~
           s(x, y, by = type, bs = "gp", k = 30) + s(som_x, som_y, bs = "gp", k = 20),
           data = land[sam,], select = TRUE, method= "REML", optimizer = "efs")
summary(mod)
AIC(mod)

X <- predict(mod, land[sam,], type = "lpmatrix")
post <- rmvn(1000, coef(mod), vcov(mod, unconditional = TRUE))
sel <- grep("type", names(coef(mod)))
post[, -sel] <- 0
lp <- X %*% t(post)
id.treat <- which(land[sam,]$type == "treatment")
mean(apply(lp[id.treat,], 2, mean))
Mode(apply(lp[id.treat,], 2, mean))
quantile2(apply(lp[id.treat,], 2, mean))

library(MatchIt)
library(lme4)

matched.cem <- matchit(type ~ z1 + z2 + z3 + z4, method = "cem", data = land)
md.cem <- match.data(matched.cem)
lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.cem) |>
# lm(response ~ type, weights = weights, data = md.cem) |>
summary()

matched.ps <- matchit(type ~ z1 + z2 + z3 + z4, method = "nearest", distance = "glm", data = land)
md.ps <- match.data(matched.ps)
lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.ps) |>
# lm(response ~ type, weights = weights, data = md.ps) |>
summary()



library(ggplot2)
ggplot(md.cem) +
  geom_raster(aes(x = x , y = y, fill = type))


data.bl <- assign_bl_som(land[sam],
                         som = som.fit,
                         cov.col = c("z1", "z2", "z3", "z4"))

post <- rmvn(1000, coef(mod), vcov(mod, unconditional = TRUE))
colnames(post) <- names(coef(mod))
post <- as_draws_matrix(post)

lp <-
  evaluate_posterior(model = mod,
                     posterior = post,
                     newdata = land[sam,],
                     id.col = "id",
                     type = "link",
                     progress = TRUE)

groups.type <-
    land[sam,] |>
    ids_by_group(id.col = "id", group.vars = "type")
id.list.type <- groups.type$ids
names(id.list.type) <- groups.type$type
yhat.type <- aggregate_variables(lp,
                    agg.fun = E,
                    ids = id.list.type,
                    )

groups.som <-
    land[sam,] |>
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

eff.arc.gam <- eff.arc


data.bl <- assign_bl_som(land,
                         som = som.fit,
                         cov.col = c("z1", "z2", "z3", "z4"))

post <- rmvn(1000, coef(mod), vcov(mod, unconditional = TRUE))
colnames(post) <- names(coef(mod))
post <- as_draws_matrix(post)

lp <-
  evaluate_posterior(model = mod,
                     posterior = post,
                     newdata = land,
                     id.col = "id",
                     type = "link",
                     progress = TRUE)

groups.id <-
    land |>
    ids_by_group(id.col = "id", group.vars = "id")
id.list.id <- groups.id$ids

groups.som <-
    land |>
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
  lapply(land$id,
         \(x) {
               extract_weights(ids.units[.(x)],
                               w.col = "som_bmu.bl.w",
                               by.col = "som_bmu.bl",
                               standardize = TRUE)
              })

eff.bl.id <- reweigh_posterior(yhat.som, w = w.points)
colnames(eff.bl.id) <- 1:40000

eff.id <- arc(lp, eff.bl.id)
lp.dt <- as_draws_df(eff.id) |> as.data.table()
lp.dt <-
  melt(lp.dt[, -c(".chain", ".iteration", ".draw")],
     variable.name = "id", value.name = "treatment.pred")
lp.dt <- lp.dt[, .(treatment.pred = mean(treatment.pred)), id]
lp.dt[,id := as.numeric(as.factor(id))]

lp.dt <-
  merge(land[,.(id, x, y, treatment)], lp.dt) |>
  melt(id.vars = c("id", "x", "y"))



library(ggplot2)

ggplot(lp.dt, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  facet_wrap(vars(variable)) +
  scale_fill_continuous_divergingx("Roma")

library(colorspace)

library(brms)

modb <- brm(response ~
           s(x, y, by = type, bs = "gp", k = 30) + s(som_x, som_y, bs = "gp", k = 20),
           data = land[sam,], cores = 4)

pr <- get_prior(response ~
           s(x, y, by = type, bs = "gp", k = 30) + s(som_x, som_y, bs = "gp", k = 20),
           data = land[sam,])

lp <- posterior_predict(modb, land[sam,])
colnames(lp) <- land[sam, id]
lp <- as_draws_matrix(lp)

groups.type <-
    land[sam,] |>
    ids_by_group(id.col = "id", group.vars = "type")
id.list.type <- groups.type$ids
names(id.list.type) <- groups.type$type
yhat.type <- aggregate_variables(lp,
                    agg.fun = E,
                    ids = id.list.type,
                    )

groups.som <-
    land[sam,] |>
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

summary(eff.arc.gam)

