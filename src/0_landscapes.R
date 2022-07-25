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

RFoptions(install="no")

source("~/projects/fcne_analysis/src/utilities.R")


equilibrate <- function(x, ...) {
  UseMethod("equilibrate", x)
}

equilibrate.matrix <- function(x, tolerance = .Machine$double.eps) {
  while(any(abs(rowSums(x) - 1) > .Machine$double.eps) |
        any(abs(colSums(x) - 1) > .Machine$double.eps)) {
    x <- apply(x, 2, \(x) x / sum(x))
    x <- apply(x, 1, \(x) x / sum(x))
  }
  return(x)
}


rescale_uint <- function(x, ...) {
  UseMethod("rescale_uint", x)
}

rescale_uint.stars <- function(x) {
  n.att <- length(x)
  for(i in 1:n.att) {
    min.att <- min(x[[i]], na.rm = TRUE)
    max.att <- max(x[[i]], na.rm = TRUE)
    x[[i]] <- (x[[i]] - min.att) / (max.att - min.att)
  }
  return(x)
}

scale_int <- function(x, ...) {
  UseMethod("scale_int", x)
}

scale_int.stars <- function(x, int = c(0, 1)) {
  n.att <- length(x)
  for(i in 1:n.att) {
    min.att <- min(x[[i]], na.rm = TRUE)
    max.att <- max(x[[i]], na.rm = TRUE)
    x[[i]] <- int[1] * (int[2] - int[1]) * (x[[i]] - min.att) / (max.att - min.att)
  }
  return(x)
}

matrix2stars <- function(x, res = 1, val = "value", ...) {
  x.stars <-
    st_as_stars(x) |>
    st_set_dimensions(names = c("x", "y")) |>
    st_set_dimensions("x", offset = 0, delta = res) |>
    st_set_dimensions("y", offset = nrow(x) * res, delta = - res) |>
    setNames(val)
  return(x.stars)
}


generate_gradient <- function(xdim, ydim, phi) {
  mat <- outer(1:x.dim, 1:y.dim, FUN = \(x,y) x * cos(phi) - y * sin(phi))
  mat.scaled <- 
    (mat - min(mat))/ (max(mat) - min(mat))
  r <- matrix2stars(mat.scaled) 
  return(r)
}

n.cov <- 5
n.t <- 1
x.dim <- 200
y.dim <- 200
# n.poly <- (x.dim * y.dim) / 100
n.poly <- 200
prop.treat <- 0.5
n.poly.treat <- round(prop.treat * n.poly, 0)
cov.alpha <- 2 * runif(n.cov)
# cov.phi <- sample((2*pi)/n.cov * 1:n.cov, n.cov)
cov.phi <- runif(5) * (2*pi)


# cov.combine <- matrix(0, ncol = n.cov, nrow = n.cov)
# diag(cov.combine) <- 0.5
# idx.bands <- 
#   (row(cov.combine) == col(cov.combine) + 1) |
#   (row(cov.combine) == col(cov.combine) - 1)
# cov.combine[idx.bands] <- 0.25
# cov.combine[1, n.cov] <- 0.25
# cov.combine[n.cov, 1] <- 0.25

# cov.combine <- matrix(runif(n.cov^2), n.cov, n.cov)
# cov.combine <- equilibrate(cov.combine)
# cov.combine <-
#   cov.combine *
#   matrix(sign(runif(n.cov^2, -1, 1)), n.cov, n.cov)

cov.combine <- diag(n.cov)

# t.combine <- matrix(1/n.cov, nrow = n.cov, ncol = n.t)
t.combine <- matrix(runif(n.cov, 0, 1), nrow = n.cov, ncol = n.t)

fields <- list()
covariates <- list()


# cov.mod <-
#   RMmatern(nu = 1.5, var = 0.5, scale = 10)

for(i in 1:n.cov){
  cov.mod <-
    RMfbm(alpha = cov.alpha[i])
  l.mat <-
    RFsimulate(model = cov.mod, x = 1:x.dim, y = 1:y.dim, grid = TRUE) |>
    as.matrix() |>
    matrix2stars()
  l.pg <- generate_gradient(x.dim, y.dim, phi = cov.phi[i])
  l.merged <-
    rescale_uint(0.8 * l.mat + 0.2 * l.pg)
    # rescale_uint(l.mat)
  fields[[i]] <- l.merged
}

fields <- 
  do.call(c, fields) |>
  setNames(paste0("field", 1:length(fields)))


for(i in 1:n.cov){
  w <- cov.combine[,i]
  fields.w <- list()
  for(j in seq_along(w)){
    fields.w[[j]] <-
      matrix2stars(fields[[j]] * w[j])
  }
  covariates[[i]] <-
    Reduce("+", fields.w) |>
    rescale_uint()
}

covariates <- 
  do.call(c, covariates) |>
  setNames(paste0("X", 1:length(covariates)))


fields.combined <-
  merge(fields)|>
  st_apply(1:2, sum) |>
  setNames("fields.combined") |>
  rescale_uint()

poly.lim <- 
  st_bbox(fields.combined) |>
  st_as_sfc()

poly.geom <- 
  poly.lim |>
  st_sample(n.poly) |>
  st_union() |>
  st_voronoi() |>
  st_cast() |>
  st_intersection(poly.lim)

poly.ex <- 
  st_extract(fields.combined, poly.geom, FUN = mean) |>
  st_as_sf()
names(poly.ex) <- c("fields.mean", "geometry")
poly.ex$id <- 1:n.poly
poly.ex$weight <- 
  with(poly.ex, fields.mean/sum(fields.mean))


Aw <-
  poly2nb(poly.ex) |>
  nb2listw(style = "B") |> 
  listw2mat() |>
  apply(2, \(x) x * poly.ex$weight)

G <- graph_from_adjacency_matrix(Aw, weighted = TRUE)

# One chunk area
path <- eulerian_path(G)$vpath
upath <- unique(path)
weights.upath <- poly.ex$weight[upath]
chunk.start <- 1:n.poly.treat
chunk.end <- (n.poly - n.poly.treat + 1):n.poly
if(sum(weights.upath[chunk.start]) > 
   sum(weights.upath[chunk.end])) {
  idx.treat <- upath[chunk.start]
} else {
  idx.treat <- upath[chunk.end]
}


# # Random walk
# walk.start <- with(poly.ex, sample(id, 1, prob = weight))
# walk <- random_walk(G, walk.start, n.poly.treat)
# while(length(unique(walk)) < n.poly.treat) {
#   walk <- c(walk, random_walk(G, walk[length(walk)], n.poly.treat))
# }
# idx.treat <- unique(walk)[1:n.poly.treat]

# Sample
# idx.treat <- with(poly.ex, sample(id, n.poly.treat, prob = weight))


poly.ex$treatment <- 0
poly.ex$treatment[idx.treat] <- 1

poly.control <- subset(poly.ex, treatment == 0) |> st_union()
poly.treat <- subset(poly.ex, treatment == 1) |> st_union()

treatment <-
  st_rasterize(poly.ex[, "treatment"], fields.combined)
plot(treatment)


effects <- runif(5, -2, 2)

fx <- list()
for(j in seq_along(effects)){
  fx[[j]] <-
    matrix2stars(covariates[[j]] * effects[j])
}

fx <-
  do.call(c, fx) |>
  setNames(paste0("fx", 1:n.cov))


covariates.dt <-
covariates |>
as.data.frame() |>
as.data.table()


fx.dt <-
fx |>
as.data.frame() |>
as.data.table()


bg.mod <-
  RMmatern(nu = 1.5, var = 0.1, scale = 10)
bg <-
  RFsimulate(model = bg.mod, x = 1:x.dim, y = 1:y.dim, grid = TRUE) |>
  as.matrix() |>
  matrix2stars()

coord <- st_centroid(c(poly.control, poly.treat)) |> st_coordinates()
phi <- atan2((coord[2,2] - coord[1,2]), (coord[2,1] - coord[1,1]))
treat.grad <- generate_gradient(x.dim, y.dim, phi = phi)
treat.mat <-
  RFsimulate(model = RMnugget(var = 0.1), x = 1:x.dim, y = 1:y.dim, grid = TRUE) |>
  as.matrix() |>
  matrix2stars()
# treat.mat <- 0

treat.c <- rescale_uint(treat.grad + treat.mat)

plot(treat.c)

c.bg <- mean(bg[poly.control]$value, na.rm = TRUE)
t.bg <- mean(bg[poly.treat]$value, na.rm = TRUE)
c.grad <- mean(treat.c[poly.control]$value, na.rm = TRUE)
t.grad <- mean(treat.c[poly.treat]$value, na.rm = TRUE)

f.treatment <-
(1 + c.bg - t.bg ) / (t.grad - c.grad) * treat.c + bg

f.treatment <-
  f.treatment - mean(f.treatment[poly.control]$value, na.rm = TRUE) + 0.5

mean(f.treatment[poly.control]$value, na.rm = TRUE)
mean(f.treatment[poly.treat]$value, na.rm = TRUE)

f.treatment.dt <-
  f.treatment |>
  setNames("f.treatment") |>
  as.data.frame() |>
  as.data.table()

e.mod <-
RMmatern(nu = 1.5, var = 0.01, scale = 10) +
RMnugget(var = 0.01)
e <-
  RFsimulate(model = e.mod, x = 1:x.dim, y = 1:y.dim, grid = TRUE) |>
  as.matrix() |>
  matrix2stars() |>
  setNames("error") |>
  as.data.frame() |>
  as.data.table()

# e$error <- e$error * (0.1 / sd(e$error))

full <- 
as.data.table(as.data.frame(treatment)) |>
merge(covariates.dt) |>
merge(fx.dt) |>
merge(f.treatment.dt)
full$e <- e$error
# full[, response := f.treatment +  fx1^2 + exp(fx2) + fx3 + fx4 + fx5 + e] 
full[, response := f.treatment +  fx1 + fx2 + fx3 + fx4 + fx5 + e] 
full[, response0 := f.treatment + e]
full[, type := factor(fifelse(treatment == 0, "control", "treatment"),
                      levels = c("control", "treatment"))]
full[, type.o := factor(fifelse(treatment == 0, "control", "treatment"),
                      levels = c("control", "treatment"), ordered = TRUE)]

# full[, pi := inv_logit(f.treatment + fx1^2 + exp(fx2) + fx3 + fx4 + fx5)]
full[, pi := inv_logit(f.treatment + fx1^2 + exp(fx2) + fx3 + fx4 + fx5)]
full[, response.bin := rbinom(pi, 1, pi)]

sam <- sample(1:nrow(full), 1500)
full.sam <- full[sam,]

library(kohonen)

cov.z <- scale(full[, .(X1, X2, X3, X4, X5)],
               center = TRUE, scale = TRUE)
grid <- somgrid(xdim = 10, ydim = 10, 
                topo = "rectangular", 
                neighbourhood.fct = "gaussian")
som.fit <- som(cov.z,
           grid = grid, 
           rlen = 200,
           radius = 100,
           init = init_som(cov.z, 10, 10),
           mode = "pbatch", 
           cores = 4,
           normalizeDataLayers = FALSE)
som.fit$scale <- list(mean = attr(cov.z, "scaled:center"),
                  sd = attr(cov.z, "scaled:scale"))
mapped <- 
    full.sam[, .(X1, X2, X3, X4, X5)] |>
    scale_data_som(som = som.fit) |>
    embed_som(som = som.fit,
              grid.coord = TRUE)
full.sam[,
         `:=`(som_bmu = mapped$bmu[,1],
              som_x = mapped$grid.coordinates$bmu.1[,"x"],
              som_y = mapped$grid.coordinates$bmu.1[,"y"])
         ]

library(mgcv)

full.sam[,mean(f.treatment), type]
full.sam[,sd(f.treatment), type]

mod <- gam(response0 ~ type,
           data = full.sam, method = "REML", select = FALSE)
summary(mod)

mod <- gam(response ~ type,
           data = full.sam, method = "REML", select = FALSE)
summary(mod)

mod <- gam(response ~ type + X1 + X2 + X3 + X4 + X5,
           data = full.sam, method = "REML", select = FALSE)
summary(mod)

mod <- gam(response ~ s(x, y),
           data = full.sam, method = "REML", select = TRUE)
summary(mod)

mod <- gam(response ~ s(x, y, by = type),
           data = full.sam, method = "REML", select = TRUE)
summary(mod)

mod <- gam(response ~ s(x, y, by = type, k = 30) + X1 + X2 + X3 + X4 + X5,
           data = full.sam, method = "REML", select = TRUE)
summary(mod)

mod <- gam(response ~
           s(x, y, k = 30) +
           s(x, y, by = type.o, k = 30) + s(som_x, som_y),
           data = full.sam, method = "REML", select = TRUE)
summary(mod)

X <- full[, .(x, y, X1, X2)]

X@RMfixed(effects)

scale <-
RMcovariate(covariates$X1, raw = TRUE)

mod <-
RMexp()

RFsimulate(model = mod, 1:x.dim, 1:y.dim)


Ffctn(RMcovariate(data=z, x=1:10), c(2, 2.1, 2.5, 3))
?RFfctn

X <- predict(mod, full.sam, type = "lpmatrix")
post <- rmvn(1000, coef(mod), vcov(mod, unconditional = TRUE))
sel <- grep("treatment", names(coef(mod)))
post[, -sel] <- 0
lp <- X %*% t(post)
id.treat <- which(full.sam$treatment == 1)
mean(apply(lp[id.treat,], 2, mean))
Mode(apply(lp[id.treat,], 2, mean))
quantile2(apply(lp[id.treat,], 2, mean))

vis.gam(mod, view = c("som_x", "som_y"), se = 1)

plot(merge(covariates))

matched <- matchit(type ~ X1 + X2 + X3 + X4 + X5, method = "cem", data = full)
md <- match.data(matched)

library(lme4)

lmer(response ~ type + (1 + treatment | subclass) , data = md) |>
summary()

lm(response ~ type, data = md) |>
# lm(response ~ type + X1 + X2 + X3 + X4 + X5, data = md, weights = weights) |>
summary()


lm(response ~ type, data = full.sam) |>
# lm(response ~ type + X1 + X2 + X3 + X4 + X5, data = md, weights = weights) |>
summary()

library(glmmTMB)

full.sam$pos <- with(full.sam, numFactor(x, y))
full.sam$group <- factor(rep(1, nrow(full.sam)))


glmmTMB(response ~ type + X1 + X2 + X3 + X4 + X5 + mat(pos + 0 | group) ,
        data = full.sam) |>
summary()

glmmTMB(response.bin ~ type + X1 + X2 + X3 + X4 + X5 + (1 + treatment | subclass) , data = md, family = binomial) |>
summary()

glm(response.bin ~ type + X1 + X2 + X3 + X4 + X5, weights = weights  , data = md, family = binomial) |>
summary()

library(ggplot2)

ggplot(full) +
  geom_raster(aes(x = x , y = y, fill = response)) +
  scale_fill_viridis_c()


ggplot(full) +
  geom_raster(aes(x = x , y = y, fill = f.treatment)) +
  scale_fill_viridis_c()

ggplot(md) +
  geom_raster(aes(x = x , y = y, fill = treatment)) +
  scale_fill_viridis_c()




full[, mean(X5), type]

cor(full[, .(treatment, X1, X2, X3, X4, X5)])

plot(mod, all = TRUE)

full.sam$fitted <- fitted(mod)

full.sam[, mean(fitted), by = type]
effects

gam.check(mod)

library(ggplot2)

ggplot(full) +
  geom_raster(aes(x = x , y = y, fill = f.treatment)) +
  scale_fill_viridis_c()



gf <-
nlm_gaussianfield(x.dim, y.dim, mean = 0, autocorr_range = 25, mag_var = 0.1, rescale = TRUE) |>
# nlm_fbm(x.dim, y.dim, fract_dim = 1.5) |>
st_as_stars() |>
setNames("f")


x <- (1:x.dim) - 0.5
y <- (1:y.dim) - 0.5

grad.center <-
  st_union(sub.control) |>
  st_sample(3) |>
  st_coordinates()

dat <- expand.grid(x = x, y = y)

dat$f <-
with(dat,
     sample(c(-1,1), 1) * runif(1, 0.5, 1) * exp(-(((x-grad.center[1,1])/50)^2
                      + ((y - grad.center[1,2])/50)^2)) +
     sample(c(-1,1), 1) * runif(1, 0.5, 1) * exp(-(((x-grad.center[2,1])/50)^2
                      + ((y - grad.center[2,2])/50)^2)) +
      sample(c(-1,1), 1) * runif(1, 0.5, 1) * exp(-(((x-grad.center[3,1])/50)^2
                      + ((y - grad.center[3,2])/50)^2))) 

contr.s <-
st_as_stars(dat)


grad.center <-
  st_union(sub.treat) |>
  st_sample(3) |>
  st_coordinates()

dat <- expand.grid(x = x, y = y)

dat$f <-
with(dat,
     sample(c(1,1), 1) * runif(1, 0.5, 1) * exp(-(((x-grad.center[1,1])/50)^2
                      + ((y - grad.center[1,2])/50)^2)) +
     sample(c(1,1), 1) * runif(1, 0.5, 1) * exp(-(((x-grad.center[2,1])/50)^2
                      + ((y - grad.center[2,2])/50)^2)) +
      sample(c(1,1), 1) * runif(1, 0.5, 1) * exp(-(((x-grad.center[3,1])/50)^2
                      + ((y - grad.center[3,2])/50)^2))) 


treat.s <-
st_as_stars(dat)


ggplot(dat) +
  geom_raster(aes(x = x, y =y, fill = f)) +
  scale_fill_viridis_c()


treat.s <-
  treat.s *
  (
   (1 +
    mean(contr.s[sub.control]$f, na.rm = TRUE) -
    mean(contr.s[sub.treat]$f, na.rm = TRUE)) /
   (mean(treat.s[sub.treat]$f, na.rm = TRUE) -
    mean(treat.s[sub.control]$f, na.rm = TRUE))
  )


new.s <-
contr.s + treat.s


mean(new.s[sub.treat]$f, na.rm = TRUE) -
mean(new.s[sub.control]$f, na.rm = TRUE)


plot(new.s)

as.data.frame(new.s) |>
ggplot() +
  geom_raster(aes(x = x, y =y, fill = f)) +
  scale_fill_viridis_c()

f.treatment <- new.s





show_landscape(covariates)

Aw[1:5, 1:5]

rescale_uint(fields.combined)


?st_extract()

plot(poly.ex)


poly.ex$weight <- 
  with(poly.ex, layer/sum(layer))
poly.ex$id <- 1:n.poly


idx <- with(poly.ex, sample(id, 100, prob = weight^2))
poly.ex[poly.ex$id %in% idx,] |> st_union() |> plot()


library(spdep)


plot(poly.nb, st_centroid(st_geometry(poly.ex)))


plot(voronoi_tess)





treatments <- list()

for(i in 1:n.t) {
  w <- t.combine[,i]
  treatments[[i]] <- 
    mapply(\(x, y) x * y, w, fields) |>
    brick() |>
    calc(sum) |>
    util_rescale() |>
    util_classify(n = 5)
}

show_landscape(treatments[[1]])


landscapetools::util_classify.RasterLayer

?util_classify
r <- treatments[[1]]

clust <- nlm_mosaictess(1000, 1000, resolution = 1, germs = 100)
show_landscape(clust)

nlm_mosaictess


??nlm_

r[r == 1] <- 0
r[r == 2] <- 1

n.fun <- function(x) {
  s <- sum(x, na.rm = TRUE)
  if(s < 4) 0 else 1
}

n.fun(c(1:10))

for(i in 1:10) {
  r <- focal(r, matrix(c(1,1,1,1,0,1,1,1,1), ncol = 3), fun = n.fun) 
}

plot(r)

      ?focal


r[r != 1] <- NA

pol <-
st_as_sf(st_as_stars(r)) |>
st_union(is_coverage = TRUE) |>
st_cast("POLYGON") |>
st_as_sf()

pol$area <- st_area(pol)

pol[pol$area >= 100,] |>
st_union(is_coverage = TRUE) |>
plot()

rasterToPolygon(r)

library(stars)
library(sf)
st_as_stars(treatments[[1]])

nlm_planargradient

# Mode 1
util_classify(fractal_landscape,
              n = 3,
              level_names = c("Land Use 1", "Land Use 2", "Land Use 3"))

# Mode 2
util_classify(fractal_landscape,
              weighting = c(0.5, 0.25, 0.25),
              level_names = c("Land Use 1", "Land Use 2", "Land Use 3"))

# Mode 3
real_land <- util_classify(gradient_landscape,
              n = 3,
              level_names = c("Land Use 1", "Land Use 2", "Land Use 3"))

fractal_landscape_real <- util_classify(fractal_landscape, real_land = real_land)
fractal_landscape_mask <- util_classify(fractal_landscape, real_land = real_land, mask_val = 1)


com <- cbind(as.data.frame(fields[[1]]), as.data.frame(fields[[2]]), as.data.frame(fields[[3]]))
names(com) <- c("field1", "field2", "field3")

library(mgcv)

mod <- gam(field1 ~ s(field2, field3, k = 200), data = com)
summary(mod)
gam.check(mod)

vis.gam(mod, theta = 210)



cluster

library(stars)

f1s <- st_as_stars(fields[[1]])
f2s <- st_as_stars(fields[[2]])

0.2 * f1s + f2s

nlm_fbm

?RMfbm


library(RandomFields)

model <-
  RMatern(nu = 2.5) + # with variance 4 and scale 10
  RMnugget(var=20) + # nugget

## define the locations:
from <- 1
to <- 100
x.seq <- seq(from, to, length=100)
y.seq <- seq(from, to, length=100)


model <-
  RMmatern(nu = 3.5, var = 0.1, scale = 10)
   # with variance 4 and scale 10
  # RMnugget(var=0) # nugget
simu <- RFsimulate(model, x=x.seq, y=y.seq, grid = TRUE) |>
as.matrix() |>
matrix2stars()

plot(simu)

nlm_planargradient(100, 100) |> st_as_stars()

phi = 1.75 * pi


generate_gradient(100,100,phi) |> plot()

nlm_planargradient(100, 100, direction = (180/pi * phi)) |>
show_landscape()

nlm_distancegradient

library(stars)
st_as_stars(cov)



new <- simu + cov

plot(new[,1:50])
mean(new[,1:50]$value)
mean(new[,51:100]$value)

plot(new)
plot(cov)


?st_as_stars
st_as_stars(cov)

image(cov, col = hcl.colors(12))

mean(cov[,1:50])
mean(cov[,51:100])

n.poly <- 100
n.poly.treat <- 50

poly.lim <- 
  st_bbox(simu) |>
  st_as_sfc()

poly.geom <- 
  poly.lim |>
  st_sample(n.poly) |>
  st_union() |>
  st_voronoi() |>
  st_cast() |>
  st_intersection(poly.lim)

poly.ex <- 
  st_extract(simu, poly.geom, FUN = mean) |>
  st_as_sf()
names(poly.ex) <- c("fields.mean", "geometry")
poly.ex$id <- 1:n.poly
poly.ex$weight <- 1


Aw <-
  poly2nb(poly.ex) |>
  nb2listw(style = "B") |> 
  listw2mat() |>
  apply(2, \(x) x * poly.ex$weight)

G <- graph_from_adjacency_matrix(Aw, weighted = TRUE)

# One chunk area
path <- eulerian_path(G)$vpath
upath <- unique(path)
weights.upath <- poly.ex$weight[upath]
treat.start <-
  rollapply(weights.upath, n.poly.treat, sum, align = "right") |>
  which.max()
treat.end <- treat.start + n.poly.treat - 1
idx.treat <- upath[treat.start:treat.end]

poly.treat <- subset(poly.ex, id %in% idx.treat) |> st_union()
poly.control <- subset(poly.ex, !id %in% idx.treat) |> st_union()

a <- mean(simu[poly.control]$value, na.rm = TRUE)
b <- mean(simu[poly.treat]$value, na.rm = TRUE)
c <- mean(cov[poly.control]$value, na.rm = TRUE)
d <- mean(cov[poly.treat]$value, na.rm = TRUE)


new <-
(1 + a - b ) / (d -c) * cov + simu

new <- new - mean(new[poly.control]$value, na.rm = TRUE) + 0.5

mean(new[poly.control]$value, na.rm = TRUE)
mean(new[poly.treat]$value, na.rm = TRUE)

c(poly.treat, poly.control)

coord <- st_centroid(c(poly.control, poly.treat)) |> st_coordinates()
phi <- atan2((coord[2,2] - coord[1,2]), (coord[2,1] - coord[1,1]))
generate_gradient(100, 100, phi = phi)

plot(new)

plot(poly.treat)

tan(atan(0))
