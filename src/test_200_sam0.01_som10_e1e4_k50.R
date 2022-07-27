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

args <- commandArgs(trailingOnly = TRUE)
n.threads <- as.integer(args[1])
# n.threads <- 4

# RFoptions(cores = n.threads)

path.base <- "../"
path.tests <- paste0(path.base, "tests/")

if(!dir.exists(path.tests)) dir.create(path.tests, recursive = TRUE)

file.results <- paste0(path.tests, "test_200_sam0.01_som10_e1e4_k50.rds")

n <- 1000
# n <- 2
ls.dim <- 200
sam.frac <- 0.01
som.dim <- 10
som.rlen <- 10000
k.max <- 50
eval.som <- FALSE
egp.basis <- "gp"
egp.select <- TRUE

set.seed(19010511) # Rose Ausländer
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
             z1.fbm.w = 0.8,
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
             z3.acc = 0.1,
             z4.seg.n = sample(10:25, n, replace = TRUE),
             z4.mat.nu = runif(n, 1, 2),
             z4.mat.var = runif(n, 0.1, 1),
             z4.mat.scale = runif(n, 0.1, 1),
             z4.mat.w = 0.8,
             # z4.grad.phi = runif(n, -1, 1) * 2 * pi,
             z4.grad.phi = runif(n, -pi/8, pi/8) +
                           sample(c(0, pi), n, replace = TRUE),
             z4.grad.w = 0.2,
             split.n = 200,
             split.prop = runif(n, 0.4, 0.6),
             e.exp.var = 0.1,
             e.exp.scale = 100,
             e.nug.var = 0.1,
             e.rand.var = 0.1
             )


results <- list()
som.quality <- list()
edf <- list()

# i<-sample(1:nrow(parameters), 1)

for(i in 1:nrow(parameters)) {

ta <- Sys.time()

system.time({
ls.par <- 
    as.list(parameters[i,]) |>
    lapply(unlist)
ls <- do.call(generate_landscape_4cov_nl, ls.par)
})

# plot_landscape_4cov_lin(ls)

sam <- sample(1:nrow(ls), round(sam.frac * nrow(ls)))
ls.sam <- ls[sam,]
# sam <- 1:10000
# som.dim <- 25
# som.dim <- 10
grid <- somgrid(xdim = som.dim, ydim = som.dim, 
                topo = "rectangular", 
                neighbourhood.fct = "gaussian")
cov.z <- scale(ls.sam[, .(z1, z2, z3, z4)],
               center = TRUE, scale = TRUE)
som.fit <- som(cov.z,
               grid = grid, 
               rlen = som.rlen,
               radius = som.dim,
               init = init_som(cov.z, som.dim, som.dim),
               mode = "pbatch", 
               cores = n.threads,
               normalizeDataLayers = FALSE)
som.fit$scale <- list(mean = attr(cov.z, "scaled:center"),
                      sd = attr(cov.z, "scaled:scale"))

mapped <- 
    ls.sam[, .(z1, z2, z3, z4)] |>
    # land[sam, .(z1, z2, z3, z4)] |>
    scale_data_som(som = som.fit) |>
    embed_som(som = som.fit,
              grid.coord = TRUE)
ls.sam[,
         `:=`(som_bmu = mapped$bmu[,1],
              som_x = mapped$grid.coordinates$bmu.1[,"x"],
              som_y = mapped$grid.coordinates$bmu.1[,"y"])
         ]

if(eval.som) {
  quality <-
    evaluate_embedding(ls.sam[, .(z1, z2, z3, z4)],
                       mapped = as.matrix(ls.sam[, .(som_x, som_y)]),
                       k.max = 400,
                       combined = FALSE)

  som.quality[[i]] <-
    data.table(id = parameters[i, id],
               rec = mean(quality$dev.expl^-1)^-1,
               ve = variance_explained(som.fit),
               te = topological_error(som.fit))

  print(som.quality[[i]])
}

if(is.null(k.max)) {
  k.max <- min(nrow(unique(ls.sam[, .(som_x, som_y)]) - 1),
               floor((nrow(ls.sam) * (0.25))))
}

# mod.egp <- bam(response ~
#            # s(x, y, by = type, bs = "gp", k = 200) + s(som_x, som_y, bs = "gp", k = 200),
#            # s(x, y, by = type, bs = "gp", k = k.max) + s(som_x, som_y, bs = "gp", k = k.max),
#            type + s(x, y, bs = "gp", k = k.max) + s(som_x, som_y, bs = "gp", k = k.max),
#            data = ls.sam,
#            select = TRUE
#           )
mod.egp <- gam(response ~
           # s(x, y, by = type, bs = "gp", k = 200) + s(som_x, som_y, bs = "gp", k = 200),
           s(x, y, by = type, bs = egp.basis, k = k.max) + s(som_x, som_y, bs = egp.basis, k = k.max),
           # type + s(som_x, som_y, bs = "gp", k = k.max),
           data = ls.sam,
           select = egp.select,
           method= "REML",
           optimizer = "efs"
          )
# summary(mod.egp)
# AIC(mod.egp)

egp.edf <- NA
labels <- NA
n.smooth <- length(mod.egp$smooth)
for (j in 1:n.smooth) egp.edf[j] <- sum(mod.egp$edf[mod.egp$smooth[[j]]$first.para:mod.egp$smooth[[j]]$last.para])

# edf[n.smooth+1] <- k.max
for (j in 1:n.smooth) labels[j] <- mod.egp$smooth[[j]]$label
names(egp.edf) <- labels

edf[[i]] <- as.data.table(as.list(c(egp.edf, "k" = k.max, "id" = parameters[i, id])))

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
                    agg.fun = E,
                    ids = id.list.type,
                    )

groups.som <-
    ls.sam |>
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

mod.lm <- lm(response ~ type, data = ls.sam)
mod.lmcov <- lm(response ~ type + z1 + z2 + z3 + z4, data = ls.sam)


matched.cem <- matchit(type ~ z1 + z2 + z3 + z4, method = "cem", data = ls.sam)
md.cem <- match.data(matched.cem)
mod.cem <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.cem)

matched.nn.ps <-
  matchit(type ~ z1 + z2 + z3 + z4, method = "nearest", distance = "glm", data = ls.sam)
md.nn.ps <- match.data(matched.nn.ps)
mod.nn.ps <- lm(response ~ type + z1 + z2 + z3 + z4, weights = weights, data = md.nn.ps)

results[[i]] <-
  data.table(
             id = parameters[i, id],
             lm = coef(mod.lm)["typetreatment"],
             lmcov = coef(mod.lmcov)["typetreatment"],
             egp = mean(extract_variable(eff.arc, "treatment")),
             cem = coef(mod.cem)["typetreatment"],
             nn = coef(mod.nn.ps)["typetreatment"]
             )

print(results[[i]])

tb <- Sys.time()
te <- tb-ta
print(te)

}

results <- rbindlist(results)
edf <- rbindlist(edf)

export <- list(results = results, edf = edf)

saveRDS(export, file.results)
