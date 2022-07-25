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

RFoptions(install="no")

source("~/projects/fcne_analysis/src/utilities.R")
source("utilities.R")

path.base <- "../"
ls.type <- "100_4cov_lin"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
path.ls.plots <- paste0(path.ls, "plots/")

if(!dir.exists(path.ls.data)) dir.create(path.ls.data, recursive = TRUE)
if(!dir.exists(path.ls.plots)) dir.create(path.ls.plots, recursive = TRUE)

n = 25
seed <- 19010511 # Rose Ausländer

set.seed(seed)
parameters <-
  data.table(
             id = 1:n,
             seed = sample(1:1e8, n),
             x.dim = 100,
             y.dim = 100,
             # Effect sizes
             treatment.effect.size.sp = 1,
             treatment.effect.size.bd = 0,
             z1.effect.size = round(runif(n, 1, 2) * sample(c(1, -1), n, replace = TRUE), 1),
             z2.effect.size = round(runif(n, 1, 2) * sample(c(1, -1), n, replace = TRUE), 1),
             z3.effect.size = round(runif(n, 1, 2) * sample(c(1, -1), n, replace = TRUE), 1),
             z4.effect.size = round(runif(n, 1, 2) * sample(c(1, -1), n, replace = TRUE), 1),
             # Parameters for generating functions
             treatment.nuclei = sample(2:7, n, replace = TRUE),
             treatment.phi.range = list(c(-pi/4, pi/4)),
             treatment.shift.range = list(c(-0.25, 0.25)),
             treatment.x.scale.range = list(c(0.5, 1)),
             treatment.y.scale.range = list(c(0.5, 1)),
             treatment.nuc.eff.range =  list(c(1,2)),
             z1.fbm.alpha = runif(n, 0.5, 1.5),
             z1.fbm.var = runif(n, 0.1, 1),
             z1.fbm.scale = runif(n, 0.1, 1),
             z1.fbm.w = 0.9,
             z1.grad.phi = runif(n, 0, 2*pi),
             z1.grad.w = 0.1,
             z2.fbm.alpha = runif(n, 0.5, 1.5),
             z2.fbm.var = runif(n, 0.1, 1),
             z2.fbm.scale = runif(n, 0.1, 1),
             z2.fbm.w = 0.1,
             z2.grad.phi = runif(n, 0, 2*pi),
             z2.grad.shift = runif(n, 0, 0.5),
             z2.grad.w = 0.9,
             z3.dist.n = sample(5:25, n, replace = TRUE),
             z3.dist.w = 0.9,
             z3.grad.phi = runif(n, 0, 2*pi),
             z3.grad.w = 0.1,
             z3.acc = 0.1,
             z4.seg.n = sample(10:50, n, replace = TRUE),
             z4.seg.nu = runif(n, 1, 2),
             z4.seg.var = runif(n, 0.1, 1),
             z4.seg.scale = runif(n, 0.1, 1),
             z4.seg.w = 0.9,
             z4.grad.phi = runif(n, 0, 2*pi),
             z4.grad.w = 0.1,
             split.n = 200,
             split.prop = 0.5,
             e.exp.var = 0.1,
             e.exp.scale = 100,
             e.nug.var = 0.1,
             e.rand.var = 0.1
             )
# parameters[, `:=`(z1.grad.phi = z2.grad.phi,
#                   z3.grad.phi = z2.grad.phi,
#                   z4.grad.phi = z2.grad.phi)]
# parameters

paste0(path.ls, "parameters.rds") |>
saveRDS(parameters, file = _)

# ls.par <- 
#   as.list(parameters[sample(1:.N, 1),]) |>
#   lapply(unlist)

system.time({
for(i in 1:nrow(parameters)) {
  ls.par <- 
    as.list(parameters[i,]) |>
    lapply(unlist)
  file.ls <- paste0(path.ls.data,
                    stri_pad_left(ls.par$id, 4, 0), ".rds")
  file.ls.plot <- paste0(path.ls.plots,
                         stri_pad_left(ls.par$id, 4, 0), ".tif")
  ls <- do.call(generate_landscape_4cov_lin, ls.par)
  ls.plot <- plot_landscape_4cov_lin(ls)
  saveRDS(ls, file.ls)
  # tiff("test", width = 6.7, height = 12, unit = "in", res = 300)
  tiff(file.ls.plot, width = 5, height = 9, unit = "in", res = 150)
  print(ls.plot)
  dev.off()
  }
})


