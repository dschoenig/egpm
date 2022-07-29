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

RFoptions(install="no")
RFoptions(spConform=FALSE)

source("utilities_fcne.R")
source("utilities.R")

task_id <- as.integer(args[1])
task_count <- as.integer(args[2])

# task_id <- 1
# task_count <- 100

path.base <- "../"
ls.type <- "1000_4cov_nl"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
if(!dir.exists(path.ls.data)) dir.create(path.ls.data, recursive = TRUE)

file.par <- paste0(path.ls, "parameters.rds")

n <- 1000
ls.dim <- 1000

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
             z3.acc = 1,
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

if(!file.exists(file.par) | task_id == 1) saveRDS(parameters, file.par)

row.chunks <- chunk_seq(1, nrow(parameters), ceiling(nrow(parameters) / task_count))
chunk <- row.chunks$from[task_id]:row.chunks$to[task_id]

for(i in chunk) {

  message(paste0("Generating landscape ", i, " / ", nrow(parameters),
                 " (chunk size ", row.chunks$size[i], ") …"))

  ta <- Sys.time()

  ls.par <- 
    as.list(parameters[i,]) |>
    lapply(unlist)
  file.ls <- paste0(path.ls.data,
                    stri_pad_left(ls.par$id, 4, 0), ".rds")
  ls <- do.call(generate_landscape_4cov_nl, ls.par)
  saveRDS(ls, file.ls)

  rm(ls)

  tb <- Sys.time()
  te <- tb-ta
  print(te)

}


