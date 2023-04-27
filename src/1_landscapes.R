args <- commandArgs(trailingOnly = TRUE)

# SETUP #############################################################

source("utilities.R")

overwrite <- TRUE

ls.name <- args[1]
n <- as.integer(args[2])
ls.dim <- as.integer(args[3])
ls.imbalance <- as.numeric(args[4])
task.id <- as.integer(args[5])
task.count <- as.integer(args[6])
parallel <- as.integer(args[7])
if(is.na(parallel)) parallel <- FALSE

# ls.name <- "imbalance_high"
# n <- 1000
# ls.dim <- 100
# ls.imbalance <- 0.3
# task.id <- 1
# task.count <- 1000
# parallel <- 4

paste0("Settings:\n",
       "N ", n, "\n",
       "Dimension ", ls.dim, "\n",
       "Imbalance ", ls.imbalance, "\n",
       "Parallel ", parallel) |>
message()

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.name, "/")
path.ls.data <- paste0(path.ls, "data/")
if(!dir.exists(path.ls.data)) dir.create(path.ls.data, recursive = TRUE)

file.par <- paste0(path.ls, "parameters.rds")
file.log <- paste0(path.ls, "simulate.log")


# PARAMETERS ########################################################

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
             areas.seg.seed = 3,
             areas.score.sam = 1e5,
             areas.imbalance.tol = 0,
             areas.area.prop = runif(n, 0.35, 0.65),
             areas.area.tol = list(c(0.35, 0.65)),
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
             areas.opt.pop = 100,
             areas.opt.prec = 1e-3,
             areas.opt.pcrossover = 0.9,
             areas.opt.pmutation = 0.5,
             areas.opt.elitism = 5,
             areas.opt.max.iter = 250,
             areas.opt.run = 100,
             areas.opt.parallel = parallel,
             areas.opt.fine = TRUE,
             areas.opt.fine.max.iter = 100,
             areas.opt.fine.constr = TRUE,
             areas.opt.fine.tol = 1e-6,
             areas.opt.cache = TRUE
             )

parameters[, file.name := stri_pad_left(id, ceiling(log10(n))+1, 0)]
parameters[, file.path := paste0(path.ls.data, file.name, ".rds")]


if(!file.exists(file.par) | task.id == 1) saveRDS(parameters, file.par)

if(!overwrite & file.exists(file.log)) {
  simulated <- as.integer(read.table(file.log)[,1])
} else {
  simulated <- integer(0)
}

row.chunks <- chunk_seq(1,
                        nrow(parameters),
                        ceiling(nrow(parameters) / task.count))
chunk <- row.chunks$from[task.id]:row.chunks$to[task.id]

chunk <- setdiff(chunk, simulated)

for(i in chunk) {

  message(paste0("Generating landscape ", i, " / ", nrow(parameters), " …"))

  ta <- Sys.time()

  ls.par <- 
    as.list(parameters[id == i,]) |>
    lapply(unlist, recursive = FALSE)
  ls.par$areas.opt.imb.agg <- ls.par$areas.opt.imb.agg[[1]]
  ls.par$verbose <- 2
  file.ls <- paste0(path.ls.data,
                    stri_pad_left(ls.par$id, ceiling(log10(n))+1, 0),
                    ".rds")

  ls <- NULL
  while(is.null(ls)) {
    try({ls <- do.call(generate_landscape_4cov_nl, ls.par)})
    if(is.null(ls)) message("Simulation failed. Trying again …")
  }

  message(paste0("Saving to ", file.ls, " …"))
  saveRDS(ls, file.ls)

  rm(ls)

  tb <- Sys.time()
  te <- tb-ta
  print(te)

  system(paste0('echo "', i, '" >> ', file.log), intern = TRUE)

}

