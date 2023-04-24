args <- commandArgs(trailingOnly = TRUE)

library(mgcv)

source("utilities.R")

ls.type <- args[1]
mod.type <- args[2]
n.threads <- as.integer(args[3])
task_id <- as.integer(args[4])
task_count <- as.integer(args[5])

# ls.type <- "imbalance_medium"
# mod.type <- "egp_som10"
# n.threads <- 4
# task_id <- 1
# task_count <- 200

sam.frac <- 0.01
som.dim <- 10
som.epochs <- 1000
egp.k.som <- 100
egp.k.geo <- 250
egp.max.knots.som <- som.dim^2
egp.max.knots.geo <- egp.k.geo*10
egp.approx <- TRUE
overwrite <- TRUE

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
path.mod <- paste0(path.base, "models/", ls.type, "/")
if(!dir.exists(path.mod)) dir.create(path.mod, recursive = TRUE)

file.par <- paste0(path.ls, "parameters.rds")
file.log <- paste0(path.mod, mod.type, ".log")

parameters <- readRDS(file.par)
ls.total <- nrow(parameters)

row.chunks <- chunk_seq(1, nrow(parameters), ceiling(nrow(parameters) / task_count))
chunk <- row.chunks$from[task_id]:row.chunks$to[task_id]

if(!overwrite & file.exists(file.log)) {
  simulated <- as.integer(read.table(file.log)[,1])
} else {
  simulated <- integer(0)
}

chunk <- setdiff(chunk, simulated)

files.tmp <- paste0(tempdir(), "/", mod.type, "_",
                    parameters[id %in% chunk, file.name],
                    ".rds")

files.res <- paste0(path.mod, mod.type, "_",
                    parameters[id %in% chunk, file.name],
                    ".rds")

# chunk <- 5
for(i in chunk) {

  ta <- Sys.time()

  i.step <- which(chunk == i)
  
  message(paste0("Fitting EGP models for landscape ", i,
                 "/", ls.total, " (",
                 i.step, "/", length(chunk),
                 " in chunk)", " …"))

  results.ls <- list()

  ls.par <- 
    as.list(parameters[id == i,]) |>
    lapply(unlist)

  file.ls <- ls.par$file.path

  ls <- readRDS(file.ls)$landscape

  set.seed(ls.par$seed) 
  sam <- sample(1:nrow(ls), round(sam.frac * nrow(ls)))
  ls.sam <- na.omit(ls[sam,])

  mod.para <-
    CJ(landscape = ls.par$id,
       sam.frac,
       type = "egp",
       som.dim,
       som.epochs,
       som.topology = c("rectangular", "hexagonal"),
       cf.nb.strategy = c("sequential", "expand"),
       egp.k.som,
       egp.k.geo,
       egp.max.knots.som,
       egp.max.knots.geo,
       egp.approx)


  results.mod <- list()

  for(j in 1:nrow(mod.para)) {

    message(paste0("Parameter combination ", j, "/", nrow(mod.para)), " …")

    ls.fit <- copy(ls.sam)

    message("Fitting SOM …")

    som.egp <-
      egp_som(ls.fit,
              topo = mod.para$som.topology[j],
              x.dim = mod.para$som.dim[j],
              y.dim = mod.para$som.dim[j],
              epochs = mod.para$som.epochs[j],
              vars = c("z1", "z2", "z3", "z4"),
              parallel = n.threads)

    ls.fit[,
           c("som.unit", "som.x", "som.y") :=
             get_bmu(som.egp, coord = TRUE, list = TRUE)]


    message("Fitting GAM …")

    if(egp.approx) {
       mod.egp <- bam(response ~
                      s(x, y, by = type, bs = "gp", k = egp.k.geo,
                        xt = list(max.knots = egp.max.knots.geo)) +
                      s(som.x, som.y, bs = "gp", k = egp.k.som,
                        xt = list(max.knots = egp.max.knots.som)),
                      data = ls.fit,
                      select = TRUE,
                      discrete = TRUE,
                      nthreads = n.threads
                      )
     } else {
       mod.egp <- gam(response ~
                      s(x, y, by = type, bs = egp.basis, k = egp.k.geo,
                        xt = list(max.knots = egp.max.knots.geo)) +
                      s(som.x, som.y, bs = egp.basis, k = egp.k.som,
                        xt = list(max.knots = egp.max.knots.som)),
                      data = ls.fit,
                      select = TRUE,
                      method= "REML",
                      optimizer = "efs"
                     )
     }

    # summary(mod.egp)
    # AIC(mod.egp)


    message("Evaluating marginal effect …")

    egp.post <-
      egp_posterior_draw(mod.egp, 1000, unconditional = TRUE, package = "mgcv")

    egp.pred <-
      egp_posterior_predict(
                            model = mod.egp,
                            data = ls.fit,
                            id.var = "cell",
                            posterior = egp.post)

    egp.def <-
      egp_define_counterfactual(
                                data = ls.fit,
                                som = som.egp,
                                fac.ids = ls.fit[type == "treatment", cell],
                                cf.ids = ls.fit[type == "reference", cell],
                                nb.strategy = mod.para$cf.nb.strategy[j],
                                group.by = list(NULL, "poly"),
                                som.var = "som.unit",
                                id.var = "cell",
                                deg.max = NULL,
                                n.min = 1)

    egp.units <- egp_summarize_units(egp.pred, egp.def)
    egp.fac <- egp_evaluate_factual(egp.pred, egp.def, agg.size = 1e2)
    egp.count <- egp_evaluate_counterfactual(egp.pred, egp.def, egp.units)
    egp.mar <- egp_marginal(egp.fac, egp.count)

    results.mod[[j]] <-
      list(som = som.egp,
           model = mod.egp,
           posterior = egp.post,
           egp.def = egp.def,
           marginal = egp.mar)

  rm(ls.fit,
     som.egp, mod.egp,
     egp.post, egp.pred,
     egp.def,
     egp.units, egp.fac, egp.count,
     egp.mar)

  }

  mars <- numeric(0)
  for(k in seq_along(results.ls$models)) {
    mars[k] <-
      results.ls$models[[k]]$marginal[group.id == 1, mean(marginal)]
  }
  
  results.ls <-
    list(parameters = mod.para,
         models = results.mod)


  saveRDS(object = results.ls, file = files.tmp[i.step])

  rm(results.mod) 
  
  tb <- Sys.time()
  te <- tb-ta
  print(te)

}

message("Copying results to final destination …")

for(i in seq_along(files.tmp)) {
  file.copy(files.tmp[i], files.res[i], overwrite = TRUE)
}

paste0('echo "',
       paste0(chunk, collapse = "\n"),
       '" >> ',
       file.log) |>
system()
