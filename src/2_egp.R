args <- commandArgs(trailingOnly = TRUE)

library(mgcv)

source("utilities.R")

ls.type <- args[1]
mod.type <- args[2]
resp.type <- args[3]
sam.frac <- as.numeric(args[4])
som.dim <- as.numeric(args[5])
som.epochs <- as.integer(args[6])
egp.k.som <- as.integer(args[7])
egp.k.geo <- as.integer(args[8])
egp.approx <- as.logical(args[9])
n.threads <- as.integer(args[10])
task_id <- as.integer(args[11])
task_count <- as.integer(args[12])
overwrite <- as.logical(args[13])
if(is.na(overwrite)) overwrite <- TRUE

# ls.type <- "binary_imbalance_high"
# mod.type <- "egp_som25"
# resp.type <- "binary"
# n.threads <- 4
# task_id <- 407
# task_count <- 1000
# sam.frac <- 0.01
# som.dim <- 25
# som.epochs <- 1000
# egp.k.som <- 250
# egp.k.geo <- 250
# egp.approx <- TRUE
# overwrite <- TRUE

egp.max.knots.som <- som.dim^2
egp.max.knots.geo <- egp.k.geo*10

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


setDTthreads(n.threads)

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
       geo = c(FALSE, TRUE),
       egp.k.som,
       egp.k.geo,
       egp.max.knots.som,
       egp.max.knots.geo,
       egp.approx)

  results.mod <- list()

  for(j in 1:nrow(mod.para)) {

    message(paste0("Parameter combination ", j, "/", nrow(mod.para)), " …")

    ls.fit <- copy(ls.sam)

    if(j > 1) {
      som.same <-
        mod.para[c(j-1, j),
                 all(unlist(lapply(.SD, \(x) x[1] == x[2]))),
                 .SDcols = c("som.topology", "som.dim", "som.epochs"),
                 nomatch = NULL]
    } else {
      som.same <- FALSE
    }

    if(som.same) {
      message("Reusing previous SOM …")
    } else {
      message("Fitting SOM …")
      som.egp <-
        egp_som(ls.fit,
                topo = mod.para$som.topology[j],
                x.dim = mod.para$som.dim[j],
                y.dim = mod.para$som.dim[j],
                epochs = mod.para$som.epochs[j],
                vars = c("z1", "z2", "z3", "z4"),
                parallel = n.threads)
    }

    ls.fit[,
           c("som.unit", "som.x", "som.y") :=
             get_bmu(som.egp, coord = TRUE, list = TRUE)]

    ls.fit[, type := factor(type, levels = levels(type), ordered = TRUE)]

    mod.fam <-
      switch(resp.type,
             normal = gaussian(link = "identity"),
             binary = binomial(link = "logit"),
             beta = betar(link = "logit", eps = .Machine$double.eps * 1e4),
             tweedie = tw(link = "log"))

    if(j > 1) {
      para.same <-
        mod.para[c(j-1, j),
                 all(unlist(lapply(.SD, \(x) x[1] == x[2]))),
                 .SDcols = c("egp.k.geo", "egp.max.knots.geo",
                             "egp.k.som", "egp.max.knots.som"),
                 nomatch = NULL]
    } else {
      para.same <- FALSE
    }

    mod.same <- all(para.same, som.same)

    if(mod.same) {
      message("Reusing previous GAM and predictions …")
    } else {

      message("Fitting GAM …")

      if(j > 1) {
        rm(mod.egp, egp.post, egp.pred)
        gc()
      }

      if(egp.approx) {

        if(resp.type != "beta") {
          mod.egp <- bam(response ~
                         s(x, y, bs = "gp", k = egp.k.geo,
                           xt = list(max.knots = egp.max.knots.geo)) +
                         s(x, y, by = type, bs = "gp", k = egp.k.geo,
                           xt = list(max.knots = egp.max.knots.geo)) +
                         s(som.x, som.y, bs = "gp", k = egp.k.som,
                           xt = list(max.knots = egp.max.knots.som)),
                         family = mod.fam,
                         data = ls.fit,
                         select = TRUE,
                         discrete = TRUE,
                         nthreads = n.threads
          )
        } else {
          # Small model for better initial smoothing parameter guess
          small.sam <- sample(ls.fit$cell, 1e3)
          mod.egp.small <- gam(response ~
                               s(x, y, bs = "gp", k = round(egp.k.geo/10),
                                 xt = list(max.knots = round(egp.max.knots.geo/10))) +
                               s(x, y, by = type, bs = "gp", k = round(egp.k.geo/10),
                                 xt = list(max.knots = round(egp.max.knots.geo/10))) +
                               s(som.x, som.y, bs = "gp", k = round(egp.k.som/10),
                                 xt = list(max.knots = round(egp.max.knots.som/10))),
                               family = mod.fam,
                               optimizer = "efs",
                               method = "REML",
                               data = ls.fit[cell %in% small.sam],
                               select = TRUE
                               )
          mod.egp <- bam(response ~
                         s(x, y, bs = "gp", k = egp.k.geo,
                           xt = list(max.knots = egp.max.knots.geo)) +
                         s(x, y, by = type, bs = "gp", k = egp.k.geo,
                           xt = list(max.knots = egp.max.knots.geo)) +
                         s(som.x, som.y, bs = "gp", k = egp.k.som,
                           xt = list(max.knots = egp.max.knots.som)),
                         family = mod.fam,
                         data = ls.fit,
                         sp = mod.egp.small$sp,
                         select = TRUE,
                         discrete = TRUE,
                         nthreads = n.threads
                         )
        }

      } else {
        mod.egp <- gam(response ~
                       s(x, y, bs = "gp", k = egp.k.geo,
                         xt = list(max.knots = egp.max.knots.geo)) +
                       s(x, y, by = type, bs = "gp", k = egp.k.geo,
                         xt = list(max.knots = egp.max.knots.geo)) +
                       s(som.x, som.y, bs = "gp", k = egp.k.som,
                         xt = list(max.knots = egp.max.knots.som)),
                       family = mod.fam,
                       data = ls.fit,
                       select = TRUE,
                       method= "REML",
                       optimizer = "efs"
        )
      }


      message("Evaluating posterior predictive distributions …")

      egp.post <-
        egp_posterior_draw(mod.egp, 1000, unconditional = TRUE, package = "mgcv")

      predict.chunk <-
        switch(resp.type,
               normal = NULL,
               binary = NULL,
               beta = NULL,
               tweedie = 2500)

      egp.pred <-
        egp_posterior_predict(model = mod.egp,
                              data = ls.fit,
                              id.var = "cell",
                              posterior = egp.post,
                              predict.chunk = predict.chunk,
                              epred = FALSE)

    }

    # summary(mod.egp)
    # AIC(mod.egp)
    # library(DHARMa)
    # mod.res <- simulateResiduals(mod.egp)
    # plot(mod.res)

    message("Evaluating marginal effect …")

    if(mod.para$geo[j]) {
      mod.geo.vars <- c("x", "y")
    } else {
      mod.geo.vars <- NULL
    }

    egp.def <-
      egp_define_counterfactual(data = ls.fit,
                                som = som.egp,
                                fac.ids = ls.fit[type == "treatment", cell],
                                cf.ids = ls.fit[type == "reference", cell],
                                nb.strategy = mod.para$cf.nb.strategy[j],
                                group.by = list(NULL, "poly"),
                                som.var = "som.unit",
                                geo.vars = mod.geo.vars,
                                geo.kernel = "matern32",
                                id.var = "cell",
                                deg.max = NULL,
                                n.min = 1)

    egp.fac <- egp_evaluate_factual(egp.pred, egp.def, agg.size = 5e3)
    egp.count <- egp_evaluate_counterfactual(egp.pred, egp.def, agg.size = 5e3)
    egp.mar <- egp_marginal(egp.fac, egp.count)

    print(egp.mar[group.id == 1, mean(marginal)])
    print(egp.mar[group.id == 1, quantile(marginal, c(0.025, 0.975))])

    results.mod[[j]] <-
      list(som = som.egp,
           model = mod.egp,
           posterior = egp.post,
           egp.def = egp.def,
           marginal = egp.mar)

    rm(ls.fit,
       egp.def,
       egp.fac, egp.count,
       egp.mar)
    gc()

  }

  # mars <- numeric(0)
  # for(k in seq_along(results.mod)) {
  #   mars[k] <-
  #     results.mod[[k]]$marginal[group.id == 1, mean(marginal)]
  # }

  results.ls <-
    list(parameters = mod.para,
         sample = ls.sam,
         models = results.mod)


  saveRDS(object = results.ls, file = files.tmp[i.step], compress = FALSE)

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
