args <- commandArgs(trailingOnly = TRUE)

library(mgcv)

source("utilities.R")

ls.type <- args[1]
task_id <- as.integer(args[2])
task_count <- as.integer(args[3])
overwrite <- as.logical(args[4])
if(is.na(overwrite)) overwrite <- TRUE

ls.type <- "imbalance_high"
task_id <- 1
task_count <- 1000

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")

file.par <- paste0(path.ls, "parameters.rds")
file.log <- paste0(path.ls, "mar.log")

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

# chunk <- 5
for(i in chunk) {

  ta <- Sys.time()

  i.step <- which(chunk == i)
  
  message(paste0("Fixing marginal for landscape ", i,
                 "/", ls.total, " (",
                 i.step, "/", length(chunk),
                 " in chunk)", " …"))

  ls.par <- 
    as.list(parameters[id == i,]) |>
    lapply(unlist)
  file.ls <- ls.par$file.path

  ls <- readRDS(file.ls)

  ls <- ls$landscape

  cov.effects <- names(ls)[names(ls) %like% "f."]
  cov.n <- length(cov.effects)


  cov.mat <-
    ls[,
       # c(cov.effects, "res.sp"),
       c(cov.effects, "res.sp"),
       with = FALSE] |>
    as.matrix()
  median(rowSums(cov.mat))


  a <- quantile(rowSums(cov.mat), 1-p.mean, names = FALSE)


  ref.form <-
    parse(text = paste(c("intercept", cov.effects, "res.sp"),
                       collapse = " + "))
  trt.form <-
    parse(text = paste(c("intercept", "treatment", cov.effects, "res.sp"),
                       collapse = " + "))

  mar.trt <-
    ls.bin[type == "treatment",
           .(p.ref = sum(eval(ref.form) > 0) / .N,
             p.trt = sum(eval(trt.form) > 0) / .N)]

  n.poly.trt <- ls.bin[type == "treatment", length(unique(poly))]

  if(n.poly.trt > 1) {
  mar.trt <-
    rbind(mar.trt,
          ls.bin[type == "treatment",
                 .(p.ref = sum(eval(ref.form) > 0) / .N,
                   p.trt = sum(eval(trt.form) > 0) / .N),
                 by = poly][order(poly)],
          fill = TRUE)
  setcolorder(mar.trt, c("poly", "p.ref", "p.trt"))
  }









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
      message("Reusing previous GAM …")
    } else {
      message("Fitting GAM …")
      if(egp.approx) {
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
    }

    # summary(mod.egp)
    # AIC(mod.egp)
    # library(DHARMa)
    # mod.res <- simulateResiduals(mod.egp)
    # plot(mod.res)

    message("Evaluating marginal effect …")

    egp.post <-
      egp_posterior_draw(mod.egp, 1000, unconditional = TRUE, package = "mgcv")

    egp.pred <-
      egp_posterior_predict(model = mod.egp,
                            data = ls.fit,
                            id.var = "cell",
                            posterior = egp.post)

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

    egp.fac <- egp_evaluate_factual(egp.pred, egp.def)
    egp.count <- egp_evaluate_counterfactual(egp.pred, egp.def)
    egp.mar <- egp_marginal(egp.fac, egp.count)

    print(egp.mar[group.id == 1, mean(marginal)])

    results.mod[[j]] <-
      list(som = som.egp,
           model = mod.egp,
           posterior = egp.post,
           egp.def = egp.def,
           marginal = egp.mar)

    rm(ls.fit,
       egp.post, egp.pred,
       egp.def,
       egp.fac, egp.count,
       egp.mar)

  }

  # mars <- numeric(0)
  # for(k in seq_along(results.mod)) {
  #   mars[k] <-
  #     results.mod[[k]]$marginal[group.id == 1, mean(marginal)]
  # }

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
