library(brms)
library(data.table)

source("utilities.R")

mod.type <- "normal_high"

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.csv")
file.model <- paste0(path.comp, "mod.brm.rds")

estimates <- fread(file.estimates, yaml = TRUE)

# mod <- readRDS(file.model)
mod <- readRDS("../results/comparison/mod.brm.28.rds")


path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.csv")

estimates <- fread(file.estimates, yaml = TRUE)

fit.dt <-
  CJ(
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("high", "low"),
     area.type = c("treatment"),
     mod.name = c("egp_som25", "match"),
     egp.som.topology = c(NA, "rectangular"),
     egp.cf.nb = c(NA, "sequential"),
     egp.geo.w = c(NA, TRUE),
     match.mod.cov = c(NA, "interact"),
     match.trt.int = c(NA, FALSE),
     sorted = FALSE)

estimates.fit <- 
  estimates[fit.dt,
            on = names(fit.dt),
            nomatch = NULL]

estimates.fit[,
              `:=`(ls.response = factor(ls.response,
                                        levels = c("normal", "tweedie", "binary")),
                   ls.imbalance = factor(ls.imbalance,
                                         levels = c("low", "high")),
                   ls.id = factor(ls.id),
                   name.short = factor(name.short,
                                       levels = c("EGP", "GLM",
                                                  "CEM-ST", "CEM-SC", "CEM-FD",
                                                  "NN-PS-NR", "NN-PS-RE",
                                                  "NN-MA-NR", "NN-MA-RE")))]

ls.type.lev <-
  with(estimates.fit,
       paste0(rep(levels(ls.response), each = 2), "_",
              rep(levels(ls.imbalance), times = 3)))

estimates.fit[,
              ls.type := factor(paste0(ls.response, "_", ls.imbalance),
                                levels = ls.type.lev)]
estimates.fit[,
              ls.uid := as.integer((as.integer(ls.type)-1) * 1000) + as.integer(as.character(ls.id))]


method.chunk <- estimates.fit[, unique(name.short)]
pred.l <- list()
for(i in seq_along(method.chunk)) {
  method.foc <- method.chunk[i]
  pred.mat <-
    t(posterior_predict(mod,
                        newdata = estimates.fit[name.short == method.foc],
                        draw_ids = draw.ids.pred))
  pred.method <- draws_dt(pred.mat, estimates.fit[name.short == method.foc])
  rm(pred.mat)
  pred.l[[i]] <-
    pred.method[,
                .(bias = bias(value, 1),
                  rmse = rmse(value, 1),
                  ser = ser(value)),
                by = c(".draw")] |>
    melt(id.vars = c(".draw"),
         measure.vars = c("bias", "rmse", "ser"),
         variable.name = "stat",
         value.name = "value")
  rm(pred.method)
  gc()
  message(i)
}

names(pred.l) <- method.chunk
pred.dt <- rbindlist(pred.l, idcol = "name.short")
pred.dt[, name.short := factor(name.short, levels = levels(estimates.fit$name.short))]
setorder(pred.dt, name.short, .draw)


est.dt <-
  estimates.fit[order(name.short),
                .(bias = bias(mar.std, 1),
                  rmse = rmse(mar.std, 1),
                  ser = ser(mar.std)),
                by = .(name.short)]


prepare_table_all <- function(estimates, predictions) {

  pred.dt <- predictions
  est.dt <- estimates

  pred.vs.egp <-
    dcast(pred.dt,
          .draw + stat ~ name.short, value.var = "value") |>
    _[,
      c(.(.draw = .draw),
        lapply(.SD, \(x) abs(x) - .SD[, abs(EGP)])),
      by = stat,
      .SDcols = levels(estimates.fit$name.short)] |>
    _[,
      lapply(.SD, \(x) sum(x > 0)/.N),
      by = stat,
      .SDcols = levels(estimates.fit$name.short)] |>
    melt(id.vars = "stat",
         variable.name = "name.short") |>
    dcast(name.short ~ stat,
          value.var = "value")
  setnames(pred.vs.egp,
           c("bias", "rmse", "ser"),
           paste0(c("bias", "rmse", "ser"), ".p.egp"))

  pred.vs.glm <-
    dcast(pred.dt,
          .draw + stat ~ name.short, value.var = "value") |>
    _[,
      c(.(.draw = .draw),
        lapply(.SD, \(x) abs(x) - .SD[, abs(GLM)])),
      by = stat,
      .SDcols = levels(estimates.fit$name.short)] |>
    _[,
      lapply(.SD, \(x) sum(x > 0)/.N),
      by = stat,
      .SDcols = levels(estimates.fit$name.short)] |>
    melt(id.vars = "stat",
         variable.name = "name.short") |>
    dcast(name.short ~ stat,
          value.var = "value")
  setnames(pred.vs.glm,
           c("bias", "rmse", "ser"),
           paste0(c("bias", "rmse", "ser"), ".p.glm"))


  est.vs.egp <-
    melt(est.dt,
         measure.vars = c("bias", "rmse", "ser"),
         variable.name = "stat",
         value.name = "value") |>
    dcast(stat ~ name.short, value.var = "value") |>
    _[,
      lapply(.SD, \(x) abs(x) - .SD[, abs(EGP)]),
      by = stat,
      .SDcols = levels(estimates.fit$name.short)] |>
    melt(id.vars = "stat",
         variable.name = "name.short") |>
    dcast(name.short ~ stat,
          value.var = "value") |>
    setnames(c("bias", "rmse", "ser"),
             paste0(c("bias", "rmse", "ser"), ".vs.egp"))

  est.vs.glm <-
    melt(est.dt,
         measure.vars = c("bias", "rmse", "ser"),
         variable.name = "stat",
         value.name = "value") |>
    dcast(stat ~ name.short, value.var = "value") |>
    _[,
      lapply(.SD, \(x) abs(x) - .SD[, abs(GLM)]),
      by = stat,
      .SDcols = levels(estimates.fit$name.short)] |>
    melt(id.vars = "stat",
         variable.name = "name.short") |>
    dcast(name.short ~ stat,
          value.var = "value") |>
    setnames(c("bias", "rmse", "ser"),
             paste0(c("bias", "rmse", "ser"), ".vs.glm"))

  est.dt[est.vs.egp,
         on = c("ls.response", "ls.imbalance", "name.short")
         ][est.vs.glm,
         on = c("ls.response", "ls.imbalance", "name.short")]


  est.comb <-
    cbind(est.dt[order(name.short)],
          est.vs.egp[order(name.short), -"name.short"],
          pred.vs.egp[order(name.short), -"name.short"],
          est.vs.glm[order(name.short), -"name.short"],
          pred.vs.glm[order(name.short), -"name.short"])

  cols.num <- c("bias", "rmse", "ser", grep(".vs.", names(est.comb), value = TRUE))
  cols.prob <- grep(".p.", names(est.comb), value = TRUE)

  est.comb.print <- copy(est.comb)
  est.comb.print[,
                 (cols) := lapply(.SD, num_format),
                 .SDcols = cols.num,
                 env = list(cols = I(cols.num))]
  est.comb.print[,
                 (cols) := lapply(.SD, prob_format),
                 .SDcols = cols.prob,
                 env = list(cols = I(cols.prob))]

  table.out <-
    list(combined = est.comb.print,
         estimates = est.dt[order(name.short)],
         vs.egp = est.vs.egp[order(name.short), -"name.short"],
         prob.egp = pred.vs.egp[order(name.short), -"name.short"],
         vs.glm = est.vs.glm[order(name.short), -"name.short"],
         prob.glm = pred.vs.glm[order(name.short), -"name.short"])

  return(table.out)
}



prepare_table_ls <- function(estimates, predictions) {

  pred.dt <- predictions
  est.dt <- estimates

  pred.vs.egp <-
    dcast(pred.dt,
          .draw + ls.response + ls.imbalance + stat ~ name.short,
          value.var = "value") |>
    _[,
      c(.(.draw = .draw),
        lapply(.SD, \(x) abs(x) - .SD[, abs(EGP)])),
      by = c("ls.response", "ls.imbalance", "stat"),
      .SDcols = levels(estimates.fit$name.short)] |>
    _[,
      lapply(.SD, \(x) sum(x > 0)/.N),
      by = c("ls.response", "ls.imbalance", "stat"),
      .SDcols = levels(estimates.fit$name.short)] |>
    melt(id.vars = c("ls.response", "ls.imbalance", "stat"),
         variable.name = "name.short") |>
    dcast(ls.response + ls.imbalance + name.short ~ stat,
          value.var = "value")
  setnames(pred.vs.egp,
           c("bias", "rmse", "ser"),
           paste0(c("bias", "rmse", "ser"), ".p.egp"))

  pred.vs.glm <-
    dcast(pred.dt,
          .draw + ls.response + ls.imbalance + stat ~ name.short,
          value.var = "value") |>
    _[,
      c(.(.draw = .draw),
        lapply(.SD, \(x) abs(x) - .SD[, abs(GLM)])),
      by = c("ls.response", "ls.imbalance", "stat"),
      .SDcols = levels(estimates.fit$name.short)] |>
    _[,
      lapply(.SD, \(x) sum(x > 0)/.N),
      by = c("ls.response", "ls.imbalance", "stat"),
      .SDcols = levels(estimates.fit$name.short)] |>
    melt(id.vars = c("ls.response", "ls.imbalance", "stat"),
         variable.name = "name.short") |>
    dcast(ls.response + ls.imbalance + name.short ~ stat,
          value.var = "value")
  setnames(pred.vs.glm,
           c("bias", "rmse", "ser"),
           paste0(c("bias", "rmse", "ser"), ".p.glm"))


  est.vs.egp <-
    melt(est.dt,
         id.vars = c("ls.response", "ls.imbalance", "name.short"),
         measure.vars = c("bias", "rmse", "ser"),
         variable.name = "stat",
         value.name = "value") |>
    dcast(ls.response + ls.imbalance + stat ~ name.short,
          value.var = "value") |>
    _[,
      lapply(.SD, \(x) abs(x) - .SD[, abs(EGP)]),
      by = c("ls.response", "ls.imbalance", "stat"),
      .SDcols = levels(estimates.fit$name.short)] |>
    melt(id.vars = c("ls.response", "ls.imbalance", "stat"),
         variable.name = "name.short") |>
    dcast(ls.response + ls.imbalance + name.short ~ stat,
          value.var = "value") |>
    setnames(c("bias", "rmse", "ser"),
             paste0(c("bias", "rmse", "ser"), ".vs.egp"))

  est.vs.glm <-
    melt(est.dt,
         id.vars = c("ls.response", "ls.imbalance", "name.short"),
         measure.vars = c("bias", "rmse", "ser"),
         variable.name = "stat",
         value.name = "value") |>
    dcast(ls.response + ls.imbalance + stat ~ name.short,
          value.var = "value") |>
    _[,
      lapply(.SD, \(x) abs(x) - .SD[, abs(GLM)]),
      by = c("ls.response", "ls.imbalance", "stat"),
      .SDcols = levels(estimates.fit$name.short)] |>
    melt(id.vars = c("ls.response", "ls.imbalance", "stat"),
         variable.name = "name.short") |>
    dcast(ls.response + ls.imbalance + name.short ~ stat,
          value.var = "value") |>
    setnames(c("bias", "rmse", "ser"),
             paste0(c("bias", "rmse", "ser"), ".vs.glm"))

  col.id <- c("name.short", "ls.response", "ls.imbalance")

  est.comb <-
    cbind(est.dt[order(ls.response, ls.imbalance, name.short)],
          est.vs.egp[order(ls.response, ls.imbalance, name.short),
                     -c("name.short", "ls.response", "ls.imbalance")],
          pred.vs.egp[order(ls.response, ls.imbalance, name.short),
                      -c("name.short", "ls.response", "ls.imbalance")],
          est.vs.glm[order(ls.response, ls.imbalance, name.short),
                     -c("name.short", "ls.response", "ls.imbalance")],
          pred.vs.glm[order(ls.response, ls.imbalance, name.short),
                      -c("name.short", "ls.response", "ls.imbalance")])

  cols.num <- c("bias", "rmse", "ser", grep(".vs.", names(est.comb), value = TRUE))
  cols.prob <- grep(".p.", names(est.comb), value = TRUE)

  est.comb.print <- copy(est.comb)
  est.comb.print[,
                 (cols) := lapply(.SD, num_format),
                 .SDcols = cols.num,
                 env = list(cols = I(cols.num))]
  est.comb.print[,
                 (cols) := lapply(.SD, prob_format),
                 .SDcols = cols.prob,
                 env = list(cols = I(cols.prob))]

  table.out <-
    list(combined = est.comb,
         combined.print,
         estimates = est.dt[order(name.short)],
         vs.egp = est.vs.egp[order(name.short), -"name.short"],
         prob.egp = pred.vs.egp[order(name.short), -"name.short"],
         vs.glm = est.vs.glm[order(name.short), -"name.short"],
         prob.glm = pred.vs.glm[order(name.short), -"name.short"])

  return(table.out)
}
















fit.dt <-
  CJ(
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("high", "low"),
     area.type = c("treatment"),
     mod.name = c("egp_som25", "match"),
     egp.som.topology = c(NA, "rectangular"),
     egp.cf.nb = c(NA, "sequential"),
     egp.geo.w = c(NA, TRUE),
     match.mod.cov = c(NA, "interact"),
     match.trt.int = c(NA, FALSE),
     sorted = FALSE)

estimates.fit <- 
  estimates[fit.dt,
            on = names(fit.dt),
            nomatch = NULL]

estimates.fit[,
              `:=`(ls.response = factor(ls.response,
                                        levels = c("normal", "tweedie", "binary")),
                   ls.imbalance = factor(ls.imbalance,
                                         levels = c("low", "high")),
                   ls.id = factor(ls.id),
                   name.short = factor(name.short,
                                       levels = c("EGP", "GLM",
                                                  "CEM-ST", "CEM-SC", "CEM-FD",
                                                  "NN-PS-NR", "NN-PS-RE",
                                                  "NN-MA-NR", "NN-MA-RE")))]

ls.type.lev <-
  with(estimates.fit,
       paste0(rep(levels(ls.response), each = 2), "_",
              rep(levels(ls.imbalance), times = 3)))

estimates.fit[,
              ls.type := factor(paste0(ls.response, "_", ls.imbalance),
                                levels = ls.type.lev)]
estimates.fit[,
              ls.uid := as.integer((as.integer(ls.type)-1) * 1000) + as.integer(as.character(ls.id))]


method.chunk <- estimates.fit[, unique(name.short)]
pred.global.l <- list()
pred.ls.l <- list()
for(i in seq_along(method.chunk)) {
  method.foc <- method.chunk[i]
  pred.mat <-
    t(posterior_predict(mod,
                        newdata = estimates.fit[name.short == method.foc]))
  pred.method <- draws_dt(pred.mat, estimates.fit[name.short == method.foc])
  rm(pred.mat)
  pred.global.l[[i]] <-
    pred.method[,
                .(bias = bias(value, 1),
                  rmse = rmse(value, 1),
                  ser = ser(value)),
                by = ".draw"] |>
    melt(id.vars = ".draw",
         measure.vars = c("bias", "rmse", "ser"),
         variable.name = "stat",
         value.name = "value")
  pred.ls.l[[i]] <-
    pred.method[,
                .(bias = bias(value, 1),
                  rmse = rmse(value, 1),
                  ser = ser(value)),
                by = c("ls.response", "ls.imbalance", ".draw")] |>
    melt(id.vars = c("ls.response", "ls.imbalance", ".draw"),
         measure.vars = c("bias", "rmse", "ser"),
         variable.name = "stat",
         value.name = "value")
  rm(pred.method)
  gc()
  message(i)
}

names(pred.global.l) <- method.chunk
pred.global.dt <- rbindlist(pred.global.l, idcol = "name.short")
pred.global.dt[, name.short := factor(name.short, levels = levels(estimates.fit$name.short))]
setorder(pred.global.dt, ls.response, ls.imbalance, name.short)

names(pred.ls.l) <- method.chunk
pred.ls.dt <- rbindlist(pred.ls.l, idcol = "name.short")
pred.ls.dt[, name.short := factor(name.short, levels = levels(estimates.fit$name.short))]
setorder(pred.ls.dt, ls.response, ls.imbalance, name.short)

est.global.dt <-
  estimates.fit[,
                .(bias = bias(mar.std, 1),
                  rmse = rmse(mar.std, 1),
                  ser = ser(mar.std)),
                by = "name.short"]
setorder(est.global.dt, name.short)

est.ls.dt <-
  estimates.fit[,
                .(bias = bias(mar.std, 1),
                  rmse = rmse(mar.std, 1),
                  ser = ser(mar.std)),
                by = c("ls.response", "ls.imbalance", "name.short")]
setorder(est.ls.dt, ls.response, ls.imbalance, name.short)


