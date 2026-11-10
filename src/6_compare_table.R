library(data.table)
library(brms)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.csv")
file.comp.global <- paste0(path.comp, "comparison.global.csv")
file.comp.ls <- paste0(path.comp, "comparison.ls.csv")

estimates <- fread(file.estimates, yaml = TRUE)

sub.dt <-
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

estimates.sub <- 
  estimates[sub.dt,
            on = names(sub.dt),
            nomatch = NULL]

estimates.sub[,
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


num.prec <- 3

est.egp.global <-
  estimates.sub[name.short == "EGP",
                .(ls.response = "all",
                  ls.imbalance = "all",
                  bias.egp = round(mean(mar.std-1), num.prec),
                  rmse.egp = round(rmse(mar.std, 1), num.prec),
                  ser.egp = round(sum(sign(mar.est) != sign(mar.true)) / .N, num.prec))]

est.sum.global <-
  estimates.sub[,
                .(ls.response = "all",
                  ls.imbalance = "all",
                  bias = round(mean(mar.std-1), num.prec),
                  rmse = round(rmse(mar.std, 1), num.prec),
                  ser = round(sum(sign(mar.est) != sign(mar.true)) / .N, num.prec)),
                  by = .(name.short)] |>
  merge(est.egp.global,
        by = c("ls.response", "ls.imbalance")) |>
  _[,
    `:=`(bias.vs = abs(bias) - abs(bias.egp),
         rmse.vs = rmse - rmse.egp,
         ser.vs = ser - ser.egp)
    ][order(ls.response, ls.imbalance, name.short),
      -c("ls.response", "ls.imbalance",
         "bias.egp", "rmse.egp", "ser.egp")]


est.egp.ls <-
  estimates.sub[name.short == "EGP",
                .(bias.egp = round(mean(mar.std-1), num.prec),
                  rmse.egp = round(rmse(mar.std, 1), num.prec),
                  ser.egp = round(sum(sign(mar.est) != sign(mar.true)) / .N, num.prec)),
                  by = .(ls.response, ls.imbalance)]
est.sum.ls <-
  estimates.sub[,
                .(bias = round(mean(mar.std-1), num.prec),
                  rmse = round(rmse(mar.std, 1), num.prec),
                  ser = round(sum(sign(mar.est) != sign(mar.true)) / .N, num.prec)),
                  by = .(ls.response, ls.imbalance, name.short)] |>
  merge(est.egp.ls,
        by = c("ls.response", "ls.imbalance")) |>
  _[,
    `:=`(bias.vs = abs(bias) - abs(bias.egp),
         rmse.vs = rmse - rmse.egp,
         ser.vs = ser - ser.egp)
    ][order(ls.response, ls.imbalance, name.short),
      -c("bias.egp", "rmse.egp", "ser.egp")]

fwrite(est.sum.global, file.comp.global, yaml = TRUE)
fwrite(est.sum.ls, file.comp.ls, yaml = TRUE)
