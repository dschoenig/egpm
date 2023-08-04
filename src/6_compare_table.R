library(data.table)
library(brms)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.csv")

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


estimates.sub[
              # ls.imbalance == "low"
              # ls.response == "tweedie" & ls.imbalance == "low"
              ,
              .(bias = mean(mar.std-1),
                cv = sd(mar.std)/mean(mar.std),
                ser = sum(sign(mar.est) != sign(mar.true)) / .N,
                rmse = rmse(mar.std, 1)),
                by = .(name.short, mod.id,
                       method,
                       match.distance,
                       match.cutpoints,
                       match.mod.cov,
                       match.trt.int,
                       match.replace
                       # )][order(-bias)]
                       )][order(name.short)]

estimates.sub[,
              `:=`(ls.response = factor(ls.response,
                                        levels = c("normal", "tweedie", "binary")),
                   ls.imbalance = factor(ls.imbalance,
                                         levels = c("low", "high")),
                   ls.id = factor(ls.id))]

fwrite(estimates.sub, file.estfit)



