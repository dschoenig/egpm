library(data.table)
library(progressr)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
path.dec <- paste0(path.results, "decision/")
file.estimates <- paste0(path.comp, "estimates.rds")
file.egp.subareas.boot <- paste0(path.comp, "egp.subareas.boot.rds")
file.egp.subareas.dec.all <- paste0(path.dec, "egp.subareas.dec.all.rds")
file.egp.subareas.dec.imb <- paste0(path.dec, "egp.subareas.dec.imb.rds")
file.egp.subareas.dec.ls <- paste0(path.dec, "egp.subareas.dec.ls.rds")
file.noeff.egp.subareas.dec.all <- paste0(path.dec, "noeff.egp.subareas.dec.all.rds")
file.noeff.egp.subareas.dec.imb <- paste0(path.dec, "noeff.egp.subareas.dec.imb.rds")
file.noeff.egp.subareas.dec.ls <- paste0(path.dec, "noeff.egp.subareas.dec.ls.rds")
main.seed <- 15930708 # Artemisia Gentileschi 

estimates <- readRDS(file.estimates)


## Subareas

sub.dt <-
  CJ(trt.effect = TRUE,
     sam.frac = c(0.01),
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     area.type = "subarea",
     mod.name = c("egp_som25"),
     egp.som.topology = c("rectangular"),
     egp.cf.nb = c("sequential"),
     egp.geo.w = c(TRUE),
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)
estimates.sub <- estimates.sub[treatment.prop >= 0.01]

# estimates.sub[order(mar.std), .(ls.response, ls.imbalance, mar.std, mar.est, subarea.id)]


## Performance measures


n.boot <- 1e4
chunk.size <- 125
handlers(global = TRUE)
options(future.globals.maxSize= 1024^3)


set.seed(main.seed+1)
egp.subareas.global <-
  compare_performance_boot(estimates.sub,
                           comparisons = NULL,
                           n.boot = n.boot,
                           chunk.size = chunk.size)

set.seed(main.seed+2)
egp.subareas.resp <-
  compare_performance_boot(estimates.sub,
                           by.landscape = "ls.response",
                           comparisons = NULL,
                           n.boot = n.boot,
                           chunk.size = chunk.size)

set.seed(main.seed+3)
egp.subareas.imb <-
  compare_performance_boot(estimates.sub,
                           by.landscape = "ls.imbalance",
                           comparisons = NULL,
                           n.boot = n.boot,
                           chunk.size = chunk.size)

set.seed(main.seed+4)
egp.subareas.resp.imb <-
  compare_performance_boot(estimates.sub,
                           by.landscape = c("ls.response", "ls.imbalance"),
                           comparisons = NULL,
                           n.boot = n.boot,
                           chunk.size = chunk.size)

estimates.sub[, .N, by = .(ls.response, ls.imbalance)]


ls.response.lev <- c("all", levels(estimates.sub$ls.response))
ls.imbalance.lev <- c("all", levels(estimates.sub$ls.imbalance))

egp.subareas.global[, 
                 `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                      ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
egp.subareas.resp[, 
               `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                    ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
egp.subareas.imb[, 
              `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                   ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]
egp.subareas.resp.imb[, 
                   `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                        ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]

egp.subareas <-
  do.call(rbind, list(egp.subareas.global, egp.subareas.resp, egp.subareas.imb, egp.subareas.resp.imb))

setcolorder(egp.subareas, c("ls.response", "ls.imbalance", "name.short"))
setorder(egp.subareas, ls.response, ls.imbalance, name.short)

saveRDS(egp.subareas, file.egp.subareas.boot)


## Decision support

sub.dt <-
  CJ(trt.effect = c(TRUE, FALSE),
     sam.frac = c(0.01),
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     area.type = "subarea",
     mod.name = c("egp_som25"),
     egp.som.topology = c("rectangular"),
     egp.cf.nb = c("sequential"),
     egp.geo.w = c(TRUE),
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)
estimates.sub <- estimates.sub[treatment.prop >= 0.01]

subareas <- unique(estimates.sub[order(ls.uid, subarea.id), .(ls.uid, subarea.id)])
subareas[, subarea.uid := 1:.N]
estimates.sub <- merge(estimates.sub, subareas, by = c("ls.uid", "subarea.id"))

# est.ci <- c("mar.q5", "mar.q95")
est.ci <- c("mar.q2.5", "mar.q97.5")
# est.ci <- c("mar.q0.5", "mar.q99.5")

egp.subareas.dec.all <-
  classify_estimates(estimates.sub[trt.effect == TRUE],
                     ls.id = "subarea.uid",
                     est.ci = est.ci,
                     group.vars = "trt.effect")

saveRDS(egp.subareas.dec.all, file.egp.subareas.dec.all)

egp.subareas.dec.imb <-
  classify_estimates(estimates.sub[trt.effect == TRUE],
                     ls.id = "subarea.uid",
                     est.ci = est.ci,
                     group.vars = "ls.imbalance")

saveRDS(egp.subareas.dec.imb, file.egp.subareas.dec.imb)

egp.subareas.dec.ls <-
  classify_estimates(estimates.sub[trt.effect == TRUE],
                     ls.id = "subarea.uid",
                     est.ci = est.ci,
                     group.vars = c("ls.response", "ls.imbalance"))

saveRDS(egp.subareas.dec.ls, file.egp.subareas.dec.ls)


noeff.egp.subareas.dec.all <-
  classify_estimates(estimates.sub[trt.effect == FALSE],
                     type = "noeffect",
                     noeff.cut = 0.5,
                     ls.id = "subarea.uid",
                     est.ci = est.ci,
                     est.std = "response.sd",
                     group.vars = NULL)

saveRDS(noeff.egp.subareas.dec.all, file.noeff.egp.subareas.dec.all)


noeff.egp.subareas.dec.imb <-
  classify_estimates(estimates.sub[trt.effect == FALSE],
                     type = "noeffect",
                     noeff.cut = 0.5,
                     ls.id = "subarea.uid",
                     est.ci = est.ci,
                     est.std = "response.sd",
                     group.vars = "ls.imbalance")


saveRDS(noeff.egp.subareas.dec.imb, file.noeff.egp.subareas.dec.imb)


noeff.egp.subareas.dec.ls <-
  classify_estimates(estimates.sub[trt.effect == FALSE],
                     type = "noeffect",
                     noeff.cut = 0.5,
                     ls.id = "subarea.uid",
                     est.ci = est.ci,
                     est.std = "response.sd",
                     group.vars = c("ls.response", "ls.imbalance"))

saveRDS(noeff.egp.subareas.dec.ls, file.noeff.egp.subareas.dec.ls)


egp.subareas.dec.ls$ls[, .N, by = .(ls.response, ls.imbalance)]
noeff.egp.subareas.dec.ls$ls[, .N, by = .(ls.response, ls.imbalance)]


