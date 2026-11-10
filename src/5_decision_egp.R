library(data.table)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
path.dec <- paste0(path.results, "decision/")
file.estimates <- paste0(path.comp, "estimates.rds")


file.egp.dec.all <- paste0(path.dec, "egp.dec.all.rds")
file.egp.dec.imb <- paste0(path.dec, "egp.dec.imb.rds")
file.egp.dec.ls <- paste0(path.dec, "egp.dec.ls.rds")

file.noeff.egp.dec.all <- paste0(path.dec, "noeff.egp.dec.all.rds")
file.noeff.egp.dec.imb <- paste0(path.dec, "noeff.egp.dec.imb.rds")
file.noeff.egp.dec.ls <- paste0(path.dec, "noeff.egp.dec.ls.rds")

estimates <- readRDS(file.estimates)

sub.dt <-
  CJ(trt.effect = c(TRUE, FALSE),
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     area.type = "treatment",
     mod.name = c("egp_som25"),
     egp.som.topology = c("rectangular"),
     egp.cf.nb = c("sequential"),
     egp.geo.w = c(TRUE),
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)

# est.ci <- c("mar.q5", "mar.q95")
est.ci <- c("mar.q2.5", "mar.q97.5")
# est.ci <- c("mar.q0.5", "mar.q99.5")

estimates.sub[trt.effect == FALSE, .N, by = ls.name]

egp.dec.all <-
  classify_estimates(estimates.sub[trt.effect == TRUE],
                     ls.id = "ls.uid",
                     est.ci = est.ci,
                     group.vars = "trt.effect")

saveRDS(egp.dec.all, file.egp.dec.all)

egp.dec.imb <-
  classify_estimates(estimates.sub[trt.effect == TRUE],
                     ls.id = "ls.uid",
                     est.ci = est.ci,
                     group.vars = "ls.imbalance")

saveRDS(egp.dec.imb, file.egp.dec.imb)

egp.dec.ls <-
  classify_estimates(estimates.sub[trt.effect == TRUE],
                     ls.id = "ls.id",
                     est.ci = est.ci,
                     group.vars = c("ls.response", "ls.imbalance"))

saveRDS(egp.dec.ls, file.egp.dec.ls)


noeff.egp.dec.all <-
  classify_estimates(estimates.sub[trt.effect == FALSE],
                     type = "noeffect",
                     noeff.cut = 0.5,
                     ls.id = "ls.uid",
                     est.ci = est.ci,
                     est.std = "response.sd",
                     group.vars = NULL)

saveRDS(noeff.egp.dec.all, file.noeff.egp.dec.all)


noeff.egp.dec.imb <-
  classify_estimates(estimates.sub[trt.effect == FALSE],
                     type = "noeffect",
                     noeff.cut = 0.5,
                     ls.id = "ls.uid",
                     est.ci = est.ci,
                     est.std = "response.sd",
                     group.vars = "ls.imbalance")


saveRDS(noeff.egp.dec.imb, file.noeff.egp.dec.imb)


noeff.egp.dec.ls <-
  classify_estimates(estimates.sub[trt.effect == FALSE],
                     type = "noeffect",
                     noeff.cut = 0.5,
                     ls.id = "ls.id",
                     est.ci = est.ci,
                     est.std = "response.sd",
                     group.vars = c("ls.response", "ls.imbalance"))

saveRDS(noeff.egp.dec.ls, file.noeff.egp.dec.ls)


# noeff.egp.dec.ls$ls[, .N, by = .(ls.imbalance, .type)]
