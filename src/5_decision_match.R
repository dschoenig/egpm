library(data.table)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
path.dec <- paste0(path.results, "decision/")
file.estimates <- paste0(path.comp, "estimates.rds")


file.match.dec.all <- paste0(path.dec, "match.dec.all.rds")
file.match.dec.imb <- paste0(path.dec, "match.dec.imb.rds")
file.match.dec.ls <- paste0(path.dec, "match.dec.ls.rds")

file.noeff.match.dec.all <- paste0(path.dec, "noeff.match.dec.all.rds")
file.noeff.match.dec.imb <- paste0(path.dec, "noeff.match.dec.imb.rds")
file.noeff.match.dec.ls <- paste0(path.dec, "noeff.match.dec.ls.rds")

estimates <- readRDS(file.estimates)

sub.dt <-
  CJ(trt.effect = c(TRUE, FALSE),
     sam.frac = c(0.01),
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     area.type = "treatment",
     mod.name = "match",
     match.mod.cov = "interact",
     match.trt.int = FALSE,
     sorted = FALSE)


estimates.sub <- subset_estimates(estimates, sub.dt)

# est.ci <- c("mar.q5", "mar.q95")
est.ci <- c("mar.q2.5", "mar.q97.5")
# est.ci <- c("mar.q0.5", "mar.q99.5")

estimates.sub[trt.effect == FALSE, .N, by = ls.name]

match.dec.all <-
  classify_estimates(estimates.sub[trt.effect == TRUE],
                     ls.id = "ls.uid",
                     est.ci = est.ci,
                     est.range = c(-3, 3),
                     group.vars = "name.short")

saveRDS(match.dec.all, file.match.dec.all)


match.dec.imb <-
  classify_estimates(estimates.sub[trt.effect == TRUE],
                     ls.id = "ls.uid",
                     est.ci = est.ci,
                     est.range = c(-3, 3),
                     group.vars = c("name.short", "ls.imbalance"))

saveRDS(match.dec.imb, file.match.dec.imb)

match.dec.ls <-
  classify_estimates(estimates.sub[trt.effect == TRUE],
                     ls.id = "ls.id",
                     est.ci = est.ci,
                     est.range = c(-3, 3),
                     group.vars = c("name.short", "ls.response", "ls.imbalance"))

saveRDS(match.dec.ls, file.match.dec.ls)


noeff.match.dec.all <-
  classify_estimates(estimates.sub[trt.effect == FALSE],
                     type = "noeffect",
                     noeff.cut = 0.5,
                     ls.id = "ls.uid",
                     est.ci = est.ci,
                     est.std = "response.sd",
                     est.range = c(-3, 3),
                     group.vars = c("name.short"))

saveRDS(noeff.match.dec.all, file.noeff.match.dec.all)


noeff.match.dec.imb <-
  classify_estimates(estimates.sub[trt.effect == FALSE],
                     type = "noeffect",
                     noeff.cut = 0.5,
                     ls.id = "ls.uid",
                     est.ci = est.ci,
                     est.std = "response.sd",
                     est.range = c(-3, 3),
                     group.vars = c("name.short", "ls.imbalance"))


saveRDS(noeff.match.dec.imb, file.noeff.match.dec.imb)


noeff.match.dec.ls <-
  classify_estimates(estimates.sub[trt.effect == FALSE],
                     type = "noeffect",
                     noeff.cut = 0.5,
                     ls.id = "ls.id",
                     est.ci = est.ci,
                     est.std = "response.sd",
                     est.range = c(-3, 3),
                     group.vars = c("name.short", "ls.response", "ls.imbalance"))

saveRDS(noeff.match.dec.ls, file.noeff.match.dec.ls)


# noeff.match.dec.all$ls[, .N, by = .(ls.imbalance, .type)]
