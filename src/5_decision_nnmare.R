library(data.table)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
path.dec <- paste0(path.results, "decision/")
file.estimates <- paste0(path.comp, "estimates.rds")


file.nnmare.dec_95.all <- paste0(path.dec, "nnmare.dec_95.all.rds")
file.nnmare.dec_95.ls <- paste0(path.dec, "nnmare.dec_95.ls.rds")
file.nnmare.dec_95.all.sum <- paste0(path.dec, "nnmare.dec_95.all.sum.csv")
file.nnmare.dec_95.ls.sum <- paste0(path.dec, "nnmare.dec_95.ls.sum.csv")


estimates <- readRDS(file.estimates)

sub.dt <-
  CJ(
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     area.type = "treatment",
     mod.name = c("match"),
     # egp.som.topology = c(NA),
     # egp.cf.nb = c(NA),
     # egp.geo.w = c(NA),
     match.replace = TRUE,
     match.distance = "mahalanobis",
     match.mod.cov = factor(c("interact")),
     match.trt.int = c(FALSE),
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)

nnmare.dec_95.all <-
  classify_estimates(estimates.sub,
                     ls.id = "ls.uid",
                     est.ci = c("mar.q2.5", "mar.q97.5"),
                     group.vars = "ls.imbalance")

saveRDS(nnmare.dec_95.all, file.nnmare.dec_95.all)

nnmare.dec_95.all.sum <-
  nnmare.dec_95.all$ls[order(ls.imbalance, .type),
                    .(count = .N),
                    by = c("ls.imbalance", ".type")]

nnmare.dec_95.all.sum[, freq := count/sum(count), by = "ls.imbalance"]


nnmare.dec_95.ls <-
  classify_estimates(estimates.sub,
                     ls.id = "ls.id",
                     est.ci = c("mar.q2.5", "mar.q97.5"),
                     group.vars = c("ls.response", "ls.imbalance"))

saveRDS(nnmare.dec_95.ls, file.nnmare.dec_95.ls)

nnmare.dec_95.ls.sum <-
  nnmare.dec_95.ls$ls[order(ls.imbalance, .type),
                   .(count = .N),
                   by = c("ls.response", "ls.imbalance", ".type")]

nnmare.dec_95.ls.sum[, freq := count/sum(count), by = c("ls.response", "ls.imbalance")]
