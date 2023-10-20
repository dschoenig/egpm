library(data.table)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
path.dec <- paste0(path.results, "decision/")
file.estimates <- paste0(path.comp, "estimates.rds")


file.cemst.dec_95.all <- paste0(path.dec, "cemst.dec_95.all.rds")
file.cemst.dec_95.ls <- paste0(path.dec, "cemst.dec_95.ls.rds")
file.cemst.dec_95.all.sum <- paste0(path.dec, "cemst.dec_95.all.sum.csv")
file.cemst.dec_95.ls.sum <- paste0(path.dec, "cemst.dec_95.ls.sum.csv")


estimates <- readRDS(file.estimates)

sub.dt <-
  CJ(
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     area.type = "treatment",
     mod.name = c("match"),
     egp.som.topology = c(NA),
     egp.cf.nb = c(NA),
     egp.geo.w = c(NA),
     match.cutpoints = "st",
     match.mod.cov = factor(c("interact")),
     match.trt.int = c(FALSE),
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)

cemst.dec_95.all <-
  classify_estimates(estimates.sub,
                     ls.id = "ls.uid",
                     est.ci = c("mar.q2.5", "mar.q97.5"),
                     group.vars = "ls.imbalance")

saveRDS(cemst.dec_95.all, file.cemst.dec_95.all)

cemst.dec_95.all.sum <-
  cemst.dec_95.all$ls[order(ls.imbalance, .type),
                    .(count = .N),
                    by = c("ls.imbalance", ".type")]

cemst.dec_95.all.sum[, freq := count/sum(count), by = "ls.imbalance"]


cemst.dec_95.ls <-
  classify_estimates(estimates.sub,
                     ls.id = "ls.id",
                     est.ci = c("mar.q2.5", "mar.q97.5"),
                     group.vars = c("ls.response", "ls.imbalance"))

saveRDS(cemst.dec_95.ls, file.cemst.dec_95.ls)

cemst.dec_95.ls.sum <-
  cemst.dec_95.ls$ls[order(ls.imbalance, .type),
                   .(count = .N),
                   by = c("ls.response", "ls.imbalance", ".type")]

cemst.dec_95.ls.sum[, freq := count/sum(count), by = c("ls.response", "ls.imbalance")]
