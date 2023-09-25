library(data.table)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.rds")
file.egp.match.perm <- paste0(path.comp, "egp.match.perm.rds")

estimates <- readRDS(file.estimates)


## Permutation p-values for differences between EGP and matching methods


sub.dt <-
  CJ(
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     area.type = "treatment",
     mod.name = c("egp_som25", "match"),
     egp.som.topology = c(NA, "rectangular"),
     egp.cf.nb = c(NA, "sequential"),
     egp.geo.w = c(NA, TRUE),
     match.mod.cov = factor(c(NA, "interact")),
     match.trt.int = c(NA, FALSE),
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)

n.mc <- 1e4
chunk.size <- 1e3

egp.match.global <-
  compare_permutation(estimates.sub,
                      n.mc = n.mc,
                      chunk.size = chunk.size)


egp.match.ls <-
  compare_permutation(estimates.sub,
                      by.landscape = c("ls.response", "ls.imbalance"),
                      n.mc = n.mc,
                      chunk.size = chunk.size)

ls.response.lev <- c("all", levels(estimates.sub$ls.response))
ls.imbalance.lev <- c("all", levels(estimates.sub$ls.imbalance))

egp.match.global[, 
                 `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                      ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
egp.match.ls[, 
             `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                  ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]

egp.match.ls2 <- egp.match.ls
egp.match.global2 <- egp.match.global

setcolorder(egp.match, c("ls.response", "ls.imbalance", "name.short"))
setorder(egp.match, ls.response, ls.imbalance, name.short)

saveRDS(egp.match, file.egp.match.perm)

