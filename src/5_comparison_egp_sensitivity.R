library(data.table)
library(progressr)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.rds")
file.egp.sompar.boot <- paste0(path.comp, "egp.sompar.boot.rds")
file.egp.sam.boot <- paste0(path.comp, "egp.sam.boot.rds")
main.seed <- 15930708 # Artemisia Gentileschi 

estimates <- readRDS(file.estimates)


n.boot <- 1e4
# n.boot <- 1e1
chunk.size <- 125
handlers(global = TRUE)
options(future.globals.maxSize= 3*1024^3)



## Different SOM configurations

sub.dt <-
  CJ(sam.frac = c(0.01),
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     area.type = "treatment",
     mod.name = c("egp_som25", "egp_som10", "egp_som50", "egp_som25_unequal"),
     egp.som.topology = c("rectangular", "hexagonal"),
     egp.som.dim = c(10, 25, 50),
     egp.cf.nb = c(NA, "sequential"),
     egp.geo.w = c(NA, TRUE),
     match.mod.cov = NA,
     match.trt.int = NA,
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)


estimates.sub[,
              name.egp := 
                fcase(mod.name == "egp_som25" & egp.som.topology == "rectangular",
                      "default",
                      mod.name == "egp_som25" & egp.som.topology == "hexagonal",
                      "hex",
                      mod.name == "egp_som10" & egp.som.topology == "rectangular",
                      "dim_10",
                      mod.name == "egp_som50" & egp.som.topology == "rectangular",
                      "dim_50",
                      mod.name == "egp_som25_unequal" & egp.som.topology == "rectangular",
                      "unequal")]

estimates.sub <- estimates.sub[!is.na(name.egp) & trt.effect == TRUE]
estimates.sub[,
              name.egp :=
                factor(name.egp,
                       levels = c("default", "hex", "unequal", "dim_10", "dim_50"))]


set.seed(main.seed+1)
egp.sompar.global <-
  compare_performance_boot(copy(estimates.sub),
                           by.method = "name.egp",
                           comparisons = "default",
                           pe.type = "data",
                           n.boot = n.boot,
                           chunk.size = chunk.size)


set.seed(main.seed+2)
egp.sompar.resp <-
  compare_performance_boot(estimates.sub,
                           by.method = "name.egp",
                           comparisons = "default",
                           by.landscape = "ls.response",
                           n.boot = n.boot,
                           chunk.size = chunk.size)

set.seed(main.seed+3)
egp.sompar.imb <-
  compare_performance_boot(estimates.sub,
                           by.method = "name.egp",
                           comparisons = "default",
                           by.landscape = "ls.imbalance",
                           n.boot = n.boot,
                           chunk.size = chunk.size)

set.seed(main.seed+4)
egp.sompar.resp.imb <-
  compare_performance_boot(estimates.sub,
                           by.method = "name.egp",
                           comparisons = "default",
                           by.landscape = c("ls.response", "ls.imbalance"),
                           n.boot = n.boot,
                           chunk.size = chunk.size)


ls.response.lev <- c("all", levels(estimates.sub$ls.response))
ls.imbalance.lev <- c("all", levels(estimates.sub$ls.imbalance))

egp.sompar.global[, 
                 `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                      ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
egp.sompar.resp[, 
               `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                    ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
egp.sompar.imb[, 
              `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                   ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]
egp.sompar.resp.imb[, 
                   `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                        ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]

egp.sompar <-
  do.call(rbind, list(egp.sompar.global, egp.sompar.resp, egp.sompar.imb, egp.sompar.resp.imb))

setcolorder(egp.sompar, c("ls.response", "ls.imbalance", "name.egp"))
setorder(egp.sompar, ls.response, ls.imbalance, name.egp)

saveRDS(egp.sompar, file.egp.sompar.boot)



## Different sample sizes

# sub.dt <-
#   CJ(sam.frac = c(0.01, 0.005, 0.02),
#      ls.response = c("normal", "tweedie", "binary"),
#      ls.imbalance = c("low", "high"),
#      area.type = "treatment",
#      egp.som.topology = c("rectangular"),
#      egp.som.dim = c(25),
#      egp.som.unequal = FALSE,
#      egp.cf.nb = c(NA, "sequential"),
#      egp.geo.w = c(NA, TRUE),
#      match.mod.cov = NA,
#      match.trt.int = NA,
#      sorted = FALSE)

# estimates.sub <- subset_estimates(estimates, sub.dt)


# estimates.sub[,
#               name.egp := 
#                 fcase(sam.frac == 0.01,
#                       "default",
#                       sam.frac == 0.005,
#                       "sam_005",
#                       sam.frac == 0.02,
#                       "sam_02")]

# estimates.sub <- estimates.sub[!is.na(name.egp)]
# estimates.sub[,
#               name.egp :=
#                 factor(name.egp,
#                        levels = c("default", "sam_005", "sam_02"))]


# set.seed(main.seed+1)
# egp.sam.global <-
#   compare_performance_boot(copy(estimates.sub),
#                by.method = "name.egp",
#                comparisons = "default",
#                pe.type = "data",
#                n.boot = n.boot,
#                chunk.size = chunk.size)


# set.seed(main.seed+2)
# egp.sam.resp <-
#   compare_performance_boot(estimates.sub,
#                by.method = "name.egp",
#                comparisons = "default",
#                by.landscape = "ls.response",
#                n.boot = n.boot,
#                chunk.size = chunk.size)

# set.seed(main.seed+3)
# egp.sam.imb <-
#   compare_performance_boot(estimates.sub,
#                by.method = "name.egp",
#                comparisons = "default",
#                by.landscape = "ls.imbalance",
#                n.boot = n.boot,
#                chunk.size = chunk.size)

# set.seed(main.seed+4)
# egp.sam.resp.imb <-
#   compare_performance_boot(estimates.sub,
#                by.method = "name.egp",
#                comparisons = "default",
#                by.landscape = c("ls.response", "ls.imbalance"),
#                n.boot = n.boot,
#                chunk.size = chunk.size)


# ls.response.lev <- c("all", levels(estimates.sub$ls.response))
# ls.imbalance.lev <- c("all", levels(estimates.sub$ls.imbalance))

# egp.sam.global[, 
#                  `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
#                       ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
# egp.sam.resp[, 
#                `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
#                     ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
# egp.sam.imb[, 
#               `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
#                    ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]
# egp.sam.resp.imb[, 
#                    `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
#                         ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]

# egp.sam <-
#   do.call(rbind, list(egp.sam.global, egp.sam.resp, egp.sam.imb, egp.sam.resp.imb))

# setcolorder(egp.sam, c("ls.response", "ls.imbalance", "name.egp"))
# setorder(egp.sam, ls.response, ls.imbalance, name.egp)

# saveRDS(egp.sam, file.egp.sam.boot)

