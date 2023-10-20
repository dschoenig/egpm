library(data.table)
library(progressr)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.rds")
file.egp.match.boot <- paste0(path.comp, "egp.match.boot.rds")
main.seed <- 15930708 # Artemisia Gentileschi 

estimates <- readRDS(file.estimates)


## Entire treatment area

sub.dt <-
  CJ(sam.frac = c(0.01),
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

n.boot <- 1e4
chunk.size <- 250
handlers(global = TRUE)
options(future.globals.maxSize= 1024^3)

set.seed(main.seed+1)
egp.match.global <-
  compare_boot(estimates.sub,
               n.boot = n.boot,
               chunk.size = chunk.size)

set.seed(main.seed+2)
egp.match.resp <-
  compare_boot(estimates.sub,
               by.landscape = "ls.response",
               n.boot = n.boot,
               chunk.size = chunk.size)

set.seed(main.seed+3)
egp.match.imb <-
  compare_boot(estimates.sub,
               by.landscape = "ls.imbalance",
               n.boot = n.boot,
               chunk.size = chunk.size)

set.seed(main.seed+4)
egp.match.resp.imb <-
  compare_boot(estimates.sub,
               by.landscape = c("ls.response", "ls.imbalance"),
               n.boot = n.boot,
               chunk.size = chunk.size)


ls.response.lev <- c("all", levels(estimates.sub$ls.response))
ls.imbalance.lev <- c("all", levels(estimates.sub$ls.imbalance))

egp.match.global[, 
                 `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                      ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
egp.match.resp[, 
               `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                    ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
egp.match.imb[, 
              `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                   ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]
egp.match.resp.imb[, 
                   `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                        ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]

egp.match <-
  do.call(rbind, list(egp.match.global, egp.match.resp, egp.match.imb, egp.match.resp.imb))

setcolorder(egp.match, c("ls.response", "ls.imbalance", "name.short"))
setorder(egp.match, ls.response, ls.imbalance, name.short)

saveRDS(egp.match, file.egp.match.boot)

