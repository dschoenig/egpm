library(data.table)
library(progressr)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.rds")
file.egp.subarea.boot <- paste0(path.comp, "egp.subarea.boot.rds")
main.seed <- 15930708 # Artemisia Gentileschi 

estimates <- readRDS(file.estimates)


## Entire treatment area

sub.dt <-
  CJ(sam.frac = c(0.01),
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     area.type = "subarea",
     mod.name = c("egp_som25"),
     egp.som.topology = c(NA, "rectangular"),
     egp.cf.nb = c(NA, "sequential"),
     egp.geo.w = c(NA, TRUE),
     match.mod.cov = factor(c(NA, "interact")),
     match.trt.int = c(NA, FALSE),
     sorted = FALSE)


estimates.sub <- subset_estimates(estimates, sub.dt)

subareas <- unique(estimates.sub[order(ls.uid, subarea.id), .(ls.uid, subarea.id)])
subareas[, subarea.uid := 1:.N]
estimates.sub <- merge(estimates.sub, subareas, by = c("ls.uid", "subarea.id"))


n.boot <- 1e4
chunk.size <- 250
handlers(global = TRUE)
options(future.globals.maxSize= 2*1024^3)

set.seed(main.seed+1)
egp.subarea.global <-
  compare_boot(estimates.sub,
               by.method = "name.short",
               comparisons = NULL,
               ls.id = "subarea.uid",
               n.boot = n.boot,
               chunk.size = chunk.size)

set.seed(main.seed+2)
egp.subarea.resp <-
  compare_boot(estimates.sub,
               by.method = "name.short",
               comparisons = NULL,
               ls.id = "subarea.uid",
               by.landscape = "ls.response",
               n.boot = n.boot,
               chunk.size = chunk.size)

set.seed(main.seed+3)
egp.subarea.imb <-
  compare_boot(estimates.sub,
               by.method = "name.short",
               comparisons = NULL,
               ls.id = "subarea.uid",
               by.landscape = "ls.imbalance",
               n.boot = n.boot,
               chunk.size = chunk.size)

set.seed(main.seed+4)
egp.subarea.resp.imb <-
  compare_boot(estimates.sub,
               by.method = "name.short",
               comparisons = NULL,
               ls.id = "subarea.uid",
               by.landscape = c("ls.response", "ls.imbalance"),
               n.boot = n.boot,
               chunk.size = chunk.size)


ls.response.lev <- c("all", levels(estimates.sub$ls.response))
ls.imbalance.lev <- c("all", levels(estimates.sub$ls.imbalance))

egp.subarea.global[, 
                   `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                        ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
egp.subarea.resp[, 
                 `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                      ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
egp.subarea.imb[, 
                `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                     ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]
egp.subarea.resp.imb[, 
                     `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                          ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]

egp.subarea <-
  do.call(rbind, list(egp.subarea.global, egp.subarea.resp, egp.subarea.imb, egp.subarea.resp.imb))

setcolorder(egp.subarea, c("ls.response", "ls.imbalance", "name.short"))
setorder(egp.subarea, ls.response, ls.imbalance, name.short)

saveRDS(egp.subarea, file.egp.subarea.boot)

