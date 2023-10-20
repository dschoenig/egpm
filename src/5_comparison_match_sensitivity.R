library(data.table)
library(progressr)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.rds")
file.match.covpar.boot <- paste0(path.comp, "match.covpar.boot.rds")
file.match.sam.boot <- paste0(path.comp, "match.sam.boot.rds")
main.seed <- 15930708 # Artemisia Gentileschi 

estimates <- readRDS(file.estimates)


n.boot <- 1e4
# n.boot <- 1e1
chunk.size <- 250
handlers(global = TRUE)
options(future.globals.maxSize= 10*1024^3)



## Different matching configurations

sub.dt <-
  rbind(
    CJ(sam.frac = c(0.01),
       ls.response = c("normal", "tweedie", "binary"),
       ls.imbalance = c("low", "high"),
       area.type = "treatment",
       # mod.name = "match",
       method = c("cem", "nearest"), 
       match.mod.cov = c("interact", "include", "exclude"),
       match.trt.int = FALSE,
       sorted = FALSE),
    CJ(sam.frac = c(0.01),
       ls.response = c("normal", "tweedie", "binary"),
       ls.imbalance = c("low", "high"),
       area.type = "treatment",
       # mod.name = "match",
       method = c("cem", "nearest"), 
       match.mod.cov = "interact",
       match.trt.int = TRUE,
       sorted = FALSE))

estimates.sub <- subset_estimates(estimates, sub.dt)

estimates.sub[,
              name.match := 
                fcase(
                      match.mod.cov == "interact" & match.trt.int == FALSE,
                      "default",
                      match.mod.cov == "exclude" & match.trt.int == FALSE,
                      "cov_exclude",
                      match.mod.cov == "include" & match.trt.int == FALSE,
                      "cov_include",
                      match.mod.cov == "interact" & match.trt.int == TRUE,
                      "trt_int")]

estimates.sub <- estimates.sub[!is.na(name.match)]
estimates.sub[,
              name.match :=
                factor(name.match,
                       levels = c("default", "cov_exclude", "cov_include", "trt_int"))]

match.methods <- as.character(unique(estimates.sub$name.short))

match.covpar.l <- list()

for(i in seq_along(match.methods)) {

  message(paste0("Processing method ", i, "/", length(match.methods), " …"))

  estimates.sub.method <- estimates.sub[name.short == match.methods[i]]

  set.seed(main.seed+1)
  match.covpar.global <-
    compare_boot(estimates.sub.method,
                 by.method = "name.match",
                 by.landscape = "name.short",
                 comparisons = "default",
                 pe.type = "data",
                 n.boot = n.boot,
                 chunk.size = chunk.size)


  set.seed(main.seed+2)
  match.covpar.resp <-
    compare_boot(estimates.sub.method,
                 by.method = "name.match",
                 comparisons = "default",
                 by.landscape = c("name.short", "ls.response"),
                 n.boot = n.boot,
                 chunk.size = chunk.size)

  set.seed(main.seed+3)
  match.covpar.imb <-
    compare_boot(estimates.sub.method,
                 by.method = "name.match",
                 comparisons = "default",
                 by.landscape = c("name.short", "ls.imbalance"),
                 n.boot = n.boot,
                 chunk.size = chunk.size)

  set.seed(main.seed+4)
  match.covpar.resp.imb <-
    compare_boot(estimates.sub.method,
                 by.method = "name.match",
                 comparisons = "default",
                 by.landscape = c("name.short", "ls.response", "ls.imbalance"),
                 n.boot = n.boot,
                 chunk.size = chunk.size)


  ls.response.lev <- c("all", levels(estimates.sub$ls.response))
  ls.imbalance.lev <- c("all", levels(estimates.sub$ls.imbalance))

  match.covpar.global[, 
                   `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                        ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
  match.covpar.resp[, 
                 `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                      ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
  match.covpar.imb[, 
                `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                     ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]
  match.covpar.resp.imb[, 
                     `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                          ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]

  match.covpar.l[[i]] <- 
    do.call(rbind,
            list(match.covpar.global, match.covpar.resp, match.covpar.imb, match.covpar.resp.imb))


}

names(match.covpar.l) <- match.methods
match.covpar <- rbindlist(match.covpar.l)
# match.covpar[, name.short := factor(name.short, levels = levels(estimates$name.short))]

setcolorder(match.covpar, c("ls.response", "ls.imbalance", "name.short", "name.match"))
setorder(match.covpar, ls.response, ls.imbalance, name.short, name.match)

saveRDS(match.covpar, file.match.covpar.boot)



## Different sample sizes

sub.dt <-
  CJ(sam.frac = c(0.01, 0.005, 0.02),
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     area.type = "treatment",
     method = c("cem", "nearest"), 
     match.mod.cov = "include",
     match.trt.int = FALSE,
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)


estimates.sub[,
              name.match := 
                fcase(sam.frac == 0.01,
                      "default",
                      sam.frac == 0.005,
                      "sam_005",
                      sam.frac == 0.02,
                      "sam_02")]

estimates.sub <- estimates.sub[!is.na(name.match)]
estimates.sub[,
              name.match :=
                factor(name.match,
                       levels = c("default", "sam_005", "sam_02"))]


match.methods <- as.character(unique(estimates.sub$name.short))

match.sam.l <- list()

for(i in seq_along(match.methods)) {

  message(paste0("Processing method ", i, "/", length(match.methods), " …"))

  estimates.sub.method <- estimates.sub[name.short == match.methods[i]]

  set.seed(main.seed+1)
  match.sam.global <-
    compare_boot(estimates.sub.method,
                 by.method = "name.match",
                 by.landscape = "name.short",
                 comparisons = "default",
                 pe.type = "data",
                 n.boot = n.boot,
                 chunk.size = chunk.size)


  set.seed(main.seed+2)
  match.sam.resp <-
    compare_boot(estimates.sub.method,
                 by.method = "name.match",
                 comparisons = "default",
                 by.landscape = c("name.short", "ls.response"), 
                 n.boot = n.boot,
                 chunk.size = chunk.size)

  set.seed(main.seed+3)
  match.sam.imb <-
    compare_boot(estimates.sub.method,
                 by.method = "name.match",
                 comparisons = "default",
                 by.landscape = c("name.short", "ls.imbalance"), 
                 n.boot = n.boot,
                 chunk.size = chunk.size)

  set.seed(main.seed+4)
  match.sam.resp.imb <-
    compare_boot(estimates.sub.method,
                 by.method = "name.match",
                 comparisons = "default",
                 by.landscape = c("name.short", "ls.response", "ls.imbalance"),
                 n.boot = n.boot,
                 chunk.size = chunk.size)


  ls.response.lev <- c("all", levels(estimates.sub$ls.response))
  ls.imbalance.lev <- c("all", levels(estimates.sub$ls.imbalance))

  match.sam.global[, 
                   `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                        ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
  match.sam.resp[, 
                 `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                      ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
  match.sam.imb[, 
                `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                     ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]
  match.sam.resp.imb[, 
                     `:=`(ls.response = factor(ls.response, levels = ls.response.lev),
                          ls.imbalance = factor(ls.imbalance, levels = ls.imbalance.lev))]

  match.sam.l[[i]] <-
    do.call(rbind, list(match.sam.global, match.sam.resp, match.sam.imb, match.sam.resp.imb))
}

match.sam <- rbindlist(match.sam.l)

setcolorder(match.sam, c("ls.response", "ls.imbalance", "name.short", "name.match"))
setorder(match.sam, ls.response, ls.imbalance, name.short, name.match)

saveRDS(match.sam, file.match.sam.boot)

