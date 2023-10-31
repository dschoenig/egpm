library(data.table)

source("utilities.R")

args <- commandArgs(trailingOnly = TRUE)
ls.type <- args[1]
resp.type <- args[2]

# ls.type <- "noeff_imbalance_high"
# resp.type <- "normal"

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")

file.results <- paste0(path.ls, "summary.rds")

file.par <- paste0(path.ls, "parameters.rds")

parameters <- readRDS(file.par)

ids <- parameters$id

cov.l <- list()
optim.l <- list()
means.l <- list()
marginal.l <- list()
obs.l <- list()


for(i in ids) {

  if(file.exists(parameters[id == i, file.path])) {

    ta <- Sys.time()

    message(paste0("Summarizing results for landscape ", i, " …"))
    
    ls <- readRDS(parameters[id == i, file.path])

    cov.names <- stri_subset(names(ls$landscape), regex = "^z\\d+")

    setindex(ls$landscape, type)
    obs.trt <- ls$landscape[.("treatment"), ..cov.names, on = "type"]
    obs.ref <- ls$landscape[.("reference"), ..cov.names, on = "type"]

    cov.j <- list()
    for(j in seq_along(cov.names)) {
      obs.trt.cov <- obs.trt[[cov.names[j]]]
      obs.ref.cov <- obs.ref[[cov.names[j]]]
      cov.j[[j]] <-
        data.table(variable = cov.names[j],
                   d_cohen = d_cohen(obs.trt.cov, obs.ref.cov,
                                     na.rm = TRUE),
                   var_ratio = var_ratio(obs.trt.cov, obs.ref.cov,
                                         na.rm = TRUE),
                   ks_stat = ks_stat(obs.trt.cov, obs.ref.cov,
                                     na.rm = TRUE),
                   js_div = js_div(obs.trt.cov, obs.ref.cov,
                                   type = "discrete", na.rm = TRUE))
    }
    cov.l[[i]] <- rbindlist(cov.j)

    imb <- ls$optim[names(ls$optim) %like% "imbalance"]
    optim.l[[i]] <- 
      data.table(imbalance = list(imb),
                 imbalance.mean = mean(imb),
                 area.prop = ls$optim["area.prop"])

    means.ls <-
      rbind(
            ls$landscape[,
                         .(response.mean = mean(response),
                           response.sd = sd(response))],
            ls$landscape[,
                         .(response.mean = mean(response),
                           response.sd = sd(response)),
                         by = "type"],
            fill = TRUE)
    if(ls$landscape[type == "treatment", length(unique(poly))] > 1) {
      means.ls <-
      rbind(means.ls,
            ls$landscape[type == "treatment",
                         .(response.mean = mean(response),
                           response.sd = sd(response)),
                         by = c("type", "poly")][order(poly)],
            fill = TRUE) 
    } else {
      means.ls[, poly := NA]
    }
    setcolorder(means.ls, c("type", "poly"))

    means.l[[i]] <- means.ls

    if(resp.type == "normal") {
      mar.ls <-
        ls$landscape[type == "treatment",
                     .(diff = mean(treatment))]
      if(ls$landscape[type == "treatment", length(unique(poly))] > 1) {
        mar.ls <-
        rbind(mar.ls,
              ls$landscape[type == "treatment",
                           .(diff = mean(treatment)),
                           by = poly],
              fill = TRUE) 
      } else {
        mar.ls[, poly := NA]
      }
      setcolorder(mar.ls, "poly")
      marginal.l[[i]] <- mar.ls
    }

    if(resp.type == "binary") {
      marginal.l[[i]] <- ls$marginal
      marginal.l[[i]][,
                      `:=`(p.diff = p.trt - p.ref,
                           p.ratio = p.trt/p.ref)]
    }

    if(resp.type == "tweedie") {
      marginal.l[[i]] <- ls$marginal
      marginal.l[[i]][, diff := trt - ref]
    }

    obs.l[[i]] <-
      rbind(
        ls$landscape[type == "treatment",
                     .(n = .N, prop = .N/nrow(ls$landscape)),
                     by = c("type", "poly")],
        ls$landscape[,
                     .(n = .N, prop = .N/nrow(ls$landscape)),
                     by = c("type")],
        fill = TRUE)

    tb <- Sys.time()
    te <- tb-ta
    print(te)

  }
}

optim <- rbindlist(optim.l, idcol = "id")
cov <- rbindlist(cov.l, idcol = "id")
marginal <- rbindlist(marginal.l, idcol = "id", fill = TRUE)
marginal[, type := fifelse(is.na(poly), "treatment", "subarea")]
marginal[, type := factor(type, levels = c("treatment", "subarea"))]
setcolorder(marginal, c("id", "type", "poly"))
setorderv(marginal, c("id", "type", "poly"))
means <- rbindlist(means.l, idcol = "id", fill = TRUE)
setorderv(means, c("id", "type", "poly"))
obs <- rbindlist(obs.l, idcol = "id", fill = TRUE)
setorderv(obs, c("id", "type", "poly"))

results <-
  list(objectives = optim,
       marginal = marginal,
       balance = cov,
       means = means,
       obs = obs)

saveRDS(results, file.results)

print(results)



