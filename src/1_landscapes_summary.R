library(data.table)

source("utilities.R")

args <- commandArgs(trailingOnly = TRUE)
ls.type <- args[1]
resp.type <- args[2]

# ls.type <- "binary_imbalance_high"
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
marginal.l <- list()
scaling.l <- list()

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


    if(resp.type == "normal") {
      marginal.l[[i]] <-
        rbind(ls$landscape[type == "treatment",
                           .(marginal = mean(treatment))],
              ls$landscape[type == "treatment",
                           .(marginal = mean(treatment)),
                           by = poly],
              fill = TRUE) |>
        setcolorder("poly")
    }

    if(resp.type == "binary") {
      marginal.l[[i]] <- ls$marginal
      marginal.l[[i]][,
                      `:=`(p.diff = p.trt - p.ref,
                           p.ratio = p.trt/p.ref)]
      scaling.l[[i]] <- as.data.table(as.list(ls$scaling))
    }

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
setcolorder(marginal, c("id", "type"))

if(resp.type == "normal") {
  results <-
    list(objectives = optim,
         marginal = marginal,
         balance = cov)
}

if(resp.type == "binary") {
  scaling <- rbindlist(scaling.l, idcol = "id")
  results <-
    list(objectives = optim,
         marginal = marginal,
         scaling = scaling,
         balance = cov)
}

print(results)

saveRDS(results, file.results)


