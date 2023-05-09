library(data.table)

source("utilities.R")

args <- commandArgs(trailingOnly = TRUE)
ls.type <- args[1]

# ls.type <- "binary_imbalance_high"

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")

file.results <- paste0(path.ls, "summary.rds")

file.par <- paste0(path.ls, "parameters.rds")

parameters <- readRDS(file.par)

ids <- parameters$id

optim.l <- list()
marginal.l <- list()
scaling.l <- list()

for(i in ids) {

  if(file.exists(parameters[id == i, file.path])) {

    ta <- Sys.time()

    message(paste0("Summarizing results for landscape ", i, " …"))
    
    ls <- readRDS(parameters[id == i, file.path])
   
    imb <- ls$optim[names(ls$optim) %like% "imbalance"]

    optim.l[[i]] <- 
      data.table(imbalance = list(imb),
                 imbalance.mean = mean(imb),
                 area.prop = ls$optim["area.prop"])

    marginal.l[[i]] <- ls$marginal
    marginal.l[[i]][,
                    `:=`(p.diff = p.trt - p.ref,
                         p.ratio = p.trt/p.ref)]

    scaling.l[[i]] <- as.data.table(as.list(ls$scaling))

    tb <- Sys.time()
    te <- tb-ta
    print(te)

  }
}

optim <- rbindlist(optim.l, idcol = "id")
marginal <- rbindlist(marginal.l, idcol = "id", fill = TRUE)
scaling <- rbindlist(scaling.l, idcol = "id")

marginal[, type := fifelse(is.na(poly), "treatment", "subarea")]
marginal[, type := factor(type, levels = c("treatment", "subarea"))]
setcolorder(marginal, c("id", "type"))

saveRDS(list(objectives = optim,
             marginal = marginal,
             scaling = scaling),
        file.results)


