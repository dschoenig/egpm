library(data.table)
library(mgcv)
library(stringi)

source("utilities.R")

args <- commandArgs(trailingOnly = TRUE)
ls.type <- args[1]
mod.type <- args[2]

# ls.type <- "binary_imbalance_high"
# mod.type <- "egp_som25"

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
path.mod <- paste0(path.base, "models/", ls.type, "/")
path.results <- paste0(path.base, "results/", ls.type, "/")
if(!dir.exists(path.results)) dir.create(path.results, recursive = TRUE)

file.results <- paste0(path.results, mod.type, ".rds")
file.results.inter <- paste0(path.results, mod.type, ".inter.rds")
file.parameters <- paste0(path.ls, "parameters.rds")


para <- readRDS(file.parameters)

ids <- as.integer(para$file.name)

files.mod <- paste0(path.mod, mod.type, "_", para$file.name, ".rds")

if(file.exists(file.results.inter)) {
  res.int <- read.rds(file.results.inter)
  params.i <- res.int$params.i
  marginals.i <- res.int$marginals.i
  terms.i <- res.int$terms.i
  dev.expl.i <- res.int$dev.expl.i
  start.i <- res.int$i.proc
  rm(res.int)
} else {
  params.i <- list()
  marginals.i <- list()
  terms.i <- list()
  dev.expl.i <- list()
  start.i <- ids[1]
  # start.i <- 963
}

ids.proc <- ids[ids >= start.i]

block <- 25

ta <- Sys.time()
for(i in ids.proc) {

  if(i == ids.proc[1]) {
    message(paste0("Starting with model ", i, " …"))
  }

  mod.res <- NULL
  n <- 0
  while(is.null(mod.res)) {
    try({mod.res <- read.rds(files.mod[i])})
    if(is.null(mod.res)) {
      n <- n + 1
      if(n <= 10) {
        message(paste0("Loading model ", i, " failed. Trying again …"))
        Sys.sleep(10)
      } else {
        message("Loading model failed. Saving intermediary output …")
        int.out <-
          list(i.proc = i,
               params.i = params.i,
               marginals.i = marginals.i,
               terms.i = terms.i,
               dev.expl.i = dev.expl.i)
        saveRDS(int.out, file.results.inter)
        stop("Aborted due to read error")
      }
    }
  }

  marginals.j <- list()
  terms.j <- list()
  dev.expl.j <- list()

  for(j in seq_along(mod.res$models)) {

    grp.vars <- c("group.id", mod.res$models[[j]]$egp.def$group.by.c)

    marginals.j[[j]] <- 
      mod.res$models[[j]]$marginal[,
                                   .(mean = mean(marginal),
                                     median = median(marginal),
                                     se = sd(marginal),
                                     q0.5 = quantile(marginal, 0.005),
                                     q2.5 = quantile(marginal, 0.025),
                                     q5 = quantile(marginal, 0.05),
                                     q25 = quantile(marginal, 0.25),
                                     q75 = quantile(marginal, 0.75),
                                     q95 = quantile(marginal, 0.95),
                                     q97.5 = quantile(marginal, 0.975),
                                     q99.5 = quantile(marginal, 0.995)),
                                   by = grp.vars
                                   ]

    mod.egp <- mod.res$models[[j]]$model

    edf <- numeric(0)
    k.mod <- numeric(0)
    labels <- character(0)
    n.smooth <- length(mod.egp$smooth)
    for (k in 1:n.smooth) {
      labels[k] <- mod.egp$smooth[[k]]$label
      edf[k] <-
        sum(mod.egp$edf[mod.egp$smooth[[k]]$first.para:mod.egp$smooth[[k]]$last.para])
      k.mod[k] <- 
        mod.egp$smooth[[k]]$bs.dim
    }

    terms.j[[j]] <-
      data.table(term = labels,
                 edf = edf,
                 k = k.mod)

    dev.expl.j[[j]] <-
      data.table(dev.expl = 
                   ((mod.egp$null.deviance - mod.egp$deviance) /
                    mod.egp$null.deviance))

  }

  params.i[[i]] <- mod.res$parameters

  rm(mod.res)
  gc()

  params.i[[i]][, mod.id := 1:.N]
  marginals.i[[i]] <- rbindlist(marginals.j, idcol = "mod.id")
  terms.i[[i]] <- rbindlist(terms.j, idcol = "mod.id")
  dev.expl.i[[i]] <- rbindlist(dev.expl.j, idcol = "mod.id")

  if(i %% block == 0) {
    message(paste0("Model ", i, "/", length(ids), "."))
    int.out <-
      list(i.proc = i+1,
           params.i = params.i,
           marginals.i = marginals.i,
           terms.i = terms.i,
           dev.expl.i = dev.expl.i)
    saveRDS(int.out, file.results.inter)
    tb <- Sys.time()
    te <- tb-ta
    print(te)
    rm(tb, te)
    ta <- Sys.time()
  }

}

params <- rbindlist(params.i)
setcolorder(params.i[[i]], c("landscape", "mod.id"))
marginals <- rbindlist(marginals.i, idcol = "landscape")
setcolorder(marginals, c("landscape", "mod.id"))
terms <- rbindlist(terms.i, idcol = "landscape")
setcolorder(terms, c("landscape", "mod.id"))
dev.expl <- rbindlist(dev.expl.i, idcol = "landscape")
setcolorder(dev.expl, c("landscape", "mod.id"))

results <-
  list(parameters = params,
       marginal = marginals,
       terms = terms,
       deviance = dev.expl)

saveRDS(results, file.results)


