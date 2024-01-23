library(data.table)
library(stringi)

source("utilities.R")

## Paths and parameters

response.types <-
  c("normal", "tweedie", "binary")

imbalance.types <- c("high", "low")



landscapes <-
  CJ(trt.effect = c(FALSE, TRUE),
     response = response.types, 
     imbalance = imbalance.types,
     sorted = FALSE)
landscapes[trt.effect == TRUE,
           name := fifelse(response == "normal",
                           paste0("imbalance_", imbalance),
                           paste0(response, "_imbalance_", imbalance))]
landscapes[trt.effect == FALSE,
           name := fifelse(response == "normal",
                           paste0("noeff_", "imbalance_", imbalance),
                           paste0("noeff_", response, "_imbalance_", imbalance))]


path.base <- "../"
path.landscapes <- paste0(path.base, "landscapes/")
path.results <- paste0(path.base, "results/")
file.ls.sum <- paste0(path.results, "landscapes.rds")


## Get landscape info

ls.objectives.l <- list()
ls.marginal.l <- list()
ls.balance.l <- list()
ls.means.l <- list()
ls.obs.l <- list()

cols.old <-
  c("id", "type", "poly", "n", "prop")
cols.new <-
  c("ls.id", "area.type", "subarea.id", "treatment.n", "treatment.prop")


for(i in 1:nrow(landscapes)) {

  ls.file <- paste0(path.landscapes, landscapes$name[i], "/summary.rds")
  ls.sum <- readRDS(ls.file)

  ls.sum <-
    lapply(ls.sum,
           \(x) x[,
                  `:=`(trt.effect = landscapes$trt.effect[i],
                       ls.name = landscapes$name[i],
                       ls.response = landscapes$response[i],
                       ls.imbalance = landscapes$imbalance[i])])
                     
  ls.objectives.l[[i]] <- ls.sum$objectives
  ls.marginal.l[[i]] <- ls.sum$marginal
  ls.balance.l[[i]] <- ls.sum$balance
  ls.means.l[[i]] <- ls.sum$means
  ls.obs.l[[i]] <- ls.sum$obs

}

ls.all <-
  list(objectives = ls.objectives.l,
       marginal = ls.marginal.l,
       balance = ls.balance.l,
       means = ls.means.l,
       obs = ls.obs.l) |>
  lapply(rbindlist, fill = TRUE) |>
  lapply(setnames,
         old = cols.old,
         new = cols.new,
         skip_absent = TRUE)



for(i in seq_along(ls.all)) {

  x <- ls.all[[i]]

  x[,
    `:=`(
         ls.response = factor(ls.response,
                              levels = c("normal", "tweedie", "binary")),
         ls.imbalance = factor(ls.imbalance,
                               levels = c("low", "high")),
         ls.id = factor(ls.id))]

  ls.type.lev <-
    with(x,
         paste0(rep(levels(ls.response), each = 2), "_",
                rep(levels(ls.imbalance), times = 3)))

  x[,
    ls.type := factor(paste0(ls.response, "_", ls.imbalance),
                      levels = ls.type.lev)]

  x[trt.effect == TRUE,
            ls.uid :=
              as.integer((as.integer(ls.type)-1) * 1000) +
              as.integer(as.character(ls.id))]
  x[trt.effect == FALSE,
            ls.uid :=
              as.integer((length(unique(ls.type))) * 1000) +
              as.integer((as.integer(ls.type)-1) * 1000) +
              as.integer(as.character(ls.id))]

  ls.all[[i]] <- x

}

ls.all$marginal[,
                area.type := factor(area.type,
                                    levels = c("treatment", "subarea"))]
ls.all$marginal[,
                mar.true := fcase(ls.response == "normal", diff,
                                  ls.response == "tweedie", diff,
                                  ls.response == "beta", diff,
                                  ls.response == "binary", p.diff)]
ls.all$means[,
             area.type := factor(area.type,
                                 levels = c("treatment", "subarea"))]
ls.all$obs[,
           area.type := factor(area.type,
                               levels = c("treatment", "subarea"))]


ls.all <-
  lapply(ls.all,
         \(x) setorder(x, -trt.effect, ls.response, ls.imbalance, ls.id))

ls.cols <- c("trt.effect",
             "ls.uid", "ls.type", "ls.name",
             "ls.response", "ls.imbalance", "ls.id")
             # "area.type", "subarea.id",
             # "response.mean", "response.sd",
             # "treatment.n", "treatment.prop")
ls.all <-
  lapply(ls.all,
         \(x) setcolorder(x, ls.cols))

saveRDS(ls.all, file.ls.sum)
