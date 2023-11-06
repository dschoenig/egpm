library(data.table)
library(stringi)

source("utilities.R")

## Paths and parameters

# models <- c(
#             "egp_som25",
#             "egp_som10",
#             "egp_som50",
#             "egp_som25_sam005",
#             "egp_som25_sam02",
#             "egp_som25_unequal",
#             "match",
#             "match_sam005",
#             "match_sam02")

models <- c("egp_som25",
            "egp_som10",
            "egp_som50",
            "egp_som25_unequal",
            "match")

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
path.comp <- paste0(path.results, "comparison/")
file.estimates.csv <- paste0(path.comp, "estimates.csv")
file.estimates.rds <- paste0(path.comp, "estimates.rds")


## Get landscape info

ls.mar.l <- list()
ls.mean.l <- list()
ls.obs.l <- list()

for(i in 1:nrow(landscapes)) {

  ls.file <- paste0(path.landscapes, landscapes$name[i], "/summary.rds")
  ls.sum <- readRDS(ls.file)

  # ls.obs <- ls.sum$obs

  # ls.obs <- ls.obs[type == "treatment"]
  # ls.obs[, poly.max := max(poly, na.rm = TRUE), by = id]
  # ls.obs <- copy(ls.obs[is.na(poly) | (!is.na(poly) & poly.max > 1)])

  # ls.obs[, type := fifelse(is.na(poly), "treatment", "subarea")]
  # ls.obs[, type := factor(type, levels = c("treatment", "subarea"))]
  # setcolorder(ls.obs, c("id", "type", "poly"))
  # setorderv(ls.obs, c("id", "type", "poly"))

  # ls.sum$obs <- ls.obs

  # saveRDS(ls.sum, ls.file)

  ls.mar.l[[i]] <-
    ls.sum$marginal[,
                    .(trt.effect = landscapes$trt.effect[i],
                      ls.name = landscapes$name[i],
                      ls.id = id,
                      area.type = type,
                      subarea.id = poly,
                      mar.true = switch(landscapes$response[i],
                                        normal = diff,
                                        tweedie = diff,
                                        binary = p.diff))]

  ls.mean.l[[i]] <-
    ls.sum$mean[is.na(type) & is.na(poly),
                .(trt.effect = landscapes$trt.effect[i],
                  ls.name = landscapes$name[i],
                  ls.id = id,
                  response.mean,
                  response.sd)]

  ls.obs.l[[i]] <-
    ls.sum$obs[,
               .(trt.effect = landscapes$trt.effect[i],
                 ls.name = landscapes$name[i],
                 ls.id = id,
                 area.type = type,
                 subarea.id = poly,
                 treatment.n = n,
                 treatment.prop = prop)]

  # print(nrow(ls.sum$marginal) == nrow(ls.sum$obs))

}

ls.mar <-
  merge(rbindlist(ls.mar.l),
        rbindlist(ls.obs.l), all = TRUE) |>
  merge(rbindlist(ls.mean.l), by = c("trt.effect", "ls.name", "ls.id")) |>
  merge(landscapes[,
                   .(trt.effect,
                     ls.name = name,
                     ls.response = response,
                     ls.imbalance = imbalance)],
        by = c("trt.effect", "ls.name")
        )


## Get model info

mod.mar.i <- list()

for(i in 1:nrow(landscapes)) {

  files.all <- list.files(paste0(path.results, landscapes$name[i]))
  mod.sel <- models[paste0(models, ".rds") %in% files.all]
  files.load <- 
    paste0(path.results,
           landscapes$name[i], "/",
           paste0(mod.sel, ".rds"))

  mod.mar.j <- list()

  for(j in seq_along(mod.sel)) {

    res <- readRDS(files.load[j])

    res.est <-
      merge(res$parameters,
            res$marginal,
            by = c("landscape", "mod.id"),
            all = TRUE)
    setnames(res.est, c("landscape"), c("ls.id"))

    res.est[,
            `:=`(trt.effect = landscapes$trt.effect[i],
                 ls.name = landscapes$name[i],
                 mod.name = mod.sel[j])]

    mod.mar.j[[j]] <- res.est

  }

  mod.mar.i[[i]] <- rbindlist(mod.mar.j, fill = TRUE)

}

mod.mar <- rbindlist(mod.mar.i)
mod.mar[,
          `:=`(area.type = ifelse(is.na(poly),
                                   "treatment",
                                   "subarea"),
               subarea.id = poly)]        
mod.mar <- mod.mar[, -c("group.id", "poly")]
setnames(mod.mar,
         c("mean", "median", "se", "q0.5", "q2.5", "q5", "q25", "q75", "q95", "q97.5", "q99.5"),
         paste0("mar.", c("est", "median", "se", "q0.5", "q2.5", "q5", "q25", "q75", "q95", "q97.5", "q99.5")))
         # c("mean", "se", "q2.5", "q97.5"),
         # paste0("mar.", c("est", "se", "q2.5", "q97.5")))
setnames(mod.mar,
         c("type", "cf.nb.strategy", "geo",
           "som.dim", "som.epochs", "som.topology",
           "mod.cov", "trt.int", "cutpoints",
           "distance", "replace"),
         c("method", "egp.cf.nb", "egp.geo.w",
           "egp.som.dim", "egp.som.epochs", "egp.som.topology",
           "match.mod.cov", "match.trt.int", "match.cutpoints",
           "match.distance", "match.replace"))


# Model abbreviations

names.short <-
  list(
       data.table(method = "egp",
                  name.short = "EGP"),
       data.table(method = "glm",
                  name.short = "GLM"),
       data.table(method = rep("cem", 3),
                  match.cutpoints = c("st", "sc", "fd"),
                  name.short = c("CEM-ST", "CEM-SC", "CEM-FD")),
       data.table(method = rep("nearest", 4),
                  match.distance = rep(c("glm", "mahalanobis"), each = 2),
                  match.replace = rep(c(FALSE, TRUE), times = 2),
                  name.short = c("NN-PS-NR", "NN-PS-RE",
                                 "NN-MA-NR", "NN-MA-RE"))) |>
  rbindlist(fill = TRUE)
names.short[,
            name.short := factor(name.short,
                                 levels = c("EGP", "GLM",
                                            "CEM-ST", "CEM-SC", "CEM-FD",
                                            "NN-PS-NR", "NN-PS-RE",
                                            "NN-MA-NR", "NN-MA-RE"))]



# Combine estimates and landscape info, add short names

estimates <-
  merge(mod.mar, ls.mar) |>
  merge(names.short,
        by = c("method", "match.distance", "match.replace", "match.cutpoints"),
        all.x = TRUE)

mar.cols <- names(estimates)
mar.cols[mar.cols %like% "mar."]

estimates[trt.effect == TRUE, mar.std := mar.est / mar.true]
estimates[trt.effect == FALSE, mar.std := mar.est / response.sd]


estimates[,
          `:=`(egp.som.topology = factor(egp.som.topology,
                                         levels = c("rectangular",
                                                    "hexagonal")),
               egp.cf.nb = factor(egp.cf.nb,
                                  levels = c("sequential",
                                             "expand")),
               match.distance = factor(match.distance,
                                       levels = c("glm",
                                                  "mahalanobis")),
               match.cutpoints = factor(match.cutpoints,
                                        levels = c("st", "sc", "fd")),
               match.mod.cov = factor(match.mod.cov,
                                      levels = c("exclude",
                                                 "include",
                                                 "interact")),
               area.type = factor(area.type,
                                  levels = c("treatment", "subarea")),
               ls.response = factor(ls.response,
                                    levels = c("normal", "tweedie", "binary")),
               ls.imbalance = factor(ls.imbalance,
                                     levels = c("low", "high")),
               ls.id = factor(ls.id))]

estimates[method == "egp",
          egp.som.unequal := fifelse(mod.name %like% "unequal", TRUE, FALSE)]


ls.type.lev <-
  with(estimates,
       paste0(rep(levels(ls.response), each = 2), "_",
              rep(levels(ls.imbalance), times = 3)))

estimates[,
          ls.type := factor(paste0(ls.response, "_", ls.imbalance),
                                levels = ls.type.lev)]

estimates[trt.effect == TRUE,
          ls.uid :=
            as.integer((as.integer(ls.type)-1) * 1000) +
            as.integer(as.character(ls.id))]

estimates[trt.effect == FALSE,
          ls.uid :=
            as.integer((length(unique(ls.type))) * 1000) +
            as.integer((as.integer(ls.type)-1) * 1000) +
            as.integer(as.character(ls.id))]

setorder(estimates, -trt.effect, ls.response, ls.imbalance, ls.id, name.short)

est.names <- names(estimates)

ls.cols <- c("trt.effect",
             "ls.uid", "ls.type", "ls.name",
             "ls.response", "ls.imbalance", "ls.id",
             "area.type", "subarea.id", "sam.frac",
             "response.mean", "response.sd",
             "treatment.n", "treatment.prop")
mod.cols <- c("name.short", "method", "mod.name", "mod.id")
egp.cols <- sort(est.names[est.names %like% "egp."])
match.cols <- est.names[est.names %like% "match."]
mar.cols <- est.names[est.names %like% "mar."]
setcolorder(estimates, c(ls.cols, mod.cols, egp.cols, match.cols, mar.cols))

setnames(estimates, c("p.pos", "p.neg"), c("prob.pos", "prob.neg"))

fwrite(estimates, file.estimates.csv, yaml = TRUE, na = "NA")
saveRDS(estimates, file.estimates.rds)

estimates[area.type == "treatment", .N, by = ls.name]

