library(data.table)
library(stringi)

source("utilities.R")

## Paths and parameters

models <- c(
            "egp_som25",
            "egp_som10",
            "egp_som50",
            "egp_som25_sam005",
            "egp_som25_sam02",
            "egp_som25_unequal",
            "match",
            "match_sam005",
            "match_sam02")

response.types <-
  c("normal", "tweedie", "binary")

imbalance.types <- c("high", "low")



landscapes <-
  CJ(response = response.types, 
     imbalance = imbalance.types,
     sorted = FALSE)
landscapes[,
           name := fifelse(response == "normal",
                           paste0("imbalance_", imbalance),
                           paste0(response, "_imbalance_", imbalance))]


path.base <- "../"
path.landscapes <- paste0(path.base, "landscapes/")
path.results <- paste0(path.base, "results/")
file.estimates <- paste0(path.results, "estimates.csv")


## Get landscape info

ls.mar.l <- list()

for(i in 1:nrow(landscapes)) {

  ls.file <- paste0(path.landscapes, landscapes$name[i], "/summary.rds")
  ls.sum <- readRDS(ls.file)

  ls.mar.l[[i]] <-
    ls.sum$marginal[,
                    .(ls.name = landscapes$name[i],
                      ls.id = id,
                      area.type = type,
                      subarea.id = poly,
                      mar.true = switch(landscapes$response[i],
                                        normal = diff,
                                        tweedie = diff,
                                        binary = p.diff))]

}

ls.mar <-
  merge(rbindlist(ls.mar.l),
        landscapes[,
                   .(ls.name = name,
                     ls.response = response,
                     ls.imbalance = imbalance)])


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
            `:=`(ls.name = landscapes$name[i],
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
         c("mean", "se", "q2.5", "q97.5"),
         paste0("mar.", c("est", "se", "q2.5", "q97.5")))
setnames(mod.mar,
         c("type", "cf.nb.strategy", "geo",
           "som.dim", "som.epochs", "som.topology",
           "mod.cov", "trt.int", "cutpoints",
           "distance", "replace"),
         c("method", "egp.cf.nb", "egp.geo.w",
           "egp.som.dim", "egp.som.epochs", "egp.som.topology",
           "match.mod.cov", "match.trt.int", "match.cutpoints",
           "match.distance", "match.replace"))


# Combine estimates and landscape info

estimates <- merge(mod.mar, ls.mar)

est.names <- names(estimates)

ls.cols <- c("ls.name", "ls.response", "ls.imbalance", "ls.id",
             "area.type", "subarea.id", "sam.frac")
mod.cols <- c("method", "mod.name", "mod.id")
egp.cols <- est.names[est.names %like% "egp."]
match.cols <- est.names[est.names %like% "match."]
mar.cols <- est.names[est.names %like% "mar."]

setcolorder(estimates, c(ls.cols, mod.cols, egp.cols, match.cols, mar.cols))

estimates[, mar.std := mar.est / mar.true]

fwrite(estimates, file.estimates, yaml = TRUE, na = "NA")




egp.comp <-
  CJ(ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("high", "low"),
     area.type = "treatment",
     mod.name = c("egp_som25"),
     sorted = FALSE)

estimates[egp.comp,
          on = names(egp.comp),
          nomatch = NULL
          ][ls.response == "binary" & ls.imbalance == "high",
            .(
              bias.std = mean(mar.std-1),
              cv = sd(mar.std),
              ser = sum(sign(mar.est) != sign(mar.true)) / .N,
              rmse = rmse(mar.std, 1)),
              by = .(
                     # ls.response, ls.imbalance,
                     mod.name,
                     egp.som.topology, egp.cf.nb, egp.geo.w)
              ]

sub.dt <-
  CJ(ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("high", "low"),
     area.type = c("treatment", "subarea"),
     mod.name = c("egp_som25", "match"),
     egp.som.topology = c(NA, "rectangular"),
     egp.cf.nb = c(NA, "sequential"),
     egp.geo.w = c(NA, TRUE),
     match.mod.cov = c(NA, "interact"),
     match.trt.int = c(NA, FALSE),
     sorted = FALSE)

estimates[sub.dt,
          on = names(sub.dt),
          nomatch = NULL
          ][
            # ls.imbalance == "low",
            # ls.response == "tweedie" & ls.imbalance == "low",
            ,
            .(bias = mean(mar.std-1),
              cv = sd(mar.std)/mean(mar.std),
              ser = sum(sign(mar.est) != sign(mar.true)) / .N,
              rmse = rmse(mar.std, 1)),
              by = .(mod.name, mod.id,
                     mod.type
                     , match.distance, match.cutpoints, match.mod.cov, match.trt.int, match.replace
                     )][order(rmse)]


estimates.sub <-
  estimates[sub.dt,
            on = names(sub.dt),
            nomatch = NULL]

lab.mod.method <- c(egp = "EGP", glm = "GLM", cem = "CEM", nearest = "NN")
lab.cutpoints <- c(fd = "FD", sc = "SC", st = "ST")
lab.distance <- c(mahalanobis = "MA", glm = "PS")
lab.replace <- c("FALSE" = "NR", "TRUE" = "RE")
lab.levels <- c("EGP", "GLM", 
                "CEM-FD", "CEM-SC", "CEM-ST",
                "NN-MA-NR", "NN-MA-RE", "NN-PS-NR", "NN-PS-RE")
estimates.sub[, method.char := as.character(method)]
estimates.sub[,
              mod.label := fcase(method.char %in% c("egp", "glm"),
                                 paste(lab.mod.method[method.char],
                                       sep = "-"),
                                 method.char == "cem",
                                 paste(lab.mod.method[method.char],
                                       lab.cutpoints[as.character(match.cutpoints)],
                                       sep = "-"),
                                 method.char == "nearest",
                                 paste(lab.mod.method[method.char],
                                       lab.distance[as.character(match.distance)],
                                       lab.replace[as.character(match.replace)],
                                       sep = "-"))]
estimates.sub[, mod.label := factor(mod.label, levels = lab.levels)]



estimates.sub[area.type == "treatment"] |>
ggplot() +
  stat_slabinterval(aes(x = mar.std, y = mod.label, fill = mod.label),
                    # slab_linewidth = 0.5,
                    p_limits = c(0.001, 0.999),
                    point_interval = "mean_qi") +
  geom_vline(xintercept = 0, colour = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = 1) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(-0.5, 2.5)) +
  theme_ggdist()
  # theme_ggdist() +
  # facet_grid(rows = vars(ls.imbalance), cols = vars(ls.response))


estimates.sub[mod.label %in% c("EGP", "NN-MA-RE")] |>
ggplot() +
  stat_slabinterval(aes(x = mar.std, y = area.type, fill = mod.label),
                    p_limits = c(0.001, 0.999),
                    point_interval = "mean_qi") +
  geom_vline(xintercept = 0, colour = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = 1) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(xlim = c(-1, 2)) +
  facet_wrap(vars(mod.label)) +
  # facet_grid(rows = vars(ls.response), cols = vars(ls.imbalance)) +
  theme_ggdist()