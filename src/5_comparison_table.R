library(data.table)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.rds")

file.egp.match.perm <- paste0(path.comp, "egp.match.perm.rds")
file.egp.match.table <- paste0(path.comp, "egp.match.csv")


estimates <- readRDS(file.estimates)


# 1. EGP vs matching (and GLM) ########################################


egp.match.perm <- readRDS(file.egp.match.perm)

sub.dt <-
  CJ(
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

egp.match.global <- 
  compare_performance(estimates.sub,
                      by.method = "name.short",
                      comparisons = c("EGP", "GLM"))
egp.match.ls <- 
  compare_performance(estimates.sub,
                      by.method = "name.short",
                      by.landscape = c("ls.response", "ls.imbalance"),
                      comparisons = c("EGP", "GLM"))


ls.response.lev <- c("all", levels(estimates.sub$ls.response))
ls.imbalance.lev <- c("all", levels(estimates.sub$ls.imbalance))

egp.match.global[, 
                 `:=`(ls.response = factor(rep("all", .N), levels = ls.response.lev),
                      ls.imbalance = factor(rep("all", .N), levels = ls.imbalance.lev))]
egp.match.ls <- set_levels_dt(egp.match.ls, egp.match.global)

egp.match <-
  rbind(egp.match.global, egp.match.ls, use.names = TRUE) |>
  merge(egp.match.perm)
setcolorder(egp.match, c("ls.response", "ls.imbalance", "name.short"))
setorder(egp.match, ls.response, ls.imbalance, name.short)


egp.match.table <- format_table(egp.match)
fwrite(egp.match.table, file.egp.match.table)




sub.dt <-
  CJ(
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     # area.type = "subarea",
     area.type = "treatment",
     mod.name = c("egp_som25", "egp_som25_unequal"),
     egp.som.topology = c("hexagonal", "rectangular"),
     egp.cf.nb = c("sequential"),
     egp.geo.w = c(TRUE),
     # egp.cf.nb = c("expand", "sequential"),
     # egp.geo.w = c(FALSE, TRUE),
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)
estimates.sub

estimates.sub[, name.comp := paste(mod.name, egp.som.topology, egp.cf.nb, egp.geo.w, sep = ".")]

estimates.sub

method.comp <-
  CJ(ls.uid = 1:6000,
     method1 = unique(estimates.sub$name.comp),
     method2 = unique(estimates.sub$name.comp))

method1 <-
  estimates.sub[,
                .(ls.response, ls.imbalance, ls.uid, 
                  method1 = name.comp,
                  est1 = mar.std)
                ]
method2 <-
      estimates.sub[,
                    .(ls.uid, 
                      method2 = name.comp,
                      est2 = mar.std)
                    ]

methods <- 
  method2[method1[method.comp, on = c("ls.uid", "method1")], on = c("ls.uid", "method2")]

ggplot(methods) +
  geom_point(aes(x = est1, y = est2, colour = ls.response),
             alpha = 0.1, shape = 16) +
  facet_grid(rows = vars(method2), cols = vars(method1))
          



compare_performance(estimates.sub,
                    by.method = "name.comp",
                    # by.landscape = c("ls.response", "ls.imbalance"),
                    # comparisons = c("rect_seq_geo"),
                    comparisons = NULL
                    ) |>
# _[order(ls.response, ls.imbalance, ser), .(ls.response, ls.imbalance, name.comp, ser)]
_[order(ser), .(name.comp, ser)]
# _[order(abs(bias)), .(name.comp, bias)]
# _[order(rmse), .(name.comp, rmse)]



sub.dt <-
  CJ(
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     # area.type = "subarea",
     area.type = "treatment",
     mod.name = c("egp_som25"),
     egp.som.topology = c("hexagonal", "rectangular"),
     egp.cf.nb = c("sequential"),
     egp.geo.w = c(TRUE),
     # egp.cf.nb = c("expand", "sequential"),
     # egp.geo.w = c(FALSE, TRUE),
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)
estimates.sub

names.type <-
  list(
       data.table(egp.som.topology = rep("hexagonal", 4),
                  egp.cf.nb = rep(c("expand", "sequential"), each = 2),
                  egp.geo.w = rep(c(FALSE, TRUE), times = 2),
                  name.type = c("hex_ex_ng", "hex_ex_geo", "hex_seq_ng", "hex_seq_geo")),
       data.table(egp.som.topology = rep("rectangular", 4),
                  egp.cf.nb = rep(c("expand", "sequential"), each = 2),
                  egp.geo.w = rep(c(FALSE, TRUE), times = 2),
                  name.type = c("rect_ex_ng", "rect_ex_geo", "rect_seq_ng", "rect_seq_geo"))) |>
  rbindlist(fill = TRUE)

names.type[, name.type := factor(name.type, levels = name.type)]
estimates.sub <- merge(estimates.sub, names.type)

compare_performance(estimates.sub,
                    by.method = "name.type",
                    # by.landscape = c("ls.response", "ls.type"),
                    comparisons = c("rect_seq_geo")) |>
_[, c(1:3, 7:9)]

unique(estimates$mod.name)

sub.dt <-
  CJ(
     ls.response = c("normal"),
     ls.imbalance = c("high"),
     area.type = c("treatment", "subarea"),
     mod.name = c("egp_som25"),
     egp.som.topology = c("rectangular"),
     egp.cf.nb = c("sequential"),
     egp.geo.w = c(TRUE),
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)
estimates.sub

compare_performance(estimates.sub,
                    by.method = "area.type",
                    by.landscape = c("ls.response", "ls.type"),
                    comparisons = c("treatment"))




sub.dt <-
  CJ(
     ls.response = c("normal"),
     ls.imbalance = c("high"),
     area.type = "treatment",
     mod.name = c("egp_som25", "egp_som25_unequal"),
     egp.som.topology = c(NA, "rectangular"),
     egp.cf.nb = c(NA, "sequential"),
     egp.geo.w = c(NA, TRUE),
     match.mod.cov = factor(c(NA, "interact")),
     match.trt.int = c(NA, FALSE),
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)
estimates.sub


compare_performance(estimates.sub,
                    by.method = "mod.name",
                    by.landscape = c("ls.response", "ls.type"),
                    comparisons = c("egp_som25"))
