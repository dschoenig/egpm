library(data.table)
library(mgcv)
library(kohonen)
library(MatchIt)
library(ranger)
library(ggdist)

source("utilities_fcne.R")
source("utilities.R")


# args <- commandArgs(trailingOnly = TRUE)
# ls.type <- args[1]
# mod.type <- args[2]

ls.type <- "1000_4cov_nl"
egp.type <- "egp_sam0.01_som50"
match.type <- "match_sam0.01"
rf.type <- "rf_sam0.01_t1000"

# ls.type <- "100_4cov_nl"
# egp.type <- "egp_sam1_som50"
# match.type <- "match_sam1"


path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
path.mod <- paste0(path.base, "models/", ls.type, "/")
path.results <- paste0(path.base, "results/", ls.type, "/")

file.egp <- paste0(path.results, egp.type, ".rds")
file.match <- paste0(path.results, match.type, ".rds")
file.rf <- paste0(path.results, rf.type, ".rds")

eff.egp <- readRDS(file.egp)
eff.match <- readRDS(file.match)
eff.rf <- readRDS(file.rf)

effects <-
  rbind(
    rbind(eff.egp$noint$effects[, .(id, interactions = FALSE, method = "egp", mean, q2.5, q97.5)],
          eff.match$noint$effects[, .(id, interactions = FALSE, method, mean, q2.5, q97.5)],
          eff.rf$noint$effects[, .(id, interactions = FALSE, method = "rf", mean, q2.5 = NA, q97.5 = NA)]),
    rbind(eff.egp$int$effects[, .(id, interactions = TRUE, method = "egp", mean, q2.5, q97.5)],
          eff.match$int$effects[, .(id, interactions = TRUE, method, mean, q2.5, q97.5)],
          eff.rf$int$effects[, .(id, interactions = TRUE, method = "rf", mean, q2.5 = NA, q97.5 = NA)]))

effects[, method := factor(method, levels = c("rf",
                                              "nn.mh.re", "nn.mh.nr",
                                              "nn.ps.re", "nn.ps.nr",
                                              "cem.st", "cem.sc", "cem.fd",
                                              "egp",
                                              "lmcov", "lm"))]

# effects[, dist := abs(1 - mean)]
# effects[mean < 0, length(mean), by = c("method", "interactions")]
# effects[, sqrt(mean((1-mean)^2)), c("method", "interactions")]
# effects[, sd(mean), c("method", "interactions")]
effects <- effects[method %in% c("lm", "egp", "cem.st", "nn.ps.re", "nn.mh.re", "rf")]

plot.noint <- 
ggplot(effects[interactions == FALSE]) +
# geom_violin(aes(x = method, y = mean))
stat_halfeye(aes(y = method, x = mean, fill_ramp = stat(cut_cdf_qi(cdf))),
             alpha = 0.8,
             point_interval = mean_qi,
             n = 1001) +
geom_vline(xintercept = 0, alpha = 0.8) +
geom_vline(xintercept = 1, colour = 2, alpha = 0.8) +
geom_vline(xintercept = effects[interactions == FALSE & method == "egp", mean(mean)],
           colour = 1, linetype = "dashed", alpha = 0.8) +
# stat_gradientinterval(aes(x = method, y = mean)) +
# stat_dots(aes(y = method, x = mean)) +
scale_fill_ramp_discrete(na.translate = FALSE, from = "grey85", range = c(1,0)) +
scale_fill_discrete_qualitative("Set2") +
coord_cartesian(xlim = c(-2, 3)) +
theme_ggdist()

plot.int <- 
ggplot(effects[interactions == TRUE]) +
# geom_violin(aes(x = method, y = mean))
stat_halfeye(aes(y = method, x = mean, fill_ramp = stat(cut_cdf_qi(cdf))),
             alpha = 0.8,
             point_interval = mean_qi,
             n = 1001) +
geom_vline(xintercept = 0, alpha = 0.8) +
geom_vline(xintercept = 1, colour = 2, alpha = 0.8) +
geom_vline(xintercept = effects[interactions == TRUE & method == "egp", mean(mean)],
           colour = 1, linetype = "dashed", alpha = 0.8) +
# stat_gradientinterval(aes(x = method, y = mean)) +
# stat_dots(aes(y = method, x = mean)) +
scale_fill_ramp_discrete(na.translate = FALSE, from = "grey85", range = c(1,0)) +
scale_fill_discrete_qualitative("Set2") +
coord_cartesian(xlim = c(-2, 3)) +
theme_ggdist()

plot.noint + plot.int

effects[interactions == TRUE, sd(mean), method]

library(brms)
library(bayesplot)

effects.br <- effects[interactions == TRUE]
effects.br[, method := relevel(method, "egp")]

priors <- c(prior(cauchy(0, 1), class = sd),
            prior(student_t(3, 0 , 1), class = b),
            prior(student_t(3, 0 , 1), class = b, dpar = sigma))

mod.brm <- brm(bf(mean ~ method + (1|id), sigma ~ method),
               prior = priors,
               data = effects.br,
               chains = 4,
               cores = 4,
               warmup = 1000,
               iter = 2000)


comp.methods <- c("egp", "rf", "nn.ps.re", "nn.mh.re", "cem.st")

pred.brm <- predict(mod.brm, summary = FALSE)
pred.dt <- as.data.table(t(pred.brm))
setnames(pred.dt, paste0("draw.", 1:ncol(pred.dt)))
pred.dt <- cbind(effects.br[, .(id, interactions, method)], pred.dt)
pred.dt <-
  melt(pred.dt, measure.vars = grep("draw", names(pred.dt)),
      variable.name = "draw", value.name = "yhat")
pred.dt[, draw := as.integer(stri_replace_first_fixed(draw, "draw.", ""))]

# Difference in relative bias
comp.bias <-
  pred.dt[, .(bias = abs(1 - mean(yhat) / 1)), by = c("method", "draw")] |>
  dcast(draw ~ method, value.var = "bias")
comp.bias <-
  cbind(comp.bias[, .(draw)],
        comp.bias[, lapply(.SD, \(x) egp - x), .SDcols = comp.methods]) |>
  melt(id.vars = "draw", variable.name = "method", value.name = "bias.diff")


# Difference in coefficent of variation
comp.sd <-
  pred.dt[, .(sd = sd(yhat)), by = c("method", "draw")] |>
  dcast(draw ~ method, value.var = "sd")
comp.sd <-
  cbind(comp.sd[, .(draw)],
        comp.sd[, lapply(.SD, \(x) (egp - x) / 1), .SDcols = comp.methods]) |>
  melt(id.vars = "draw", variable.name = "method", value.name = "sd.diff")

# Difference in false negative rate
comp.neg <-
  pred.dt[,
          .(neg.prop = sum(ifelse(yhat < 0, 1, 0)) / .N),
          by = c("method", "draw")] |>
  dcast(draw ~ method, value.var = "neg.prop")
comp.neg <-
  cbind(comp.neg[, .(draw)],
        comp.neg[, lapply(.SD, \(x) egp -x), .SDcols = comp.methods]) |>
  melt(id.vars = "draw", variable.name = "method", value.name = "neg.prop")

comp.senst

ggplot(comp.bias[method != "egp"]) +
# geom_violin(aes(x = method, y = mean))
stat_halfeye(aes(y = method, x = bias.diff, fill_ramp = stat(cut_cdf_qi(cdf))),
             alpha = 0.8,
             point_interval = mean_qi,
             n = 1001) +
geom_vline(xintercept = 0, colour = 2, alpha = 0.8) +
scale_fill_ramp_discrete(na.translate = FALSE, from = "grey85", range = c(1,0)) +
scale_fill_discrete_qualitative("Set2") +
scale_x_continuous(limits = c(-0.55, 0.25), breaks = seq(-0.5, 0.5, 0.1)) +
# theme_default()
theme_ggdist()

ggplot(comp.sd[method != "egp"]) +
# geom_violin(aes(x = method, y = mean))
stat_halfeye(aes(y = method, x = sd.diff, fill_ramp = stat(cut_cdf_qi(cdf))),
             alpha = 0.8,
             point_interval = mean_qi,
             n = 1001) +
geom_vline(xintercept = 0, colour = 2, alpha = 0.8) +
scale_fill_ramp_discrete(na.translate = FALSE, from = "grey85", range = c(1,0)) +
scale_fill_discrete_qualitative("Set2") +
scale_x_continuous(limits = c(-0.25, 0.55), breaks = seq(-0.5, 1.5, 0.1)) +
# coord_cartesian(xlim = c(-1, 1)) +
# theme_default()
theme_ggdist()

ggplot(comp.neg[method != "egp"]) +
# geom_violin(aes(x = method, y = mean))
stat_halfeye(aes(y = method, x = neg.prop, fill_ramp = stat(cut_cdf_qi(cdf))),
             alpha = 0.8,
             point_interval = mean_qi,
             n = 1001) +
geom_vline(xintercept = 0, colour = 2, alpha = 0.8) +
scale_fill_ramp_discrete(na.translate = FALSE, from = "grey85", range = c(1,0)) +
scale_fill_discrete_qualitative("Set2") +
# coord_cartesian(xlim = c(-1, 1)) +
# theme_default()
theme_ggdist()

