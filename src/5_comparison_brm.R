library(data.table)
library(brms)

source("utilities.R")

# overwrite <- TRUE

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.csv")

file.mod <- paste0(path.comp, "mod.brm..rds")

estimates <- fread(file.estimates, yaml = TRUE)


fit.dt <-
  CJ(
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("high", "low"),
     area.type = c("treatment"),
     mod.name = c("egp_som25", "match"),
     egp.som.topology = c(NA, "rectangular"),
     egp.cf.nb = c(NA, "sequential"),
     egp.geo.w = c(NA, TRUE),
     match.mod.cov = c(NA, "interact"),
     match.trt.int = c(NA, FALSE),
     sorted = FALSE)

estimates.fit <- 
  estimates[fit.dt,
            on = names(fit.dt),
            nomatch = NULL]

estimates.fit[,
              `:=`(ls.response = factor(ls.response,
                                        levels = c("normal", "tweedie", "binary")),
                   ls.imbalance = factor(ls.imbalance,
                                         levels = c("low", "high")),
                   ls.id = factor(ls.id),
                   name.short = factor(name.short,
                                       levels = c("EGP", "GLM",
                                                  "CEM-ST", "CEM-SC", "CEM-FD",
                                                  "NN-PS-NR", "NN-PS-RE",
                                                  "NN-MA-NR", "NN-MA-RE")))]

ls.type.lev <-
  with(estimates.fit,
       paste0(rep(levels(ls.response), each = 2), "_",
              rep(levels(ls.imbalance), times = 3)))

estimates.fit[,
              ls.type := factor(paste0(ls.response, "_", ls.imbalance),
                                levels = ls.type.lev)]
estimates.fit[,
              ls.uid := ((as.integer(ls.type)-1) * 1000) + as.integer(as.character(ls.id))]



priors <- c(
            prior(student_t(3, 1, 1), class = Intercept),
            prior(student_t(3, 0, 1), class = sd),
            prior(student_t(3, 0, 1), class = Intercept, dpar = nu),
            prior(student_t(3, 0, 1), class = sd, dpar = nu),
            prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
            prior(student_t(3, 0, 1), class = sd, dpar = sigma))

mod.form <-
  bf(mar.std ~
       1 +
       (1 |m| name.short) +
       (1 |r| ls.response:name.short) +
       (1 |i| ls.imbalance:ls.response:name.short),
     nu ~
       1 +
       (1 |m| name.short) +
       (1 |r| ls.response:name.short) +
       (1 |i| ls.imbalance:ls.response:name.short),
     sigma ~
       1 +
       (1 |m| name.short) +
       (1 |r| ls.response:name.short) +
       (1 |i| ls.imbalance:ls.response:name.short))

mod.mar <- brm(mod.form,
               family = student(),
               prior = priors,
               data = estimates.fit,
               chains = 4,
               cores = 4,
               threads = 8,
               warmup = 5000,
               iter = 10000,
               # save_pars = save_pars(all = TRUE),
               init = 0,
               thin = 2,
               control = list(adapt_delta = 0.95),
               refresh = 25,
               empty = FALSE)


saveRDS(mod.mar, file.mod)


# salloc --account=def-cricrime --cpus-per-task=32 --mem=12G --time=0:30:00


