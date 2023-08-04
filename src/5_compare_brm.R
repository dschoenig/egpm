library(data.table)
library(brms)

source("utilities.R")

overwrite <- TRUE

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.csv")

file.mod <- paste0(path.comp, "mod.brm.rds")

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

priors <- c(prior(cauchy(0, 1), class = sd),
            prior(student_t(3, 0 , 1), class = b),
            prior(gamma(2, 1),  class = nu),
            prior(student_t(3, 0 , 1), class = b, dpar = sigma))

mod.form <-
  bf(mar.std ~ name.short * ls.response * ls.imbalance +
               (1|ls.id:ls.imbalance:ls.response),
     sigma ~ name.short * ls.response * ls.imbalance)

if(file.exists(file.mod) & overwrite == FALSE) {
  mod.mar <- readRDS(file.mod)
} else {
  mod.mar <- brm(bf(mar.std ~ name.short * ls.response * ls.imbalance +
                    (1|ls.id:ls.imbalance:ls.response),
                    sigma ~ name.short * ls.response * ls.imbalance),
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 8,
                 warmup = 8000,
                 iter = 10000,
                 control = list(max_treedepth = 15),
                 refresh = 10,
                 empty = FALSE)
  saveRDS(mod.mar, file.mod)
}

# salloc --account=def-cricrime --cpus-per-task=32 --mem=12G --time=0:30:00

