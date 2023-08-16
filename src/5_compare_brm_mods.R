args <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(brms)

source("utilities.R")

# overwrite <- TRUE
mod.id <- as.integer(args[1])
# mod.id <- 2

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.csv")

file.mod <- paste0(path.comp, "mod.brm.", mod.id, ".rds")

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

# TEST
estimates.fit <- estimates.fit[ls.type == "normal_high"]




if(mod.id == 1) {

  priors <- c(
              # prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(student_t(3, 0, 1),  class = sigma))
              # prior(student_t(3, 0 , 1), class = Intercept, dpar = sigma),
              # prior(student_t(3, 0, 1), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~  name.short)

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 5000,
                 iter = 7500,
                 save_pars = save_pars(all = TRUE),
                 # init = 0,
                 # thin = 2,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 2) {

  priors <- c(
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(student_t(3, 0, 1),  class = sigma))
              # prior(student_t(3, 0 , 1), class = Intercept, dpar = sigma),
              # prior(student_t(3, 0, 1), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 5000,
                 iter = 7500,
                 save_pars = save_pars(all = TRUE),
                 # init = 0,
                 # thin = 2,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)

  }


if(mod.id == 3) {

  priors <- c(
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(student_t(3, 0 , 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 | ls.uid),
       sigma ~ name.short)

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 5000,
                 iter = 7500,
                 save_pars = save_pars(all = TRUE),
                 # init = 0,
                 # thin = 2,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 4) {

  priors <- c(
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 | ls.uid),
       sigma ~ name.short + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 4,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 5) {

  priors <- c(
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(student_t(3, 0 , 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 | ls.uid),
       sigma ~ name.short)

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 4,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 6) {

  priors <- c(
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma),
              prior(normal(0, 10), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(normal(0, 10), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 | ls.uid),
       sigma ~ name.short + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 4,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)

  

}


if(mod.id == 7) {

  priors <- c(
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma),
              prior(student_t(3, 1, 1), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(student_t(3, 0, 1), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~ 0 + name.short + (1 | ls.uid),
       sigma ~ 0 + name.short + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 4,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 8) {

  estimates.fit[,
                `:=`(mod.off.mu = 1)]

  priors <- c(
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma),
              prior(student_t(3, 0, 1), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(student_t(3, 0, 1), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~ 0 + offset(mod.off.mu) + name.short + (1 | ls.uid),
       sigma ~ 0 + name.short + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 4,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)

}

if(mod.id == 9) {

  priors <- c(
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma),
              prior(student_t(3, 1, 1), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(student_t(3, 0, 1), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~ 0 + name.short:ls.response:ls.imbalance + (1 | ls.uid),
       sigma ~ 0 + name.short:ls.response:ls.imbalance + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 4,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)


}


if(mod.id == 10) {

  priors <- c(
              prior(normal(0, 1), class = sd),
              prior(normal(0, 1), class = sd, dpar = sigma),
              prior(normal(1, 1), class = Intercept),
              prior(normal(0, 1), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(normal(0, 1), class = Intercept, dpar = sigma),
              prior(normal(0, 1), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 | ls.uid),
       sigma ~ name.short + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 4,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)

}




if(mod.id == 11) {

  priors <- c(
              prior(normal(0, 1), class = sd),
              prior(normal(0, 1), class = sd, dpar = sigma),
              prior(normal(1, 1), class = Intercept),
              prior(normal(0, 1), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(normal(0, 1), class = Intercept, dpar = sigma),
              prior(normal(0, 1), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short*ls.response*ls.imbalance + (1 | ls.uid),
       sigma ~ name.short*ls.response*ls.imbalance + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 4,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)

}

if(mod.id == 12) {

  priors <- c(
              prior(normal(0, 0.5), class = sd),
              prior(normal(0, 0.5), class = sd, dpar = sigma),
              prior(normal(1, 0.5), class = Intercept),
              prior(normal(0, 0.5), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(normal(0, 0.5), class = Intercept, dpar = sigma),
              prior(normal(0, 0.5), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 | ls.uid),
       sigma ~ name.short + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 4,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 13) {

  priors <- c(
              prior(student_t(3, 0, 0.5), class = sd),
              prior(student_t(3, 0, 0.5), class = sd, dpar = sigma),
              prior(student_t(3, 1, 0.5), class = Intercept),
              prior(student_t(3, 0, 0.5), class = b),
              prior(gamma(2, 0.1),  class = nu),
              prior(student_t(3, 0, 0.5), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 0.5), class = b, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short*ls.response*ls.imbalance + (1 | ls.uid),
       sigma ~ name.short*ls.response*ls.imbalance + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 4,
                 # control = list(max_treedepth = 15),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 14) {

  priors <- c(
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 2, 1), class = Intercept, dpar = nu),
              prior(student_t(3, 0, 1), class = b, dpar = nu),
              prior(student_t(3, 0, 1), class = sd, dpar = nu),
              prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 | ls.uid),
       nu ~ name.short + (1 | ls.uid),
       sigma ~ name.short + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 4,
                 # control = list(adapt_delta = 0.7),
                 refresh = 25,
                 empty = FALSE)

}

if(mod.id == 15) {

  priors <- c(
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 2, 1), class = Intercept, dpar = nu),
              prior(student_t(3, 0, 1), class = b, dpar = nu),
              prior(student_t(3, 0, 1), class = sd, dpar = nu),
              prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 | ls.uid),
       nu ~ name.short + (1 | ls.uid),
       sigma ~ name.short + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 2,
                 control = list(adapt_delta = 0.7),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 16) {

  priors <- c(
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 2, 1), class = Intercept, dpar = nu),
              prior(student_t(3, 0, 1), class = b, dpar = nu),
              prior(student_t(3, 0, 1), class = sd, dpar = nu),
              prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 | ls.uid),
       nu ~ name.short + (1 | ls.uid),
       sigma ~ name.short + (1 | ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 8,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 2,
                 control = list(adapt_delta = 0.7),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 17) {

  priors <- c(
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 2, 1), class = Intercept, dpar = nu),
              prior(student_t(3, 0, 1), class = b, dpar = nu),
              prior(student_t(3, 0, 1), class = sd, dpar = nu),
              prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 |l| ls.uid),
       nu ~ name.short + (1 |l| ls.uid),
       sigma ~ name.short + (1 |l| ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 8,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 init = 0,
                 thin = 2,
                 control = list(adapt_delta = 0.7),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 18) {

  priors <- c(
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(student_t(3, 0, 1), class = sd),
              prior(gamma(2, 0.1),  class = nu),
              prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 |l| ls.uid),
       sigma ~ name.short + (1 |l| ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 4,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 2,
                 # control = list(adapt_delta = 0.7),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 19) {

  priors <- c(
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 2, 1), class = Intercept, dpar = nu),
              prior(student_t(3, 0, 1), class = b, dpar = nu),
              prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 |l| ls.uid),
       nu ~ ls.response*ls.imbalance,
       sigma ~ name.short + (1 |l| ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 8,
                 warmup = 10000,
                 iter = 20000,
                 # save_pars = save_pars(all = TRUE),
                 # init = 0,
                 thin = 4,
                 # control = list(adapt_delta = 0.7),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 20) {

  priors <- c(
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 2, 1), class = Intercept, dpar = nu),
              prior(student_t(3, 0, 1), class = b, dpar = nu),
              prior(student_t(3, 0, 1), class = sd, dpar = nu),
              prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short*ls.response*ls.imbalance + (1 |l| ls.uid),
       nu ~ name.short*ls.response*ls.imbalance + (1 |l| ls.uid),
       sigma ~ name.short*ls.response*ls.imbalance + (1 |l| ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 8,
                 cores = 8,
                 threads = 4,
                 warmup = 5000,
                 iter = 10000,
                 # save_pars = save_pars(all = TRUE),
                 init = 0,
                 thin = 2,
                 # control = list(adapt_delta = 0.7),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 21) {

  priors <- c(
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(student_t(3, 0, 1), class = sd),
              prior(student_t(3, 2, 1), class = Intercept, dpar = nu),
              prior(student_t(3, 0, 1), class = b, dpar = nu),
              prior(student_t(3, 0, 1), class = sd, dpar = nu),
              prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short*ls.response*ls.imbalance + (1 |l| ls.uid),
       nu ~ name.short*ls.response*ls.imbalance + (1 |l| ls.uid),
       sigma ~ name.short*ls.response*ls.imbalance + (1 |l| ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 8,
                 warmup = 5000,
                 iter = 25000,
                 # save_pars = save_pars(all = TRUE),
                 init = 0,
                 thin = 4,
                 # control = list(adapt_delta = 0.7),
                 refresh = 25,
                 empty = FALSE)

}



if(mod.id == 22) {

  priors <- c(
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(student_t(3, 0, 1), class = sd),
              prior(gamma(2, 0.1),  class = nu),
              prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma))

  mod.form <-
    bf(mar.std ~ name.short + (1 |l| ls.uid),
       sigma ~ name.short + (1 |l| ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 8,
                 warmup = 5000,
                 iter = 25000,
                 # save_pars = save_pars(all = TRUE),
                 init = 0,
                 thin = 4,
                 # control = list(adapt_delta = 0.7),
                 refresh = 25,
                 empty = FALSE)

}

if(mod.id == 23) {

  priors <- c(
              prior(student_t(3, 1, 1), class = Intercept),
              prior(student_t(3, 0, 1), class = b),
              prior(student_t(3, 0, 1), class = sd),
              prior(gamma(2, 0.1),  class = nu),
              prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
              prior(student_t(3, 0, 1), class = b, dpar = sigma),
              prior(student_t(3, 0, 1), class = sd, dpar = sigma))

  mod.form <-
    bf(mar.std ~ 1 + (1 |m| name.short) + (1 |l| ls.uid),
       sigma ~ 1 + (1 |m| name.short) + (1 |l| ls.uid))

  mod.mar <- brm(mod.form,
                 family = student(),
                 prior = priors,
                 data = estimates.fit,
                 chains = 4,
                 cores = 4,
                 threads = 8,
                 warmup = 5000,
                 iter = 25000,
                 # save_pars = save_pars(all = TRUE),
                 init = 0,
                 thin = 4,
                 # control = list(adapt_delta = 0.7),
                 refresh = 25,
                 empty = FALSE)

}


if(mod.id == 24) {

priors <- c(
            prior(student_t(3, 1, 1), class = Intercept),
            prior(student_t(3, 0, 1), class = b),
            prior(student_t(3, 0, 1), class = sd),
            prior(gamma(2, 0.1),  class = nu),
            prior(student_t(3, 0, 1), class = Intercept, dpar = sigma),
            prior(student_t(3, 0, 1), class = b, dpar = sigma),
            prior(student_t(3, 0, 1), class = sd, dpar = sigma))

mod.form <-
  bf(mar.std ~ 1 + (1 |i| name.short) + (1 |i| ls.uid),
     sigma ~ 1 + (1 |i| name.short) + (1 |i| ls.uid))

mod.mar <- brm(mod.form,
               family = student(),
               prior = priors,
               data = estimates.fit,
               chains = 4,
               cores = 4,
               threads = 8,
               warmup = 5000,
               iter = 25000,
               # save_pars = save_pars(all = TRUE),
               init = 0,
               thin = 4,
               # control = list(adapt_delta = 0.7),
               refresh = 25,
               empty = FALSE)

}

saveRDS(mod.mar, file.mod)


# salloc --account=def-cricrime --cpus-per-task=32 --mem=12G --time=0:30:00


