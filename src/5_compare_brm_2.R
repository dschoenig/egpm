library(data.table)
library(brms)

source("utilities.R")

overwrite <- TRUE

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.csv")

file.mod <- paste0(path.comp, "mod.brm_2.rds")

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


# Test section #################

# priors <- c(
#             prior(cauchy(0, 1), class = sd),
#             # prior(student_t(3, 1 , 1), class = Intercept),
#             prior(student_t(3, 1 , 1), class = b),
#             prior(gamma(2, 1),  class = nu),
#             prior(student_t(3, 0, 1), class = b, dpar = sigma))

# ls.sam <- sample(1:1000, 10)
# estimates.fit2 <- estimates.fit[ls.id %in% ls.sam]

# mod.form <-
#   bf(mar.std ~ 0 + name.short + name.short:ls.type + (1|ls.uid),
#      sigma ~ 0 + name.short + name.short:ls.type)

# mod.mar <- brm(mod.form,
#                family = student(),
#                prior = priors,
#                data = estimates.fit2,
#                chains = 4,
#                cores = 4,
#                warmup = 4000,
#                iter = 5000,
#                control = list(max_treedepth = 15),
#                refresh = 10,
#                empty = FALSE)


# priors <- c(prior(cauchy(0, 1), class = sd),
#             prior(student_t(3, 1 , 1), class = Intercept),
#             prior(student_t(3, 0 , 1), class = b),
#             prior(gamma(2, 1),  class = nu),
#             prior(student_t(3, 0 , 1), class = b, dpar = sigma))

# estimates.fit2 <- estimates.fit

# mod.mar <- brm(bf(mar.std ~ name.short + (1|ls.id:ls.imbalance:ls.response),
#                   sigma ~ name.short),
#                family = student(),
#                prior = priors,
#                data = estimates.fit2,
#                chains = 4,
#                cores = 4,
#                warmup = 4000,
#                iter = 5000,
#                control = list(max_treedepth = 12),
#                refresh = 10,
#                empty = FALSE)

# End test section #####################


# priors <- c(prior(cauchy(0, 1), class = sd),
#             prior(student_t(3, 0 , 1), class = b),
#             prior(gamma(2, 1),  class = nu),
#             prior(student_t(3, 0 , 1), class = b, dpar = sigma))


priors <- c(
            prior(student_t(3, 0, 1), class = sd),
            prior(student_t(3, 1 , 1), class = Intercept),
            prior(student_t(3, 0 , 1), class = b),
            prior(gamma(2, 1),  class = nu),
            prior(student_t(3, 0 , 1), class = Intercept, dpar = sigma),
            prior(student_t(3, 0, 1), class = b, dpar = sigma))

# priors <- c(
#             prior(student_t(3, 0, 5), class = sd),
#             prior(normal(1 , 2), class = Intercept),
#             prior(normal(0 , 2), class = b),
#             # prior(gamma(4, 1),  class = nu),
#             prior(normal(0, 2), class = b, dpar = sigma))


mod.form <-
  bf(mar.std ~ name.short + (1|ls.uid),
     sigma ~ name.short)

# mod.form <-
#   bf(mar.std ~ name.short + name.short:ls.type
#                (1|ls.id:ls.imbalance:ls.response),
#      # nu = 4,
#      sigma ~ name.short * ls.response * ls.imbalance)


if(file.exists(file.mod) & overwrite == FALSE) {
  mod.mar <- readRDS(file.mod)
} else {

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

  print(summary(mod.mar))

  saveRDS(mod.mar, file.mod)
}

# salloc --account=def-cricrime --cpus-per-task=32 --mem=12G --time=0:30:00

