library(data.table)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.estimates <- paste0(path.comp, "estimates.csv")

file.mod <- paste0(path.comp, "mod.brm.rds")

estimates <- fread(file.estimates)

sub.dt <-
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

estimates.sub <- 
  estimates[sub.dt,
            on = names(sub.dt),
            nomatch = NULL] |>
  merge(names.short, by = c("method", "match.distance", "match.replace", "match.cutpoints"))


estimates.sub[
              # ls.imbalance == "low"
              # ls.response == "tweedie" & ls.imbalance == "low"
              ,
              .(bias = mean(mar.std-1),
                cv = sd(mar.std)/mean(mar.std),
                cv2 = sd(mar.std),
                ser = sum(sign(mar.est) != sign(mar.true)) / .N,
                rmse = rmse(mar.std, 1)),
                by = .(name.short, mod.id,
                       method,
                       match.distance,
                       match.cutpoints,
                       match.mod.cov,
                       match.trt.int,
                       match.replace
                       # )][order(-bias)]
                       )][order(name.short)]



library(brms)
library(bayesplot)

priors <- c(prior(cauchy(0, 1), class = sd),
            prior(student_t(3, 0 , 1), class = b),
            prior(gamma(2, 1),  class = nu),
            prior(student_t(3, 0 , 1), class = b, dpar = sigma))

# ids <- sample(1:1000, 10)
# estimates.fit <- estimates.sub[ls.id %in% ids]
estimates.fit <- estimates.sub

mod.form <-
  bf(mar.std ~ name.short * ls.response * ls.imbalance +
               (1|ls.id:ls.imbalance:ls.response),
     sigma ~ name.short * ls.response * ls.imbalance)

mod.mar <- brm(bf(mar.std ~ name.short * ls.response * ls.imbalance +
                  (1|ls.id:ls.imbalance:ls.response),
                  # (1|ls.id:ls.imbalance:ls.response),
                  sigma ~ name.short * ls.response * ls.imbalance),
               family = student(),
               prior = priors,
               data = estimates.fit,
               chains = 4,
               threads = 8,
               cores = 4,
               warmup = 8000,
               iter = 10000,
               control = list(max_treedepth = 15),
               refresh = 10,
               empty = FALSE)

saveRDS(mod.mar, file.mod)

