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
mod.type <- "egp_sam0.0"

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
path.results <- paste0(path.base, "results/", ls.type, "/")

# files.res <- paste0(path.results, list.files(path.results, pattern = mod.type))
# som.sizes <- c(10, 100, 25, "50_200", "50_2500", "50_500")

# effects <- list()

# for(i in seq_along(files.res)){
#   eff.egp <- readRDS(files.res[i])
#   effects[[i]] <- 
#     rbind(eff.egp$noint$effects[, .(id, interactions = FALSE, method = "egp", som = som.sizes[i],
#                                     mean, q2.5, q97.5)],
#           eff.egp$int$effects[, .(id, interactions = TRUE, method = "egp", som = som.sizes[i], 
#                                   mean, q2.5, q97.5)])
# }

# effects <- rbindlist(effects)
# effects[, som := factor(som)]



files.res <- paste0(path.results, "egp_sam0.0", c("05", "1", "2"), c("_som50.rds"))
sam.sizes <- c(0.005, 0.01, 0.02)
som.sizes <- c(50, 50, 50)

effects <- list()

for(i in seq_along(files.res)){
  eff.egp <- readRDS(files.res[i])
  effects[[i]] <- 
    rbind(eff.egp$noint$effects[, .(id, interactions = FALSE, method = "egp",
                                    sam = sam.sizes[i], som = som.sizes[i],
                                    mean, q2.5, q97.5)],
          eff.egp$int$effects[, .(id, interactions = TRUE, method = "egp",
                                  sam = sam.sizes[i], som = som.sizes[i],
                                  mean, q2.5, q97.5)])
}

effects <- rbindlist(effects)
effects[, `:=`(som = factor(som),
               sam = factor(sam))]

effects[interactions == TRUE, mean(mean), c("som", "sam")]

ggplot(effects[interactions == FALSE]) +
# geom_violin(aes(x = method, y = mean))
stat_halfeye(aes(y = sam, x = mean, fill_ramp = stat(cut_cdf_qi(cdf))),
             alpha = 0.8,
             point_interval = mean_qi,
             n = 1001) +
geom_vline(xintercept = 0, alpha = 0.8) +
geom_vline(xintercept = 1, colour = 2, alpha = 0.8) +
geom_vline(xintercept = effects[interactions == TRUE & method == "egp" & som == 50, mean(mean)],
           colour = 1, linetype = "dashed", alpha = 0.8) +
# stat_gradientinterval(aes(x = method, y = mean)) +
# stat_dots(aes(y = method, x = mean)) +
scale_fill_ramp_discrete(na.translate = FALSE, from = "grey85", range = c(1,0)) +
scale_fill_discrete_qualitative("Set2") +
coord_cartesian(xlim = c(-2, 3)) +
theme_ggdist()

