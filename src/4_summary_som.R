library(data.table)
library(stringi)

source("utilities.R")

## Paths and parameters

models <- c("egp_som25",
            "egp_som10",
            "egp_som50",
            "egp_som25_unequal")

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
file.som_eval <- paste0(path.comp, "som_eval.rds")

## Get model info

som.eval.i <- list()

for(i in 1:nrow(landscapes)) {

  files.all <- list.files(paste0(path.results, landscapes$name[i]))
  som.sel <- models[paste0(models, "_embedding.rds") %in% files.all]
  files.load <- 
    paste0(path.results,
           landscapes$name[i], "/",
           paste0(som.sel, "_embedding.rds"))

  som.eval.j <- list()

  for(j in seq_along(som.sel)) {

    res <- readRDS(files.load[j])

    res[,
        `:=`(trt.effect = landscapes$trt.effect[i],
             ls.name = landscapes$name[i],
             mod.name = som.sel[j])]

    som.eval.j[[j]] <- res

  }

  som.eval.i[[i]] <- rbindlist(som.eval.j, fill = TRUE)

}

som.eval <- rbindlist(som.eval.i)

som.eval[, ls.id := factor(ls.id)]

som.cols <-
  c("quant.error", "topo.error", "var.expl", "units.ov_frac", "units.js_div")
som.newcols <- paste0("som.", som.cols)

setnames(som.eval, som.cols, som.newcols)
setnames(som.eval, c("landscape"), c("ls.id"))
setcolorder(som.eval, c("trt.effect", "ls.id", "ls.name", "mod.name", "mod.id"))

saveRDS(som.eval, file.som_eval)


file.estimates.rds <- paste0(path.comp, "estimates.rds")

estimates <- readRDS(file.estimates.rds)

est.som <- merge(estimates, som.eval)
est.som.sub <- est.som[trt.effect == TRUE & area.type == "treatment" & mod.id == 8 & mod.name == "egp_som25"]
est.som.sub[, error := abs(mar.std -1)]
est.som.sub[, ls_name := factor(ls.name)]


est.som.sub |>
  ggplot(aes(x = som.var.expl, y = error)) +
  geom_point() +
  facet_grid(cols = vars(ls.response), rows = vars(ls.imbalance))

mod.var.expl <-
  gam(error ~ ls.name + s(som.var.expl, by = ls_name), data = est.som.sub, method = "REML", select = TRUE)
summary(mod.var.expl)
plot(mod.var.expl, pages = 1)


est.som.sub |>
  ggplot(aes(x = som.units.ov_frac, y = error)) +
  geom_point() +
  facet_grid(cols = vars(ls.response), rows = vars(ls.imbalance))

mod.units.ov_frac <-
  gam(error ~ ls.name + s(som.units.ov_frac, by = ls_name), data = est.som.sub, method = "REML", select = TRUE)
summary(mod.units.ov_frac)
plot(mod.units.ov_frac, pages = 1)


est.som.sub |>
  ggplot(aes(x = som.units.js_div, y = error)) +
  geom_point() +
  facet_grid(cols = vars(ls.response), rows = vars(ls.imbalance))

mod.units.js_div <-
  gam(error ~ ls.name + s(som.units.js_div, by = ls_name), data = est.som.sub, method = "REML", select = TRUE)
summary(mod.units.js_div)
plot(mod.units.js_div, pages = 1)

est.som.sub |>
ggplot(aes(x = som.topo.error, y = error)) +
geom_point() +
facet_grid(cols = vars(ls.response), rows = vars(ls.imbalance))

mod.topo.error <-
gam(error ~ ls.name + s(som.topo.error, by = ls_name), data = est.som.sub, method = "REML", select = TRUE)
summary(mod.topo.error)
plot(mod.topo.error, pages = 1)

est.som.sub |>
ggplot(aes(x = som.quant.error, y = error)) +
geom_point() +
facet_grid(cols = vars(ls.response), rows = vars(ls.imbalance))

mod.quant.error <-
gam(error ~ ls.name + s(som.quant.error, by = ls_name), data = est.som.sub, method = "REML", select = TRUE)
summary(mod.quant.error)
plot(mod.quant.error, pages = 1)
