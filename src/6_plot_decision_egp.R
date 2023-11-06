library(data.table)
library(ggplot2)
library(ggdist)
library(scico)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
path.fig <- paste0(path.results, "figures/")
path.dec <- paste0(path.results, "decision/")

file.egp.dec.all <- paste0(path.dec, "egp.dec.all.rds")
file.egp.dec.imb <- paste0(path.dec, "egp.dec.imb.rds")
file.egp.dec.ls <- paste0(path.dec, "egp.dec.ls.rds")

file.noeff.egp.dec.all <- paste0(path.dec, "noeff.egp.dec.all.rds")
file.noeff.egp.dec.imb <- paste0(path.dec, "noeff.egp.dec.imb.rds")
file.noeff.egp.dec.ls <- paste0(path.dec, "noeff.egp.dec.ls.rds")

egp.dec.all <- readRDS(file.egp.dec.all)
egp.dec.imb <- readRDS(file.egp.dec.imb)
egp.dec.ls <- readRDS(file.egp.dec.ls)

noeff.egp.dec.all <- readRDS(file.noeff.egp.dec.all)
noeff.egp.dec.imb <- readRDS(file.noeff.egp.dec.imb)
noeff.egp.dec.ls <- readRDS(file.noeff.egp.dec.ls)


# PLOT SETUP ########################################################


base.size <- 9 
base.family <- "IBMPlexSansCondensed"

plot_theme <-
  theme_light(base_family = base.family,
              base_size = base.size) +
  theme(
        plot.title = element_text(hjust = 0,
                                  # face = "bold",
                                  face = "plain",
                                  size = rel(1.25),
                                  margin = margin(l = 0, b = base.size, t = base.size/3)),
        plot.tag = element_text(face = "bold"),
        axis.line.x = element_line(color = "black",
                                   linewidth = rel(0.5)),
        axis.line.y = element_line(color = "black",
                                   linewidth = rel(0.5)),
        axis.title.x = element_text(margin = margin(t = base.size/2)),
        axis.title.y = element_text(margin = margin(r = base.size/2)),
        axis.text.y = element_text(color = "black", size = rel(1.2)),
        axis.text.x = element_text(color = "black"),
        axis.ticks = element_line(colour = "grey30"),
        legend.title = element_text(margin = margin(b = base.size/3)),
        legend.position = "right",
        legend.justification = "center",
        legend.key.width = unit(base.size, "pt"),
        legend.key.height = unit(1.75 * base.size, "pt"),
        legend.spacing.y = unit(0.25 * base.size, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(base.size*2, "pt"),
        panel.spacing.y = unit(base.size, "pt"),
        plot.margin = margin(3, 3, 9, 3),
        strip.text = element_text(size = rel(1),
                                  # face = "bold",
                                  hjust = 0.5,
                                  color = "black",
                                  # margin = margin(base.size,
                                  #                 base.size,
                                  #                 base.size,
                                  #                 base.size)),
                                  margin = margin(base.size/2,
                                                  base.size/2,
                                                  base.size/2,
                                                  base.size/2)),
        strip.background = element_rect(fill = "gray93", colour = NA)
        # strip.background = element_rect(fill = NA, colour = NA)
  )



# Labels and colours

ls.imbalance.labs <- c(low = "Moderate imbalance", high = "High imbalance", all = "All")
ls.imbalance.labs <- factor(ls.imbalance.labs, levels = ls.imbalance.labs)
ls.response.labs <- c(normal = "Gaussian", tweedie = "Tweedie", binary = "Binary")
ls.response.labs <- factor(ls.response.labs, levels = ls.response.labs)
type.labs <- 
  c("correct" =  "Correct sign,\nincludes true effect",
    "over" = "Correct sign,\noverestimation",
    "under" = "Correct sign,\nunderestimation",
    "non" = "Non-detection",
    "sign" = "Sign error")
type.labs <- factor(type.labs, levels = type.labs)
type.cols <-
  # c("#4daf4a", "#984ea3", "#377eb8", "#737373", "#e41a1c")
  c("#41ab5d", "#807dba", "#4292c6", "#727272", "#ef3b2c")
names(type.cols) <- type.labs
noeff.type.labs <- 
  c("correct" =  "Correct,\nnon-detection",
    "minor" = "False effect,\nsmall error",
    "moderate" = "False effect,\nmoderate error",
    "severe" = "False effect,\nsevere error")
noeff.type.labs <- factor(noeff.type.labs, levels = noeff.type.labs)
noeff.type.cols <-
  c("#41ab5d", scico(5, palette = "hawai")[3:1])
names(noeff.type.cols) <- noeff.type.labs

plot_dec <- function(x,
                     type = "effect",
                     plot.theme,
                     type.cols,
                     est.lim = NULL,
                     title = waiver(),
                     tag = waiver()){
  if(type == "effect") {
    p.dec <-
      ggplot(x) +
      geom_raster(aes(x = .est, y = .rank, fill = type.lab)) +
      geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.3) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = 1, linewidth = 0.3)
  }
  if(type == "noeffect") {
    p.dec <-
      ggplot(x) +
      geom_raster(aes(x = .est, y = .rank, fill = type.lab)) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = 1, linewidth = 0.3) +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.3)
  }
  p.dec <-
    p.dec +
    scale_fill_manual(values = type.cols, drop = FALSE) +
    guides(fill = guide_legend(byrow = TRUE)) +
    coord_cartesian(xlim = est.lim) +
    labs(title = title,
         tag = tag,
         y = "Landscape identity (sorted)",
         x = "Estimated treatment effect (standardized)",
         fill = NULL) +
    plot_theme
  return(p.dec)
}


egp.dec.all <-
  lapply(egp.dec.all,
         \(x) {
           x[, type.lab := type.labs[as.character(.type)]]})
egp.dec.imb <-
  lapply(egp.dec.imb,
         \(x) {
           x[, ls.imbalance.lab := ls.imbalance.labs[as.character(ls.imbalance)]]
           x[, type.lab := type.labs[as.character(.type)]]})
egp.dec.ls <-
  lapply(egp.dec.ls,
         \(x) {
           x[, ls.imbalance.lab := ls.imbalance.labs[as.character(ls.imbalance)]]
           x[, ls.response.lab := ls.response.labs[as.character(ls.response)]]
           x[, type.lab := type.labs[as.character(.type)]]})
noeff.egp.dec.all <-
  lapply(noeff.egp.dec.all,
         \(x) {
           x[, type.lab := noeff.type.labs[as.character(.type)]]})
noeff.egp.dec.imb <-
  lapply(noeff.egp.dec.imb,
         \(x) {
           x[, ls.imbalance.lab := ls.imbalance.labs[as.character(ls.imbalance)]]
           x[, type.lab := noeff.type.labs[as.character(.type)]]})
noeff.egp.dec.ls <-
  lapply(noeff.egp.dec.ls,
         \(x) {
           x[, ls.imbalance.lab := ls.imbalance.labs[as.character(ls.imbalance)]]
           x[, ls.response.lab := ls.response.labs[as.character(ls.response)]]
           x[, type.lab := noeff.type.labs[as.character(.type)]]})


file.dec.all <- paste0(path.fig, "dec.egp.all.png")
png(file.dec.all, width = 7, height = 3.5, unit = "in", res = 600)
p.dec.all.eff <-
  egp.dec.all$rast[!is.na(.type)] |>
  plot_dec(plot.theme = plot_theme,
           type.cols = type.cols,
           title = "Treatment effect present",
           tag = "(A)")
p.dec.all.noeff <-
  noeff.egp.dec.all$rast[!is.na(.type)] |>
  plot_dec(type = "noeffect",
           plot.theme = plot_theme,
           type.cols = noeff.type.cols,
           title = "No treatment effect",
           tag = "(B)")
p.dec.all.eff+p.dec.all.noeff
dev.off()


file.dec.imb <- paste0(path.fig, "dec.egp.imb.png")
png(file.dec.imb, width = 7, height = 7, unit = "in", res = 600)
p.dec.imb.eff <-
  egp.dec.imb$rast[!is.na(.type)] |>
  plot_dec(plot.theme = plot_theme,
           type.cols = type.cols,
           title = "Treatment effect present",
           tag = "(A)") +
  facet_wrap(vars(ls.imbalance.lab), nrow = 1)
p.dec.imb.noeff <-
  noeff.egp.dec.imb$rast[!is.na(.type)] |>
  plot_dec(type = "noeffect",
           plot.theme = plot_theme,
           type.cols = noeff.type.cols,
           title = "No treatment effect",
           tag = "(B)") +
  facet_wrap(vars(ls.imbalance.lab), nrow = 1)
p.dec.imb.eff/p.dec.imb.noeff
dev.off()


file.dec.ls <- paste0(path.fig, "dec.egp.ls.png")
png(file.dec.ls, width = 7, height = 8.5, unit = "in", res = 600)
p.dec.ls.eff <-
  egp.dec.ls$rast[!is.na(.type)] |>
  plot_dec(plot.theme = plot_theme,
           type.cols = type.cols,
           title = "Treatment effect present",
           tag = "(A)") +
  facet_grid(cols = vars(ls.response.lab),
             rows = vars(ls.imbalance.lab))
p.dec.ls.noeff <-
  noeff.egp.dec.ls$rast[!is.na(.type)] |>
  plot_dec(type = "noeffect",
           plot.theme = plot_theme,
           type.cols = noeff.type.cols,
           title = "No treatment effect",
           tag = "(B)") +
  facet_grid(cols = vars(ls.response.lab),
             rows = vars(ls.imbalance.lab))
p.dec.ls.eff / p.dec.ls.noeff &
  theme(legend.position = "bottom",
        legend.key.height = unit(base.size, "pt"),
        legend.text = element_text(margin = margin(0, base.size, 0, 0,)))
dev.off()

