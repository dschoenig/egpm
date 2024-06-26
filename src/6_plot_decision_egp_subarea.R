library(data.table)
library(ggplot2)
library(ggdist)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
path.fig <- paste0(path.results, "figures/")
path.dec <- paste0(path.results, "decision/")

file.egp_subarea.dec_95.all <- paste0(path.dec, "egp_subarea.dec_95.all.rds")
file.egp_subarea.dec_95.ls <- paste0(path.dec, "egp_subarea.dec_95.ls.rds")

egp_subarea.dec_95.all <- readRDS(file.egp_subarea.dec_95.all)
egp_subarea.dec_95.ls <- readRDS(file.egp_subarea.dec_95.ls)



# PLOT SETUP ########################################################


base.size <- 9 
base.family <- "IBMPlexSansCondensed"

plot_theme <-
  theme_light(base_family = base.family,
              base_size = base.size) +
  theme(
        plot.title = element_text(hjust = 0,
                                  face = "bold",
                                  margin = margin(l = 0, b = base.size/3, t = base.size/3)),
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
        legend.position = "bottom",
        legend.justification = "center",
        legend.key.size = unit(base.size, "pt"),
        legend.text = element_text(margin = margin(0,
                                                   base.size,
                                                   0,
                                                   0,)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(base.size*2, "pt"),
        panel.spacing.y = unit(base.size*3, "pt"),
        plot.margin = margin(3, 3, 3, 3),
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
type.cols <-
  c("#41ab5d", "#807dba", "#4292c6", "#737373", "#ef3b2c")
names(type.cols) <- type.labs



egp_subarea.dec_95.all <-
  lapply(egp_subarea.dec_95.all,
         \(x) {
           x[, ls.imbalance.lab := ls.imbalance.labs[as.character(ls.imbalance)]]
           x[, type.lab := type.labs[as.character(.type)]]})

egp_subarea.dec_95.ls <-
  lapply(egp_subarea.dec_95.ls,
         \(x) {
           x[, ls.imbalance.lab := ls.imbalance.labs[as.character(ls.imbalance)]]
           x[, ls.response.lab := ls.response.labs[as.character(ls.response)]]
           x[, type.lab := type.labs[as.character(.type)]]})




file.dec_95.all <- paste0(path.fig, "dec_95.egp_subarea.all.png")
png(file.dec_95.all, width = 7, height = 5, unit = "in", res = 600)
xlims <-
  egp_subarea.dec_95.all$rast[!is.na(.type), quantile(c(mar.q2.5, mar.q97.5), c(0.01, 0.99))]
egp_subarea.dec_95.all$rast[!is.na(.type)] |>
ggplot() +
  geom_raster(aes(x = .est, y = .rank, fill = type.lab)) +
  geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 1, linewidth = 0.3) +
  scale_fill_manual(values = type.cols) +
  labs(y = "Landscape identity (sorted)",
       x = "Estimated average treatment effect (ATT)",
       fill = NULL) +
  facet_wrap(vars(ls.imbalance.lab), nrow = 1) +
  coord_cartesian(xlim = xlims) +
  plot_theme
dev.off()


file.dec_95.ls <- paste0(path.fig, "dec_95.egp_subarea.ls.png")
png(file.dec_95.ls, width = 7, height = 5, unit = "in", res = 600)
xlims <-
  egp_subarea.dec_95.ls$rast[!is.na(.type), quantile(c(mar.q2.5, mar.q97.5), c(0.01, 0.99))]
egp_subarea.dec_95.ls$rast[!is.na(.type)] |>
ggplot() +
  geom_raster(aes(x = .est, y = .rank, fill = type.lab)) +
  geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = 1, linewidth = 0.3) +
  scale_fill_manual(values = type.cols) +
  labs(y = "Landscape identity (sorted)",
       x = "Estimated average treatment effect (ATT)",
       fill = NULL) +
  facet_grid(
             rows = vars(ls.imbalance.lab),
             cols = vars(ls.response.lab),
             scales = "free_x") +
  coord_cartesian(xlim = xlims) +
  plot_theme
dev.off()


