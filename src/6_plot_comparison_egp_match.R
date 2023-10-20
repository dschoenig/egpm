library(data.table)
library(ggplot2)
library(ggdist)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
path.fig <- paste0(path.results, "figures/")
file.estimates <- paste0(path.comp, "estimates.rds")
file.egp.match.boot <- paste0(path.comp, "egp.match.boot.rds")

estimates <- readRDS(file.estimates)
egp.match.boot <- readRDS(file.egp.match.boot)


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
        legend.position = "right",
        legend.justification = "center",
        legend.key.size = unit(base.size, "pt"),
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


# Function to plot comparisons

plot_ls_comp <- function(x, plot.theme, method.cols) {

  ref.line <-
    list(
      x[,
        .(type = "Theoretical maximum",
          xint = 0,
          ymin = -Inf,
          ymax = Inf),
         # ymax = 1 + length(unique(x$name.short))),
         by = c("ls.imbalance.lab", ".stat.lab")],
      x[name.short == "EGP",
        .(type = "EGP",
          xint = mean,
          ymin = -Inf,
          ymax = 0.5 + length(unique(x$name.short))),
         by = c("ls.imbalance.lab", ".stat.lab")],
      x[name.short == "GLM",
        .(type = "GLM",
          xint = mean,
          ymin = -Inf,
          ymax = 0.5 + length(unique(x$name.short))),
         by = c("ls.imbalance.lab", ".stat.lab")]) |>
    rbindlist()

  ref.ind <-
    rbind(
          CJ(x = 0,
             y = length(unique(x$name.short))+1.25,
             side = "left",
             full = TRUE,
             hjust = 1,
             label =  "better ▶ ",
             .stat.lab = unique(x$.stat.lab)),
          CJ(x = 0,
             y = length(unique(x$name.short))+1.25,
             side = "right",
             full = TRUE,
             hjust = 0,
             label =  " ◀ better",
             .stat.lab = unique(x$.stat.lab)),
          CJ(x = 0,
             y = length(unique(x$name.short))+1.25,
             side = "left",
             full = FALSE,
             hjust = 1,
             label =  "▶ ",
             .stat.lab = unique(x$.stat.lab)),
          CJ(x = 0,
             y = length(unique(x$name.short))+1.25,
             side = "right",
             full = FALSE,
             hjust = 0,
             label =  " ◀",
             .stat.lab = unique(x$.stat.lab)))

  ref.ind.sel <-
    x[,
      .(min = min(mean),
        max = max(mean),
        closest.left = ifelse(any(mean < 0),
                              max(mean[mean < 0]),
                              0),
        closest.right = ifelse(any(mean > 0),
                               min(mean[mean > 0]),
                               0),
        lim.left = min(c(mean-se, 0)),
        lim.right = max(c(mean+se, 0))),
      by = c(".stat.lab")] |>
    _[,
      .(.stat.lab,
        left = min < 0,
        right = max > 0,
        left.full = abs(lim.left) > 0.3 * abs(lim.right - lim.left),
        right.full = abs(lim.right) > 0.3 * abs(lim.right - lim.left))] |>
    melt(measure.vars = c("left", "right"),
         variable.name = "side",
         value.name = "side.sel") |>
    _[side.sel == TRUE,
      .(.stat.lab,
        side,
        full = fcase(side == "left" & left.full == TRUE, TRUE,
                     side == "right" & right.full == TRUE, TRUE,
                     default = FALSE))]

  ref.ind <- ref.ind[ref.ind.sel, on = names(ref.ind.sel)]

  ggplot(x) +
    geom_linerange(data = ref.line,
                   mapping = aes(x = xint, ymin = ymin, ymax = ymax, linetype = type),
                   linewidth = 0.3) +
    geom_text(data = ref.ind,
              mapping = aes(x = x, y = y, label = label, hjust = hjust),
              vjust = 0.5,
              family = base.family,
              # fontface = "italic",
              size = base.size/4,
              colour = "grey30") +
    geom_pointrange(aes(x = mean, y = name.short,
                        xmin = mean - se,
                        xmax = mean + se,
                        colour = method.lab),
                        # xmin = ci_l,
                        # xmax = ci_u),
                    size = 2/base.size,
                    fill = NA,
                    linewidth = 5/base.size) +
   # geom_linerange(aes(y = name.short,
   #                    xmin = mean - se,
   #                    xmax = mean + se),
   #                    # xmin = ci_l,
   #                    # xmax = ci_u),
   #                linewidth = 0.1) +
    scale_y_discrete(limits = rev, expand = expansion(add = c(0.5, 1.75))) +
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    scale_linetype_manual(values = c("Theoretical maximum" = "dotted",
                                     "EGP" = "solid",
                                     "GLM" = "dashed"),
                          drop = FALSE) +
    scale_colour_manual(values = c(method.cols,
                                   GLM = unname(method.cols[1]),
                                   EGP = unname(method.cols[2])),
                          drop = FALSE) +
    guides(linetype = "none", colour = "none") +
    labs(colour = NULL, x = "Performance", y = NULL) +
    facet_grid(rows = vars(ls.imbalance.lab), cols = vars(.stat.lab),
               scales = "free_x") +
    plot.theme

}

# Labels and colours

stat.labs <- c(bias = "Bias", rmse = "RMSE", ser = "SER")
stat.labs <- factor(stat.labs, levels = stat.labs)
ls.imbalance.labs <- c(low = "Moderate imbalance", high = "High imbalance", all = "All")
ls.imbalance.labs <- factor(ls.imbalance.labs, levels = ls.imbalance.labs)
ls.response.labs <- c(normal = "Gaussian", tweedie = "Tweedie", binary = "Binary")
ls.response.labs <- factor(ls.response.labs, levels = ls.response.labs)
method.labs = c("Generalized linear model",
                "Embedded gaussian process",
                "Coarsened exact matching",
                "Nearest neighbour matching")
method.cols <- c("grey50", "#984ea3", "#377eb8", "#4daf4a")
names(method.cols) <- method.labs
# method.cols.light <- c("grey80", "#984ea399", "#377eb899", "#4daf4a99")
# method.cols.light <- c("grey80", "#984ea38c", "#377eb88c", "#4daf4a8c")
method.cols.light <- c("grey80", "#984ea380", "#377eb880", "#4daf4a80")
names(method.cols.light) <- method.labs



# DATA PREPARATION ##################################################


# Prepare data for distribution comparison

sub.dt <-
  CJ(
     ls.response = c("normal", "tweedie", "binary"),
     # ls.response = c("binary"),
     ls.imbalance = c("low", "high"),
     area.type = "treatment",
     mod.name = c("egp_som25", "match"),
     egp.som.topology = c(NA, "rectangular"),
     egp.cf.nb = c(NA, "sequential"),
     egp.geo.w = c(NA, TRUE),
     match.mod.cov = factor(c(NA, "interact")),
     match.trt.int = c(NA, FALSE),
     sorted = FALSE)

estimates.sub.p <- subset_estimates(estimates, sub.dt)

estimates.sub.p[, ls.imbalance.lab := ls.imbalance.labs[as.character(ls.imbalance)]]
estimates.sub.p[, ls.response.lab := ls.response.labs[as.character(ls.response)]]
estimates.sub.p[, name.short := relevel(name.short, "GLM")]
estimates.sub.p[name.short %like% "NN-", method.lab := method.labs[4]]
estimates.sub.p[name.short %like% "CEM-", method.lab := method.labs[3]]
estimates.sub.p[name.short == "EGP", method.lab := method.labs[2]]
estimates.sub.p[name.short == "GLM", method.lab := method.labs[1]]
estimates.sub.p[, method.lab := factor(method.lab, levels = method.labs)]


# Prepare data for summary comparison

egp.match.l <- comparisons_cast(egp.match.boot)
egp.match.l[, .stat.lab := stat.labs[.stat]]
egp.match.l[, ls.imbalance.lab := ls.imbalance.labs[as.character(ls.imbalance)]]
egp.match.l[, name.short := relevel(name.short, "GLM")]

method.labs = c("Generalized linear model",
                "Embedded gaussian process",
                "Coarsened exact matching",
                "Nearest neighbour matching")
egp.match.l[name.short %like% "NN-", method.lab := method.labs[4]]
egp.match.l[name.short %like% "CEM-", method.lab := method.labs[3]]
egp.match.l[name.short == "EGP", method.lab := method.labs[2]]
egp.match.l[name.short == "GLM", method.lab := method.labs[1]]
egp.match.l[, method.lab := factor(method.lab, levels = method.labs)]


# PLOT DISTRIBUTIONS ################################################

file.dist <- paste0(path.fig, "comp.egp_match.dist.all.png")
png(file.dist, width = 5, height = 4, unit = "in", res = 600)
ggplot(estimates.sub.p) +
  stat_halfeye(aes(y = name.short, x = mar.std, fill = method.lab),
               point_interval = "mean_qi",
               # normalize = "panels",
               scale = 0.9,
               interval_size_range = c(0.35, 0.9),
               fatten_point = 1.35,
               p_limits = c(0.001, 0.999),
               n = 501,
               .width = c(0.5, 0.95)) +
  geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", colour = 2, linewidth = 0.5) +
  scale_fill_manual(values = method.cols.light,
                        drop = FALSE) +
  guides(linetype = "none", fill = "none") +
  scale_y_discrete(limits = rev) +
  coord_cartesian(xlim = c(-0.5, 2.5)) +
  labs(x = "Estimated average treatment effect (ATT)",
       y = NULL) +
  facet_wrap(vars(ls.imbalance.lab), nrow = 1) +
  plot_theme
dev.off()


file.dist <- paste0(path.fig, "comp.egp_match.dist.ls.png")
png(file.dist, width = 7, height = 7, unit = "in", res = 600)
ggplot(estimates.sub.p) +
  stat_halfeye(aes(y = name.short, x = mar.std, fill = method.lab),
               point_interval = "mean_qi",
               # normalize = "panels",
               scale = 0.9,
               interval_size_range = c(0.35, 0.9),
               fatten_point = 1.35,
               p_limits = c(0.001, 0.999),
               n = 501,
               .width = c(0.5, 0.95)) +
  geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", colour = 2, linewidth = 0.5) +
  scale_fill_manual(values = method.cols.light,
                        drop = FALSE) +
  guides(linetype = "none", fill = "none") +
  scale_y_discrete(limits = rev) +
  coord_cartesian(xlim = c(-0.5, 2.5)) +
  labs(x = "Estimated average treatment effect (ATT)",
       y = NULL) +
  facet_grid(rows = vars(ls.imbalance.lab),
             cols = vars(ls.response.lab),
             scales = "free_x") +
  plot_theme
dev.off()



# PLOT SUMMARIES ###################################################

egp.match.all <-
  egp.match.l[ls.response == "all" &
              ls.imbalance == "all" &
              is.na(.comp)]

egp.match.all[ls.imbalance == "all",
              ls.imbalance.lab := "Moderate or high imbalance (N=6000)"]

paste0(path.fig, "comp.egp_match.sum.all_all.png") |>
png(width = 7, height = 3, unit = "in", res = 600)
egp.match.all |>
plot_ls_comp(plot.theme = plot_theme, method.cols = method.cols) |>
print()
dev.off()


ls.sum <- levels(egp.match.l$ls.response)

for(i in seq_along(ls.sum)) {
  file.ls.sum <- paste0(path.fig, "comp.egp_match.sum.", ls.sum[i], ".png")
  png(file.ls.sum, width = 7, height = 5, unit = "in", res = 600)
  egp.match.l[ls.response == ls.sum[i] &
              ls.imbalance != "all" &
              is.na(.comp)] |>
  plot_ls_comp(plot.theme = plot_theme, method.cols = method.cols) |>
  print()
  dev.off()
}

