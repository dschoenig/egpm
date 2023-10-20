library(data.table)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
path.dec <- paste0(path.results, "decision/")
file.estimates <- paste0(path.comp, "estimates.rds")


file.egp.dec_95.all <- paste0(path.dec, "egp.dec_95.all.rds")
file.egp.dec_95.ls <- paste0(path.dec, "egp.dec_95.ls.rds")
file.egp.dec_95.all.sum <- paste0(path.dec, "egp.dec_95.all.sum.csv")
file.egp.dec_95.ls.sum <- paste0(path.dec, "egp.dec_95.ls.sum.csv")


estimates <- readRDS(file.estimates)

sub.dt <-
  CJ(
     ls.response = c("normal", "tweedie", "binary"),
     ls.imbalance = c("low", "high"),
     area.type = "treatment",
     mod.name = c("egp_som25"),
     egp.som.topology = c(NA, "rectangular"),
     egp.cf.nb = c(NA, "sequential"),
     egp.geo.w = c(NA, TRUE),
     match.mod.cov = factor(c(NA, "interact")),
     match.trt.int = c(NA, FALSE),
     sorted = FALSE)

estimates.sub <- subset_estimates(estimates, sub.dt)

egp.dec_95.all <-
  classify_estimates(estimates.sub,
                     ls.id = "ls.uid",
                     est.ci = c("mar.q2.5", "mar.q97.5"),
                     group.vars = "ls.imbalance")

saveRDS(egp.dec_95.all, file.egp.dec_95.all)

egp.dec_95.all.sum <-
  egp.dec_95.all$ls[order(ls.imbalance, .type),
                    .(count = .N),
                    by = c("ls.imbalance", ".type")]

egp.dec_95.all.sum[, freq := count/sum(count), by = "ls.imbalance"]


egp.dec_95.ls <-
  classify_estimates(estimates.sub,
                     ls.id = "ls.id",
                     est.ci = c("mar.q2.5", "mar.q97.5"),
                     group.vars = c("ls.response", "ls.imbalance"))

saveRDS(egp.dec_95.ls, file.egp.dec_95.ls)

egp.dec_95.ls.sum <-
  egp.dec_95.ls$ls[order(ls.imbalance, .type),
                   .(count = .N),
                   by = c("ls.response", "ls.imbalance", ".type")]

egp.dec_95.ls.sum[, freq := count/sum(count), by = c("ls.response", "ls.imbalance")]


# ggplot(xc[!is.na(.type) & ls.response == "normal" & ls.imbalance == "low"]) +
#   geom_raster(aes(x = .est, y = .rank, fill = type.int))



# # classify_estimates(estimates.sub, ls.id = "ls.uid", group.vars = "ls.imbalance")$rast |>
# classify_estimates(estimates.sub, ls.id = "ls.id", group.vars = c("ls.response", "ls.imbalance"))$rast |>
# # classify_estimates(estimates.sub, ls.id = "ls.uid", group.vars = NULL)$rast |>
# _[!is.na(.type)] |>


# classify_estimates(estimates.sub,
#                    ls.id = "ls.uid",
#                    # est.ci = c("mar.q2.5", "mar.q97.5"),
#                    est.ci = c("mar.q0.5", "mar.q99.5"),
#                    # est.ci = c("mar.q25", "mar.q75"),
#                    group.vars = "ls.imbalance") |>
# _$rast[!is.na(.type)] |>
# ggplot() +
#   geom_raster(aes(x = .est, y = .rank, fill = .type)) +
#   # geom_point(data = xc[p.within == TRUE],
#   #            mapping = aes(x = .p, y = .rank),
#   #            shape = ".",
#   #            colour = "grey20") +
#   scale_fill_manual(values = type.cols) +
#   # coord_cartesian(xlim = c(-0.75, 2.75)) +
#   # facet_grid(rows = vars(ls.imbalance),
#   #            cols = vars(ls.response),
#   #            scales = "free_x") +
#   facet_wrap(vars(ls.imbalance), nrow = 1) +
#   plot_theme




# base.size <- 9 
# base.family <- "IBMPlexSansCondensed"

# plot_theme <-
#   theme_light(base_family = base.family,
#               base_size = base.size) +
#   theme(
#         plot.title = element_text(hjust = 0,
#                                   face = "bold",
#                                   margin = margin(l = 0, b = base.size/3, t = base.size/3)),
#         plot.tag = element_text(face = "bold"),
#         axis.line.x = element_line(color = "black",
#                                    linewidth = rel(0.5)),
#         axis.line.y = element_line(color = "black",
#                                    linewidth = rel(0.5)),
#         axis.title.x = element_text(margin = margin(t = base.size/2)),
#         axis.title.y = element_text(margin = margin(r = base.size/2)),
#         axis.text.y = element_text(color = "black", size = rel(1.2)),
#         axis.text.x = element_text(color = "black"),
#         axis.ticks = element_line(colour = "grey30"),
#         legend.title = element_text(margin = margin(b = base.size/3)),
#         legend.position = "right",
#         legend.justification = "center",
#         legend.key.size = unit(base.size, "pt"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.spacing.x = unit(base.size*2, "pt"),
#         panel.spacing.y = unit(base.size*3, "pt"),
#         plot.margin = margin(3, 3, 3, 3),
#         strip.text = element_text(size = rel(1),
#                                   # face = "bold",
#                                   hjust = 0.5,
#                                   color = "black",
#                                   # margin = margin(base.size,
#                                   #                 base.size,
#                                   #                 base.size,
#                                   #                 base.size)),
#                                   margin = margin(base.size/2,
#                                                   base.size/2,
#                                                   base.size/2,
#                                                   base.size/2)),
#         strip.background = element_rect(fill = "gray93", colour = NA)
#         # strip.background = element_rect(fill = NA, colour = NA)
#   )




# type.cols <-
#   c(over = "#807dba",
#     correct = "#41ab5d",
#     under = "#4292c6",
#     non = "#737373",
#     sign = "#ef3b2c")


# if(is.null(int.labels)) {
#   type.int.labs <- 
#     c("correct" =  "Correct"
#       "over" = "Correct sign,\noverestimation",
#       "under" = "Correct sign,\nunderestimation",
#       "non" = "Non-detection",
#       "sign" = "Sign error")
# } else {
#   type.int.labs <- int.labels
# }
