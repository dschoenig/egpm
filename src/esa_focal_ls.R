library(data.table)
library(mgcv)
library(kohonen)
library(ggdist)

source("utilities_fcne.R")
source("utilities.R")

draw_quad <- function(x.dim, y.dim, increment = 1){
  full.grid <- expand.grid(x = 1:x.dim,
                           y = 1:y.dim)
  full.grid$id <- 1:nrow(full.grid)
  starting.points <-
    expand.grid(x = seq(1, x.dim-increment, increment),
                y = seq(1, y.dim-increment, increment))
  quads <- matrix(NA, nrow = 4, ncol = nrow(starting.points))
  for(i in 1:nrow(starting.points)){
    quads[,i] <-  
      c(full.grid$id[which(full.grid$x == starting.points$x[i] &
                           full.grid$y == starting.points$y[i])],
        full.grid$id[which(full.grid$x == starting.points$x[i] + increment &
                           full.grid$y == starting.points$y[i])],
        full.grid$id[which(full.grid$x == starting.points$x[i] + increment &
                           full.grid$y == starting.points$y[i] + increment)],
        full.grid$id[which(full.grid$x == starting.points$x[i] &
                           full.grid$y == starting.points$y[i] + increment)])
  }
  return(quads)
}



ls.type <- "1000_4cov_nl"
egp.type <- "egp_sam0.01_som25"

path.base <- "../"
path.ls <- paste0(path.base, "landscapes/", ls.type, "/")
path.ls.data <- paste0(path.ls, "data/")
path.mod <- paste0(path.base, "models/", ls.type, "/")
path.results <- paste0(path.base, "results/esa/")
if(!dir.exists(path.results)) dir.create(path.results, recursive = TRUE)

ls.id <- 249

file.ls <- paste0(path.ls.data,
                  stri_pad_left(ls.id, 4, 0), ".rds")
file.egp <- paste0(path.mod, egp.type, "_",
                   stri_pad_left(ls.id, 4, 0), ".rds")

ls <- readRDS(file.ls)
mod.egp <- readRDS(file.egp)

str(mod.egp, max.level = 1)

ls.sam <- mod.egp$sample
som.fit <- mod.egp$som


## ILUSTRATE SOM FIT ###########################################################


library(rgl)
  
som.dim = 25

pts <- scale(as.matrix(ls.sam[1:1000,.(z1, z2, z3)]),
             center = TRUE, scale = TRUE)
quads <- draw_quad(som.dim, som.dim, 1)

som.init <- mesh3d(init_som(pts, som.dim, som.dim), quads = quads)


grid <- somgrid(xdim = som.dim, ydim = som.dim, 
                topo = "rectangular", 
                neighbourhood.fct = "gaussian")
cov.z <- scale(ls.sam[, .(z1, z2, z3, z4)],
               center = TRUE, scale = TRUE)

head(init_som(na.omit(cov.z), som.dim, som.dim))


soms <- list()
meshes <- list()
r.len <- c(10, 20)
for(i in seq_along(r.len)) {
  soms[[i]] <- som(cov.z,
                 grid = grid, 
                 rlen = r.len[i],
                 radius = som.dim,
                 init = init_som(na.omit(cov.z), som.dim, som.dim),
                 mode = "pbatch", 
                 cores = 4,
                 normalizeDataLayers = FALSE)
  meshes[[i]] <- mesh3d(soms[[i]]$codes[[1]][,1:3], quads = quads)
}

mesh.som.fit <- mesh3d(som.fit$codes[[1]][,1:3], quads = quads)

# lims <- apply(rbind(cov.z, init_som(cov.z, som.dim, som.dim)), 2, range)
# lims <- apply(cov.z, 2, range)
lims <- matrix(rep(c(-3.25, 3.25), 3), ncol = 3, byrow = FALSE)

meshes.plot <- c(list(som.init), meshes, list(mesh.som.fit))

# xlim <- c(-3, 3)
# ylim <- c(-3, 3)
# zlim <- c(-3, 3)
xlim <- c(-3.25, 3.25)
ylim <- c(-3.25, 3.25)
zlim <- c(-3.25, 3.25)
# thetas <- c(10, 95, 175, 260, 45, 45)
# phis <- c(25, 25, 25, 25, 60, -60)
# thetas <- c(seq(0, 359, 90), rep(0, 2))
# phis <- c(rep(0, 4), 90, -90)
thetas <- c(-30, -110)
phis <- c(30, 30)
point.alpha <- 0.15
wire.lwd <- 2
shade.alpha <- 0.2

font.paths <- paste0("/home/danielschoenig/.local/share/fonts/IBMPlexSans-",
                     c("Regular", "Bold", "Italic", "BoldItalic"), ".ttf")
rglFonts(plex = font.paths)
# library(extrafont)
# choose_font("Open Sans")
# fonts()
# rglExtrafonts(plex = "IBM Plex Sans")
# rglExtrafonts(dj = "DejaVu Sans Mono")
# fonts()

for(j in seq_along(thetas)) {
  open3d()
  pars <- par3d() 
  par3d(family = "plex", cex = 1.5)
  plot3d(pts, col = "darkred", type = "s", size = 0.75, alpha = point.alpha, lit = FALSE,
         xlim = lims[,1], ylim = lims[,2], zlim = lims[,3], axes = FALSE,
         xlab = "Cov 1", ylab = "Cov 2", zlab = "Cov 3")
  box3d(lwd = 1.5, xlim = lims[,1], ylim = lims[,2], zlim = lims[,3])
  view3d(theta = thetas[j], phi = phis[j], zoom = 0.9)
  snapshot3d(paste0(path.results, "som_fit_0_", j, ".png"),
             width = 750, height = 750, webshot = FALSE)
  close3d()
}

for(i in seq_along(meshes.plot)) {
  for(j in seq_along(thetas)) {
    open3d()
    pars <- par3d() 
    par3d(family = "plex", cex = 1.5)
    plot3d(pts, col = "darkred", type = "s", size = 0.75, alpha = point.alpha, lit = FALSE,
           xlim = lims[,1], ylim = lims[,2], zlim = lims[,3], axes = FALSE,
           xlab = "Cov 1", ylab = "Cov 2", zlab = "Cov 3")
    box3d(lwd = 1.5, xlim = lims[,1], ylim = lims[,2], zlim = lims[,3])
    wire3d(meshes.plot[[i]], col = "grey20", lit = FALSE, add = TRUE,
           xlim = lims[,1], ylim = lims[,2], zlim = lims[,3], alpha = 1, lwd = 1)
    shade3d(meshes.plot[[i]], col = "darkblue", lit = TRUE, add = TRUE,
            xlim = lims[,1], ylim = lims[,2], zlim = lims[,3], alpha = 0.2)
    view3d(theta = thetas[j], phi = phis[j], zoom = 0.9)
    snapshot3d(paste0(path.results, "som_fit_", i, "_", j, ".png"),
               width = 750, height = 750, webshot = FALSE)
    close3d()
  }
}

som2d <- as.data.table(cbind(som.fit$grid$pts, som.fit$codes[[1]]))

ggplot(som2d) +
  geom_raster(aes(x = x, y = y, fill = z4)) +
  scale_fill_continuous_sequential("viridis", rev = FALSE) +
  labs(x = "SOM 1", y = "SOM 2")


ls.sam[, .(n = .N), by = c("som_x", "som_y", "type")]|>
  ggplot() +
  geom_raster(aes(x = som_x, y = som_y, fill = n)) +
  scale_fill_continuous_sequential("viridis", rev = FALSE) +
  labs(x = "SOM 1", y = "SOM 2") +
  coord_fixed() +
  facet_wrap(vars(type))

## MODEL PREDICTIONS ###########################################################

mapped <-
  ls$landscape[, .(z1, z2, z3, z4)] |>
  scale_data_som(som = som.fit) |>
  embed_som(som = som.fit,
            grid.coord = TRUE)

ls$landscape[,
             `:=`(som_bmu = mapped$bmu[,1],
                  som_x = mapped$grid.coordinates$bmu.1[,"x"],
                  som_y = mapped$grid.coordinates$bmu.1[,"y"])
             ]
    
ls$landscape[, id := cell]
ls$landscape <- ls$landscape[, !"cell"]
setcolorder(ls$landscape, "id")

mean(mod.egp$estimates.noint$effects[,2])

egp.pred <- predict(mod.egp$estimates.noint$gam, ls$landscape, nthreads = 4)

ls$landscape$response.pred <- egp.pred

ls$landscape |>
ggplot(aes(x = x, y = y, fill = response.pred)) +
  geom_raster() +
  scale_fill_continuous_divergingx("Roma", rev = TRUE)

ls$landscape |>
ggplot(aes(x = x, y = y, fill = response)) +
  geom_raster() +
  scale_fill_continuous_divergingx("Roma", rev = TRUE)

som.eval <-
  ls$landscape[, .(z1, z2, z3, z4)] |>
  scale_data_som(som = som.fit) |>
  evaluate_embedding(mapped = mapped$grid.coordinates[[1]], combined = TRUE, bs = "gp", k.max = 250)

quantization_error(som.fit)
topological_error(som.fit)
variance_explained(som.fit, n.cores = 4)

## LS plots


ls_theme <-
  theme_void(base_family = "IBMPlexSans", base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(l = 3, b = 3, t = 3)),
        plot.margin = margin(3, 3, 3, 3),
        legend.position = "bottom",
        legend.justification = c(0.5,1),
        legend.title = element_text(size = rel(0.9), hjust = 0.5,
                                    margin = margin(t = 3, b = 6)),
        legend.text = element_text(size = rel(0.8)),
        legend.spacing.y = unit(2, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.box.spacing = unit(12, "pt")
  )

guide_fill <-
  guides(fill = guide_colorbar(
                               ticks.colour = "grey5",
                               ticks.linewidth = 1,
                               frame.colour = "grey5",
                               frame.linewidth = 0.5,
                               barwidth = 10,
                               barheight = 0.5,
                               title.position = "top",
                               label.position = "bottom",
                               # label.hjust = 0.5,
                               draw.ulim = TRUE,
                               draw.llim = TRUE
                               ))


lab.coord <- ls$landscape[, lapply(.SD, mean), .SDcols = c("x", "y"), by = type]

    
png(paste0(path.results, "ls_response.png"),
    width = 300, height = 300)
  ggplot(ls$landscape, aes(x = x, y = y, fill = response)) +
  geom_raster() +
  # scale_fill_viridis_c()
  scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                   limits = range(ls$landscape$response)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() +
  guide_fill +
  labs(title = "Outcome", fill = "Observed value") +
  ls_theme
dev.off()

 
png(paste0(path.results, "ls_type.png"),
   width = 250, height = 250)
  ggplot(ls$landscape, aes(x = x, y = y)) +
  geom_raster(aes(fill = type)) +
  geom_text(data = lab.coord,
            aes(x = x, y = y,
                label = c("A\n(reference)", "B"),
                color = type), size = 6)  +
  scale_fill_manual(breaks = levels(ls$landscape$type),
                    values = c("grey20", "grey80"),
                    # labels = capwords(levels(ls$landscape$type))) +
                    labels = c("A (reference)", "B")) +
  scale_colour_manual(breaks = levels(ls$landscape$type),
                      values = c("grey95", "grey5")) +
  coord_fixed(expand = FALSE) +
  labs(title = "Area type") +
  guides(fill = "none", color = "none") +
  ls_theme
dev.off()

cov.var <- paste0("z", 1:4)
plots.z <- list()
  for(i in 1:4) {
    plots.z[[cov.var[i]]] <-
      ls$landscape[, .(x, y, value = eval(parse(text = cov.var[i])))] |>
      ggplot() +
      geom_raster(aes(x = x, y = y, fill = value)) +
      scale_fill_continuous_sequential(palette = "Viridis", rev = FALSE, limits = c(0, 1)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed() +
      guide_fill +
      # labs(title = cov.var[i], fill = "Covariate value") +
      labs(title = paste0("Cov ", i), fill = "Covariate value") +
      ls_theme
  }


png(paste0(path.results, "ls_cov.png"),
   width = 450, height = 450)
  plots.z$z1 + plots.z$z2 + plots.z$z3 + plots.z$z4 +
  plot_layout(design = layout1, guides = "collect") &
  ls_theme
dev.off()

  lim.effects <-
    unlist(ls$landscape[,
             .(min(as.matrix(.SD)), max(as.matrix(.SD))),
             .SDcols = c("treatment", paste0("f.", cov.var))],
           use.names = FALSE)
  lim.effects <- 
    (floor(lim.effects / 0.25) + c(0, 1)) * 0.25

png(paste0(path.results, "ls_treatment.png"),
    width = 300, height = 300)
  ggplot(ls$landscape, aes(x = x, y = y, fill = treatment)) +
    geom_raster() +
    scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                     limits = lim.effects) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    guide_fill +
    labs(title = "Effect of area type", fill = "Marginal effect") +
    ls_theme
dev.off()

eff.var <-  paste0("f.", cov.var)
plots.f <- list()
  for(i in seq_along(eff.var)) {
    plots.f[[eff.var[i]]] <-
      ls$landscape[, .(x, y, value = eval(parse(text = eff.var[i])))] |>
      ggplot(aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                       limits = lim.effects) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed() +
      guide_fill +
      labs(title = paste0("Cov ", i), fill = "Effect size") +
      ls_theme
  }

png(paste0(path.results, "ls_cov_effects.png"),
   width = 450, height = 450)
  plots.f$f.z1 + plots.f$f.z2 + plots.f$f.z3 + plots.f$f.z4 +
  plot_layout(design = layout1, guides = "collect") &
  ls_theme
dev.off()


png(paste0(path.results, "ls_response_sam.png"),
    width = 300, height = 300)
  ls$landscape |>
  ggplot(aes(x = x, y = y, fill = response)) +
  geom_raster(alpha = 0.2) +
  geom_raster(data = ls$landscape[.(ls.sam$id), on = "id"], fill = "black") +
  # geom_point(pch = 15, colour = "grey5", size = 0.01) +
  # scale_fill_viridis_c()
  scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                   limits = range(ls$landscape$response)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() +
  guide_fill +
  labs(title = "Outcome", fill = "Observed value") +
  ls_theme
dev.off()

som_theme <-
  theme_void(base_family = "IBMPlexSans", base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(l = 3, b = 3, t = 3)),
        plot.margin = margin(3, 3, 3, 3),
        legend.position = "bottom",
        legend.justification = c(0.5,1),
        legend.title = element_text(size = rel(0.9), hjust = 0.5,
                                    margin = margin(t = 3, b = 6)),
        axis.title = element_text(),
        legend.text = element_text(size = rel(0.8)),
        legend.spacing.y = unit(2, "pt"),
        legend.key.size = unit(12, "pt"),
        legend.box.spacing = unit(12, "pt")
  )


## SOM visualization 2D

som2d <- as.data.table(cbind(som.fit$grid$pts, som.fit$codes[[1]]))

png(paste0(path.results, "som2d.png"), width = 400, height = 400)
expand.grid(x = 1:25,
            y = 1:25) |>
ggplot() +
  geom_rect(aes(xmin = x - 1, ymin = y - 1, xmax = x, ymax = y),
            fill = "darkblue", colour = "grey35", alpha = 0.2) +
  labs(x = "SOM 1", y = "SOM 2", title = "") +
  guides(fill = "none") +
  som_theme
dev.off()

png(paste0(path.results, "som2d_cov1.png"), width = 400, height = 400)
ggplot(som2d) +
  geom_raster(aes(x = x-0.5, y = y-0.5, fill = z1)) +
  scale_fill_continuous_sequential("viridis", rev = FALSE) +
  labs(x = "SOM 1", y = "SOM 2", title = "Cov 1") +
  guides(fill = "none") +
  som_theme
dev.off()

png(paste0(path.results, "som2d_cov2.png"), width = 400, height = 400)
ggplot(som2d) +
  geom_raster(aes(x = x-0.5, y = y-0.5, fill = z2)) +
  scale_fill_continuous_sequential("viridis", rev = FALSE) +
  labs(x = "SOM 1", y = "SOM 2", title = "Cov 2") +
  guides(fill = "none") +
  som_theme
dev.off()

png(paste0(path.results, "ls_response_pred.png"),
    width = 300, height = 300)
  ggplot(ls$landscape, aes(x = x, y = y, fill = response.pred)) +
  geom_raster() +
  # scale_fill_viridis_c()
  scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                   limits = range(ls$landscape$response)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() +
  guide_fill +
  labs(title = "Model (GAM)", fill = "Predicted value") +
  ls_theme
dev.off()

png(paste0(path.results, "ls_response_pred_geo.png"),
    width = 300, height = 300)
  ggplot(ls$landscape, aes(x = x, y = y, fill = response.pred)) +
  geom_raster() +
  # scale_fill_viridis_c()
  scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                   limits = range(ls$landscape$response)) +
  # scale_x_continuous(expand = c(0, 0)) +
  # scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() +
  guide_fill +
  labs(title = "Predictions over landscape", fill = "Posterior mean",
       x = "", y = "") +
  ls_theme
dev.off()

png(paste0(path.results, "ls_response_pred_som.png"),
    width = 308, height = 308)
  ls$landscape[, .(som.pred = mean(response.pred)), c("som_x", "som_y")] |>
    ggplot(aes(x = som_x, y = som_y, fill = som.pred)) +
    geom_raster() +
    # scale_fill_viridis_c()
    scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                     limits = range(ls$landscape$response)) +
    coord_fixed() +
    guide_fill +
    labs(title = "Predictions over SOM", fill = "Posterior mean",
         x = "SOM 1", y = "SOM 2") +
    som_theme
dev.off()

candidates.pt <-
merge(
  dcast(ls.sam[, .(n = .N), c("som_x", "som_y", "type")],
        som_x + som_y ~ type
        )[control >= 3 & !is.na(treatment)],
  ls.sam)[type == "treatment",
         .(id, control, som_x, som_y, x, y,
           diff = sqrt((som_x/25 - x/1000)^2 + (som_y/25 - y/1000)^2))][order(-diff)][1:20]

candidates.pt

ls$landscape[.(ls.sam[.(candidates.pt$id), on = "id"][,.(id)]), on = "id"][, .(id, eff.mar, treatment)]


ls.sam[type == "treatment" & som_x == 15 & som_y == 18]
ls.sam[id == 935853]
ls.sam[type == "control" & som_x == 11 & som_y == 24]

ls$landscape[id == focal.id]

focal.id <- 935853

which(ls.sam$id == focal.id)
which(ls.sam$id %in% ls.sam[type == "control" & som_x == 11 & som_y == 24, id])


i.draws <- c(7, which.max(lp[,as.character(focal.id)]), which.min(lp[,as.character(123611)]))

lp[,as.character(focal.id)][i.draws,]
lp[, as.character(ls.sam[type == "control" & som_x == 11 & som_y == 24, id])][i.draws,]
rowMeans(lp[, as.character(ls.sam[type == "control" & som_x == 11 & som_y == 24, id])][i.draws,])
colMeans(lp[, as.character(ls.sam[type == "control" & som_x == 11 & som_y == 24, id])])

as.vector(lp[,as.character(focal.id)][i.draws,]) - 
rowMeans(lp[, as.character(ls.sam[type == "control" & som_x == 11 & som_y == 24, id])][i.draws,])

# \mu_{1811} &= -1.2
# \mu_{1811} &= -1.2 \ldots -0.8 \ldots -1.0 \ldots

# \begin{align*}
# \mu_{1462} &= -2.2 \\
# \mu_{2763} &= -2.5 \\
# \mu_{3358} &= -2.4 \\
# \vdots
# \end{align*}


# \mu^*_{1811} = -2.4
# \mu^*_{1811} = -2.4 \ldots -2.1, \ldots -2.5 \ldots

# \tilde{\mu_{1811}} = \mu_{1811} - \mu^*_{1811} &= 1.2
# \tilde{\mu_{1811}} = -1.2 \ldots -1.3 \ldots -1.5 \ldots


png(paste0(path.results, "ls_response_pred_geo_pt.png"),
    width = 300, height = 300)
  ggplot(ls$landscape, aes(x = x, y = y, fill = response.pred)) +
  geom_raster() +
  geom_point(data = ls.sam[id == focal.id],
             pch = 1, fill = NA, colour = "grey5", size = 3, stroke = 2) +
  # scale_fill_viridis_c()
  scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                   limits = range(ls$landscape$response)) +
  # scale_x_continuous(expand = c(0, 0)) +
  # scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() +
  guide_fill +
  labs(title = "Predictions over landscape", fill = "Posterior mean",
       x = "", y = "") +
  ls_theme
dev.off()

png(paste0(path.results, "ls_response_pred_som_pt.png"),
    width = 308, height = 308)
  ls$landscape[, .(som.pred = mean(response.pred)), c("som_x", "som_y")] |>
  ggplot(aes(x = som_x, y = som_y, fill = som.pred)) +
  geom_raster() +
  geom_point(data = ls.sam[id == focal.id],
             pch = 1, fill = NA, colour = "grey5", size = 3, stroke = 2) +
  # scale_fill_viridis_c()
  scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                   limits = range(ls$landscape$response)) +
  coord_fixed() +
  guide_fill +
  labs(title = "Predictions over SOM", fill = "Posterior mean",
       x = "SOM 1", y = "SOM 2") +
  som_theme
dev.off()

png(paste0(path.results, "ls_eff.png"), width = 450, height = 450)
data.table(type = "treatment",
           effect = extract_variable(mod.egp$estimates.noint$effects,
                                     "treatment")) |>
ggplot() +
  stat_halfeye(aes(x = effect, y = type), alpha = 0.75) +
  scale_x_continuous(limits = c(0.7, 1)) +
  scale_y_discrete(labels = c(treatment = "Treatment")) +
  coord_cartesian(ylim = c(3, NA)) +
  labs(x = "Marginal effect", y = "") +
  theme_minimal(base_family = "IBMPlexSans", base_size = 16) +
  theme(
        axis.line.x = element_line(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 14, colour = "grey5"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 12)),
        panel.grid = element_blank())
dev.off()

poly_treat <-
  st_as_stars(ls$landscape[type == "treatment", .(x, y, type)]) |>
  st_as_sf() |>
  st_union(is_coverage = TRUE) |>
  st_as_sf()

plot(poly_treat)

png(paste0(path.results, "ls_response_pred_geo_treat.png"),
    width = 300, height = 300)
  ggplot(ls$landscape) +
  geom_raster(aes(x = x, y = y, fill = response.pred)) +
  geom_sf(data = poly_treat, fill = NA, colour = "grey5", lwd = 1.5) +
  # scale_fill_viridis_c()
  scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                   limits = range(ls$landscape$response)) +
  # scale_x_continuous(expand = c(0, 0)) +
  # scale_y_continuous(expand = c(0, 0)) +
  # coord_fixed() +
  guide_fill +
  labs(title = "Predictions over landscape", fill = "Posterior mean",
       x = "", y = "") +
  ls_theme
dev.off()



## Marginal effect for individual observations


# system.time({
# data.bl <- assign_bl_som(ls$landscape,
#                          som = som.fit,
#                          cov.col = c("z1", "z2", "z3", "z4"),
#                          id.col = "id")
# }) # Takes around 10 minutes
file.bl <- paste0(path.results, "data.bl_", ls.id, ".rds")
# saveRDS(data.bl, file.bl)
data.bl <- readRDS(file.bl)

post <- mod.egp$estimates.noint$post[1:100,]

lp <-
  evaluate_posterior(model = mod.egp$estimates.noint$gam,
                     posterior = post,
                     newdata = ls$landscape,
                     id.col = "id",
                     type = "response",
                     progress = TRUE,
                     predict.chunk = 1e5)

id.list <- as.list(ls$landscape$id)
names(id.list) <- colnames(lp)

groups.som <-
    ls$landscape[type == "control"] |>
    ids_by_group(id.col = "id", group.vars = "som_bmu")
id.list.som <- groups.som$ids
names(id.list.som) <- groups.som$som_bmu
yhat.som <- aggregate_variables(lp,
                    agg.fun = mean,
                    ids = id.list.som,
                    agg.size = 1e5
                    )

ids.units <- data.bl[,
                     .(id,
                       som_bmu.bl,
                       som_bmu.bl.w)
                     ][,
                       lapply(.SD, unlist), c("id")
                       ][,
                         .(id,
                           som_bmu.bl,
                           som_bmu.bl.w)]
setkey(ids.units, id)

# Reweigh baseline SOM units for each group, based on what points they where
# assigned to
system.time({
w.points <-
  lapply(id.list,
         \(x) {
               extract_weights(ids.units[.(x)],
                               w.col = "som_bmu.bl.w",
                               by.col = "som_bmu.bl",
                               standardize = TRUE)
              })
yhat.bl <- reweigh_posterior(yhat.som, w = w.points)
})



eff.mar <- arc(lp, yhat.bl)

eff.mar.means <- apply(eff.mar, 2, mean)

# export.l <- list(data.bl = data.bl,
#                  lp = lp, id.list = id.list,
#                  groups.som = groups.som, id.list.som = id.list.som,
#                  yhat.som = yhat.som, ids.units = ids.units,
#                  w.points = w.points, yhat.bl = yhat.bl,
#                  eff.mar = eff.mar, eff.mar.means = eff.mar.means)
# saveRDS(export.l, paste0(path.results, "mar.pred.rds"))
# export.l <- readRDS(paste0(path.results, "mar.pred.rds"))

ls$landscape$eff.mar <- eff.mar.means


ls$landscape |>
ggplot(aes(x = x, y = y, fill = eff.mar)) +
  geom_raster() +
  scale_fill_continuous_divergingx("Roma", rev = TRUE)

ls$landscape |>
ggplot(aes(x = x, y = y, fill = treatment)) +
  geom_raster() +
  scale_fill_continuous_divergingx("Roma", rev = TRUE)

ls.sam2 <- copy(ls.sam)

ls.sam2[, type := ordered(type)]

    mod2 <- bam(response ~
                   s(x, y, bs = "gp", m = 3, k = 250) +
                   s(x, y, by = type, bs = "gp", k = 250) +
                   s(som_x, som_y, bs = "gp", k = 250,
                     xt = list(max.knots = 625)),
                   data = ls.sam2,
                   select = TRUE
                  )

sd(fitted(mod2) - fitted(mod.egp$estimates.noint$gam))

summary(mod2)




## Marginal effect for quadrants

data.bl <- assign_bl_som(ls.sam,
                         som = som.fit,
                         cov.col = c("z1", "z2", "z3", "z4"),
                         id.col = "id")

post <- mod.egp$estimates.noint$post

lp <-
  evaluate_posterior(model = mod.egp$estimates.noint$gam,
                     posterior = post,
                     newdata = ls.sam,
                     id.col = "id",
                     type = "response",
                     progress = TRUE,
                     predict.chunk = 1e5)


binned.9 <- bin_cols(ls.sam, c("x", "y"), c(1000/3, 1000/3), c(0, 0), bin.names = c("x9", "y9")) |>
lapply(\(x) round(x, 2))
binned.16 <- bin_cols(ls.sam, c("x", "y"), c(1000/4, 1000/4), c(0, 0), bin.names = c("x16", "y16")) |>
lapply(\(x) round(x, 2))
ls.sam <- cbind(ls.sam, do.call(cbind, binned.9), do.call(cbind, binned.16))

groups.bin9 <-
    ls.sam |>
    ids_by_group(id.col = "id", group.vars = c("x9", "y9"),
                 # group.labels =  list(c("166.67" = "A", "500" = "B", "833.33" = "C"),
                 #                      c("166.67" = "A", "500" = "B", "833.33" = "C")),
                 expand.label = TRUE)
groups.bin16 <-
    ls.sam |>
    ids_by_group(id.col = "id", group.vars = c("x16", "y16"),
                 # group.labels =  list(c("125" = "A", "375" = "B", "625" = "C", "875" = "D"),
                 #                      c("125" = "A", "375" = "B", "625" = "C", "875" = "D")),
                 expand.label = TRUE)
groups.bin <- rbind(groups.bin9, groups.bin16, fill = TRUE)
id.list.bin <- groups.bin$ids
names(id.list.bin) <- groups.bin$group.label

yhat.bin <- aggregate_variables(lp,
                    agg.fun = mean,
                    ids = id.list.bin,
                    agg.size = 1e5
                    )

groups.som <-
    ls.sam |>
    ids_by_group(id.col = "id", group.vars = "som_bmu")
id.list.som <- groups.som$ids
names(id.list.som) <- groups.som$som_bmu
yhat.som <- aggregate_variables(lp,
                    agg.fun = mean,
                    ids = id.list.som,
                    agg.size = 1e5
                    )

ids.units <- data.bl[,
                     .(id,
                       som_bmu.bl,
                       som_bmu.bl.w)
                     ][,
                       lapply(.SD, unlist), c("id")
                       ][,
                         .(id,
                           som_bmu.bl,
                           som_bmu.bl.w)]
setkey(ids.units, id)

# Reweigh baseline SOM units for each group, based on what points they where
# assigned to
system.time({
w.points <-
  lapply(id.list.bin,
         \(x) {
               extract_weights(ids.units[.(x)],
                               w.col = "som_bmu.bl.w",
                               by.col = "som_bmu.bl",
                               standardize = TRUE)
              })
yhat.bl <- reweigh_posterior(yhat.som, w = w.points)
})



eff.mar <- arc(yhat.bin, yhat.bl)

comp.quads <-
merge(
  rbind(ls.sam[order(x9, y9), .(obs.mean = mean(treatment)), by = c("x9", "y9")],
        ls.sam[order(x16, y16), .(obs.mean = mean(treatment)), by = c("x16", "y16")],
        fill = TRUE),
  cbind(groups.bin[, .(x9, y9, x16, y16)], mar.mean = apply(eff.mar, 2, mean)),
  sort = FALSE
  )

comp.quads.l <-
  melt(comp.quads, measure.vars = c("obs.mean", "mar.mean"))

# comp.quads[, relbias := (mar.mean)/obs.mean]
# comp.quads[,mean(relbias)]

ls$landscape[, treatment0 := treatment]
ls$landscape[type == "control", treatment0 := 0]


ls$landscape |>
ggplot(aes(x = x, y = y, fill = treatment0)) +
  geom_raster() +
  scale_fill_continuous_divergingx("Roma", rev = TRUE)


comp.quads.l |>
ggplot(aes(x = x9, y = y9, fill = value)) +
  geom_raster() +
  scale_fill_continuous_divergingx("Roma", rev = TRUE) +
  coord_fixed() +
 facet_wrap(vars(variable))

comp.quads.l |>
ggplot(aes(x = x16, y = y16, fill = value)) +
  geom_raster() +
  scale_fill_continuous_divergingx("Roma", rev = TRUE) +
  coord_fixed() +
  facet_wrap(vars(variable))

