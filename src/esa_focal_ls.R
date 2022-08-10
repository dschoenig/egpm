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
egp.type <- "egp_sam0.01_som50"

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
  
som.dim = 50

pts <- scale(as.matrix(ls.sam[1:1000,.(z1, z2, z3)]),
             center = TRUE, scale = TRUE)
quads <- draw_quad(som.dim, som.dim, 4)

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


egp.pred <- predict(mod.egp$estimates.noint$gam, ls$landscape, nthreads = 4)

ls$landscape$response.pred <- egp.pred

ls$landscape |>
ggplot(aes(x = x, y = y, fill = response.pred)) +
  geom_raster() +
  scale_fill_continuous_divergingx("Roma", rev = TRUE)




# system.time({
# data.bl <- assign_bl_som(ls$landscape,
#                          som = som.fit,
#                          cov.col = c("z1", "z2", "z3", "z4"),
#                          id.col = "id")
# }) Takes around 10 minutes
file.bl <- paste0(path.results, "data.bl_", ls.id, ".rds")
# saveRDS(data.bl, file.bl)
data.bl <- readRDS(file.bl)

post <- mod.egp$estimates.noint$post

lp <-
  evaluate_posterior(model = mod.egp$estimates.noint$gam,
                     posterior = post,
                     newdata = ls$landscape,
                     id.col = "id",
                     type = "response",
                     progress = TRUE,
                     predict.chunk = 1e5)

groups.type <-
    ls$landscape |>
    ids_by_group(id.col = "id", group.vars = c("x", "y"))

id.list.type <- groups.type$ids
names(id.list.type) <- groups.type$type
yhat.type <- aggregate_variables(lp,
                    agg.fun = mean,
                    ids = id.list.type,
                    )

groups.som <-
    ls.sam[type == "control"] |>
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
w.points <-
  lapply(id.list.type,
         \(x) {
               extract_weights(ids.units[.(x)],
                               w.col = "som_bmu.bl.w",
                               by.col = "som_bmu.bl",
                               standardize = TRUE)
              })
yhat.bl.type <- reweigh_posterior(yhat.som, w = w.points)

eff.mar <- arc(yhat.type, yhat.bl.type)




