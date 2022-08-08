library(mgcv)
library(posterior)
library(RandomFields)
library(RandomFieldsUtils)
library(raster)
library(sf)
library(stars)
library(spdep)
library(igraph)
library(data.table)
library(mvnfast)
library(colorspace)
library(stringi)
library(ggplot2)
library(patchwork)
library(kohonen)

RFoptions(install="no")

equilibrate <- function(x, ...) {
  UseMethod("equilibrate", x)
}

equilibrate.matrix <- function(x, tolerance = .Machine$double.eps) {
  while(any(abs(rowSums(x) - 1) > .Machine$double.eps) |
        any(abs(colSums(x) - 1) > .Machine$double.eps)) {
    x <- apply(x, 2, \(x) x / sum(x))
    x <- apply(x, 1, \(x) x / sum(x))
  }
  return(x)
}

point_on_square <- function(phi, prec = 6){
  phi <- phi %% (2*pi)
  # tangens
  if((phi >= 0 & phi <= pi/4) |
     (phi >= 7*pi/4 & phi < 2*pi)) {
    crd <-c(1, tan(phi))
  }
  if(phi >= 3*pi/4 & phi <= 5*pi/4) {
    crd <- c(-1, -tan(phi))
  }
  # cotangens
  if(phi > pi/4 & phi < 3*pi/4) {
    crd <- c(1/tan(phi), 1)
  }
  if(phi > 5*pi/4 & phi < 7*pi/4) {
    crd <- c(-1/tan(phi), -1)
  }
  crd <- (crd + 1) / 2
  return(round(crd, 6))
}


scale_int <- function(x, ...) {
  UseMethod("scale_int", x)
}

scale_int.stars <- function(x, int = c(0, 1)) {
  if(!is.null(int)) {
    n.att <- length(x)
    for(i in 1:n.att) {
      min.att <- min(x[[i]], na.rm = TRUE)
      max.att <- max(x[[i]], na.rm = TRUE)
      x[[i]] <- int[1] + (int[2] - int[1]) * (x[[i]] - min.att) / (max.att - min.att)
    }
  }
  return(x)
}


interpolate_na <- function(x, ...) {
  UseMethod("interpolate_na", x)
}


interpolate_na.stars <- function(x, method = "average", nh = "moore", range = 1) {
  # Default: Moore neighbourhood
  window <- expand.grid(seq(-range, range, 1), seq(-range, range, 1))
  window <- unname(as.matrix(window[-which(window[,1] == 0 & window[,2] == 0),]))
  if(nh == "neumann") {
    dist <- apply(abs(window), 1, sum)
    window <- window[which(dist <= range),]
  }
  n.att <- length(x)
  for(i in 1:n.att) {
    ind.na <- which(is.na(x[[i]]), arr.ind = TRUE)
    dims <- dim(x[[i]])
    inter <- NA
    if(length(ind.na) > 0) {
      for(j in 1:nrow(ind.na)) {
        foc <- ind.na[j,]
        ipol <- t(apply(window, 1, \(x) x + foc))
        rem <- c(which(ipol[,1] > dims[1] | ipol[, 1] < 1),
                 which(ipol[,2] > dims[2] | ipol[, 2] < 1))
        if(length(rem) > 0){
          ipol <- ipol[-rem,]
        }
        if(method == "average"){
          inter[j] <- mean(x[[i]][ipol], na.rm = TRUE)
        }
        if(method == "gam") {
          mod.df <- na.omit(data.frame(x = ipol[,1], y = ipol[,2], z = x[[i]][ipol]))
          pred.df <- data.frame(x = unname(foc[1]), y = unname(foc[2]))
          mod <- gam(z ~ s(x, y, k = nrow(mod.df)), data = mod.df)
          inter[j] <- predict(object = mod, newdata = pred.df)
        }
      }
      x[[i]][ind.na] <- inter
    }
  }
  return(x)
}


scale_int.numeric <- function(x, int = c(0, 1)) {
  if(!is.null(int)) {
      min.x <- min(x, na.rm = TRUE)
      max.x <- max(x, na.rm = TRUE)
      y <- int[1] + (int[2] - int[1]) * (x - min.x) / (max.x - min.x)
  }
  return(y)
}


matrix2stars <- function(x, res = 1, name = "value", ...) {
  x.stars <-
    st_as_stars(x) |>
    st_set_dimensions(names = c("x", "y")) |>
    st_set_dimensions("x", offset = 0, delta = res) |>
    st_set_dimensions("y", offset = nrow(x) * res, delta = - res) |>
    setNames(name)
  return(x.stars)
}

reset_dim <- function(x, ...) {
  UseMethod("reset_dim", x)
}

reset_dim.stars <- function(x, y = NULL, ...) {
  if(is.null(y)) {
    dims <- dim(x)
    y <- generate_empty(x.dim = dims[1], y.dim = dims[2])
  }
  st_dimensions(x) <- st_dimensions(y)
  return(x)
}


limit_polygon <- function(x.dim,
                          y.dim) {
  limit <-
    list(matrix(c(0, 0,
                x.dim, 0,
                x.dim, y.dim,
                0, y.dim,
                0, 0),
              nrow = 5, ncol = 2,
              byrow = TRUE)) |>
    st_polygon() |>
    st_sfc()
  return(limit)
}

generate_empty <- function(x.dim,
                           y.dim,
                           name = "value") {
  mat <- matrix(0, nrow = x.dim, ncol = y.dim)
  r <- matrix2stars(mat, name = name)
  return(r)
}

generate_linear_gradient <- function(x.dim,
                                     y.dim,
                                     phi,
                                     rescale = c(0,1),
                                     name = "value") {
  grad.mat <- outer(1:x.dim, 1:y.dim, FUN = \(x,y) x * cos(phi) - y * sin(phi))
  grad <-
    matrix2stars(grad.mat, name = name) |> 
    scale_int(int = rescale)
  return(grad)
}

generate_sigmoid_gradient <- function(x.dim,
                                     y.dim,
                                     phi,
                                     rescale = c(0,1),
                                     name = "value") {
  lin.grad <-
    generate_linear_gradient(x.dim = x.dim,
                             y.dim = y.dim,
                             phi = phi,
                             rescale = c(0, 1),
                             name = name)
  sig.grad <- (1/ (1 + exp(-lin.grad))) |>
    scale_int(int = rescale)
  return(sig.grad)
}



generate_exp_gradient <- function(x.dim,
                                  y.dim,
                                  phi = NULL,
                                  shift = NULL,
                                  centre = NULL,
                                  x.scale = 1,
                                  y.scale = 1,
                                  rescale = c(0,1),
                                  name = "value") {
  if(is.null(centre)) {
    centre <- point_on_square(phi)
    centre <- centre - (shift * (centre - 0.5))
  } else {
    centre <- centre / c(x.dim, y.dim)
  }
  scale <- pi / sqrt(x.dim^2 + y.dim^2)
  x.seq <- (1:x.dim) * scale / x.scale
  y.seq <- (y.dim:1) * scale / y.scale
  centre <- centre * c(max(x.seq), max(y.seq))
  grad.mat <-
    outer(x.seq, y.seq,
          # FUN = \(x,y) exp(-((x - x.off)^2 + (y - y.off)^2)))
          FUN = \(x,y) exp(-((x - centre[1])^2 + (y - centre[2])^2)))
  matrix2stars(grad.mat)
  grad <-
      matrix2stars(grad.mat, name = name) |> 
      scale_int(int = rescale)
  return(grad)
}

generate_fbm <- function(x.dim,
                         y.dim,
                         alpha = 1,
                         var = 1,
                         scale = 1,
                         rescale = c(0,1),
                         name = "value") {
  mod <- RMfbm(alpha = alpha, var = var, scale = scale)
  fbm <- 
    RFsimulate(model = mod, x = 1:x.dim, y = 1:y.dim) |>
    as.matrix() |>
    t() |>
    matrix2stars(name = name) |>
    scale_int(int = rescale)
  return(fbm)
}

generate_matern <- function(x.dim,
                            y.dim,
                            nu = 1,
                            var = 1,
                            scale = 1,
                            mean = 0,
                            rescale = c(0,1),
                            name = "value") {
  mod <- RMmatern(nu = nu, var = var, scale = scale) + RMtrend(mean = mean)
  matern <- 
    RFsimulate(model = mod, x = 1:x.dim, y = 1:y.dim) |>
    as.matrix() |>
    t() |>
    matrix2stars(name = name) |>
    scale_int(int = rescale)
  return(matern)
}

generate_mpd <- function(x.dim,
                         y.dim,
                         r,
                         rdev,
                         rescale = c(0, 1),
                         name = "value") {
  x.dim.mpd <- ifelse(x.dim %% 2 == 0, x.dim + 3, x.dim + 2)
  y.dim.mpd <- ifelse(y.dim %% 2 == 0, x.dim + 3, x.dim + 2)
  mpd.mat <-
    nlm_mpd(ncol = x.dim.mpd,
            nrow = y.dim.mpd,
            roughness = r,
            rand_dev = rdev,
            verbose = FALSE,
            rescale = FALSE) |>
    as.matrix() |>
    t()
  mpd <-
    matrix2stars(mpd.mat[1:x.dim, 1:y.dim]) |>
    scale_int(int = rescale) |>
    setNames(name)
  return(mpd)
}

# generate_distance <- function(x.dim,
#                               y.dim,
#                               n,
#                               dist.acc = 0.1,
#                               rescale = c(0,1),
#                               name = "value") {
#   poly.lim <- limit_polygon(x.dim, y.dim)
#   target <- st_sample(poly.lim, n) |> st_union()
#   x.crd <- seq(1, x.dim, 1 / dist.acc) + (1/(2*dist.acc))
#   y.crd <- seq(1, y.dim, 1 / dist.acc) + (1/(2*dist.acc))
#   points <-
#     cbind(x = rep(x.crd, times = length(y.crd)),
#           y = rep(y.crd, each = length(y.crd))) |>
#     st_multipoint() |>
#     st_sfc() |>
#     st_cast("POINT") |>
#     st_as_sf()
#   points$dist <- st_distance(points, target)
#   r <-
#     st_rasterize(points) |>
#     st_warp(use_gdal = TRUE, cellsize = 1,
#             method = "cubicspline", no_data_value = 0) |>
#     setNames(name) |>
#     scale_int(int = rescale)
#   return(r)
# }

generate_distance <- function(x.dim,
                              y.dim,
                              n,
                              grad.phi,
                              grad.prop = 1,
                              dist.acc = 0.1,
                              rescale = c(0,1),
                              name = "value") {
  poly.lim <- limit_polygon(x.dim, y.dim)
  target <- st_sample(poly.lim, n) |> st_union()
  grad.lin <- 
    generate_linear_gradient(x.dim = x.dim,
                             y.dim = y.dim,
                             phi = grad.phi,
                             rescale = c(1 - grad.prop, 1))
  grad.lin <-
    as.data.table(as.data.frame(grad.lin))
  nuclei.id <- sample(1:nrow(grad.lin), size = n, prob = grad.lin$value)
  target <-
    grad.lin[nuclei.id, .(x, y)] |>
    as.matrix() |>
    st_multipoint() |>
    st_sfc()
  x.crd <- seq(1, x.dim, 1 / dist.acc) + (1/(2*dist.acc))
  y.crd <- seq(1, y.dim, 1 / dist.acc) + (1/(2*dist.acc))
  points <-
    cbind(x = rep(x.crd, times = length(y.crd)),
          y = rep(y.crd, each = length(y.crd))) |>
    st_multipoint() |>
    st_sfc() |>
    st_cast("POINT") |>
    st_as_sf()
  points$dist <- st_distance(points, target)
  dist <-
    st_rasterize(points) |>
    st_warp(use_gdal = TRUE, cellsize = 1,
            method = "cubicspline", no_data_value = 0)
  dist <- reset_dim(dist)
  # st_dimensions(dist) <-
  #   st_dimensions(generate_empty(x.dim = x.dim, y.dim = y.dim))
  dist <-
    setNames(dist, name) |>
    scale_int(int = rescale)
  return(dist)
}

generate_segments <- function(x.dim,
                              y.dim,
                              n,
                              rescale = c(0,1),
                              name = "value") {
  poly.lim <- limit_polygon(x.dim, y.dim)
  poly.sam <- st_sample(poly.lim, n)
  seg <- 
    poly.sam |>
    st_union() |>
    st_voronoi() |>
    st_cast() |>
    st_intersection(poly.lim) |>
    st_as_sf() |>
    st_set_geometry("geometry")
  seg$segment <- 1:n
  return(seg)
}


segment_field <- function(field,
                          segments,
                          fun = mean) {
  x.dim <- dim(field)[1]
  y.dim <- dim(field)[2]
  seg <-
    aggregate(field, segments, FUN = fun) |>
    st_as_sf() |>
    st_rasterize(template = field)
  seg <- reset_dim(seg)
  # seg <- generate_empty(x.dim = x.dim, y.dim = y.dim) + seg
  return(seg)
}


generate_split <- function(x.dim,
                           y.dim,
                           n,
                           prop = 0.5,
                           name = "split") {
  poly.lim <- limit_polygon(x.dim, y.dim)
  poly <- 
    poly.lim |>
    st_sample(n) |>
    st_union() |>
    st_voronoi() |>
    st_cast() |>
    st_intersection(poly.lim) |>
    st_as_sf()
  poly$id <- 1:n
  poly$area <- st_area(poly) / prod(x.dim, y.dim)
  Aw <-
    poly2nb(poly) |>
    nb2listw(style = "B") |> 
    listw2mat()
  G <- graph_from_adjacency_matrix(Aw)
  path <- eulerian_path(G)$vpath
  upath <- unique(path)
  idx.split <- upath[cumsum(poly$area[upath]) > prop]
  poly$split <- 0
  poly$split[idx.split] <- 1
  poly0 <- subset(poly, split == 0) |> st_union()
  poly1 <- subset(poly, split == 1) |> st_union()
  split <-
    st_as_sf(c(poly0, poly1)) |>
    st_set_geometry("geometry")
  split[[name]] <- c(0, 1)
  return(split)
}


# generate_z1 <- function(x.dim,
#                         y.dim,
#                         fbm1.alpha,
#                         fbm1.var,
#                         fbm1.scale,
#                         fbm2.alpha,
#                         fbm2.var,
#                         fbm2.scale,
#                         fbm.ratio,
#                         grad.phi,
#                         rescale = c(0, 1),
#                         name = "z1") {
#   fbm1 <- generate_fbm(x.dim = x.dim,
#                        y.dim = y.dim,
#                        alpha = fbm1.alpha,
#                        var = fbm1.var,
#                        scale = fbm1.scale,
#                        rescale = c(0, 1))
#   fbm2 <- generate_fbm(x.dim = x.dim,
#                        y.dim = y.dim,
#                        alpha = fbm2.alpha,
#                        var = fbm2.var,
#                        scale = fbm2.scale,
#                        rescale = c(0, 1))
#   grad.lin <- generate_linear_gradient(x.dim = x.dim,
#                                        y.dim = y.dim,
#                                        phi = grad.phi,
#                                        rescale = c(0, 1))
#   z <-
#      (fbm1 * fbm.ratio * grad.lin) + 
#      (fbm2  * (1 - grad.lin))
#   z <-
#     scale_int(z, int = rescale) |>
#     setNames(name)
# }


generate_z1 <- function(x.dim,
                        y.dim,
                        fbm.alpha,
                        fbm.var,
                        fbm.scale,
                        fbm.w,
                        grad.phi,
                        grad.w,
                        rescale = c(0, 1),
                        name = "z1") {
  fbm <- generate_fbm(x.dim = x.dim,
                       y.dim = y.dim,
                       alpha = fbm.alpha,
                       var = fbm.var,
                       scale = fbm.scale,
                       rescale = c(0, 1))
  grad.lin <- generate_linear_gradient(x.dim = x.dim,
                                       y.dim = y.dim,
                                       phi = grad.phi,
                                       rescale = c(0, 1))
  z <- fbm.w * fbm + grad.w * grad.lin
  z <-
    scale_int(z, int = rescale) |>
    setNames(name)
}



generate_z2 <- function(x.dim,
                        y.dim,
                        fbm.alpha,
                        fbm.var,
                        fbm.scale,
                        fbm.w,
                        grad.phi,
                        grad.shift,
                        grad.w,
                        rescale = c(0, 1),
                        name = "z2") {
  fbm <- generate_fbm(x.dim = x.dim,
                    y.dim = y.dim,
                    alpha = fbm.alpha,
                    var = fbm.var,
                    scale = fbm.scale,
                    rescale = c(0, 1))
  # grad <- generate_linear_gradient(x.dim = x.dim,
  #                                  y.dim = y.dim,
  #                                  phi = grad.phi,
  #                                  rescale = c(0, 1))
  grad <- generate_exp_gradient(x.dim = x.dim,
                                   y.dim = y.dim,
                                   phi = grad.phi,
                                   shift = grad.shift,
                                   rescale = c(0, 1))
  z <- fbm.w * fbm + grad.w * grad
  z <-
    scale_int(z, int = rescale) |>
    setNames(name)
}

generate_z3 <- function(
                        x.dim,
                        y.dim,
                        dist.n,
                        grad.phi,
                        grad.prop,
                        acc = 0.1,
                        rescale = c(0,1),
                        name = "z3"
                        ) {
  dist <- generate_distance(x.dim = x.dim,
                            y.dim = y.dim,
                            n = dist.n,
                            grad.phi = grad.phi,
                            grad.prop = grad.prop,
                            dist.acc = acc,
                            rescale = c(0,1))
  z <- interpolate_na(dist, method = "gam", nh = "neumann", range = 2)
  z <-
    scale_int(z, int = rescale) |>
    setNames(name)
  return(z)
}

generate_z4 <- function(x.dim,
                        y.dim,
                        seg.n,
                        mat.nu,
                        mat.var,
                        mat.scale,
                        mat.w,
                        grad.phi,
                        grad.w,
                        rescale = c(0,1),
                        name = "z4") {
  seg <- generate_segments(x.dim = x.dim,
                           y.dim = y.dim,
                           n = seg.n)
  mat <-
    generate_matern(x.dim = x.dim,
                    y.dim = y.dim,
                    nu = mat.nu,
                    scale = mat.scale,
                    var = mat.var,
                    rescale = c(0, 1)) |>
    segment_field(seg) |>
    scale_int()
  grad.lin <-
    generate_linear_gradient(x.dim = x.dim,
                             y.dim = y.dim,
                             phi = grad.phi,
                             rescale = c(0, 1)) |>
    segment_field(seg) |>
    scale_int()
  z <- mat.w * mat + grad.w * grad.lin
  z <-
    scale_int(z, int = rescale) |>
    setNames(name)
  return(z)
}


x.dim = 100
y.dim = 100
effect.size.sp = 1
effect.size.bd = 0
nuclei = 10
poly = generate_split(100, 100, 200, 0.5)
phi.range = c(-1, 1)
shift.range = c(0, 0.5)
x.scale.range = c(0.5, 1)
y.scale.range = c(0.5, 1)
nuc.eff.range = c(0.5, 2)
name = "treatment"

generate_treatment <- function(x.dim,
                               y.dim,
                               effect.size.sp,
                               effect.size.bd,
                               nuclei,
                               poly,
                               phi.range,
                               shift.range,
                               x.scale.range,
                               y.scale.range,
                               nuc.eff.range,
                               name = "treatment") {
  treatment <- generate_empty(x.dim = x.dim,
                              y.dim = y.dim)
  centroids <-
    st_centroid(poly[1:2,]) |>
    st_coordinates() |>
    suppressWarnings()
  phi.shift <-
    atan2((centroids[2,2] - centroids[1,2]),
          (centroids[2,1] - centroids[1,1]))
  phi.range.trt <- phi.range + phi.shift
  phi.range.ctr <- phi.range + phi.shift + pi
  phi.trt <-
    sample(seq(phi.range.trt[1], phi.range.trt[2], length.out = 2 * nuclei), nuclei)
  phi.ctr <-
    sample(seq(phi.range.ctr[1], phi.range.ctr[2], length.out = 2 * nuclei), nuclei)
  # phi.trt <- runif(n,
  #                  phi.range[1] + phi.shift,
  #                  phi.range[2] + phi.shift)
  # phi.ctr <- runif(n,
  #                  phi.range[1] + phi.shift + pi,
  #                  phi.range[2] + phi.shift + pi)
  shift.trt <- runif(nuclei,
                     shift.range[1],
                     shift.range[2])
  shift.ctr <- runif(nuclei,
                     shift.range[1],
                     shift.range[2])
  x.scale.trt <- runif(nuclei,
                     x.scale.range[1],
                     x.scale.range[2])
  x.scale.ctr <- runif(nuclei,
                     x.scale.range[1],
                     x.scale.range[2])
  y.scale.trt <- runif(nuclei,
                     y.scale.range[1],
                     y.scale.range[2])
  y.scale.ctr <- runif(nuclei,
                     y.scale.range[1],
                     y.scale.range[2])
  trt.eff <- runif(nuclei,
                   nuc.eff.range[1] * effect.size.sp,
                   nuc.eff.range[2] * effect.size.sp)
  ctr.eff <- runif(nuclei,
                   nuc.eff.range[1] * effect.size.sp,
                   nuc.eff.range[2] * effect.size.sp)
  for(i in 1:nuclei) {
    trt.grad <-
      generate_exp_gradient(x.dim = x.dim,
                            y.dim = y.dim,
                            phi = phi.trt[i],
                            shift = shift.trt[i],
                            x.scale = x.scale.trt[i],
                            y.scale = y.scale.trt[i])
    ctr.grad <-
      generate_exp_gradient(x.dim = x.dim,
                            y.dim = y.dim,
                            phi = phi.ctr[i],
                            shift = shift.ctr[i],
                            x.scale = x.scale.ctr[i],
                            y.scale = y.scale.ctr[i])
    treatment <- trt.eff[i] * trt.grad - ctr.eff[i] * ctr.grad + treatment
  }
  eff.scale <-
    effect.size.sp /
    (mean(treatment[poly[2,]][[1]], na.rm = TRUE) -
     mean(treatment[poly[1,]][[1]], na.rm = TRUE))
  treatment <- treatment * eff.scale
  eff.shift <- mean(treatment[poly[1,]][[1]], na.rm = TRUE)
  treatment <- treatment - eff.shift
  treatment <- setNames(treatment, name)
  poly.bd <-
    st_geometry(poly[2,]) |>
    st_as_sf()
  poly.bd$eff.bd <- effect.size.bd
  treatment.bd <- st_rasterize(poly.bd, template = generate_empty(x.dim = x.dim, y.dim = y.dim))
  treatment <- treatment + reset_dim(treatment.bd, treatment)
  return(treatment)
}


# generate_treatment <- function(x.dim,
#                                y.dim,
#                                effect.size.sp,
#                                effect.size.bd,
#                                nuclei,
#                                poly,
#                                phi.range,
#                                shift.range,
#                                x.scale.range,
#                                y.scale.range,
#                                nuc.eff.range,
#                                grad.prop,
#                                mean,
#                                name = "treatment") {
#   treatment <- generate_empty(x.dim = x.dim,
#                               y.dim = y.dim)
#   centroids <-
#     st_centroid(poly[1:2,]) |>
#     st_coordinates() |>
#     suppressWarnings()
#   phi.trt <-
#     atan2((centroids[2,2] - centroids[1,2]),
#           (centroids[2,1] - centroids[1,1]))
#   # prob.grad <-
#   #   generate_linear_gradient(x.dim = x.dim,
#   #                            y.dim = y.dim,
#   #                            phi = phi.trt,
#   #                            rescale = c(1-grad.prop, 1)) |>
#   # as.data.table()
#   prob.grad <-
#     generate_linear_gradient(x.dim = x.dim,
#                              y.dim = y.dim,
#                              phi = phi.trt,
#                              rescale = c(0, 1))
#   prob.grad <-
#     (1/ (1 + exp(prob.grad))) |>
#     scale_int(int = c(0, 1)) |>
#     as.data.table()
#   centres.trt <-
#     prob.grad[sample(1:nrow(prob.grad), nuclei, prob = prob.grad$value),
#               .(x, y)]
#   centres.ctr <-
#     prob.grad[sample(1:nrow(prob.grad), nuclei, prob = 1 - prob.grad$value),
#               .(x, y)]
#   x.scale.trt <- runif(nuclei,
#                      x.scale.range[1],
#                      x.scale.range[2])
#   x.scale.ctr <- runif(nuclei,
#                      x.scale.range[1],
#                      x.scale.range[2])
#   y.scale.trt <- runif(nuclei,
#                      y.scale.range[1],
#                      y.scale.range[2])
#   y.scale.ctr <- runif(nuclei,
#                      y.scale.range[1],
#                      y.scale.range[2])
#   trt.eff <- runif(nuclei,
#                    nuc.eff.range[1] * effect.size.sp,
#                    nuc.eff.range[2] * effect.size.sp)
#   ctr.eff <- runif(nuclei,
#                    nuc.eff.range[1] * effect.size.sp,
#                    nuc.eff.range[2] * effect.size.sp)
#   for(i in 1:nuclei) {
#     trt.grad <-
#       generate_exp_gradient(x.dim = x.dim,
#                             y.dim = y.dim,
#                             centre = centres.trt[i, c(x, y)],
#                             x.scale = x.scale.trt[i],
#                             y.scale = y.scale.trt[i])
#     ctr.grad <-
#       generate_exp_gradient(x.dim = x.dim,
#                             y.dim = y.dim,
#                             centre = centres.ctr[i, c(x, y)],
#                             x.scale = x.scale.ctr[i],
#                             y.scale = y.scale.ctr[i])
#     treatment <- trt.eff[i] * trt.grad - ctr.eff[i] * ctr.grad + treatment
#   }
#   eff.scale <-
#     effect.size.sp /
#     (mean(treatment[poly[2,]][[1]], na.rm = TRUE) -
#      mean(treatment[poly[1,]][[1]], na.rm = TRUE))
#   treatment <- treatment * eff.scale
#   # eff.shift <- mean(treatment[poly[1,]][[1]], na.rm = TRUE)
#   # treatment <- treatment - eff.shift
#   poly.bd <-
#     st_geometry(poly[2,]) |>
#     st_as_sf()
#   poly.bd$eff.bd <- effect.size.bd
#   treatment.bd <- st_rasterize(poly.bd, template = generate_empty(x.dim = x.dim, y.dim = y.dim))
#   treatment <- treatment + treatment.bd
#   treatment <- treatment + (mean - mean(treatment[[1]], na.rm = TRUE))
#   treatment <- setNames(treatment, name)
#   return(treatment)
# }


# generate_treatment <- function(x.dim,
#                                y.dim,
#                                effect.size.sp,
#                                effect.size.bd,
#                                nuclei,
#                                poly,
#                                phi.range,
#                                shift.range,
#                                x.scale.range,
#                                y.scale.range,
#                                nuc.eff.range,
#                                grad.prop,
#                                grad.w,
#                                mat.nu,
#                                mat.var,
#                                mat.scale,
#                                mat.w,
#                                global.mean,
#                                name = "treatment") {
#   treatment <- generate_empty(x.dim = x.dim,
#                               y.dim = y.dim)
#   centroids <-
#     st_centroid(poly[1:2,]) |>
#     st_coordinates() |>
#     suppressWarnings()
#   phi.trt <-
#     atan2((centroids[2,2] - centroids[1,2]),
#           (centroids[2,1] - centroids[1,1]))
#   # prob.grad <-
#   #   generate_linear_gradient(x.dim = x.dim,
#   #                            y.dim = y.dim,
#   #                            phi = phi.trt,
#   #                            rescale = c(1-grad.prop, 1)) |>
#   # as.data.table()
#   grad.sig <-
#       generate_sigmoid_gradient(x.dim = x.dim,
#                                 y.dim = y.dim,
#                                 phi = phi.trt,
#                                 rescale = c(0, 1))
#   mat <-
#     generate_matern(x.dim = x.dim,
#                     y.dim = y.dim,
#                     nu = mat.nu,
#                     var = mat.var,
#                     scale = mat.scale,
#                     rescale = c(0,1))
#   treatment <- grad.w * grad.sig + mat.w * mat
#   eff.scale <-
#     effect.size.sp /
#     (mean(treatment[poly[2,]][[1]], na.rm = TRUE) -
#      mean(treatment[poly[1,]][[1]], na.rm = TRUE))
#   treatment <- treatment * eff.scale
#   # eff.shift <- mean(treatment[poly[1,]][[1]], na.rm = TRUE)
#   # treatment <- treatment - eff.shift
#   poly.bd <-
#     st_geometry(poly[2,]) |>
#     st_as_sf()
#   poly.bd$eff.bd <- effect.size.bd
#   treatment.bd <- st_rasterize(poly.bd, template = generate_empty(x.dim = x.dim, y.dim = y.dim))
#   treatment <- treatment + treatment.bd
#   treatment <- treatment + (global.mean - mean(treatment[[1]], na.rm = TRUE))
#   treatment <- setNames(treatment, name)
#   return(treatment)
# 


generate_nonlinear_effect <- function(field,
                                      eff.type = "random",
                                      eff.range = c(-1, 1),
                                      mu = 0,
                                      f.acc = 0.01) {
  field.rs <- scale_int(field, int = c(-1, 1))
  if(eff.type == "random") {
    eff.type <- sample(c("sigmoid", "minimum", "unimodal", "bimodal"), 1)
  }
  if(eff.type == "sigmoid") {
    par <- runif(1, 1, 10)
    fn <- function(x) 1/(1+exp(par*(-x)))
  }
  if(eff.type == "minimum") {
    par <- runif(1, -3/4, -1/4) * pi
    # par <- c(0.25, 1)
    fn <- function(x) sin(0.5*pi * (x) + par[1])
  }
  if(eff.type == "unimodal") {
    par <- c(runif(1, exp(1)/4, exp(1)), runif(1, -0.5, 0.5))
    fn <- function(x) exp(-((par[1] * (x - par[2]))^2))
  }
  if(eff.type == "bimodal") {
    par <- c(
             runif(1, 1, 2),
             runif(1, 1, 2),
             runif(1, exp(1), 2*exp(1)),
             runif(1, exp(1), 2*exp(1)),
             runif(1, -0.75, -0.5),
             runif(1, 0.5, 0.75),
             runif(1, 0.5, 1) * sample(c(-1, 1), 1))
    fn <- function(x) {
      par[1] * exp(-((par[3] * (x - par[5]))^2)) +
      par[2] * exp(-((par[4] * (x - par[6]))^2)) +
      par[7] * x
    }
  }
  effect <-
    fn(field.rs) |>
    scale_int(int = c(0, 1))
  effect <- eff.range * effect
  effect <- effect - mean(effect[[1]], na.rm = TRUE) + mu
  fun <- data.table(x = seq(0, 1, f.acc))
  fun[, f.x := eff.range * scale_int(fn((2 * x) - 1), int = c(0, 1))]
  fun[, f.x := f.x - mean(f.x, na.rm = TRUE) + mu]
  setnames(fun, c("x", "f.x"), c(names(field), paste0("f.", names(field)))) 
  return(list(effect = effect,
              fun = fun))
}


generate_interaction_effect <- function(field1,
                                        field2,
                                        range,
                                        mu,
                                        nuclei,
                                        f.acc = 0.01) {
  field1 <- scale_int(field1, int = c(0, 1))
  field2 <- scale_int(field2, int = c(0, 1))
  centres <- data.table(x = runif(nuclei, 0, 1),
                        y = runif(nuclei, 0, 1),
                        max = runif(nuclei, 0.5, 2),
                        sign = sample(c(-1, 1), nuclei, replace = TRUE),
                        scale.x = runif(nuclei, 0.5 * pi, 2 * pi))
  centres[, scale.y := scale.x * runif(nuclei, 0.5, 2)]
  fn <- function(x, y, centre.x, centre.y, sign, max, scale.x, scale.y) {
    sign * max * exp(-((scale.x * (x - centre.x))^2 + (scale.y * (y - centre.y))^2))
  }
  int.effect <- generate_empty(x.dim = dim(field1)[1],
                               y.dim = dim(field2)[2],
                               name = paste0("f.", names(field1), names(field2)))
  int.fun <- CJ(x = seq(0, 1, f.acc), y = seq(0, 1, f.acc))
  int.fun.com <- matrix(NA, nrow = nrow(int.fun), ncol = nuclei)
  for(i in 1:nrow(centres)){
    int.effect <-
      int.effect +
      fn(field1, field2,
         centre.x = centres[i, x],
         centre.y = centres[i, y],
         max = centres[i, max],
         sign = centres[i, sign],
         scale.x = centres[i, scale.x],
         scale.y = centres[i, scale.y])
    int.fun.com[,i] <-
      fn(int.fun$x, int.fun$y,
         centre.x = centres[i, x],
         centre.y = centres[i, y],
         max = centres[i, max],
         sign = centres[i, sign],
         scale.x = centres[i, scale.x],
         scale.y = centres[i, scale.y])
  }
  int.effect <- scale_int(int.effect, c(0, 1)) * range
  int.effect <- int.effect - mean(int.effect[[1]]) + mu
  int.fun[, f.xy := apply(int.fun.com, 1, sum)]
  setnames(int.fun, c(names(field1), names(field2),
                      paste0("f.", names(field1), names(field2))))
  return(list(effect = int.effect, fun = int.fun))
}


generate_landscape_4cov_nl <-
  function(seed = NULL,
           x.dim,
           y.dim,
           treatment.effect.size.sp,
           treatment.effect.size.bd,
           z1.effect.type,
           z1.effect.range,
           z1.effect.mu,
           z2.effect.type,
           z2.effect.range,
           z2.effect.mu,
           z3.effect.type,
           z3.effect.range,
           z3.effect.mu,
           z4.effect.type,
           z4.effect.range,
           z4.effect.mu,
           int12.effect.range,
           int12.effect.mu,
           int12.effect.nuclei,
           int13.effect.range,
           int13.effect.mu,
           int13.effect.nuclei,
           int14.effect.range,
           int14.effect.mu,
           int14.effect.nuclei,
           int23.effect.range,
           int23.effect.mu,
           int23.effect.nuclei,
           int24.effect.range,
           int24.effect.mu,
           int24.effect.nuclei,
           int34.effect.range,
           int34.effect.mu,
           int34.effect.nuclei,
           treatment.nuclei,
           treatment.phi.range,
           treatment.shift.range,
           treatment.x.scale.range,
           treatment.y.scale.range,
           treatment.nuc.eff.range,
           z1.fbm.alpha,
           z1.fbm.var,
           z1.fbm.scale,
           z1.fbm.w,
           z1.grad.phi,
           z1.grad.w,
           z2.fbm.alpha,
           z2.fbm.var,
           z2.fbm.scale,
           z2.fbm.w,
           z2.grad.phi,
           z2.grad.shift,
           z2.grad.w,
           z3.dist.n,
           z3.grad.phi,
           z3.grad.prop,
           z3.acc,
           z4.seg.n,
           z4.mat.nu,
           z4.mat.var,
           z4.mat.scale,
           z4.mat.w,
           z4.grad.phi,
           z4.grad.w,
           split.n,
           split.prop,
           e.exp.var,
           e.exp.scale,
           e.nug.var,
           e.rand.var,
           ...) {
  ls.par <- as.list(match.call(expand.dots=FALSE))
  set.seed(seed)
  z1 <- generate_z1(x.dim = x.dim,
                    y.dim = y.dim,
                    fbm.alpha = z1.fbm.alpha,
                    fbm.var = z1.fbm.var,
                    fbm.scale = z1.fbm.scale,
                    fbm.w = z1.fbm.w,
                    grad.phi = z1.grad.phi,
                    grad.w = z1.grad.w,
                    rescale = c(0,1),
                    name = "z1")

  z2 <- generate_z2(
                    x.dim = x.dim,
                    y.dim = y.dim,
                    fbm.alpha = z2.fbm.alpha,
                    fbm.var = z2.fbm.var,
                    fbm.scale = z2.fbm.scale,
                    fbm.w = z2.fbm.w,
                    grad.phi = z2.grad.phi,
                    grad.shift = z2.grad.shift,
                    grad.w = z2.grad.w,
                    rescale = c(0,1),
                    name = "z2")

  z3 <- generate_z3(x.dim = x.dim,
                    y.dim = y.dim,
                    dist.n = z3.dist.n,
                    grad.phi = z3.grad.phi,
                    grad.prop = z3.grad.prop,
                    acc = z3.acc,
                    rescale = c(0,1),
                    name = "z3")

  z4 <- generate_z4(x.dim = x.dim,
                    y.dim = y.dim,
                    seg.n = z4.seg.n,
                    mat.nu = z4.mat.nu,
                    mat.var = z4.mat.var,
                    mat.scale = z4.mat.scale,
                    mat.w = z4.mat.w,
                    grad.phi = z4.grad.phi,
                    grad.w = z4.grad.w,
                    rescale = c(0,1),
                    name = "z4")

  split <- generate_split(x.dim = x.dim, 
                          y.dim = y.dim,
                          n = split.n,
                          prop = split.prop,
  )
  type <- st_rasterize(split, template = generate_empty(x.dim = x.dim, y.dim = y.dim))
  type <- reset_dim(type)
  type <- setNames(type, "type")

  treatment <-
    generate_treatment(x.dim = x.dim,
                       y.dim = y.dim,
                       effect.size.sp = treatment.effect.size.sp,
                       effect.size.bd = treatment.effect.size.bd,
                       nuclei = treatment.nuclei,
                       poly = split,
                       phi.range = treatment.phi.range,
                       shift.range = treatment.shift.range,
                       x.scale.range = treatment.x.scale.range,
                       y.scale.range = treatment.y.scale.range,
                       nuc.eff.range = treatment.nuc.eff.range,
                       name = "treatment")

  covariates <- c(z1, z2, z3, z4)
  
  cov.main.names <- paste0("z", 1:4)
  cov.main.effect.type <- ls.par[paste0(cov.main.names, ".effect.type")]
  cov.main.effect.range <- ls.par[paste0(cov.main.names, ".effect.range")]
  cov.main.effect.mu <- ls.par[paste0(cov.main.names, ".effect.mu")]
  cov.main.effects <- list()
  cov.main.funs <- list()
  for(i in seq_along(cov.main.names)){
    nl.effect <-
      generate_nonlinear_effect(covariates[i],
                                eff.type = cov.main.effect.type[[i]],
                                eff.range = cov.main.effect.range[[i]],
                                mu = cov.main.effect.mu[[i]])
    cov.main.effects[[i]] <- nl.effect$effect
    setnames(nl.effect$fun, c("val", "f"))
    nl.effect$fun[, cov := cov.main.names[i]]
    cov.main.funs[[i]] <- nl.effect$fun
  }
  cov.main.effects <-
    do.call(c, cov.main.effects) |>
    setNames(paste0("f.", cov.main.names))
  cov.main.funs <- rbindlist(cov.main.funs)
  setcolorder(cov.main.funs, c("cov", "val", "f"))

  cov.comb <- t(combn(1:4, 2))
  cov.int.names <- apply(cov.comb, 1, \(x) paste0("int", x[1], x[2]))
  cov.int.effect.range <- ls.par[paste0(cov.int.names, ".effect.range")]
  cov.int.effect.mu <- ls.par[paste0(cov.int.names, ".effect.mu")]
  cov.int.effect.nuclei <- ls.par[paste0(cov.int.names, ".effect.nuclei")]
  cov.int.effects <-list()
  cov.int.funs <- list()
  for(i in seq_along(cov.int.names)){
    int.effect <-
      generate_interaction_effect(field1 = get(cov.main.names[cov.comb[i, 1]]),
                                  field2 = get(cov.main.names[cov.comb[i, 2]]),
                                  range = cov.int.effect.range[[i]],
                                  mu = cov.int.effect.mu[[i]],
                                  nuclei = cov.int.effect.nuclei[[i]])
    cov.int.effects[[i]] <- int.effect$effect
    setnames(int.effect$fun, c("val1", "val2", "f"))
    int.effect$fun[, `:=`(int = cov.int.names[i],
                          cov1 = cov.main.names[cov.comb[i, 1]],
                          cov2 = cov.main.names[cov.comb[i, 2]])]
    cov.int.funs[[i]] <- int.effect$fun
  }
  cov.int.effects <-
    do.call(c, cov.int.effects) |>
    setNames(paste0("f.", cov.int.names))
  cov.int.funs <- rbindlist(cov.int.funs)
  setcolorder(cov.int.funs, c("int", "cov1", "cov2", "val1", "val2", "f"))

  mod.error <-
    RMexp(var = e.exp.var, scale = e.exp.scale) +
    RMnugget(var = e.nug.var)
  error.sp <- RFsimulate(mod.error, x = 1:x.dim, y = 1:y.dim)
  error.sp <- reset_dim(st_as_stars(error.sp))
  error.sp <- setNames(error.sp, "error")
  error.rand <-
    matrix(rnorm(prod(x.dim, y.dim), 0, e.rand.var), nrow = x.dim, ncol = y.dim) |>
    matrix2stars()
  error <- error.sp + error.rand

  response <-
    c(cov.main.effects, treatment, error) |>
    merge() |>
    st_apply(1:2, sum) |>
    setNames("response")

  response.int <-
    c(cov.main.effects, cov.int.effects, treatment, error) |>
    merge() |>
    st_apply(1:2, sum) |>
    setNames("response.int")

  landscape <- c(response, response.int, type, covariates, treatment, cov.main.effects, cov.int.effects, error)

  landscape.dt <- 
    as.data.frame(landscape)
  setDT(landscape.dt)
  landscape.dt[, type := factor(ifelse(type == 0, "control", "treatment"),
                        levels = c("control", "treatment"))]
  landscape.dt$cell <- 1:nrow(landscape.dt)
  setcolorder(landscape.dt, c("x", "y", "cell"))

  fun <- list(main = cov.main.funs, interactions = cov.int.funs)
  
  return(list(landscape = landscape.dt, fun = fun))
  }


zplot <- function(z, n.colours = 128) {
  plot(z, nbreaks = n.colours + 1, breaks = "equal", col = hcl.colors(n.colours))
}


assign_bl_som <- function(data,
                          som,
                          cov.col,
                          id.col = "id",
                          bl.col = "type",
                          bl.level = "control",
                          som.col = "som_bmu",
                          min.obs = 1,
                          n.threads = 4,
                          scale = TRUE
                          ){
  data <- setnames(copy(data), c(id.col, bl.col, som.col), c("id", "type", "som_bmu"))
  # Baseline units
  bl.bmu <- data[type == bl.level, .(n = .N), som_bmu]
  # If BMU contains a sufficient number of baseline observations determined by
  # `min.obs`, it is set to be the baseline BMU.
  data.bmu.bl <-
    copy(data)[som_bmu %in% bl.bmu[n >= min.obs, som_bmu]
                    ][, som_bmu.bl := as.list(som_bmu)]
  # Observations to which a different (or multiple, if min.obs > 1) baseline
  # BMU has to be assigned because there are not enough observations in the
  # first-order BMU.
  data.bmu.bl.mult <-
    data[!som_bmu %in% bl.bmu[n >= min.obs, som_bmu]]
  if(scale) {
    data.bmu.bl.mult$som_bmu.bl.mult <-
      scale_data_som(data.bmu.bl.mult[, ..cov.col],
                     som) |>
      bmu_match_reference(som, bl.bmu[, .(som_bmu, n)], min.obs, n.threads)
  } else {
    data.bmu.bl.mult$som_bmu.bl.mult <-
      as.matrix(data.bmu.bl.mult[, ..cov.col]) |>
      bmu_match_reference(som, bl.bmu[, .(som_bmu, n)], min.obs, n.threads)
  }
  data.bmu <-
    data |>
    merge(data.bmu.bl[, .(id, som_bmu.bl)], all = TRUE, sort = FALSE) |>
    merge(data.bmu.bl.mult[, .(id, som_bmu.bl.mult)], all = TRUE, sort = FALSE)
  data.bmu[, som_bmu.bl.mult := as.list(som_bmu.bl.mult)]
  data.bmu[, `:=`(som_bmu.bl = fifelse(unlist(lapply(som_bmu.bl, is.null)),
                                       som_bmu.bl.mult, som_bmu.bl))]
  # Weights for each BMU
  id.bl <-
    data.bmu[,
             .(id, som_bmu.bl)
             ][, lapply(.SD, \(x) as.numeric(unlist(x))), by = id] |>
    merge(bl.bmu[, .(som_bmu.bl = som_bmu, n.bmu = n)], by = "som_bmu.bl", all = TRUE, sort = FALSE)
  id.bl[, n.id := sum(n.bmu), by = id][, w := n.bmu/n.id]
  data.bl <- id.bl[,.(som_bmu.bl = list(som_bmu.bl), som_bmu.bl.w = list(w)), id]
  return(data.bl)
}




generate_landscape_4cov_lin <-
  function(seed = NULL,
           x.dim,
           y.dim,
           treatment.effect.size.sp,
           treatment.effect.size.bd,
           z1.slope,
           z1.mean,
           z2.slope,
           z2.mean,
           z3.slope,
           z3.mean,
           z4.slope,
           z4.mean,
           # z1.effect.size,
           # z2.effect.size,
           # z3.effect.size,
           # z4.effect.size,
           treatment.nuclei,
           treatment.phi.range,
           treatment.shift.range,
           treatment.x.scale.range,
           treatment.y.scale.range,
           treatment.nuc.eff.range,
           treatment.grad.prop,
           treatment.grad.w,
           treatment.mat.nu,
           treatment.mat.var,
           treatment.mat.scale,
           treatment.mat.w,
           treatment.global.mean,
           z1.fbm.alpha,
           z1.fbm.var,
           z1.fbm.scale,
           # z1.fbm1.alpha,
           # z1.fbm1.var,
           # z1.fbm1.scale,
           # z1.fbm2.alpha,
           # z1.fbm2.var,
           # z1.fbm2.scale,
           # z1.fbm.ratio,
           z1.fbm.w,
           z1.grad.phi,
           z1.grad.w,
           z2.fbm.alpha,
           z2.fbm.var,
           z2.fbm.scale,
           z2.fbm.w,
           z2.grad.phi,
           z2.grad.shift,
           z2.grad.w,
           z3.dist.n,
           z3.grad.phi,
           z3.grad.prop,
           z3.acc,
           z4.seg.n,
           z4.mat.nu,
           z4.mat.var,
           z4.mat.scale,
           z4.mat.w,
           z4.grad.phi,
           z4.grad.w,
           split.n,
           split.prop,
           e.exp.var,
           e.exp.scale,
           e.nug.var,
           e.rand.var,
           ...) {
    set.seed(seed)

  z1 <- generate_z1(x.dim = x.dim,
                    y.dim = y.dim,
                    fbm.alpha = z1.fbm.alpha,
                    fbm.var = z1.fbm.var,
                    fbm.scale = z1.fbm.scale,
                    fbm.w = z1.fbm.w,
                    grad.phi = z1.grad.phi,
                    grad.w = z1.grad.w,
                    rescale = c(0,1),
                    name = "z1")

  z2 <- generate_z2(
                    x.dim = x.dim,
                    y.dim = y.dim,
                    fbm.alpha = z2.fbm.alpha,
                    fbm.var = z2.fbm.var,
                    fbm.scale = z2.fbm.scale,
                    fbm.w = z2.fbm.w,
                    grad.phi = z2.grad.phi,
                    grad.shift = z2.grad.shift,
                    grad.w = z2.grad.w,
                    rescale = c(0,1),
                    name = "z2")

  z3 <- generate_z3(x.dim = x.dim,
                    y.dim = y.dim,
                    dist.n = z3.dist.n,
                    grad.phi = z3.grad.phi,
                    grad.prop = z3.grad.prop,
                    acc = z3.acc,
                    rescale = c(0,1),
                    name = "z3")

  z4 <- generate_z4(x.dim = x.dim,
                    y.dim = y.dim,
                    seg.n = z4.seg.n,
                    mat.nu = z4.mat.nu,
                    mat.var = z4.mat.var,
                    mat.scale = z4.mat.scale,
                    mat.w = z4.mat.w,
                    grad.phi = z4.grad.phi,
                    grad.w = z4.grad.w,
                    rescale = c(0,1),
                    name = "z4")

  split <- generate_split(x.dim = x.dim, 
                          y.dim = y.dim,
                          n = split.n,
                          prop = split.prop,
  )
  type <- st_rasterize(split, template = generate_empty(x.dim = x.dim, y.dim = y.dim))
  type <- reset_dim(type)
  type <- setNames(type, "type")

  treatment <-
    generate_treatment(x.dim = x.dim,
                       y.dim = y.dim,
                       effect.size.sp = treatment.effect.size.sp,
                       effect.size.bd = treatment.effect.size.bd,
                       nuclei = treatment.nuclei,
                       poly = split,
                       phi.range = treatment.phi.range,
                       shift.range = treatment.shift.range,
                       x.scale.range = treatment.x.scale.range,
                       y.scale.range = treatment.y.scale.range,
                       nuc.eff.range = treatment.nuc.eff.range,
                       name = "treatment")

  # treatment <-
  #   generate_treatment(x.dim = x.dim, y.dim = y.dim,
  #                      effect.size.sp = treatment.effect.size.sp,
  #                      effect.size.bd = treatment.effect.size.bd,
  #                      nuclei = treatment.nuclei,
  #                      poly = split,
  #                      phi.range = treatment.phi.range,
  #                      shift.range = treatment.shift.range,
  #                      x.scale.range = treatment.x.scale.range,
  #                      y.scale.range = treatment.y.scale.range,
  #                      nuc.eff.range = treatment.nuc.eff.range,
  #                      grad.prop = treatment.grad.prop,
  #                      grad.w = treatment.grad.w,
  #                      mat.nu = treatment.mat.nu,
  #                      mat.var = treatment.mat.var,
  #                      mat.scale = treatment.mat.scale,
  #                      mat.w = treatment.mat.w,
  #                      global.mean = treatment.global.mean,
  #                      name = "treatment")


  covariates <- c(z1, z2, z3, z4)

  # cov.effect.size <- c(z1.effect.size, z2.effect.size, z3.effect.size, z4.effect.size)
  # cov.effects <- list()
  # for(i in seq_along(cov.effect.size)){
  #   cov.effects[[i]] <-
  #     matrix2stars(covariates[[i]] * cov.effect.size[i])
  # }

  cov.mean <- c(z1.mean, z2.mean, z3.mean, z4.mean)
  cov.slope <- c(z1.slope, z2.slope, z3.slope, z4.slope)
  cov.effects <- list()
  for(i in seq_along(cov.mean)){
    mu <- mean(covariates[[i]], na.rm = TRUE)
    cov.effects[[i]] <-
      matrix2stars((covariates[[i]] * cov.slope[i]) + (cov.mean[i] - mu))
  }

  cov.effects <-
    do.call(c, cov.effects) |>
    setNames(paste0("f", 1:4))

  mod.error <-
    RMexp(var = e.exp.var, scale = e.exp.scale) +
    RMnugget(var = e.nug.var)
  error.sp <- RFsimulate(mod.error, x = 1:x.dim, y = 1:y.dim)
  error.sp <- reset_dim(st_as_stars(error.sp))
  error.sp <- setNames(error.sp, "error")
  error.rand <-
    matrix(rnorm(prod(x.dim, y.dim), 0, e.rand.var), nrow = x.dim, ncol = y.dim) |>
    matrix2stars()
  error <- error.sp + error.rand

  response <-
    c(cov.effects, treatment, error) |>
    merge() |>
    st_apply(1:2, sum) |>
    setNames("response")

  landscape <- c(response, type, covariates, treatment, cov.effects, error)

  landscape.dt <- 
    as.data.frame(landscape)
  setDT(landscape.dt)
  landscape.dt[, type := factor(ifelse(type == 0, "control", "treatment"),
                        levels = c("control", "treatment"))]
  landscape.dt$cell <- 1:nrow(landscape.dt)
  setcolorder(landscape.dt, c("x", "y", "cell"))
  
  return(landscape.dt)
  }


zplot <- function(z, n.colours = 128) {
  plot(z, nbreaks = n.colours + 1, breaks = "equal", col = hcl.colors(n.colours))
}


plot_landscape_4cov_lin <- function(x) {

  ls_theme <-
    theme_void(base_family = "IBMPlexSans", base_size = 7) +
    theme(plot.title = element_text(hjust = 0.5, margin = margin(l = 3, b = 3, t = 3)),
          plot.margin = margin(3, 3, 3, 3),
          legend.position = "right",
          legend.justification = c(0,1),
          legend.title = element_text(size = rel(0.9), hjust = 0,
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
                                 barwidth = 0.5,
                                 barheight = 5,
                                 label.position = "left",
                                 label.hjust = 1,
                                 draw.ulim = TRUE,
                                 draw.llim = TRUE
                                 ))

  
  lab.coord <- x[, lapply(.SD, mean), .SDcols = c("x", "y"), by = type]

  plots <- list()

  plots[["type"]] <- 
    ggplot(x, aes(x = x, y = y)) +
    geom_raster(aes(fill = type)) +
    geom_text(data = lab.coord,
              aes(x = x, y = y,
                  label = capwords(as.character(type)), color = type), size = 3)  +
    scale_fill_manual(breaks = levels(x$type),
                      values = c("grey20", "grey80"),
                      labels = capwords(levels(x$type))) +
    scale_colour_manual(breaks = levels(x$type),
                        values = c("grey95", "grey5")) +
    # scale_x_continuous(expand = c(0, 0)) +
    # scale_y_continuous(expand = c(0, 0)) +
    coord_fixed(expand = FALSE) +
    guide_fill +
    labs(title = "Area type") +
    guides(fill = "none", color = "none") +
    ls_theme

  cov.var <- paste0("z", 1:4)

  for(i in 1:4) {
    plots[[cov.var[i]]] <-
      x[, .(x, y, value = eval(parse(text = cov.var[i])))] |>
      ggplot() +
      geom_raster(aes(x = x, y = y, fill = value)) +
      scale_fill_continuous_sequential(palette = "Viridis", rev = FALSE, limits = c(0, 1)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed() +
      guide_fill +
      labs(title = cov.var[i], fill = "Covariate\nvalue") +
      ls_theme
  }

  eff.var <-  paste0("f.", cov.var)

  lim.effects <-
    unlist(x[,
             .(min(as.matrix(.SD)), max(as.matrix(.SD))),
             .SDcols = c("treatment", "error", paste0("f.", cov.var))],
           use.names = FALSE)
  lim.effects <- 
    (floor(lim.effects / 0.25) + c(0, 1)) * 0.25

  for(i in seq_along(eff.var)) {
    plots[[eff.var[i]]] <-
      x[, .(x, y, value = eval(parse(text = eff.var[i])))] |>
      ggplot(aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                       limits = lim.effects) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed() +
      guide_fill +
      labs(title = eff.var[i], fill = "Effect size") +
      ls_theme
  }

  plots[["treatment"]] <-
    ggplot(x, aes(x = x, y = y, fill = treatment)) +
    geom_raster() +
    scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                     limits = lim.effects) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    guide_fill +
    labs(title = "Treatment effect", fill = "Effect size") +
    ls_theme

  plots[["error"]] <-
    ggplot(x, aes(x = x, y = y, fill = error)) +
    geom_raster() +
    scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                     limits = lim.effects) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    guide_fill +
    labs(title = "Residual variation (error)", fill = "Effect size") +
    ls_theme

  plots[["response"]] <-
    ggplot(x, aes(x = x, y = y, fill = response)) +
    geom_raster() +
    # scale_fill_viridis_c()
    scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    guide_fill +
    labs(title = "Response", fill = "Observed\nvalue") +
    ls_theme

  layout1 <-
    "ABC
     #DE"

  layout2 <-
    "ABC
     FDE"



  p1 <- 
  plots$type + plots$z1 + plots$z2 + plots$z3 + plots$z4 +
  plot_layout(design = layout1, guides = "collect") &
  ls_theme

  p2 <-
  plots$treatment + plots$f.z1 + plots$f.z2 + plots$f.z3 + plots$f.z4 + plots$error +
  plot_layout(design = layout2, guides = "collect") &
  ls_theme
  
  p3 <-
  plots$response + plot_layout(design = "A##", guides = "keep") +
  ls_theme

  pw <- wrap_plots(p1, p2, p3, ncol = 1, heights = c(2.06, 2.06, 1))

  return(pw)
}



capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
                  {s <- substring(s, 2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}



