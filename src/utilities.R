library(mgcv)
library(posterior)
library(RandomFields)
library(RandomFieldsUtils)
library(raster)
library(sf)
library(stars)
library(lwgeom)
library(spdep)
library(igraph)
library(data.table)
library(mvnfast)
library(colorspace)
library(stringi)
library(ggplot2)
library(patchwork)
library(kohonen)
library(GA)
library(memoise)
library(Matrix)

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
    # For loop is *much* faster than merge, st_apply, split
    for(i in 1:n.att) {
      min.att <- min(x[[i]], na.rm = TRUE)
      max.att <- max(x[[i]], na.rm = TRUE)
      x[[i]] <- int[1] + (int[2] - int[1]) * (x[[i]] - min.att) / (max.att - min.att)
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


scale.stars <- function(x, center = TRUE, scale = TRUE) {
  if(!is.null(int)) {
    x.d <- as.data.table(as.data.frame(x))
    att <- names(x)
    x <-
      x.d[, (att) := lapply(.SD, \(x) scale(x, center = center, scale = scale)),
          .SDcols = att] |>
      st_as_stars()
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


inverse <- function(x, ...) {
  UseMethod("inverse", x)
}


inverse.stars <- function(x, attributes = NULL) {
  if(is.null(attributes)) {
    inv.att <- 1:length(x)
  } else {
    if(is.logical(attributes)) {
      inv.att <- which(attributes)
    }
    if(is.character(attributes)) {
      inv.att <- which(names(x) %in% attributes)
    }
  }
  for(i in inv.att) {
    x[[i]] <- x[[i]] * -1
  }
  return(x)
}


score <- function(x, ...) {
  UseMethod("score", x)
}

score.stars <- function(x,
                        type = "euclidean") {
  att <- names(x)
  x.d <- as.data.table(as.data.frame(x))
  att.mat <- as.matrix(x.d[, ..att])
  att.min <- apply(att.mat, 2, min)
  if(type == "qss")
    x.d[, score := apply(apply(att.mat, 2, \(x) ecdf(x)(x)),
                         1,  \(x) sum(x^2))]
  if(type == "qprod")
    x.d[, score := apply(apply(att.mat, 2, \(x) ecdf(x)(x)),
                         1, \(x) prod(x)^(1/length(att)))]
  if(type == "mahalanobis")
    x.d[, score := mahalanobis(att.mat, center = att.min, cov = cov(att.mat))]
  if(type == "sumofsquares")
    x.d[, score := apply(att.mat, 1, \(x) sum((x - att.min)^2))]
  if(type == "euclidean")
    x.d[, score := apply(att.mat, 1, \(x) sqrt(sum((x - att.min)^2)))]
  if(type == "prod")
    x.d[, score := apply(att.mat, 1, \(x) prod(x - att.min))]
  x.s <- st_as_stars(x.d[, -..att])
  return(x.s)
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

resize_polygon <- function(x,
                           limit,
                           dist,
                             ...) {
  st_agr(x) <- "constant"
  x.buf <-
    st_buffer(x, dist = dist, ...) |>
    st_intersection(limit) |>
    st_as_sf() |>
    st_set_agr("constant")
  # poly.new <-
    # res.grid[unlist(st_contains(poly.buf, st_centroid(grid))),]
  # return(poly.new)
  return(x.buf[, names(x)])
}


resize_polygons <- function(x,
                            limit,
                            dist,
                             ...) {
  st_agr(x) <- "constant"
  x.buf.l <- list()
  for(i in 1:nrow(x)) {
    x.buf.l[[i]] <-
      resize_polygon(x[i,],
                     limit = limit[i,],
                     dist = dist,
                     ...)
  }
  x.buf <-
    do.call(rbind, x.buf.l) |>
    st_as_sf()
  return(x.buf[, names(x)])
}


d_cohen <- function(x, y) {
  mu_x <- mean(x)
  mu_y <- mean(y)
  s <- sqrt((sum((x-mu_x)^2) + sum((y-mu_y)^2)) / (length(x) + length(y) - 2))
  d <- (mu_x - mu_y) / s
  return(d)
}

var_ratio <- function(x, y) {
  var(x) / var(y)
}

ks_stat <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  n <- length(x)
  y <- y[!is.na(y)]
  n.x <- as.double(n)
  n.y <- length(y)
  w <- c(x, y)
  z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
  z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)] #exclude ties
  ks <- max(abs(z))
  return(ks)
}

get_breaks <- function(x, breaks = "FD", ...) {
  x.range <- range(x)
  if(is.character(breaks)) {
    breaks <- match.arg(tolower(breaks),
                        c("fd", "scott", "sturges"))
    breaks <- switch(breaks,
                     fd = nclass.FD(x),
                     scott = nclass.scott(x),
                     sturges = nclass.Sturges(x),
                     stop(paste("Unknown algorithm to calculate breaks.",
                                "Possible values are: `FD`, `Scott`, or `Sturges`.")))
  }
  if(is.numeric(breaks) & length(breaks) == 1) {
    breaks <- seq(x.range[1], x.range[2], length.out = breaks+1)
  }
  return(breaks)
}

assign_breaks <- function(x, breaks = "FD", ...) {
  breaks <- get_breaks(x, breaks = breaks)
  breaks.assigned = cut(x, breaks, include.lowest = TRUE, labels = FALSE)
  return(breaks.assigned)
}


to_pdist <- function(x, ...) {
    p <- density(x, ...)$y
    p <- p/sum(p)
    return(p)
}


to_pmass <- function(x, breaks = "FD", ...) {
  # p.range <- range(x)
  # if(is.character(breaks)) {
  #   breaks <- match.arg(tolower(breaks),
  #                       c("fd", "scott", "sturges"))
  #   breaks <- switch(breaks,
  #                    fd = nclass.FD(x),
  #                    scott = nclass.scott(x),
  #                    sturges = nclass.Sturges(x),
  #                    stop(paste("Unknown algorithm to calculate breaks.",
  #                               "Possible values are: `FD`, `Scott`, or `Sturges`.")))
  # }
  # if(is.numeric(breaks) & length(breaks) == 1) {
  #   breaks <- seq(p.range[1], p.range[2], length.out = breaks+1)
  # }
  x.dt <- data.table(x, h = assign_breaks(x, breaks = breaks))
  x.dt <- x.dt[, .(x = .N/nrow(x.dt)), h]
  if(nrow(x.dt) < length(breaks)) {
    i.br <- 1:length(breaks)
    x.dt <- rbind(x.dt, data.table(h = i.br[!(i.br %in% x.dt$h)], x = 0))
  }
  return(x.dt[is.na(x), x := 0][order(h), x])
}



# kl_div <- function(x, y, base = 2, symmetric = FALSE, to.pdist = FALSE, p.range = NULL, p.n = 512, p.bw = "SJ-dpi") {
#   if(to.pdist) {
#     if(is.null(p.range)) {
#       p.range <- range(c(range(x), range(y)))
#     }
#     if(is.null(p.n)) {
#       p.n <- 512
#     }
#     x <- to_pdist(x, bw = p.bw, from = p.range[1], to = p.range[2], n = 512)
#     y <- to_pdist(y, bw = p.bw, from = p.range[1], to = p.range[2], n = 512)
#   }
#   if(symmetric) {
#     d <-
#       0.5 *
#       (kl_dist(x, y, base = base, symmetric == FALSE, to.pdist = FALSE) /
#        kl_dist(y, x, base = base, symmetric == FALSE, to.pdist = FALSE))
#   } else {
#     d <- sum(x * log(x / y, base = base))
#   }
#   return(d)
# }


kl_div <- function(
                   x,
                   y,
                   base = 2,
                   symmetric = FALSE,
                   type = "raw",
                   disc.breaks = "FD",
                   dens.range = NULL,
                   dens.n = 512,
                   dens.bw = "SJ") {
  if(type == "density") {
    if(is.null(dens.range)) {
      p.range <- range(c(range(x), range(y)))
    }
    if(is.null(dens.n)) {
      dens.n <- 512
    }
    x <- to_pdist(x, bw = dens.bw, from = p.range[1], to = p.range[2], n = dens.n)
    y <- to_pdist(y, bw = dens.bw, from = p.range[1], to = p.range[2], n = dens.n)
  }
  if(type == "discrete") {
    disc.breaks <- get_breaks(c(x, y), disc.breaks = disc.breaks)
    x <- to_pmass(x, breaks = disc.breaks)
    y <- to_pmass(y, breaks = disc.breaks)
  }
  x[x == 0] <- .Machine$double.eps
  y[y == 0] <- .Machine$double.eps
  if(length(x) != length(y)) {
    print(x)
    print(y)
    # stop("Wrong length.")
  }
  if(symmetric) {
    d <-
      0.5 *
      (kl_dist(x, y, base = base, symmetric == FALSE, type = "raw") /
       kl_dist(y, x, base = base, symmetric == FALSE, type = "raw"))
  } else {
    d <- sum(x * log(x / y, base = base))
  }
  return(d)
}


# js_div <- function(x, y, base = 2, to.pdist = FALSE, p.range = NULL, p.n = 512, p.bw = "SJ") {
#   if(type == "density") {
#     if(is.null(dens.range)) {
#       p.range <- range(c(range(x), range(y)))
#     }
#     if(is.null(dens.n)) {
#       dens.n <- 512
#     }
#     x <- to_pdist(x, bw = dens.bw, from = p.range[1], to = p.range[2], n = dens.n)
#     y <- to_pdist(y, bw = dens.bw, from = p.range[1], to = p.range[2], n = dens.n)
#   }
#   if(type == "discrete") {
#     p.range <- range(c(range(x), range(y)))
#     if(is.character(disc.breaks)) {
#       disc.breaks <- match.arg(tolower(disc.breaks),
#                                c("fd", "scott", "sturges"))
#       disc.breaks <- switch(disc.breaks,
#                             fd = nclass.FD(x),
#                             scott = nclass.scott(x),
#                             sturges = nclass.Sturges(x),
#                             stop(paste("Unknown algorithm to calculate breaks.",
#                                        "Possible values are: `FD`, `Scott`, or `Sturges`.")))
#     }
#     if(is.numeric(disc.breaks) & length(disc.breaks) == 1) {
#       disc.breaks <- seq(p.range[1], p.range[2], length.out = disc.breaks+1)
#     }
#     x <- to_pmass(x, breaks = disc.breaks)
#     y <- to_pmass(y, breaks = disc.breaks)
#   }
#   x[x == 0] <- .Machine$double.eps
#   y[y == 0] <- .Machine$double.eps
#   z <- 0.5 * (x + y)
#   d <-
#     0.5 *
#     (kl_div(x, z, base = base, symmetric = FALSE, type = "raw") +
#      kl_div(y, z, base = base, symmetric = FALSE, type = "raw"))
#   return(d)
# }

js_div <- function(x,
                   y,
                   base = 2,
                   type = "raw",
                   disc.breaks = "FD",
                   dens.range = NULL,
                   dens.n = 512,
                   dens.bw = "SJ") {
  if(type == "density") {
    if(is.null(dens.range)) {
      p.range <- range(c(range(x), range(y)))
    }
    if(is.null(dens.n)) {
      dens.n <- 512
    }
    x <- to_pdist(x, bw = dens.bw, from = p.range[1], to = p.range[2], n = dens.n)
    y <- to_pdist(y, bw = dens.bw, from = p.range[1], to = p.range[2], n = dens.n)
  }
  if(type == "discrete") {
    disc.breaks <- get_breaks(c(x, y), breaks = disc.breaks)
    x <- to_pmass(x, breaks = disc.breaks)
    y <- to_pmass(y, breaks = disc.breaks)
  }
  x[x == 0] <- .Machine$double.eps
  y[y == 0] <- .Machine$double.eps
  z <- 0.5 * (x + y)
  d <-
    0.5 *
    (kl_div(x, z, base = base, symmetric = FALSE, type = "raw") +
     kl_div(y, z, base = base, symmetric = FALSE, type = "raw"))
  return(d)
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


generate_areas_poly <- function(x.dim,
                                y.dim,
                                area.prop = 0.5,
                                area.exact = FALSE,
                                shape.constraints = TRUE,
                                shape.compactness = 2,
                                score = NULL,
                                imbalance = 0.1,
                                score.prop = 0.5,
                                score.sam = 1e4,
                                seg.n,
                                seg.res,
                                seg.min.dist,
                                min.bound.dist = 0,
                                verbose = FALSE,
                                iter.max = 25,
                                opt.tol = 1e-4,
                                ...
                                ) {

# for(m in 1:20){
  if(verbose) message("Setting up candidate polygons ...")

  if(is.null(score)) {
    score <- generate_empty(x.dim = x.dim, y.dim = y.dim) + 1
    score.prop <- 0
  } else {
    score <- score[1]
  }
  names(score)[1] <- "score"
  score <- scale_int(score, int = c(0, 1))

  poly.lim <- limit_polygon(x.dim, y.dim)

  # set up segment grid
  x.crd.seg <- seq(1, x.dim - (0.5 * seg.min.dist), seg.min.dist) + (0.5 * seg.min.dist)
  y.crd.seg <- seq(1, y.dim - (0.5 * seg.min.dist), seg.min.dist) + (0.5 * seg.min.dist)
  seg.grid <-
    cbind(x = rep(x.crd.seg, times = length(y.crd.seg)),
          y = rep(y.crd.seg, each = length(y.crd.seg))) |>
    st_multipoint() |>
    st_sfc() |>
    st_cast("POINT") |>
    st_as_sf()
  seg.sam <- sample(1:nrow(seg.grid), seg.n)
  seg.coarse <- 
      seg.grid[seg.sam,] |>
      st_union() |>
      st_voronoi() |>
      st_cast() |>
      st_intersection(poly.lim) |>
      st_as_sf() |>
      st_set_geometry("geometry")

  # set up segment grid
  x.crd.res <- seq(1, x.dim - (0.5 * seg.res), seg.res) + (0.5 * seg.res)
  x.crd.res <- x.crd.res + runif(n = length(x.crd.res), -seg.res/2, seg.res/2)
  y.crd.res <- seq(1, y.dim - (0.5 * seg.res), seg.res) + (0.5 * seg.res)
  y.crd.res <- y.crd.res + runif(n = length(y.crd.res), -seg.res/2, seg.res/2)

  res.grid <-
    (cbind(x = rep(x.crd.res, times = length(y.crd.res)),
           y = rep(y.crd.res, each = length(y.crd.res))) +
     matrix(runif(n = 2 * prod(length(x.crd.res), length(y.crd.res)),
                  -seg.res/2, seg.res/2),
            ncol = 2)) |>
    st_multipoint() |>
    st_sfc() |>
    st_cast("POINT") |>
    st_as_sf() |>
    st_union() |>
    st_voronoi() |>
    st_cast() |>
    st_intersection(poly.lim) |>
    st_as_sf() |>
    st_set_geometry("geometry")

  res.grid$grid <- as.integer(1:nrow(res.grid))
  res.grid <- st_set_agr(res.grid, "constant")
  res.grid$segment <- as.integer(st_within(st_centroid(res.grid), seg.coarse))
  res.grid$dist.lim <- st_distance(res.grid, st_cast(poly.lim, "LINESTRING"))
  res.grid$eligible <- res.grid$dist.lim >= min.bound.dist

  seg.fine <-
    aggregate(res.grid[res.grid$eligible,],
              by = list(seg = res.grid$segment[res.grid$eligible]), identity) |>
    st_geometry() |>
    st_as_sf() |>
    st_set_geometry("geometry")
  seg.fine$segment <- as.integer(1:nrow(seg.fine))

  seg.score <- 
    st_rasterize(res.grid, template = generate_empty(x.dim = x.dim, y.dim = y.dim)) |>
    as.data.table() |>
    merge(as.data.table(score), all.x = FALSE, all.y = TRUE)
  if(score.sam < nrow(seg.score)) {
    seg.score <- seg.score[sample(1:nrow(seg.score), score.sam)][order(x, y)]
  }
  setindex(seg.score, grid)
  setindex(seg.score, segment)
  score.breaks <- get_breaks(seg.score$score, breaks = "FD")
  seg.score[, score.bin := assign_breaks(score, breaks = score.breaks)]

  if(verbose) message("Building areas (coarse selection) …")

  infl <- prod(x.dim, y.dim) / sum(st_area(seg.fine))
  n.tot <- round(area.prop * nrow(seg.fine) * infl)
  n.score <- round(score.prop * n.tot)
  n.rand <- n.tot - n.score

  seg.score.bin <-
    merge(seg.score[, .(n.bin = .N), by = c("segment", "score.bin")],
          as.data.table(expand.grid(segment = unique(seg.score$segment),
                        score.bin = 1:length(score.breaks))),
          by = c("segment", "score.bin"), all = TRUE)
  seg.score.bin[is.na(n.bin), n.bin := 0]
  setindex(seg.score.bin, segment)
  setindex(seg.score.bin, score.bin)

  idx.trt.score <- integer(0)
  if(n.score > 0) {
    for(i in 1:n.score) {
      seg.ids <- with(seg.fine, segment[!(segment %in% idx.trt.score)])
      dist <- list()
      for(j in seq_along(seg.ids)) {
        idx <- c(idx.trt.score, seg.ids[j])
        foc.obs <- seg.score[.(idx), na.omit(score), on = "segment"]
        rem.obs <- seg.score[!.(idx), na.omit(score), on = "segment"]
        foc.freq <-
          seg.score.bin[.(idx),
                        .(n = sum(n.bin)),
                        by = "score.bin", on = "segment"
                        ][!is.na(score.bin)
                          ][order(score.bin), n/sum(n)]
        rem.freq <-
          seg.score.bin[!.(idx),
                        .(n = sum(n.bin)),
                        by = "score.bin", on = "segment"
                        ][!is.na(score.bin)
                          ][order(score.bin), n/sum(n)]
        if(length(foc.freq) > 0 & length(rem.freq) > 0) {
          dist[[j]] <- data.table(segment = seg.ids[j],
                                  dist = js_div(foc.freq, rem.freq,
                                                type = "raw"),
                                  smd = d_cohen(foc.obs, rem.obs))
        }
      }
      dist <- rbindlist(dist)
      idx.trt.score <-
        c(idx.trt.score,
          dist[smd >= 0][which.min(abs(dist - imbalance)), segment])
      }
  }

  # # Old algorithm
  # idx.trt.score <- integer(0)
  # if(n.score > 0) {
  #   for(i in 1:n.score) {
  #     seg.ids <- with(seg.fine, segment[!segment %in% idx.trt.score])
  #     dist <- list()
  #     for(j in seq_along(seg.ids)) {
  #       idx <- c(idx.trt.score, seg.ids[j])
  #       foc.obs <- seg.score[.(idx), na.omit(score), on = "segment"]
  #       rem.obs <- seg.score[!.(idx), na.omit(score), on = "segment"]
  #       if(length(foc.obs) > 0 & length(rem.obs) > 0) {
  #         dist[[j]] <- data.table(segment = seg.ids[j],
  #                                 dist = js_div(foc.obs, rem.obs,
  #                                               type = "discrete"),
  #                                 smd = d_cohen(foc.obs, rem.obs))
  #       }
  #     }
  #     dist <- rbindlist(dist)
  #     idx.trt.score <-
  #       c(idx.trt.score,
  #         dist[smd >= 0][which.min(abs(dist - imbalance)), segment])
  #     }
  # }

  idx.trt.rand <-
    sample(seg.score$segment[-which(seg.score$segment %in% idx.trt.score)], n.rand)
  
  idx.trt <- c(idx.trt.score, idx.trt.rand)
  idx.ref <- seg.fine$segment[!(seg.fine$segment %in% idx.trt)]

  # seg.trt <-
  #   seg.fine[which(seg.fine$segment %in% idx.trt),] |>
  #   st_union(is_coverage = TRUE) |>
  #   st_as_sf() |>
  #   st_set_geometry("geometry")


  if(verbose) message("Adjusting areas (fine tuning) …")

   A.nb <-
      poly2nb(res.grid, queen = FALSE) |>
      nb2listw(style = "B") |> 
      listw2mat()

  grid.trt <- with(res.grid, grid[eligible & segment %in% idx.trt])
  grid.ref <- setdiff(res.grid$grid, grid.trt)
  trt.prop <- seg.score[.(grid.trt), length(segment) / nrow(seg.score), on = "grid"]
  trt.imb <-
    js_div(x = seg.score[.(grid.trt), na.omit(score), on = "grid"],
           y = seg.score[!.(grid.trt), na.omit(score), on = "grid"],
           type = "discrete")
  if(n.rand > 0) {
    # Maintain deviation from target imbalance if randomness was introduced earlier.
    imbalance <- trt.imb
  }
  opt.null <- log((trt.imb - imbalance)^2 + (trt.prop - area.prop)^2)


  seg.score.bin.f <-
    merge(seg.score[, .(n.bin = .N), by = c("grid", "score.bin")],
          expand.grid(grid = unique(seg.score$grid),
                      score.bin = 1:length(score.breaks)),
          by = c("grid", "score.bin"), all = TRUE)
  seg.score.bin.f[is.na(n.bin), n.bin := 0]
  setindex(seg.score.bin.f, grid)
  setindex(seg.score.bin.f, score.bin)
    

  for(i in 1:iter.max) {
    if(trt.prop >= area.prop) {
      grid.adj <- grid.trt[rowSums(A.nb[grid.trt, -grid.trt]) > 0]
      # if(shape.constraints) {
      #     ex.bridge <- exclude_bridges(A.nb, grid.ref, grid.trt)
      #     grid.adj <- grid.adj[!(grid.adj %in% ex.bridge)]
      # }
      adj <- list()
      for(j in seq_along(grid.adj)) {
        idx.trt <- grid.trt[!(grid.trt %in% grid.adj[j])]
        idx.ref <- with(res.grid, grid[!(grid %in% idx.trt)])
        if(shape.constraints) {
          ex.islands.ref <- exclude_islands(A.nb, idx.ref, shape.compactness)
          ex.islands.trt <- exclude_islands(A.nb, idx.trt, shape.compactness)
          idx.trt <- c(idx.trt[!(idx.trt %in% ex.islands.trt)], ex.islands.ref)
        }
        idx.trt.l[[j]] <- idx.trt
        foc.obs <- seg.score[.(idx.trt), na.omit(score), on = "grid"]
        rem.obs <- seg.score[!.(idx.trt), na.omit(score), on = "grid"]
        foc.freq <-
          seg.score.bin.f[.(idx.trt),
                          .(n = sum(n.bin)),
                          by = "score.bin", on = "grid"
                          ][!is.na(score.bin)
                            ][order(score.bin), n/sum(n)]
        rem.freq <-
          seg.score.bin.f[!.(idx.trt),
                          .(n = sum(n.bin)),
                          by = "score.bin", on = "grid"
                          ][!is.na(score.bin)
                            ][order(score.bin), n/sum(n)]
        if(length(foc.obs) > 0 & length(rem.obs) > 0) {
          adj[[j]] <- data.table(grid = grid.adj[j],
                                 dist = js_div(foc.freq, rem.freq,
                                               type = "raw"),
                                 prop = seg.score[.(idx.trt),
                                                  length(grid) / nrow(seg.score),
                                                  on = "grid"],
                                 smd = d_cohen(foc.obs, rem.obs))
        }
      }
    }
    if(trt.prop < area.prop) {
      grid.adj <-
        grid.ref[rowSums(A.nb[grid.ref, -grid.ref]) > 0 &
                 grid.ref %in% with(res.grid, grid[eligible])]
      if(shape.constraints) {
          ex.bridge <- exclude_bridges(A.nb, grid.trt, grid.ref)
          grid.adj <- grid.adj[!(grid.adj %in% ex.bridge)]
      }
      adj <- list()
      idx.trt.l <- list()
      for(j in seq_along(grid.adj)) {
        idx.trt <- c(grid.trt, grid.adj[j])
        idx.ref <- with(res.grid, grid[!(grid %in% idx.trt)])
        if(shape.constraints) {
          ex.islands.ref <- exclude_islands(A.nb, idx.ref, shape.compactness)
          ex.islands.trt <- exclude_islands(A.nb, idx.trt, shape.compactness)
          idx.trt <- c(idx.trt[!(idx.trt %in% ex.islands.trt)], ex.islands.ref)
        }
        idx.trt.l[[j]] <- idx.trt
        foc.obs <- seg.score[.(idx.trt), na.omit(score), on = "grid"]
        rem.obs <- seg.score[!.(idx.trt), na.omit(score), on = "grid"]
        foc.freq <-
          seg.score.bin.f[.(idx.trt),
                          .(n = sum(n.bin)),
                          by = "score.bin", on = "grid"
                          ][!is.na(score.bin)
                            ][order(score.bin), n/sum(n)]
        rem.freq <-
          seg.score.bin.f[!.(idx.trt),
                          .(n = sum(n.bin)),
                          by = "score.bin", on = "grid"
                          ][!is.na(score.bin)
                            ][order(score.bin), n/sum(n)]
        if(length(foc.obs) > 0 & length(rem.obs) > 0) {
          adj[[j]] <- data.table(grid = grid.adj[j],
                                 dist = js_div(foc.freq, rem.freq,
                                               type = "raw"),
                                 prop = seg.score[.(idx.trt),
                                                  length(grid) / nrow(seg.score),
                                                  on = "grid"],
                                 smd = d_cohen(foc.obs, rem.obs))
        }
      }
    } # fi area <
    adj.dt <- rbindlist(adj)
    adj.dt[, `:=`(opt = log((dist - imbalance)^2 + (prop - area.prop)^2))]
    grid.change.j <- adj.dt[smd > 0 & opt < opt.null, which.min(opt)]
    if(length(grid.change.j) > 0) {
      grid.trt <- idx.trt.l[[grid.change.j]]
      grid.ref <- setdiff(res.grid$grid, grid.trt)
      trt.prop <- adj.dt[grid.change.j, prop] 
      trt.imb <- adj.dt[grid.change.j, dist] 
      grid.opt <- adj.dt[grid.change.j, opt]
      if(abs(opt.null - grid.opt) > opt.tol) {
        opt.null <- grid.opt
      } else {
        break
      }
    } else {
      break
    }
    if(verbose & i == iter.max) {
      message(paste0("Maximum number of iterations reached (", iter.max,
                     "). Further optimization is possible."))
    }
  }

# }

  poly.trt <- 
    res.grid[res.grid$grid %in% grid.trt,] |>
    st_union() |>
    st_as_sf() |>
    st_set_geometry("geometry")


  if(trt.prop != area.prop & area.exact == TRUE) {
    if(verbose) message("Resizing area … ")
    opt.bounds <- c(-2, 2) * (abs(sum(trt.prop) - area.prop) * sqrt(prod(x.dim, y.dim)))
    opt.res <-
      optimize(opt_resize_poly, interval = opt.bounds,
               poly = poly.trt, limit = poly.lim, prop.target = area.prop)
    poly.trt <-
      resize_polygon(poly.trt, limit = poly.lim,
                       dist = opt.res$minimum ) |>
      st_union() |>
      st_as_sf() |>
      st_set_geometry("geometry")
  }

  poly.ref <-
    st_difference(poly.lim, poly.trt) |>
    st_as_sf() |>
    st_set_geometry("geometry")

  poly.areas <- st_as_sf(rbind(poly.ref, poly.trt))
  poly.areas$type <- factor(c("reference", "treatment"),
                            levels = c("reference", "treatment"))

  if(verbose) {
    areas.dt <-
      st_rasterize(poly.areas, template = score) |>
      as.data.table() |>
      merge(score, all.x = FALSE, all.y = TRUE)
    trt.prop <- areas.dt[type == "treatment", length(type) / nrow(areas.dt)]
    trt.imb <-
      js_div(x = areas.dt[type == "treatment", na.omit(score)],
             y = areas.dt[type == "reference", na.omit(score)],
             type = "discrete")
    message(paste0("Imbalance: ", round(trt.imb, getOption("digits")),
                   "\nArea proportion: ", round(trt.prop, getOption("digits"))))
  }
  return(poly.areas)
}


match_to_segments_old <- function(segments, pw.dist, center, buffer, collapse = TRUE) {
  stopifnot(length(center) == length(buffer))
  center.inb <- center %in% 1:nrow(pw.dist)
  matched.l <- list()
  for(i in seq_along(center)) {
    if(center.inb[i]) {
      seg.sel <- which(segments == i)
      grid.sel <- unname(which(pw.dist[, center[i]] <= buffer[i]))
      matched.l[[i]] <- grid.sel[grid.sel %in% seg.sel]
    } else {
      matched.l[[i]] <- integer(0)
    }
  }
  if(collapse == TRUE) matched.l <- do.call(c, matched.l)
  return(matched.l)
}


match_to_segments <- function(pw.dist,
                              center,
                              buffer,
                              min.dist = 0,
                              include.unmatched = TRUE) {
    stopifnot(length(center) == length(buffer))
    center.inb <- center %in% 1:nrow(pw.dist)
    matched.l <- list()
    if(any(center.inb)) {
      grid.seg <-
        data.table(grid = 1:nrow(pw.dist),
                   segment =
                     apply(pw.dist[, center[center.inb], drop = FALSE],
                           1, which.min) |>
                     lapply(\(x) ifelse(length(x) > 0, x, 0)) |>
                     unlist())
      grid.seg[, dist := pw.dist[grid, center[center.inb][segment]], by = .I]
      for(i in seq_along(center)) {
        if(center.inb[i]) {
          matched.l[[i]] <- grid.seg[segment == i & dist <= buffer[i], grid]
        } else {
          matched.l[[i]] <- integer(0)
        }
      }
      matched.grid <- do.call(c, matched.l)
      matched.segment <- do.call(c, mapply(\(x, y) rep(y, length(x)),
                                           matched.l, 1:sum(center.inb)))
    } else {
      matched.grid <- integer(0)
      matched.segment <- integer(0)
    }
    if(include.unmatched) {
      unmatched.grid <-  (1:nrow(pw.dist))[!1:nrow(pw.dist) %in% matched.grid]
      unmatched.segment <- rep(0, length(unmatched.grid))
      matched.grid <- c(matched.grid, unmatched.grid)
      matched.segment <- c(matched.segment, unmatched.segment)
    }
    return(list(grid = matched.grid, segment = matched.segment))
  }


segment_min_distance <- function(pw.dist,
                                 grid.cols = c("grid1", "grid2"),
                                 grid.sel,
                                 segments) {
  grid.seg <- data.table(grid = grid.sel, segment = segments)
  grid.dist.dt <- as.data.table(pw.dist[grid.sel, grid.sel])
  names(grid.dist.dt) <- as.character(grid.sel)
  grid.dist.dt[, grid1 := grid.sel]
  grid.dist.dt <-
    melt(grid.dist.dt,
         id.vars = "grid1",
         variable.name = "grid2", value.name = "dist",
         variable.factor = FALSE) |>
    DT(, grid2 := as.integer(grid2)) |>
    merge(grid.seg[, .(grid1 = grid, segment1 = segment)],
          by = "grid1") |>
    merge(grid.seg[, .(grid2 = grid, segment2 = segment)],
          by = "grid2") |>
  DT(segment1 != segment2,
     .(dist.min = min(dist)),
     by = c("segment1", "segment2")) |>
  DT(order(segment1, segment2))
  return(grid.dist.dt)
}


encode_grid_buffer <- function(grid, buffer, order.grid, order.buffer) {
  grid.gray <- lapply(grid, \(x) binary2gray(decimal2binary(x, order.grid)))
  buffer.gray <- lapply(buffer, \(x) binary2gray(decimal2binary(x, order.buffer)))
  string <- c(do.call(c, grid.gray), do.call(c, buffer.gray))
  return(string)
}

decode_grid_buffer <- function(string, bit.orders) {
  string <- split(string, rep.int(seq.int(bit.orders), times = bit.orders))
  orders <- sapply(string, function(x) { binary2decimal(gray2binary(x)) })
  dec <- unname(orders)
  grid.i <- 1:round(0.5*length(dec))
  return(list(grid = dec[grid.i],
              buffer = dec[-grid.i]))
}

exclude_islands <- function(x, a, min.vertices = 1) {
  x.a <- x[a, a, drop = FALSE] 
  diag(x.a) <- 0
  to_exclude <- a[rowSums(x.a) < min.vertices]
  return(to_exclude)
}

exclude_bridges <- function(x, a, b) {
  filter.a <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  filter.a[a, a] <- 1
  diag(filter.a) <- 0
  filter.b <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  filter.b[b, b] <- 1
  diag(filter.b) <- 0
  filter.ab <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  filter.ab[a, b] <- 1
  filter.ab[b, a] <- 1
  diag(filter.ab) <- 0
  # Which nodes of type `a` are seperated by more than 2 other nodes of type `a`
  x.a <- x * filter.a
  x.a2 <- x.a %*% x.a
  x.a3 <- x.a2 %*% x.a
  x.asep2 <- x.a + x.a2 + x.a3
  x.asep2 <- (x.asep2 == 0) * filter.a
  # Identify type `a` nodes that are connected to a type `b` nodes; either
  # directly, or through one other `a` node.
  x.ab <- x * filter.ab
  x.b <- x * filter.b
  x.b2a <- x.ab + (x.ab %*% x.b)
  # Which nodes of type `b` are on a path that connects two type `a` nodes
  # that are seperated as defined above ?
  x.bridge <- x.b2a %*% x.asep2 %*% x.b2a
  to_exclude <- b[diag(x.bridge)[b] > 0]
  return(to_exclude)
}


generate_areas_poly2 <- function(x.dim,
                                y.dim,
                                imbalance = 0.1,
                                imbalance.tol = NULL,
                                area.prop = 0.5,
                                area.tol = NULL,
                                area.exact = FALSE,
                                score,
                                score.sam = 1e4,
                                seg.seed,
                                seg.res = 0.1 * min(x.dim, y.dim),
                                seg.min.dist = 0,
                                seg.min.area = 0,
                                seg.even = 2,
                                seg.prec = 0.1 * seg.res,
                                min.bound.dist = 0,
                                verbose = FALSE,
                                opt.imp.imb = 1,
                                opt.imp.area = 1,
                                opt.pop = 50,
                                opt.prec = 1e-4,
                                opt.pcrossover = 0.8,
                                opt.pmutation = 0.1,
                                opt.max.iter = 100,
                                opt.run = opt.max.iter,
                                opt.parallel = FALSE,
                                opt.fine = TRUE,
                                opt.fine.max.iter = 25,
                                opt.fine.constr = TRUE,
                                opt.fine.tol = 1e-4,
                                ...
                                ) {

# for(m in 1:20){
  if(verbose) message("Setting up candidate polygons ...")

  if(is.null(area.tol)) {
    area.lim <- rep(area.prop, 2)
  } else{
    if(length(area.tol) == 1) {
      area.lim <- area.prop + c(-area.tol, area.tol)
    } else {
      area.lim <- area.prop + area.tol
    }
  }

  if(is.null(imbalance.tol)) {
    imbalance.lim <- rep(imbalance, 2)
  } else{
    if(length(imbalance.tol) == 1) {
      imbalance.lim <- imbalance + c(-imbalance.tol, imbalance.tol)
    } else {
      imbalance.lim <- imbalance + imbalance.tol
    }
  }

  # if(is.null(score)) {
  #   score <- generate_empty(x.dim = x.dim, y.dim = y.dim) + 1
  #   score.prop <- 0
  # } else {
  # }
  score <- score[1]
  names(score)[1] <- "score"
  score <- scale_int(score, int = c(0, 1))

  poly.lim <- limit_polygon(x.dim, y.dim)
  
  # if(seg.seed > 1) {
  #   seg.sp <- max(seg.min.dist * 2, seg.res)
  #   x.crd.seg <- seq(0, x.dim - (0.5 * seg.sp), seg.sp) + (0.5 * seg.sp)
  #   y.crd.seg <- seq(0, y.dim - (0.5 * seg.sp), seg.sp) + (0.5 * seg.sp)
  #   seg.grid <-
  #     cbind(x = rep(x.crd.seg, times = length(y.crd.seg)),
  #           y = rep(y.crd.seg, each = length(y.crd.seg))) |>
  #     st_multipoint() |>
  #     st_sfc() |>
  #     st_cast("POINT") |>
  #     st_as_sf()
  #   seg.sam <- sample(1:nrow(seg.grid), seg.seed)
  #   seg.coarse <- 
  #       seg.grid[seg.sam,] |>
  #       st_union() |>
  #       st_voronoi() |>
  #       st_cast() |>
  #       st_intersection(poly.lim) |>
  #       st_as_sf() |>
  #       st_set_geometry("geometry")
  # } else {
  #   seg.coarse <-
  #     st_as_sf(poly.lim) |>
  #     st_set_geometry("geometry")
  # }

  # set up segment grid

  res.grid.coord <- 
    expand.grid(x = seq(1, x.dim - (0.5 * seg.res), seg.res) + (0.5 * seg.res),
                y = seq(1, y.dim - (0.5 * seg.res), seg.res) + (0.5 * seg.res)) |>
    as.data.table()
  res.grid.coord


  # seg.even <- 0.75
  x.crd.res <- seq(- seg.res, x.dim + seg.res, seg.res) + (0.5 * seg.res)
  y.crd.res <- seq(- seg.res, y.dim + seg.res, seg.res) + (0.5 * seg.res)
  # x.crd.res <- seq(-0.5 * seg.res, x.dim + (0.5 * seg.res), seg.res) + (0.5 * seg.res)
  # y.crd.res <- seq(-0.5 * seg.res, y.dim + (0.5 * seg.res), seg.res) + (0.5 * seg.res)
  grid.per <- 
    matrix((rbeta(n = 2 * length(x.crd.res) * length(y.crd.res),
                  seg.even, seg.even) * seg.res) -
           (0.5 * seg.res),
           ncol = 2)
  grid.per <- 
    round(grid.per / seg.prec) * seg.prec
  res.grid <-
    (cbind(x = rep(x.crd.res, times = length(y.crd.res)),
           y = rep(y.crd.res, each = length(y.crd.res))) +
     grid.per) |>
    st_multipoint() |>
    st_sfc() |>
    st_cast("POINT") |>
    st_as_sf() |>
    st_union() |>
    st_voronoi() |>
    # st_triangulate() |>
    st_cast() |>
    st_intersection(poly.lim) |>
    st_as_sf() |>
    # st_simplify(preserveTopology = TRUE, dTolerance = seg.prec) |> 
    st_snap_to_grid(seg.prec) |>
    st_make_valid() |>
    st_set_geometry("geometry")
  st_agr(res.grid) <- "constant"
  res.grid <- res.grid[!st_is_empty(res.grid) & st_geometry_type(res.grid) == "POLYGON",]
  # plot(res.grid)

  res.grid$area <- st_area(res.grid)
  res.grid$dist.lim <- st_distance(res.grid, st_cast(poly.lim, "LINESTRING"))
  res.grid <- res.grid[res.grid$dist.lim >= min.bound.dist,]
  res.grid$grid <- as.integer(1:nrow(res.grid))
  res.grid$selected <- FALSE
  res.grid <- st_set_agr(res.grid, "constant")

  res.grid.dt <- as.data.table(st_drop_geometry(res.grid))
  setindex(res.grid.dt, grid)
  setindex(res.grid.dt, selected)
  setorder(res.grid.dt, grid)

  grid.dist.pw <- suppressWarnings(st_distance(res.grid))


  seg.score <- 
    st_rasterize(res.grid[, "grid"],
                 template = generate_empty(x.dim = x.dim, y.dim = y.dim)) |>
    as.data.table() |>
    merge(as.data.table(score), all.x = FALSE, all.y = TRUE)
  if(score.sam < nrow(seg.score)) {
    seg.score <- seg.score[sample(1:nrow(seg.score), score.sam)][order(x, y)]
  }
  setindex(seg.score, grid)
  score.breaks <- get_breaks(seg.score$score, breaks = "FD")
  seg.score[, score.bin := assign_breaks(score, breaks = score.breaks)]

  seg.score.bin <-
    merge(seg.score[, .(n.bin = .N), by = c("grid", "score.bin")],
          expand.grid(grid = unique(seg.score$grid),
                      score.bin = 1:length(score.breaks)),
          by = c("grid", "score.bin"), all = TRUE)
  seg.score.bin[is.na(n.bin), n.bin := 0]
  setindex(seg.score.bin, grid)
  setindex(seg.score.bin, score.bin)


  if(verbose) message("Building areas (optimisation) …")

  order.grid <- ceiling(log2(nrow(res.grid.dt)))
  order.buffer <- ceiling(log2(sqrt(x.dim^2 + y.dim^2)))
  order.c <- rep(c(order.grid, order.buffer), each = seg.seed)
  n.bits <- sum(order.c)

  f_opt_imb_area <- function(x) {
    x <- decode_grid_buffer(x, order.c)
    n.seg <- length(x$grid)
    trt.proposed <-
      match_to_segments.m(pw.dist = grid.dist.pw,
                          center = x$grid,
                          buffer = x$buffer,
                          include.unmatched = FALSE) |>
      as.data.table()
    if(length(trt.proposed$grid) == 0) {
      f <- length(x$grid)^2 * sqrt(.Machine$double.xmax)
    } else {
      foc.obs <- seg.score[.(trt.proposed$grid), na.omit(score), on = "grid"]
      rem.obs <- seg.score[!.(trt.proposed$grid), na.omit(score), on = "grid"]
      foc.freq <-
        seg.score.bin[.(trt.proposed$grid),
                      .(n = sum(n.bin)),
                      by = "score.bin", on = "grid"
                      ][!is.na(score.bin)
                        ][order(score.bin), n/sum(n)]
      rem.freq <-
        seg.score.bin[!.(trt.proposed$grid),
                      .(n = sum(n.bin)),
                      by = "score.bin", on = "grid"
                      ][!is.na(score.bin)
                       ][order(score.bin), n/sum(n)]
      trt.imb <- js_div.m(foc.freq, rem.freq, type = "raw")
      trt.a <- 
        merge(res.grid.dt[.(trt.proposed$grid),
                          on = "grid"],
              trt.proposed) |>
        DT(, .(area = sum(area)), by = "segment")
      # Any remaining segments without assigned grid polygons?
      seg.rem <- setdiff(1:n.seg, trt.a$segment)
      if(length(seg.rem) > 0) {
        trt.a <-
          rbind(trt.a,
                data.table(segment = seg.rem,
                           area = rep(0, length(seg.rem))))
      }
      trt.area <- sum(trt.a) / prod(x.dim, y.dim)
      imbalance.loss <- (trt.imb - imbalance)^2
      area.loss <- (trt.area - area.prop)^2
      penalty.min.area <- pmax(-trt.a$area + seg.min.area, 0) * sqrt(.Machine$double.xmax)
      penalty.min.dist <-
        segment_min_distance.m(pw.dist = grid.dist.pw,
                               grid.cols = c("grid1", "grid2"),
                               grid.sel = trt.proposed$grid,
                               segments = trt.proposed$segment) |>
        DT(, as.numeric(dist.min < seg.min.dist) * sqrt(.Machine$double.xmax))
      penalty.oob <- as.numeric(!(x$grid %in% res.grid.dt$grid)) * sqrt(.Machine$double.xmax)
      f <-
        opt.imp.imb * imbalance.loss +
        opt.imp.area * area.loss +
        sum(penalty.min.area) + sum(penalty.oob) + sum(penalty.min.dist)
    }
    return(-f)
  }

  # Cache function evaluations
  fm_opt_imb_area <- memoise(f_opt_imb_area)
  match_to_segments.m <- memoise(match_to_segments)
  segment_min_distance.m <- memoise(segment_min_distance)
  js_div.m <- memoise(js_div)
  
  pop.start.mat <- matrix(NA, nrow = opt.pop, ncol = sum(order.c))
  for(i in 1:opt.pop) {
    grid <- sample(res.grid.dt$grid, seg.seed)
    buffer <- round(runif(n = 3,
                          sqrt(seg.min.area),
                          sqrt(area.prop * x.dim * y.dim / seg.seed)))
    pop.start.mat[i,] <- encode_grid_buffer(grid, buffer, order.grid, order.buffer)
  }

  f_max <- -((opt.imp.imb + opt.imp.area) * opt.prec)

  f_monitor <- ifelse(verbose, gaMonitor, FALSE)

  buffer.opt <- NULL
  while(is.null(buffer.opt)) {
    try({
      buffer.opt <-
          ga(type =  "binary",
             fitness = fm_opt_imb_area,
             maxFitness = f_max,
             # lower = rep(0, seg.seed), upper = rep(max(x.dim, y.dim), seg.seed),
             nBits = sum(order.c), suggestions = pop.start.mat,
             run = opt.run, popSize = opt.pop, maxiter = opt.max.iter,
             optim = TRUE, optimArgs = list(poptim = 0.25),
             parallel = opt.parallel,
             pcrossover = opt.pcrossover, pmutation = opt.pmutation,
             monitor = f_monitor)
    })
    if(verbose & is.null(buffer.opt)) message("Optimization failed. Trying again …")
  }

  imb_area.sol <- decode_grid_buffer(buffer.opt@solution[1,], order.c)

  trt.proposed <-
    match_to_segments.m(pw.dist = grid.dist.pw,
                      center = imb_area.sol$grid,
                      buffer = imb_area.sol$buffer,
                      include.unmatched = TRUE) |>
  as.data.table()

  # Update segment and selection information for grid polygons

  res.grid <-
    as.data.table(res.grid) |>
    merge(trt.proposed) |>
    DT(segment != 0, selected := TRUE) |>
    st_as_sf()

  res.grid.dt <-
    merge(res.grid.dt,
          trt.proposed) |>
    DT(segment != 0, selected := TRUE)
  setindex(res.grid.dt, grid)
  setindex(res.grid.dt, segment)
  setindex(res.grid.dt, selected)
  setorder(res.grid.dt, grid)

  grid.dist.dt <- as.data.table(grid.dist.pw)
  names(grid.dist.dt) <- as.character(1:nrow(res.grid))
  grid.dist.dt[, grid1 := 1:.N]
  grid.dist.dt <-
    melt(grid.dist.dt, id.vars = "grid1",
         variable.name = "grid2", value.name = "dist",
         variable.factor = FALSE)
  grid.dist.dt[,
              c("grid1", "grid2") := lapply(.SD, \(x) as.integer(as.character(x))),
              .SDcols = c("grid1", "grid2")]
  grid.dist.dt <-
    merge(grid.dist.dt,
          res.grid.dt[, .(grid1 = grid, segment1 = segment)], by = "grid1") |>
    merge(res.grid.dt[, .(grid2 = grid, segment2 = segment)], by = "grid2")

  # plot(res.grid[, "selected"])

  # Extract results of coarse optimisation

  grid.trt <- with(res.grid, grid[selected])
  grid.ref <- setdiff(res.grid$grid, grid.trt)
  trt.a <- 
    res.grid.dt[.(grid.trt),
                .(area = sum(area)),
                by = "segment",
                on = "grid"
                ][, area]

  trt.area <- seg.score[.(grid.trt), length(grid) / nrow(seg.score), on = "grid"]
  trt.imb <-
    js_div(x = seg.score[.(grid.trt), na.omit(score), on = "grid"],
           y = seg.score[!.(grid.trt), na.omit(score), on = "grid"],
           type = "discrete")


  # remove cached functions
  rm(fm_opt_imb_area, match_to_segments.m, segment_min_distance.m, js_div.m)

  if(verbose) {
    message(paste0("Imbalance: ", round(trt.imb, getOption("digits")),
                   ". Area proportion: ", round(trt.area, getOption("digits"))))
  }

  if(opt.fine) {

    if(verbose) message("Attempt to further adjust areas (fine tuning) …")

    penalty <- sqrt(opt.imp.imb + opt.imp.area)
    pen.imb = ifelse(trt.imb < area.lim[1] | trt.imb > area.lim[2],
                     min(abs(trt.imb - area.lim)) * sqrt(penalty), 0)
    pen.area = ifelse(trt.area < area.lim[1] | trt.area > area.lim[2],
                     min(abs(trt.area - area.lim)) * sqrt(penalty), 0)
    pen.min.area <- sum(pmax(-trt.a + seg.min.area, 0) * sqrt(penalty))
    pen.min.dist <-
      grid.dist.dt[grid1 %in% grid.trt & grid2 %in% grid.trt,
                   .(dist.min = min(dist)),
                   by = c("segment1", "segment2")
                   ][segment1 != segment2,
                     sum(0.5 * as.numeric(dist.min < seg.min.dist) * sqrt(penalty))]
    imb.loss <- (trt.imb - imbalance)^2
    area.loss <- (trt.area - area.prop)^2
    loss.null <- opt.imp.imb * imb.loss +
                 opt.imp.area * area.loss +
                 pen.min.area + pen.min.dist + pen.imb + pen.area

    A.nb <- st_relate(res.grid, res.grid, pattern = "F***1****", sparse = FALSE) * 1
    diag(A.nb) <- 0
    A.nb <- as(A.nb, "sparseMatrix")

#     A.nb <-
#         poly2nb(res.grid, queen = FALSE) |>
#         nb2listw(style = "B") |> 
#         listw2mat()
    loss.hist <- list()

    for(i in 1:opt.fine.max.iter) {
      # if(trt.area <= area.prop) {
        grid.adj.add <-
          grid.ref[rowSums(A.nb[grid.ref, -grid.ref]) > 0 &
                   grid.ref %in% with(res.grid, grid[segment != 0])]
        # if(opt.fine.constr) {
        #     ex.bridge <- exclude_bridges(A.nb, grid.trt, grid.ref)
        #     grid.adj.add <- grid.adj.add[!(grid.adj.add %in% ex.bridge)]
        # }
        adj.add <- list()
        idx.trt.add.l <- list()
        for(j in seq_along(grid.adj.add)) {
          idx.trt <- c(grid.trt, grid.adj.add[j])
          idx.ref <- with(res.grid, grid[!(grid %in% idx.trt)])
          if(opt.fine.constr) {
            ex.islands.ref <- exclude_islands(A.nb, idx.ref, 2)
            ex.islands.trt <- exclude_islands(A.nb, idx.trt, 2)
            idx.trt <- c(idx.trt[!(idx.trt %in% ex.islands.trt)], ex.islands.ref)
            idx.trt <- idx.trt[idx.trt %in% with(res.grid, grid[segment != 0])]
          }
          idx.trt.add.l[[j]] <- idx.trt
          foc.obs <- seg.score[.(idx.trt), na.omit(score), on = "grid"]
          rem.obs <- seg.score[!.(idx.trt), na.omit(score), on = "grid"]
          foc.freq <-
            seg.score.bin[.(idx.trt),
                          .(n = sum(n.bin)),
                          by = "score.bin", on = "grid"
                          ][!is.na(score.bin)
                            ][order(score.bin), n/sum(n)]
          rem.freq <-
            seg.score.bin[!.(idx.trt),
                          .(n = sum(n.bin)),
                          by = "score.bin", on = "grid"
                          ][!is.na(score.bin)
                            ][order(score.bin), n/sum(n)]
         foc.a <- 
            res.grid.dt[.(idx.trt),
                        .(area = sum(area)),
                        by = "segment",
                        on = "grid"
                        ][, area]
        pen.min.dist <-
          grid.dist.dt[grid1 %in% idx.trt & grid2 %in% idx.trt,
                       .(dist.min = min(dist)),
                       by = c("segment1", "segment2")
                       ][segment1 != segment2,
                         0.5 * as.numeric(dist.min < seg.min.dist) * sqrt(penalty)]
          if(length(foc.obs) > 0 & length(rem.obs) > 0) {
            adj.add[[j]] <- data.table(grid = grid.adj.add[j],
                                       # smd = d_cohen(foc.obs, rem.obs),
                                       imb = js_div(foc.freq, rem.freq,
                                                    type = "raw"),
                                       area = sum(foc.a) / prod(x.dim, y.dim),
                                       pen.min.area = sum(pmax(-foc.a + seg.min.area, 0) *
                                                          sqrt(penalty)),
                                       pen.min.dist = sum(pen.min.dist))
          }
        }
      # } # fi area <
      # if(trt.area > area.prop) {
        grid.adj.rem <- grid.trt[rowSums(A.nb[grid.trt, -grid.trt]) > 0]
        # if(opt.fine.constr) {
        #     ex.bridge <- exclude_bridges(A.nb, grid.ref, grid.trt)
        #     grid.adj.rem <- grid.adj.rem[!(grid.adj.rem %in% ex.bridge)]
        # }
        adj.rem <- list()
        idx.trt.rem.l <- list()
        for(j in seq_along(grid.adj.rem)) {
          idx.trt <- grid.trt[!(grid.trt %in% grid.adj.rem[j])]
          idx.ref <- with(res.grid, grid[!(grid %in% idx.trt)])
          if(opt.fine.constr) {
            ex.islands.ref <- exclude_islands(A.nb, idx.ref, 2)
            ex.islands.trt <- exclude_islands(A.nb, idx.trt, 2)
            idx.trt <- c(idx.trt[!(idx.trt %in% ex.islands.trt)], ex.islands.ref)
            idx.trt <- idx.trt[idx.trt %in% with(res.grid, grid[segment != 0])]
          }
          idx.trt.rem.l[[j]] <- idx.trt
          foc.obs <- seg.score[.(idx.trt), na.omit(score), on = "grid"]
          rem.obs <- seg.score[!.(idx.trt), na.omit(score), on = "grid"]
          foc.freq <-
            seg.score.bin[.(idx.trt),
                          .(n = sum(n.bin)),
                          by = "score.bin", on = "grid"
                          ][!is.na(score.bin)
                            ][order(score.bin), n/sum(n)]
          rem.freq <-
            seg.score.bin[!.(idx.trt),
                          .(n = sum(n.bin)),
                          by = "score.bin", on = "grid"
                          ][!is.na(score.bin)
                            ][order(score.bin), n/sum(n)]
          foc.a <- 
            res.grid.dt[.(idx.trt),
                        .(area = sum(area)),
                        by = "segment",
                        on = "grid"
                        ][, area]
          if(length(foc.obs) > 0 & length(rem.obs) > 0) {
            adj.rem[[j]] <- data.table(grid = grid.adj.rem[j],
                                       # smd = d_cohen(foc.obs, rem.obs),
                                       imb = js_div(foc.freq, rem.freq,
                                                    type = "raw"),
                                       area = sum(foc.a) / prod(x.dim, y.dim),
                                       pen.min.area = sum(pmax(-foc.a + seg.min.area, 0) *
                                                          sqrt(penalty)),
                                       pen.min.dist = sum(pen.min.dist))
          }
        }
      # }
      idx.trt.l <- c(idx.trt.add.l, idx.trt.rem.l)
      adj.dt <- rbind(rbindlist(adj.add), rbindlist(adj.rem))
      adj.dt[, `:=`(pen.imb = 0, pen.area = 0)]
      adj.dt[imb < imbalance.lim[1] | imb > imbalance.lim[2],
             `:=`(
                  pen.imb = fifelse(imb < imbalance.lim[1],
                                    abs(imb - imbalance.lim[1]) * sqrt(penalty),
                                    abs(imb - imbalance.lim[2]) * sqrt(penalty)
                                    )
                  )]
      adj.dt[area < area.lim[1] | area > area.lim[2],
             `:=`(
                  pen.area = fifelse(area < area.lim[1],
                                     abs(area - area.lim[1]) * sqrt(penalty),
                                     abs(area - area.lim[2]) * sqrt(penalty)
                                     )
                  )]
      adj.dt[, `:=`(imbalance.loss = opt.imp.imb * (imb - imbalance)^2,
                    area.loss = opt.imp.area * (area - area.prop)^2
                 )]
      adj.dt[, `:=`(loss = opt.imp.imb * imbalance.loss +
                           opt.imp.area * area.loss +
                           pen.min.area + pen.imb + pen.area)]
      grid.change.j <- adj.dt[loss < loss.null, which.min(loss)]
      if(length(grid.change.j) > 0) {
        grid.trt <- idx.trt.l[[grid.change.j]]
        grid.ref <- setdiff(res.grid$grid, grid.trt)
        trt.area <- adj.dt[grid.change.j, area] 
        trt.imb <- adj.dt[grid.change.j, imb] 
        grid.loss <- adj.dt[grid.change.j, loss]
        loss.hist[[i]] <- grid.loss
        if((loss.null - grid.loss) > opt.fine.tol) {
          loss.null <- grid.loss
          if(verbose) {
            message(paste0("Iteration ", i,
                           ". Imbalance: ", round(trt.imb, getOption("digits")),
                           ". Area proportion: ", round(trt.area, getOption("digits")),
                           ". Loss: ", round(grid.loss, getOption("digits"))))
          }
        } else {
          message(paste0("Iteration ", i, ". No further improvement."))
          break
        }
      } else {
        message(paste0("Iteration ", i, ". No further improvement."))
        break
      }
      if(verbose & i == opt.fine.max.iter) {
        message(paste0("Maximum number of iterations reached (", opt.fine.max.iter,
                       "). Further optimization may be possible."))
      }
    }

  } # End fine-scale optimisation

  res.grid$selected[-grid.trt] <- FALSE
  res.grid$selected[grid.trt] <- TRUE

  res.grid.sel <- res.grid[res.grid$selected,]

  poly.trt <-
    with(res.grid.sel, aggregate(geometry, by = list(segment = segment), st_union)) |>
    st_as_sf()

  if(trt.area != area.prop & area.exact == TRUE) {
    if(verbose) message("Resizing area … ")

    res.grid.seg <- res.grid[res.grid$segment != 0,]
    poly.seg.ex <-
      with(res.grid.seg, aggregate(geometry, by = list(segment = segment), st_union)) |>
      st_as_sf()
    st_agr(poly.seg.ex) <- "constant"

    f_opt_area_exact <- function(x) {
      poly.new <- resize_polygons(poly.trt,
                                  limit = poly.seg.ex,
                                  dist = x,
                                  endCapStyle  = "SQUARE",
                                  joinStyle = "MITRE",
                                  mitreLimit = 1)
      prop.new <- sum(st_area(poly.new)) / prod(x.dim, y.dim)
      return(-abs(prop.new - area.prop))
    }

    opt.bounds <- c(-2, 2) * (abs(sum(trt.area) - area.prop) * sqrt(prod(x.dim, y.dim)))

    opt.res <-
      ga(type = "real-valued", fitness = \(x) suppressWarnings(f_opt_area_exact(x)),
         lower = opt.bounds[1], upper = opt.bounds[2],
         popSize = 20, run = 10)

    poly.trt <-
      resize_polygons(poly.trt,
                      limit = poly.seg.ex,
                      dist = opt.res@solution[1,],
                      endCapStyle  = "SQUARE",
                      joinStyle = "MITRE",
                      mitreLimit = 1)
  }

  poly.trt <-
    st_union(poly.trt) |>
    st_as_sf() |>
    st_set_geometry("geometry")

  poly.ref <-
    st_difference(poly.lim, poly.trt) |>
    st_as_sf() |>
    st_set_geometry("geometry")

  poly.areas <- st_as_sf(rbind(poly.ref, poly.trt))
  poly.areas$type <- factor(c("reference", "treatment"),
                            levels = c("reference", "treatment"))

  areas.dt <-
    st_rasterize(poly.areas, template = score) |>
    as.data.table() |>
    merge(score, all.x = FALSE, all.y = TRUE)
  trt.area <- areas.dt[type == "treatment", length(type) / nrow(areas.dt)]
  trt.imb <-
    js_div(x = areas.dt[type == "treatment", na.omit(score)],
           y = areas.dt[type == "reference", na.omit(score)],
           type = "discrete")

  if(verbose) {
    message(paste0("Final result: Imbalance: ", round(trt.imb, getOption("digits")),
                   ". Area proportion: ", round(trt.area, getOption("digits"))))
  }

  return(list(shape = poly.areas, table = areas.dt[, .(x, y, type)],
              performance = c(imbalance = trt.imb, area.prop = trt.area)))
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
                          prop = split.prop)
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

