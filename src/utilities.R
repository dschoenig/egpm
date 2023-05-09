library(mgcv)
library(RandomFields)
library(RandomFieldsUtils)
RFoptions(install="no")
RFoptions(spConform=FALSE)
# library(raster)
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


cloglog <- function(mu) {
  log(-log(1 - mu))
}

inv_cloglog <- function(eta) {
  pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
}

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
  if(all(int == FALSE)) {
    return(x)
  } else if(!is.null(int)) {
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
  if(all(int == FALSE)) {
    return(x)
  } else if(!is.null(int)) {
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


invert <- function(x, ...) {
  UseMethod("invert", x)
}


invert.stars <- function(x, attributes = NULL) {
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
                           value = 0,
                           name = "value") {
  mat <- matrix(value, nrow = x.dim, ncol = y.dim)
  r <- matrix2stars(mat, name = name)
  return(r)
}

generate_nugget <- function(x.dim,
                            y.dim,
                            var = 1,
                            name = "value") {
  mat <- matrix(rnorm(prod(x.dim, y.dim), 0, var),
                nrow = x.dim,
                ncol = y.dim)
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


cos2_2d_trunc <- function(x, y, a = 1, xd = 1, yd = xd, x0 = 0, y0 = 0) {
  bx <- pi/(2*xd)
  by <- pi/(2*yd)
  x.t <- bx * (x - x0)
  y.t <- by * (y - y0)
  # x.t <- 1/b * x - x0
  # y.t <- 1/b * y - y0
  z <- a * ((cos(x.t))^2 * (cos(y.t))^2)
  z[x.t < -pi/2 | x.t > pi/2] <- 0
  z[y.t < -pi/2 | y.t > pi/2] <- 0
  return(z)
}

generate_sin_gradient <- function(x.dim,
                                  y.dim,
                                  phi = NULL,
                                  shift = NULL,
                                  centre = NULL,
                                  x.scale = x.dim/2,
                                  y.scale = y.dim/2,
                                  rescale = c(0,1),
                                  name = "value") {
  if(is.null(centre)) {
    centre <- point_on_square(phi)
    centre <- centre - (shift * (centre - 0.5))
    centre <- centre * c(x.dim, y.dim)
  }
  bx <- pi/(2*x.scale)
  by <- pi/(2*y.scale)
  x.seq <- bx * (1:x.dim - centre[1])
  y.seq <- by * (1:y.dim - centre[2])
  x.seq[x.seq < -pi/2 | x.seq > pi/2] <- pi/2
  y.seq[y.seq < -pi/2 | y.seq > pi/2] <- pi/2
  grad.mat <-
    outer(x.seq, y.seq,
          # FUN = \(x,y) exp(-((x - x.off)^2 + (y - y.off)^2)))
          FUN = \(x,y) ((cos(x))^2 * (cos(y))^2))
  matrix2stars(grad.mat)
  grad <-
      matrix2stars(grad.mat, name = name) |> 
      scale_int(int = rescale)
  return(grad)
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

generate_matern <-
  function(x.dim,
           y.dim,
           nu = 1,
           var = 1,
           scale = 1,
           mean = 0,
           rescale = c(0,1),
           name = "value") {
  mod <- RMmatern(nu = nu, var = var, scale = scale) + mean
  matern <- 
    RFsimulate(model = mod, x = 1:x.dim, y = 1:y.dim) |>
    as.matrix() |>
    t() |>
    matrix2stars(name = name)
  if(is.numeric(rescale)) {
    matern <- scale_int(matern, int = rescale)
  }
  return(matern)
}



# # Backup if RandomFields eventually breaks
# generate_matern_fields <-
#   function(x.dim,
#            y.dim,
#            nu = 1,
#            scale = max(x.dim, y.dim)/10,
#            var = 1,
#            mean = 0,
#            rescale = c(0,1),
#            name = "value",
#            ...) {
#   ce.setup <-
#     find_ce_domain(x.dim = x.dim,
#                    y.dim = y.dim,
#                    cov.args = list(Covariance= "Matern",
#                                             aRange = scale,
#                                             smoothness = nu),
#                    cov.function = "stationary.cov",
#                    return.setup = TRUE,
#                    ...)
#   if(is.null(ce.setup)) stop("Domain to small for chosen range.")
#   field <-
#     (circulantEmbedding(ce.setup) * var) |>
#     matrix2stars(name = name)
#   field <- field + mean
#   if(is.numeric(rescale)) {
#     field <- scale_int(field, int = rescale)
#   }
#   return(field)
#   }

# find_ce_domain <-
#   function(x.dim,
#            y.dim,
#            ex.max = 5,
#            cov.args,
#            cov.function = "stationary.cov",
#            return.setup = TRUE
#            ) {
#     p2 <- 2^(1:ex.max)
#     for(i in seq_along(p2)) {
#       M <- p2[i] * c(x.dim, y.dim)
#       dx <- c(1,1)
#       sim.grid <- as.matrix(expand.grid(1:M[1], 1:M[2]))
#       sim.center <- rbind(round(M/2))
#       cov.out <- do.call(cov.function, c(cov.args, list(x1 = sim.grid, x2 = sim.center)))
#       cov.out <- array(c(cov.out), M)
#       arr.center <- array(0, M)
#       arr.center[rbind(sim.center)] <- 1
#       weights <- fft(cov.out)/(fft(arr.center) * prod(M))
#       if(!any(Re(weights) < 0)) break
#     }
#     if(!(any(Re(weights) < 0))) {
#       if(return.setup) {
#         ce.setup <- list(m = c(x.dim, y.dim),
#                          grid = list(1:x.dim, 1:y.dim),
#                          dx = c(1, 1),
#                          M = M,
#                          wght = weights,
#                          call = NULL)
#         return(ce.setup)
#       } else {
#       return(M)
#       }
#     } else {
#       return(NULL)
#     }
#   }


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
                              dist.n,
                              dist.acc = 1,
                              mix = NULL,
                              mix.prop = 1,
                              rescale = c(0,1),
                              name = "value") {
  poly.lim <- limit_polygon(x.dim, y.dim)
  target <- st_sample(poly.lim, dist.n) |> st_union()
  if(is.null(mix)) {
    mix <- generate_empty(x.dim = x.dim, y.dim = y.dim) + 1
  }
  mix.dt <- as.data.table(as.data.frame(mix))
  nuclei.id <- sample(1:nrow(mix.dt), size = dist.n, prob = mix.dt$value)
  target <-
    mix.dt[nuclei.id, .(x, y)] |>
    as.matrix() |>
    st_multipoint() |>
    st_sfc()
  x.crd <- seq(0 + dist.acc/2, x.dim - dist.acc/2, dist.acc)
  y.crd <- seq(0 + dist.acc/2, y.dim - dist.acc/2, dist.acc)
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
  dist[[1]][which(is.na(dist[[1]]))] <- 0
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


generate_areas_poly.old <- function(x.dim,
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
  if(verbose) message("Setting up candidate polygons …")

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



decode_seg <- function(x, n = 1, digits = 0) {
  seg.list <- split(round(x, digits = digits),
                    f = rep(1:5, times = n))
  names(seg.list) <- c("x.seg", "y.seg", "x.cen", "y.cen", "buffer")
  return(seg.list)
}


# REMOVE
# x <- decode_seg(x)
# x.dt <- as.data.table(x)
# x.dt[,
#      `:=`(grid.seg = res.grid$grid[st_intersects(res.grid,
#                                                  st_point(c(x.seg, y.seg)),
#                                                  sparse = FALSE)],
#           grid.cen = res.grid$grid[st_intersects(res.grid,
#                                                  st_point(c(x.cen, y.cen)),
#                                                  sparse = FALSE)]),
#      by = .I]


# st_as_sf(st

# res.grid$grid[st_intersects(res.grid, st_point(c(5, 25)), sparse = FALSE)]

# REMOVE
# pw.dist = grid.dist.pw
# segments = x.dt$grid.seg
# centers = x.dt$grid.cen
# buffer = x.dt$buffer
# include.unmatched = FALSE
# plot(res.grid[res.grid$grid %in% c(segments, grid),])

match_to_segments <- function(pw.dist,
                              segments,
                              centers,
                              buffer,
                              min.dist = 0,
                              include.unmatched = TRUE) {
    stopifnot(length(segments) == length(buffer) &
              length(centers) == length(buffer))
    n.seg <- length(segments)
    matched.l <- list()
    grid.seg <-
      data.table(grid = 1:nrow(pw.dist),
                 segment =
                   apply(pw.dist[, segments, drop = FALSE],
                         1, which.min) |>
                   lapply(\(x) ifelse(length(x) > 0, x, 0)) |>
                   unlist())
    grid.seg[, dist := pw.dist[grid, centers[segment]], by = .I]
    matched.l <- list()
    for(i in 1:n.seg) {
      matched.l[[i]] <- grid.seg[segment == i & dist <= buffer[i], grid]
    }
    matched.grid <- do.call(c, matched.l)
    matched.segment <-
      do.call(c, mapply(\(x, y) rep(y, length(x)),
                                         matched.l,
                                         1:n.seg,
                                         SIMPLIFY = FALSE))
    if(include.unmatched) {
      unmatched.grid <-  (1:nrow(pw.dist))[!1:nrow(pw.dist) %in% matched.grid]
      unmatched.segment <- rep(0, length(unmatched.grid))
      matched.grid <- c(matched.grid, unmatched.grid)
      matched.segment <- c(matched.segment, unmatched.segment)
    }
    return(list(grid = matched.grid, segment = matched.segment))
  }

#3 REMOVE
# matched <-
# match_to_segments(pw.dist = grid.dist.pw,
# segments = x.dt$grid.seg,
# centers = x.dt$grid.cen,
# buffer = x$buffer,
# include.unmatched = TRUE) |>
# as.data.table()

# plot(cbind(res.grid, merge(res.grid.dt, matched)[order(grid), segment]))

# system.time(segment_min_distance(grid.dist.pw, grid.sel = grid.sel, segments = segments))

segment_min_distance.old <- function(pw.dist,
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

segment_min_distance <- function(pw.dist,
                                 grid.sel,
                                 segments
                                 ) {
  grid.seg <- data.table(grid = grid.sel, segment = segments)
  grid.seg[, `:=`(grid1 = grid,
                  grid2 = grid,
                  segment1 = segment,
                  segment2 = segment)]
  if(is.data.table(pw.dist)) {
    min.dist <-
      pw.dist[CJ(grid1 = grid.sel, grid2 = grid.sel),
              on = c("grid1", "grid2")
                ][grid.seg[, .(grid1, segment1)],
                  on = "grid1"
                  ][grid.seg[, .(grid2, segment2)],
                    on = "grid2"] |>
      DT(segment1 != segment2,
         .(dist.min = min(dist)),
         by = c("segment1", "segment2")) |>
      DT(order(segment1, segment2))
  }
  if(is.matrix(pw.dist)) {
    grid.dist.dt <- as.data.table(pw.dist[grid.sel, grid.sel])
    names(grid.dist.dt) <- as.character(grid.sel)
    grid.dist.dt[, grid1 := grid.sel]
    min.dist <-
      melt(grid.dist.dt,
           id.vars = "grid1",
           variable.name = "grid2", value.name = "dist",
           variable.factor = FALSE) |>
      DT(, grid2 := as.integer(grid2)) |>
      merge(grid.seg[, .(grid1, segment1)],
            by = "grid1") |>
      merge(grid.seg[, .(grid2, segment2)],
            by = "grid2") |>
      DT(segment1 != segment2,
         .(dist.min = min(dist)),
         by = c("segment1", "segment2")) |>
      DT(order(segment1, segment2))
  }
  return(min.dist)
}


encode_grid_buffer.old <- function(grid, buffer, order.grid, order.buffer) {
  grid.gray <- lapply(grid, \(x) binary2gray(decimal2binary(x, order.grid)))
  buffer.gray <- lapply(buffer, \(x) binary2gray(decimal2binary(x, order.buffer)))
  string <- c(do.call(c, grid.gray), do.call(c, buffer.gray))
  return(string)
}

encode_grid_buffer <- function(grid.seg, grid.cen, buffer, order.grid, order.buffer) {
  grid.gray <- lapply(c(grid.seg, grid.cen), \(x) binary2gray(decimal2binary(x, order.grid)))
  buffer.gray <- lapply(buffer, \(x) binary2gray(decimal2binary(x, order.buffer)))
  string <- c(do.call(c, grid.gray), do.call(c, buffer.gray))
  return(string)
}

decode_grid_buffer <- function(string, bit.orders) {
  string <- split(string, rep.int(seq.int(bit.orders), times = bit.orders))
  orders <- sapply(string, function(x) { binary2decimal(gray2binary(x)) })
  decoded <- split(unname(orders), f = rep(1:3, each = length(bit.orders)/3))
  return(list(grid.seg = decoded[[1]],
              grid.cen = decoded[[2]],
              buffer = decoded[[3]]))
}

decode_grid_buffer.old <- function(string, bit.orders) {
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

gmean <- function(x, na.rm = FALSE, zero.allow = TRUE) {
  if(any(x == 0)) {
    if(zero.allow) {
      return(0)
    } else {
      stop("Values of `x` must be strictly positive when `zero.allow = FALSE`.")
    }
  }
  gmean <- exp(mean(log(x), na.rm = na.rm))
  return(gmean)
}



generate_areas_poly <- function(x.dim,
                                y.dim,
                                seg.seed,
                                score,
                                score.sam = 1e4,
                                imbalance = 0.1,
                                imbalance.tol = NULL,
                                area.prop = 0.5,
                                area.tol = NULL,
                                area.exact = FALSE,
                                seg.res = 0.1 * min(x.dim, y.dim),
                                seg.min.dist = 0,
                                seg.min.area = 0,
                                seg.even = 2,
                                seg.prec = 0.1 * seg.res,
                                min.bound.dist = 0,
                                verbose = 2,
                                opt.imp.imb = 1,
                                opt.imp.even = 1,
                                opt.imp.area = 1,
                                opt.imb.agg = gmean,
                                opt.pop.n = 50,
                                opt.pop.rand = TRUE,
                                opt.prec = 1e-3,
                                opt.pcrossover = 0.8,
                                opt.pmutation = 0.1,
                                opt.elitism = max(1, round(opt.pop.n*0.1)),
                                opt.max.iter = 100,
                                opt.run = opt.max.iter,
                                opt.parallel = FALSE,
                                opt.fine = TRUE,
                                opt.fine.max.iter = 25,
                                opt.fine.constr = TRUE,
                                opt.fine.tol = 1e-4,
                                opt.cache = TRUE,
                                ...
                                ) {

  if(verbose > 0) message("Setting up candidate polygons …")

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
  # score <- score[1]
  # names(score)[1] <- "score"
  # score <- scale_int(score, int = c(0, 1))

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
  grid.dist.dt[grid1 == grid2, dist := -1e-6]
  setindex(grid.dist.dt, grid1)
  setindex(grid.dist.dt, grid2)


  cov.names <- paste0("cov", 1:length(score))
  score <- setNames(score, cov.names)
  cov.dt <- as.data.table(score)

  seg.cov <- 
    st_rasterize(res.grid[, "grid"],
                 template = generate_empty(x.dim = x.dim, y.dim = y.dim)) |>
    as.data.table() |>
    merge(as.data.table(cov.dt), all.x = FALSE, all.y = TRUE)
  if(score.sam < nrow(seg.cov)) {
    seg.cov <- seg.cov[sample(1:nrow(seg.cov), score.sam)][order(x, y)]
  }
  seg.cov <-
    melt(seg.cov,
         measure.vars = cov.names,
         variable.name = "var", value.name = "value")
  setindex(seg.cov, grid)

  seg.cov.obs <- seg.cov[var == cov.names[1], .(x, y, grid)]
  setindex(seg.cov.obs, grid)

  cov.bins <- list()
  for(i in seq_along(cov.names)) {
    cov.breaks <- get_breaks(seg.cov[var == cov.names[i], value], breaks = "FD")
    seg.cov[var == cov.names[i],
            value.bin := assign_breaks(value, breaks = cov.breaks)]
    cov.bins[[i]] <- 
      expand.grid(var = cov.names[i],
                  grid = unique(seg.cov$grid),
                  value.bin = 1:length(cov.breaks)) |>
      as.data.table()
  }

  cov.bins <- rbindlist(cov.bins)

  seg.cov.bin <-
    merge(seg.cov[, .(n.bin = .N), by = c("var", "grid", "value.bin")],
          cov.bins,
          by = c("var", "grid", "value.bin"),
          all = TRUE)
  seg.cov.bin[is.na(n.bin), n.bin := 0]
  setindex(seg.cov.bin, var)
  setindex(seg.cov.bin, grid)
  setindex(seg.cov.bin, value.bin)

  get_freq <- function(x) {
    foc.freq <- 
      seg.cov.bin[.(x),
                    .(n = sum(n.bin)),
                    by = c("value.bin", "var"),
                    on = "grid"
                    ][!is.na(value.bin)
                      ][order(var, value.bin),
                        .(value.bin,
                          foc = n/sum(n)),
                        by = "var"]
    ref.freq <- 
      seg.cov.bin[!.(x),
                    .(n = sum(n.bin)),
                    by = c("value.bin", "var"),
                    on = "grid"
                    ][!is.na(value.bin)
                      ][order(var, value.bin),
                        .(value.bin,
                          ref = n/sum(n)),
                        by = "var"]
    return(list(foc.freq = foc.freq, ref.freq = ref.freq))
  }

  if(verbose > 0) message("Building areas (optimization) …")

  order.grid <- ceiling(log2(nrow(res.grid.dt)))
  order.buffer <- ceiling(log2(sqrt(x.dim^2 + y.dim^2)))
  order.c <- rep(c(order.grid, order.grid, order.buffer), each = seg.seed)
  n.bits <- sum(order.c)

  f_opt_imb_area <- function(x) {
    x.decode <- decode_grid_buffer(x, order.c)
    x.dt <- as.data.table(x.decode)
    x.dt[!grid.seg %in% res.grid$grid, grid.seg := NA]
    x.dt[!grid.cen %in% res.grid$grid, grid.cen := NA]
    # If coordinates out of bounds, return penalty
    if(any(is.na(x.dt$grid.seg)) | any(is.na(x.dt$grid.cen))) {
      f <-
        sum(c(is.na(x.dt[, c(grid.seg, grid.cen)])))^2 *
        sqrt(.Machine$double.xmax)
    } else {
      trt.proposed <-
        match_to_segments.m(pw.dist = grid.dist.pw,
                            segments = x.dt$grid.seg,
                            centers = x.dt$grid.cen,
                            buffer = x.dt$buffer,
                            include.unmatched = FALSE) |>
        as.data.table()
      # If empty treatment area, return penalty
      if(nrow(trt.proposed) == 0) {
        f <- length(x.decode$grid.seg)^2 * sqrt(.Machine$double.xmax)
      # Otherwise evaluate function
      } else {
        trt.a <- 
          merge(res.grid.dt[.(trt.proposed$grid),
                            on = "grid"],
                trt.proposed) |>
          DT(, .(area = sum(area)), by = "segment")
        # Any remaining segments without assigned grid polygons?
        seg.rem <- setdiff(1:seg.seed, trt.a$segment)
        if(length(seg.rem) > 0) {
          trt.a <-
            rbind(trt.a,
                  data.table(segment = seg.rem,
                             area = rep(0, length(seg.rem))))
        }
        trt.area <- sum(trt.a) / prod(x.dim, y.dim)
        penalty.min.area <- pmax(-trt.a$area + seg.min.area, 0) * sqrt(.Machine$double.xmax)
        if(length(unique(trt.proposed$segment)) > 1) {
          min.dist <-
            segment_min_distance.m(pw.dist = grid.dist.dt,
                                   grid.sel = trt.proposed$grid,
                                   segments = trt.proposed$segment)
          penalty.min.dist <-
            min.dist[,
                     pen.dist := fifelse(dist.min < seg.min.dist,
                                         sqrt(.Machine$double.xmax),
                                         0)
                     ][, pen.dist]
        } else {
          penalty.min.dist <- 0
        }
        penalties <- sum(penalty.min.area) + sum(penalty.min.dist)
        if(penalties > 0) {
          f <- penalties
        } else {
          freqs <- get_freq.m(trt.proposed$grid)
          if(nrow(freqs$ref.freq) > 0) {
            js.dt <-
              merge(freqs$foc.freq, freqs$ref.freq, by = c("var", "value.bin"))
            trt.imb <-
              js.dt[, 
                    .(js_div = js_div(foc, ref, type = "raw")),
                    by = "var"
                    ][,js_div]
          } else {
            trt.imb <- rep(0, length(unique(freqs$foc.freq$var)))
          }
          imbalance.loss <- (opt.imb.agg(trt.imb) - imbalance)^2
          even.loss <- var(trt.imb)
          area.loss <- (trt.area - area.prop)^2
          # REMOVE
          # penalty.oob <- as.numeric(!(x$grid %in% res.grid.dt$grid)) * sqrt(.Machine$double.xmax)
          f <-
            opt.imp.imb * imbalance.loss +
            opt.imp.even * even.loss +
            opt.imp.area * area.loss +
            penalties
        }
      }
    }
    return(-f)
  }
  
  if(opt.cache) {
    # Cache function evaluations
    f_opt_imb_area.m <- memoise(f_opt_imb_area)
    match_to_segments.m <- memoise(match_to_segments)
    segment_min_distance.m <- memoise(segment_min_distance)
    get_freq.m <- memoise(get_freq)
  } else {
    f_opt_imb_area.m <- f_opt_imb_area
    match_to_segments.m <- match_to_segments
    segment_min_distance.m <- segment_min_distance
    get_freq.m <- get_freq
  }

  # Starting population

  if(verbose > 0) message("Creating starting population for genetic optimization …")

  pop.start.mat <- matrix(NA, nrow = opt.pop.n, ncol = n.bits)

  if(!opt.pop.rand & verbose > 1) {
    prog <- txtProgressBar(min = 0, max = opt.pop.n, initial = 0,
                           char = "=", width = NA, style = 3)
  }
  n.pop <- 1
  while(n.pop <= opt.pop.n) {
    grid.seg.sam <- sample(res.grid.dt$grid, seg.seed)
    grid.cen.sam <- sample(res.grid.dt$grid, seg.seed, replace = TRUE)
    buffer.sam <- round(runif(n = seg.seed, 0, sqrt(x.dim^2 + y.dim^2)))
    pop.try <- encode_grid_buffer(grid.seg.sam,
                                  grid.cen.sam,
                                  buffer.sam,
                                  order.grid,
                                  order.buffer)
    if(opt.pop.rand) {
      fit.try = 0
    } else {
      fit.try <- f_opt_imb_area.m(pop.try)
    }
    if(fit.try > -sum(c(opt.imp.imb, opt.imp.even, opt.imp.area))) {
      pop.start.mat[n.pop,] <- pop.try
      if(!opt.pop.rand & verbose > 1) {
        setTxtProgressBar(prog, n.pop)
      }
      n.pop <- n.pop + 1
    }
  }
  if(!opt.pop.rand & verbose > 1) close(prog)

  # pop.start.mat <- matrix(NA, nrow = opt.pop.n, ncol = n.bits)
  # for(i in 1:opt.pop.n) {
  #   grid.seg.sam <- sample(res.grid.dt$grid, seg.seed)
  #   grid.cen.sam <-
  #     grid.dist.dt[grid2 %in% grid.seg.sam,
  #                  .SD[which.min(dist)][1],
  #                  by = grid1,
  #                  .SDcols = c("grid2", "dist")
  #                  ][,
  #                    .SD[sample(1:.N, 1)],
  #                    by = grid2
  #                    ][.(grid.seg.sam), grid1, on = "grid2"]
  #   buffer.sam <- round(runif(n = seg.seed,
  #                             sqrt(seg.min.area),
  #                             sqrt(area.prop * x.dim * y.dim / seg.seed)))
  #   pop.start.mat[i,] <-
  #     encode_grid_buffer(grid.seg.sam,
  #                        grid.cen.sam,
  #                        buffer.sam,
  #                        order.grid,
  #                        order.buffer)
  # }

  f_max <- -((opt.imp.imb + opt.imp.even + opt.imp.area) * opt.prec)
  f_monitor <- ifelse(verbose > 1, gaMonitor, FALSE)

  if(verbose > 0) message("Genetic optimization …")

  buffer.opt <- NULL
  while(is.null(buffer.opt)) {
    try({
      buffer.opt <-
          ga(type =  "binary",
             fitness = f_opt_imb_area.m,
             maxFitness = f_max,
             # lower = f_lower,
             # upper = f_upper,
             nBits = n.bits,
             suggestions = pop.start.mat,
             run = opt.run,
             popSize = opt.pop.n,
             maxiter = opt.max.iter,
             optim = TRUE,
             optimArgs = list(poptim = 0.25),
             parallel = opt.parallel,
             pcrossover = opt.pcrossover,
             pmutation = opt.pmutation,
             elitism = opt.elitism,
             monitor = f_monitor,
             ...)
    })
    if(verbose > 0 & is.null(buffer.opt)) message("Optimization failed. Trying again …")
  }

  imb_area.sol <- decode_grid_buffer(buffer.opt@solution[1,], order.c)

  trt.proposed <-
    match_to_segments.m(pw.dist = grid.dist.pw,
                        segments = imb_area.sol$grid.seg,
                        centers = imb_area.sol$grid.cen,
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

  # Extract results of coarse optimization

  grid.trt <- with(res.grid, grid[selected])
  grid.ref <- setdiff(res.grid$grid, grid.trt)
  trt.a <- 
    res.grid.dt[.(grid.trt),
                .(area = sum(area)),
                by = "segment",
                on = "grid"
                ][, area]

  trt.area <-
    res.grid.dt[selected == TRUE, sum(area)] / res.grid.dt[, sum(area)]
  trt.imb <-
    seg.cov[,
            .(js = js_div(x = .SD[grid %in% grid.trt, value],
                          y = .SD[!grid %in% grid.trt, value],
                          type = "discrete")),
            by = "var",
            .SDcols = c("grid", "value")
            ][, js]

  # remove cached functions
  rm(f_opt_imb_area.m, match_to_segments.m, segment_min_distance.m)

  if(verbose > 0) {
    message(paste0("Imbalance (approx.) ",
                   paste(round(trt.imb, getOption("digits")), collapse = ", "),
                   "; aggregate ",
                   round(opt.imb.agg(trt.imb), getOption("digits")),
                   ".\nArea proportion (approx.) ",
                   round(trt.area, getOption("digits")),
                   "."))
  }

  if(opt.fine) {

    if(verbose > 0) message("Attempt to further adjust areas (fine tuning) …")
    penalty <- sqrt(opt.imp.imb + opt.imp.even + opt.imp.area)
    pen.imb = ifelse(opt.imb.agg(trt.imb) < imbalance.lim[1] |
                     opt.imb.agg(trt.imb) > imbalance.lim[2],
                     abs(opt.imb.agg(trt.imb) - imbalance.lim) * sqrt(penalty), 0)
    pen.area = ifelse(trt.area < area.lim[1] | trt.area > area.lim[2],
                      abs(trt.area - area.lim) * sqrt(penalty), 0)
    pen.min.area <- sum(pmax(-trt.a + seg.min.area, 0) * sqrt(penalty))
    pen.min.dist <-
      grid.dist.dt[grid1 %in% grid.trt & grid2 %in% grid.trt,
                   .(dist.min = min(dist)),
                   by = c("segment1", "segment2")
                   ][segment1 != segment2,
                     sum(0.5 * as.numeric(dist.min < seg.min.dist) * sqrt(penalty))]
    imb.loss <- (opt.imb.agg(trt.imb) - imbalance)^2
    even.loss <- var(trt.imb)
    area.loss <- (trt.area - area.prop)^2
    loss.null <- opt.imp.imb * imb.loss +
                 opt.imp.even * even.loss +
                 opt.imp.area * area.loss +
                 pen.min.area + pen.min.dist +
                 pen.imb + pen.area

    A.nb <-
      poly2nb(res.grid, queen = FALSE) |>
      nb2listw(style = "B") |> 
      listw2mat() |>
      Matrix()

    # A.nb <- st_relate(res.grid, res.grid, pattern = "F***1****", sparse = FALSE) * 1
    # diag(A.nb) <- 0
    # A.nb <- as(A.nb, "sparseMatrix")

#     A.nb <-
#         poly2nb(res.grid, queen = FALSE) |>
#         nb2listw(style = "B") |> 
#         listw2mat()
    loss.hist <- list()

    for(i in 1:opt.fine.max.iter) {
      # if(trt.area <= area.prop) {
        grid.adj.add <-
          grid.ref[rowSums(A.nb[grid.ref, -grid.ref]) > 0]
          # grid.ref[rowSums(A.nb[grid.ref, -grid.ref]) > 0 &
          #          grid.ref %in% with(res.grid, grid[segment != 0])]
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
          foc.obs <- seg.cov.obs[.(idx.trt), .N, on = "grid"]
          rem.obs <- seg.cov.obs[!.(idx.trt), .N, on = "grid"]
          foc.a <- 
            res.grid.dt[.(idx.trt),
                        .(area = sum(area)),
                        by = "segment",
                        on = "grid"
                        ][, area]
          idx.area <- sum(foc.a) / res.grid.dt[, sum(area)]
          idx.js <-
            seg.cov[,
                    .(js = js_div(x = .SD[grid %in% idx.trt, value],
                                  y = .SD[!grid %in% idx.trt, value],
                                  type = "discrete")),
                    by = "var",
                    .SDcols = c("grid", "value")
                    ][, js]
          idx.imb <- opt.imb.agg(idx.js)
          idx.even <- sd(idx.js)
          pen.min.dist <-
            grid.dist.dt[grid1 %in% idx.trt & grid2 %in% idx.trt,
                         .(dist.min = min(dist)),
                         by = c("segment1", "segment2")
                         ][segment1 != segment2,
                           0.5 * as.numeric(dist.min < seg.min.dist) * sqrt(penalty)]
          if(foc.obs > 0 & rem.obs > 0) {
            adj.add[[j]] <- data.table(grid = grid.adj.add[j],
                                       # smd = d_cohen(foc.obs, rem.obs),
                                       imb = idx.imb,
                                       even = idx.even,
                                       area = idx.area,
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
          foc.obs <- seg.cov.obs[.(idx.trt), .N, on = "grid"]
          rem.obs <- seg.cov.obs[!.(idx.trt), .N, on = "grid"]
          foc.a <- 
            res.grid.dt[.(idx.trt),
                        .(area = sum(area)),
                        by = "segment",
                        on = "grid"
                        ][, area]
          idx.area <- sum(foc.a) / res.grid.dt[, sum(area)]
          idx.js <-
            seg.cov[,
                    .(js = js_div(x = .SD[grid %in% idx.trt, value],
                                  y = .SD[!grid %in% idx.trt, value],
                                  type = "discrete")),
                    by = "var",
                    .SDcols = c("grid", "value")
                    ][, js]
          idx.imb <- opt.imb.agg(idx.js)
          idx.even <- sd(idx.js)
          pen.min.dist <-
            grid.dist.dt[grid1 %in% idx.trt & grid2 %in% idx.trt,
                         .(dist.min = min(dist)),
                         by = c("segment1", "segment2")
                         ][segment1 != segment2,
                           0.5 * as.numeric(dist.min < seg.min.dist) * sqrt(penalty)]
          if(foc.obs > 0 & rem.obs > 0) {
            adj.rem[[j]] <- data.table(grid = grid.adj.rem[j],
                                       # smd = d_cohen(foc.obs, rem.obs),
                                       imb = idx.imb,
                                       even = idx.even,
                                       area = idx.area,
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
      adj.dt[, `:=`(loss = opt.imp.imb * (imb - imbalance)^2 +
                           opt.imp.even * even^2 +
                           opt.imp.area * (area - area.prop)^2 +
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
          if(verbose > 0) {
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
      if(verbose > 0 & i == opt.fine.max.iter) {
        message(paste0("Maximum number of iterations reached (", opt.fine.max.iter,
                       "). Further optimization may be possible."))
      }
    }

  } # End fine-scale optimization


  res.grid$selected[-grid.trt] <- FALSE
  res.grid$selected[grid.trt] <- TRUE

  res.grid.sel <- res.grid[res.grid$selected,]

  poly.trt <-
    with(res.grid.sel, aggregate(geometry, by = list(segment = segment), st_union)) |>
    st_as_sf()

  if(trt.area != area.prop & area.exact == TRUE) {
    if(verbose > 0) message("Resizing area … ")

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
  st_agr(poly.areas) <- "constant"

  areas.dt <-
    st_geometry(poly.areas[poly.areas$type == "treatment",]) |>
    st_cast("POLYGON") |>
    st_as_sf() |>
    st_set_geometry("geometry") |>
    st_rasterize(template = generate_empty(x.dim = x.dim,
                                           y.dim = y.dim)) |>
    as.data.table() |>
    setnames(c("x", "y", "poly"))
  areas.dt[, type := fifelse(poly == 0, "reference", "treatment")]
  areas.dt[, type := factor(type, c("reference", "treatment"))]

  trt.area <- areas.dt[type == "treatment", length(type) / nrow(areas.dt)]

  trt.loc <- areas.dt[type == "treatment", .(x, y)]

  trt.obs <-
    cov.dt[trt.loc,
           ..cov.names,
           on = c("x", "y")] |>
    as.list()
  ref.obs <-
    cov.dt[!trt.loc,
           ..cov.names,
           on = c("x", "y")] |>
    as.list()
  trt.imb <-
    mapply(\(x, y) js_div(x, y, type = "discrete"), trt.obs, ref.obs) |>
    unname()


  if(verbose > 0) {
    message(paste0("Final result (exact):\nImbalance ",
                   paste(round(trt.imb, getOption("digits")), collapse = ", "),
                   "; aggregate ",
                   round(opt.imb.agg(trt.imb), getOption("digits")),
                   ".\nArea proportion ",
                   round(trt.area, getOption("digits")),
                   "."))
  }

  return(list(shape = poly.areas,
              table = areas.dt[, .(x, y, poly, type)],
              performance = c(imbalance = trt.imb,
                              area.prop = trt.area)))
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
                        fbm.alpha = 1,
                        fbm.var = 1,
                        fbm.scale = 1,
                        fbm.w = 1,
                        fbm.rescale = c(0, 1),
                        mix = NULL,
                        mix.w = 0,
                        mix.rescale = c(0, 1),
                        rescale = c(0, 1),
                        name = "z1") {
  fbm <- generate_fbm(x.dim = x.dim,
                       y.dim = y.dim,
                       alpha = fbm.alpha,
                       var = fbm.var,
                       scale = fbm.scale,
                       rescale = fbm.rescale)
  if(is.null(mix)) {
    z <- fbm
  } else {
  z <- fbm.w * fbm + mix.w * scale_int(mix, mix.rescale)
  }
  z <-
    scale_int(z, int = rescale) |>
    setNames(name)
}



generate_z2 <- function(x.dim,
                        y.dim,
                        grad.phi = 0,
                        grad.shift = 0,
                        grad.w = 0,
                        grad.x.scale = 1,
                        grad.y.scale = 1,
                        grad.rescale = c(0, 1),
                        mix = NULL,
                        mix.w = 1,
                        mix.rescale = c(0, 1),
                        rescale = c(0, 1),
                        name = "z2") {
  grad <- generate_exp_gradient(x.dim = x.dim,
                                y.dim = y.dim,
                                phi = grad.phi,
                                shift = grad.shift,
                                x.scale = grad.x.scale,
                                y.scale = grad.y.scale,
                                rescale = c(0, 1))
  if(is.null(mix)) {
    z <- grad
  } else {
  z <- grad.w * grad + mix.w * scale_int(mix, mix.rescale)
  }
  z <-
    scale_int(z, int = rescale) |>
    setNames(name)
}

generate_z3 <- function(
                        x.dim,
                        y.dim,
                        dist.n = 10,
                        dist.acc = 1,
                        mix = NULL,
                        mix.prop = 1,
                        rescale = c(0,1),
                        name = "z3"
                        ) {
  dist <- generate_distance(x.dim = x.dim,
                            y.dim = y.dim,
                            dist.n = dist.n,
                            dist.acc = dist.acc,
                            mix = mix,
                            mix.prop = mix.prop,
                            rescale = FALSE)
  # z <- interpolate_na(dist, method = "gam", nh = "neumann", range = 2)
  z <-
    scale_int(dist, int = rescale) |>
    setNames(name)
  return(z)
}

generate_z4 <- function(x.dim,
                        y.dim,
                        seg.n = 3,
                        mat.nu = 1,
                        mat.var = 1,
                        mat.scale = 1,
                        mat.mean = 0,
                        mat.w = 1,
                        mat.rescale = c(0, 1),
                        mix = NULL,
                        mix.w = 0,
                        mix.rescale = c(0, 1),
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
                    mean = mat.mean,
                    rescale = mat.rescale) |>
    segment_field(seg) |>
    scale_int()
    if(is.null(mix)) {
      z <- mat
    } else {
    mix.seg <-
      mix |>
      segment_field(seg) |>
      scale_int(mix.rescale)
    z <- mat.w * mat + mix.w * mix.seg
    }
    z <-
      scale_int(z, int = rescale) |>
      setNames(name)
  return(z)
}


# REMOVE
# generate_treatment.old <- function(x.dim,
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
#                                name = "treatment") {
#   treatment <- generate_empty(x.dim = x.dim,
#                               y.dim = y.dim)
#   centroids <-
#     st_centroid(poly[1:2,]) |>
#     st_coordinates() |>
#     suppressWarnings()
#   phi.shift <-
#     atan2((centroids[2,2] - centroids[1,2]),
#           (centroids[2,1] - centroids[1,1]))
#   phi.range.trt <- phi.range + phi.shift
#   phi.range.ctr <- phi.range + phi.shift + pi
#   phi.trt <-
#     sample(seq(phi.range.trt[1], phi.range.trt[2], length.out = 2 * nuclei), nuclei)
#   phi.ctr <-
#     sample(seq(phi.range.ctr[1], phi.range.ctr[2], length.out = 2 * nuclei), nuclei)
#   # phi.trt <- runif(n,
#   #                  phi.range[1] + phi.shift,
#   #                  phi.range[2] + phi.shift)
#   # phi.ctr <- runif(n,
#   #                  phi.range[1] + phi.shift + pi,
#   #                  phi.range[2] + phi.shift + pi)
#   shift.trt <- runif(nuclei,
#                      shift.range[1],
#                      shift.range[2])
#   shift.ctr <- runif(nuclei,
#                      shift.range[1],
#                      shift.range[2])
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
#                             phi = phi.trt[i],
#                             shift = shift.trt[i],
#                             x.scale = x.scale.trt[i],
#                             y.scale = y.scale.trt[i])
#     ctr.grad <-
#       generate_exp_gradient(x.dim = x.dim,
#                             y.dim = y.dim,
#                             phi = phi.ctr[i],
#                             shift = shift.ctr[i],
#                             x.scale = x.scale.ctr[i],
#                             y.scale = y.scale.ctr[i])
#     treatment <- trt.eff[i] * trt.grad - ctr.eff[i] * ctr.grad + treatment
#   }
#   eff.scale <-
#     effect.size.sp /
#     (mean(treatment[poly[2,]][[1]], na.rm = TRUE) -
#      mean(treatment[poly[1,]][[1]], na.rm = TRUE))
#   treatment <- treatment * eff.scale
#   eff.shift <- mean(treatment[poly[1,]][[1]], na.rm = TRUE)
#   treatment <- treatment - eff.shift
#   treatment <- setNames(treatment, name)
#   poly.bd <-
#     st_geometry(poly[2,]) |>
#     st_as_sf()
#   poly.bd$eff.bd <- effect.size.bd
#   treatment.bd <- st_rasterize(poly.bd, template = generate_empty(x.dim = x.dim, y.dim = y.dim))
#   treatment <- treatment + reset_dim(treatment.bd, treatment)
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


generate_treatment <- function(x.dim,
                               y.dim,
                               trt.poly,
                               eff.mean = 1,
                               mat.nu = 1,
                               mat.scale = max(x.dim, y.dim) / 10,
                               mat.var = 1,
                               damp.type = "asymmetric",
                               damp.infl = 0.2,
                               damp.scale = 1/(5*pi),
                               name = "treatment"
                               ) {
 
  # Prepare treatment areas

  trt.poly <-
    st_geometry(trt.poly) |>
    st_as_sf() |>
    st_set_geometry("geometry")

  trt.dt <-
    trt.poly |>
    st_rasterize(template = generate_empty(x.dim = x.dim,
                                           y.dim = y.dim)) |>
    as.data.table() |>
    setnames(c("x", "y", "poly"))

  # Set up grid and calculate distances to polygon boundaries and
  # centroids.

  grid.coord <-
    expand.grid(x = (1:x.dim)-0.5,
                y = (1:y.dim)-0.5)

  grid.pts <-
    st_as_sf(grid.coord, coords = c(1, 2))

  dist.cen <-
    st_distance(st_centroid(trt.poly), grid.pts) |>
    apply(2, min)

  dist.lim <-
    st_cast(trt.poly, "POLYGON") |>
    st_cast("LINESTRING") |>
    st_distance(grid.pts) |>
    apply(2, min)

  dist.dt <-
    cbind(grid.coord, dist.cen, dist.lim) |>
    as.data.table()

  # Calculate dampening due to edge effects

  if(damp.type == "asymmetric") {
    dist.dt[, damp := inv_cloglog(((dist.lim / (dist.cen + dist.lim)) - damp.infl) / damp.scale)]
  }
  if(damp.type == "symmetric") {
    dist.dt[, damp := plogis(dist.lim / (dist.cen + dist.lim),
                             location = damp.infl,
                             scale = damp.scale)]
  }

  dist.dt <- merge(dist.dt, trt.dt)

  dist.dt[poly == 0, damp := 0]

  # Generate raw treatment surface

  matern <-
    generate_matern(x.dim = x.dim,
                    y.dim = y.dim,
                    mean = eff.mean,
                    nu = mat.nu,
                    scale = mat.scale,
                    var = mat.var,
                    name = "matern",
                    rescale = FALSE)

  matern <- matern - min(matern$matern, na.rm = TRUE)
  matern.dt <- as.data.table(matern)

  dist.mat.dt <- merge(dist.dt, matern.dt)

  # Dampen and scale

  eff.raw <- dist.mat.dt[poly != 0, mean(matern * damp)]
  dist.mat.dt[, treatment := (eff.mean/eff.raw) * matern * damp]

  treat <- st_as_stars(dist.mat.dt[, .(x, y, treatment)])
  names(treat) <- name

  return(treat)
}


generate_nonlinear_effect <- function(field,
                                      type = "random",
                                      range = 1,
                                      mu = 0,
                                      f.acc = 0.01) {
  f.min <- min(field[[1]], na.rm = TRUE)
  f.max <- max(field[[1]], na.rm = TRUE)
  f.range <- f.max - f.min
  field.sc <- (field - f.min) / f.range
  if(type == "random") {
    type <- sample(c("sigmoid", "minimum", "unimodal", "bimodal"), 1)
  }
  if(type == "sigmoid") {
    ex <- runif(1, 0, 2)
    sign <- sample(c(-1, 1), 1)
    par <- sign * 10^ex
    fn <- function(x) 1/(1+exp(-par*(x-0.5)))
  }
  if(type == "minimum") {
    par <- runif(1, -3/4, -1/4) * pi
    # par <- c(0.25, 1)
    fn <- function(x) sin(0.5*pi * (x-0.5) + par[1])
  }
  if(type == "unimodal") {
    par <- c(runif(1, exp(1)/4, exp(1)), runif(1, 0, 1))
    fn <- function(x) exp(-((par[1] * (x - par[2]))^2))
  }
  if(type == "bimodal") {
    par <- c(
             runif(1, 1, 2),
             runif(1, 1, 2),
             runif(1, exp(1), 2*exp(1)),
             runif(1, exp(1), 2*exp(1)),
             runif(1, 0, 0.25),
             runif(1, 0.75, 1),
             runif(1, 0.5, 1) * sample(c(-1, 1), 1))
    fn <- function(x) {
      par[1] * exp(-((par[3] * (x - par[5]))^2)) +
      par[2] * exp(-((par[4] * (x - par[6]))^2)) +
      par[7] * x
    }
  }
  effect.sc <- field.sc
  effect.sc[[1]] <- fn(effect.sc[[1]])
  eff.min <- min(effect.sc[[1]], na.rm = TRUE)
  eff.max <- max(effect.sc[[1]], na.rm = TRUE)
  eff.range <- eff.max - eff.min
  effect <- (range / eff.range) * (effect.sc - eff.min)
  eff.mean <- mean(effect[[1]], na.rm = TRUE)
  effect <- effect - eff.mean + mu
  fun <- data.table(x = seq(0, 1, f.acc))
  fun[, f.x := (range / eff.range) * (fn(x) - eff.min) - eff.mean + mu]
  fun[, x := f.min + (f.range * x)]
  setnames(fun, c("x", "f.x"), c(names(field), paste0("f.", names(field)))) 
  return(list(effect = effect, fun = fun))
}


generate_interaction_effect <- function(field1,
                                        field2,
                                        range,
                                        mu,
                                        nuclei,
                                        f.acc = 0.01) {
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
  f1.min <- min(field1[[1]], na.rm = TRUE)
  f1.max <- max(field1[[1]], na.rm = TRUE)
  f1.range <- f1.max - f1.min
  f2.min <- min(field2[[1]], na.rm = TRUE)
  f2.max <- max(field2[[1]], na.rm = TRUE)
  f2.range <- f2.max - f2.min
  field1.sc <- (field1 - f1.min) / f1.range
  field2.sc <- (field2 - f2.min) / f2.range
  fields.c <- c(field1.sc[1], field2.sc[1])
  names(fields.c) <- c("field1", "field2")
  fields.dt <- as.data.table(fields.c)
  fields.dt[,
            `:=`(x = f.acc * round(field1 / f.acc, 0),
                 y = f.acc * round(field2 / f.acc, 0))]
  int.ex <- unique(fields.dt[, .(x, y, exists = TRUE)])
  int.fun <- CJ(x = seq(min(fields.c[[1]], na.rm = TRUE),
                        max(fields.c[[1]], na.rm = TRUE),
                        f.acc),
                y = seq(min(fields.c[[2]], na.rm = TRUE),
                        max(fields.c[[2]], na.rm = TRUE),
                        f.acc))
  int.fun <-
    merge(int.fun, int.ex, by = c("x", "y"), all.x = TRUE)
  int.fun[is.na(exists), exists := FALSE]
  int.fun.com <- matrix(NA, nrow = nrow(int.fun), ncol = nuclei)
  for(i in 1:nrow(centres)){
    int.effect[[1]] <-
      int.effect[[1]] +
      fn(field1.sc[[1]], field2.sc[[1]],
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
  f.min <- min(int.effect[[1]], na.rm = TRUE)
  f.max <- max(int.effect[[1]], na.rm = TRUE)
  f.range <- f.max - f.min
  int.effect <- (range / f.range) * (int.effect - f.min)
  f.mean <- mean(int.effect[[1]])
  int.effect <- int.effect - f.mean + mu
  int.fun[, f.xy := apply(int.fun.com, 1, sum)]
  int.fun[, f.xy := (range / f.range) * (f.xy - f.min) - f.mean + mu]
  int.fun[,
          `:=`(x = f1.min + (f1.range * x),
               y = f2.min + (f2.range * y))]
  setnames(int.fun, c(names(field1[1]), names(field2[1]),
                      "exists",
                      paste0("f.", names(field1[1]), names(field2[1]))))
  return(list(effect = int.effect, fun = int.fun))
}


generate_landscape_4cov_nl <-
  function(seed = NULL,
           x.dim,
           y.dim,
           mix.alpha,
           z1.fbm.alpha,
           z1.mix.w,
           z2.grad.phi,
           z2.mix.w,
           z3.dist.n,
           z3.dist.acc,
           z3.mix.prop,
           z4.seg.n,
           z4.mat.nu,
           z4.mat.scale,
           z4.mix.w,
           areas.seg.seed,
           areas.score.sam,
           areas.imbalance,
           areas.imbalance.tol,
           areas.area.prop,
           areas.area.tol,
           areas.area.exact,
           areas.seg.res,
           areas.seg.min.dist,
           areas.seg.min.area,
           areas.seg.even,
           areas.seg.prec,
           areas.min.bound.dist,
           areas.opt.imp.imb,
           areas.opt.imp.even,
           areas.opt.imp.area,
           areas.opt.imb.agg,
           areas.opt.pop.n,
           areas.opt.pop.rand,
           areas.opt.prec,
           areas.opt.pcrossover,
           areas.opt.pmutation,
           areas.opt.elitism,
           areas.opt.max.iter,
           areas.opt.run,
           areas.opt.parallel,
           areas.opt.fine,
           areas.opt.fine.max.iter,
           areas.opt.fine.constr,
           areas.opt.fine.tol,
           areas.opt.cache,
           treatment.eff.mean,
           treatment.mat.nu,
           treatment.mat.scale,
           treatment.mat.var,
           treatment.damp.type,
           treatment.damp.infl,
           treatment.damp.scale,
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
           e.exp.var,
           e.exp.scale,
           e.nug.var,
           verbose = 2,
           ...) {

  ls.par <- as.list(match.call(expand.dots=FALSE))
  set.seed(seed)
 
  if(verbose > 0) message("Simulating covariates …") 
 
  fbm.mix <- generate_fbm(x.dim = x.dim,
                          y.dim = y.dim,
                          alpha = mix.alpha,
                          var = 1,
                          scale = 1,
                          rescale = c(0,1),
                          name = "mix")

  z1 <- generate_z1(x.dim = x.dim,
                    y.dim = y.dim,
                    fbm.alpha = z1.fbm.alpha,
                    fbm.var = 1,
                    fbm.scale = 1,
                    fbm.w = 1-z1.mix.w,
                    fbm.rescale = FALSE,
                    mix = scale_int((sample(c(-1, 1), 1) * fbm.mix)),
                    mix.w = z1.mix.w,
                    rescale = c(0,1),
                    name = "z1")

  z2 <- generate_z2(
                    x.dim = x.dim,
                    y.dim = y.dim,
                    grad.phi = z2.grad.phi,
                    grad.shift = 0,
                    grad.x.scale = 2,
                    grad.y.scale = 2,
                    grad.w = 1-z2.mix.w,
                    mix = scale_int((sample(c(-1, 1), 1) * fbm.mix)),
                    mix.w = z2.mix.w,
                    rescale = c(0,1),
                    name = "z2")

  z3 <- generate_z3(x.dim = x.dim,
                    y.dim = y.dim,
                    dist.n = z3.dist.n,
                    dist.acc = z3.dist.acc,
                    mix = scale_int((sample(c(-1, 1), 1) * fbm.mix)),
                    mix.prop = z3.mix.prop,
                    rescale = c(0,1),
                    name = "z3")

  z4 <- generate_z4(x.dim = x.dim,
                    y.dim = y.dim,
                    seg.n = z4.seg.n,
                    mat.nu = z4.mat.nu,
                    mat.scale = z4.mat.scale,
                    mat.var = 1,
                    mat.mean = 0,
                    mat.rescale = FALSE,
                    mat.w = 1-z4.mix.w,
                    mix = scale_int((sample(c(-1, 1), 1) * fbm.mix)),
                    mix.w = z4.mix.w,
                    rescale = c(0,1),
                    name = "z4")

  covariates <- c(z1, z2, z3, z4)


  if(verbose > 0) message("Simulating treatment areas …") 

  areas.poly <-
    generate_areas_poly(x.dim = x.dim,
                        y.dim = y.dim,
                        seg.seed = areas.seg.seed,
                        # score = score.trt,
                        score = covariates,
                        score.sam = areas.score.sam,
                        imbalance = areas.imbalance,
                        imbalance.tol = areas.imbalance.tol,
                        area.prop = areas.area.prop,
                        area.tol = areas.area.tol,
                        area.exact = areas.area.exact,
                        seg.res = areas.seg.res,
                        seg.min.dist = areas.seg.min.dist,
                        seg.min.area = areas.seg.min.area,
                        seg.even = areas.seg.even,
                        seg.prec = areas.seg.prec,
                        min.bound.dist = areas.min.bound.dist,
                        verbose = verbose,
                        opt.imp.imb = areas.opt.imp.imb,
                        opt.imp.even = areas.opt.imp.even,
                        opt.imp.area = areas.opt.imp.area,
                        opt.imb.agg = areas.opt.imb.agg,
                        opt.pop.n = areas.opt.pop.n,
                        opt.pop.rand = areas.opt.pop.rand,
                        opt.prec = areas.opt.prec,
                        opt.pcrossover = areas.opt.pcrossover,
                        opt.pmutation = areas.opt.pmutation,
                        opt.elitism = areas.opt.elitism,
                        opt.max.iter = areas.opt.max.iter,
                        opt.run = areas.opt.run,
                        opt.parallel = areas.opt.parallel,
                        opt.fine = areas.opt.fine,
                        opt.fine.max.iter = areas.opt.fine.max.iter,
                        opt.fine.constr = areas.opt.fine.constr,
                        opt.fine.tol = areas.opt.fine.tol,
                        opt.cache = areas.opt.cache)

  trt.poly <- 
    areas.poly$shape[areas.poly$shape$type == "treatment",] |>
    st_cast("POLYGON")


  if(verbose > 0) message("Simulating treatment effect …") 

  treatment <-
    generate_treatment(x.dim = x.dim,
                       y.dim = y.dim,
                       trt.poly = trt.poly,
                       eff.mean = treatment.eff.mean,
                       mat.nu = treatment.mat.nu,
                       mat.scale = treatment.mat.scale,
                       mat.var = treatment.mat.var,
                       damp.type = treatment.damp.type,
                       damp.infl = treatment.damp.infl,
                       damp.scale = treatment.damp.scale,
                       name = "treatment"
                       )

  if(verbose > 0) message("Simulating covariate effects …") 

  cov.main.names <- paste0("z", 1:4)
  cov.main.effect.type <- ls.par[paste0(cov.main.names, ".effect.type")]
  cov.main.effect.range <- ls.par[paste0(cov.main.names, ".effect.range")]
  cov.main.effect.mu <- ls.par[paste0(cov.main.names, ".effect.mu")]
  cov.main.effects <- list()
  cov.main.funs <- list()
  for(i in seq_along(cov.main.names)){
    nl.effect <-
      generate_nonlinear_effect(field = covariates[i],
                                type = cov.main.effect.type[[i]],
                                range = cov.main.effect.range[[i]],
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
    setnames(int.effect$fun, c("val1", "val2", "exists", "f"))
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


  if(verbose > 0) message("Simulating residual variation …") 

  mod.residual <-
    RMexp(var = e.exp.var, scale = e.exp.scale) +
    RMnugget(var = e.nug.var)
  residual <- RFsimulate(mod.residual, x = 1:x.dim, y = 1:y.dim)
  residual <- reset_dim(st_as_stars(residual))
  residual <- setNames(residual, "residual")


  if(verbose > 0) message("Building landscape …") 

  # response <-
  #   c(cov.main.effects, treatment, residual) |>
  #   merge() |>
  #   st_apply(1:2, sum) |>
  #   setNames("response")

  response <-
    c(cov.main.effects, cov.int.effects, treatment, residual) |>
    merge() |>
    st_apply(1:2, sum) |>
    setNames("response")

  landscape.dt <-
    c(response,
      covariates,
      treatment,
      cov.main.effects,
      cov.int.effects,
      residual) |>
    as.data.table()

  landscape.dt <-
    merge(areas.poly$table,
          landscape.dt,
          by = c("x", "y"))

  landscape.dt$cell <- 1:nrow(landscape.dt)

  setcolorder(landscape.dt, c("x", "y", "cell"))

  fun <- list(main = cov.main.funs, interactions = cov.int.funs)
 
  return(list(landscape = landscape.dt, fun = fun, optim = areas.poly$performance))
  }



generate_landscape_4cov_nl_binary <- function(ls,
                                              p.ref,
                                              p.trt,
                                              e.var,
                                              e.df,
                                              opt.grid,
                                              verbose = 2,
                                              ...) {

  if(verbose > 0) message("Scaling treatment effect …") 

  cov.effects <- names(ls)[names(ls) %like% "f."]
  cov.n <- length(cov.effects)

  cov.mat.trt <-
    ls[type == "treatment",
       c("treatment", cov.effects),
       with = FALSE] |>
    as.matrix()

  f_opt_b <- function(b, a, n, p.ref, p.trt) {
    trt.scale <- c(b, rep(1, n))
    trt.raw <- as.vector(cov.mat.trt %*% trt.scale)
    p.trt.s <- sum(trt.raw >= a) / length(trt.raw)
    f <- (p.trt.s - p.trt)^2
    return(f)
  }

  a <- quantile(cov.mat.trt %*% c(0, rep(1, cov.n)), 1-p.ref, names = FALSE)

  resp.r <- range(ls$response)
  resp.m <- median(ls$response)
  opt.lim <- (max(abs(resp.r)) * c(-1, 1))

  b.seq <- seq(opt.lim[1], opt.lim[2], length.out = opt.grid)
  f.seq <- sapply(b.seq, f_opt_b, a = a, n = cov.n, p.ref = p.ref, p.trt = p.trt)

  f_opt_sp <- splinefun(b.seq, f.seq)

  df_opt_sp <- function(x, ...) {
    f_opt_sp(x, deriv = 1)
  }

  if(verbose > 1) opt.trace <- TRUE else opt.trace <- TRUE

  opt.b <-
    optim(par = 1,
          fn = f_opt_sp,
          gr = df_opt_sp,
          method = "L-BFGS-B",
          lower = opt.lim[1],
          upper = opt.lim[2],
          control = list(factr = 1,
                         trace = opt.trace))

  if(verbose > 0) message("Building landscape …") 

  b <- opt.b$par

  pars <- c(a = a, b = b)

  ls.bin <- copy(ls)
  ls.bin[, `:=`(intercept = -pars[1],
                treatment = pars[2] * treatment,
                residual = qt(pnorm(residual,
                                    mean = 0,
                                    sd = sqrt(e.var)),
                              df = e.df))]

  lat.form <-
    parse(text = paste(c("intercept", "treatment", cov.effects, "residual"),
                       collapse = " + "))

  ls.bin[, latent := eval(lat.form)]
  ls.bin[, response := fifelse(latent > 0, 1, 0)]

  trt.scale <- c(b, rep(1, cov.n))
  trt.raw <- as.vector(cov.mat.trt %*% trt.scale)
  mar.trt <- sum(trt.raw >= a) / length(trt.raw)

  ref.form <-
    parse(text = paste(c("intercept", cov.effects),
                       collapse = " + "))
  trt.form <-
    parse(text = paste(c("intercept", "treatment", cov.effects),
                       collapse = " + "))

  mar.trt <-
    ls.bin[type == "treatment",
           .(p.ref = sum(eval(ref.form) > 0) / .N,
             p.trt = sum(eval(trt.form) > 0) / .N)]

  n.poly.trt <- ls.bin[type == "treatment", length(unique(poly))]

  if(n.poly.trt > 1) {
  mar.trt <-
    rbind(mar.trt,
          ls.bin[type == "treatment",
                 .(p.ref = sum(eval(ref.form) > 0) / .N,
                   p.trt = sum(eval(trt.form) > 0) / .N),
                 by = poly][order(poly)],
          fill = TRUE)
  setcolorder(mar.trt, c("poly", "p.ref", "p.trt"))
  }

  return(list(landscape = ls.bin,
              scaling = pars,
              marginal = mar.trt))

}


plot_landscape_4cov_nl <- function(x,
                                   select = "overview",
                                   interactions = TRUE,
                                   ls_theme = NULL,
                                   ls_guide = NULL,
                                   title.type = "Area type",
                                   legend.type = "Area type",
                                   title.cov = "Covariates",
                                   legend.cov = "Covariate\nvalue",
                                   title.resp = "Response",
                                   legend.resp = "Response",
                                   title.trt = "Treatment effect",
                                   legend.trt = "Treatment\neffect",
                                   title.resid = "Residual variation",
                                   legend.resid = "Residual\nvariation",
                                   title.eff.main = "Covariate effects (main)",
                                   legend.eff.main = "Covariate\neffect",
                                   title.eff.int = "Covariate effects (interactions)",
                                   legend.eff.int = "Covariate\neffect",
                                   limits.common = TRUE
                                   ) {
  # select options for landscape: overview, effects, all

  if(is.null(ls_theme)) {
    ls_theme <-
      theme_minimal(base_family = "IBMPlexSansCondensed", base_size = 11) +
      theme(plot.title = element_text(hjust = 0,
                                      face = "bold",
                                      margin = margin(l = 0, b = 5, t = 11)),
            plot.margin = margin(3, 3, 3, 3),
            strip.text = element_text(size = rel(0.8),
                                      hjust = 0.5),
            strip.background = element_rect(fill = "grey90", color = NA),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "right",
            legend.justification = "top")
  }

  if(is.null(ls_guide)) {
    guide_fill <-
      guides(fill = guide_colorbar(
                                   ticks.colour = "grey5",
                                   ticks.linewidth = 0.2,
                                   frame.colour = "grey5",
                                   frame.linewidth = 0.2,
                                   barwidth = 1,
                                   barheight = 5,
                                   label.position = "right",
                                   label.hjust = 1,
                                   draw.ulim = TRUE,
                                   draw.llim = TRUE
                                   ))
  } else {
    guide_fill <- ls_guide
  }

  if(select == "interactions") interactions <- TRUE
  
  plots <- list()

  # Response
  
  plots[["resp"]] <-
    ggplot(x, aes(x = x, y = y, fill = response)) +
    geom_raster() +
    # scale_fill_viridis_c()
    scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    guide_fill +
    labs(title = title.resp, fill = legend.resp) +
    ls_theme

  # Variables: area type

  plots[["type"]] <- 
    ggplot(x, aes(x = x, y = y)) +
    geom_raster(aes(fill = type)) +
    # geom_text(data = lab.coord,
    #           aes(x = x, y = y,
    #               # label = capwords(as.character(type)),
    #               color = type), size = 3)  +
    scale_fill_manual(breaks = levels(x$type),
                      values = c("grey20", "grey80"),
                      labels = capwords(levels(x$type))) +
    # scale_colour_manual(breaks = levels(x$landscape$type),
    #                     values = c("grey95", "grey5")) +
    # scale_x_continuous(expand = c(0, 0)) +
    # scale_y_continuous(expand = c(0, 0)) +
    coord_fixed(expand = FALSE) +
    labs(title = title.type, fill = legend.type) +
    # guides(color = "none") +
    ls_theme

  # Variables: covariates

  cov.var <- paste0("z", 1:4)
  cov.lab <- toupper(cov.var)
  names(cov.lab) <- cov.var

  cov.dt <-
    x[,
      c("x", "y", cov.var),
      with = FALSE] |>
    melt(id.vars = c("x", "y"),
         measure.vars = cov.var,
         variable.name = "covariate",
         value.name = "value")

  cov.dt[, label := cov.lab[covariate]]

  # for(i in 1:4) {
  #   p.title <- switch(i, `1` = "Covariates", " ")
  #   plots[[cov.var[i]]] <-
  #     cov.dt[covariate == cov.var[i]] |>
  #     ggplot() +
  #     geom_raster(aes(x = x, y = y, fill = value)) +
  #     scale_fill_continuous_sequential(palette = "Viridis", rev = FALSE, limits = c(0, 1)) +
  #     scale_x_continuous(expand = c(0, 0)) +
  #     scale_y_continuous(expand = c(0, 0)) +
  #     coord_fixed() +
  #     guide_fill +
  #     facet_wrap(vars(label)) +
  #     labs(title = p.title, fill = "Covariate\nvalue") +
  #     ls_theme
  # }

  plots[["cov"]] <-
    ggplot(cov.dt) +
    geom_raster(aes(x = x, y = y, fill = value)) +
    scale_fill_continuous_sequential(palette = "Viridis", rev = FALSE, limits = c(0, 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    guide_fill +
    facet_wrap(vars(label), ncol = 2, nrow = 2) +
    labs(title = title.cov, fill = legend.cov) +
    ls_theme


  # Effects: covariates

  eff.main <-  paste0("f.", cov.var)
  eff.int <- 
    paste0("f.int", c("12", "13", "14", "23", "24", "34"))
  eff.var <- c(eff.main, eff.int)

  eff.main.lab <- toupper(cov.var)
  names(eff.main.lab) <- eff.main
  eff.int.lab <-
    c("z1:z2", "z1:z3", "z1:z4",
      "z2:z3", "z2:z4", "z3:z4") |>
    toupper()
  names(eff.int.lab) <- eff.int
  eff.lab <- c(eff.main.lab, eff.int.lab)

  eff.dt <-
    x[,
      c("x", "y", eff.var),
      with = FALSE] |>
    melt(id.vars = c("x", "y"),
         measure.vars = eff.var,
         variable.name = "effect",
         value.name = "value")

  eff.dt[, label := eff.lab[effect]]

  if(limits.common) {
    lim.effects <-
      unlist(x[,
               .(min(as.matrix(.SD)),
                 max(as.matrix(.SD))),
               .SDcols = c("treatment",
                           "residual",
                           eff.var)],
             use.names = FALSE)
    lim.effects <- 
      (floor(lim.effects / 0.25) + c(0, 1)) * 0.25
  } else {
    lim.effects <- c(NA, NA)
  }

    # for(i in seq_along(eff.main)) {
    #   p.title <- switch(i, `1` = "Covariate effects (main)", " ")
    #   plots[[eff.main[i]]] <-
    #     eff.dt[effect == eff.main[i]] |>
    #     ggplot(aes(x = x, y = y, fill = value)) +
    #     geom_raster() +
    #     scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
    #                                      limits = lim.effects) +
    #     scale_x_continuous(expand = c(0, 0)) +
    #     scale_y_continuous(expand = c(0, 0)) +
    #     coord_fixed() +
    #     guide_fill +
    #     facet_wrap(vars(label)) +
    #     labs(title = p.title,
    #          fill = "Covariate\neffects") +
    #     ls_theme
    # }

    plots[["eff.main"]] <-
      eff.dt[effect %in% eff.main] |>
      ggplot(aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                       limits = lim.effects) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed() +
      guide_fill +
      facet_wrap(vars(label)) +
      labs(title = title.eff.main,
           fill = legend.eff.main) +
      ls_theme


    plots[["eff.int"]] <-
      eff.dt[effect %in% eff.int] |>
      ggplot(aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_continuous_divergingx(palette = "Roma",
                                       rev = TRUE, mid = 0,
                                       limits = lim.effects) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed() +
      guide_fill +
      facet_wrap(vars(label)) +
      labs(title = title.eff.int,
           fill = legend.eff.int) +
      ls_theme


    # Effects: treatment

    plots[["trt"]] <-
      ggplot(x, aes(x = x, y = y, fill = treatment)) +
      geom_raster() +
      scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                       limits = lim.effects) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed() +
      guide_fill +
      labs(title = title.trt, fill = legend.trt) +
      ls_theme

    # Effects: residual

    plots[["resid"]] <-
      ggplot(x, aes(x = x, y = y, fill = residual)) +
      geom_raster() +
      scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                       limits = lim.effects) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed() +
      guide_fill +
      labs(title = title.resid, fill = legend.resid) +
      ls_theme


    # for(i in seq_along(eff.int)) {
    #   p.title <- switch(i, `1` = "Covariate effects (interactions)", " ")
    #   plots[[eff.int[i]]] <-
    #     eff.dt[effect == eff.int[i]] |>
    #     ggplot(aes(x = x, y = y, fill = value)) +
    #     geom_raster() +
    #     scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
    #                                      limits = lim.effects) +
    #     scale_x_continuous(expand = c(0, 0)) +
    #     scale_y_continuous(expand = c(0, 0)) +
    #     coord_fixed() +
    #     guide_fill +
    #     facet_wrap(vars(label)) +
    #     labs(title = p.title,
    #          fill = "Covariate\neffects") +
    #     ls_theme
    # }


  if(select == "overview") {
    layout.overview <-
      "ACC
       BCC"
    p.select <-
      (plots$type +
      plots$resp +
      plots$cov +
      plot_layout(design = layout.overview, widths = 1, heights = 1, guides = "collect")) &
      ls_theme
  }

  if(select == "effects" & interactions) {
    layout.eff <-
      "ABB
       DBB
       CCC
       CCC"
    p.select <-
      (plots$trt +
      plots$eff.main +
      plots$eff.int +
      plots$resid +
      plot_layout(design = layout.eff, widths = 1, guides = "collect")) &
      ls_theme
  }
  

  if(select == "effects" & !interactions) {
    layout.eff.noint <-
      "ABB
       CBB"
    p.select <-
      (plots$trt +
      plots$eff.main +
      plots$resid +
      plot_layout(design = layout.eff.noint, widths = 1, guides = "collect")) &
      ls_theme
  }

  if(select == "all") {
    layout.all <-
      "ACC
       BCC
       DEE
       GEE
       FFF
       FFF"
    p.select <-
      (plots$type +
      plots$resp +
      plots$cov +
      plots$trt +
      plots$eff.main +
      plots$eff.int +
      plots$resid +
      plot_layout(design = layout.all, widths = 1, guides = "collect")) &
      ls_theme
  }

  if(select == "all" & !interactions) {
    layout.all.noint <-
      "ACC
       BCC
       DEE
       FEE"
    p.select <-
      (plots$type +
      plots$resp +
      plots$cov +
      plots$trt +
      plots$eff.main +
      plots$resid +
      plot_layout(design = layout.all.noint, widths = 1, guides = "collect")) &
      ls_theme
  }
  
  if(select == "list") {
    p.select <- plots
  }

  return(p.select)
}



plot_landscape_4cov_nl_binary <- function(x,
                                          select = "overview",
                                          interactions = TRUE,
                                          ls_theme = NULL,
                                          ls_guide = NULL,
                                          title.type = "Area type",
                                          legend.type = "Area type",
                                          title.cov = "Covariates",
                                          legend.cov = "Covariate\nvalue",
                                          title.resp = "Response",
                                          legend.resp = "Response",
                                          title.trt = "Treatment effect",
                                          legend.trt = "Treatment\neffect",
                                          title.resid = "Residual variation",
                                          legend.resid = "Residual\nvariation",
                                          title.eff.main = "Covariate effects (main)",
                                          legend.eff.main = "Covariate\neffect",
                                          title.eff.int = "Covariate effects (interactions)",
                                          legend.eff.int = "Covariate\neffect",
                                          limits.common = TRUE
                                          ) {
  # select options for landscape: overview, effects, all

  if(is.null(ls_theme)) {
    ls_theme <-
      theme_minimal(base_family = "IBMPlexSansCondensed", base_size = 11) +
      theme(plot.title = element_text(hjust = 0,
                                      face = "bold",
                                      margin = margin(l = 0, b = 5, t = 11)),
            plot.margin = margin(3, 3, 3, 3),
            strip.text = element_text(size = rel(0.8),
                                      hjust = 0.5),
            strip.background = element_rect(fill = "grey90", color = NA),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "right",
            legend.justification = "top")
  }

  if(is.null(ls_guide)) {
    guide_fill <-
      guides(fill = guide_colorbar(
                                   ticks.colour = "grey5",
                                   ticks.linewidth = 0.2,
                                   frame.colour = "grey5",
                                   frame.linewidth = 0.2,
                                   barwidth = 1,
                                   barheight = 5,
                                   label.position = "right",
                                   label.hjust = 1,
                                   draw.ulim = TRUE,
                                   draw.llim = TRUE
                                   ))
  } else {
    guide_fill <- ls_guide
  }

  if(select == "interactions") interactions <- TRUE
  
  plots <- list()

  # Response
 
  resp.col <- divergingx_hcl(11, "Roma", rev = TRUE)[c(6, 11)]
  names(resp.col) <- c("0", "1")

  plots[["resp"]] <-
    ggplot(x, aes(x = x, y = y,
                  fill = factor(as.character(response), levels = c("0", "1")))) +
    geom_raster() +
    scale_fill_manual(values = resp.col) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    labs(title = title.resp, fill = legend.resp) +
    ls_theme

  # Variables: area type

  plots[["type"]] <- 
    ggplot(x, aes(x = x, y = y)) +
    geom_raster(aes(fill = type)) +
    scale_fill_manual(breaks = levels(x$type),
                      values = c("grey20", "grey80"),
                      labels = capwords(levels(x$type))) +
    coord_fixed(expand = FALSE) +
    labs(title = title.type, fill = legend.type) +
    ls_theme

  # Variables: covariates

  cov.var <- paste0("z", 1:4)
  cov.lab <- toupper(cov.var)
  names(cov.lab) <- cov.var

  cov.dt <-
    x[,
      c("x", "y", cov.var),
      with = FALSE] |>
    melt(id.vars = c("x", "y"),
         measure.vars = cov.var,
         variable.name = "covariate",
         value.name = "value")

  cov.dt[, label := cov.lab[covariate]]

  plots[["cov"]] <-
    ggplot(cov.dt) +
    geom_raster(aes(x = x, y = y, fill = value)) +
    scale_fill_continuous_sequential(palette = "Viridis", rev = FALSE, limits = c(0, 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed() +
    guide_fill +
    facet_wrap(vars(label), ncol = 2, nrow = 2) +
    labs(title = title.cov, fill = legend.cov) +
    ls_theme

  # Effects: covariates

  eff.main <-  paste0("f.", cov.var)
  eff.int <- 
    paste0("f.int", c("12", "13", "14", "23", "24", "34"))
  eff.var <- c(eff.main, eff.int)

  eff.main.lab <- toupper(cov.var)
  names(eff.main.lab) <- eff.main
  eff.int.lab <-
    c("z1:z2", "z1:z3", "z1:z4",
      "z2:z3", "z2:z4", "z3:z4") |>
    toupper()
  names(eff.int.lab) <- eff.int
  eff.lab <- c(eff.main.lab, eff.int.lab)

  eff.dt <-
    x[,
      c("x", "y", eff.var),
      with = FALSE] |>
    melt(id.vars = c("x", "y"),
         measure.vars = eff.var,
         variable.name = "effect",
         value.name = "value")

  eff.dt[, label := eff.lab[effect]]

  if(limits.common) {
    lim.effects <-
      unlist(x[,
               .(min(as.matrix(.SD)),
                 max(as.matrix(.SD))),
               .SDcols = c("treatment",
                           eff.var)],
             use.names = FALSE)
    lim.effects <- 
      (floor(lim.effects / 0.25) + c(0, 1)) * 0.25
  } else {
    lim.effects <- c(NA, NA)
  }

    plots[["eff.main"]] <-
      eff.dt[effect %in% eff.main] |>
      ggplot(aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                       limits = lim.effects) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed() +
      guide_fill +
      facet_wrap(vars(label)) +
      labs(title = title.eff.main,
           fill = legend.eff.main) +
      ls_theme


    plots[["eff.int"]] <-
      eff.dt[effect %in% eff.int] |>
      ggplot(aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_continuous_divergingx(palette = "Roma",
                                       rev = TRUE, mid = 0,
                                       limits = lim.effects) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed() +
      guide_fill +
      facet_wrap(vars(label)) +
      labs(title = title.eff.int,
           fill = legend.eff.int) +
      ls_theme

    # Effects: treatment

    plots[["trt"]] <-
      ggplot(x, aes(x = x, y = y, fill = treatment)) +
      geom_raster() +
      scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                       limits = lim.effects) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed() +
      guide_fill +
      labs(title = title.trt, fill = legend.trt) +
      ls_theme

    # Effects: residual

    res.lim <- quantile(x$residual, c(0.01, 0.99))

    plots[["resid"]] <-
      ggplot(x, aes(x = x, y = y, fill = residual)) +
      geom_raster() +
      scale_fill_continuous_divergingx(palette = "Roma", rev = TRUE, mid = 0,
                                       limits = res.lim,
                                       oob = scales::squish
                                       ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed() +
      guide_fill +
      labs(title = title.resid, fill = legend.resid) +
      ls_theme

  if(select == "overview") {
    layout.overview <-
      "ACC
       BCC"
    p.select <-
      (plots$type +
      plots$resp +
      plots$cov +
      plot_layout(design = layout.overview, widths = 1, heights = 1, guides = "collect")) &
      ls_theme
  }

  if(select == "effects" & interactions) {
    layout.eff <-
      "ABB
       DBB
       CCC
       CCC"
    p.select <-
      (plots$trt +
      plots$eff.main +
      plots$eff.int +
      plots$resid +
      plot_layout(design = layout.eff, widths = 1, guides = "collect")) &
      ls_theme
  }
  

  if(select == "effects" & !interactions) {
    layout.eff.noint <-
      "ABB
       CBB"
    p.select <-
      (plots$trt +
      plots$eff.main +
      plots$resid +
      plot_layout(design = layout.eff.noint, widths = 1, guides = "collect")) &
      ls_theme
  }

  if(select == "all") {
    layout.all <-
      "ACC
       BCC
       DEE
       GEE
       FFF
       FFF"
    p.select <-
      (plots$type +
      plots$resp +
      plots$cov +
      plots$trt +
      plots$eff.main +
      plots$eff.int +
      plots$resid +
      plot_layout(design = layout.all, widths = 1, guides = "collect")) &
      ls_theme
  }

  if(select == "all" & !interactions) {
    layout.all.noint <-
      "ACC
       BCC
       DEE
       FEE"
    p.select <-
      (plots$type +
      plots$resp +
      plots$cov +
      plots$trt +
      plots$eff.main +
      plots$resid +
      plot_layout(design = layout.all.noint, widths = 1, guides = "collect")) &
      ls_theme
  }
  
  if(select == "list") {
    p.select <- plots
  }

  return(p.select)
}




plot_functions_4cov_nl <- function(x = NULL,
                                   main = NULL,
                                   select = "all",
                                   interactions = NULL,
                                   exist.only = FALSE,
                                   title.main = "Covariate effects (main)",
                                   title.int = "Covariate effects (interactions)",
                                   legend.int = "Partial\neffect",
                                   contour.bins = 11,
                                   fun_theme = NULL,
                                   fun_guide = NULL) {

  if(!is.null(x)) {
    main <- x[[1]]
    interactions <- x[[2]]
  }

  if(select != "list") {
    if(!is.null(main) & is.null(interactions))
      select <- "main"
    if(is.null(main) & !is.null(interactions))
      select <- "interactions"
  }


  lim.f <- range(c(main$f, interactions$f), na.rm = TRUE)
  lim.f <- (floor(lim.f / 0.25) + c(0, 1)) * 0.25

  if(is.null(fun_theme)) {
    fun_theme <-
      theme_light(base_family = "IBMPlexSansCondensed", base_size = 11) +
        theme(
              plot.title = element_text(hjust = 0,
                                        face = "bold",
                                        margin = margin(l = 0, b = 11, t = 11)),
              axis.line.x = element_line(color = "black",
                                         linewidth = rel(0.5)),
              axis.line.y = element_line(color = "black",
                                         linewidth = rel(0.5)),
              axis.title.x = element_text(margin = margin(t = 7)),
              axis.title.y = element_text(margin = margin(r = 7)),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              plot.margin = margin(3, 3, 3, 3),
              strip.text = element_text(size = rel(0.8),
                                        hjust = 0.5,
                                        color = "black",
                                        margin = margin(6, 6, 6, 6)),
              strip.background = element_rect(fill = "gray90"))
  }

  if(is.null(fun_guide)) {
    guide_fill <-
      guides(fill = guide_colorbar(ticks.colour = "grey5",
                                   ticks.linewidth = 0.2,
                                   frame.colour = "grey5",
                                   frame.linewidth = 0.2,
                                   barwidth = 1,
                                   barheight = 5,
                                   label.position = "right",
                                   label.hjust = 1,
                                   draw.ulim = TRUE,
                                   draw.llim = TRUE))
  } else {
    guide_fill <- fun_guide
  }


  plots <- list()

  if(!is.null(main)) {
    lim.main <- range(main$f, na.rm = TRUE)
    lim.main <- (floor(lim.main / 0.25) + c(0, 1)) * 0.25
    cov.lab <- toupper(main$cov)
    plots[["main"]] <-
      cbind(main, cov.lab) |>
      ggplot() +
        geom_line(aes(x = val, y = f)) +
        facet_wrap(vars(cov.lab), nrow = 2) +
        scale_y_continuous(lim = lim.main) +
        labs(title = title.main,
             x = "Covariate value",
             y = "Partial effect") +
        guide_fill +
        fun_theme
  }

  if(!is.null(interactions)) {
    if(exist.only) {
      interactions <- interactions[exists == TRUE]
    }
    lim.int <- range(interactions$f, na.rm = TRUE)
    lim.int <- (floor(lim.int / 0.25) + c(0, 1)) * 0.25
    cov1.lab <- toupper(interactions$cov1)
    cov2.lab <- toupper(interactions$cov2)
    plots[["int"]] <-
      cbind(interactions, cov1.lab, cov2.lab) |>
      ggplot() +
        aes(x = val1, y = val2) +
        geom_raster(aes(fill = f)) +
        geom_contour(aes(z = f),
                     bins = contour.bins,
                     colour = "gray35",
                     linewidth = 0.2) +
        facet_grid(cols = vars(cov1.lab),
                   rows = vars(cov2.lab)) +
        scale_fill_continuous_diverging(palette = "Tropic",
                                        lim = lim.int) +
        scale_x_continuous(n.breaks = 3) +
        scale_y_continuous(n.breaks = 3) +
        labs(title = title.int,
             x = "First covariate",
             y = "Second covariate",
             fill = legend.int) +
        guide_fill +
        fun_theme
  }

  if(select == "list") {
    return(plots)
  } else {
    if(select == "all") {
      p.select <- plots$main + plots$int + plot_layout(nrow = 2)
    }
    if(select == "main") {
      p.select <- plots$main
    }
    if(select == "interactions") {
      p.select <- plots$int
    }
    return(p.select)
  }
}



zplot <- function(z, n.colours = 128) {
  if(is.null(n.colours)) {
    n.colours <- min(length(unique(z[[1]])), 128)
  }
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



## EGP functions

som_init_pc <- function(grid, data) {
  # Calculate principal components
  init.pca <- prcomp(x = data, center = FALSE, scale = FALSE)
  init.max <- apply(init.pca$x[, 1:2], 2, max)
  init.min <- apply(init.pca$x[, 1:2], 2, min)
  grid.min <- apply(grid$pts, 2, min)
  grid.max <- apply(grid$pts, 2, max)
  grid.scale <- grid.max - grid.min
  init.scale <- (init.max - init.min) / grid.scale
  # Scale coordinates
  init.coord.pc <- t(apply(grid$pts, 1, \(x) init.scale * (x - grid.min) + init.min))
  # Map to covariate space
  init.coord.cov <- init.coord.pc %*% t(init.pca$rotation[,1:2])
  attr(init.coord.cov, "som.dim") <- c(grid$xdim, grid$ydim)
  return(init.coord.cov)
}

egp_som <- function(x,
                    x.dim,
                    y.dim = x.dim,
                    vars = NULL,
                    scale = TRUE,
                    topo = "hexagonal",
                    nb.fun = "gaussian",
                    dist.fun = "sumofsquares",
                    radius = max(x.dim, y.dim),
                    epochs = 1000,
                    init = som_init_pc,
                    mode = "pbatch",
                    parallel = 1) {

  if(inherits(x, "matrix")) {
    if(!is.null(vars)) {
      x <- x[,vars]
    }
  } else {
    if(is.null(vars)) {
      x <- as.matrix(as.data.frame(x))
    } else {
      x <- as.matrix(as.data.frame(x)[,vars])
    }
  }

  x <- na.omit(x)

  if(scale) {
    x.scaled <- scale(x, center = TRUE, scale = TRUE)
  }

  scale.att <- list(mean = apply(x, 2, mean, na.rm = TRUE),
                    sd = apply(x, 2, sd, na.rm = TRUE))

  som.grid <- somgrid(xdim = x.dim,
                      ydim = y.dim, 
                      topo = topo)

  if(topo == "hexagonal") {
    grid.poly <- grid_hex_poly(centroids = som.grid$pts)
  } else {
    grid.poly <- grid_rect_poly(centroids = som.grid$pts)
  }

  grid.nb <-
    poly2nb(grid.poly, queen = FALSE) |>
    nb2listw(style = "B") |> 
    listw2mat() |>
    Matrix()

  if(epochs * nrow(x.scaled) < 1e6)
    warning("Number of epochs may be too small for the number of observations.")

  if(is.function(init)) {
    init <- init(som.grid, x.scaled)
  }
 
  som.fit <- som(x.scaled,
                 grid = som.grid, 
                 rlen = epochs,
                 radius = radius,
                 init = init,
                 dist.fcts = dist.fun,
                 mode = "pbatch", 
                 cores = parallel,
                 normalizeDataLayers = FALSE)
  
  if(dist.fun == "sumofsquares") {
    dist.mat <-
      dist(som.fit$codes[[1]], method = "euclidean") |>
      as.matrix()
    dimnames(dist.mat) <- list(NULL, NULL)
    dist.mat <- dist.mat^2
  }

  if(dist.fun %in% c("euclidean", "manhattan")) {
    dist.mat <-
      dist(som.fit$codes[[1]], method = dist.fun) |>
      as.matrix()
    dimnames(dist.mat) <- list(NULL, NULL)
  }

  if(dist.fun == "tanimoto") {
    dist.mat <-
      dist(som.fit$codes[[1]], method = "binary") |>
      as.matrix()
    dimnames(dist.mat) <- list(NULL, NULL)
    dist.mat <- dist.mat^2
  }

  som.fit$unit.dist <-
    apply(dist.mat, 1, \(x) data.table(id.to = 1:length(x),
                                       dist = x)) |>
    rbindlist(idcol = "id.from")

  som.fit$grid.poly <- grid.poly
  som.fit$grid.nb <- grid.nb
  som.fit$init <- init
  som.fit$scale <- scale.att

  class(som.fit) <- c(class(som.fit), "egp_som")

  return(som.fit)
}



grid_hex_poly <- function(centroids = NULL,
                          x.dim = NULL,
                          y.dim = NULL,
                          d = 1,
                          phi = 0) {
  if(phi != 0) {
    rot <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), ncol = 2)
  }
  if(is.null(centroids)) {
    x <- 1L:x.dim
    y <- 1L:y.dim
    centroids <- as.matrix(expand.grid(x = x, y = y))
    centroids[, 1L] <- centroids[, 1L] + 0.5 * (centroids[, 2L]%%2)
    centroids[, 2L] <- sqrt(3)/2 * centroids[, 2L]
    if(phi != 0) {
      centroids <- centroids %*% t(rot)
    }
  }
  R <- d/sqrt(3)
  c30 <- cos(pi/6)
  s30 <- sin(pi/6)
  corners <-
    matrix(c(c30, s30,
             0, 1,
             -c30, s30,
             -c30, -s30,
             0, -1,
             c30, -s30,
             c30, s30),
          byrow = TRUE, nrow = 7) * R
  if(phi != 0) {
    corners <- corners %*% t(rot)
  }
  hex <- 
    apply(centroids, 1,
          \(x) 
          st_polygon(list(matrix(rep(x, 7), nrow = 7, byrow = TRUE) + corners)),
          simplify = FALSE) |>
    st_sfc() |>
    st_sf() |>
    st_set_geometry("geometry")
  hex$id <- 1:nrow(hex)
  return(hex)
}

grid_rect_poly <- function(centroids = NULL,
                           x.dim = NULL,
                           y.dim = NULL,
                           d = 1) {
  if(is.null(centroids)) {
    x <- 1L:x.dim
    y <- 1L:y.dim
    centroids <- as.matrix(expand.grid(x = x, y = y))
  }
  corners <-
    matrix(c(-0.5, -0.5,
              0.5, -0.5,
              0.5,  0.5,
             -0.5,  0.5,
             -0.5, -0.5),
           byrow = TRUE, nrow = 5) * d
  square <- 
    apply(centroids, 1,
          \(x) 
          st_polygon(list(matrix(rep(x, 5), nrow = 5, byrow = TRUE) + corners)),
          simplify = FALSE) |>
    st_sfc() |>
    st_sf() |>
    st_set_geometry("geometry")
  square$id <- 1:nrow(square)
  return(square)
}


egp_embed <- function(x,
                      som,
                      vars = NULL,
                      scale = TRUE,
                      bmu.name = "som.unit",
                      coord = TRUE,
                      coord.names = c("som.x", "som.y"),
                      dist = FALSE,
                      dist.name = "som.dist",
                      append = FALSE,
                      list = FALSE) {
  if(inherits(x, "matrix")) {
    append.mode <- "mat"
    if(!is.null(vars)) {
      x <- x[,vars]
    }
  } else {
    append.mode <- ifelse(is.data.table(x), "dt", "df")
    if(is.null(vars)) {
      x <- as.matrix(as.data.frame(x))
    } else {
      x <- as.matrix(as.data.frame(x)[,vars])
    }
  }
  x <- na.omit(x)
  if(scale) {
    x.scaled <- t(apply(x, 1, \(x) (x - som$scale$mean) / som$scale$sd))
  }
  mapped <- map(som, x.scaled)
  embedding <- mapped$unit.classif
  if(coord) {
    coord.mapped <- som$grid$pts[embedding,]
    colnames(coord.mapped) <- coord.names
    embedding <- cbind(embedding, coord.mapped)
  }
  if(dist) {
    dist.mapped <- mapped$distances
    embedding <- cbind(embedding, dist.mapped)
  }
  if(is.matrix(embedding)) {
    colnames(embedding)[1] <- bmu.name
  }
  if(append) {
    embedding <- cbind(x, embedding)
    if(append.mode %in% c("dt", "df")) {
      embedding <- as.data.frame(embedding)
    }
    if(append.mode == "dt") {
      setDT(embedding)
    }
  }
  if(list) {
    embedding <- as.list(as.data.frame(embedding))
  }
  return(embedding)
}

get_bmu <- function(x,
                    coord = TRUE,
                    bmu.name = "som.unit",
                    coord.names = c("som.x", "som.y"),
                    list = FALSE) {
  embedding <- x$unit.classif
  if(coord) {
    coord.mapped <- x$grid$pts[embedding,]
    colnames(coord.mapped) <- coord.names
    embedding <- cbind(embedding, coord.mapped)
  }
  if(is.matrix(embedding)) {
    colnames(embedding)[1] <- bmu.name
  }
  if(list) {
    embedding <- as.list(as.data.frame(embedding))
    if(coord == FALSE) names(embedding) <- bmu.name
  }
  return(embedding)
}

get_coord <- function(x) {x$grid$pts[x$unit.classif,]}
get_codes <- function(x) {x$codes}
get_grid <- function(x) {x$grid}
get_poly <- function(x) {x$grid.poly}
get_nb <- function(x) {x$grid.nb}


.ids_by_group <- function(data,
                          id.var,
                          group.vars = NULL,
                          group.labels = NULL,
                          add.label = TRUE,
                          expand.label = TRUE,
                          label.prefix = NULL,
                          ...){
  data.dt <- as.data.table(data)
  if(is.null(group.vars)) {
    groups <-
        data.dt[,
                .(id.col = list(id.col)),
                  env = list(id.col = id.var)]
  } else {
    groups <-
        data.dt[,
                .(id.col = list(id.col)),
                  by = group.vars,
                  env = list(id.col = id.var)]
    setorderv(groups, group.vars)
    if(add.label) {
      labels <- groups[, ..group.vars]
    }
  }
  if(!is.null(group.labels)) {
    for(i in 1:length(group.vars)) {
      old.labels <- as.character(groups[[group.vars[i]]])
      new.labels <- factor(group.labels[[i]][old.labels], 
                           levels = group.labels[[i]])
      groups[, (group.vars[i]) := new.labels]
    }
    setorderv(groups, group.vars)
  }
  if(add.label == TRUE & !is.null(group.vars)) {
    groups <- .add_group_label(data = groups,
                               cols = group.vars,
                               expand.label = expand.label,
                               ...)
  }
  if(is.character(add.label) & !is.null(group.vars)) {
    groups <- .add_group_label(data = groups,
                               cols = group.vars,
                               label.name = add.label,
                               expand.label = expand.label,
                               ...)
  }
  if(!is.null(label.prefix)) {
    groups[[add.label]] <- paste0(label.prefix, groups[[add.label]])
  }
  return(groups)
}


.add_group_label <- function(data,
                             cols,
                             label.name = "group.label",
                             expand.label = TRUE,
                             col.sep = ":",
                             cat.sep = "."
                             ) {
    data <- copy(data)
    labels <- data[, ..cols]
    if(expand.label) {
      for(i in 1:length(cols)) {
        labels[[cols[i]]] <- paste(cols[i], labels[[cols[i]]], sep = cat.sep)
      }
    }
    labels.c <- do.call(paste, c(labels, sep = col.sep))
    data[, label.col := labels.c, env = list(label.col = label.name)]
    setorderv(data, c(cols, label.name))
    setcolorder(data, c(label.name, cols))
    return(data)
}


.nb_sequential <- function(dist,
                           n,
                           n.min = 1,
                           deg.max = NULL) {

  if(is.null(deg.max)) {
    deg.max <- length(n)
  }

  unit.n <- data.table(id.to = 1:length(n),
                       n = n)

  dist <-
    copy(dist) |>
    merge(unit.n, by = "id.to")
  setkey(dist, id.from)

  dist.ord <-
    dist[order(id.from, dist, n)
         ][n > 0,
           .(id.to,
             n,
             deg = 1:length(dist),
             nc = cumsum(n),
             dist),
           by = id.from
           ]

  nb.dt <-
    merge(dist.ord,
          dist.ord[nc >= n.min,
                   .(n.min.u = min(nc)),
                   by = id.from]) |>
    DT(nc <= n.min.u,
       .(neighbourhood = list(id.to),
         n = unique(n.min.u),
         degree = max(deg)),
       by = id.from)

  return(as.list(nb.dt[, -"id.from"]))
}


.nb_expand <- function(A,
                       n,
                       n.min = 1,
                       deg.max = NULL) {

  if(is.null(deg.max)) {
    deg.max <- length(n)
  }

  nbh <- list() 

  n.nbh <- n
  deg.nbh <- integer(length(n))

  # deg 0 nbh

  bmu.sel <- which(n.nbh >= n.min)
  nbh[bmu.sel] <- bmu.sel

  deg <- 0
  while(any(n.nbh < n.min)) {
    deg <- deg + 1
    # deg 1 nbh
    bmu.sel <- which(n.nbh < n.min)
    if(deg < 2) {
      A.nbh <- A
    } else {
      A.nbh <- A.nbh %*% A
    }

    nbh.sel <- apply(A.nbh[bmu.sel,,drop = FALSE], 1,
                    \(x) which(x > 0),
                    simplify = FALSE)
    n.nbh[bmu.sel] <- unlist(lapply(nbh.sel, \(x) sum(n[x])))
    deg.nbh[bmu.sel] <- deg
    nbh[bmu.sel] <- nbh.sel
    if(deg == deg.max) {
      nbh[which(n.nbh < n.min)] <- NA
      break
    }
  }

  neighbourhood <-
    list(neighbourhood = nbh,
         n = n.nbh,
         degree = deg.nbh)

  return(neighbourhood)
}


egp_define_counterfactual <-
  function(data,
           som,
           cf.ids = NULL,
           fac.ids = NULL,
           compare.by = NULL,
           group.by = NULL,
           som.var = "som.unit",
           id.var = "id",
           group.name = "group.id",
           unit.name = "cf.unit",
           assign.name = "assigned",
           assign.cat = c("counterfactual", "factual"),
           n.min = 1,
           nb.strategy = "sequential",
           deg.max = NULL) {

  data.dt <- as.data.table(data)

  if(!id.var %in% names(data.dt)) {
    data.dt[, id.col := 1:.N, env = list(id.col = id.var)]
  }

  if(is.list(group.by)) {
    group.by.c <- unique(do.call(c, group.by))
  }
  if(!is.list(group.by)) {
    group.by.c <- group.by
    group.by <- list(group.by)
  }

  vars.sel <- c(id.var, compare.by, group.by.c, som.var)

  data.dt <- data.dt[, ..vars.sel]

  if(is.null(fac.ids) & is.null(cf.ids)) {
    stop("Either `fac.ids` or `cf.ids` must be provided")
  }
  if(!is.null(fac.ids) & is.null(cf.ids)) {
    cf.ids <-
      data.dt[!id.col %in% fac.ids,
              id.col,
              env = list(id.col = id.var)]
  }
  if(is.null(fac.ids) & !is.null(cf.ids)) {
    fac.ids <-
      data.dt[!id.col %in% cf.ids,
              id.col,
              env = list(id.col = id.var)]
  }

  # Reference SOM units

  data.dt <- data.dt[.(id.col = unique(c(fac.ids, cf.ids))),
                     on = id.var,
                     env = list(id.col = id.var)]
  setkeyv(data.dt, id.var)

  bmu.ref <- data.dt[.(cf.ids), on = id.var]

  if(is.null(compare.by)) {

    n.bmu <- integer(nrow(som$grid.nb))
    n.ref <- bmu.ref[order(som.col),
                     .(n = .N),
                     by = som.col,
                     env = list(som.col = som.var)]
    n.bmu[n.ref[, som.col, env = list(som.col = som.var)]] <- n.ref$n
    n.bmu[is.na(n.bmu)] <- 0


    if(nb.strategy == "sequential") {
      nbh <-
        .nb_sequential(dist = som$unit.dist,
                       n = n.bmu,
                       n.min = n.min,
                       deg.max = deg.max)
    }
    if(nb.strategy == "expand") {
      nbh <-
        .nb_expand(A = som$grid.nb,
                   n = n.bmu,
                   n.min = n.min,
                   deg.max = deg.max)
    }

    cf.bmu.dt <- 
      data.dt[.(fac.ids),
              on = id.var
              ][,
                .(id.col,
                  cf.col = nbh$neighbourhood[som.col]),
                env = list(id.col = id.var,
                           cf.col = unit.name,
                           som.col = som.var)]

    cf.ids.dt <-
      data.dt[.(cf.ids),
              on = id.var
              ][order(som.col, id.col),
                .(id.col = list(id.col)),
                by = (cf.col = som.col),
                env = list(id.col = id.var,
                           som.col = som.var,
                           cf.col = unit.name)]

  } else {


    bmus <- 1:nrow(som$grid.nb)

    bmu.by <-
      data.dt[,
              .(som.col = 1:nrow(som$grid.nb)),
              by = compare.by,
              env = list(som.col = som.var)]

    comp <-
      data.dt[, .(.n = .N), by = c(compare.by)]
    comp[, .cfgrp := 1:.N] 

    comp.bmu.cf <-
      data.dt[.(cf.ids), on = id.var
              ][, .(.n = .N), by = c(som.var, compare.by)] |>
      merge(bmu.by, by = c(som.var, compare.by), all = TRUE)
 
    nbh.l <- list()
    cf.ids.l <- list()
    for(i in 1:nrow(comp)) {
      
      n.ref <- comp.bmu.cf[comp[i, compare.by, with = FALSE],
                           on = compare.by
                           ][order(som.col),
                             .(som.col, .n),
                             env = list(som.col = som.var)]

      n.bmu <- integer(nrow(som$grid.nb))
      n.bmu[n.ref[, som.col, env = list(som.col = som.var)]] <- n.ref$.n
      n.bmu[is.na(n.bmu)] <- 0

      if(nb.strategy == "expand") {
        nbh.l[[i]] <-
          .nb_expand(A = som$grid.nb,
                     n = n.bmu,
                     n.min = n.min,
                     deg.max = deg.max)$neighbourhood
      }
      if(nb.strategy == "expand") {
        nbh.l[[i]] <-
          .nb_sequential(dist = som$unit.dist,
                         n = n.bmu,
                         n.min = n.min,
                         deg.max = deg.max)$neighbourhood
      }
      
      cf.ids.l[[i]] <-
        data.dt[.(cf.ids),
                on = id.var
                ][comp[i, compare.by, with = FALSE],
                  on = compare.by
                  ][order(som.col, id.col),
                    .(id.col = list(id.col)),
                    by = .(cf.col = som.col),
                    env = list(id.col = id.var,
                               som.col = som.var,
                               cf.col = unit.name)]

    }

    comp[, `:=`(.nbh = nbh.l)]

    comp.cf <-
      comp[,
           .(cf.col = unlist(.nbh, recursive = FALSE),
             som.col = bmus),
           by = c(".cfgrp", compare.by),
           env = list(som.col = som.var,
                      cf.col = unit.name)
           ]

    cf.bmu.dt <-
      cbind(
            data.dt[.(fac.ids),
                    on = id.var
                    ][order(id.col),
                      env = list(id.col = id.var)
                      ][,
                        c(id.var, compare.by), with = FALSE],
            comp.cf[data.dt[.(fac.ids),
                            on = id.var
                            ][order(id.col),
                              env = list(id.col = id.var)
                              ][,
                                c(compare.by, som.var),
                                with = FALSE],
                                .(cf.col),
                                on = c(compare.by, som.var),
                                env = list(cf.col = unit.name)])
 
    cf.ids.dt <- 
      rbindlist(cf.ids.l, idcol = ".cfgrp") |>
      merge(comp[, -c(".n", ".nbh")], by = ".cfgrp") |>
      setcolorder(c(compare.by, unit.name, id.var)) |>
      DT(, -".cfgrp")

  }

  data.cf <- data.dt[.(cf.ids), on = id.var ]
  data.cf[,
          assign.col := assign.cat[1],
          env = list(assign.col = assign.name)]

  data.fac <- data.dt[.(fac.ids), on = id.var]
  data.fac[,
           assign.col := assign.cat[2],
           env = list(assign.col = assign.name)]
  
  data.ret <- rbind(data.fac, data.cf)

  # Groups

  groups.l <- list()
  for(i in seq_along(group.by)){
    groups.l[[i]] <-
      .ids_by_group(data.dt[.(fac.ids)],
                    id.var = id.var,
                    group.vars = group.by[[i]],
                    add.label = FALSE)
  }

  groups <- rbindlist(groups.l, fill = TRUE)

  groups[,
         group.col := 1:nrow(groups),
         env = list(group.col = group.name)]
  setcolorder(groups, c(group.name, group.by.c, id.var)) 

  counterfactual <- list(data = data.ret,
                         assignment = cf.bmu.dt,
                         groups = groups,
                         units = cf.ids.dt,
                         compare.by = compare.by,
                         group.by = group.by,
                         group.by.c = group.by.c,
                         id.var = id.var,
                         som.var = som.var,
                         group.var = group.name,
                         unit.var = unit.name,
                         assign.var = assign.name,
                         assign.cat = assign.cat)

  class(counterfactual) <- c(class(counterfactual), "egp_cf_def")

  return(counterfactual)
}


get_influence <- function(x) {
  counts <-
    merge(
          x$units[,
                  .(.w = 1/unlist(lapply(id.col, length)),
                  id.col = unlist(id.col)),
                  by = c(x$unit.var, x$compare.by),
                  env = list(id.col = x$id.var)] ,
          x$assignment[,
                       .(cf.col = unlist(cf.col)),
                       by = eval(x$compare.by),
                       env = list(cf.col = x$unit.var)
                       ][,
                         .(.n = .N),
                         by = c(x$unit.var, x$compare.by)],
          all.x = TRUE)
  counts[is.na(.n), .n := 0]
  counts[, .influence := .n * .w]
  influence <-
    counts[,
           c(x$compare.by, x$unit.var, x$id.var, ".influence"),
           with = FALSE] |>
      setorderv(c(x$compare.by, ".influence", x$id.var),
                order = c(rep(1, length(x$compare.by)), -1, 1))
  return(influence)
}


chunk_seq <- function(from, to, size = to) {
  chunk.from <- seq(from, to, size)
  if(length(chunk.from) > 1) {
    chunk.to <- c(chunk.from[2:length(chunk.from)]-1, to)
  } else {
    chunk.to <- to
  }
  chunk.size <- chunk.to - c(0, chunk.to[-length(chunk.to)])
  return(list(from = chunk.from,
              to = chunk.to,
              size = chunk.size))
}


egp_posterior_draw <- function(model,
                               n = 1000,
                               unconditional = TRUE,
                               package = "mgcv",
                               parallel = 1) {
  if(package == "mgcv") {
    post <-
        mgcv::rmvn(n = n,
                   mu = coef(model),
                   V = vcov(model, unconditional = unconditional))
  }
  if(package == "mvnfast") {
    post <-
      mvnfast::rmvn(n = n,
                 mu = coef(model),
                 sigma = vcov(model, unconditional = unconditional),
                 ncores = parallel)
  }
  colnames(post) <- names(coef(model))
  post <- Matrix(post)
  return(post)
}


.evaluate_posterior <- 
  function(
           model, 
           posterior,
           data,
           id.var = "id",
           draw.name = ".draw",
           pred.name = "predicted",
           type = "link",
           obs = NULL,
           coef = NULL,
           marginals = NULL,
           marginal.ids = NULL,
           predict.chunk = NULL,
           post.chunk = NULL,
           progress = TRUE
           ) {
  posterior <- Matrix(posterior)
  data.dt <- as.data.table(copy(data))
  if(!id.var %in% names(data.dt)) {
    data.dt[, id.col := 1:.N, env = list(id.col = id.var)]
  }
  if(!is.null(obs)) {
    no.obs <-
      data.dt[id.col %in% obs, length(id.col) < 1, env = list(id.col = id.var)]
    if(no.obs) return(NULL)
    data.dt <- copy(data.dt[.(obs), on = id.var])
  }
  setkeyv(data.dt, id.var)
  data.dt
  n <- nrow(data.dt)
  m <- nrow(posterior)
  if(is.null(predict.chunk)) predict.chunk <- n
  if(is.null(post.chunk)) post.chunk <- m
  predict.chunks <- chunk_seq(1, n, predict.chunk)
  post.chunks <- chunk_seq(1, m, post.chunk)
  # Set excluded coefficients to 0
  if(!is.null(coef)) {
    posterior[,-coef] <- 0
  }
  if(is.null(marginals)) {
    marginals <- list(marginal = 1:ncol(posterior))
  }
  if(is.null(marginal.ids)) {
    marginal.ids <- list(marginal = data.dt[[id.var]])
  }
  if(length(marginals) != length(marginal.ids))
    stop("Number of elements in `marginals` and `marginal.ids` must match.")
  id.lu <- cbind(data.dt[[id.var]], 1:nrow(data.dt))
  colnames(id.lu) <- c(id.var, ".row")
  id.lu <- as.data.table(id.lu)
  setorder(id.lu, .row)
  mar.lu <- list()
  for(i in seq_along(marginal.ids)) {
    mar.lu[[i]] <- id.lu[id.col %in% marginal.ids[[i]],
                         env = list(id.col = id.var)]
    mar.lu[[i]][, marginal := names(marginals)[i]]
  }
  mar.lu <- rbindlist(mar.lu)
  mar.lu$chunk <-
    cut(mar.lu$.row,
        with(predict.chunks, c(to[1] - size[1], to)),
        labels = FALSE)
  mar.lu <-
    mar.lu[,
           .(n = .N,
             row.ids = list(mar.lu$.row[.I]),
             ids = list(mar.lu[[id.var]][.I])),
           by = c("marginal", "chunk")]
  mar.lu[, `:=`(eval.id.from = 1 + cumsum(n) - n, eval.id.to = cumsum(n)), marginal]
  evaluated <- list()
  for(i in seq_along(marginals)) {
      # idx <- id.lu[eval(parse(text = paste(id.col, "%in% marginal.ids[[i]]")))]
      idx <- id.lu[id.col %in% marginal.ids[[i]],
                   env = list(id.col = id.var)]
      evaluated[[i]] <- Matrix(numeric(0), nrow = nrow(idx), ncol = m)
      dimnames(evaluated[[i]])[1] <- list(paste0(id.var, ":", idx[[id.var]]))
  }
  if(length(predict.chunks$from) < 2) progress <- FALSE
  if(progress) {
    prog <- txtProgressBar(min = 0, max = length(predict.chunks$from), initial = 0,
                           char = "=", width = NA, title = "Progress", style = 3)
  }
  for(i in 1:length(predict.chunks$from)) {
    Xp <-
      predict(model,
              newdata = data.dt[predict.chunks$from[i]:predict.chunks$to[i],],
              type = "lpmatrix",
              block.size = predict.chunks$size[i],
              newdata.guaranteed = TRUE,
              discrete = FALSE) |>
      Matrix()
    for(j in 1:length(marginals)) {
      chunk.lu <- mar.lu[chunk == i & marginal == names(marginals)[j]]
      if(nrow(chunk.lu) < 1) next
      pc.rows <- unlist(chunk.lu$row.ids) - with(predict.chunks, to[i] - size[i])
      eval.rows <- with(chunk.lu, eval.id.from:eval.id.to)
      m.predict.chunk <- Matrix(numeric(0), nrow = chunk.lu$n, ncol = m)
      for(k in 1:length(post.chunks$from)) {
        lp <-
          Xp[pc.rows, marginals[[j]]] %*%
          t(posterior[post.chunks$from[k]:post.chunks$to[k], marginals[[j]]])
        if(type == "response") {
          fam <- model$family
          m.predict.chunk[, post.chunks$from[k]:post.chunks$to[k]] <-
           fam$linkinv(as.matrix(lp))
        } else {
          m.predict.chunk[, post.chunks$from[k]:post.chunks$to[k]] <- lp
        }
        rm(lp)
      }
      evaluated[[j]][eval.rows,] <- m.predict.chunk
      rm(m.predict.chunk)
    }
    rm(Xp)
    gc()
    if(progress) {
      setTxtProgressBar(prog, i)
    }
  }
  if(progress) close(prog)
  evaluated.dt <- list()
  for(i in seq_along(evaluated)) {
    eval.dt <- as.data.table(as.matrix(t(evaluated[[i]])))
    eval.dt[, draw.col := 1:nrow(eval.dt), env = list(draw.col = draw.name)]
    toid <- list(as.integer)
    names(toid) <- id.var
    evaluated.dt[[i]] <-
      melt(eval.dt,
           measure.vars = measurev(toid, pattern = paste0(id.var, ":(.*)")),
           value.name = pred.name)
    setkeyv(evaluated.dt[[i]], id.var)
    setindexv(evaluated.dt[[i]], draw.name)
    setorderv(evaluated.dt[[i]], c(draw.name, id.var))
  }
  if(length(marginals) == 1) {
    return(evaluated.dt[[1]])
  } else {
    names(evaluated.dt) <- names(marginals)
    return(evaluated.dt)
  }
}


egp_posterior_predict <- function(model,
                                  posterior,
                                  data = NULL,
                                  id.var = "id",
                                  type = "response",
                                  ids = NULL,
                                  pred.name = NULL,
                                  predict.chunk = NULL,
                                  post.chunk = NULL,
                                  progress = TRUE,
                                  ...
                                  ) {
  if(is.null(pred.name)) {
    if(type == "response") {
      pred.name <- all.vars(update(model$formula, . ~ 1))
    } else {
      pred.name <- ".eta"
    }
  }
  predicted <-
    .evaluate_posterior(model = model,
                        posterior = posterior,
                        data = data,
                        id.var = id.var,
                        pred.name = pred.name,
                        type = type,
                        obs = ids,
                        predict.chunk = predict.chunk,
                        post.chunk = post.chunk,
                        progress = progress,
                        ...)
  return(predicted)
}

# predictions = pred
# ids <- cf$groups$cell
# trans.fun = NULL
# agg.fun = E
# pred.var = "response.int"
# draw.var = ".draw"
# id.var = "cell"
# # ids = NULL
# draw.ids = NULL
# draw.chunk = NULL
# agg.size = 3500
# clamp = NULL
# parallel = NULL
# progress = TRUE

.aggregate_variables <- function(predictions, ...) {
  UseMethod(".aggregate_variables", predictions)
}



.aggregate_variables.data.table <- 
  function(predictions,
           agg.fun = mean,
           trans.fun = NULL,
           ids = NULL,
           pred.var = "predicted",
           draw.var = ".draw",
           id.var = "id",
           agg.name = "aggregated",
           group.name = "group.id",
           draw.ids = NULL,
           draw.chunk = NULL,
           agg.size = NULL,
           parallel = NULL,
           progress = TRUE,
           ...) {
  predictions.dt <- as.data.table(predictions)
  setkeyv(predictions.dt, id.var)
  setindexv(predictions.dt, draw.var)
  if(is.numeric(parallel)) {
    dt.threads.old <- getDTthreads()
    arrow.threads.old <- cpu_count()
    setDTthreads(parallel)
    set_cpu_count(parallel)
  }
  if(is.null(draw.ids)) {
    draw.ids <- predictions.dt[, unique(draw.col), env = list(draw.col = draw.var)]
  }
  if(is.null(draw.chunk)) {
    draw.chunk <- length(draw.ids)
  }
  if(is.null(ids)) {
    ids <- list(predictions.dt[, unique(id.col), env = list(id.col = id.var)])
  }
  draw.chunks <- chunk_seq(1, length(draw.ids), draw.chunk)
  if(is.null(agg.size))
    agg.size <- min(unlist(lapply(ids, length)))
  ids.dt <-
    list(1:length(ids), ids) |>
    as.data.table() |>
    setnames(c(group.name, id.var))
  setorderv(ids.dt, group.name)
  ids.dt[,N := as.numeric(unname(unlist(lapply(ids, length))))]
  ids.dt[order(-N), Nc := cumsum(N)]
  ids.dt[,
         `:=`(aggregate = ifelse(N < agg.size, TRUE, FALSE))
         ][aggregate == TRUE, 
           agg.id := ceiling(Nc / agg.size)]
  setkeyv(ids.dt, group.name)
  agg.ids <- na.omit(unique(ids.dt[order(agg.id), agg.id]))
  single.ids <- ids.dt[aggregate == FALSE,
                       group.col,
                       env = list(group.col = group.name)]
  if(all(ids.dt$N == 0)) return(NULL)
  ids.dt.l <-
    ids.dt[order(-N),
           .(id.col = unlist(id.col)),
           by = c(group.name, "agg.id"),
           env = list(id.col = id.var)]
  setkeyv(ids.dt.l, id.var)
  setindexv(ids.dt.l, list(group.name, "agg.id"))
  if(!is.null(trans.fun)) {
    predictions.dt[,
                   pred.col := trans(pred.col),
                   env = list(pred.col = pred.var,
                              trans = trans.fun)]
  }
  if(progress) {
    prog <- txtProgressBar(min = 0, max = max(ids.dt$Nc) * length(draw.chunks$from), initial = 0,
                           char = "=", width = NA, title = "Progress", style = 3)
    prog.counter <- 0
  }
  draws.agg.l <- list()
  for(i in seq_along(draw.chunks$from)) {
    draws.sum <- draw.ids[draw.chunks$from[i]:draw.chunks$to[i]]
    predictions.draws <- predictions.dt[.(draws.sum),
                                        on = draw.var,
                                        nomatch = NULL]
    setkeyv(predictions.draws, id.var)
    single.l <- list()
    for(j in seq_along(single.ids)) {
      ids.sum <- ids.dt.l[.(single.ids[j]), on = group.name]
      single.l[[j]] <-
        predictions.draws[.(ids.sum),
                          nomatch = NULL,
                          on = id.var
                          ][order(draw.col, group.col),
                            .(agg.col = agg(pred.col)),
                            by = c(draw.var, group.name),
                            env = list(group.col = group.name,
                                       draw.col = draw.var,
                                       pred.col = pred.var,
                                       agg.col = agg.name,
                                       agg = agg.fun)]
      if(progress) {
        prog.counter <- prog.counter + ids.dt[group.col == single.ids[j],
                                              N,
                                              env = list(group.col = group.name)]
        setTxtProgressBar(prog, prog.counter)
      }
    }
    agg.l <- list()
    for(k in seq_along(agg.ids)) {
      match.ids.groups <- ids.dt.l[.(agg.ids[k]), on = "agg.id"]
      agg.l[[k]] <-
        predictions.draws[match.ids.groups,
                          nomatch = NULL,
                          on = id.var,
                          allow.cartesian = TRUE
                          ][order(draw.col, group.col),
                            .(agg.col = agg(pred.col)),
                            by = c(draw.var, group.name),
                            env = list(group.col = group.name,
                                       draw.col = draw.var,
                                       pred.col = pred.var,
                                       agg.col = agg.name,
                                       agg = agg.fun)]
      if(progress) {
        prog.counter <- prog.counter + sum(ids.dt[agg.id == agg.ids[k], N])
        setTxtProgressBar(prog, prog.counter)
      }
    }
    draws.agg.l[[i]] <- rbindlist(c(single.l, agg.l), fill = TRUE)
  }
  draws.agg <- rbindlist(draws.agg.l, fill = TRUE)
  setorderv(draws.agg, c(draw.var, group.name))
  if(progress) close(prog)
  if(is.numeric(parallel)) {
    setDTthreads(dt.threads.old)
  }
  return(draws.agg)
}


egp_summarize_units <- function(predictions,
                                cf.def,
                                pred.var = NULL,
                                draw.chunk = NULL,
                                agg.size = NULL,
                                parallel = NULL,
                                progress = TRUE,
                                ...
                                ){

  predictions.dt <- as.data.table(predictions)

  id.var <- cf.def$id.var
  unit.var <- cf.def$unit.var
  compare.by <- cf.def$compare.by
  units <- copy(cf.def$units)

  if(is.null(pred.var)) {
    pred.sel <- which(!names(predictions.dt) %in% c(id.var, ".draw"))[1]
    pred.var <- names(predictions.dt)[pred.sel]
  }

  unit.sum <-
    .aggregate_variables(predictions.dt,
               agg.fun = mean,
               ids = units[[id.var]],
               pred.var = pred.var,
               draw.var = ".draw",
               id.var = id.var,
               agg.name = pred.var,
               group.name = ".uid",
               draw.chunk = draw.chunk,
               agg.size = agg.size,
               parallel = parallel,
               progress = progress)

  units[, .uid := 1:.N]

  unit.sum <-
    merge(units[, -id.var, with = FALSE], unit.sum, all.x = FALSE) |>
    DT(, -".uid", with = FALSE)
  setindexv(unit.sum, ".draw")
  setorderv(unit.sum, c(compare.by, unit.var, ".draw"))
  return(unit.sum)
}


egp_evaluate_factual <- function(predictions,
                                 cf.def,
                                 name = "factual",
                                 group.eval = NULL,
                                 pred.var = NULL,
                                 draw.chunk = NULL,
                                 agg.size = NULL,
                                 parallel = NULL,
                                 progress = TRUE,
                                 ...
                                 ){

  predictions.dt <- as.data.table(predictions)

  id.var <- cf.def$id.var
  group.var <- cf.def$group.var

  if(is.null(pred.var)) {
    pred.sel <- which(!names(predictions.dt) %in% c(id.var, ".draw"))[1]
    pred.var <- names(predictions.dt)[pred.sel]
  }

  if(is.null(group.eval)) {
    group.eval <- cf.def$groups[[group.var]]
  }

  factual <-
    .aggregate_variables(predictions.dt,
               agg.fun = mean,
               ids = cf.def$groups[.(group.eval),
                                   id.col,
                                   on = group.var,
                                   env = list(id.col = id.var)],
               pred.var = pred.var,
               draw.var = ".draw",
               id.var = id.var,
               agg.name = name,
               group.name = group.var,
               draw.chunk = draw.chunk,
               agg.size = agg.size,
               parallel = parallel,
               progress = progress)

  factual[,
          group.col := group.eval[group.col],
          env = list(group.col = group.var)]
  factual.aug <-
    merge(cf.def$groups[, -id.var, with = FALSE], factual, all.x = FALSE)

  setkeyv(factual.aug, group.var)
  setindexv(factual.aug, ".draw")
  return(factual.aug)
}


egp_evaluate_counterfactual <- function(predictions,
                                        cf.def,
                                        units = NULL,
                                        name = "counterfactual",
                                        group.eval = NULL,
                                        pred.var = NULL,
                                        draw.chunk = NULL,
                                        progress = TRUE,
                                        ...) {
  predictions.dt <- copy(predictions)

  id.var <- cf.def$id.var
  group.var <- cf.def$group.var
  unit.var <- cf.def$unit.var
  compare.by <- cf.def$compare.by

  if(is.null(units)) {
    units <-
      egp_summarize_units(predictions = predictions,
                          cf.def = cf.def,,
                          pred.var = pred.var,
                          draw.chunk = draw.chunk,
                          progress = progress,
                          ...)
  }
   
  draw.ids <- predictions.dt[, unique(.draw)]
  if(is.null(draw.chunk)) {
    draw.chunk <- length(draw.ids)
  }
  draw.chunks <- chunk_seq(1, length(draw.ids), draw.chunk)

  if(is.null(group.eval)) {
    group.eval <- cf.def$groups[[group.var]]
  }

  if(is.null(pred.var)) {
    pred.sel <- which(!names(predictions.dt) %in% c(id.var, ".draw"))[1]
    pred.var <- names(predictions.dt)[pred.sel]
  }

  assigned <- copy(cf.def$assignment)
  groups <- copy(cf.def$groups[.(group.eval), on = group.var])
  unit.def <- copy(cf.def$units)
  unit.sum <- copy(units)
  setkeyv(unit.sum, unit.var)
  setindexv(unit.sum, ".draw")

  assigned.l <-
    assigned[,
             .(unit.col = unlist(unit.col)),
             by = c(id.var, compare.by),
             env = list(id.col = id.var,
                        unit.col = unit.var)]
  setindexv(assigned.l, unit.var)

  unit.def[,
           .n := unlist(lapply(id.col, length)),
           env = list(id.col = id.var)]
  setkeyv(unit.def, unit.var)

  assigned.l <-
    merge(assigned.l,
          unit.def[, c(compare.by, unit.var, ".n"), with = FALSE],
          by = c(compare.by, unit.var))

  assigned.l[,
             .ntotal := sum(.n),
             by = id.var]
  assigned.l <-
    assigned.l[,
               .w := .n/.ntotal
               ][, -c(".n", ".ntotal")]
  setkeyv(assigned.l, id.var)
  setindexv(assigned.l, unit.var)

  groups.n <- unlist(lapply(groups[[id.var]], length))
  draws.n <- predictions[, length(unique(.draw))]
  total.prog <- sum(groups.n * draws.n)

  if(progress) {
    prog <- txtProgressBar(min = 0,
                           max = total.prog,
                           char = "=", width = NA, title = "Progress", style = 3)
    prog.counter <- 0
  }
  
  group.eval.l <- list()
  for(i in 1:nrow(groups)) {
    ass.gr <-
      assigned.l[.(unlist(groups[[id.var]][i])),
                 on = id.var,
                 .(.w = sum(.w)),
                 by = c(compare.by, unit.var)
                 ][, .w := .w/sum(.w)]
    setkeyv(ass.gr, unit.var)
    draws.eval.l <- list()
    for(j in seq_along(draw.chunks$from)) {
      units.sum.draw <-
        unit.sum[.(draw.chunks$from[j]:draw.chunks$to[j]),
                 on = ".draw"]
      setindexv(units.sum.draw, c(compare.by, unit.var))
      draw.chunk <-
        merge(ass.gr,
              units.sum.draw,
              by = c(compare.by, unit.var),
              allow.cartesian = TRUE)
      draws.eval.l[[j]] <-
        draw.chunk[,
                   .(name.col = sum(.w * pred.col)),
                   by = c(".draw"),
                   env = list(pred.col = pred.var,
                              name.col = name)]
      if(progress) {
        prog.counter <- prog.counter + groups.n[i] * draw.chunks$size[j]
        setTxtProgressBar(prog, prog.counter)
      }
    }
    group.eval.l[[i]] <- rbindlist(draws.eval.l)
  } 

  if(progress) {
    close(prog)
  }

  counterfactual <- rbindlist(group.eval.l, idcol = group.var)
  counterfactual[,
                 group.col := group.eval[group.col],
                 env = list(group.col = group.var)]
  counterfactual.aug <-
      merge(cf.def$groups[, -id.var, with = FALSE], counterfactual, all.x = FALSE)
  setkeyv(counterfactual.aug, group.var)
  setindexv(counterfactual.aug, ".draw")


  return(counterfactual.aug)
}


egp_marginal <- function(factual,
                         counterfactual,
                         type = "absolute",
                         marginal.name = "marginal",
                         fac.var = NULL,
                         cf.var = NULL) {

  names.com <- intersect(names(factual),
                         names(counterfactual))
  if(is.null(fac.var)) {
    fac.var <-
      names(factual)[which(!names(factual) %in% names.com)[1]]
  }
  if(is.null(cf.var)) {
    cf.var <-
      names(counterfactual)[which(!names(counterfactual) %in% names.com)[1]]
  }

  marginal <-
    merge(factual, counterfactual, by = names.com)

  if(type == "absolute") {
    marginal[,
             marginal.col := fac.col - cf.col,
             env = list(marginal.col = marginal.name,
                        fac.col = fac.var,
                        cf.col = cf.var)]
  }
  if(type == "relative") {
    marginal[,
             marginal.col := fac.col / cf.col,
             env = list(marginal.col = marginal.name,
                        fac.col = fac.var,
                        cf.col = cf.var)]
  }
  
  return(marginal)
}
