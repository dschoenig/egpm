# Truncated cos^2

cos2_trunc <- function(x, a = 1, b = 1, c = 0) {
  x.t <- pi/(2*b) * x -c
  y <- a * (cos(x.t))^2
  y[x.t < -pi/2 | x.t > pi/2] <- 0
  return(y)
}

plot(seq(-5, 5, length.out = 1000), cos2_trunc(seq(-5, 5, length.out = 1000), b = 2))

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

n=100
x <- rep(seq(-5, 5, length.out = n), each = n)
y <- rep(seq(-5, 5, length.out = n), times = n)

z <- cos(x)^2 * cos(y)^2

z1 <- cos2_2d_trunc(x, y, xd = 3, y0 = 0.5)
z2 <- cos2_2d_trunc(x, y, xd = 5, x0 = 1, y0 = -1.5)
z3 <- cos2_2d_trunc(x, y, xd = 5, x0 = -2.5, y0 = 1)
z <- z1+z2+z3
# z <- cos2_2d_trunc(x, y, xd = 5, yd = 5)
ggplot() +
  geom_raster(aes(x=x, y=y, fill = z)) +
  scale_fill_viridis_c()

y <- cos2_trunc(x, b = 2*10)
plot(x, y, type = "l")


generate_matern(x.dim = 100, y.dim = 100, scale = 10) |>
zplot()
