#| out.width: 60%
#| fig.cap: "Two-dimensional polar (red) and Cartesian (blue) coordinates"
#| code-fold: true
par(mar = rep(0,4))
plot(3, 4, xlim = c(-6,6), ylim = c(-6,6), asp = 1)
axis(1, pos = 0, at = 0:6)
axis(2, pos = 0, at = -6:6)
xd <- seq(-5, 5, by = .1)
lines(xd, sqrt(25 - xd^2), col = 'grey')
lines(xd, -sqrt(25 - xd^2), col = 'grey')
arrows(0, 0, 3, 4, col = 'red', length = .15, angle = 20)
text(1.5, 2.7, label = "r", col = 'red')
xd <- seq(3/5, 1, by = .1)
lines(xd, sqrt(1 - xd^2), col = 'red')
text(1.2, 0.5, label = parse(text = "phi"), col = 'red')
lines(c(3,3), c(0,4), lty = 2, col = 'blue')
lines(c(0,3), c(4,4), lty = 2, col = 'blue')
text(3.3, 0.3, label = "x", col = 'blue')
text(0.3, 4.3, label = "y", col = 'blue')


#| fig.cap: "Cartesian geocentric coordinates (left) measure three distances, ellipsoidal coordinates (right) measure two angles, and possibly an ellipsoidal height"
#| code-fold: true
library(sf) |> suppressPackageStartupMessages()
e <- cbind(-90:90,0) # equator
f1 <- rbind(cbind(0, -90:90)) # 0/antimerid.
f2 <- rbind(cbind(90, -90:90), cbind(270, 90:-90))# +/- 90
eq <- st_sfc(st_linestring(e), st_linestring(f1), st_linestring(f2), crs='OGC:CRS84')

geoc <- st_transform(eq, "+proj=geocent")
cc <- rbind(geoc[[1]], NA, geoc[[2]], NA, geoc[[3]])
from3d <- function(x, offset, maxz, minz) {
	x = x[,c(2,3,1)] + offset # move to y right, x up, z backw
	x[,2] = x[,2] - maxz      # shift y to left
	d = maxz
	z = x[,3] - minz + offset
	x[,1] = x[,1] * (d/z)
	x[,2] = x[,2] * (d/z)
	x[,1:2]
}
maxz <- max(cc[,3], na.rm = TRUE)
minz <- min(cc[,3], na.rm = TRUE)
offset <- 3e7
circ <- from3d(cc, offset, maxz, minz)
mx <- max(cc, na.rm = TRUE) * 1.1
x <- rbind(c(0, 0, 0), c(mx, 0, 0))
y <- rbind(c(0, 0, 0), c(0, mx, 0))
z <- rbind(c(0, 0, 0), c(0, 0, mx))
ll <- rbind(x, NA, y, NA, z)
l0 <-  from3d(ll, offset, maxz, minz)
mx <- max(cc, na.rm = TRUE) * 1.2
x <- rbind(c(0, 0, 0), c(mx, 0, 0))
y <- rbind(c(0, 0, 0), c(0, mx, 0))
z <- rbind(c(0, 0, 0), c(0, 0, mx))
ll <- rbind(x, NA, y, NA, z)
l <-  from3d(ll, offset, maxz, minz)

par(mfrow = c(1, 2))
par(mar = rep(0,4))
plot.new()
plot.window(xlim = c(min(circ[,1],na.rm = TRUE), 3607103*1.02), 
						ylim = c(min(circ[,2],na.rm = TRUE), 2873898*1.1), asp = 1)
lines(circ)
lines(l0)
text(l[c(2,5,8),], labels = c("x", "y", "z"), col = 'red')
# add POINT(60 47)
p <- st_as_sfc("POINT(60 47)", crs = 'OGC:CRS84') |> st_transform("+proj=geocent")
p <- p[[1]]
pts <- rbind(c(0,0,0), c(p[1],0,0), c(p[1],p[2],0), c(p[1],p[2],p[2]))
ptsl <- from3d(pts, offset, maxz, minz)
lines(ptsl, col = 'blue', lty = 2, lwd = 2)
points(ptsl[4,1], ptsl[4,2], col = 'blue', cex = 1, pch = 16)

plot.new()
plot.window(xlim = c(min(circ[,1],na.rm = TRUE), 3607103*1.02), 
						ylim = c(min(circ[,2],na.rm = TRUE), 2873898*1.1), asp = 1)
lines(circ)

p <- st_as_sfc("POINT(60 47)", crs = 'OGC:CRS84') |> st_transform("+proj=geocent")
p <- p[[1]]
pts <- rbind(c(0,0,0), c(p[1],p[2],p[3]))
pt <-  from3d(pts, offset, maxz, minz)
lines(pt)
points(pt[2,1], pt[2,2], col = 'blue', cex = 1, pch = 16)

p0 <- st_as_sfc("POINT(60 0)", crs = 'OGC:CRS84') |> st_transform("+proj=geocent")
p0 <- p0[[1]]
pts <- rbind(c(0,0,0), c(p0[1],p0[2],p0[3]))
pt <-  from3d(pts, offset, maxz, minz)
lines(pt)

p0 <- st_as_sfc("POINT(0 0)", crs = 'OGC:CRS84') |> st_transform("+proj=geocent")
p0 <- p0[[1]]
pts <- rbind(c(0,0,0), c(p0[1],p0[2],p0[3]))
pt <-  from3d(pts, offset, maxz, minz)
lines(pt)

p0 <- st_as_sfc("POINT(0 90)", crs = 'OGC:CRS84') |> st_transform("+proj=geocent")
p0 <- p0[[1]]
pts <- rbind(c(0,0,0), c(p0[1],p0[2],p0[3]))
pt <-  from3d(pts, offset, maxz, minz)
lines(pt, lty = 2)

p0 <- st_as_sfc("POINT(90 0)", crs = 'OGC:CRS84') |> st_transform("+proj=geocent")
p0 <- p0[[1]]
pts <- rbind(c(0,0,0), c(p0[1],p0[2],p0[3]))
pt <-  from3d(pts, offset, maxz, minz)
lines(pt, lty = 2)

f1 <- rbind(cbind(0:60, 0))
arc <- st_sfc(st_linestring(f1), crs='OGC:CRS84')
geoc <- st_transform(arc, "+proj=geocent")
cc <- rbind(geoc[[1]])
circ <- from3d(cc, offset, maxz, minz)
lines(circ, col = 'red', lwd = 2, lty = 2)

f1 <- rbind(cbind(60, 0:47))
arc <- st_sfc(st_linestring(f1), crs='OGC:CRS84')
geoc <- st_transform(arc, "+proj=geocent")
cc <- rbind(geoc[[1]])
circ <- from3d(cc, offset, maxz, minz)
lines(circ, col = 'blue', lwd = 2, lty = 2)

text(pt[1,1]+100000, pt[1,2]+50000, labels = expression(phi), col = 'blue') # lat
text(pt[1,1]+20000, pt[1,2]-50000, labels = expression(lambda), col = 'red') # lng


#| code-fold: true
#| collapse: false
p <- st_as_sfc("POINT(60 47)", crs = 'OGC:CRS84')
p[[1]]


#| code-fold: true
#| collapse: false
p <- st_as_sfc("POINT(60 47)", crs = 'OGC:CRS84') |> st_transform("+proj=geocent")
p[[1]]


#| out.width: 60%
#| fig.cap: "Angles on an ellipse: geodetic (blue) and geocentric (red) latitude"
#| code-fold: true
par(mar = rep(0,4))
x <- 4
y <- 5/8 * sqrt(48)
plot(x, y, xlim = c(-6,6), ylim = c(-8,8), asp = 1)
axis(1, pos = 0, at = 0:9)
axis(2, pos = 0, at = -5:5)
xd <- seq(-8, 8, by = .1)
lines(xd, 5/8 * sqrt(64 - xd^2), col = 'grey')
lines(xd, 5/8 * -sqrt(64 - xd^2), col = 'grey')
arrows(0, 0, x, y, col = 'red', length = .15, angle = 20)
b <- (x * 25) / (-y * 64)
a <- y - x * b
abline(a, b, col = 'grey')
b <- -1/b
x0 <- x - y / b
arrows(x0, 0, x, y, col = 'blue', length = .15, angle = 20)
text(1.2, 0.5, label = parse(text = "psi"), col = 'red')
text(3, 0.5, label = parse(text = "phi"), col = 'blue')


#| code-fold: true
#| collapse: false
pts <- st_sfc(st_point(c(13.4050, 52.5200)), st_point(c(2.3522, 48.8566)), crs = 'OGC:CRS84')
s2_orig <- sf_use_s2(FALSE)
d1 <- c(gc_ellipse = st_distance(pts)[1,2])
sf_use_s2(TRUE)
# or, without using s2, use st_distance(st_transform(pts, "+proj=cart +ellps=sphere"))
d2 <- c(gc_sphere = st_distance(pts)[1,2])
p <- st_transform(pts, "+proj=cart +ellps=WGS84")
d3 <- c(str_ellipse = units::set_units(sqrt(sum(apply(do.call(cbind, p), 1, diff)^2)), m))
p2 <- st_transform(pts, "+proj=cart +ellps=sphere")
d4 <- c(str_sphere = units::set_units(sqrt(sum(apply(do.call(cbind, p2), 1, diff)^2)), m))
res <- c(d1, d3, d2, d4) # note order
# print as km, re-add names:
sf_use_s2(s2_orig) # back to what it was before changing
res |> units::set_units(km) |> setNames(names(res)) |> print(digits = 5)


#| out.width: 70%
#| fig.cap: "UK horizontal datum grid, from datum OSGB 1936 (EPSG:4277) to datum ETRS89 (EPSG:4258); units arc-seconds"
#| code-fold: true
library(stars)
library(rnaturalearth)
library(dplyr) |> suppressPackageStartupMessages()
countries110 |> 
    st_as_sf() |>
    filter(ADMIN == "United Kingdom") |>
    st_geometry() -> uk
filename = "data/uk_os_OSTN15_NTv2_OSGBtoETRS.tif"
r <- if (file.exists(filename)) {
    r <- read_stars(filename)
  } else {
    read_stars("/vsicurl/https://cdn.proj.org/uk_os_OSTN15_NTv2_OSGBtoETRS.tif")
  }
hook <- function() {
		plot(uk, border = "orange", col = NA, add = TRUE)
}
plot(r[,,,1:2], axes = TRUE, hook = hook, key.pos = 4)


#| out.width: 70%
#| fig.cap: "UK vertical datum grid, from ETRS89 (EPSG:4937) to ODN height (EPSG:5701), units m"
#| code-fold: true
filename = "data/uk_os_OSGM15_GB.tif"
h <- if (file.exists(filename)) {
    read_stars(filename)
  } else {
    read_stars("/vsicurl/https://cdn.proj.org/uk_os_OSGM15_GB.tif")
  } 
plot(h, axes = TRUE, reset = FALSE)
plot(uk, border = "orange", col = NA, add = TRUE)

