#| fig.cap: "Sketches of the main simple feature geometry types"
#| code-fold: true
library(sf) |> suppressPackageStartupMessages()
par(mfrow = c(2,4))
par(mar = c(1,1,1.2,1))

# 1
p <- st_point(0:1)
plot(p, pch = 16)
title("point")
box(col = 'grey')

# 2
mp <- st_multipoint(rbind(c(1,1), c(2, 2), c(4, 1), c(2, 3), c(1,4)))
plot(mp, pch = 16)
title("multipoint")
box(col = 'grey')

# 3
ls <- st_linestring(rbind(c(1,1), c(5,5), c(5, 6), c(4, 6), c(3, 4), c(2, 3)))
plot(ls, lwd = 2)
title("linestring")
box(col = 'grey')

# 4
mls <- st_multilinestring(list(
  rbind(c(1,1), c(5,5), c(5, 6), c(4, 6), c(3, 4), c(2, 3)),
  rbind(c(3,0), c(4,1), c(2,1))))
plot(mls, lwd = 2)
title("multilinestring")
box(col = 'grey')

# 5 polygon
po <- st_polygon(list(rbind(c(2,1), c(3,1), c(5,2), c(6,3), c(5,3), c(4,4), c(3,4), c(1,3), c(2,1)),
	rbind(c(2,2), c(3,3), c(4,3), c(4,2), c(2,2))))
plot(po, border = 'black', col = '#ff8888', lwd = 2)
title("polygon")
box(col = 'grey')

# 6 multipolygon
mpo <- st_multipolygon(list(
	list(rbind(c(2,1), c(3,1), c(5,2), c(6,3), c(5,3), c(4,4), c(3,4), c(1,3), c(2,1)),
		rbind(c(2,2), c(3,3), c(4,3), c(4,2), c(2,2))),
	list(rbind(c(3,7), c(4,7), c(5,8), c(3,9), c(2,8), c(3,7)))))
plot(mpo, border = 'black', col = '#ff8888', lwd = 2)
title("multipolygon")
box(col = 'grey')

# 7 geometrycollection
gc <- st_geometrycollection(list(po, ls + c(0,5), st_point(c(2,5)), st_point(c(5,4))))
plot(gc, border = 'black', col = '#ff6666', pch = 16, lwd = 2)
title("geometrycollection")
box(col = 'grey')


#| code-fold: true
## p
## mp
## ls
## mls
## po
## mpo
## gc


#| code-fold: true
#| collapse: false
(ls <- st_linestring(rbind(c(0,0), c(1,1), c(2,2), c(0,2), c(1,1), c(2,0))))
c(is_simple = st_is_simple(ls))


#| code-fold: true
#| collapse: false
st_point(c(1,3,2))
st_point(c(1,3,2), dim = "XYM")
st_linestring(rbind(c(3,1,2,4), c(4,4,2,2)))


#| code-fold: true
(e <- st_intersection(st_point(c(0,0)), st_point(c(1,1))))


#| code-fold: true
#| collapse: false
st_point()
st_linestring(matrix(1,1,3)[0,], dim = "XYM")


#| code-fold: true
library(sf)
polygon <- po <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
p0 <- st_polygon(list(rbind(c(-1,-1), c(2,-1), c(2,2), c(-1,2), c(-1,-1))))
line <- li <- st_linestring(rbind(c(.5, -.5), c(.5, 0.5)))
s <- st_sfc(po, li)

par(mfrow = c(3,3))
par(mar = c(1,1,1,1))

# "1020F1102"
# 1: 1
plot(s, col = c(NA, 'darkgreen'), border = 'blue', main = expression(paste("I(pol)",intersect(),"I(line) = 1")))
lines(rbind(c(.5,0), c(.5,.495)), col = 'red', lwd = 2)
points(0.5, 0.5, pch = 1)

# 2: 0
plot(s, col = c(NA, 'darkgreen'), border = 'blue', main = expression(paste("I(pol)",intersect(),"B(line) = 0")))
points(0.5, 0.5, col = 'red', pch = 16)

# 3: 2
plot(s, col = c(NA, 'darkgreen'), border = 'blue', main = expression(paste("I(pol)",intersect(),"E(line) = 2")))
plot(po, col = '#ff8888', add = TRUE)
plot(s, col = c(NA, 'darkgreen'), border = 'blue', add = TRUE)

# 4: 0
plot(s, col = c(NA, 'darkgreen'), border = 'blue', main = expression(paste("B(pol)",intersect(),"I(line) = 0")))
points(.5, 0, col = 'red', pch = 16)

# 5: F
plot(s, col = c(NA, 'darkgreen'), border = 'blue', main = expression(paste("B(pol)",intersect(),"B(line) = F")))

# 6: 1
plot(s, col = c(NA, 'darkgreen'), border = 'blue', main = expression(paste("B(pol)",intersect(),"E(line) = 1")))
plot(po, border = 'red', col = NA, add = TRUE, lwd = 2)

# 7: 1
plot(s, col = c(NA, 'darkgreen'), border = 'blue', main = expression(paste("E(pol)",intersect(),"I(line) = 1")))
lines(rbind(c(.5, -.5), c(.5, 0)), col = 'red', lwd = 2)

# 8: 0
plot(s, col = c(NA, 'darkgreen'), border = 'blue', main = expression(paste("E(pol)",intersect(),"B(line) = 0")))
points(.5, -.5, col = 'red', pch = 16)

# 9: 2
plot(s, col = c(NA, 'darkgreen'), border = 'blue', main = expression(paste("E(pol)",intersect(),"E(line) = 2")))
plot(p0 / po, col = '#ff8888', add = TRUE)
plot(s, col = c(NA, 'darkgreen'), border = 'blue', add = TRUE)


#| code-fold: true
#| collapse: false
st_relate(polygon, line)


#| fig.cap: "For a set of points, left: convex hull (red); middle: Voronoi polygons; right: Delauney triangulation"
#| code-fold: true
#| out.width: 60%
par(mar = rep(0,4), mfrow = c(1, 3))
set.seed(133331)
mp <- st_multipoint(matrix(runif(20), 10))
plot(mp, cex = 2)
plot(st_convex_hull(mp), add = TRUE, col = NA, border = 'red')
box()
plot(mp, cex = 2)
plot(st_voronoi(mp), add = TRUE, col = NA, border = 'red')
box()
plot(mp, cex = 2)
plot(st_triangulate(mp), add = TRUE, col = NA, border = 'darkgreen')
box()


#| code-fold: true
#| out.width: 50%
#| fig.cap: "Left: three overlapping squares -- how do we identify the small box where all three overlap? Right: unique, non-overlapping n-ary intersections"
par(mar = rep(.1, 4), mfrow = c(1, 2))
sq <- function(pt, sz = 1) st_polygon(list(rbind(c(pt - sz), 
  c(pt[1] + sz, pt[2] - sz), c(pt + sz), c(pt[1] - sz, pt[2] + sz), c(pt - sz))))
x <- st_sf(box = 1:3, st_sfc(sq(c(0, 0)), sq(c(1.7, -0.5)), sq(c(0.5, 1))))
plot(st_geometry(x), col = NA, border = sf.colors(3, categorical = TRUE), lwd = 3)
plot(st_intersection(st_geometry(x)), col = sf.colors(7, categorical=TRUE, alpha = .5))


#| code-fold: true
#| out.width: 50%
#| fig.cap: "Difference between subsequent boxes, left: in original order; right: in reverse order"
par(mar = rep(.1, 4), mfrow = c(1, 2)) 
xg <- st_geometry(x)
plot(st_difference(xg), col = sf.colors(3, alpha = .5, categorical=TRUE))
plot(st_difference(xg[3:1]), col = sf.colors(3, alpha = .5, categorical=TRUE))


#| out.width: 50%
#| fig.cap: "Rasterization artifact: as only top-left corners are part of the raster cell, only cells touching the red line below the diagonal line are rasterized"
#| code-fold: true
library(stars) |> suppressPackageStartupMessages()
par(mar = rep(1, 4))
ls <- st_sf(a = 2, st_sfc(st_linestring(rbind(c(0.1, 0), c(1, .9)))))
grd <- st_as_stars(st_bbox(ls), nx = 10, ny = 10, xlim = c(0, 1.0), ylim = c(0, 1),
   values = -1)
r <- st_rasterize(ls, grd, options = "ALL_TOUCHED=TRUE")
r[r == -1] <- NA
plot(st_geometry(st_as_sf(grd)), border = 'orange', col = NA, 
	 reset = FALSE, key.pos = NULL)
plot(r, axes = FALSE, add = TRUE, breaks = "equal", main = NA) # ALL_TOUCHED=FALSE;
plot(ls, add = TRUE, col = "red", lwd = 2)

