
library(sf)
n <- 30
set.seed(13531) # remove this to create another random sequence
xy <- data.frame(x = runif(n), y = runif(n)) |> 
    st_as_sf(coords = c("x", "y"))



w1 <- st_bbox(c(xmin = 0, ymin = 0, xmax = 1, ymax = 1)) |> 
		st_as_sfc() 
w2 <- st_sfc(st_point(c(1, 0.5))) |> st_buffer(1.2)


#| fig.height: 3.3
#| fig.cap: "Depending on the observation window (grey), the same point pattern can appear completely spatially random (left), or clustered (right)"
#| code-fold: true
par(mfrow = c(1, 2), mar = c(2.1, 2.1, 0.1, 0.5), xaxs = "i", yaxs = "i")
plot(w1, axes = TRUE, col = 'grey')
plot(xy, add = TRUE)
plot(w2, axes = TRUE, col = 'grey')
plot(xy, add = TRUE, cex = .5)



library(spatstat) |> suppressPackageStartupMessages()
as.ppp(xy)



(pp1 <- c(w1, st_geometry(xy)) |> as.ppp())
c1 <- st_buffer(st_centroid(w2), 1.2)
(pp2 <- c(c1, st_geometry(xy)) |> as.ppp())


#| fig.height: 3.5
#| fig.cap: "3 $\\times$ 3 quadrat counts for the two point patterns"
par(mfrow = c(1, 2), mar = rep(0, 4))
q1 <- quadratcount(pp1, nx=3, ny=3)
q2 <- quadratcount(pp2, nx=3, ny=3)
plot(q1, main = "")
plot(xy, add = TRUE)
plot(q2, main = "")
plot(xy, add = TRUE)



quadrat.test(pp1, nx=3, ny=3)
quadrat.test(pp2, nx=3, ny=3)



den1 <- density(pp1, sigma = bw.diggle)
den2 <- density(pp2, sigma = bw.diggle)

#| fig.height: 3.0
#| fig.cap: "Kernel densities for both point patterns"
#| code-fold: true
par(mfrow = c(1, 2), mar = c(0,0,1.1,2))
plot(den1)
plot(pp1, add=TRUE)
plot(den2)
plot(pp1, add=TRUE)



library(stars)
s1 <- st_as_stars(den1)
(s2 <- st_as_stars(den2))



s1$a <- st_area(s1) |> suppressMessages()
s2$a <- st_area(s2) |> suppressMessages()
with(s1, sum(v * a, na.rm = TRUE))
with(s2, sum(v * a, na.rm = TRUE))



pt <- st_sfc(st_point(c(0.5, 0.5)))
st_as_sf(s2, as_points = TRUE, na.rm = FALSE) |>
  st_distance(pt) -> s2$dist



(m <- ppm(pp2 ~ dist, data = list(dist = as.im(s2["dist"]))))


#| fig.cap: "Predicted densities of a ppm model"
#| out.width: 50%
#| code-fold: true
plot(m, se = FALSE)



predict(m, covariates = list(dist = as.im(s2["dist"]))) |>
    st_as_stars()



system.file("gpkg/nc.gpkg", package = "sf") |> 
    read_sf() |>
	st_geometry() |>
    st_centroid() |>
    as.ppp()



longleaf
ll <- st_as_sf(longleaf)
print(ll, n = 3)



as.ppp(ll)



print(st_as_sf(copper$SouthLines), n = 5)



print(st_as_sf(chicago), n = 5)



table(st_as_sf(chicago)$label)



kappa <- 30 / st_area(w2) # intensity
th <- st_sample(w2, kappa = kappa, mu = 3, scale = 0.05, 
    type = "Thomas")
nrow(th)

#| fig.height: 4
#| out.width: 60%
#| fig.cap: "Thomas process with mu = 3 and scale = 0.05"
#| code-fold: true
par(mar = rep(0, 4))
plot(w2)
plot(th, add = TRUE)


#| fig.cap: "Points sampled on the globe over the oceans: randomly (left) and approximately regular (Fibonacci; right), shown in an orthographic projection"
#| out.width: 100%
#| code-fold: true
par(mar = rep(0, 4), mfrow = c(1, 2))
library(s2)
g <- as_s2_geography(TRUE) # Earth
co <- s2_data_countries()
oc <- s2_difference(g, s2_union_agg(co)) # oceans
b <- s2_buffer_cells(as_s2_geography("POINT(-30 -10)"), 9800000) # visible half
i <- s2_intersection(b, oc) # visible ocean
co <- s2_intersection(b, co)
ortho = st_crs("+proj=ortho +lat_0=-10 +lon_0=-30")
# background:
st_transform(st_as_sfc(i), ortho) |> plot(col = 'lightblue')
st_transform(st_as_sfc(co), ortho) |> plot(col = NA, add = TRUE, border = 'grey')
# sampling randomly from globe:
sf_use_s2(FALSE) |> suppressMessages()
st_as_stars() |> st_bbox() |> st_as_sfc() |> st_sample(1000, exact = FALSE) |>
  suppressMessages() -> pts
sf_use_s2(TRUE) |> suppressMessages()
pts |> s2_intersection(i) |> st_as_sfc() -> pts
# add:
st_transform(pts, ortho) |> plot(add = TRUE, pch = 3, cex = .5)
# right: background:
st_transform(st_as_sfc(i), ortho) |> plot(col = 'lightblue')
st_transform(st_as_sfc(co), ortho) |> plot(col = NA, add = TRUE, border = 'grey')
# Fibonacci:
sf_use_s2(FALSE) |> suppressMessages()
st_as_stars() |> st_bbox() |> st_as_sfc() |> 
  st_sample(1000, type = "Fibonacci", exact = FALSE) |> suppressMessages() -> pts
sf_use_s2(TRUE) |> suppressMessages()
pts |> s2_intersection(i) |> st_as_sfc() -> pts
st_transform(pts, ortho) |> plot(add = TRUE, pch = 3, cex = .5)



# we need to detach, to avoid name clashes of using idw() and diagnose() later on
library(spatstat)
w <- function(x) which(x == search())
detach(pos = w("package:spatstat"))
detach(pos = w("package:spatstat.linnet"))
pc <- w("package:spatstat.core")
if (length(pc)) {
	detach(pos = w("package:spatstat.core"))
} else { # > 3.0.0
	detach(pos = w("package:spatstat.model"))
	detach(pos = w("package:spatstat.explore"))
}
detach(pos = w("package:spatstat.random"))
detach(pos = w("package:spatstat.geom"))
detach(pos = w("package:spatstat.data"))

