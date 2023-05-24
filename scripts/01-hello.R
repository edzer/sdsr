#| fig.cap: "A first map: birth counts 1974-78, North Carolina counties"
#| code-fold: true
#| fig.height: 4
library(tidyverse)
library(sf)
system.file("gpkg/nc.gpkg", package="sf") |>
	read_sf() -> nc
nc.32119 <- st_transform(nc, 'EPSG:32119')
nc.32119 |>
	select(BIR74) |>
	plot(graticule = TRUE, axes = TRUE)


#| code-fold: true
#| collapse: false
nc |> select(AREA, BIR74, SID74) |> print(n = 3)


#| code-fold: true
year_labels <- c("SID74" = "1974 - 1978", "SID79" = "1979 - 1984")
nc.32119 |> select(SID74, SID79) |>
	pivot_longer(starts_with("SID")) -> nc_longer
ggplot() + geom_sf(data = nc_longer, aes(fill = value), linewidth = 0.4) + 
  facet_wrap(~ name, ncol = 1, labeller = labeller(name = year_labels)) +
  scale_y_continuous(breaks = 34:36) +
  scale_fill_gradientn(colors = sf.colors(20)) +
  theme(panel.grid.major = element_line(color = "white"))


#| code-fold: true
#| fig.cap: "Interactive map created with **mapview**: pan and zoom move the map and change scale; clicking a county pops up window with the available county properties."
library(mapview) |> suppressPackageStartupMessages()
mapviewOptions(fgb = FALSE)
nc.32119 |> mapview(zcol = "BIR74", legend = TRUE, col.regions = sf.colors)


#| fig.cap: "Interactive map created with **mapview**, showing feature attributes for a selected county in a popup window."
knitr::include_graphics("images/mapview.png")


#| code-fold: true
#| fig.cap: "Raster maps (Olinda, Atlantic coast of Brazil): Landsat-7 blue band, with colour values derived from data values (a), the top-left $10 \\times 10$ sub-image from (a) with numeric values shown (b), and overlayed by two different types of vector data: three sample points (c), and a 500 m radius around the points represented as polygons (d)"
#| fig.height: 5
library(stars)
par(mfrow = c(2, 2))
par(mar = rep(1, 4))
tif <- system.file("tif/L7_ETMs.tif", package = "stars")
x <- read_stars(tif)[,,,1]
image(x, main = "(a)")
image(x[,1:10,1:10], text_values = TRUE, border = 'grey', main = "(b)")
image(x, main = "(c)")
set.seed(131)
pts <- st_sample(st_as_sfc(st_bbox(x)), 3)
plot(pts, add = TRUE, pch = 3, col = 'blue')
image(x, main = "(d)")
plot(st_buffer(pts, 500), add = TRUE, pch = 3, border = 'blue', col = NA, lwd = 2)


#| code-fold: true
st_extract(x, pts) # query at points
aggregate(x, st_buffer(pts, 500), FUN = mean) |> st_as_sf() # aggregate over circles


#| fig.cap: "Map obtained by rasterizing county births counts for the period 1974-78 shown in 1.1"
#| code-fold: true
plot(st_rasterize(nc["BIR74"], dx = 0.1), col = sf.colors(), breaks = "equal")


#| code-fold: true
#| fig.cap: "Various raster geometry types"
x <- 1:5
y <- 1:4
d <- st_dimensions(x = x, y = y, .raster = c("x", "y"))
m <- matrix(runif(20),5,4)
r1 <- st_as_stars(r = m, dimensions = d)

r <- attr(d, "raster")
r$affine <- c(0.2, -0.2)
attr(d, "raster") = r
r2 <- st_as_stars(r = m, dimensions = d)

r <- attr(d, "raster")
r$affine <- c(0.1, -0.3)
attr(d, "raster") = r
r3 = st_as_stars(r = m, dimensions = d)

x <- c(1, 2, 3.5, 5, 6)
y <- c(1, 1.5, 3, 3.5)
d <- st_dimensions(x = x, y = y, .raster = c("x", "y"))
r4 <- st_as_stars(r = m, dimensions = d)

grd <- st_make_grid(cellsize = c(10,10), offset = c(-130,10), n = c(8,5), crs = st_crs('OGC:CRS84'))
r5 <- st_transform(grd, "+proj=laea +lon_0=-70 +lat_0=35")

par(mfrow = c(2,3), mar = c(0.1, 1, 1.1, 1))
r1 <- st_make_grid(cellsize = c(1,1), n = c(5,4), offset = c(0,0))
plot(r1, main = "regular")
plot(st_geometry(st_as_sf(r2)), main = "rotated")
plot(st_geometry(st_as_sf(r3)), main = "sheared")
plot(st_geometry(st_as_sf(r4, as_points = FALSE)), main = "rectilinear")
plot(st_geometry((r5)), main = "curvilinear")


#| code-fold: true
#| out.width: '100%'
#| fig.cap: "**sf** and its dependencies; arrows indicate strong dependency, dashed arrows weak dependency"
knitr::include_graphics("images/sf_deps.png")

