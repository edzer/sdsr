
library(sf)
p1 <- st_point(c(7.35, 52.42))
p2 <- st_point(c(7.22, 52.18))
p3 <- st_point(c(7.44, 52.19))
sfc <- st_sfc(list(p1, p2, p3), crs = 'OGC:CRS84')
st_sf(elev = c(33.2, 52.1, 81.2), 
	  marker = c("Id01", "Id02", "Id03"), geom = sfc)


#| fig-cap: components of an `sf` object
#| out.width: 100%
knitr::include_graphics("images/sf_obj.png")



library(sf)
(file <- system.file("gpkg/nc.gpkg", package = "sf"))
nc <- st_read(file)



st_layers(file)



(file = tempfile(fileext = ".gpkg"))
st_write(nc, file, layer = "layer_nc")



## nc[2:5, 3:7]



nc5 <- nc[1:5, ]
nc7 <- nc[1:7, ]
(i <- st_intersects(nc5, nc7))


#| fig.cap: "First seven North Carolina counties"
#| code-fold: true
#| fig.height: 2.5
plot(st_geometry(nc7))
plot(st_geometry(nc5), add = TRUE, border = "brown")
cc = st_coordinates(st_centroid(st_geometry(nc7)))
text(cc, labels = 1:nrow(nc7), col = "blue")



as.matrix(i)



lengths(i)



lengths(t(i))



methods(class = "sgbp")



library(tidyverse) |> suppressPackageStartupMessages()
nc |> as_tibble() |> select(BIR74) |> head(3)



orange <- nc |> dplyr::filter(NAME == "Orange")
wd <- st_is_within_distance(nc, orange, 
							units::set_units(50, km))
o50 <- nc |> dplyr::filter(lengths(wd) > 0)
nrow(o50)


#| fig.cap: "Orange County (orange), counties within a 50 km radius (black), a 50~km buffer around Orange County (brown), and remaining counties (grey)"
#| code-fold: true
og <- st_geometry(orange)
buf50 <- st_buffer(og, units::set_units(50, km))
all <- c(buf50, st_geometry(o50))
plot(st_geometry(o50), lwd = 2, extent = all)
plot(og, col = 'orange', add = TRUE)
plot(buf50, add = TRUE, col = NA, border = 'brown')
plot(st_geometry(nc), add = TRUE, border = 'grey')


#| fig-cap: "Example of `st_join` with `largest = TRUE`: the label of the polygon in the top figure with the largest intersection with polygons in the bottom figure is assigned to the polygons of the bottom figure."
#| code-fold: true
# example of largest = TRUE:
system.file("shape/nc.shp", package="sf") |> 
    read_sf() |>
    st_transform('EPSG:2264') -> nc
gr <- st_sf(
         label = apply(expand.grid(1:10, LETTERS[10:1])[,2:1], 1, paste0, collapse = ""),
         geom = st_make_grid(nc))
gr$col <- sf.colors(10, categorical = TRUE, alpha = .3)
# cut, to verify that NA's work out:
gr <- gr[-(1:30),]
suppressWarnings(nc_j <- st_join(nc, gr, largest = TRUE))
par(mfrow = c(2,1), mar = rep(0,4))
plot(st_geometry(nc_j), border = 'grey')
plot(st_geometry(gr), add = TRUE, col = gr$col)
text(st_coordinates(st_centroid(st_geometry(gr))), labels = gr$label, cex = .85)
# the joined dataset:
plot(st_geometry(nc_j), border = 'grey', col = nc_j$col)
text(st_coordinates(st_centroid(st_geometry(nc_j))), labels = nc_j$label, cex = .7)
plot(st_geometry(gr), border = '#88ff88aa', add = TRUE)



"POINT(50 50.1)" |> st_as_sfc(crs = "OGC:CRS84") -> pt



"POLYGON((40 40, 60 40, 60 50, 40 50, 40 40))" |>
  st_as_sfc(crs = "OGC:CRS84") -> pol
st_intersects(pt, pol)


#| fig.cap: "Intersection depends on whether we use geodesics/great circle arcs (left: s2) or Cartesian coordinates (right)"
#| code-fold: true
par(mfrow = c(1, 2))
par(mar = c(2.1, 2.1, 1.2, .5))
ortho <- st_crs("+proj=ortho +lon_0=50 +lat_0=45")
pol |> st_transform(ortho) |> plot(axes = TRUE, graticule = TRUE, 
								   main = 's2geometry')
pt |> st_transform(ortho) |> plot(add = TRUE, pch = 16, col = 'red')
# second plot:
plot(pol, axes = TRUE, graticule = TRUE, main = 'GEOS')
plot(pt, add = TRUE, pch = 16, col = 'red')



old <- sf_use_s2(FALSE)
st_intersects(pol, pt)
sf_use_s2(old) # restore



tif <- system.file("tif/L7_ETMs.tif", package = "stars")
library(stars)
(r <- read_stars(tif))



length(r)
class(r[[1]])
dim(r[[1]])



## st_dimensions(r)



st_bbox(r)



tf <- tempfile(fileext = ".tif")
write_stars(r, tf)



## st_drivers("raster")



r[,1:100, seq(1, 250, 5), 4] |> dim()
r[,1:100, seq(1, 250, 5), 4, drop = TRUE] |> dim()



library(dplyr, warn.conflicts = FALSE)
filter(r, x > 289000, x < 290000)



slice(r, band, 3)



b <- st_bbox(r) |>
	st_as_sfc() |>
	st_centroid() |>
	st_buffer(units::set_units(500, m))
r[b]


#| fig.cap: "Circular centre region of the Landsat 7 scene (band 1)"
#| code-fold: true
#| out.width: 70%
plot(r[b][,,,1], reset = FALSE)
plot(b, border = 'brown', lwd = 2, col = NA, add = TRUE)



r[b] |> st_normalize() |> st_dimensions()



r[b, crop = FALSE]



## st_crop(r, b)



aperm(r, c(3, 1, 2))



(rs <- split(r))
merge(rs, name = "band") |> setNames("L7_ETMs")



st_redimension(r, c(x = 349, y = 352, b1 = 3, b2 = 2))



c(r, r, along = "new_dim")



set.seed(115517)
pts <- st_bbox(r) |> st_as_sfc() |> st_sample(20)
(e <- st_extract(r, pts))



set.seed(131)


circles <- st_sample(st_as_sfc(st_bbox(r)), 3) |>
    st_buffer(500)
aggregate(r, circles, max)


#| fig.cap: "Randomly chosen sample locations for training data; red: water, yellow: land"
#| out.width: 80%
#| code-fold: true
plot(r[,,,1], reset = FALSE)
col <- rep("yellow", 20)
col[c(8, 14, 15, 18, 19)] = "red"
st_as_sf(e) |> st_coordinates() |> text(labels = 1:20, col = col)



rs <- split(r)
trn <- st_extract(rs, pts)
trn$cls <- rep("land", 20)
trn$cls[c(8, 14, 15, 18, 19)] <- "water"
model <- MASS::lda(cls ~ ., st_drop_geometry(trn))
pr <- predict(rs, model)


#| fig-cap: "Linear discriminant classifier for land/water, based on training data of @fig-rsample"
#| out.width: 70%
#| code-fold: true
plot(pr[1], key.pos = 4, key.width = lcm(3.5), key.length = lcm(2))


#| fig.cap: "Six 30 m Landsat bands downsampled to 90m for Olinda, Br."
#| code-fold: true
plot(r)


#| fig.cap: "Two colour composites"
#| code-fold: true
par(mfrow = c(1, 2))
plot(r, rgb = c(3,2,1), reset = FALSE, main = "RGB")    # rgb
plot(r, rgb = c(4,3,2), main = "False colour (NIR-R-G)") # false colour



log(r)
r + 2 * log(r)



r2 <- r
r2[r < 50] <- NA
r2



r2[is.na(r2)] <- 0
r2



st_apply(r, c("x", "y"), mean)



ndvi <- function(b1, b2, b3, b4, b5, b6) (b4 - b3)/(b4 + b3)
st_apply(r, c("x", "y"), ndvi)



ndvi2 <- function(x) (x[4]-x[3])/(x[4]+x[3])



st_apply(r, c("band"), mean) |> as.data.frame()



st_apply(r, c("band"), quantile, c(.25, .5, .75))



st_apply(r, c("x", "y"), quantile, c(.25, .5, .75))



options("rgdal_show_exportToProj4_warnings"="none") # led to:
# Warning in sp::proj4string(obj): CRS object has comment, which is
# lost in output



load("data/air.rda") # this loads several datasets in .GlobalEnv
dim(air)
stations |>
	st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
	st_geometry() -> st
d <- st_dimensions(station = st, time = dates)
(aq <- st_as_stars(list(PM10 = air), dimensions = d))


#| fig.cap: "Space time diagram of PM$_{10}$ measurements by time and station"
#| code-fold: true
par(mar = c(5.1, 4.1, 0.3, 0.1))
image(aperm(log(aq), 2:1), main = NULL)


#| fig.cap: "Locations of PM$_{10}$ measurement stations, showing mean values"
#| code-fold: true
de_nuts1 <- read_sf("data/de_nuts1.gpkg")
st_as_sf(st_apply(aq, 1, mean, na.rm = TRUE)) |>
    plot(reset = FALSE, pch = 16, extent = de_nuts1)
st_union(de_nuts1) |> plot(add = TRUE)



(a <- aggregate(aq, de_nuts1, mean, na.rm = TRUE))


#| fig.cap: "Areal mean PM$_{10}$ values, for six days"
library(tidyverse)
a |> filter(time >= "2008-01-01", time < "2008-01-07") |> 
	plot(key.pos = 4)


#| fig.cap: "Areal mean PM$_{10}$ time series for a single state"
library(xts) |> suppressPackageStartupMessages()
plot(as.xts(a)[,4], main = de_nuts1$NAME_1[4])


#| fig.height: 5
#| fig.cap: "Origin destination data zones for Bristol, UK, with zone 33 (E02003043) coloured red"
#| code-fold: true
library(spDataLarge)
plot(st_geometry(bristol_zones), axes = TRUE, graticule = TRUE)
plot(st_geometry(bristol_zones)[33], col = 'red', add = TRUE)



head(bristol_od)



nrow(bristol_zones)^2 - nrow(bristol_od) 



# create O-D-mode array:
bristol_tidy <- bristol_od |> 
	select(-all) |> 
	pivot_longer(3:6, names_to = "mode", values_to = "n")
head(bristol_tidy)



od <- bristol_tidy |> pull("o") |> unique()
nod <- length(od)
mode <- bristol_tidy |> pull("mode") |> unique()
nmode = length(mode)
a = array(0L,  c(nod, nod, nmode), 
	dimnames = list(o = od, d = od, mode = mode))
dim(a)



a[as.matrix(bristol_tidy[c("o", "d", "mode")])] <- 
		bristol_tidy$n



order <- match(od, bristol_zones$geo_code)
zones <- st_geometry(bristol_zones)[order]



library(stars)
(d <- st_dimensions(o = zones, d = zones, mode = mode))



(odm <- st_as_stars(list(N = a), dimensions = d))



plot(adrop(odm[,,33]) + 1, logz = TRUE)



d <- st_apply(odm, 2, sum)
which.max(d[[1]])



st_apply(odm, 1:2, sum)



st_apply(odm, c(1,3), sum)



st_apply(odm, c(2,3), sum)



o <- st_apply(odm, 1, sum)



d <- st_apply(odm, 2, sum)



x <- (c(o, d, along = list(od = c("origin", "destination"))))
plot(x, logz = TRUE)



library(units)
a <- set_units(st_area(st_as_sf(o)), km^2)
o$sum_km <- o$sum / a
d$sum_km <- d$sum / a
od <- c(o["sum_km"], d["sum_km"], along = 
		list(od = c("origin", "destination")))
plot(od, logz = TRUE)


#| fig.cap: "Rasterising vector geometry using `st_as_stars`"
file <- system.file("gpkg/nc.gpkg", package="sf")
read_sf(file) |> 
	st_geometry() |>
	st_as_stars() |>
	plot(key.pos = 4)



library(dplyr)
read_sf(file) |>
	mutate(name = as.factor(NAME)) |>
	select(SID74, SID79, name) |>
	st_rasterize()


#| fig.cap: "Rasterising the North Carolina county boundaries"
read_sf(file) |>
	st_cast("MULTILINESTRING") |>
	select(CNTY_ID) |>
	st_rasterize() |>
	plot(key.pos = 4)



x <- st_crs("OGC:CRS84")
x$proj4string



sf_proj_search_paths()



paste0(tail(sf_proj_search_paths(), 1), .Platform$file.sep, 
	   "proj.db")



sf_extSoftVersion()["PROJ"]



sf_proj_network()



sf_proj_network(TRUE)



list.files(sf_proj_search_paths()[1], full.names = TRUE)



(p <- sf_proj_pipelines("OGC:CRS84", "EPSG:22525"))



sf_proj_network(FALSE)
sf_proj_pipelines("OGC:CRS84", "EPSG:22525")



names(p)



p |> pull(accuracy)



p |> filter(is.na(accuracy))



## st_axis_order(TRUE)



st_crs(4326)$axes
st_crs(4326)$ud_unit
st_crs("EPSG:2053")$axes
st_crs("EPSG:2053")$ud_unit



tif <- system.file("tif/L7_ETMs.tif", package = "stars")
read_stars(tif) |>
	st_transform('OGC:CRS84')



read_stars(tif) |>
	st_warp(crs = st_crs('OGC:CRS84')) |>
	st_dimensions()



r <- read_stars(tif)
grd <- st_bbox(r) |>
		st_as_sfc() |>
		st_transform('OGC:CRS84') |>
		st_bbox() |>
		st_as_stars(nx = dim(r)["x"], ny = dim(r)["y"])
st_warp(r, grd)

