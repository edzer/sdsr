#| fig.cap: "Earth country boundaries; left: mapping long/lat linearly to $x$ and $y$ (plate carrée); right: as seen from an infinite distance (orthographic)"
#| out.width: 90%
library(sf)
library(rnaturalearth)
w <- ne_countries(scale = "medium", returnclass = "sf")
suppressWarnings(st_crs(w) <- st_crs('OGC:CRS84'))
layout(matrix(1:2, 1, 2), c(2,1))
par(mar = rep(0, 4))
plot(st_geometry(w))

# sphere:
library(s2)
g <- as_s2_geography(TRUE) # Earth
co <- s2_data_countries()
oc <- s2_difference(g, s2_union_agg(co)) # oceans
b <- s2_buffer_cells(as_s2_geography("POINT(-30 -10)"), 9800000) # visible half
i <- s2_intersection(b, oc) # visible ocean
co <- s2_intersection(b, co)
plot(st_transform(st_as_sfc(i), "+proj=ortho +lat_0=-10 +lon_0=-30"), col = 'lightblue')
plot(st_transform(st_as_sfc(co), "+proj=ortho +lat_0=-10 +lon_0=-30"), col = NA, add = TRUE)



## library(sf)
## library(rnaturalearth)
## w <- ne_countries(scale = "medium", returnclass = "sf")
## plot(st_geometry(w))



st_is_longlat(w)



DE <- st_geometry(ne_countries(country = "germany",
							  returnclass = "sf"))
DE |> st_transform("+proj=eqc +lat_ts=51.14 +lon_0=90w") ->
    DE.eqc


#| code-fold: true
print(mean(st_bbox(DE)[c("ymin", "ymax")]), digits = 4)


#| fig.height: 4.5
#| out.width: 70%
#| code-fold: true
#| fig.cap: "Germany in equirectangular projection: with axis units degrees (left) and metres in the equidistant cylindrical projection (right)"
par(mfrow = c(1, 2), mar = c(2.2, 2.2, 0.3, 0.5))
plot(DE, axes = TRUE)
plot(DE.eqc, axes = TRUE)



library(classInt)
# set.seed(1) if needed ?
r <- rnorm(100)
(cI <- classIntervals(r))
cI$brks


#| fig.cap: "Annotating base plots with a legend"
library(sf)
nc <- read_sf(system.file("gpkg/nc.gpkg", package = "sf"))
plot(nc["BIR74"], reset = FALSE, key.pos = 4)
plot(st_buffer(nc[1,1], units::set_units(10, km)), col = 'NA', 
	 border = 'red', lwd = 2, add = TRUE)


#| fig.cap: "Annotated multi-slice stars plot"
library(stars)
system.file("tif/L7_ETMs.tif", package = "stars") |>
	read_stars() -> r
st_bbox(r) |> st_as_sfc() |> st_sample(5) |> 
    st_buffer(300) -> circ
hook <- function() { 
    plot(circ, col = NA, border = 'yellow', add = TRUE)
}
plot(r, hook = hook, key.pos = 4)



library(tidyverse) |> suppressPackageStartupMessages()
nc.32119 <- st_transform(nc, 32119) 
year_labels <- 
	c("SID74" = "1974 - 1978", "SID79" = "1979 - 1984")
nc.32119 |> select(SID74, SID79) |> 
	pivot_longer(starts_with("SID")) -> nc_longer


## ggplot() + geom_sf(data = nc_longer, aes(fill = value), linewidth = 0.4) +
##   facet_wrap(~ name, ncol = 1,
## 			 labeller = labeller(name = year_labels)) +
##   scale_y_continuous(breaks = 34:36) +
##   scale_fill_gradientn(colours = sf.colors(20)) +
##   theme(panel.grid.major = element_line(colour = "white"))


#| fig.cap: "Simple facet raster plot with `ggplot2` and `geom_stars`"
library(ggplot2)
library(stars)
r <- read_stars(system.file("tif/L7_ETMs.tif", package = "stars"))
ggplot() + geom_stars(data = r) +
		facet_wrap(~band) + coord_equal() +
		theme_void() +
        scale_x_discrete(expand = c(0,0)) + 
        scale_y_discrete(expand = c(0,0)) +
		scale_fill_viridis_c()



## library(tmap)
## system.file("gpkg/nc.gpkg", package = "sf") |>
##     read_sf() |> st_transform('EPSG:32119') -> nc.32119
## tm_shape(nc.32119) +
## 	tm_polygons(c("SID74", "SID79"), title = "SIDS") +
##     tm_layout(legend.outside = TRUE,
## 			  panel.labels = c("1974-78", "1979-84")) +
##     tm_facets(free.scales=FALSE)


#| fig.cap: "**tmap**: using `tm_polygons()` with two attribute names"
#| code-fold: true
library(tmap)
system.file("gpkg/nc.gpkg", package = "sf") |>
    read_sf() |>
    st_transform('EPSG:32119') -> nc.32119
tm_shape(nc.32119) + tm_polygons(c("SID74", "SID79"), title="SIDS") + 
    tm_layout(legend.outside=TRUE, panel.labels=c("1974-78", "1979-84")) + 
    tm_facets(free.scales=FALSE)



## nc_longer <- nc.32119 |> select(SID74, SID79) |>
## 	pivot_longer(starts_with("SID"), values_to = "SID")
## tm_shape(nc_longer) + tm_polygons("SID") +
##     tm_facets(by = "name")

#| fig.cap: "**tmap**: Using `tm_facets()` on a long table"
nc.32119 |> select(SID74, SID79) |> 
	pivot_longer(starts_with("SID"), values_to = "SID") -> nc_longer
tm_shape(nc_longer) + tm_polygons("SID") + tm_facets(by = "name")



## tm_shape(r) + tm_raster()


#| fig.cap: "Simple raster plot with tmap"
tm_shape(r) + tm_raster()



## tmap_mode("view")



## tmap_mode("plot")



## tmap_mode("view")



## tmap_mode("plot")

