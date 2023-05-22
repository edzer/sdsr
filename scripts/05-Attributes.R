#| code-fold: true
#| collapse: false
library(sf) |> suppressPackageStartupMessages()
library(dplyr) |> suppressPackageStartupMessages()
system.file("gpkg/nc.gpkg", package="sf") |>
	read_sf() |>
	st_transform(32119) |>
	select(BIR74, SID74, NAME) |>
	st_centroid() -> x


#| code-fold: true
#| fig.cap: "SID74 total incidences aggregated to four areas"
#| fig.height: 4
nc <- read_sf(system.file("gpkg/nc.gpkg", package = "sf"))
# encode quadrant by two logicals:
nc$lng <- st_coordinates(st_centroid(st_geometry(nc)))[,1] > -79
nc$lat <- st_coordinates(st_centroid(st_geometry(nc)))[,2] > 35.5
nc.grp <- aggregate(nc["SID74"], list(nc$lng, nc$lat), sum)
nc.grp["SID74"] |> st_transform('EPSG:32119') |>
  plot(graticule = TRUE, axes = TRUE)


#| fig.cap: "Example target blocks plotted over North Carolina counties"
#| code-fold: true
nc <- st_transform(nc, 2264)
gr <- st_sf(
   label = apply(expand.grid(1:10, LETTERS[10:1])[,2:1], 1, paste0, collapse = " "),
   geom = st_make_grid(nc))
plot(st_geometry(nc), reset = FALSE, border = 'grey')
plot(st_geometry(gr), add = TRUE)


#| code-fold: true
#| collapse: false
a <- aggregate(nc["SID74"], gr, sum)
c(sid74_sum_counties = sum(nc$SID74),
  sid74_sum_rectangles = sum(a$SID74, na.rm = TRUE))


#| out.width: 40%
#| fig.cap: "Example data for area-weighted interpolation"
#| code-fold: true
g <- st_make_grid(st_bbox(st_as_sfc("LINESTRING(0 0,1 1)")), n = c(2,2))
par(mar = rep(0,4))
plot(g)
plot(g[1] * diag(c(3/4, 1)) + c(0.25, 0.125), add = TRUE, lty = 2)
text(c(.2, .8, .2, .8), c(.2, .2, .8, .8), c(1,2,4,8), col = 'red')

