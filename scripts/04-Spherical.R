#| code-fold: true
#| collapse: false
#| out.width: 100%
library(sf) |> suppressPackageStartupMessages()
library(maps) |> suppressPackageStartupMessages()
library(dplyr) |> suppressPackageStartupMessages()
map(fill = TRUE, plot = FALSE) |>
  st_as_sf() |>
  filter(ID == "Antarctica") -> a
st_bbox(a)


#| code-fold: true
#| collapse: false
library(s2)
s2_bounds_cap(a)


#| code-fold: true
#| collapse: false
s2_bounds_rect(a)


#| code-fold: true
#| collapse: false
map(fill = TRUE, plot = FALSE) |>
  st_as_sf() |>
  filter(ID == "Fiji") -> Fiji
st_bbox(Fiji)


#| code-fold: true
#| collapse: false
s2_bounds_rect(Fiji)


#| code-fold: true
#| out.height: 70%
#| fig.cap: "Antarctica polygon, (a, c): _not_ passing through `POINT(-180 -90)`; (b, d): passing through `POINT(-180 -90)` and `POINT(180 -90)`"
# maps:
par(mfrow = c(2,2))
par(mar = c(1,1.2,1,1))
m <- st_as_sf(map(fill=TRUE, plot=FALSE))
m <- m[m$ID == "Antarctica", ]
plot(st_geometry(m), asp = 2)
title("a (not valid)")
# ne:
library(rnaturalearth)
ne <- ne_countries(returnclass = "sf")
ne <- ne[ne$region_un == "Antarctica", "region_un"]
plot(st_geometry(ne), asp = 2)
title("b (valid)")
# 3031
m |>
  st_geometry() |>
  st_transform(3031) |>
  plot()
title("c (valid)")
ne |>
  st_geometry() |>
  st_transform(3031) |>
  plot()
title("d (not valid)")

