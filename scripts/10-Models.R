#| code-fold: true
library(tidyverse)
library(sf)
system.file("gpkg/nc.gpkg", package="sf") |>
	read_sf() -> nc



nc |> mutate(SID = SID74/BIR74, NWB = NWBIR74/BIR74) -> nc1
lm(SID ~ NWB, nc1) |>
  predict(nc1, interval = "prediction") -> pr
bind_cols(nc, pr) |> names()

