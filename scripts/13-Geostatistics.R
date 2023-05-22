
# try(load("data/ch16.rda"))
set.seed(123)



files <- list.files("aq", pattern = "*.csv", full.names = TRUE)
files <- setdiff(files, "aq/AirBase_v8_stations.csv") # metadata file
r <- lapply(files, function(f) read.csv(f))



Sys.setenv(TZ = "UTC") # don't use local time zone
r <- lapply(r, function(f) {
		f$t = as.POSIXct(f$DatetimeBegin) 
		f[order(f$t), ] 
	}
)



r <- r[sapply(r, nrow) > 1000]
names(r) <- sapply(r,
			   function(f) unique(f$AirQualityStationEoICode))
length(r) == length(unique(names(r)))



library(xts) |> suppressPackageStartupMessages()
r <- lapply(r, function(f) xts(f$Concentration, f$t))
aq <- do.call(cbind, r)



sel <- apply(aq, 2, function(x) mean(is.na(x)) < 0.25)
aqsel <- aq[, sel]



library(tidyverse) |> suppressPackageStartupMessages()
read.csv("aq/AirBase_v8_stations.csv", sep = "\t") |>
	as_tibble() |> 
	filter(country_iso_code == "DE",
		   station_type_of_area == "rural",
		   type_of_station == "Background") -> a2



library(sf) |> suppressPackageStartupMessages()
a2.sf <- st_as_sf(a2, crs = 'OGC:CRS84',
  coords = c("station_longitude_deg", "station_latitude_deg"))



sel <- colnames(aqsel) %in% a2$station_european_code
aqsel <- aqsel[, sel]
dim(aqsel)



tb <- tibble(NO2 = apply(aqsel, 2, mean, na.rm = TRUE), 
			station_european_code = colnames(aqsel))
crs <- st_crs('EPSG:32632')
right_join(a2.sf, tb) |> st_transform(crs) -> no2.sf 
read_sf("data/de_nuts1.gpkg") |> st_transform(crs) -> de



## library(gstat)
## demo(cokriging)
## demo(cosimulation)



sfc <- st_geometry(a2.sf)[match(colnames(aqsel),
						   a2.sf$station_european_code)] |>
  st_transform(crs)



library(stars)
st_as_stars(NO2 = as.matrix(aqsel)) |>
	st_set_dimensions(names = c("time", "station")) |>
	st_set_dimensions("time", index(aqsel)) |>
	st_set_dimensions("station", sfc) -> no2.st
no2.st



load(file = "data/vst.RData")


library(gstat)


v.st <- variogramST(NO2~1, no2.st[,1:(24*31)], tlags = 0:48, 
	cores = getOption("mc.cores", 2))


save(list = "v.st", file = "data/vst.RData")


#| fig.cap: "Spatiotemporal sample variogram for hourly NO$_2$ concentrations at rural background stations in Germany over 2027; in the right-hand side plot colour corresponds to time lag (yellow is later); distance in m"
#| code-fold: true
#| out.width: 100%
v1 <- plot(v.st)
v2 <- plot(v.st, map = FALSE, legend = list())
print(v1, split = c(1,1,2,1), more = TRUE)
print(v2, split = c(2,1,2,1), more = FALSE)



# product-sum
prodSumModel <- vgmST("productSum",
	space = vgm(150, "Exp", 200000, 0),
	time = vgm(20, "Sph", 6, 0),
	k = 2)
#v.st$dist = v.st$dist / 1000
StAni <- estiStAni(v.st, c(0,200000))
(fitProdSumModel <- fit.StVariogram(v.st, prodSumModel,
	fit.method = 7, stAni = StAni, method = "L-BFGS-B",
	control = list(parscale = c(1,100000,1,1,0.1,1,10)),
	lower = rep(0.0001, 7)))


#| fig.cap:  "Product-sum model, fitted to the spatiotemporal sample variogram"
#| code-fold: true
plot(v.st, fitProdSumModel, wireframe = FALSE, all = TRUE, 
     scales = list(arrows = FALSE), zlim = c(0, 150))

#| fig.cap: "Wireframe plot of the fitted spatiotemporal variogram model"
#| code-fold: true
plot(v.st, model = fitProdSumModel, wireframe = TRUE, all = TRUE, 
	 scales = list(arrows = FALSE), zlim = c(0, 185))



set.seed(1331)
pt <- st_sample(de, 2)
t <- st_get_dimension_values(no2.st, 1)
st_as_stars(list(pts = matrix(1, length(t), length(pt)))) |>
	st_set_dimensions(names = c("time", "station")) |>
	st_set_dimensions("time", t) |>
	st_set_dimensions("station", pt) -> new_pt



load("data/new_ts.RData")


no2.st <- st_transform(no2.st, crs)
new_ts <- krigeST(NO2~1, data = no2.st["NO2"], newdata = new_pt,
	     nmax = 50, stAni = StAni, modelList = fitProdSumModel,
		 progress = FALSE)


save(list = "new_ts", file = "data/new_ts.RData")


#| fig.cap: "Time series plot of spatiotemporal predictions for two points"
#| code-fold: true
plot(as.xts(new_ts[2]))



st_bbox(de) |>
  st_as_stars(dx = 10000) |>
  st_crop(de) -> grd
d <- dim(grd)
t4 <- t[(1:4 - 0.5) * (3*24*30)]
st_as_stars(pts = array(1, c(d[1], d[2], time = length(t4)))) |>
	st_set_dimensions("time", t4) |>
	st_set_dimensions("x", st_get_dimension_values(grd, "x")) |>
	st_set_dimensions("y", st_get_dimension_values(grd, "y")) |>
	st_set_crs(crs) -> grd.st



load("data/new_int.RData")


new_int <- krigeST(NO2~1, data = no2.st["NO2"], newdata = grd.st,
         nmax = 200, stAni = StAni, modelList = fitProdSumModel,
         progress = FALSE)
names(new_int)[2] = "NO2"


save(list = "new_int", file = "data/new_int.RData")


#| fig.cap: "Spatiotemporal predictions for four selected time slices"
#| code-fold: true
library(viridis)
library(viridisLite)
library(ggplot2)
g <- ggplot() + coord_equal() +
	scale_fill_viridis() +
    theme_void() +
	scale_x_discrete(expand=c(0,0)) +
	scale_y_discrete(expand=c(0,0))
g + geom_stars(data = new_int, aes(fill = NO2, x = x, y = y)) + 
    facet_wrap(~as.Date(time)) +
    geom_sf(data = st_cast(de, "MULTILINESTRING")) + 
    geom_sf(data = no2.sf, col = 'grey', cex = .5) + 
	coord_sf(lims_method = "geometry_bbox")



save(list = ls(), file = "ch13.RData")

