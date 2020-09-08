## ----load,eval=TRUE,echo=FALSE,nodetails=TRUE---------------------------------
not_loaded = TRUE
#not_loaded = !file.exists("data/ch16.rda")
#if (!not_loaded)
#	load("data/ch16.rda")


















## ----computestationmeans------------------------------------------------------
tb = tibble(NO2 = apply(aqsel, 2, mean, na.rm = TRUE), station_european_code = colnames(aqsel))
crs = 32632
right_join(a2.sf, tb) %>% st_transform(crs) -> no2.sf 
# load German boundaries
data(air, package = "spacetime")
de <- st_transform(st_as_sf(DE_NUTS1), crs)
ggplot() + geom_sf(data = de) +  geom_sf(data = no2.sf, mapping = aes(col = NO2))


## ----computevariogram---------------------------------------------------------
library(gstat)
v = variogram(NO2~1, no2.sf)
plot(v, plot.numbers = TRUE)


## ----controlcutoff------------------------------------------------------------
library(gstat)
v0 = variogram(NO2~1, no2.sf, cutoff = 100000, width = 10000)
plot(v0, plot.numbers = TRUE)


## ----fitvariogrammodel--------------------------------------------------------
v.m = fit.variogram(v, vgm(1, "Exp", 50000, 1))
plot(v, v.m, plot.numbers = TRUE)


## ----buildgridgermany---------------------------------------------------------
# build a grid over Germany:
st_bbox(de) %>%
  st_as_stars(dx = 10000) %>%
  st_set_crs(crs) %>%
  st_crop(de) -> grd
grd


## ----krigeovergermany---------------------------------------------------------
k = krige(NO2~1, no2.sf, grd, v.m)
ggplot() + geom_stars(data = k, aes(fill = var1.pred, x = x, y = y)) + 
	geom_sf(data = st_cast(de, "MULTILINESTRING")) + 
	geom_sf(data = no2.sf)


## -----------------------------------------------------------------------------
a = aggregate(no2.sf["NO2"], by = de, FUN = mean)


## ----krigeNO2regionmeans------------------------------------------------------
b = krige(NO2~1, no2.sf, de, v.m)


## -----------------------------------------------------------------------------
b$sample = a$NO2
b$kriging = b$var1.pred
b %>% select(sample, kriging) %>% gather(var, NO2, -geometry) -> b2
ggplot() + geom_sf(data = b2, mapping = aes(fill = NO2)) + facet_wrap(~var) +
	 scale_fill_gradientn(colors = sf.colors(20))




## -----------------------------------------------------------------------------
b$sample = a$NO2
b$kriging = sqrt(b$var1.var)
b %>% select(sample, kriging) %>% gather(var, NO2, -geometry) -> b2
ggplot() + geom_sf(data = b2, mapping = aes(fill = NO2)) + facet_wrap(~var) +
	 scale_fill_gradientn(colors = sf.colors(20))


## ----plotkrigingvalues--------------------------------------------------------
library(viridis)
s = krige(NO2~1, no2.sf, grd, v.m, nmax = 30, nsim = 10)
g = ggplot() + coord_equal() +
	scale_fill_viridis() +
    theme_void() +
	scale_x_discrete(expand=c(0,0)) +
	scale_y_discrete(expand=c(0,0))
g + geom_stars(data = s[,,,1:6]) + facet_wrap(~sample)






## -----------------------------------------------------------------------------
grd$pop_dens = a$Einwohner / grd$area
sum(grd$pop_dens * grd$area, na.rm = TRUE) # verify
sum(b$Einwohner)
g + geom_stars(data = grd, aes(fill = sqrt(pop_dens), x = x, y = y))


## -----------------------------------------------------------------------------
(a = aggregate(grd["pop_dens"], no2.sf, mean))
no2.sf$pop_dens = st_as_sf(a)[[1]]
summary(lm(NO2~sqrt(pop_dens), no2.sf))


## ----no2scat, out.width = '50%', fig.cap="Scatter plot of 2017 annual mean NO2 concentration against population density, for rural background air quality stations", eval=TRUE, fig=TRUE, echo=FALSE----
plot(NO2~sqrt(pop_dens), no2.sf)
abline(lm(NO2~sqrt(pop_dens), no2.sf))


## ----predictusingpopulationdensity--------------------------------------------
no2.sf = no2.sf[!is.na(no2.sf$pop_dens),]
vr = variogram(NO2~sqrt(pop_dens), no2.sf)
vr.m = fit.variogram(vr, vgm(1, "Exp", 50000, 1))
plot(vr, vr.m, plot.numbers = TRUE)


## -----------------------------------------------------------------------------
kr = krige(NO2~sqrt(pop_dens), no2.sf, grd["pop_dens"], vr.m)
k$kr1 = k$var1.pred
k$kr2 = kr$var1.pred
st_redimension(k[c("kr1", "kr2")], 
	along = list(what = c("kriging", "residual kriging"))) %>%
	setNames("NO2") -> km
g + geom_stars(data = km, aes(fill = NO2, x = x, y = y)) + 
	geom_sf(data = st_cast(de, "MULTILINESTRING")) + 
	geom_sf(data = no2.sf) + facet_wrap(~what)


## ----eval=FALSE---------------------------------------------------------------
## library(gstat)
## demo(cokriging)
## demo(cosimulation)




## ----plotvariograms-----------------------------------------------------------
v1 = plot(v.st)
v2 = plot(v.st, map = FALSE)
print(v1, split = c(1,1,2,1), more = TRUE)
print(v2, split = c(2,1,2,1), more = FALSE)


## ----fitprodsummodel----------------------------------------------------------
# product-sum
prodSumModel <- vgmST("productSum",
	space=vgm(150, "Exp", 200, 0),
	time= vgm(20, "Sph",   40, 0),
	k=2)
StAni = estiStAni(v.st, c(0,200000))
(fitProdSumModel <- fit.StVariogram(v.st, prodSumModel, fit.method = 7,
	stAni = StAni, method = "L-BFGS-B",
	control = list(parscale = c(1,10,1,1,0.1,1,10)),
	lower = rep(0.0001, 7)))
plot(v.st, fitProdSumModel, wireframe=FALSE, all=TRUE, scales=list(arrows=FALSE), zlim=c(0,150))


## ----plotprodsummodel---------------------------------------------------------
plot(v.st, model=fitProdSumModel, wireframe=TRUE, all=TRUE, scales=list(arrows=FALSE), zlim=c(0,185))


## ----createstarsobject--------------------------------------------------------
set.seed(123)
pt = st_sample(de, 2)
t = st_get_dimension_values(no2.st, 1)
st_as_stars(list(pts = matrix(1, length(t), length(pt)))) %>%
	st_set_dimensions(names = c("time", "station")) %>%
	st_set_dimensions("time", t) %>%
	st_set_dimensions("station", pt) -> new_pt



## -----------------------------------------------------------------------------
plot(xts(t(new_ts[[2]]), t), type = 'l')


## ----spatiotemporalpredictionsgrid--------------------------------------------
t4 = t[(1:4 - 0.5) * (3*24*30)]

d = dim(grd)
st_as_stars(pts = array(1, c(d[1], d[2], time=length(t4)))) %>%
	st_set_dimensions("time", t4) %>%
	st_set_dimensions("x", st_get_dimension_values(grd, "x")) %>%
	st_set_dimensions("y", st_get_dimension_values(grd, "y")) %>%
	st_set_crs(crs) -> grd.st



## -----------------------------------------------------------------------------
g + geom_stars(data = new_int, aes(fill = NO2, x = x, y = y)) + 
	facet_wrap(~time) +
	geom_sf(data = st_cast(de, "MULTILINESTRING")) + 
	geom_sf(data = no2.sf, col = 'grey', cex = .5)

