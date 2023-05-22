
## # after running the "Spatiotemporal Geostatistics" chapter, write NO2 means by
## write_csv(right_join(a2, tb), "no2.csv")



library(tidyverse) |> suppressPackageStartupMessages()
no2 <- read_csv(system.file("external/no2.csv", 
    package = "gstat"), show_col_types = FALSE)



library(sf)
crs <- st_crs("EPSG:32632")
st_as_sf(no2, crs = "OGC:CRS84", coords = 
	c("station_longitude_deg", "station_latitude_deg")) |>
	st_transform(crs) -> no2.sf



read_sf("data/de_nuts1.gpkg") |> st_transform(crs) -> de


#| fig.cap: "Mean NO$_2$ concentrations in air for rural background stations in Germany, in 2017"
#| code-fold: true
ggplot() + geom_sf(data = de) + 
	geom_sf(data = no2.sf, mapping = aes(col = NO2))



library(stars) |> suppressPackageStartupMessages()
st_bbox(de) |>
  st_as_stars(dx = 10000) |>
  st_crop(de) -> grd
grd



library(gstat)
i <- idw(NO2~1, no2.sf, grd)


#| fig.cap: "Inverse distance weighted interpolated values for NO$_2$ over Germany"
#| code-fold: true
ggplot() + geom_stars(data = i, 
					  aes(fill = var1.pred, x = x, y = y)) + 
    xlab(NULL) + ylab(NULL) +
	geom_sf(data = st_cast(de, "MULTILINESTRING")) + 
	geom_sf(data = no2.sf)



v <- variogram(NO2~1, no2.sf)


#| fig.cap: "Sample variogram plot"
#| code-fold: true
plot(v, plot.numbers = TRUE, xlab = "distance h [m]",
	 ylab = expression(gamma(h)),
	 xlim = c(0, 1.055 * max(v$dist)))



v0 <- variogram(NO2~1, no2.sf, cutoff = 100000, width = 10000)


#| fig.cap: "Sample variogram plot with adjusted cutoff and lag width"
#| code-fold: true
plot(v0, plot.numbers = TRUE, xlab = "distance h [m]",
	 ylab = expression(gamma(h)),
	 xlim = c(0, 1.055 * max(v0$dist)))



v.m <- fit.variogram(v, vgm(1, "Exp", 50000, 1))


#| fig.cap: "Sample variogram (circles) with models fitted using weighted least squares (solid line) and maximum likelihood estimation (dashed line)"
#| code-fold: true
fit.variogram_ml <- function(formula, data, init, ...) {
  stopifnot(nrow(init) <= 2, inherits(data, "sf"), inherits(formula, "formula"),
    inherits(init, "variogramModel"))
  if (nrow(init) == 2)
    stopifnot("Nug" %in% init$model)

  # convert from parameter vector to "variogramModel" class:
  # x is c(sill, range) or: c(sill, range, nugget)
  get_model <- function(x, model, min_range = 1e-10) {
    sill <- x[1]
    range <- max(x[2], min_range)
    nugget <- if (length(x) == 3)
		    x[3]
      else
		    0.
    m <- vgm(sill, model, range, nugget)
  }

  # with A <- chol(Q), solve Q x = b for x:
  ch_solve <- function(A, b) {
    backsolve(A, forwardsolve(A, b, upper.tri = TRUE, transpose = TRUE))
  }

  # negative log likelihood, excluding the constant:
  nll <- function(x, d, res, model, ...) {
    m <- get_model(x, model, ...)
    Q <- variogramLine(m, dist_vector = d, covariance = TRUE)
    Qc <- chol(Q)
    det <- 2 * sum(log(diag(Qc)))
    det + t(res) %*% ch_solve(Qc, res)
  }

  # distance matrix, for optim stability rescaled to [0,1] range
  d <- st_distance(data) |> units::drop_units()
  max_d <- max(d)
  d <- d / max_d
  # residuals y - X beta: scale to sd 1
  res <- residuals(lm(formula, data))
  v <- var(res)
  res <- res/sqrt(v)
  if (nrow(init) == 2) {
    o.init <- c(init$psill[2], init$range[2], init$psill[0])
    model <- as.character(init$model[2])
  } else {
    o.init <- c(init$psill[1], init$range[1])
    model <- as.character(init$model[1])
  }
  o.init[2] <- o.init[2] / max_d # scale to [0,1]
  o.init[-2] <- o.init[-2] / v   # scale to sd 1
  o <- optim(o.init, nll, d = d, res = res, model = model, 
			lower = rep(0, length(o.init)), method = "L-BFGS-B", ...)$par
  o[2] <- o[2] * max_d # scale back to distance units
  o[-2] <- o[-2] * v # scale back to variance v
  get_model(o, model)
}
# use WLS fit model v.m for initial parameters:
v.ml <- fit.variogram_ml(NO2~1, no2.sf, v.m)
# plot(v, v.ml, plot.numbers = TRUE)
# plot(v, v.m, plot.numbers = TRUE) ## draws a single model; draw 2 models in single plot:
par(xaxs = "i", yaxs = "i")
plot(gamma ~ dist, v, 
	 xlim = c(0, 1.075 * max(v$dist)), ylim = c(0, 1.05 * max(v$gamma)),
     xlab = "distance h [m]", ylab = expression(gamma(h)))
lines(variogramLine(v.m, 1.075 * max(v$dist)), lty = 1, col = 'blue')
lines(variogramLine(v.ml, 1.075 * max(v$dist)), lty = 2, col = 'blue')
text(v$dist, v$gamma, v$np, pos = 4)



k <- krige(NO2~1, no2.sf, grd, v.m)

#| fig.cap: "Kriged NO$_2$ concentrations over Germany"
ggplot() + geom_stars(data = k, aes(fill = var1.pred, x = x, y = y)) + 
    xlab(NULL) + ylab(NULL) +
	geom_sf(data = st_cast(de, "MULTILINESTRING")) + 
	geom_sf(data = no2.sf) +
	coord_sf(lims_method = "geometry_bbox")



a <- aggregate(no2.sf["NO2"], by = de, FUN = mean)



b <- krige(NO2~1, no2.sf, de, v.m)



b$sample <- a$NO2
b$kriging <- b$var1.pred


#| fig.cap: "Aggregated NO$_2$ values from simple averaging (left) and block kriging (right)"
b |> select(sample, kriging) |> 
		pivot_longer(1:2, names_to = "var", values_to = "NO2") -> b2
b2$var <- factor(b2$var, levels = c("sample", "kriging"))
ggplot() + geom_sf(data = b2, mapping = aes(fill = NO2)) + facet_wrap(~var) +
	 scale_fill_gradientn(colors = sf.colors(20))



SE <- function(x) sqrt(var(x)/length(x))
a <- aggregate(no2.sf["NO2"], de, SE)


#| fig.cap: "Standard errors for mean NO$_2$ values obtained by simple averaging (left) and block kriging (right)"
#| code-fold: true
b$sample <- a$NO2
b$kriging <- sqrt(b$var1.var)
b |> select(sample, kriging) |> 
		pivot_longer(1:2, names_to = "var", 
					 values_to = "Standard_error") -> b2
b2$var <- factor(b2$var, levels = c("sample", "kriging"))
ggplot() +
    geom_sf(data = b2, mapping = aes(fill = Standard_error)) +
    facet_wrap(~var, as.table = FALSE) + 
    scale_fill_gradientn(colors = sf.colors(20))



set.seed(13341)
(s <- krige(NO2~1, no2.sf, grd, v.m, nmax = 30, nsim = 6))


#| fig.cap: "Six conditional simulations for NO$_2$ values"
#| code-fold: true
library(viridis)
g <- ggplot() + coord_equal() +
	scale_fill_viridis() +
    theme_void() +
	scale_x_discrete(expand=c(0,0)) +
	scale_y_discrete(expand=c(0,0))
g + geom_stars(data = s[,,,1:6]) + facet_wrap(~sample)



v <- vroom::vroom("aq/pop/Zensus_Bevoelkerung_100m-Gitter.csv")
v |> filter(Einwohner > 0) |> 
	select(-Gitter_ID_100m) |>
	st_as_sf(coords = c("x_mp_100m", "y_mp_100m"), crs = 3035) |>
	st_transform(st_crs(grd)) -> b
a <- aggregate(b, st_as_sf(grd, na.rm = FALSE), sum)



grd$ID <- 1:prod(dim(grd)) # to identify grid cells
ii <- st_intersects(grd["ID"],
  st_cast(st_union(de), "MULTILINESTRING"), as_points = FALSE)
grd_sf <- st_as_sf(grd["ID"], na.rm = FALSE)[lengths(ii) > 0,]
st_agr(grd_sf) = "identity"
iii <- st_intersection(grd_sf, st_union(de))
grd$area <- st_area(grd)[[1]] + 
    units::set_units(grd$values, m^2)
grd$area[iii$ID] <- st_area(iii)



grd$pop_dens <- a$Einwohner / grd$area
sum(grd$pop_dens * grd$area, na.rm = TRUE) # verify
sum(b$Einwohner)


#| fig.cap: "Population density for 100 m $\\times$ 100 m grid cells"
#| code-fold: true
g + geom_stars(data = grd, aes(fill = sqrt(pop_dens), x = x, y = y))



grd |>
  select("pop_dens") |>
  st_extract(no2.sf) |>
  pull("pop_dens") |> 
  mutate(no2.sf, pop_dens = _) -> no2.sf



summary(lm(NO2~sqrt(pop_dens), no2.sf))


#| fig.cap: "Scatter plot of 2017 annual mean NO$_2$ concentration against population density, for rural background air quality stations"
#| code-fold: true
plot(NO2 ~ sqrt(pop_dens), no2.sf)
abline(lm(NO2 ~ sqrt(pop_dens), no2.sf))



no2.sf <- no2.sf[!is.na(no2.sf$pop_dens),]
vr <- variogram(NO2~sqrt(pop_dens), no2.sf)
vr.m <- fit.variogram(vr, vgm(1, "Exp", 50000, 1))

#| fig.cap: "Residual variogram after subtracting population density trend"
#| code-fold: true
plot(vr, vr.m, plot.numbers = TRUE)



kr <- krige(NO2 ~ sqrt(pop_dens), no2.sf, 
			grd["pop_dens"], vr.m)


#| fig.cap: "Kriging NO$_2$ values using population density as a trend variable"
#| code-fold: true
k$kr1 <- k$var1.pred
k$kr2 <- kr$var1.pred
st_redimension(k[c("kr1", "kr2")], 
	along = list(what = c("kriging", "residual kriging"))) |>
	setNames("NO2") -> km
g + geom_stars(data = km, aes(fill = NO2, x = x, y = y)) + 
	geom_sf(data = st_cast(de, "MULTILINESTRING")) + 
	geom_sf(data = no2.sf) + facet_wrap(~what) +
	coord_sf(lims_method = "geometry_bbox")



save(list = ls(), file = "ch12.RData")

