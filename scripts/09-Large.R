
library(sf)
file <- system.file("gpkg/nc.gpkg", package = "sf")
c(xmin = -82,ymin = 36, xmax = -80, ymax = 37) |>
    st_bbox() |> st_as_sfc() |> st_as_text() -> bb
read_sf(file, wkt_filter = bb) |> nrow() # out of 100



q <- "select BIR74,SID74,geom from 'nc.gpkg' where BIR74 > 1500"
read_sf(file, query = q) |> nrow()



q <- "select BIR74,SID74,geom from 'nc.gpkg' LIMIT 10 OFFSET 50"
read_sf(file, query = q) |> nrow()



has_PG <- any("PostgreSQL" %in% st_drivers()$name) &&
	!inherits(try(DBI::dbConnect( 
    	RPostgres::Postgres(),
    	host = "localhost",
    	dbname = "postgis"), silent = TRUE), "try-error")













## download.file(paste0("https://openstreetmap.org/api/0.6/map?",
## 	   "bbox=7.595,51.969,7.598,51.970"),
##     "data/ms.osm", method = "auto")


#| fig.cap: "OpenStreetMap vector data"
#| fig.height: 5
o <- read_sf("data/ms.osm", "lines")
p <- read_sf("data/ms.osm", "multipolygons")
st_bbox(c(xmin = 7.595, ymin = 51.969, 
		  xmax = 7.598, ymax = 51.970), crs = 'OGC:CRS84') |>
    st_as_sfc() |>
    plot(axes = TRUE, lwd = 2, lty = 2, cex.axis = .5)
plot(o[,1], lwd = 2, add = TRUE)
plot(st_geometry(p), border = NA, col = '#88888888', add = TRUE)



## options(timeout = 600) # or large in case of slow network
## install.packages("starsdata",
##     repos = "http://pebesma.staff.ifgi.de", type = "source")



library(stars) |> suppressPackageStartupMessages()
f <- paste0("sentinel/S2A_MSIL1C_20180220T105051_N0206",
		   "_R051_T32ULE_20180221T134037.zip")
granule <- system.file(file = f, package = "starsdata")
file.size(granule)
base_name <- strsplit(basename(granule), ".zip")[[1]]
s2 <- paste0("SENTINEL2_L1C:/vsizip/", granule, "/", base_name, 
	".SAFE/MTD_MSIL1C.xml:10m:EPSG_32632")
(p <- read_stars(s2, proxy = TRUE))
object.size(p)



p2 <- p * 2


#| fig.cap: "Downsampled 10 m bands of a Sentinel-2 scene"
#| code-fold: true
plot(p)


#| code-fold: true
(ds <- floor(sqrt(prod(dim(p)) / prod(dev.size("px")))))



methods(class = "stars_proxy")



dsn = paste0('ZARR:"/vsicurl/https://ncsa.osn.xsede.org',
	   '/Pangeo/pangeo-forge/gpcp-feedstock/gpcp.zarr"')



library(stars)
bounds = c(longitude = "lon_bounds", latitude = "lat_bounds")
(r = read_mdim(dsn, bounds = bounds, count = c(NA, NA, 10)))
st_bbox(r)

