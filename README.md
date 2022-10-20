# quarto sources the "Spatial Data Science" book.

A rendered (html) version of this book is available [here](https://r-spatial.org/book).
The pdf version has been submitted to CRC/Chapman and Hall, for hardcopy publication.

To recreate/reproduce this book:

* git clone this repository
* download the data used in Ch 13 from https://uni-muenster.sciebo.de/s/8mEbeHPOX9GdAYn, and extract the contents of the `aq` subdirectory into `sdsr/aq`
* install all R packages listed under [Dependencies](#dependencies)
* install [quarto](https://quarto.org/) 
* run `quarto render --to html`

See also the [Dockerfile](https://github.com/edzer/sdsr/tree/main/docker); building the image with
```
docker build . -t sdsr --build-arg TZ=`timedatectl show --property=Timezone | awk -F = '{print $2}'`
```
and running it with
```
docker run -d -p 80:80 sdsr:latest
```
will serve the html book on http://localhost:80

## Dependencies

To locally process the book, install the following R packages from CRAN:

```
install.packages(c(
  "dbscan",
  "gstat",
  "hglm",
  "igraph",
  "lme4",
  "lmtest",
  "maps",
  "mapview",
  "matrixStats",
  "mgcv",
  "R2BayesX",
  "rgeoda",
  "rnaturalearth",
  "rnaturalearthdata",
  "sf",
  "spatialreg",
  "spdep",
  "spData",
  "stars",
  "tidyverse",
  "tmap"))
```

Install `INLA`:
```
install.packages("INLA", repos = c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"))
```

Install `spDataLarge`:
```
options(timeout = 600); install.packages("spDataLarge", repos = "https://nowosad.github.io/drat/",type = "source")
```
Install `starsdata`:
```
options(timeout = 600); install.packages("starsdata", repos = "http://pebesma.staff.ifgi.de", type = "source")
```

Install `sf` and `stars` from source from github (not needed after sf 1.0-9 and stars 0.5-7 are available from CRAN):
```
# apt-get install -y  libudunits2-dev libgdal-dev libgeos-dev libproj-dev
install.packages("remotes")
remotes::install_github("r-spatial/sf")
remotes::install_github("r-spatial/stars")
```
or as binary from `r-universe`:
```
options(repos = c(
  rspatial = "https://r-spatial.r-universe.dev",
  CRAN = "https://cloud.r-project.org"))
install.packages(c("sf", "stars"))
```
