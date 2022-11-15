# quarto sources the "Spatial Data Science" book.

A rendered (html) version of this book is available [here](https://r-spatial.org/book).
The pdf version has been submitted to CRC/Chapman and Hall, for hardcopy publication.

To recreate/reproduce this book:

* git clone this repository
* download the data used in Ch 13 from https://uni-muenster.sciebo.de/s/8mEbeHPOX9GdAYn, and extract the contents of the `aq` subdirectory into `sdsr/aq`
* install all R packages listed under [Dependencies](#dependencies)
* install [quarto](https://quarto.org/) 
* run `quarto render --to html`

See also the [Dockerfile](https://github.com/edzer/sdsr/tree/main/docker); building the (18 Gb) image with
```
docker build . -t sdsr --build-arg TZ=`timedatectl show --property=Timezone | awk -F = '{print $2}'`
```
and running it with
```
docker run -p 8787:8787 -e DISABLE_AUTH=true -ti --rm sdsr
```
will serve an Rstudio server instance on <http://localhost:8787/>, without authentication.

After running the docker image and opening `rstudio` in the browser:

* click on `01-hello.qmd` in the bottom-right pane
* click on the `Render` button of the top-left pane to compile the whole book
* this should open a new browser window with the full book rendered (switch off popup blocker for localhost)
* to run an individual code section, possibly after modification, click the small green arrow symbols on the top-left corner of code blocks:
    * to prepare, first `Run All Chunks Above`,
	* to run the actual section: `Run Current Chunk`

## Dependencies

To locally process the book, download (clone) this repository and install the following R packages from CRAN:

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
options(timeout = 1200); install.packages("starsdata", repos = "http://pebesma.staff.ifgi.de", type = "source")
```

Install `stars` from source from github (not needed after stars >= 0.6-0 is available from CRAN), either from source:
```
install.packages("remotes")
remotes::install_github("r-spatial/stars")
```
or as binary from `r-universe`:
```
options(repos = c(
  rspatial = "https://r-spatial.r-universe.dev",
  CRAN = "https://cloud.r-project.org"))
install.packages(c("stars"))
```
