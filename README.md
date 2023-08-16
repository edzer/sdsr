# [quarto](https://quarto.org/) sources for ["Spatial Data Science: With Applications in R"](https://r-spatial.org/book)

The print version of this book is available from [CRC/Chapman and Hall](https://www.routledge.com/Spatial-Data-Science-With-Applications-in-R/Pebesma-Bivand/p/book/9781138311183). A complete online version of this book [is available](https://r-spatial.org/book).

To recreate/reproduce this book:

* git clone this repository
* download the [data used in Ch 13](https://uni-muenster.sciebo.de/s/8mEbeHPOX9GdAYn), and extract the contents of the `aq` subdirectory into `sdsr/aq`
* install R package dependencies [listed below](#dependencies)
* install [quarto](https://quarto.org/) 
* run `quarto render --to html`

See also the [Dockerfile](https://github.com/edzer/sdsr/tree/main/docker); building the (18 Gb) image with
```
docker build . -t sdsr
```
and running it with
```
docker run -p 8787:8787 -e DISABLE_AUTH=true -ti --rm sdsr
```
will serve an Rstudio server instance on <http://localhost:8787/>, without authentication.


## Compiling the whole book

After running the docker image and opening `rstudio` in the browser:

* click on `01-hello.qmd` in the bottom-right pane
* click on the `Render` button of the top-left pane to compile the whole book

this should open a new browser window with the full book rendered (you may need to switch off popup blockers for localhost)

## Running selected chunks 

To run a selected code section, possibly after modification, find the selected code section in the corresponding `.qmd` file, and click the small green arrow symbols on the top-right corner of the code blocks:

* to prepare, first click `Run All Chunks Above`,
* to run a selected code chunk: click `Run Current Chunk`

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
options(timeout = 1200); install.packages("starsdata", repos = "http://cran.uni-muenster.de/pebesma", type = "source")
```

Install `spatialreg` from source from github, either from source:
```
install.packages("remotes")
remotes::install_github("r-spatial/spatialreg")
```
or as binary from `r-universe`:
```
options(repos = c(
  rspatial = "https://r-spatial.r-universe.dev",
  CRAN = "https://cloud.r-project.org"))
install.packages(c("spatialreg"))
```

## Daily rendered version on GA

The entire book is recreated from source nightly with the latest released R and all updated [CRAN]() packages by a [Github Action](https://github.com/edzer/sdsr/actions) using this [script](https://github.com/edzer/sdsr/blob/main/.github/workflows/publish.yml). The online version thus rendered is found [here](https://edzer.github.io/sdsr/). As this output is not checked daily it is not automatically copied to the "official" online version, at https://r-spatial.org/book/ ; it also lacks portions that require a running PostGIS instance.

## Python version

A version "With Applications in R and Python" is under construction; the sources are in the python branch of this repository, a rendered online version is found at https://r-spatial.org/python/ .
