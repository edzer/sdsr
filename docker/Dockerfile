FROM rocker/geospatial:latest

MAINTAINER "edzerpebesma" edzer.pebesma@uni-muenster.de

<<<<<<< HEAD
RUN apt-get update && apt-get install -y python3-pip
=======
## Add CRAN packages:
RUN R -q -e 'install.packages(c("R2BayesX", "cubelyr", "dbscan", "hglm", "igraph", "lme4", "lmtest", "matrixStats", "mgcv", "rgeoda", "rnaturalearth", "rnaturalearthdata", "spData", "stars"), repos = "https://cloud.r-project.org/")'
>>>>>>> main

RUN chown -R rstudio:rstudio /home/rstudio
USER rstudio
WORKDIR /home/rstudio

## Add R packages:
RUN R -q -e 'install.packages(c("R2BayesX", "cubelyr", "dbscan", "hglm", "igraph", "lme4", "lmtest", "matrixStats", "mgcv", "reticulate", "rgeoda", "rnaturalearth", "rnaturalearthdata", "spData", "stars"))'

# Add Python packages:
RUN pip3 install geopandas descartes rasterio folium matplotlib mapclassify

# 0. install book source files, use python branch:
RUN git clone  https://github.com/edzer/sdsr.git
RUN cd sdsr && git checkout python && cd ..
RUN mv sdsr /home/rstudio

## 1. download & install sdsr data:
RUN wget -q https://uni-muenster.sciebo.de/s/8mEbeHPOX9GdAYn/download -O sds.zip
RUN unzip -q sds.zip
RUN cp -r sds/aq /home/rstudio/sdsr
RUN rm -fr sds.zip sds
RUN mv sdsr/* .

## 2. spDataLarge:
RUN R -q -e 'options(timeout = 600); install.packages("spDataLarge", repos = "https://nowosad.github.io/drat/", type = "source")'

## 3. starsdata:
RUN R -q -e 'options(timeout = 3600); install.packages("starsdata", repos = "http://cran.uni-muenster.de/pebesma", type = "source")'

## 4. INLA:
RUN R -q -e 'options(timeout = 600)' -e 'install.packages("INLA", repos = c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"))'

USER root

## 5. dev version of spatialreg:
RUN R -q -e 'remotes::install_github("r-spatial/spatialreg")'

EXPOSE 8787
