FROM rocker/geospatial:latest

MAINTAINER "edzerpebesma" edzer.pebesma@uni-muenster.de
ARG TZ # used by R package lubridate

## Add some CRAN packages:
RUN R -q -e 'install.packages(c("R2BayesX", "cubelyr", "dbscan", "hglm", "igraph", "lme4", "lmtest", "matrixStats", "mgcv", "rgeoda", "rnaturalearth", "rnaturalearthdata", "spData"))'

# 0. install book source files:
RUN git clone https://github.com/edzer/sdsr.git
RUN mv sdsr/* /home/rstudio
RUN rm -fr sdsr

## 1. download & install sdsr data:
RUN wget -q https://uni-muenster.sciebo.de/s/8mEbeHPOX9GdAYn/download -O sds.zip
RUN unzip -q sds.zip
RUN mv sds/aq/* /home/rstudio/aq
RUN rm -fr sds.zip sds
RUN chown -R rstudio:rstudio /home/rstudio/*

## 2. spDataLarge:
RUN R -q -e 'options(timeout = 600); install.packages("spDataLarge", repos = "https://nowosad.github.io/drat/", type = "source")'

## 3. starsdata:
RUN R -q -e 'options(timeout = 3600); install.packages("starsdata", repos = "http://pebesma.staff.ifgi.de", type = "source")'

## 4. INLA:
RUN R -q -e 'options(timeout = 600)' -e 'install.packages("INLA", repos = c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"))'

# 5. Install stars dev (can be removed when stars 0.5-7 hits CRAN):
RUN R -q -e 'remotes::install_github("r-spatial/stars")'

ARG EVAL_INLA # if false, do not run INLA code

EXPOSE 8787
