FROM rocker/geospatial:latest

MAINTAINER "edzerpebesma" edzer.pebesma@uni-muenster.de

RUN apt-get update && apt-get install -y python3-pip python3.10-venv

RUN chown -R rstudio:rstudio /home/rstudio
USER rstudio
WORKDIR /home/rstudio

# 0. install book source files, use python branch:
RUN git clone  https://github.com/edzer/sdsr.git
RUN cd sdsr && git checkout python && cd ..
# RUN mv sdsr /home/rstudio

# 0.1 install python dependencies in venv:
RUN python3 -m venv venv
RUN ./venv/bin/pip3 install -r /home/rstudio/sdsr/requirements.txt

## Add R packages:
RUN R -q -e 'install.packages(c("R2BayesX", "cubelyr", "dbscan", "hglm", "igraph", "lme4", "lmtest", "matrixStats", "mgcv", "reticulate", "rgeoda", "rnaturalearth", "rnaturalearthdata", "spData", "stars"), repos = "https://cloud.r-project.org/")'

# install sf from source, if needed
Run git clone https://github.com/r-spatial/sf.git
Run R CMD INSTALL sf

## 1. download & install sdsr data:
RUN wget -q https://uni-muenster.sciebo.de/s/8mEbeHPOX9GdAYn/download -O sds.zip
RUN unzip -q sds.zip
RUN cp -r sds/aq /home/rstudio/sdsr
RUN rm -fr sds.zip sds
RUN mv sdsr/* .

## 2. spDataLarge:
RUN R -q -e 'options(timeout = 600); install.packages("spDataLarge", repos = "https://nowosad.github.io/drat/", type = "source")'

## 3. starsdata:
RUN R -q -e 'options(timeout = 3600); install.packages("starsdata", repos = "https://cran.uni-muenster.de/pebesma", type = "source")'

## 4. INLA:
RUN R -q -e 'options(timeout = 600)' -e 'install.packages("INLA", repos = c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"))'

RUN quarto render --to html

USER root

EXPOSE 8787

