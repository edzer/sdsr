FROM ubuntu:22.04
# minimal docker file to recreate the sdsr book as html
# uses c2d4u for R binary packages to limit install time

MAINTAINER "edzerpebesma" edzer.pebesma@uni-muenster.de

RUN apt-get update
RUN apt-get -y install --no-install-recommends software-properties-common dirmngr wget gpg-agent
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN add-apt-repository ppa:c2d4u.team/c2d4u4.0+

RUN apt-get update
RUN apt-get upgrade -y

# install r-base + r-cran-* packages
RUN export DEBIAN_FRONTEND=noninteractive; apt-get install -y \
	r-base \
	r-cran-dbscan \
	r-cran-gstat \
	r-cran-igraph \
	r-cran-lme4 \
	r-cran-maps \
	r-cran-mapview \
	r-cran-mgcv \
	r-cran-rgeoda \
	r-cran-rnaturalearth \
	r-cran-rnaturalearthdata \
	r-cran-sf \
	r-cran-spatialreg \
	r-cran-spdep \
	r-cran-stars \
	r-cran-tidyverse \
	r-cran-tmap
 
RUN apt-get install -y \
 	git \
	nginx

# install quarto:
RUN wget -q https://github.com/quarto-dev/quarto-cli/releases/download/v1.1.251/quarto-1.1.251-linux-amd64.deb
RUN dpkg -i quarto-1.1.251-linux-amd64.deb

# install book source files:
RUN git clone https://github.com/edzer/sdsr.git

# download & install sdsr data:
RUN wget -q https://uni-muenster.sciebo.de/s/8mEbeHPOX9GdAYn/download -O sds.zip
RUN unzip -q sds.zip
RUN mv sds/aq/* sdsr/aq
RUN rm -fr sds.zip sds

# install R packages not on c2d4u:
## 1. a few CRAN packages:
RUN apt-get -y install libgsl-dev # for R2BayesX
RUN R -e 'install.packages(c("spData", "hglm", "R2BayesX"))'
## 2. INLA:
RUN R -e 'install.packages("INLA", repos = c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"))'
## 3. spDataLarge:
RUN R -e 'options(timeout = 600); install.packages("spDataLarge", repos = "https://nowosad.github.io/drat/",type = "source")'
## 4. starsdata:
RUN R -e 'options(timeout = 600); install.packages("starsdata", repos = "http://pebesma.staff.ifgi.de", type = "source")'

# To be removed when sf 1.0-9 and stars 0.5-7 hit CRAN:
RUN apt-get install -y  libudunits2-dev libgdal-dev libgeos-dev libproj-dev
RUN export DEBIAN_FRONTEND=noninteractive; apt-get install -y r-cran-remotes
RUN R -e 'remotes::install_github("r-spatial/sf")'
RUN R -e 'remotes::install_github("r-spatial/stars")'

RUN R -e 'install.packages(c("matrixStats"))'
RUN export DEBIAN_FRONTEND=noninteractive; apt-get install -y \
	r-cran-lmtest
## generate html book; result found in sdsr/_book
RUN (cd sdsr; quarto render --to html)

## serve the book on http://localhost:80
EXPOSE 80
STOPSIGNAL SIGQUIT
COPY default /etc/nginx/sites-available/default
EXPOSE 80/tcp
CMD ["nginx", "-g", "daemon off;"]
