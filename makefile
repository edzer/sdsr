all:
	quarto render --to pdf
	make sed
	make 3
	awk -f repl.awk sds.toc > x.toc
	mv x.toc sds.toc
	make 1

sed:
	sed -f repl.sed Spatial*tex > sds.tex

3:
	pdflatex sds.tex
	pdflatex sds.tex
	pdflatex sds.tex

1:
	pdflatex sds.tex

# bookdown::gitbook

cp:
	make zip
	cp sds.zip ~/sciebo/sds_pdf/sds_2023_02_15.zip
	cp sds.pdf ~/sciebo/sds_pdf/sds_2023_02_15.pdf

view:
	google-chrome _book/index.html

push:
	# make
	#git add _book/*html _book/sds_files/figure-html/* _book/libs _book/search_index.json _book/style.css _book/images/*png
	#git commit -a -m 'commit'
	#git push || true

clean:
	rm -fr *_cache *_files _book/*

zip:
	rm -fr sds.zip _book/Spatial-Data-Science.pdf
	zip -r sds.zip sds* krantz.cls images/ *_files/ _book


veryclean: 
	rm -fr data/new_int.RData data/new_ts.RData data/vst.RData

install:
	sudo apt install r-cran-tidyverse r-cran-r2bayesx r-cran-cubelyr r-cran-dbscan r-cran-igraph r-cran-lme4 r-cran-lmtest r-cran-matrixstats r-cran-mgcv r-cran-rgeoda r-cran-rnaturalearth r-cran-rnaturalearthdata r-cran-spdata r-cran-stars r-cran-tmap r-cran-spatstat r-cran-spdep r-cran-spatialreg
	R -q -e 'options(timeout = 600); install.packages("spDataLarge", repos = "https://nowosad.github.io/drat/", type = "source")'
	R -q -e 'options(timeout = 3600); install.packages("starsdata", repos = "http://pebesma.staff.ifgi.de", type = "source")'
	R -q -e 'options(timeout = 600)' -e 'install.packages("INLA", repos = c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"))'
	R -q -e 'install.packages("hglm")'
	R -q -e 'remotes::install_github("r-spatial/gstat")' # for external/no2.csv
	wget https://uni-muenster.sciebo.de/s/8mEbeHPOX9GdAYn/download -O sds.zip
	unzip sds.zip
	mv sds/aq .
	quarto render --to html
