book:
	./_build.sh

clean:
	rm -fr _book/*
	rm -fr _bookdown_files
	Rscript -e "bookdown::clean_book(TRUE)"

pdf0:
	Rscript -e 'bookdown::render_book("index.Rmd", "bookdown::pdf_book")'

# bookdown::gitbook

view:
	google-chrome _book/index.html

tangle:
	Rscript -e "knitr::purl('16-Geostatistics.Rmd', documentation=0)"

pdf:
	Rscript --quiet _render.R "bookdown::pdf_book"
	# evince _book/sds.pdf

push:
	# make
	git add _book/*html _book/sds_files/figure-html/* _book/libs _book/search_index.json _book/style.css _book/images/*png
	git commit -a -m 'commit'
	git push || true

purl:
	#Rscript -e 'knitr::purl("16-Geostatistics.Rmd")'
	Rscript -e 'knitr::purl("17-Areal.Rmd")'

objects = 01-hello.Rmd

details:
	#Rscript -e 'source("f.R"); f("01-hello.Rmd")'
	#Rscript -e 'source("f.R"); f("02-Spaces.Rmd")'
	#Rscript -e 'source("f.R"); f("03-Geometries.Rmd")'
	#Rscript -e 'source("f.R"); f("04-Spherical.Rmd")'
	Rscript -e 'source("f.R"); f("05-Attributes.Rmd")'
