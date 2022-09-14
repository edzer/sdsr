book:
	./_build.sh

clean:
	rm -fr _book/*
	rm -fr _bookdown_files
	Rscript -e "bookdown::clean_book(TRUE)"

pdf0:
	# quarto render --to pdf
	cp -rp images krantz.cls all.bib 15*_files 16*_files _book
	vi +55 Spatial-Data-Science.tex # edit

pdf1:
	(cd _book; pdflatex ../Spatial-Data-Science.tex; pdflatex ../Spatial-Data-Science.tex; makeindex Spatial-Data-Science.idx; pdflatex ../Spatial-Data-Science.tex)

# bookdown::gitbook

view:
	google-chrome _book/index.html

tangle:
	Rscript -e "knitr::purl('16-Geostatistics.Rmd', documentation=0)"

pdf:
	Rscript --quiet _render.R "bookdown::pdf_book"
	# evince _book/sds.pdf

small:
	mv 0[5-9]*Rmd 1*Rmd tmp

quarto_book:
	#./to_quarto_fold 0[0-5]*Rmd
	#./to_quarto 0[7-9]*Rmd [1-9]*Rmd
	cat book.bib packages.bib > quarto/book.bib
	cp -rp images quarto
	(cd quarto; quarto render)
	# echo '\@ref(fig:foo)' | sed 's/fig:/fig-/g' | sed 's/\\@ref//g'
	# mv quarto/SpatialDataScience.tex xx.tex 
	#sed 's/\DefineVerbatimEnvironment{Highlighting}{Verbatim}{commandchars=\\\{\}}/\DefineVerbatimEnvironment{Highlighting}{Verbatim}{commandchars=\\\{\},fontsize=\small}/ xx.tex > quarto/SpatialDataScience.tex


push:
	# make
	git add _book/*html _book/sds_files/figure-html/* _book/libs _book/search_index.json _book/style.css _book/images/*png
	git commit -a -m 'commit'
	git push || true

purl:
	Rscript -e 'knitr::purl("16-Geostatistics.Rmd")'
	#Rscript -e 'knitr::purl("17-Areal.Rmd")'

objects = 01-hello.Rmd

details:
	#Rscript -e 'source("f.R"); f("01-hello.Rmd")'
	#Rscript -e 'source("f.R"); f("02-Spaces.Rmd")'
	#Rscript -e 'source("f.R"); f("03-Geometries.Rmd")'
	#Rscript -e 'source("f.R"); f("04-Spherical.Rmd")'
	Rscript -e 'source("f.R"); f("05-Attributes.Rmd")'
