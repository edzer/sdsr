book:
	./_build.sh

clean:
	rm -fr _book/*
	rm -fr _bookdown_files
	Rscript -e "bookdown::clean_book(TRUE)"

view:
	google-chrome _book/index.html

tangle:
	Rscript -e "knitr::purl('06-Attributes.Rmd', documentation=0)"

pdf:
	Rscript --quiet _render.R "bookdown::pdf_book"

push:
	make
	git add _book/*html _book/sdsr_files/figure-html/* #_book/images/*png
	git commit -a -m 'commit'
	git push
