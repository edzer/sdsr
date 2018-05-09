book:
	./_build.sh

clean:
	rm -fr _book/*

view:
	google-chrome _book/index.html

tangle:
	Rscript -e "knitr::purl('sdsr.Rmd', documentation=0)"

pdf:
	Rscript --quiet _render.R "bookdown::pdf_book"

push:
	make
	git add _book/*html
	git commit -a -m 'commit'
	git push
