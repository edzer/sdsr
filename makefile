3:
	pdflatex Spatial*tex
	pdflatex Spatial*tex
	pdflatex Spatial*tex

# bookdown::gitbook

view:
	google-chrome _book/index.html

push:
	# make
	git add _book/*html _book/sds_files/figure-html/* _book/libs _book/search_index.json _book/style.css _book/images/*png
	git commit -a -m 'commit'
	git push || true

purl:
	Rscript -e 'knitr::purl("16-Geostatistics.Rmd")'

clean:
	rm -fr *_cache *_files _book/*

zip:
	rm -fr sds.zip
	zip -r sds.zip Spatial-Data-Science* *_files/ _book
