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

clean:
	rm -fr *_cache *_files _book/*

zip:
	rm -fr sds.zip
	zip -r sds.zip Spatial-Data-Science* krantz.cls images/ *_files/ _book


veryclean: 
	rm -fr data/new_int.RData data/new_ts.RData data/vst.RData
