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
