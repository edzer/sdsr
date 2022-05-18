#!/bin/sh
Rscript -e "bookdown::render_book(input = 'index.Rmd', output_format = bookdown::gitbook(pandoc_args = '--syntax-definition=r.xml')); warnings()"
#Rscript -e "bookdown::render_book(input = 'index.Rmd', output_format = 'bookdown::gitbook')
