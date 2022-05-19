#!/bin/sh
Rscript -e "bookdown::render_book(input = 'index.Rmd', output_format = bookdown::gitbook(pandoc_args = c('--syntax-definition=r.xml', '--include-in-header=hide_code.html'))); warnings()"
#Rscript -e "bookdown::render_book(input = 'index.Rmd', output_format = 'bookdown::gitbook')
