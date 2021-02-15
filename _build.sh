#!/bin/sh
Rscript -e "bookdown::render_book(input = 'index.Rmd', output_format = 'bookdown::gitbook')"
