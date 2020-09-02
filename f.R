# add/replace <details>...</details> sections to show hidden code blocks below output
f <- function(file, out = file) { 
	r = readLines(file)
	fout = file(out, "wt")
	inside = FALSE
	ignore = FALSE
	for (i in seq_along(r)) {
		if (grepl("<details>", r[i])) { 
			ignore = TRUE
		}
		if (!ignore)
			writeLines(r[i], fout)
		if (grepl("</details>", r[i])) { 
			ignore = FALSE
		}
		if (grepl("echo=FALSE", r[i]) && (!grepl("eval=FALSE", r[i])) && 
				(!grepl("nodetails", r[i]))) { 
			# start of a echo=FALSE section:
			s = character(0) # start
			j = 1
			s[j] = r[i]
			inside = TRUE
		} else { 
			if (inside) {
				j = j + 1
				s[j] = r[i]
			}
			# print & process:
			if (grepl("```", r[i]) && inside) { # end of a section
				inside = FALSE
				writeLines("<details>", fout)
				#writeLines("<summary>Click for Code</summary>", fout)
				writeLines("<summary style=\"color:blue;\">Click for Code</summary>", fout)
				writeLines("```{r eval=FALSE}", fout)
				for (k in 2:j)
					writeLines(s[k], fout)
				writeLines("</details>", fout)
			}
		}
	}
	close(fout)
}
