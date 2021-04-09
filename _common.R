set.seed(1014)
options(digits = 3)

knitr::opts_chunk$set(
  comment = "#",
  collapse = knitr::is_latex_output(),
  cache = TRUE,
#  out.width = "70%",
  fig.align = 'center'
#  fig.width = 6,
#  fig.asp = 0.618,  # 1 / phi
#  fig.show = "hold"
)

options(dplyr.print_min = 6, dplyr.print_max = 6)
options(stars.crs = 17)
mapview::mapviewOptions(fgb = FALSE)

# for units: (?)
Sys.setenv(UDUNITS2_XML_PATH="")
if (knitr::is_latex_output())
	options(width = 66)
	#options(width = 72)
