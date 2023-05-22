
knitr::opts_chunk$set(echo = TRUE, paged.print=FALSE)
owidth <- getOption("width")
xargs <- function(x) {
    o <- capture.output(args(x))
    oo <- gsub("  *", " ", paste(o[-length(o)], collapse=""))
    ooo <- strwrap(oo, width=getOption("width"), indent=1, exdent=3)
    cat(paste(ooo, collapse="\n"), "\n")
}



library(sf)



data(pol_pres15, package = "spDataLarge")
pol_pres15 |>
    subset(select = c(TERYT, name, types)) |>
    head()

#| out.width: 100%
#| fig.cap: "Polish municipality types 2015"
library(tmap, warn.conflicts = FALSE)
tm_shape(pol_pres15) + tm_fill("types")



if (!all(st_is_valid(pol_pres15)))
		pol_pres15 <- st_make_valid(pol_pres15)



library(spdep) |> suppressPackageStartupMessages()



## args(poly2nb)


xargs(poly2nb)



pol_pres15 |> poly2nb(queen = TRUE) -> nb_q



nb_q



old_use_s2 <- sf_use_s2()



sf_use_s2(TRUE)



(pol_pres15 |> st_transform("OGC:CRS84") -> pol_pres15_ll) |> 
    poly2nb(queen = TRUE) -> nb_q_s2



all.equal(nb_q, nb_q_s2, check.attributes=FALSE)



(nb_q |> n.comp.nb())$nc



library(Matrix, warn.conflicts = FALSE)
library(spatialreg, warn.conflicts = FALSE)
nb_q |> 
    nb2listw(style = "B") |> 
    as("CsparseMatrix") -> smat
library(igraph, warn.conflicts = FALSE)
(smat |> graph.adjacency() -> g1) |> 
    count_components()



tf <- tempfile(fileext = ".gal")
write.nb.gal(nb_q, tf)



pol_pres15 |> 
    st_geometry() |> 
    st_centroid(of_largest_polygon = TRUE) -> coords 
(coords |> tri2nb() -> nb_tri)



nb_tri |> 
    nbdists(coords) |> 
    unlist() |> 
    summary()



(nb_tri |> n.comp.nb())$nc



(nb_tri |> 
        soi.graph(coords) |> 
        graph2nb() -> nb_soi)



(nb_soi |> n.comp.nb() -> n_comp)$nc



table(n_comp$comp.id)


#| fig.cap: "Triangulated (orange + black) and sphere of influence neighbours (black); apparent holes appear for sphere of influence neighbours where an urban municipality is surrounded by a dominant rural municipality (see @fig-plotpolpres15)"
#| code-fold: true
#| out.width: 100%
opar <- par(mar = c(0,0,0,0)+0.5)
pol_pres15 |> 
    st_geometry() |> 
    plot(border = "grey", lwd = 0.5)
nb_soi |> plot(coords = coords, add = TRUE, 
			   points = FALSE, lwd = 0.5)
nb_tri |> 
    diffnb(nb_soi) |> 
    plot(coords = coords, col = "orange", add = TRUE,
		 points = FALSE, lwd = 0.5)
par(opar)



coords |> 
    knearneigh(k = 1) |> 
    knn2nb() |> 
    nbdists(coords) |> 
    unlist() |> 
    summary()



coords |> dnearneigh(0, 18000) -> nb_d18



coords |> dnearneigh(0, 18000, use_kd_tree = FALSE) -> nb_d18a



all.equal(nb_d18, nb_d18a, check.attributes = FALSE)


nb_d18



(nb_d18 |> n.comp.nb() -> n_comp)$nc


table(n_comp$comp.id)



(coords |> dnearneigh(0, 18300) -> nb_d183)


(nb_d183 |> n.comp.nb())$nc



(coords |> dnearneigh(0, 16000) -> nb_d16)



((coords |> knearneigh(k = 6) -> knn_k6) |> knn2nb() -> nb_k6)



(knn_k6 |> knn2nb(sym = TRUE) -> nb_k6s)



(nb_k6s |> n.comp.nb())$nc



old_use_s2 <- sf_use_s2()



sf_use_s2(TRUE)



pol_pres15_ll |> 
    st_geometry() |> 
    st_centroid(of_largest_polygon = TRUE) -> coords_ll



(coords_ll |> dnearneigh(0, 18.3, use_s2 = TRUE, 
						 dwithin = TRUE) -> nb_d183_ll)



isTRUE(all.equal(nb_d183, nb_d183_ll, check.attributes = FALSE))



(coords_ll |> dnearneigh(0, 18.3) -> nb_d183_llce)



isTRUE(all.equal(nb_d183_llce, nb_d183_ll,
				 check.attributes = FALSE))



(coords_ll |> knearneigh(k = 6) |> knn2nb() -> nb_k6_ll)



isTRUE(all.equal(nb_k6, nb_k6_ll, check.attributes = FALSE))



nb_q |> nbdists(coords_ll) |> unlist() |> summary()



nb_q |> nbdists(coords) |> unlist() |> summary()


sf_use_s2(old_use_s2)



## args(nb2listw)


xargs(nb2listw)



## args(spweights.constants)


xargs(spweights.constants)



(nb_q |> 
    nb2listw(style = "B") -> lw_q_B) |> 
    spweights.constants() |> 
    data.frame() |> 
    subset(select = c(n, S0, S1, S2))



(nb_q |> 
        nb2listw(style = "W") -> lw_q_W) |> 
    spweights.constants() |> 
    data.frame() |> 
    subset(select = c(n, S0, S1, S2))



nb_d183 |> 
    nbdists(coords) |> 
    lapply(function(x) 1/(x/1000)) -> gwts
(nb_d183 |> nb2listw(glist=gwts, style="B") -> lw_d183_idw_B) |> 
    spweights.constants() |> 
    data.frame() |> 
    subset(select=c(n, S0, S1, S2))



try(nb_d16 |> nb2listw(style="B") -> lw_d16_B)



nb_d16 |> 
    nb2listw(style="B", zero.policy=TRUE) |> 
    spweights.constants(zero.policy=TRUE) |> 
    data.frame() |> 
    subset(select=c(n, S0, S1, S2))



nb_q



(nb_q |> nblag(2) -> nb_q2)[[2]]



nblag_cumul(nb_q2)



union.nb(nb_q2[[2]], nb_q2[[1]])



diameter(g1)



g1 |> shortest.paths() -> sps
(sps |> apply(2, max) -> spmax) |> max()



mr <- which.max(spmax)
pol_pres15$name0[mr]


#| out.width: 100%
#| code-fold: true
#| fig.cap: "Relationship of shortest paths to distance for Lutowiska; left panel: shortest path counts from Lutowiska; right panel: plot of shortest paths from Lutowiska to other observations, and distances from Lutowiska to other observations"
pol_pres15$sps1 <- sps[,mr]
tm1 <- tm_shape(pol_pres15) +
          tm_fill("sps1", title = "Shortest path\ncount")
coords[mr] |> 
    st_distance(coords) |> 
    c() |> 
    (function(x) x/1000)() |> 
    units::set_units(NULL) -> pol_pres15$dist_52
library(ggplot2)
g1 <- ggplot(pol_pres15, aes(x = sps1, y = dist_52)) +
		geom_point() +
		xlab("Shortest path count") +
		ylab("km distance")
gridExtra::grid.arrange(tmap_grob(tm1), g1, nrow=1)



save(list = ls(), file = "ch14.RData")

