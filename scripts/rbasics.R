
## a |> b() |> c() |> d(n = 10)



## d(c(b(a)), n = 10)



## tmp1 <- b(a)
## tmp2 <- c(tmp1)
## tmp3 <- d(tmp2, n = 10)



typeof(1:10)
length(1:10)
typeof(1.0)
length(1.0)
typeof(c("foo", "bar"))
length(c("foo", "bar"))
typeof(c(TRUE, FALSE))



i <- integer(0)
typeof(i)
i
length(i)



a <- c(1,2,3)
a[2]
a[[2]]
a[2:3]
a[2:3] <- c(5,6)
a
a[[3]] <- 10
a



l <- list(3, TRUE, "foo")
typeof(l)
length(l)



l[1]
l[[1]]



l[1:2] <- list(4, FALSE)
l
l[[3]] <- "bar"
l



l <- list(first = 3, second = TRUE, third = "foo")
l



l$second
l$second <- FALSE
l



3 == NULL # not FALSE!
NULL == NULL # not even TRUE!



is.null(NULL)



l <- l[c(1,3)] # remove second, implicitly
l



l$second <- NULL
l



a <- 1:3
attr(a, "some_meta_data") = "foo"
a



attr(a, "some_meta_data")
attr(a, "some_meta_data") <- "bar"
attr(a, "some_meta_data")



attributes(a)
attributes(a) = list(some_meta_data = "foo")
attributes(a)



class(1:3)
class(c(TRUE, FALSE))
class(c("TRUE", "FALSE"))



a <- 1:3
class(a) <- "foo"
a
class(a)
attributes(a)



print.foo <- function(x, ...) { 
	print(paste("an object of class foo with length", length(x)))
}
print(a)



unclass(a)



library(sf) |> suppressPackageStartupMessages()
p <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,0))))
p



unclass(p)



a <- 1:8
class(a)
attr(a, "dim") <- c(2,4) # or: dim(a) = c(2,4)
class(a)
a
attr(a, "dim") <- c(2,2,2) # or: dim(a) = c(2,2,2)
class(a)
a



a <- c(first = 3, second = 4, last = 5)
a["second"]
attributes(a)



a <- matrix(1:4, 2, 2)
dimnames(a) <- list(rows = c("row1", "row2"),
				    cols = c("col1", "col2"))
a
attributes(a)



df <- data.frame(a = 1:3, b = c(TRUE, FALSE, TRUE))
attributes(df)



## f <- function(x) {
##    a <- create_obj(x) # call some other function
##    attributes(a) <- list(class = "foo", meta = 33)
##    a
## }



## f <- function(x) {
##    a <- create_obj(x) # call some other function
##    structure(a, class = "foo", meta = 33)
## }



system.file("gpkg/nc.gpkg", package = "sf") %>%
	read_sf() -> nc



attributes(nc)



nc$geom



attributes(nc$geom)



nc$geom[[4]] |> format(width = 60, digits = 5)



typeof(nc$geom[[4]])



attributes(nc$geom[[4]])



length(nc$geom[[4]])



lengths(nc$geom)



typeof(nc$geom[[4]])



typeof(nc$geom[[4]][[1]])



length(nc$geom[[4]][[1]])



typeof(nc$geom[[4]][[1]][[1]])
dim(nc$geom[[4]][[1]][[1]])
head(nc$geom[[4]][[1]][[1]])



nc$geom[[4]][[1]][[1]][3,2] <- 36.5

