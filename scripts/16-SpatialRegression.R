
eval_inla = Sys.getenv("EVAL_INLA") != "false"



knitr::opts_chunk$set(echo = TRUE, paged.print = FALSE)



library(sf)
library(spData)
boston_506 <- st_read(system.file("shapes/boston_tracts.shp",
					  package = "spData")[1], quiet = TRUE)


nb_q <- spdep::poly2nb(boston_506)
lw_q <- spdep::nb2listw(nb_q, style = "W")



table(boston_506$censored)


summary(boston_506$median)



boston_506$CHAS <- as.factor(boston_506$CHAS)
boston_489 <- boston_506[!is.na(boston_506$median),]
nb_q_489 <- spdep::poly2nb(boston_489)
lw_q_489 <- spdep::nb2listw(nb_q_489, style = "W",
							zero.policy = TRUE)



agg_96 <- list(as.character(boston_506$NOX_ID))
boston_96 <- aggregate(boston_506[, "NOX_ID"], by = agg_96,
					   unique)
nb_q_96 <- spdep::poly2nb(boston_96)
lw_q_96 <- spdep::nb2listw(nb_q_96)
boston_96$NOX <- aggregate(boston_506$NOX, agg_96, mean)$x
boston_96$CHAS <-
    aggregate(as.integer(boston_506$CHAS)-1, agg_96, max)$x



nms <- names(boston_506)
ccounts <- 23:31
for (nm in nms[c(22, ccounts, 36)]) {
  boston_96[[nm]] <- aggregate(boston_506[[nm]], agg_96, sum)$x
}
br2 <- 
  c(3.50, 6.25, 8.75, 12.5, 17.5, 22.5, 30, 42.5, 60) * 1000
counts <- as.data.frame(boston_96)[, nms[ccounts]]
f <- function(x) matrixStats::weightedMedian(x = br2, w = x,
									 interpolate = TRUE)
boston_96$median <- apply(counts, 1, f)
is.na(boston_96$median) <- boston_96$median > 50000
summary(boston_96$median)



POP <- boston_506$POP
f <- function(x) matrixStats::weightedMean(x[,1], x[,2])
for (nm in nms[c(9:11, 14:19, 21, 33)]) {
  s0 <- split(data.frame(boston_506[[nm]], POP), agg_96)
  boston_96[[nm]] <- sapply(s0, f)
}
boston_94 <- boston_96[!is.na(boston_96$median),]
nb_q_94 <- spdep::subset.nb(nb_q_96, !is.na(boston_96$median))
lw_q_94 <- spdep::nb2listw(nb_q_94, style="W")



boston_94a <- aggregate(boston_489[,"NOX_ID"], 
						list(boston_489$NOX_ID), unique)
nb_q_94a <- spdep::poly2nb(boston_94a)
NOX_ID_no_neighs <-
		boston_94a$NOX_ID[which(spdep::card(nb_q_94a) == 0)]
boston_487 <- boston_489[is.na(match(boston_489$NOX_ID,
									 NOX_ID_no_neighs)),]
boston_93 <- aggregate(boston_487[, "NOX_ID"],
					   list(ids = boston_487$NOX_ID), unique)
row.names(boston_93) <- as.character(boston_93$NOX_ID)
nb_q_93 <- spdep::poly2nb(boston_93,
        row.names = unique(as.character(boston_93$NOX_ID)))



form <- formula(log(median) ~ CRIM + ZN + INDUS + CHAS + 
				I((NOX*10)^2) + I(RM^2) + AGE + log(DIS) +
				log(RAD) + TAX + PTRATIO + I(BB/100) + 
				log(I(LSTAT/100)))



library(Matrix)
library(lme4)
MLM <- lmer(update(form, . ~ . + (1 | NOX_ID)), 
			data = boston_487, REML = FALSE)



boston_93$MLM_re <- ranef(MLM)[[1]][,1]



library(hglm) |> suppressPackageStartupMessages()
suppressWarnings(HGLM_iid <- hglm(fixed = form,
								  random = ~1 | NOX_ID,
								  data = boston_487,
								  family = gaussian()))
boston_93$HGLM_re <- unname(HGLM_iid$ranef)



library(spatialreg)
W <- as(spdep::nb2listw(nb_q_93, style = "B"), "CsparseMatrix")



suppressWarnings(HGLM_car <- hglm(fixed = form,
								  random = ~ 1 | NOX_ID,
								  data = boston_487,
								  family = gaussian(),
								  rand.family = CAR(D=W)))
boston_93$HGLM_ss <- HGLM_car$ranef[,1]



library(R2BayesX) |> suppressPackageStartupMessages()



BX_iid <- bayesx(update(form, . ~ . + sx(NOX_ID, bs = "re")),
				 family = "gaussian", data = boston_487,
				 method = "MCMC", iterations = 12000,
				 burnin = 2000, step = 2, seed = 123)



boston_93$BX_re <- BX_iid$effects["sx(NOX_ID):re"][[1]]$Mean



RBX_gra <- nb2gra(nb_q_93)
all.equal(row.names(RBX_gra), attr(nb_q_93, "region.id"))



all.equal(unname(diag(RBX_gra)), spdep::card(nb_q_93))



BX_mrf <- bayesx(update(form, . ~ . + sx(NOX_ID, bs = "mrf",
										 map = RBX_gra)), 
                 family = "gaussian", data = boston_487,
				 method = "MCMC", iterations = 12000,
				 burnin = 2000, step = 2, seed = 123)



boston_93$BX_ss <- BX_mrf$effects["sx(NOX_ID):mrf"][[1]]$Mean



library(INLA) |> suppressPackageStartupMessages()







ID2 <- as.integer(as.factor(boston_487$NOX_ID))







M <- Diagonal(nrow(W), rowSums(W)) - W
Cmatrix <- Diagonal(nrow(M), 1) -  M







library(mgcv)
names(nb_q_93) <- attr(nb_q_93, "region.id")
boston_487$NOX_ID <- as.factor(boston_487$NOX_ID)



GAM_MRF <- gam(update(form, . ~ . + s(NOX_ID, bs = "mrf",
									  xt = list(nb = nb_q_93))),
			   data = boston_487, method = "REML")



ssre <- predict(GAM_MRF, type = "terms", 
				se = FALSE)[, "s(NOX_ID)"]
all(sapply(tapply(ssre, list(boston_487$NOX_ID), c),
		   function(x) length(unique(round(x, 8))) == 1))



boston_93$GAM_ss <- aggregate(ssre, list(boston_487$NOX_ID), 
							  head, n=1)$x







save(list = ls(), file = "ch16.RData")

