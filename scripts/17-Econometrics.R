
eval_inla = Sys.getenv("EVAL_INLA") != "false"


knitr::opts_chunk$set(echo = TRUE, paged.print = FALSE)


load("ch16.RData")



library(spatialreg)
eigs_489 <- eigenw(lw_q_489)
SDEM_489 <- errorsarlm(form, data = boston_489, 
      listw = lw_q_489, Durbin = TRUE, zero.policy = TRUE,
      control = list(pre_eig = eigs_489))
SEM_489 <- errorsarlm(form, data = boston_489, 
      listw = lw_q_489, zero.policy = TRUE,
      control = list(pre_eig = eigs_489))



cbind(data.frame(model=c("SEM", "SDEM")), 
      rbind(broom::tidy(Hausman.test(SEM_489)), 
            broom::tidy(Hausman.test(SDEM_489))))[,1:4]



eigs_94 <- eigenw(lw_q_94)
SDEM_94 <- errorsarlm(form, data=boston_94, listw=lw_q_94,
					  Durbin = TRUE,
					  control = list(pre_eig=eigs_94))
SEM_94 <- errorsarlm(form, data = boston_94, listw = lw_q_94,
					 control = list(pre_eig = eigs_94))



cbind(data.frame(model=c("SEM", "SDEM")), 
      rbind(broom::tidy(Hausman.test(SEM_94)), 
            broom::tidy(Hausman.test(SDEM_94))))[, 1:4]



cbind(data.frame(model=c("SEM", "SDEM")),
	  rbind(broom::tidy(LR1.Sarlm(SEM_94)),
			broom::tidy(LR1.Sarlm(SDEM_94))))[,c(1, 4:6)]



o <- lmtest::lrtest(SEM_489, SDEM_489)
attr(o, "heading")[2] <- "Model 1: SEM_489\nModel 2: SDEM_489"
o


o <- lmtest::lrtest(SEM_94, SDEM_94)
attr(o, "heading")[2] <- "Model 1: SEM_94\nModel 2: SDEM_94"
o



SLX_489 <- lmSLX(form, data = boston_489, listw = lw_q_489,
				 zero.policy = TRUE)
o <- lmtest::lrtest(SLX_489, SDEM_489)
attr(o, "heading")[2] <- "Model 1: SLX_489\nModel 2: SDEM_489"
o



SLX_94 <- lmSLX(form, data = boston_94, listw = lw_q_94)
o <- lmtest::lrtest(SLX_94, SDEM_94)
attr(o, "heading")[2] <- "Model 1: SLX_94\nModel 2: SDEM_94"
o



SLX_489w <- lmSLX(form, data = boston_489, listw = lw_q_489,
				  weights = units, zero.policy = TRUE)
SDEM_489w <- errorsarlm(form, data = boston_489,
						listw = lw_q_489, Durbin = TRUE,
						weights = units, zero.policy = TRUE,
						control = list(pre_eig = eigs_489))
o <- lmtest::lrtest(SLX_489w, SDEM_489w)
attr(o, "heading")[2] <- "Model 1: SLX_489w\nModel 2: SDEM_489w"
o



SLX_94w <- lmSLX(form, data = boston_94, listw = lw_q_94,
				 weights = units)
SDEM_94w <- errorsarlm(form, data = boston_94, listw = lw_q_94,
					   Durbin = TRUE, weights = units,
                       control = list(pre_eig = eigs_94))
o <- lmtest::lrtest(SLX_94w, SDEM_94w)
attr(o, "heading")[2] <- "Model 1: SLX_94w\nModel 2: SDEM_94w"
o



sum_imp_94_SDEM <- summary(impacts(SDEM_94))
rbind(Impacts = sum_imp_94_SDEM$mat[5,], 
	  SE = sum_imp_94_SDEM$semat[5,])



sum_imp_94_SLX <- summary(impacts(SLX_94))
rbind(Impacts = sum_imp_94_SLX$mat[5,], 
	  SE = sum_imp_94_SLX$semat[5,])



sum_imp_94_SDEMw <- summary(impacts(SDEM_94w))
rbind(Impacts = sum_imp_94_SDEMw$mat[5,], 
	  SE = sum_imp_94_SDEMw$semat[5,])



sum_imp_94_SLXw <- summary(impacts(SLX_94w))
rbind(Impacts = sum_imp_94_SLXw$mat[5,], 
	  SE = sum_imp_94_SLXw$semat[5,])



nd <- boston_506[is.na(boston_506$median),]
t0 <- exp(predict(SDEM_489, newdata = nd, listw = lw_q,
				  pred.type = "TS", zero.policy  =TRUE))
suppressWarnings(t1 <- exp(predict(SDEM_489, newdata = nd,
									listw = lw_q,
									pred.type = "KP2",
                                    zero.policy = TRUE)))
suppressWarnings(t2 <- exp(predict(SDEM_489, newdata = nd,
									listw = lw_q,
									pred.type = "KP5",
                                    zero.policy = TRUE)))

