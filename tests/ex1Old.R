#### OOps! Running this in 'CMD check' or in *R* __for the first time__
#### ===== gives a wrong result (at the end) than when run a 2nd time
####-- problem disappears with introduction of   if (psw) call ... in Fortran

library(cobs)
options(digits = 6)
postscript("ex1Old.ps")

cobs <- cobsOld
###############

summaryCobs <- function(x, ...)
{
    ## Purpose: something like print(summary( cobs.result ))
    ## ----------------------------------------------------------------------
    ## Arguments: x: result of cobs()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 15 Feb 2002, 16:55
    str(x, ...)
    print(as.data.frame(x[c("cb.lo", "ci.lo", "fit", "ci.up", "cb.up")]), ...)
    cat("knots :\n"); print(x$knots, ...)
    cat("coef  :\n"); print(x$coef, ...)
    if(!is.null(x$sic)) {
        print(cbind(lambda = x$pp.lambda, SIC = x$sic), ...)
    }
}


## Simple example from  example(cobs)
set.seed(908)
x <- seq(-1,1, len = 50)
f.true <- pnorm(2*x)
y <- f.true + rnorm(50)/10
## specify constraints (boundary conditions)
con <- rbind(c( 1,min(x),0),
             c(-1,max(x),1),
             c( 0, 0,  0.5))
## obtain the median regression B-spline using automatically selected knots
coR <- cobs(x,y,constraint = "increase", pointwise = con)
summaryCobs(coR)

## compute the median smoothing B-spline using automatically chosen lambda
coS <- cobs(x,y,constraint = "increase", pointwise = con,
            lambda = -1, lstart = 7872)
summaryCobs(coS)

##-- real data example (still n = 50)
data(cars)
attach(cars)
summaryCobs(co1   <- cobs(speed, dist, "increase"))
summaryCobs(co1.1 <- cobs(speed, dist, "increase", knots.add = TRUE))
1 - sum(co1 $ resid ^2) / sum((dist - mean(dist))^2) # R^2 = 64.2%

summaryCobs(co2 <- cobs(speed, dist, "increase", lambda = -1, lstart = 7872))# warning
1 - sum(co2 $ resid ^2) / sum((dist - mean(dist))^2) # R^2 = 65.8%

summaryCobs(co3 <- cobs(speed, dist, "convex", lambda = -1, lstart = 7872))# warning
1 - sum(co3 $ resid ^2) / sum((dist - mean(dist))^2) # R^2 = 65.2%

detach(cars)

##-- another larger example using "random" x
x <- round(sort(rnorm(500)), 3) # rounding -> multiple values
sum(duplicated(x)) # 32
y <- (fx <- exp(-x)) + rt(500,4)/4
summaryCobs(cxy  <- cobs(x,y, "decrease"))
1 - sum(cxy $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 95.9%

## Interpolation
if(FALSE) { ##-- since it takes too long here!
   system.time(cxyI  <- cobs(x,y, "decrease", knots = unique(x)))[1]
   ## takes very long : 1864.46 sec. (Pent. III, 700 MHz)
   summaryCobs(cxyI)# only 8 knots remaining!
}
ux <- unique(x)
xx <- seq(ux[1], x[length(ux)], len = 201)
system.time(cxyI  <- cobs(x,y, "decrease", knots = ux, nknots = length(ux)))[1]
## 17.3 sec. (Pent. III, 700 MHz)
summaryCobs(cxyI)
1 - sum(cxyI $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 96.56%
plot(x,y, cex = 3/4, main = "Artificial (x,y), N=500 : `interpolating' cobs()")
lines(xx, exp(-xx), type = "l", col = "gray40")
lines(cxyI$z, cxyI$fit, col = "red")

system.time(cxyS <- cobs(x,y, "decrease", lambda = -1, lstart = 7872))
## somewhat <  2 sec. (Pent. III, 700 MHz)
summaryCobs(cxyS)
R2 <- 1 - sum(cxyS $ resid ^ 2) / sum((y - mean(y))^2)
R2 # R^2 = 96.3%

plot(x,y, cex = 3/4, main = "Artificial (x,y), N=500 : cobs(*, lambda = -1)")
mtext(substitute(R^2 == r2 * "%", list(r2 = round(100*R2,1))))
lines(xx, exp(-xx), type = "l", col = "gray40")
lines(cxyS$z, cxyS$fit, col = "red")

## Show print-monitoring :

cxyS <- cobs(x,y, "decrease", lambda = -1, trace = 2, lstart = 7872)
cxyS <- cobs(x,y, "none",     lambda = -1, trace = 3, lstart = 7872)

## this does NOT converge (and "trace = 3" does *not* show it -- improve!)

cxyC <- cobs(x,y, "concave", lambda = -1, lstart = 7872)
summaryCobs(cxyC)


dev.off()
