library(cobs)
data(DublinWind)
attach(DublinWind)##-> speed & day (instead of "wind.x" & "DUB.")

if(!interactive()) postscript("wind.ps", horizontal = TRUE)

stopifnot(identical(day,c(rep(c(rep(1:365,3),1:366),4),
                          rep(1:365,2))))
nknots <- 13
## Compute the quadratic median smoothing B-spline using SIC
## lambda selection
co.o50 <-
 cobs(day,speed,knots.add = TRUE,constraint="periodic", nknots=nknots,
      tau = .5,lambda = -1,factor = 3,method = "uniform")
summary(co.o50)
## Since SIC chooses a lambda that corresponds to the smoothest
## possible fit, rerun cobs with a larger lstart value
(lstart <- log(.Machine$double.xmax)^3) # 3.57 e9

co.o5 <-
    cobs(day,speed,knots.add = TRUE,constraint = "periodic", nknots = nknots,
         tau = .5,lambda = -1,factor = 3,method = "uniform",lstart = lstart)
summary(co.o5)
(pc.5 <- predict(co.o5, interval = "both"))

co.o9 <- ## Compute the .9 quantile smoothing B-spline
    cobs(day,speed,knots.add = TRUE,constraint = "periodic",nknots = nknots,
         tau = .9,lambda = -1,factor = 3,method = "uniform")
summary(co.o9)
(pc.9 <- predict(co.o9,interval = "both"))

co.o1 <- ## Compute the .1 quantile smoothing B-spline
    cobs(day,speed,knots.add = TRUE,constraint = "periodic",nknots = nknots,
         tau = .1,lambda = -1,factor = 3,method = "uniform")
summary(co.o1)
(pc.1 <- predict(co.o1, interval = "both"))

par(mfrow = c(1,2))
plot(day,speed, pch = 3, cex=0.6,  xlab = "DAYS", ylab = "SPEED (knots)")
lines(pc.5, lwd = 2.5, col = 2)
lines(pc.9, lwd = 1.5, col = 3)
lines(pc.1, lwd = 1.5, col = 3)
plot(day,speed,type = "n",xlab = "DAYS", ylab = "SPEED (knots)")
lines(pc.5, lwd = 1.5)
lines(pc.9, col = 3)
lines(pc.1, col = 3)
abline(v = co.o5$knots, lty = 3, col = "gray70")
## rather rug(co.o5$knots, lty = 2)
