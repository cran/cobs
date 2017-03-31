#### -*- mode: R; kept-new-versions: 25; kept-old-versions: 20 -*-
####
#### MM: port of Lutz Duembgen's matlab code to R

## concaveRegr <-
## concaveRegression <-
## 'conreg' is analogue to the standard R  'isoreg()'
## also, it now can be convex or concave regression
conreg <- function(x, y = NULL, w = NULL, convex = FALSE,
		   method = c("Duembgen06_R", "SR"),
		   tol = c(1e-10, 1e-7),
                   maxit = c(500, 20),
		   adjTol = TRUE, verbose = FALSE)
{
  ## Ported to R and enhanced:
  ## - work for unordered, even duplicated 'x'
  ## - made slightly faster; more options;
  ## - detect infinite loop; auto-adjust tol;
  ## - changed tol (prec) to 1e-7
  ## - new arg.	 convex == FALSE   <==>	 concave == TRUE
  ## - define "class" and many methods, plot, predict, ...

  ## Martin Maechler, 27.-28. Apr 2007

  ##  2)	use 'call'  in plot, etc...

  ##  4) ------> refind the "shape preserving splines" from numerical analysis!
  ##     and do these as a last step

  ##  5) I have upped tol to 1e-7; and added an 'adjTol' trick.
  ##   Can anyone find an example where  tol = 1e-7 is "too small"
  ##   insofar as it takes less steps than  tol = 1e-10 erronously ?

  ## berechnet zu gegebenen Vektoren x, y, w im R^n mit
  ##	x(1) < x(2) < ... < x(n)
  ## einen Spaltenvektor M im R^n mit minimaler Quadratsumme
  ##	sum(w .* (M - y)^2 )
  ## unter der Nebenbedingung
  ##	     (M[i] - M[i-1]) / (x[i] - x[i-1]) >=
  ##	     (M[i+1] - M[i]) / (x[i+1] - x[i])
  ##	fuer 1 <= i <= n ,
  ## wobei M[0] := M[n+1] := - inf.
  ##
  ## Default fuer w ist w[i] = 1.
  ##
  ## Die Berechnungsmethode ist ein Active-Set-Verfahren,
  ## TODO {this was 'matlab'}:
  ## --- dessen Zwischenschritte auch graphisch illustriert werden.
  ##
  ## iKnots gibt jene Indizes  i  an, wo Vektor M die strikte
  ## Ungleichung
  ##	 (M[i] - M[i-1]) / (x[i] - x[i-1]) >
  ##	 (M[i+1] - M[i]) / (x[i+1] - x[i])
  ## erfuellt.
  ##
  ## Lutz Duembgen, 3. Juli 2006


  ## --- use the same code as in smooth.spline() {!}
  xy <- xy.coords(x, y)
  y <- xy$y
  x <- xy$x
  n <- length(x)
  method <- match.arg(method)
  isSR <- method == "SR"
  if(isSR) {
    if(!is.null(w))
      stop("weights are not yet supported for method \"SR\"")
    ## cobs() has all the weights; handling of ties, etc etc
    ## Here, we only sort (if needed)
    if(doSort <- is.unsorted(x)) {
      i. <- sort.list(x, method="quick")
      x. <- x[i.]
      y. <- y[i.]
    } else {
      x. <- x
      y. <- y
    }
    w. <- NULL
    nx <- n

  } else { ## (method == "Duembgen06_R")
    w <-
      if(is.null(w)) rep(1, n)
      else {
        if(n != length(w)) stop("lengths of 'x' and 'w' must match")
        if(any(w < 0)) stop("all weights should be non-negative")
        if(all(w == 0)) stop("some weights should be positive")
        (w * sum(w > 0))/sum(w)
      }# now sum(w) == #{obs. with weight > 0} == sum(w > 0)

    ## Replace y[] (and w[]) for same x[] (to 6 digits precision) by their mean :
    x <- signif(x, 6)
    x. <- unique(sort(x))
    nx <- length(x.)
    if(nx == 0) stop("need at least one x value")
    ## FIXME: for large n, the following is *very* slow
    ox <- match(x, x.)
    ## Faster, much simplified version of tapply()
    tapply1 <- function (X, INDEX, FUN = NULL, ..., simplify = TRUE) {
      sapply(unname(split(X, INDEX)), FUN, ...,
             simplify = simplify, USE.NAMES = FALSE)
    }
    tmp <- matrix(unlist(tapply1(seq_along(y), ox,
                                 function(i)
                                   c(sum(w[i]), sum(w[i]*y[i]))#,sum(w[i]*y[i]^2))
                                 ), use.names=FALSE),
                  ncol = 2, #3,
                  byrow=TRUE)
    w. <- tmp[, 1]
    y. <- tmp[, 2]/ifelse(w. > 0, w., 1)
    ## yssw <- sum(tmp[, 3] - w.*y.^2) # will be added to RSS for GCV
  }

  ## sign switch: different for the two methods ! (?!)
  convSwitch <- (isSR && !convex) || (!isSR && convex)

  if(convSwitch) y. <- - y. # and will revert at the end

  stopifnot(length(maxit) >= 1, maxit == round(maxit))
  if(length(maxit) == 1) maxit <- rep.int(maxit, 2)
  stopifnot(length(tol) >= 1, tol >= 0)
  if(length(tol) == 1) tol <- rep.int(tol, 2)

  ## rescale 'tol' to be x-scale equivariant:
  if(nx > 1) tol <- tol / sd(x.)

  if(method == "SR") {

    r <- .C(SR_R, # --> ../src/SR.c
            n=  as.integer(n),
            cc= as.double(x.[n] - x.[1L]), # x-scale: diff(range(x.))
            m1= integer(1),
            ind= integer(n+1),
            x = as.double(x.),
            y = as.double(y.),
            r = double(n),
            ## MM: what are these? -- surely some are interesting
            ## ---> ../convexreg-PietGr/convexregres/main.c  outputs (x[],Y[]),  (x, H) and (x, D)
            ##   but the
            ##      ../convexreg-PietGr/convexregresnew/main.c  no longer does ..
            ##  (if not, create and free in C !!)
            R = double(n),
            H = double(n),
            S = double(n),
            Y = double(n),
            D = double(n),
            tol = as.double(tol),
            maxit = as.integer(maxit),
            verbose = as.integer(verbose),
            phiBL = double(1),
            numIt = integer(1))

    ## the vector r contains the convex LS estimates f[i] at the points x[i]
    ## the data are in the vector y and phi is the criterion function
    ## the vector of indices ind contains the indices of the points where the
    ## solution has a kink (only for method 1) and m is the number of kinks

    M <- r$ r
    H <- r$ H
    ## TODO ?  'D' , 'S', ... ???

    if(doSort) {
      ## FIXME -- revert or ?????
    }
    iK <- r$ind[1L + seq_len(r$m1)]
    iter <- r$numIt # FIXME? c(iter, innerIt)
  }
  else if(method == "Duembgen06_R") {
    ## auxiliaries in locEstimate
    i.n <- seq_len(nx)
    rtW <- sqrt(w.)
    rWy <- rtW * y.

    ## rather work with knot *indices* instead of logical (or 0/1):
    iK <- if(nx > 1) c(1,nx) else 1

    ## M = "current y.hat" :
    M  <- locEstimate(x., rWy, rtW, iK)
    wR <- (M-y.)* w.
    H  <- locDirDeriv(x., wR, iK)

    ## Precision parameter: NB  length(tol) == 2
    prec <- tol * mean(abs(wR))
    eConv <- prec[2] # should maybe be different
    prec <- prec[1]
    ## could also use "adaptive" prec:
    ## while (max(H) > (prec <- tol * mean(abs(wR),trim=.1)) && iter < maxit) {

    iter <- innerIt <- 0
    while (max(H) > prec && iter < maxit[1]) {
      ## extend the set of knots - at  max_H location:
      ## isKnot[(im <- which.max(H))] <- TRUE
      iK <- sort(c(iK, im <- which.max(H)))
      if(verbose) cat("#{knots}=",length(iK),"; new knot at", im, "")
      ## compute new candidate:
      M_new <- locEstimate(x., rWy, rtW, iK)

      ## check new candidate's concavity; local convexities, desired <= 0
      Conv	<- locConvexities(x., M)
      Conv_new	<- locConvexities(x., M_new)
      iit <- 0
      while (max(Conv_new) > eConv && iit < maxit[2]) {
        if(verbose) cat(".")
        ## modify M_new and (typically!) reduce the set of knots:
        JJ <- (Conv_new > eConv) # non-empty
        t <- min(- Conv[JJ] / (Conv_new[JJ] - Conv[JJ]))
        ## if(adjTol) ## FIXME: can restore 'M' from 'oM' -- or return 'oM' ??    oM <- M
        M <- (1 - t)*M + t*M_new
        Conv <- locConvexities(x.,M)
        oiK <- iK
        iK <- c(1, i.n[Conv < -eConv], nx) ## = locKnInd(Conv,eConv)
        if(adjTol && identical(iK, oiK)) { ## may not converge ..
          if(verbose)
            cat("inner knot set unchanged; increasing eConv\n\t")
          eConv <- 8 * eConv
          next # while
        }

        M_new	 <- locEstimate(x., rWy, rtW, iK)
        Conv_new  <- locConvexities(x.,M_new)
        iit <- iit + 1
      }
      innerIt <- innerIt + iit
      if(verbose) cat("\n")

      if(iit == maxit[2])
        warning("** inner iterations did *not* converge in ",
                maxit[2], " steps")
      M <- M_new
      Conv <- Conv_new
      iK <- c(1L, i.n[Conv < -eConv], nx) ## = locKnInd(Conv,eConv)
      wR <- (M-y.)* w.
      H	 <- locDirDeriv(x., wR, iK)
      iter <- iter+1
    } ## while (outer iter)

    iter <- c(iter, innerIt)

  } ## method "Duembgen06"


  if(iter[1] == maxit[1])
    warning("*** iterations did *not* converge in ", maxit[1], " steps")
  else if(method != "Duembgen06" || !iter) {
    ## did not need any iteration;  e.g. for n = 2
    Conv <- locConvexities(x., M)
  }

  if(convSwitch) {## switch the signs back
    y. <- -y.
    M <- -M
    H <- -H
    Conv <- -Conv
  }
  ## return
  structure(list(x = x., y = y., w = w., yf = M,
		 convex = convex, call = match.call(),
		 iKnots = iK, deriv.loc = H, conv.loc = Conv,
		 iter = iter),
	    class = "conreg")
}

## auxiliary routines: Note that we could gain speed if
## ------------------
##  1) these were *local* to the main function
##  2) they would use *updating* {convexities; DD matrix, even qr(rtW * DD) !}
##

locEstimate <- function(x, rWy, rtW, iK)
{
    n  <- length(x)
    ## Note that iK == JJ <- which(isKnot)
    p  <- length(iK)
    stopifnot(p >= 1, n >= 1)
    DD	<- matrix(0, n, p)
    DD[1,1] <- 1
    for (i in seq_len(p-1)) {
	if(iK[i]+1 <= iK[i+1]) { ## <- safety check
	    kn.i <- x[iK[i]]
	    ii <- (iK[i]+1):iK[i+1]
	    d2 <- (x[ii] - kn.i) / (x[iK[i+1]] - kn.i)
	    DD[ii, i:(i+1)] <- c(1 - d2, d2)
	}
	else ## does it happen? [--> message]  MM has never seen it - proof?
	    message("locEst: empty between knots ", iK[i]," and ", iK[i+1])
    }
    ## Solve the weighted L.S. problem
    ##	   beta = argmin_b[ sum_i w_i (y_i - (DD %*% b)_i)^2]

    ## qr.coef(.) here equivalent to solve()  {and slightly faster}
    as.vector(DD %*% qr.coef(qr(rtW * DD), rWy))
}


locDirDeriv <- function(x, wRes, iK) {
  n <- length(x)
  H  <- numeric(n)
  ## Note that iK == JJ <- which(isKnot)
  ## p  <- length(iK)
  cs <- cumsum(wRes)
  for (i in seq_len(length(iK)-1)) {
    if((ik1 <- iK[i]) <= (ik2 <- iK[i+1L] - 2L)) {
      ii <- ik1:ik2
      H[ii+1L] <- cumsum((x[ii+1L] - x[ii]) * cs[ii])
    }
    ## else ## does happen occasionally
    ## message("DirDeriv: empty between knots ", iK[i]," and ", iK[i+1])
  }
  H
}


locConvexities <- function(x,y) {
  n <- length(x)
  if(n > 2) {
      tmp  <- (y[-1L] - y[-n]) / (x[-1L] - x[-n])
      c(0, tmp[-1L] - tmp[-(n-1L)], 0)
  } else numeric(n)
}



if(FALSE) # now unused
locKnots <- function(conv, prec) {
    n <- length(conv)
    isKnot <- c(TRUE, rep(FALSE, n-2), TRUE)
    isKnot[conv < -prec] <- TRUE
    isKnot
}

if(FALSE) # now unused
locKnInd <- function(conv, prec) {
    n <- length(conv)
    c(1, seq_len(n)[conv < -prec], n)
}



###-------------- Methods for the result object/class : ----------------------

print.conreg <-
    function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall: ", deparse(x$call), "\n")
    n <- length(x$x)
    nK <- length(iK <- x$iKnots)
    cat(if(x$convex) "Convex" else "Concave",
	"regression: From", n, "separated x-values,",
	"using", nK - 2, "inner knots,\n  ")
    ## FIXME:  mention "original n"
    r <- resid(x)
    RSS <- sum(r^2)
    cat(paste(formatC(x$x[iK[-c(1,nK)]], digits = max(1, digits - 1)),
	      collapse = ", "), ".\n",
	"RSS = ", formatC(RSS, digits=digits),"; R^2 = ",
	formatC(1 - RSS / sum((x$y - mean(x$y))^2), digits = digits),
	";\n needed (",paste(x$iter, collapse=","),") iterations\n",
	sep = "")
    ## TODO: mention 'tol'; summarize 'deriv.loc' and 'conv.loc'
    invisible(x)
}

knots.conreg <- function(Fn, ...) Fn$x[Fn$iKnots]
residuals.conreg <- function(object, ...) object$y - object$yf
fitted.conreg <- function(object, ...) object$yf

predict.conreg <- function (object, x, deriv = 0, ...)
{
    if (missing(x) || is.null(x))
	return(fitted(object))
    ## else
    x <- as.numeric(x)
    stopifnot(deriv == round(deriv), deriv >= 0, deriv <= 1)
    iK <- object$iKnots
    xK <- object$x[iK]
    if(deriv == 1) {
	slopes <- diff(object$yf[iK]) / diff(xK)
	approx(xK[-length(xK)], slopes, xout = x, rule = 2,
	       method = "constant")$y
    } else {
	approx(xK, object$yf[iK], xout = x, rule = 2)$y
    }
}

plot.conreg <-
    function(x, type = "l", col = 2, lwd = 1.5, show.knots = TRUE,
	     add.iSpline = TRUE, force.iSpl = FALSE,
	     xlab = "x", ylab = expression(s[c](x)),
	     sub = "simple concave regression", col.sub = col, ...)
{
    plot(x$x, x$yf, type=type, col=col, col.sub=col.sub, lwd=lwd,
	 xlab = xlab, ylab = ylab, sub = sub, ...)
    x.kn <- x$x [x$iKnots]
    y.kn <- x$yf[x$iKnots]
    if(show.knots)
	points(x.kn, y.kn, col = "orange", pch = 3, lwd=lwd)
    invisible(if(add.iSpline)
	      .add.ispline(x.kn, y.kn, lwd, force.iSpl))
}

##' Get  splines :: interpSpline() from a "conreg" object
interpSplineCon <- function(object, ...) {
    nk <- object$iKnots
    interpSpline(object$x [nk],
		 object$yf[nk], ...)
}

##' Check if the cubic interpolation spline fulfills the same convexity / concavity
##' property as the as "linear spline" [which is the solution of conreg()]
isIsplineCon <- function(object, isp, ...) {
    if(missing(isp))
	isp <- interpSplineCon(object, ...)
    ## quadratic coef -> 2nd derivative <= 0 (for concave;  >= 0 for convex):
    ## if(object$convex) all(coef(isp)[,3] >= 0) else all(coef(isp)[,3] <= 0):
    sgn <- if(object$convex) -1 else 1
    all(sgn*coef(isp)[,3] <= 0)
}

.add.ispline <- function(x, y, lwd, force) {
    isp <- interpSpline(x, y)
    if(force || all(coef(isp)[,3] <= 0))
	## quadratic coef -> 2nd derivative <= 0
	lines(predict(isp), col = "gray", lwd=lwd)
    else warning("the cubic spline through the knots is not convex and not drawn",
		 call.=FALSE)
    invisible(isp)
}

lines.conreg <-
    function(x, type = "l", col = 2, lwd = 1.5,
	     show.knots = TRUE, add.iSpline = TRUE,
	     force.iSpl = FALSE, ...)
{
    lines(x$x, x$yf, type=type, col=col, lwd=lwd, ...)
    x.kn <- x$x [x$iKnots]
    y.kn <- x$yf[x$iKnots]
    if(show.knots)
	points(x.kn, y.kn, col = "orange", pch = 3, lwd=lwd)
    invisible(if(add.iSpline)
	      .add.ispline(x.kn, y.kn, lwd, force.iSpl))
}
