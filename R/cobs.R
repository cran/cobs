####-*- mode: R; kept-old-versions: 12;  kept-new-versions: 20; -*-

n1000cut <- function(n) floor(ifelse(n > 1000, 671+log(n)^3, n))
## original had  "670 +" but that was *not* monotone

if(!exists("is.R", mode="function"))
    is.R <- function()
    exists("version") && !is.null(vl <- version$language) && vl == "R"

## S+ does not allow "cut(*, labels = FALSE)" -- use cut00() for compatibility:
if(is.R()) {
    .First.lib <- function(lib, pkg) library.dynam("cobs",pkg,lib)
    cut00 <- function(x, breaks)
        cut.default(x, breaks, labels = FALSE, include.lowest = TRUE)
} else { ## S-plus  (tested only with S+ 6.0):
    cut00 <- function(x, breaks)
        as.integer(cut.default(x, breaks, include.lowest = TRUE))
}

cobs <- function(x, y, constraint = c("none", "increase", "decrease",
                       "convex", "concave", "periodic"),
                 knots, nknots, method = "quantile",
                 degree = 2, tau = 0.5, lambda = 0, ic = "aic",
                 n.sub = n1000cut(n),
                 knots.add = FALSE, pointwise = NULL,
                 print.warn = TRUE, print.mesg = TRUE, trace = print.mesg,
                 coef = rep(0,nvar), w = rep(1,n),
                 maxiter = 20*n, lstart = 7872, toler.kn = 1e-6,
                 eps = .Machine$double.eps, factor = 1)
{
    ##=########################################################################
    ##
    ## S interface for He and Ng (1997), ``COBS-Qualitatively Constrained
    ##   Smoothing via Linear Programming''.
    ##
    ##=########################################################################

    ## lstart had default = log(big)^2 , where big = .Machine$single.xmax
    ##        this is 7871.74216945922  for IEEE
    ## preamble
    ##
    cl <- match.call()
    if(!is.R() && !is.loaded(symbol.For("drqssbc")))# keep former S+ setup
        dyn.load("cobs_l.o")
    constraint <- match.arg(constraint)
    na.idx <- is.na(x) | is.na(y)
    x <- x[!na.idx]
    y <- y[!na.idx]
    minx <- min(x)
    maxx <- max(x)
    n <- nrq <- length(x)
    tmin <- 1/n
    ox <- order(x)
    xo <- x[ox]
    nj0 <- 1
    lam <- 1
    Tlambda <- lambda
    if(length(unique(y)) == 2 && print.warn)
        ## warn(7)
        cat("\n It looks like you are fitting a binary choice model.",
            "We recommend pre-smoothing your data using smooth.spline(),",
            "loess() or ksmooth() before calling COBS\n", sep="\n ")

    if(is.null(pointwise)) {
        equal <- greater <- smaller <- gradient <- NULL
        n.equal <- n.greater <- n.smaller <- n.gradient <- 0
    }
    else { ## `pointwise'
        if(!is.matrix(pointwise)|| dim(pointwise)[2] != 3)
            stop("  Argument `pointwise' has to be a three-column matrix.")
        kind <- pointwise[,1] # .Alias
        equal   <- pointwise[kind ==  0, , drop = FALSE]
        greater <- pointwise[kind ==  1, , drop = FALSE]
        smaller <- pointwise[kind == -1, , drop = FALSE]
        gradient<- pointwise[kind ==  2, , drop = FALSE]
        n.equal   <- nrow(equal) # maybe 0 in R and Sv4 (i.e. S+ (>=5)
        n.greater <- nrow(greater)
        n.smaller <- nrow(smaller)
        n.gradient<- nrow(gradient)
    }
    ##
    ## generate default knots sequence
    ##
    Tnknots <- if(lambda == 0) 6 else 20
    if(missing(knots)) {
        mk.flag <- TRUE
        if(missing(nknots))
            nknots <- Tnknots
        if(method == "quantile") {
            lux <- length(ux <- unique(xo))
            if(lux <= nknots) {
                nknots <- lux
                knots <- ux
            } else { # nknots < lux: take ``rounded'' quantiles
                knots  <- ux[seq(1, lux, len = nknots)]
            }
        }
        else ## "equidistant" :
            knots <- seq(xo[1],xo[n], len = nknots)
    }
    else {
        knots <- sort(knots)
        mk.flag <- missing(nknots) || nknots != length(knots)
        names(knots) <- NULL
        nknots <- length(knots)
        if(knots[1] > minx || knots[nknots] < maxx)
            stop("  The range of knots should cover the range of x.")
    }
    if(nknots < 2) stop("  A minimum of two knots is needed.")
    if(nknots == 2) {
        if(lambda == 0 && mk.flag)
            stop("  Can't perform automatic knot selection with nknots == 2.")
        else if(degree == 1)
            stop("  You need at least 3 knots when lambda!=0 and degree==1")
    }
    if(nknots != Tnknots && print.warn)
        cat("\n You are using",nknots,
            "knots instead of the default number of",Tnknots,"knots.\n")
    knots[1] <- knots[1]-toler.kn
    knots[nknots] <- knots[nknots]+toler.kn
    ## make sure that there is at least one observation between any pair of
    ## adjacent knots
    if(length(unique(cut00(x, knots))) != nknots - 1)
        stop(" There is at least one pair of adjacent knots that contains no observation.")
    kmax <- nknots
    ##
    ## set up proper dimension for the pseudo design matrix
    ##
    dim.o <- getdim(degree,nknots,constraint)
    neqc <- dim.o$n.eqc + n.equal + n.gradient
    nl1 <- 0
    ks <- dim.o$ks
    n.iqc <- dim.o$n.iqc
    niqc  <- n.iqc + n.greater + n.smaller
    if(lambda == 0) { ## quantile B-splines without penalty
        pw <- 0
        nvar <- dim.o$nvar
        ##
        ## compute B-spline coefficients for quantile B-spline with stepwise
        ## knots selection, quantile B-spline with fixed knots
        ##
        rr <- qbsks(x,y,w,pw, knots,nknots, degree,Tlambda,constraint, n.sub,
                    equal,smaller,greater,gradient, coef,maxiter,
                    trace, n.equal,n.smaller,n.greater,n.gradient,
                    nrq,nl1, neqc,nj0, tau,lam,tmin,kmax,lstart,
                    ks,mk.flag, knots.add, ic,
                    print.mesg = print.mesg, factor=factor,
                    tol.kn = toler.kn, eps = eps, print.warn = print.warn)
        knots <- rr$knots
        nknots <- rr$nknots
    }
    else { ## lambda !=0 : quantile smoothing B-Splines with penalty
        if(degree == 1) {
            nl1 <-  nknots - 2
            nvar <- dim.o$nvar
            pw <- rep(1, nl1)
        }
        else {
            nl1 <- 1
            nvar <- dim.o$nvar + 1      #one more parameter for sigma
            niqc <- niqc + 2 * (nknots - 1)
            pw <- rep(1, nknots-1)
        }
        if(lambda < 0) {                # lambda is chosen by sic
            lam <- -1
            Tlambda <- 1
            nj0 <- 50*n ## = nsol <<-- determines size of sol[,] !
        }
        ##
        ## compute B-spline coefficients for quantile smoothing B-spline
        ##
        if(lambda < 0 && print.mesg)
            cat("\n Searching for optimal lambda. This may take a while.\n",
                "  While you are waiting, here is something you can consider\n",
                "  to speed up the process:\n",
                "      (a) Use a smaller number of knots;\n",
                "      (b) Increase `factor' to 2 or 3; \n",
                "      (c) Set lambda==0 to exclude the penalty term.\n")# 3

        rr <- drqssbc(x,y,w,pw, knots, degree,Tlambda,constraint, n.sub,
                      equal,smaller,greater,gradient, coef,maxiter,
                      trace, n.equal,n.smaller,n.greater,n.gradient,
                      nrq,nl1, neqc,niqc,nvar,nj0, tau,lam,tmin,kmax,lstart,
                      factor, eps, print.warn)
    }
    if(rr$ifl != 1) { # had problem
        if(rr$ifl < 1)
            warning(" ifl = ",rr$ifl," < 1 -- should NOT happen !!")
        else
            switch(rr$ifl,

        1, # 1
        stop("At least one of the additional pointwise constraints is not feasible.\n Please check your `pointwise' argument and rerun cobs."), # 2
        ## warn(6,maxiter)
        cat("\n WARNING! The algorithm has not converged after",maxiter,
            "iterations.\n",
            "Increase the `maxiter' counter and restart cobs with both\n",
            "`coef' and `knots' set to the values at the last iteration.\n"), # 3
        stop("  The problem is ill-conditioned."), # ifl = 4
        stop(" ifl = 5 -- shold not have happened"), # 5
        warning(" ifl = 6 --- nj0 = nsol = 50*n  is too small!"), # 6
        stop(" ifl = 7 : lambda < 0  *and*  tau outside [0,1]")  # 7
                   )
    }
    nvar <- rr$nvar
    Tcoef <- rr$coef[1:nvar]
    ##
    ## compute the residual
    ##
    y.hat <- .splValue(degree, knots, Tcoef, xo)
    y.hat <- y.hat[order(ox)]# original (unsorted) ordering

    r <- list(call = cl,
              tau = tau, degree = degree, constraint = constraint,
              pointwise = pointwise,
              coef = Tcoef, knots = knots, ifl = rr$ifl, icyc = rr$icyc,
              k = min(rr$k, nknots-2+ks), k0 = rr$k,
              x.ps = rr$pseudo.x,
              resid = y - y.hat, fitted = y.hat,
              SSy = sum((y - mean(y))^2),
              lambda = rr$lambda,
              pp.lambda = if(lambda < 0) rr$pp.lambda,
              sic       = if(lambda < 0) log(rr$sic))
    class(r) <- "cobs"
    r
}## cobs()

print.cobs <- function(x, digits = getOption("digits"), ...) {
    if(!is.numeric(lam <- x$lambda))
        stop("`x' is not a valid \"cobs\" object")
    cat("COBS ", if(lam == 0) "regression" else "smoothing",
        " spline (degree = ", x$degree, ") from call:\n  ", deparse(x$call),
        "\n{tau=",format(x$tau,digits),"}-quantile",
        ";  dimensionality of fit: ",x$k," (",x$k0,")\n", sep="")
    nkn <- length(x$knots)
    cat("knots[1 .. ", nkn,"]: ", sep = "")
    chk <- format(x$knots[if(nkn <= 5) 1:nkn else c(1:4, nkn)], digits=digits)
    if(nkn > 5) chk[4] <- "... "
    cat(chk, sep = ", "); cat("\n")
    if(lam != 0) {
        cat("lambda =", format(lam, digits = digits))
        if((nlam <- length(x$pp.lambda))) {
            cat(", selected via SIC, out of", nlam, "ones.")
        }
        cat("\n")
    }
    invisible(x)
} # print

summary.cobs <- function(object, digits = getOption("digits"), ...) {
    if(!is.numeric(lam <- object$lambda))
        stop("`object' is not a valid \"cobs\" object")
    print(object, digits = digits, ...)# includes knots
    if(!is.null(pw <- object$pointwise)) {
        cat("with",nrow(pw),"pointwise constraints\n")
    }
    cat("coef  :\n"); print(object$coef, digits = digits, ...)
    tau <- object$tau
    if(abs(tau - 0.50) < 1e-6)
        cat("R^2 = ", round(100 * (1 - sum(object$resid^2) / object$SSy), 2),
            "% ;  ", sep="")
    k <- sum((r <- resid(object)) <= 0)
    n <- length(r)
    cat("empirical tau (over all): ",k,"/",n,"  = ", format(k/n,digits),
        " (target tau : ",tau,")\n", sep="")

    ## add more -- maybe finally *return* an object and define
    ## print.summary.cobs <- function(x, ...)

} # summary

residuals.cobs <- function (object, ...) object$resid
fitted.cobs <- function (object, ...) object$fitted

predict.cobs <-
    function(object, z, minz = knots[1], maxz = knots[nknots], nz = 100,
             interval = c("none", "confidence", "simultaneous", "both"),
             level = 0.95, ...)
{
    if(is.null(knots <- object$knots) ||
       is.null(coef  <- object$coef)  ||
       is.null(degree<- object$degree)  ||
       is.null(tau   <- object$tau)) stop("not a valid `cobs' object")

    interval <- match.arg(interval)

    big        <- if(is.R()) 3.4028234663852886e+38 else .Machine$single.xmax
    ##IN
    single.eps <- if(is.R())1.1920928955078125e-07 else .Machine$single.eps

    nknots <- length(knots)
    ord <- as.integer(degree + 1)
    nvar   <- length(coef)

    ##DBG cat("pr..cobs(): (ord, nknots, nvar) =",ord, nknots, nvar, "\n")

    ##
    ## compute fitted value at z
    ##
    ## MM: why should z be *inside* (even strictly) the knots interval? _FIXME_
    if(missing(z)) {
        if(minz >= maxz) stop("minz >= maxz")
        ##NOT YET (for "R CMD check" compatibility):
        ## zo <- seq(minz, maxz, len = nz)
        zo <- seq(max(minz,knots[1]     + single.eps),
                  min(maxz,knots[nknots]- single.eps), len = nz)
    }
    else {
        zo <- sort(z)
        ##IN zo <- zo[zo > knots[1] & zo < knots[nknots]]
        nz <- length(zo)
    }

    fit <- .splValue(degree, knots, coef, zo)

    if(interval != "none") {
        ##
        ## compute confidence bands
        ## both (pointwise and simultaneous : as cheap as only one !
        ##
        z3 <- .splBasis(ord = ord, knots, ncoef = nknots + degree - 1, xo = zo)
        idx <- cbind(rep(1:nz, rep(ord, nz)),
                     c(outer(1:ord, z3$offsets, "+")))
        X <- matrix(0, nz, nvar)
        X[idx] <- z3$design
        if(any(ibig <- abs(X) > big)) { ## MM: no sense here!
            X[ibig] <- X[ibig] / big^0.5
            warning("re-scaling ", sum(ibig), "values of spline basis `X'")
        }

        Tqr <- qr(crossprod(object$x.ps))
        if(Tqr$rank != dim(object$x.ps)[2])
            stop("The pseudo design matrix is singular; this can most likely be solved by using a smaller lambda")# when obtaining qsbs.out
        ## Improved way of computing
        ## xQx = diag(X %*%  solve(Tqr)  %*% t(X)) :
        tX <- t(X)
        xQx <- colSums(qr.coef(Tqr, tX) * tX)

        res <- object$resid
        n <- length(res)

        s <- shat(res, tau, 1 - level, hs = TRUE)
        cn <- sqrt(xQx * tau * (1 - tau))
        sde <- cn * s
        an <- sqrt(qchisq(level, object$k0))   * sde
        bn <- qt((1 + level)/2, n - object$k0) * sde

        cbind(z = zo, fit = fit,
              cb.lo = fit - an, cb.up = fit + an,
              ci.lo = fit - bn, ci.up = fit + bn)
    }# interval
    else
        cbind(z = zo, fit = fit)
} # predict



getdim <- function(degree, nknots,
                   constraint = c("none", "increase", "decrease",
                   "convex", "concave", "periodic"))
{
    ##=########################################################################
    ##
    ## Compute the appropriate dimension for the pseudo design matrix
    ##
    ##=########################################################################
    if(degree == 1)	ks <- 2
    else if(degree == 2)ks <- 3
    else stop("degree has to be either 1 or 2")
    constraint <- match.arg(constraint)
    nvar <- nknots - 2 + ks
    n.eqc <- # the number of EQuality Constraints
        if(constraint == "periodic") 2 else 0
    n.iqc <- # the number of IneQuality Constraints
    if(constraint == "increase" || constraint == "decrease")
        nknots - as.integer(degree == 1)
    else if(constraint == "concave" || constraint == "convex")
        nknots - 1 - as.integer(degree == 1)
    else if(constraint == "periodic" || constraint == "none")
         0
    return(n.iqc = n.iqc, n.eqc = n.eqc, ks = ks, nvar = nvar)
} ## getdim()

### These are (only) used for confidence intervals :

shat <- function(residual, tau, alpha, hs)
{
    ##=########################################################################
    ##
    ## sparsity estimate from empirical quantile function using residuals
    ##
    ##=########################################################################
    residual <- sort(residual)
    n <- length(residual)
    residual <- c(residual[1], residual, residual[n])
    grid <- c(0, seq(0.5/n, 1 - 0.5/n, 1/n), 1)
    hn <- dn(tau, n, hs = hs, alpha)
    ## for small n,  tau +/- hn might be outside [0,1]
    bound <- pmax(0, pmin(1, c(tau - hn, tau + hn)))
    idx <- cut00(bound, grid)
    lambda <- bound * n - (idx - 1) + 0.5
    return(diff(lambda * residual[idx + 1] +
                (1 - lambda) * residual[idx])/(2 * hn))
}

dn <- function(p, n, hs = FALSE, alpha)
{
    ##=########################################################################
    ##
    ## compute window width for sparsity estimator
    ## at quantile p, level = 1-alpha,  n observations
    ## according to
    ##   Hall and Sheather (1988),  <<-  hs=TRUE,   or
    ##   Bofinger (1975),           <<-  hs=FALSE
    ##=########################################################################
    x0 <- qnorm(p)
    f0 <- dnorm(x0)
    if(as.logical(hs))
        n^(-1/3) * qnorm(1 - alpha/2)^(2/3) *
            ((1.5 * f0^2)/(2 * x0^2 + 1))^(1/3)
    else n^-0.2 * ((4.5 * f0^4)/(2 * x0^2 + 1)^2)^ 0.2
}
