cobsOld <- function(x, y, constraint = c("none", "increase", "decrease",
                       "convex", "concave", "periodic"),
                 z, minz = knots[1], maxz = knots[nknots], nz = 100,
                 knots, nknots, method = "quantile",
                 degree = 2, tau = 0.5, lambda = 0, ic = "aic",
                 knots.add = FALSE, alpha = 0.1, pointwise,
                 print.warn = TRUE, print.mesg = TRUE, trace = print.mesg,
                 coef = rep(0,nvar), w = rep(1,n),
                 maxiter = 20*n, lstart = log(big)^2, toler = 1e-6,
                 factor = 1)
{
    ##=########################################################################
    ##
    ## S interface for He and Ng (1997), ``COBS-Qualitatively Constrained
    ##   Smoothing via Linear Programming''.
    ##
    ##=########################################################################

mesg <- function(number,...) {
    ##=########################################################################
    ##
    ## S function to print intermediate messages
    ##
    ##=########################################################################
    switch(number,
           cat("\n"),
           cat("\n Now we are fitting cobs() ...\n"), # 2
           cat("\n Searching for optimal lambda. This may take a while.\n",
               "  While you are waiting, here is something you can consider\n",
               "  to speed up the process:\n",
               "      (a) Use a smaller number of knots;\n",
               "      (b) Increase `factor' to 2 or 3; \n",
               "      (c) Set lambda==0 to exclude the penalty term.\n") # 3
           )
}
    ## preamble
    ##
    if(!is.R() && !is.loaded(symbol.For("drqssbc")))# keep former S+ setup
        dyn.load("cobs_l.o")
    constraint <- match.arg(constraint)
    na.idx <- is.na(x) | is.na(y)
    x <- x[!na.idx]
    y <- y[!na.idx]
    big        <- if(is.R()) 3.4028234663852886e+38 else .Machine$single.xmax
    single.eps <- if(is.R()) 1.1920928955078125e-07 else .Machine$single.eps
    ## FIXME: `single.eps' is only used for cutting off log(<very small>)
    ## double.eps <- .Machine$double.eps
    minx <- min(x)
    maxx <- max(x)
    n <- nrq <- length(x)
    tmin <- 1/n
    ox <- order(x)
    xo <- x[ox]
    yo <- y[ox]
    nj0 <- 1
    lam <- 1
    Tlambda <- lambda
    if(length(unique(y)) == 2 && print.warn)
        ## warn(7)
        cat("\n It looks like you are fitting a binary choice model.",
            "We recommend pre-smoothing your data using smooth.spline(),",
            "loess() or ksmooth() before calling COBS\n", sep="\n ")

    if(missing(pointwise)) {
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
        n.equal   <- nrow(equal) # maybe 0 in R!
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
            if(lux < nknots) {
                nknots <- lux
            }
            knots <- ux[seq(1, lux, len = nknots)]
            names(knots) <- NULL
        }
        else
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
        cat("\n WARNING! It looks like you are using",nknots,
            "knots instead of the\n   default number of",Tnknots,"knots.\n")
    knots[1] <- knots[1]-toler
    knots[nknots] <- knots[nknots]+toler
    ## make sure that there is at least one observation between any pair of
    ## adjacent knots
    if(length(unique(cut(x, knots, labels=FALSE, include.lowest = TRUE))) != nknots - 1)
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
    if(lambda == 0) { ## quantile B-splines without penalty
        pw <- 0
        niqc <- n.iqc + n.greater + n.smaller
        nvar <- dim.o$nvar

        ##
        ## compute B-spline coefficients for quantile B-spline with stepwise
        ## knots selection, quantile B-spline with fixed knots
        ##
        qsbs.o <- qbsks(x,y,w,pw, knots,nknots, degree,Tlambda,constraint,
                        n.sub = n1000cut(n),
                        equal,smaller,greater,gradient, coef,maxiter,
                        trace, n.equal,n.smaller,n.greater,n.gradient,
                        nrq,nl1, neqc,nj0, tau,lam,tmin,kmax,lstart,
                        ks,mk.flag, knots.add, ic, print.mesg,## method,
                        factor = factor, print.warn = print.warn)
        knots <- qsbs.o$knots
        nknots <- qsbs.o$nknots
    }
    else { ## lambda !=0 : quantile smoothing B-Splines with penalty
        if(degree == 1) {
            nl1 <-  nknots - 2
            nvar <- dim.o$nvar
            niqc <- n.iqc + n.greater + n.smaller
            pw <- rep(1, nl1)
        }
        else {
            nl1 <- 1
            nvar <- dim.o$nvar + 1      #one more parameter for sigma
            niqc <- 2 * (nknots - 1) + n.iqc + n.greater + n.smaller
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
        if(lambda < 0 && print.mesg) mesg(3)
        qsbs.o <- drqssbc(x,y,w,pw, knots, degree,Tlambda,constraint,
                          n.sub = n1000cut(n),
                          equal,smaller,greater,gradient, coef,maxiter,
                          trace, n.equal,n.smaller,n.greater,n.gradient,
                          nrq,nl1, neqc,niqc,nvar,nj0, tau,lam,tmin,kmax,lstart,
                          factor, print.warn = print.warn)
    }
    if(qsbs.o$ifl != 1) { # had problem
        if(qsbs.o$ifl < 1)
            warning(" ifl = ",qsbs.o$ifl," < 1 -- should NOT happen !!")
        else
            switch(qsbs.o$ifl,

        1, # 1
        stop("  At least one of the additional pointwise constraints is not feasible. Please check your `pointwise' argument and rerun cobs."), # 2
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
    nvar <- qsbs.o$nvar
    ##
    ## compute fitted value at z
    ##
    if(missing(z)) {
        zo <- seq(max(minz,knots[1]     + single.eps),
                  min(maxz,knots[nknots]- single.eps), len = nz)
    }
    else {
        zo <- z[order(z)]
        zo <- zo[zo > knots[1] & zo < knots[nknots]]
        nz <- length(zo)
    }
    oz <- order(zo)

    new.knots <- c(rep(knots[1], ks-1), knots, rep(knots[nknots], ks-1))
    nk <- length(new.knots)
    derivs <- as.integer(0)
    Tcoef <- qsbs.o$coef[1:nvar]
    Tncoef <- length(Tcoef) ## ==!== nvar
    fit <- .C("spline_value",
              as.double(new.knots),
              as.double(Tcoef),
              Tncoef,
              as.integer(ks),
              as.double(zo),
              as.integer(nz),
              derivs,
              y = double(nz)) $ y
    ##
    ## compute the residual
    ##
    z2 <- .C("spline_value",
             as.double(new.knots),
             as.double(Tcoef),
             Tncoef,
             as.integer(ks),
             as.double(xo),
             as.integer(n),
             derivs,
             y = double(n))
    resid <- (y[ox] - z2$y)[order(ox)]
    ##
    ## compute the confidence band
    ##
    X <- matrix(0, nz, nvar)
    z3 <- .C("spline_basis",
             as.double(new.knots),
             ncoef = as.integer(nk - ks),
             as.integer(ks),
             as.double(zo),
             derivs = integer(nz),# 0
             as.integer(nz),
             design = array(0, c(ks, nz)),
             offsets = integer(nz))
    idx <- cbind(rep(oz, rep(ks, nz)),
                 c(outer(1:ks, z3$offsets, "+")))
    X[idx] <- z3$design
    if(any(ibig <- abs(X) > big)) { ## make sure that drqssbc won't overflow
        X[ibig] <- X[ibig] / big^0.5
        warning("re-scaling ", sum(ibig), "values of spline basis `X'")
    }
    chisq.alpha <- qchisq(1 - alpha, qsbs.o$k)
    z.alpha <- qt(1 - alpha/2, n - qsbs.o$k)
    ## MM: not clear if this is good: QR( X' X ) ; later inverse
    Tqr <- qr(t(qsbs.o$pseudo.x) %*% qsbs.o$pseudo.x)
    if(Tqr$rank != dim(qsbs.o$pseudo.x)[2])
        stop("The pseudo design matrix is singular; this can most likely be solved by using a smaller lambda when obtaining qsbs.out"
             )
    Q <- solve(Tqr)                     #note: no n here
    xQx <- diag(X %*% Q %*% t(X))

##Dbg cat("after ` Q <- solve(Tqr) ' and `xQx <- diag(X %*% Q %*% t(X))'")
##Dbg browser()

    s <- shat(resid, tau, alpha, hs = TRUE)
    cn <- sqrt(xQx * tau * (1 - tau))
    an <- sqrt(chisq.alpha) * cn * s
    bn <- z.alpha * cn * s

    return(coef = Tcoef, fit = fit, resid = resid,
           z = zo, knots = knots, ifl = qsbs.o$ifl, icyc = qsbs.o$icyc,
           k = min(qsbs.o$k, nknots-2+ks), lambda = qsbs.o$lambda,
           pp.lambda = if(lambda < 0) qsbs.o$pp.lambda,
           sic       = if(lambda < 0) log(qsbs.o$sic),
           cb.lo = fit - an, cb.up = fit + an,
           ci.lo = fit - bn, ci.up = fit + bn)
}## cobsOld()
