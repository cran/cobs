### used to be part of ./cobs.R

qbsks <- function(x,y,w,pw, knots,nknots, degree,Tlambda, constraint,
                  n.sub = n1000cut(n),
                  equal,smaller, greater,gradient, coef,maxiter,
                  trace, n.equal,n.smaller,n.greater,n.gradient,
                  nrq,nl1, neqc, nj0, tau,lam,tmin,kmax,lstart,
                  ks,mk.flag, knots.add, ic, print.mesg,
                  factor, tol.kn = 1e-6, eps = .Machine$double.eps, print.warn)
{
    ##=########################################################################
    ##
    ## Compute B-spline coefficients for quantile B-spline with stepwise knots
    ## selection, quantile B-spline with fixed knots (REGRESSION SPLINE), using
    ##		Ng (1996)  `An Algorithm for Quantile Smoothing Splines',
    ## 		Computational Statistics & Data Analysis, 22, 99-118.
    ##
    ##=########################################################################

    ## single.eps <- if(is.R()) 1.1920928955078125e-07 else .Machine$single.eps
    ## double.eps <- .Machine$double.eps
    smll.log <- 50*floor(log(.Machine$double.xmin)/50) # heavily rounding down
    finite.log <- function(x) {
        r <- log(x)
        if(is.na(r) || r > -Inf) r else smll.log # = -750 for IEEE arithmetic
    }
    n <- n.old <- nrq # "n <-" : for default of n.sub
    n <- nrq <- n.sub <- as.integer(n.sub)
    if(n != n.old) {
        ##
        ## select a sub-sample of size n.sub
        ##
        sub.idx <- seq(1,n.old, length = n)
	x.old <- x; x <- x[sub.idx]
	y.old <- y; y <- y[sub.idx]
	w.old <- w; w <- w[sub.idx]
    }
    xo <- x[order(x)]
    logn <- log(n)
    const <- switch(ic, aic = 2/n, sic = logn/n)
    constraint.old <- constraint
    if(mk.flag) { ##  perform first step knots selection
        if(print.mesg) cat("\n Performing general knot selection ...\n")# 4
        Tic <- Tifl <- double(nknots-1)
        for(i in 1:(nknots-1)) {
            Tknots <- knots[seq(1,nknots, len = i+1)]
            Tnknots <- length(Tknots)
            if(Tnknots == 2 && degree == 1 &&
               (constraint == "convex" || constraint == "concave"))
                ## guard against convex fit when only 2 knots are used
                constraint <- "none"
            dim.o <- getdim(degree, Tnknots, constraint)
            ks <- dim.o$ks
            n.iqc <- dim.o$n.iqc
            Tnvar <- dim.o$nvar
            niqc <- n.iqc + n.greater + n.smaller
            rqss <- drqssbc(x,y,w,pw,Tknots,degree,Tlambda,constraint, n.sub,
                            equal,smaller,greater,gradient,coef,maxiter,
                            trace,n.equal,n.smaller,n.greater,n.gradient,
                            nrq,nl1,neqc,niqc,Tnvar,nj0,
                            tau,lam,tmin,kmax,lstart, factor,eps, print.warn)
            constraint <- constraint.old
            Tic[i] <- finite.log(rqss$fidel) -logn + (i-1+ks)*const
            Tifl[i] <- rqss$ifl
        }
        Tic.min <- min(Tic)
        Tifl.final <- Tifl[Tic == Tic.min]
        nknots.min <- min((1:(nknots-1))[Tic == Tic.min])
        if(nknots.min == 1 && Tifl.final != 1) {
            ##
            ## when the chosen nknots=2, guard against anomaly of ifl=5 when
            ## constraint=='periodic', or ifl=2 when the chosen model is infeasible.
            ##
            Tic.min <- min(Tic[2:length(Tic)])
            Tifl.final <- Tifl[Tic == Tic.min]
            nknots.min <- min((1:(nknots-1))[Tic == Tic.min])
        }
        if(Tifl.final == 2 || Tifl.final == 3 || Tifl.final == 4)
            return(list(ifl = Tifl.final))
        if((nknots.min + 1) == nknots && print.warn)
            ## warn(5,nknots,ic)
            cat("\n WARNING! Since the number of ",nknots," knots selected by ",
                ic," reached the\n",
                "  upper bound during general knot selection, you might want to rerun\n",
                "  cobs with a larger number of knots. \n")

        knots <- knots[seq(1,nknots, len = nknots.min+1)]
        names(knots) <- NULL
        ##
        ## perform knots deletion
        ##
        delete <- TRUE
        if(print.mesg) cat("\n Deleting unnecessary knots ...\n")# 5
        while(delete && nknots.min > 1) {
            Tnknots <- length(knots)
            Tic1 <- rep(0,(Tnknots-2))
            Tnknots.1 <- Tnknots - 1
            if(Tnknots.1 == 2 && degree == 1  &&
               (constraint == "convex" || constraint == "concave"))
                ## guard against convex fit when only 2 knots are used
                constraint <- "none"
            dim.o <- getdim(degree,Tnknots.1,constraint)
            ks <- dim.o$ks
            n.iqc <- dim.o$n.iqc
            Tnvar <- dim.o$nvar
            niqc <- n.iqc + n.greater + n.smaller
            Tcoef <- rep(0,Tnvar)
            for(i in 2:(Tnknots-1)) {
                Tknots <- knots[-i]
                rqss <- drqssbc(x,y,w,pw,Tknots, degree,Tlambda,constraint, n.sub,
                                equal,smaller,greater,gradient,
                                Tcoef, maxiter, trace,
                                n.equal,n.smaller,n.greater,n.gradient,
                                nrq,nl1,neqc,niqc,Tnvar,nj0,
                                tau,lam,tmin,kmax,lstart,factor, eps,print.warn)
                constraint <- constraint.old
                Tic1[i-1] <- finite.log(rqss$fidel)-logn+(Tnknots.1-2+ks)*const
                Tcoef <- rqss$coef
            }
            Tic1.min <- min(Tic1)
            idx.del <- min((2:(Tnknots-1))[Tic1 == Tic1.min])
            if((delete <- Tic1.min <= Tic.min)) {
                Tic.min <- Tic1.min
                if(print.mesg >= 3)
                    cat("\n A knot at ",signif(knots[idx.del]),
                        " is deleted.\n")# 6
                knots <- knots[-idx.del]
                nknots.min <- length(knots)-1
            }
        }
        if(print.mesg >= 2) cat("\n No more knot to be deleted.\n") # 7
        ##
        ## perform knots addition
        ##
        if(knots.add) {
            add <- TRUE
            Tnknots <- length(knots)
            if(print.mesg) cat("\n Searching for missing knots ...\n")# 8
            while(add && Tnknots < nknots) {
                Tic2 <- double(Tnknots-1)
                knots.add <- (knots[1:(Tnknots-1)]+knots[2:Tnknots])/2
                Tnknots.1 <- Tnknots + 1
                dim.o <- getdim(degree,Tnknots.1,constraint)
                ks <- dim.o$ks
                n.iqc <- dim.o$n.iqc
                Tnvar <- dim.o$nvar
                niqc <- n.iqc + n.greater + n.smaller
                Tcoef <- double(Tnvar)
                for(i in 1:(Tnknots-1)) {
                    Tknots <- sort(c(knots,knots.add[i]))
                    if(length(unique(cut00(x, Tknots))) != Tnknots)
                        Tic2[i] <- Tic.min+1
                    else {
                        rqss <-
                            drqssbc(x,y,w,pw,Tknots,degree,Tlambda, constraint, n.sub,
                                    equal,smaller,greater,gradient,
                                    Tcoef,maxiter, trace,
                                    n.equal,n.smaller,n.greater,n.gradient,
                                    nrq,nl1,neqc,niqc,Tnvar,nj0,
                                    tau,lam,tmin,kmax,lstart, factor, eps,print.warn)
                        Tic2[i] <- finite.log(rqss$fidel) -logn +
                            (Tnknots.1-2+ks)*const
                        Tcoef <- rqss$coef
                    }
                }
                Tic2.min <- min(Tic2)
                idx.add <- min((1:(Tnknots-1))[Tic2 == Tic2.min])
                if((add <- Tic2.min <= Tic.min)) {
                    Tic.min <- Tic2.min
                    knots <- sort(c(knots,knots.add[idx.add]))
                    if(print.mesg >= 2)
                        cat("\n A knot at ",signif(knots.add[idx.add]),
                            " is added.\n")# 9
                }
                Tnknots <- length(knots)
            }# end while(add ..)
            if(print.mesg >= 2) cat("\n No more knot to be added.\n")# 10
        } # (knots.add)
        if(print.mesg) cat("\n Computing the final fit ...\n")# 11
    } # end if(mk.flag)

    ##
    ## compute the B-spline coefficients for the full sample
    ##
    nknots <- length(knots)
    rk <- diff(range(knots))
    knots[1] <- knots[1] - tol.kn*rk
    knots[nknots] <- knots[nknots] + tol.kn*rk
    if(n != n.old) {
        x <- x.old
        y <- y.old
        w <- w.old
        nrq <- n.old
    }
    if(nknots == 2 && (constraint == "convex" || constraint == "concave") &&
       degree == 1) # guard against convex fit when only 2 knots are used
        constraint <- "none"
    dim.o <- getdim(degree,nknots,constraint)
    ks <- dim.o$ks
    n.iqc <- dim.o$n.iqc
    Tnvar <- dim.o$nvar
    niqc <- n.iqc + n.greater + n.smaller
    rqss <- drqssbc(x,y,w,pw, knots, degree,Tlambda,constraint, n.sub,
                    equal,smaller,greater,gradient, coef,maxiter,
                    trace, n.equal,n.smaller,n.greater,n.gradient,
                    nrq,nl1, neqc,niqc, Tnvar,nj0,
                    tau,lam,tmin,kmax,lstart,
                    factor, eps,print.warn)
    constraint <- constraint.old
    return(coef = rqss$coef, fidel = rqss$fidel,
           k = nknots-2+ks, ifl = rqss$ifl, icyc = rqss$icyc,
           knots = knots, nknots = nknots, nvar = Tnvar, lambda = Tlambda,
           pseudo.x = rqss$pseudo.x)
}# end qbsks()
