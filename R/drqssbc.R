### used to be part of ./cobs.R

drqssbc <- function(x,y, w = rep(1,n), pw, knots, degree,Tlambda, constraint,
                    n.sub = n1000cut(nrq),
		    equal,smaller, greater,gradient, coef, maxiter = 20*n,
		    trace = 1,
                    n.equal = nrow(equal), n.smaller = nrow(smaller),
                    n.greater = nrow(greater), n.gradient = nrow(gradient),
		    nrq = length(x), nl1, neqc,niqc, nvar,nj0,
		    tau = 0.50, lam, tmin, kmax, lstart, factor = 1,
                    eps = .Machine$double.eps, print.warn = TRUE)
{
    ##=########################################################################
    ##
    ## Estimate the B-spline coefficients for quantile *smoothing* spline, using
    ##		Ng (1996)  `An Algorithm for Quantile Smoothing Splines',
    ##		Computational Statistics & Data Analysis, 22, 99-118.
    ##
    ##=########################################################################
    big	       <- if(is.R()) 3.4028234663852886e+38 else .Machine$single.xmax
    single.eps <- if(is.R()) 1.1920928955078125e-07 else .Machine$single.eps
    toler.kn <- 1e-6
    ## Note   nrq != length(x) , e.g., in case of sub.sampling+ fit full
    if(lam >= 0)
        n <- nrq
    else if((n.old <- nrq) != (n <- as.integer(n.sub))) {
	##
	## sub-sampling for smoothing B-splines parametric programming
        ## select a sub-sample of size n.sub
	##
	sub.idx <- seq(1,n.old, length = n)
	x.old <- x; x <- x[sub.idx]
	y.old <- y; y <- y[sub.idx]
	w.old <- w; w <- w[sub.idx]
    }
    if(degree == 1) {
	X <- l1.design(x,w,constraint,equal,smaller,greater,gradient,knots,
		       pw,n.equal,n.smaller,n.greater,n.gradient,
		       nrq=n, nl1,neqc,niqc,nvar,Tlambda)
	niqc1 <- 0
    }
    else {
	X <- loo.design(x,w,constraint,equal,smaller,greater,gradient,knots,
			pw,n.equal,n.smaller,n.greater,n.gradient,
			nrq=n, nl1,neqc,niqc,nvar,Tlambda)
	niqc1 <- if(Tlambda == 0) 0 else 2*(length(knots) - 1)
    }
    if(any(ibig <- abs(X) > big)) { ## make sure that drqssbc won't overflow
	X[ibig] <- X[ibig] / big^0.5
	warning("re-scaling ", sum(ibig), "values of Lp-design `X'")
    }
    Tnobs <- nrow(X)
    Tequal   <- if(n.equal    > 0)    equal[,3] # else NULL
    Tsmaller <- if(n.smaller  > 0) -smaller[,3]
    Tgreater <- if(n.greater  > 0)  greater[,3]
    Tgradient<- if(n.gradient > 0) gradient[,3]

    Y <- c(y*w, rep(0,nl1), Tequal, Tgradient,
	   rep(0,Tnobs-n-nl1-n.equal-n.gradient-n.smaller-n.greater),
	   Tsmaller,Tgreater)
    ##storage.mode(X) <- "single" # would round to ~ 7 digits in S+, not in R
    d <- matrix(0.0, Tnobs + 5, nvar + 2) # double
    sol <- matrix(0.0, nvar + 6, nj0) # double -- to contain "sol"ution
    z0 <- .Fortran("drqssbc",
		   as.integer(n),
		   as.integer(nl1),
		   as.integer(neqc),
		   as.integer(niqc),
		   as.integer(niqc1),
		   as.integer(nvar),
		   integer(1),		# nact
		   ifl = integer(1),
		   as.integer(maxiter), # mxs
		   as.integer(trace),
		   X = as.double(t(X)), # e
		   as.integer(nvar),	# ner
		   coef = as.double(coef),# x
		   as.double(Y),	# f
		   obj = double(1),	# erql1n
		   resid = double(Tnobs),
		   integer(Tnobs),	# indx
		   double(((3 * nvar + 13) * nvar + 2)/2 + 2 * Tnobs),# w
		   nt = integer(1),	# nt
		   as.integer(nj0),	# nsol
		   sol = sol,
		   as.double(c(tau, lam)),# == tl[1:2] == (t, lam)
		   as.double(toler.kn),
		   as.double(big),
		   as.double(eps),
		   icyc = integer(2),
		   as.double(tmin),
		   k = integer(1),
		   as.integer(kmax),	# k0
		   as.double(lstart),
		   as.double(factor),
                   PACKAGE = "cobs")
    sol <- z0$sol[,1:z0$nt]
    names(z0$icyc) <- c("icyc", "tot.cyc")
    if(lam < 0) {
        ##
        ## search for optimal lambda
        ##
        ifl.idx <-
            if(maxiter > 20*n) sol[3,] != 3 # trap ifl=3
            else rep(TRUE, ncol(sol))
        fidel <- sol[4,ifl.idx]
        eff.k <- sol[6,ifl.idx]
        sic <- fidel/n * n^(eff.k/(2*n))
        mlam.idx <- min((1:ncol(sol[,ifl.idx]))[sic == min(sic)])
        ##
        ## trap infeasible solution when lam<0
        ##
        if(sol[,ifl.idx][3,mlam.idx] == 0)
            return(list(ifl = 2))

        Tcoef <- sol[,ifl.idx][7:nrow(sol),mlam.idx]
        if(print.warn && sol[3,mlam.idx] != 3) {
            if(sol[,ifl.idx][6,mlam.idx] >= length(knots))
                ## warn(2)
                cat("\n WARNING! Since the optimal lambda chosen by SIC corresponds to the",
                    "   roughest possible fit, you might want to consider doing one of",
                    "   the following:",
                    "   (1) plot the components $sic against $pp.lambda of cobs to see",
                    "   if a bigger lambda value at another local minimum of $sic will",
                    "   yield a more reasonable fit;",
                    "   (2) increase the number of knots.\n", sep="\n ")
            else if(abs(sic[mlam.idx]-sic[length(sic)]) < single.eps * sic[length(sic)])
                ## warn(3)
                cat("\n WARNING! Since the optimal lambda chosen by SIC corresponds to the",
                    "   roughest possible fit, you might want to plot the components",
                    "   $sic against $pp.lambda of cobs to see if a bigger lambda value",
                    "   at another local minimum of $sic will yield a more reasonable fit.\n",
                    sep="\n ")

            if(sol[,ifl.idx][2,mlam.idx] == lstart)
                ## warn(4)
                cat("\n WARNING!  Since the optimal lambda chosen by SIC reached the smoothest",
                    "   possible fit allowed by `lstart', you might want to rerun cobs with",
                    "   a larger `lstart' value to see if it makes a difference if you haven't",
                    "   done so.\n", sep="\n ")
        }
        if(n.old != n) { ## was sub sampling; refit full sample for one Tlambda
            Tlambda <- sol[,ifl.idx][2,mlam.idx]
	    rqss <- drqssbc(x.old,y.old,w.old, pw,knots,degree, Tlambda,
			    constraint, n.sub,
			    equal,smaller,greater,gradient,
			    Tcoef,maxiter,trace,
			    n.equal,n.smaller,n.greater,n.gradient,
			    nrq = n.old, nl1,neqc,niqc,nvar, nj0 = 1,
			    tau,lam = 1, tmin,kmax,lstart,
			    factor, eps, print.warn)
	    list(coef = rqss$coef, ifl = sol[,ifl.idx][3,mlam.idx],
		 icyc = rqss$icyc, nvar = nvar, lambda = Tlambda,
		 pp.lambda = sol[,ifl.idx][2,], sic = sic,
		 k = min(rqss$k, length(knots)-2+degree+1), pseudo.x = X)
	}
	else
	    list(coef = Tcoef, ifl = sol[,ifl.idx][3,mlam.idx],
		 icyc = z0$icyc, nvar = nvar,
		 lambda	   = sol[,ifl.idx][2,mlam.idx],
		 pp.lambda = sol[,ifl.idx][2,], sic = sic,
		 k = min(sol[,ifl.idx][6,mlam.idx], length(knots)-2+degree+1),
		 pseudo.x = X)
    }
    else # `lam >= 0'
	list(coef = z0$coef,
	     fidel = sum((tau-(-z0$resid[1:n] < 0))*(-z0$resid[1:n]))*2,
	     k = min(z0$k, length(knots)-2+degree+1), ifl = z0$ifl,
	     icyc = z0$icyc, nvar = nvar, lambda = Tlambda, pseudo.x = X)
}## drqssbc
