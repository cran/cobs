####  Create B-Spline Design matrices for COBS :
####
####  l1.design()  -- L_1         penalty
####  loo.design() -- L_Infinity  penalty

### used to be part of ./cobs.R

l1.design <-
function(x, w, constraint, equal, smaller, greater, gradient,
         knots, pw, n.equal, n.smaller, n.greater, n.gradient,
         nrq, nl1, neqc, niqc, nvar, lambda)
{
    ##=########################################################################
    ##
    ## Generate the pseudo design matrix for L1 penalty
    ##
    ##=########################################################################

    ## create the pseudo design
    ##
    ##x <- as.matrix(x, ncol = 1)

    ks <- 2 ## <==> degree == 1
    nk <- length(knots) + 2 # *(ks - 1)
    nkm3 <- nk - 3
    nkm4 <- nk - 4
    ncoef <- nk - ks
    ox <- order(x)
    sortx <- x[ox]
    nobs <- nrq + nl1 + neqc + niqc
    X <- matrix(0, nobs, nvar)

    z1 <- .splBasis(ord = ks, knots, ncoef, xo = sortx)
    idx1 <- cbind(rep(ox, rep(ks, nrq)),
                  c(outer(1:ks, z1$offsets, "+")))
    X[idx1] <- t(t(z1$design) * rep(w[ox], ks))
    ##
    ## formulate inequality constraints for the pseudo X
    ##
    if(lambda != 0 || constraint != "none") {
        z2 <- .splBasis(ord = ks, knots, ncoef, xo = knots[1:nkm3],
                        derivs = rep(1, nkm3))
        if(lambda != 0) {
            idx2 <- array(rep(1:nkm4, 2), c(nkm4, 2))
            idx2[, 1] <- idx2[, 1] + nrq
            X[idx2]   <-  - z2$design[1, 1:nkm4] * lambda
            idx2[, 2] <- idx2[, 2] + 2
            X[idx2]   <- z2$design[2, 2:nkm3] * lambda
            idx2[, 2] <- idx2[, 2] - 1
            X[idx2] <- (z2$design[1, 2:nkm3] -
                        z2$design[2, 1:nkm4]) * lambda
            ##
            ## assign different weight to roughness penalty
            ##
            X[nrq + 1:nkm4, ] <- X[nrq + 1:nkm4, ] * pw
        }
        if(constraint == "increase" || constraint == "decrease") {
            niqc1 <- nkm3
            idx3 <- cbind(rep(nrq+nl1+neqc + 1:niqc1, rep(ks, niqc1)),
                          c(outer(1:ks, z2$offsets, "+")))
            X[idx3] <- if(constraint == "increase") z2$design else -z2$design
        }
        else if (constraint == "periodic") {
            ##
            ## this portion corresponds to equality constraints
            ##
            niqc1 <- 0
            neqc3 <- 2
            z1.3 <- .splBasis(ord = ks, knots, ncoef,
                              xo = sortx[c(1,nrq)], derivs = c(1,1))
            idx1.3 <- cbind(rep(1:neqc3, rep(ks, neqc3)),
                            c(outer(1:ks, z1.3$offsets,"+")))
            X.temp <- matrix(0,neqc3,nvar)
            X.temp[idx1.3] <- z1.3$design
            X[nrq+nl1+neqc,  ] <- X.temp[2,] - X.temp[1,]
            X[nrq+nl1+neqc-1,] <- X[ox[nrq],]- X[ox[1],]
        }
        else if(constraint == "convex" || constraint == "concave") {
            niqc1 <- nkm4
            idx3 <- array(rep(1:nkm4, 2), c(nkm4, 2))
            idx3[, 1] <- idx3[, 1] + nrq + nl1 + neqc
            sgn <- if(constraint == "convex") +1 else -1
            X[idx3] <- -sgn* z2$design[1, 1:nkm4]
            idx3[,2] <- idx3[,2] + 2
            X[idx3] <- sgn * z2$design[2, 2:nkm3]
            idx3[,2] <- idx3[,2] - 1
            X[idx3] <- sgn * (z2$design[1, 2:nkm3] - z2$design[2, 1:nkm4])
        }
    }
    niqc2 <- n.smaller
    niqc3 <- n.greater
    if(n.smaller > 0) {
        o.smaller <- order(smaller[,2])
        smaller.o <- smaller[,2][o.smaller]
        z3.1 <- .splBasis(ord = ks, knots, ncoef, xo = smaller.o)
        idx3.1 <- cbind(rep(nrq+nl1+neqc+niqc1 + o.smaller,rep(ks,niqc2)),
                        c(outer(1:ks,z3.1$offsets,"+")))
        X[idx3.1] <- -z3.1$design
    }
    if(n.greater > 0) {
        o.greater <- order(greater[,2])
        greater.o <- greater[,2][o.greater]
        z3.2 <- .splBasis(ord = ks, knots, ncoef, xo = greater.o)
        idx3.2 <- cbind(rep(nrq+nl1+neqc+niqc1+niqc2 + o.greater,rep(ks,niqc3)),
                        c(outer(1:ks,z3.2$offsets,"+")))
        X[idx3.2] <- z3.2$design
    }
    neqc1 <- n.equal
    neqc2 <- n.gradient

    if(n.equal > 0) { ## formulate equality constraints for the pseudo X

        o.equal <- order(equal[,2])
        equal.o <- equal[,2][o.equal]
        z1.1 <- .splBasis(ord = ks, knots, ncoef, xo = equal.o)
        idx1.1 <- cbind(rep(nrq+nl1+ o.equal, rep(ks, neqc1)),
                        c(outer(1:ks, z1.1$offsets,"+")))
        X[idx1.1] <- z1.1$design
    }

    if(n.gradient > 0) { ## gradient constraints for the pseudo X

        o.gradient <- order(gradient[,2])
        gradient.o <- gradient[,2][o.gradient]
        z1.2 <- .splBasis(ord = ks, knots, ncoef, xo = gradient.o,
                          derivs = rep(1, neqc2))
        idx1.2 <- cbind(rep(nrq+nl1+neqc1+ o.gradient, rep(ks, neqc2)),
                        c(outer(1:ks, z1.2$offsets,"+")))
        X[idx1.2] <- z1.2$design
    }
    return(X)
} ## l1.design

loo.design <- function(x, w, constraint, equal, smaller, greater, gradient,
                       knots, pw, n.equal, n.smaller, n.greater, n.gradient,
                       nrq, nl1, neqc, niqc, nvar, lambda)
{
    ##=########################################################################
    ##
    ## Generate the pseudo design matrix for L_oo penalty
    ##
    ##=########################################################################

    ##x <- as.matrix(x, ncol = 1)

    ks <- 3 ## <==> degree == 2
    nk <- length(knots) + 2*(ks - 1)# = length(new.knots)
    ncoef <- nk - ks
    nd <- nk - 5
    ox <- order(x)
    sortx <- x[ox]
    nrql1 <- nrq + nl1
    nrleq <- nrql1 + neqc
    nobs <- nrleq + niqc
    X <- matrix(0, nobs, nvar)
    z1 <- .splBasis(ord = ks, knots, ncoef, xo = sortx)
    idx1 <- cbind(rep(ox, rep(ks, nrq)),
                  c(outer(1:ks, z1$offsets, "+")))
    X[idx1] <- t(t(z1$design) * rep(w[ox], ks))
    ##
    ## formulate the inequality constraints -- s''()
    ##
    if(lambda != 0) {
        niqc1 <- 2*nd
        z2 <- .splBasis(ord = ks, knots, ncoef, xo = knots[1:nd],
                        derivs = rep(2, nd))
        X[(nrq+1):nrql1,] <- cbind(c(rep(0,nvar-1), lambda))
        idx2 <- cbind(rep(nrleq + 1:nd, rep(ks,nd)),
                      c(outer(1:ks,z2$offsets,"+")))
        X[idx2] <- z2$design
        idx2[,1] <- idx2[,1]+nd
        X[idx2] <- -z2$design
        ## assign different weight to roughness penalty
        X[nrleq + 1:niqc1, ] <-
            X[nrleq + 1:niqc1, ] * rep(pw, 2)
        X[nrleq + 1:niqc1, nvar ] <- rep(1,niqc1)
    } else niqc1 <- 0
    niqc2 <- n.smaller
    niqc3 <- n.greater
    niqc4 <- niqc - niqc1 -niqc2 -niqc3
    if(constraint != "none") {
        if(constraint == "convex" || constraint == "concave") {
            if(lambda == 0) # z2 not yet above
                z2 <- .splBasis(ord = ks, knots, ncoef, xo = knots[1:niqc4],
                                derivs = rep(2, niqc4))
            idx3 <- cbind(rep(nrleq+niqc1 + 1:niqc4, rep(ks,niqc4)),
                          c(outer(1:ks,z2$offsets,"+")))
            X[idx3] <- if(constraint == "convex") z2$design else -z2$design
        }
        else if(constraint == "periodic") {
            neqc3 <- 2
            z2 <- .splBasis(ord = ks, knots, ncoef,
                            xo = sortx[c(1,nrq)], derivs = c(1,1))
            idx2 <- cbind(rep(1:neqc3, rep(ks,neqc3)),
                          c(outer(1:ks,z2$offsets,"+")))
            X.temp <- matrix(0,neqc3,nvar)
            X.temp[idx2] <- z2$design
            X[nrleq,  ] <- X.temp[2,]-X.temp[1,]
            X[nrleq-1,] <- X[ox[nrq],]-X[ox[1],]
        }
        else {
            z3 <- .splBasis(ord = ks, knots, ncoef, xo = knots[1:niqc4],
                            derivs = rep(1, niqc4))
            idx3 <- cbind(rep(nrleq+niqc1 + 1:niqc4, rep(ks,niqc4)),
                          c(outer(1:ks,z3$offsets,"+")))
            X[idx3] <- if(constraint == "increase") z3$design else -z3$design
        }
    }
    if(n.smaller > 0) {
        o.smaller <- order(smaller[,2])
        smaller.o <- smaller[,2][o.smaller]
        z3.1 <- .splBasis(ord = ks, knots, ncoef, xo = smaller.o)
        idx3.1 <- cbind(rep(nrleq+niqc1+niqc4 + o.smaller,rep(ks,niqc2)),
                        c(outer(1:ks,z3.1$offsets,"+")))
        X[idx3.1] <- -z3.1$design
    }
    if(n.greater > 0) {
        o.greater <- order(greater[,2])
        greater.o <- greater[,2][o.greater]
        z3.2 <- .splBasis(ord = ks, knots, ncoef, xo = greater.o)
        idx3.2 <- cbind(rep(nrleq+niqc1+niqc4+niqc2 + o.greater, rep(ks,niqc3)),
                        c(outer(1:ks,z3.2$offsets,"+")))
        X[idx3.2] <- z3.2$design
    }
    ##
    ## formulate the equality constraints
    ##
    neqc1 <- n.equal
    neqc2 <- n.gradient
    if(n.equal > 0) {
        derivs <- rep(0, neqc1)
        o.equal <- order(equal[,2])
        equal.o <- equal[,2][o.equal]
        z1.1 <- .splBasis(ord = ks, knots, ncoef, xo = equal.o)
        idx1.1 <- cbind(rep(nrql1 + o.equal, rep(ks, neqc1)),
                        c(outer(1:ks, z1.1$offsets,"+")))
        X[idx1.1] <- z1.1$design
    }
    if(n.gradient > 0) {
        o.gradient <- order(gradient[,2])
        gradient.o <- gradient[,2][o.gradient]
        z1.1 <- .splBasis(ord = ks, knots, ncoef, xo = gradient.o,
                          derivs = rep(1, neqc2))
        idx1.2 <- cbind(rep(nrql1+neqc1 + o.gradient, rep(ks, neqc2)),
                        c(outer(1:ks, z1.2$offsets,"+")))
        X[idx1.2] <- z1.2$design
    }
    return(X)
} ## loo.design
