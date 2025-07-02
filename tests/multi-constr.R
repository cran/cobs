#### Examples which use the new feature of more than one 'constraint'.

suppressMessages(library(cobs))

## do *not* show platform info here (as have *.Rout.save), but in 0_pt-ex.R
options(digits = 6)

if(!dev.interactive(orNone=TRUE)) pdf("multi-constr.pdf")

source(system.file("util.R", package = "cobs"))
source(system.file(package="Matrix", "test-tools-1.R", mustWork=TRUE))
##--> tryCatch.W.E(), showProc.time(), assertError(), relErrV(), ...
 Lnx  <- Sys.info()[["sysname"]] == "Linux"
isMac <- Sys.info()[["sysname"]] == "Darwin"
x86 <- (arch <- Sys.info()[["machine"]]) == "x86_64"
noLdbl <- (.Machine$sizeof.longdouble <= 8) ## TRUE when --disable-long-double
## IGNORE_RDIFF_BEGIN
Sys.info()
noLdbl
## IGNORE_RDIFF_END


Rsq <- function(obj) {
    stopifnot(inherits(obj, "cobs"), is.numeric(res <- obj$resid))
    1 - sum(res^2)/obj$SSy
}
list_ <- function (...) `names<-`(list(...), vapply(sys.call()[-1L], as.character, ""))
is.cobs <- function(x) inherits(x, "cobs")

set.seed(908)
x <- seq(-1,2, len = 50)
f.true <- pnorm(2*x)
y <- f.true + rnorm(50)/10
plot(x,y); lines(x, f.true, col="gray", lwd=2, lty=3)

## constraint on derivative at right end:
(con <- rbind(c(2 , max(x), 0))) # f'(x_n) == 0

## Using 'trace = 3' --> 'trace = 2' inside drqssbc2()
##
## Regression splines (lambda = 0)
c2   <- cobs(x,y,                               trace = 3) # gave warnings (till early 2025)
c2i  <- cobs(x,y, constraint = "increase",      trace = 3)
c2c  <- cobs(x,y, constraint = "concave" ,      trace = 3)
c2IC <- cobs(x,y, constraint=c("inc","concav"), trace = 3)
## here, it *was* the same as just "i":
## IGNORE_RDIFF_BEGIN
all.equal(fitted(c2i), fitted(c2IC)) ## (2024-12) no longer ?!????
## IGNORE_RDIFF_END

c1   <- cobs(x,y, degree = 1,                          trace = 3) # gave warnings (no longer in 2025-06)
c1i  <- cobs(x,y, degree = 1, constraint = "increase", trace = 3)
c1c  <- cobs(x,y, degree = 1, constraint = "concave" , trace = 3) # no warnings (2025-06)
## now gives warning (not error):
c1IC <- cobs(x,y, degree = 1, constraint=c("inc","concav"), trace = 3)

plot(c1)
lines(predict(c1i), col="forest green")
## IGNORE_RDIFF_BEGIN
all.equal(fitted(c1), fitted(c1i), tol = 1e-9)# but not 1e-10
## (2024-12:-- now mean rel.diff. 0.0215671 <--> IGNORE)
## IGNORE_RDIFF_END


cp2  <- cobs(x,y,                          pointwise = con, trace = 3)

## Here, warning ".. 'ifl'.. " on *some* platforms (e.g. Windows 32bit) :
r2i <- tryCatch.W.E( cobs(x,y, constraint = "increase", pointwise = con) )
cp2i <- r2i$value
## IGNORE_RDIFF_BEGIN
r2i$warning
## IGNORE_RDIFF_END
## when plotting it, we see that it gave a trivial constant!!
cp2c  <- cobs(x,y, constraint = "concave",  pointwise = con, trace = 3)

## now gives warning (not error): but no warning on M1 mac -> IGNORE
## IGNORE_RDIFF_BEGIN
cp2IC <- cobs(x,y, constraint = c("inc", "concave"), pointwise = con, trace = 3)
## IGNORE_RDIFF_END
cp1   <- cobs(x,y, degree = 1,                            pointwise = con, trace = 3)
cp1i  <- cobs(x,y, degree = 1, constraint = "increase",   pointwise = con, trace = 3)
cp1c  <- cobs(x,y, degree = 1, constraint = "concave",    pointwise = con, trace = 3)

cp1IC <- cobs(x,y, degree = 1, constraint=c("inc","concav"), pointwise = con, trace = 3)

## Named list of all cobs() results above -- sort() collation order matters for ls() !
(curLC <- Sys.getlocale("LC_COLLATE"))
Sys.setlocale("LC_COLLATE", "C")
cobsL <- mget(Filter(\(nm) is.cobs(.GlobalEnv[[nm]]), ls(patt="c[12p]")),
              envir = .GlobalEnv)
Sys.setlocale("LC_COLLATE", curLC) # reverting

knL <- lapply(cobsL, `[[`, "knots")
str(knL[order(lengths(knL))])

gotRsqrs <- sapply(cobsL, Rsq)
Rsqrs <- c(c1  = 0.95079126, c1IC = 0.92974549, c1c  = 0.92974549, c1i  = 0.95079126,
           c2  = 0.94637437, c2IC = 0.91375404, c2c  = 0.92505977, c2i  = 0.95022829,
           cp1 = 0.9426453, cp1IC = 0.92223149, cp1c = 0.92223149, cp1i = 0.9426453,
           cp2 = 0.94988863, cp2IC= 0.90051964, cp2c = 0.91375409, cp2i = 0.93611487)
## M1 mac   "  =     "     , cp2IC= 0.91704726,  "   =      "    , cp2i = 0.94620178
## noLD     "  =     "     , cp2IC=-0.08244284,  "   =      "    , cp2i = 0.94636815
## ATLAS    "  =     "     , cp2IC= 0.91471729,  "   =      "    , cp2i = 0.94506339
## openBLAS "  =     "     , cp2IC= 0.91738019,  "   =      "    , cp2i = 0.93589404
## MKL      "  =     "     , cp2IC= 0.91765403,  "   =      "    , cp2i = 0.94501205
## Intel    "  =     "     , cp2IC= 0.91765403,  "   =      "    , cp2i = 0.94501205
##                                  ^^^^^^^^^^                            ^^^^^^^^^^
## remove these two from testing, notably for the M1 Mac & noLD .. :
##iR2 <- if(!x86 || noLdbl) setdiff(names(cobsL), c("cp2IC", "cp2i")) else TRUE
## actually everywhere, because of ATLAS, openBLAS, MKL, Intel... :
iR2 <- setdiff(names(cobsL), nR2 <- c("cp2IC", "cp2i"))
## IGNORE_RDIFF_BEGIN
dput(signif(gotRsqrs, digits=8))
all.equal(Rsqrs[iR2], gotRsqrs[iR2], tolerance=0)# 2.6277e-9 (Lnx F 38); 2.6898e-9 (M1 mac)
## c1 and c2 changed (Fedora new gcc/clang, quantreg 5.99.1 Date/Publication: 2024-11-22)
## BUG ?? (in quantreg / cobs / Fortran / C/ .... ?) ____ TODO: *Find* out (c1 R^2 is *higher*, 2nd one is lower ..)
##            vvvvvv
## c(c1 = 0.95341697, c1IC = 0.92974549, c1c = 0.92974549, c1i = 0.95079126,
##   c2 = 0.94864721, c2IC = 0.91375404, c2c = 0.92505977, c2i = 0.95022829,
##            ^^^^^^
all.equal(Rsqrs[nR2], gotRsqrs[nR2], tolerance=0)# differ; drastically only for 'noLD'
## IGNORE_RDIFF_END
stopifnot(exprs = {
    all.equal(Rsqrs[iR2], gotRsqrs[iR2], tolerance = 0.0006)
    identical(c(5L, 3L, 3L, 5L,
                3L, 2L, 3L, 4L,
                5L, 3L, 3L, 5L,
                4L, 2L, 2L, 4L), unname(lengths(knL)))
})

plot(x,y, main = "cobs(*, degree= 1, constraint = *, pointwise= *)")
matlines(x,cbind(fitted(c1),
                 fitted(c1i),
                 fitted(c1c),
                 fitted(cp1),
                 fitted(cp1i),
                 fitted(cp1c)),
        col = 1:6, lty=1)
legend("bottomright", inset = .02, col = 1:6, lty=1,
       legend = c("none", "increase","concave",
       "pt", "pt + incr.", "pt + conc."))

if(dev.interactive()) x11() # cheap way to get plot in new window, when testing

plot(x,y, main = "cobs(*, degree= 2, constraint = *, pointwise= *)")
matlines(x,cbind(fitted(c2),
                 fitted(c2i),
                 fitted(c2c),
                 fitted(cp2),
                 fitted(cp2i),
                 fitted(cp2c)),
        col = 1:6, lty=1)
legend("bottomright", inset = .02, col = 1:6, lty=1,
       legend = c("none", "increase","concave",
       "pt", "pt + incr.", "pt + conc."))

##--> "increase + pointwise" gives constant which seems plain wrong  <<<< BUG ???
