## The concept of namespace [in many programming
## languages, not just R] is used to "protect" the
## package functions from users who (accidentally) define a
## function with the name as one of the "system" functions.

## In the case of package *DEVELOPMENT*,
## the protection of course is undesirable.

## I'd do the following (execute the following R code):

library(cobs) # load the unchanged packge and its C/Fortran...

## export to of the NAMESPACE_hidden symbols that you later need:
spline_value <- cobs:::spline_value
spline_basis <- cobs:::spline_basis

example('source') # defines the sourceDir() utility function
if(identical(c(USER = "maechler"), Sys.getenv("USER"))) {
    sourceDir("~/R/Pkgs/cobs/R")
} else sourceDir("..../cobs/R")
##                ^^^^ replace by correct path

## now test to see that things work:
example(cobs, package="cobs")

## Now edit cobs function, and they will be used instead of the "cobs ones"
