library(cobs)
options(digits = 6)
postscript("ex3.ps")

data(women) # 15 obs.
attach(women)

## Interpolation! very easy problem,   BUT :
try( ## gives  "ifl = 5" !!!!!
cobw <- cobs(weight, height, knots = weight, nknots = length(weight))
)
