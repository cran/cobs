#### --> read ../src/dunif01-tst.c

library(cobs)
postscript("dunif01.ps")

n <- as.integer(100000)
ii <- as.integer(.C("COBS_dunif01",
                    iseed = as.integer(2), n, u = double(n))$u * 99730)
w1 <- which(ii == 78104)
w1
diff(w1) # all = 9972 --> periodicity = 9972 :
all( ii[-(1:9972)] == ii[-((n-9972+1):n)])# proof

i0 <- ii[1:(9972 + 20)] # only first period + 20 obs
i0[1:100]
i0[9971 + 1:20]
library(ts)
lag.plot(i0)# looks fine (they said it passed the spectral test...)
plot(ts(i0[1:200]))# fine
acf(ts(i0[1:1000]))# fine
acf(ts(abs(i0[1:1000] - mean(range(i0)))))# fine
