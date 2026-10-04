### ############################################################################
### Solution du laboratoire informatique #######################################
### Théorie du risque - 4 octobre 2026 #########################################
################################################################################

n <- 2^20; k <- 0:(n - 1); v <- c(10, 50)
r <- c(0.5, 2.5, 5); q <- c(1/11, 1/3, 1/2)

EspX <- r * (1 - q)/q

fx <- sapply(1:3, function(i) dnbinom(k, r[i], q[i]))
fs <- Re(fft(apply(mvfft(fx), 1, prod), TRUE))/n

# a)
Fs <- cumsum(fs)
Fs[v + 1]

# b)
re <- 10; la <- 40
prime <- sum(pmin(pmax(k - re, 0), la) * fs)

# c)
VaR <- function(u) k[min(which(Fs >= u))]
VaR(0.99)

TVaR <- function(u) sum(pmax(k - VaR(u), 0) * fs)/(1 - u) + VaR(u)
TVaR(0.99)

# d)
EA <- sapply(1:3, function(i) Re(
    fft(fft(k * fx[, i]) * apply(mvfft(fx[, -i]), 1, prod), TRUE)/n
))
round(EA[v + 1, ], 4)
colSums(EA)

# e)
EC <- sapply(1:3, function(i) zapsmall(EA[, i])/fs)
round(EC[v + 1, ], 4)
rowSums(round(EC[v + 1, ], 4))

###
### méthode alternative en utilisant directement les ogfs
###

e <- exp(-2i * pi * k/n)
fgp <- function(t) prod((q/(1 - (1 - q) * t))^r)

fs.alt <- Re(fft(sapply(e, fgp), TRUE))/n
all.equal(fs, fs.alt)

ogf <- function(i, t) t * r[i] * (1 - q[i])/(1 - (1 - q[i]) * t) * fgp(t)
EA.ogf <- sapply(
    1:3, function(i) Re(fft(sapply(e, function(t) ogf(i, t)), TRUE))/n
)
all.equal(EA[v + 1, ], EA.ogf[v + 1, ])


###
### méthode alternative (v2) en identifiant des fgps dans l'ogf
###

# fmp d'une géométrique avec translation
fyi <- sapply(1:3, function(i) c(0, dgeom(0:(n - 2), q[i])))

EA.fgp <- sapply(
    1:3, function(i) EspX[i] * Re(fft(fft(fyi[, i]) * fft(fs), TRUE)/n)
)

all.equal(EA[v + 1, ], EA.fgp[v + 1, ])
