################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen partiel informatique 2017 ################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())

VaR1 <- function(u) 1000 * (-log(u))^(-1/2)
VaR2 <- function(u) -2000 * log(2 - 2^u)

VaRcx <- function(u, VaR.f) (
    1/4 * VaR.f(1 - (1 - u)/4) +
    1/2 * VaR.f(1 - (1 - u)/2) +
    1/4 * VaR.f(1 - 3 * (1 - u)/4)
)
VaRS <- function(u) VaR1(u) + VaR2(u)

## e)
VaRcx.res <- c(VaRcx(0.9, VaR1), VaRcx(0.9, VaR2))
c(VaRcx.res, sum(VaRcx.res))

## f)
Fs <- function(x) optimize(function(u) abs(VaRS(u) - x), c(0, 1))$minimum
Fs(5000)

###
### Question 2      À TERMINER!!
###

rm(list = ls())

m <- 1E6; set.seed(2017)
sig <- c(1, 0.8); mu <- log(c(100, 200)) - sig^2/2; rho <- 0.5

## a)
EspX12 <- exp(sum(mu) + (sum(sig^2) + 2 * rho * prod(sig))/2)
CovX12 <- EspX12 - prod(exp(mu + sig^2/2))
rhoX12 <- CovX12/sqrt(prod((exp(2 * mu + sig^2) * (exp(sig^2) - 1))))

## b)
EspS <- sum(exp(mu + sig^2/2))
VarS <- sum(exp(2 * mu + sig^2) * (exp(sig^2) - 1)) + 2 * CovX12

## d)
U <- matrix(runif(2 * m), ncol = 2, byrow = TRUE)
Z <- matrix(qnorm(U), ncol = 2)

X1 <- exp(mu[1] + sig[1] * Z[, 1])
X2 <- exp(mu[2] + sig[2] * (rho * Z[, 1] + sqrt(1 - rho^2) * Z[, 2]))

# i)
c(U[2, ], Z[2, ], X1[2], X2[2])

# ii) et iii)
S <- X1 + X2

VaR <- function(u, X) sort(X)[u * m]
sapply(list(X1, X2, S), VaR, u = 0.99)

VaRcx <- function(u, X) (
    1/4 * VaR(u + (1 - u)/4, X) +
    1/2 * VaR(u + (1 - u)/2, X) +
    1/4 * VaR(u + 3 * (1 - u)/4, X)
)
sapply(list(X1, X2, S), VaRcx, u = 0.99)

TVaR <- function(u, X) mean(X[X > VaR(u, X)])
sapply(list(X1, X2, S), TVaR, u = 0.99)

# iv) voir comment approximer la VaR à 0 ?

# v)
BM <- function(u) sum(sapply(list(X1, X2), TVaR, u = u)) - TVaR(u, S)
BM(0.999)

###
### Question 3
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
g0 <- 0.2; al <- c(0.8, 0.3); m <- 10; q <- c(0.2, 0.3)

## d)
fgpB1 <- function(t) (1 - q[1] + q[1] * t)^m
fgpB2 <- function(t) (1 - q[2] + q[2] * t)^m

fgpM1M2 <- function(t1, t2) exp(
    g0 * (t1 * t2 - 1) + al[1] * (t1 - 1) + al[2] * (t2 - 1)
)
fgpS <- function(t) fgpM1M2(fgpB1(t), fgpB2(t))

# i)
fs <- Re(fft(sapply(e, fgpS)))/n
fs[c(0, 1, 5, 10) + 1]

# ii)
VaR <- function(u) k[min(which(cumsum(fs) >= u))]
sapply(c(0.9, 0.99, 0.999), VaR)


###
### Question 4      À FAIRE!!
###

rm(list = ls())


###
### Question 5
###

rm(list = ls())

al <- 1.2; lam <- 2
VaRX <- function(u) lam * ((1 - u)^(-1/al) - 1)
VaRS <- function(u) VaRX((1 - u)/2) + VaRX((1 + u)/2)

# attention, je n'ai pas exactement la même réponse que dans le corrigé
Fs <- function(x) optimize(function(u) abs(VaRS(u) - x), c(0, 1))$minimum
Fs(100)

###
### Question 6
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
lam <- c(-log(0.01)/3, -log(0.01)/2); lamN <- sum(lam); q <- c(0.2, 0.4)

## b)
1/lam

## c)
fgpC <- function(t) sum(lam/lamN * q * t/(1 - (1 - q) * t))
fc <- Re(fft(sapply(e, fgpC)))/n
fc[c(1, 2) + 1]

## d)
1/lamN

## e)
fgpS <- function(t) exp(4 * lamN * (fgpC(t) - 1))
fs <- Re(fft(sapply(e, fgpS)))/n
fs[c(50, 60, 70) + 1]

Fs <- cumsum(fs)
1 - Fs[105 + 1]

## f)
VaR <- function(u) k[min(which(Fs >= u))]
VaR(0.9)
