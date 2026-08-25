################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen final informatique 2016 ##################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 3
###

rm(list = ls())

n <- 2^7; k <- 0:(n - 1)
al <- 5; lam <- c(2, 3)

gumbel <- function(u1, u2) exp(-((-log(u1))^al + (-log(u2))^al)^(1/al))
Fm1m2 <- function(k1, k2) gumbel(ppois(k1, lam[1]), ppois(k2, lam[2]))

## b)
fm1m2 <- function(k1, k2) (
    Fm1m2(k1, k2) - Fm1m2(k1 - 1, k2) -
    Fm1m2(k1, k2 - 1) + Fm1m2(k1 - 1, k2 - 1)
)
c(fm1m2(2, 3), fm1m2(3, 2))
probs <- outer(k, k, fm1m2)

## d)
CovM1M2 <- sum(outer(k, k) * probs) - lam[1] * lam[2]
CovM1M2

RhoM1M2 <- CovM1M2/sqrt(lam[1] * lam[2])
RhoM1M2

## f)
sum(outer(pmax(k - 4, 0), I(k > 5)) * probs)


###
### Question 4
###

rm(list = ls())


###
### Question 5
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
lam <- c(-log(0.01)/3, -log(0.01)/2); lamN <- sum(lam); q <- c(0.2, 0.4)

## d)
fgpC <- function(t) sum(lam/lamN * q * t/(1 - (1 - q) * t))
fc <- Re(fft(sapply(e, fgpC), TRUE))/n

## f)
fgpS <- function(t) exp(4 * lamN * (fgpC(t) - 1))
fs <- Re(fft(sapply(e, fgpS), TRUE))/n

Fs <- cumsum(fs)
1 - Fs[105 + 1]

## g)
VaR <- function(u) k[min(which(Fs >= u))]
VaR(0.9)


###
### Question 6
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)

## a)
fgp_n <- function(t) (0.8 + 0.16 * t + 0.04 * t^2)^12
fn <- Re(fft(fgp_n(e), TRUE))/n
fn[5 + 1]

## b)
psi1 <- function(t)
    t * 12 * (0.8 + 0.1 * t + 0.06 * t + 0.04 * t^2)^11 * (0.1 + 0.04 * t)

em1 <- Re(fft(psi1(e), TRUE))/n
em1[5 + 1]

## c)
psi2 <- function(t)
    t * 12 * (0.8 + 0.1 * t + 0.06 * t + 0.04 * t^2)^11 * (0.06 + 0.04 * t)

em2 <- Re(fft(psi2(e), TRUE))/n
em2[5 + 1]

## d)
VaR <- function(u) k[min(which(cumsum(fn) >= u))]
VaR(0.9)

CVaR1 <- function(u) em1[VaR(u) + 1]/fn[VaR(u) + 1]
CVaR1(0.9)

CVaR2 <- function(u) em2[VaR(u) + 1]/fn[VaR(u) + 1]
CVaR2(0.9)
