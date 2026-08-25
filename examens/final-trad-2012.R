################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen final traditionnel 2012 ##################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); i <- 1:5
lam <- 0.06 - 0.01 * i; lamN <- sum(lam)

## b)
fb <- c(0, lam/lamN, rep(0, n - length(lam) - 1))

## c)
fs <- Re(fft(exp(lamN * (fft(fb) - 1)), TRUE))/n


###
### Question 2
###

rm(list = ls())

mu <- c(0.06, 0.05); sig <- c(0.2, 0.15)
al <- 4; tau <- 0.5; be <- 1/1000

Fri <- function(x, i) 1/2 * (
    1 + (x - mu[i])/(sig[i] * sqrt(2 + ((x - mu[i])/sig[i])^2))
)
gumbel <- function(u1, u2) exp(-((-log(u1))^al + (-log(u2))^al)^(1/al))

## a)
Fr12 <- function(x1, x2) gumbel(Fri(x1, 1), Fri(x2, 2))
Fr12(-0.1, -0.1)

## b)
q <- Fri(-0.1, 1) + Fri(-0.1, 2) - Fr12(-0.1, -0.1)
q

## d)
Fx <- function(x) 1 - q + q * (1 - exp(-(be * x)^tau))
sapply(c(0, 10000), Fx)

## e)
VaR <- function(u) ifelse(
    u < 1 - q, 0, optimize(function(x) abs(Fx(x) - u), c(0, 1e5))$min
)
sapply(c(0.5, 0.995), VaR)


###
### Question 3
###

rm(list = ls())

fx12 <- rbind(c(0.05, 0.25, 0.4), c(0.01, 0.1, 0.09), c(0.04, 0.05, 0.01))
k <- 0:2 * 1000; fx1 <- rowSums(fx12); fx2 <- colSums(fx12)

## a)
values <- outer(0:2 * 1000, 0:2 * 1000, '+')
fs <- as.vector(tapply(fx12, values, sum))

## a)
EspX1 <- sum(k * fx1)
EspX2 <- sum(k * fx2)
EspX12 <- sum(outer(k, k) * fx12)
CovX12 <- EspX12 - EspX1 * EspX2

s <- 0:4 * 1000
EspS <- sum(s * fs)
VarS <- sum(s^2 * fs) - EspS^2

c(CovX12, EspS, VarS)

## b)
VarX1 <- sum(k^2 * fx1) - EspX1^2
VarX2 <- sum(k^2 * fx2) - EspX2^2

## c)
VaRS <- function(u) s[min(which(cumsum(fs) >= u))]
VaRS(0.9)

TVaRS <- function(u) sum(pmax(s - VaRS(u), 0) * fs)/(1 - u) + VaRS(u)
TVaRS(0.9)

## d)
EspX1 + (VarX1 + CovX12)/VarS * (TVaRS(0.9) - EspS)
EspX2 + (VarX2 + CovX12)/VarS * (TVaRS(0.9) - EspS)

## e) à terminer


###
### Question 9      À FAIRE!!
###

rm(list = ls())

b <- c(1000, 2000); be <- c(0.05, 0.06); al <- 1; d <- 0.04

fgm.exp <- function(t, be) be/(be - t)
fgmT12 <- function(t1, t2) (
    (1 + al) * fgm.exp(t1, be[1]) * fgm.exp(t2, be[2]) +
    al * fgm.exp(t1, 2 * be[1]) * fgm.exp(t2, 2 * be[2]) -
    al * fgm.exp(t1, 2 * be[1]) * fgm.exp(t2, be[2]) -
    al * fgm.exp(t1, be[1]) * fgm.exp(t2, 2 * be[2])
)
Ft12 <- function(x1, x2) prod(pexp(c(x1, x2), be)) * (
    1 + al * prod(pexp(c(x1, x2), be, lower = FALSE))
)

## b)
EspZ12 <- prod(b) * fgmT12(-d, -d)
CovZ12 <- EspZ12 - 555.56 * 1200

## c) à vérifier, on dirait que le solutionnaire a fait une erreur?
Fz12 <- function(x1, x2)
{
    gam <- c(-1/d * log(c(x1, x2)/b))
    1 - Ft12(gam[1], Inf) - Ft12(Inf, gam[2]) + Ft12(gam[1], gam[2])
}
Fz12(700, 1500)