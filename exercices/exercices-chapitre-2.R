################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions exercices - Chapitre 2 ##########################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

#===============================================================================
#== Exercices informatiques ====================================================
#===============================================================================

###
### Exercice 1
###

rm(list = ls())

be <- c(1/10, 1/20)

VaR1 <- function(u) qgamma(u, 2, be[1])
VaR1(0.99)

ci <- sapply(1:2, function(i) be[-i]/(be[-i] - be[i]))
Fx2 <- function(x) sum(ci * pexp(x, be))

VaR2 <- function(u) optimize(function(x) abs(Fx2(x) - u), c(0, 1000))$minimum
VaR2(0.99)


###
### Exercice 2
###

rm(list = ls())

be <- c(1/2, 1/3, 1/5)

ci <- sapply(1:3, function(i) prod(be[-i]/(be[-i] - be[i])))
Fr <- function(x) sum(ci * pexp(x, be))

VaRr <- function(u) optimize(function(x) abs(Fr(x) - u), c(0, 100))$minimum
VaRr(0.999)

VaRs <- function(u) qgamma(u, 3, be[1])
VaRs(0.999)

VaRt <- function(u) qgamma(u, 3, be[3])
VaRt(0.999)


###
### Exercice 3
###

rm(list = ls())

sig <- c(0.5, 1); a <- 10; mu <- log(a) - sig^2/2

# i)
VarX1etX2 <- (exp(sig^2) - 1) * exp(2 * mu + sig^2)
VarX1etX2

# ii)
c <- a * exp(prod(sig)/2)
u <- c(0.01, 0.5, plnorm(c, mu[1], sig[1]), 0.99)

# iii)
VaR1 <- function(u) qlnorm(u, mu[1], sig[1])
sapply(u, VaR1)

VaR2 <- function(u) qlnorm(u, mu[2], sig[2])
sapply(u, VaR2)

# iv)
TVaR1 <- function(u) exp(mu[1] + sig[1]^2/2)/(1 - u) * pnorm(
    (mu[1] + sig[1]^2 - log(VaR1(u)))/sig[1]
)
sapply(u, TVaR1)

TVaR2 <- function(u) exp(mu[2] + sig[2]^2/2)/(1 - u) * pnorm(
    (mu[2] + sig[2]^2 - log(VaR2(u)))/sig[2]
)
sapply(u, TVaR2)


###
### Exercice 4
###

rm(list = ls())

al <- 0.5; be <- al/10; n <- c(1, 10, 100, 1000)
alW <- n * al; beW <- n * be
u <- c(0.01, 0.5, 0.99)

# i)
EspW <- alW/beW
EspW

# ii)
VarW <- alW/beW^2
VarW

## iii)
VaRW <- function(k) qgamma(k, alW, beW)
sapply(u, VaRW)

# iv)
TVaRW <- function(k) alW/beW * pgamma(VaRW(k), alW + 1, beW, low = F)/(1 - k)
sapply(u, TVaRW)


###
### Exercice 5
###

rm(list = ls())

n <- 10; lam <- a <- 2; q <- a/n
k <- 0:100; u <- c(0.01, 0.5, 0.99)

# ii)
VarM <- n * q * (1 - q)
VarN <- a

c(VarM, VarN)

# iii)
VaRM <- function(u) qbinom(u, n, q)
sapply(u, VaRM)

VaRN <- function(u) qpois(u, a)
sapply(u, VaRN)

# iv)
TVaRM <- function(u) VaRM(u) + 1/(1 - u) * sum(
    pmax(k - VaRM(u), 0) * dpois(k, a)
)
sapply(u, TVaRM)

TVaRN <- function(u) VaRN(u) + 1/(1 - u) * sum(
    pmax(k - VaRN(u), 0) * dbinom(k, n, q)
)
sapply(u, TVaRN)


###
### Exercice 6
###

rm(list = ls())

lam <- c(2, 5); u <- c(0.01, 0.5, 0.99)

sapply(lam, function(l) qpois(u, l))
