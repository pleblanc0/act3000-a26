################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen final informatique 2015 ##################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())


###
### Question 4
###

## a)
rm(list = ls())
al <- 2.5; lam <- 150

# i)
2 * qpareto(1/2, al, lam)

# iii)
VaR <- function(u) qpareto((1 - u)/2, al, lam) + qpareto((1 + u)/2, al, lam)

Fs <- function(x) optimize(function(u) abs(VaR(u) - x), c(0, 1))$minimum
Fs(1000)

## b)
rm(list = ls())

## c)
rm(list = ls())

al <- 1.5; lam <- 50; be <- 1/100; mu <- log(100) - 1/2; sig <- 1

# ii)
VaR <- function(u) qpareto(u, al, lam) + qexp(u, be) + qlnorm(u, mu, sig)
VaR(0.99)

# iii)
Fs <- function(x) optimize(function(u) abs(VaR(u) - x), c(0, 1))$minimum
Fs(2000)
