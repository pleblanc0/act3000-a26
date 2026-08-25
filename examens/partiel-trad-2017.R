################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen partiel traditionnel 2017 ################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1      À FAIRE!!
###

rm(list = ls())


###
### Question 2      À FAIRE!!
###

rm(list = ls())


###
### Question 4
###

rm(list = ls())


###
### Question 5
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
fgp <- function(t1, t2) (
    0.4 + 0.3 * (0.9 + 0.1 * t1)^2 + 0.2 * (0.8 + 0.2 * t2)^2 +
    0.1 * (0.9 + 0.1 * t1)^2 * (0.8 + 0.2 * t2)^2
)

## a)
fx1 <- Re(fft(sapply(e, function(t) fgp(t, 1))))/n
fx1[0:2 + 1]

## b)
fx2 <- Re(fft(sapply(e, function(t) fgp(1, t))))/n
fx2[0:2 + 1]

## c)
fs <- Re(fft(sapply(e, function(t) fgp(t, t))))/n
fs[0:4 + 1]


###
### Question 7
###

rm(list = ls())

al <- 3; lam <- 1
EspX <- lam/(al - 1)

VaRX <- function(u) lam * ((1 - u)^(-1/al) - 1)
TVaRX <- function(u) lam * (al/(al - 1) * (1 - u)^(-1/al) - 1)
LTVaRX <- function(u) 1/u * (EspX - (1 - u) * TVaRX(u))

VaRS <- function(u) VaRX((1 - u)/2) + VaRX((1 + u)/2)
VaRS(0.99)

TVaRS <- function(u) LTVaRX((1 - u)/2) + TVaRX((1 + u)/2)
TVaRS(0.99)


###
### Question 10
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
p <- 0.4; r <- 0.2; q <- 1/2; gam <- 1/4

## c)
fgpY <- function(t) (q/(1 - (1 - q) * gam * t/(1 - (1 - gam) * t)))^r
fy <- Re(fft(sapply(e, fgpY)))/n
fy[c(0, 3) + 1]

## d)
fgpX <- function(t) 1 - p + p * fgpY(t)
fx <- Re(fft(sapply(e, fgpX)))/n
fx[c(0, 3) + 1]


###
### Question 11
###

rm(list = ls())
