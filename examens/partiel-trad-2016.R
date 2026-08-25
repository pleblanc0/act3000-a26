################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen partiel traditionnel 2016 ################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################


###
### Question 1
###

## a)
rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
r <- 0.2; q <- 1/2; gam <- 1/4

fgp <- function(t) (q/(1 - (1 - q) * gam * t/(1 - (1 - gam) * t)))^r
fx <- Re(fft(sapply(e, fgp)))/n
fx[0:3 + 1]

## b)
rm(list = ls())

be <- c(1/2, 1); X1 <- 10
U <- pexp(X1, be[1])
X2 <- qexp(1 - U, be[2])

## c)
rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
fgpM12 <- function(t1, t2) (0.56 + 0.14 * t1 + 0.24 * t2 + 0.06 * t1 * t2)^6

fn <- Re(fft(sapply(e, function(t) fgpM12(t, t))))/n
fn[0:1 + 1]


###
### Question 4      À FAIRE!!
###

rm(list = ls())


###
### Question 5
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
q <- 0.7; lam <- 0.6; gam <- c(4, 6); p <- 1/gam
p00 <- 0.1; p10 <- 0.2; p01 <- 0.3; p11 <- 0.4

fgpJ <- function(t) 1 - q + q * t
fgpK <- function(t) exp(lam * (t - 1))

fgpB <- function(t, i) p[i] * t/(1 - (1 - p[i]) * t)
fgpI1I2 <- function(t1, t2) p00 + p10 * t1 + p01 * t2 + p11 * t1 * t2

## b)
EspX <- q * lam * (p10 * gam[1] + p01 * gam[2] + p11 * sum(gam))
EspX

## c)
fgpX <- function(t) fgpJ(fgpK(fgpI1I2(fgpB(t, 1), fgpB(t, 2))))
fx <- Re(fft(sapply(e, fgpX)))/n
fx[0:3 + 1]


###
### Question 6      À FAIRE!!
###

rm(list = ls())


###
### Question 7
###

rm(list = ls())

be <- c(1/2, 1/3, 1/5)
Fxi <- function(x, i) pexp(x, be[i])

Fx123 <- function(x1, x2, x3) {
    max(min(Fxi(x1, 1), Fxi(x3, 3)) + Fxi(x2, 2) - 1, 0)
}
Fx123(4, 4, 4)

###
### Question 9      À FAIRE!!
###

## a)
rm(list = ls())

## b)
rm(list = ls())
