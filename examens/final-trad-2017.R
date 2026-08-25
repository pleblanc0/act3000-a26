################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen final traditionnel 2017 ##################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)

fgp <- function(t1, t2) (
    0.2 + 0.4 * (0.7 + 0.3 * t1)^2 + 0.3 * (0.8 + 0.2 * t2)^2 +
    0.1 * (0.65 + 0.15 * t1 + 0.05 * t2 + 0.15 * t1 * t2)^2
)

## a)
fx1 <- Re(fft(sapply(e, function(t) fgp(t, 1)), TRUE))/n
fx1[0:2 + 1]

## b)
fx2 <- Re(fft(sapply(e, function(t) fgp(1, t)), TRUE))/n
fx2[0:2 + 1]

## c)
fs <- Re(fft(sapply(e, function(t) fgp(t, t)), TRUE))/n
fs[0:4 + 1]


###
### Question 2
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
r <- 2; m <- 20; q <- c(0.9, 0.05); p <- c(1/2, 1/3); gam <- 3

## a)
fgpX1 <- function(t) (q[1]/(1 - (1 - q[1]) * p[1] * t/(1 - (1 - p[1]) * t)))^r
fx1 <- Re(fft(sapply(e, fgpX1), TRUE))/n
fx1[c(0, 2) + 1]

## b)
fgpX2 <- function(t) (1 - q[2] + q[2] * p[2] * t/(1 - (1 - p[2]) * t))^m
fx2 <- Re(fft(sapply(e, fgpX2), TRUE))/n
fx2[c(0, 2) + 1]

## c)
gumbel <- function(u1, u2) exp(-((-log(u1))^gam + (-log(u2))^gam)^(1/gam))
Fx1x2 <- function(k1, k2) gumbel(cumsum(fx1)[k1 + 1], cumsum(fx2)[k2 + 1])

fx1x2 <- function(k1, k2) (
    Fx1x2(k1, k2) - Fx1x2(k1 - 1, k2) -
    Fx1x2(k1, k2 - 1) + Fx1x2(k1 - 1, k2 - 1)
)
fx1x2(1, 2)


###
### Question 3      À TERMINER!!
###

rm(list = ls())

al <- cbind(c(1, 3, 2), c(1.5, 2.5, 3.5))
be <- cbind(c(1/2, 1/4, 1/3), c(1/5, 1/6, 1/7))
gam <- c(0.1, 0.2, 0.3, 0.4)
vk <- 3.3674

EspX3j <- 1/be[3, ] * (al[3, ] + sum(gam[c(0, 3) + 1]))
VarX3j <- 1/be[3, ]^2 * (al[3, ] + sum(gam[c(0, 3) + 1]))

EspS3 <- sum(EspX3j)
VarS3 <- sum(VarX3j) + 2/prod(be[3, ]) * sum(gam[c(0, 3) + 1])

BMP3 <- vk * (sum(sqrt(VarX3j)) - sqrt(VarS3))


###
### Question 4      À FAIRE!!
###

rm(list = ls())

mu <- c(5, 10, 15); sig <- c(1, 2, 3); rho <- matrix(
    c(1, 0.2, 0.6, 0.2, 1, 0.4, 0.6, 0.4, 1), nrow = 3, byrow = TRUE
)

## b)
muS <- sum(mu)
sigS <- sqrt(t(sig) %*% rho %*% sig)

# VÉRIFIER LA SUITE
# ## c)
# TVaRS <- function(u) muS + sigS * dnorm(qnorm(u))/(1 - u)
# TVaRS(0.9999)

# ## d)
# CTVaRS <- function(u) sapply(1:3, function(i)
#     mu[i] + sig[i] * rho[i, ] %*% sig/sigS * dnorm(qnorm(u))/(1 - u)
# )
# CTVaR(0.99)


###
### Question 5      À FAIRE!!
###

rm(list = ls())


###
### Question 7      À FAIRE!!
###

rm(list = ls())

q <- c(0.3, 0.4); be <- c(1/2, 1)

## a)
Fx12 <- function(x1, x2) max(sum(1 - q + q * pexp(be * c(x1, x2))) - 1, 0)

## b) à faire


###
### Question 10     À FAIRE!!
###

rm(list = ls())
