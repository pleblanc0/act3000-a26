################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen final traditionnel 2016 ##################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())

q <- c(0.2, 0.1); be <- c(1, 1/2); U <- 0.85

Fxi <- function(x, i) 1 - q[i] + q[i] * pexp(be[i] * x)
VaRi <- function(u, i) ifelse(
    u < 1 - q[i], 0, optimize(function(x) abs(Fxi(x, i) - u), c(0, 10))$min
)

## a)
VaRi(U, 1)

## b)
VaRi(U, 2)

## c)
VaRi(U, 1) + VaRi(U, 2)

## d)
VaRi(U, 1) + VaRi(1 - U, 2)


###
### Question 2
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
gam <- c(0.3, 0.25, 0.2, 0.15, 0.1)

lamS <- sum(gam); v <- 1:3
tau <- c(sum(gam[1:3])/lamS, gam[4]/lamS, gam[5]/lamS)

# en utilisant une fgp composée
fgpC <- function(t) sum(tau * t^v)
fgpS <- function(t) exp(lamS * (fgpC(t) - 1))

fs <- Re(fft(sapply(e, fgpS), TRUE))/n
fs[c(0, 3) + 1]

# directement la fgp de S
fgpM <- function(t1, t2, t3) exp(
    gam[1] * (t1 - 1) + gam[2] * (t2 - 1) + gam[3] * (t3 - 1) +
    gam[4] * (t2 * t3 - 1) + gam[5] * (t1 * t2 * t3 - 1)
)

fs.direct <- Re(fft(sapply(e, function(t) fgpM(t, t, t)), TRUE))/n
fs.direct[c(0, 3) + 1]


###
### Question 5
###

## a)
rm(list = ls())
gam <- c(0.5, 0.3, 0.1); be <- 1

VarS <- gam[1] + gam[2] + 4 * gam[3]
sqrt(VarS)

## b)
rm(list = ls())
n <- 2; q <- 0.2; al <- 2; lam <- 10

Fs <- function(x) 1 - q + q * ppareto(x/n, al, lam)
Fs(60)

## c)
rm(list = ls())
fgp <- function(t1, t2) (0.8 + 0.1 * t1 + 0.06 * t2 + 0.04 * t1 * t2)^20
fgp(0, 0)


###
### Question 6
###

rm(list = ls())

mu <- c(10, 5, 15); sig <- c(2, 12, 7); rho <- matrix(
    c(1, -0.9, -0.8, -0.9, 1, -0.7, -0.8, -0.7, 1), nrow = 3, byrow = TRUE
)

## b)
muS <- sum(mu)
sigS <- sqrt(t(sig) %*% rho %*% sig)

## c)
TVaRS <- function(u) muS + sigS * dnorm(qnorm(u))/(1 - u)
TVaRS(0.99)

## d)
CTVaRS <- function(u) sapply(1:3, function(i)
    mu[i] + sig[i] * rho[i, ] %*% sig/sigS * dnorm(qnorm(u))/(1 - u)
)
CTVaR(0.99)


###
### Question 7
###

rm(list = ls())


###
### Question 8      À TERMINER!!
###

## a)
rm(list = ls())
a <- c(0.3, 0.4); b <- c(1.9, 1.7)

Fs <- function(x) (x - (a[1] + b[2]))/((b[1] - a[1]) - (b[2] - a[2]))
Fs(2.1)

## b) à faire
rm(list = ls())


###
### Question 9
###

rm(list = ls())

al <- 2; lam <- 1; tau <- 1/2; be <- 1/2
copule <- function(u1, u2) (
    1 - ((1 - u1)^5 + (1 - u2)^5 - ((1 - u1) * (1 - u2))^5)^(1/5)
)

Fx1 <- function(x) 1 - (lam/(lam + x))^al
Fx2 <- function(x) 1 - exp(-(be * x)^tau)

Fx12 <- function(x1, x2) Fx1(x1) - copule(Fx1(x1), 1 - Fx2(x2))
Fx12(2, 2)
