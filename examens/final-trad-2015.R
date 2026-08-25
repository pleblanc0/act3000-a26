################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen final traditionnel 2015 ##################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
g0 <- 0.2; gam <- c(0.4, 0.3)

# avec une fgp composée
lamS <- g0 + sum(gam); v <- 1:3
eta <- c(gam/lamS, g0/lamS)

fgp <- function(t) exp(lamS * (sum(eta * t^v) - 1))
fs <- Re(fft(sapply(e, fgp), TRUE))/n
fs[0:2 + 1]

# directement la fgp de S
fgp.direct <- function(t1, t2) exp(
    g0 * (t1 * t2 - 1) + gam[1] * (t1 - 1) + gam[2] * (t2 - 1)
)
fs.direct <- Re(fft(sapply(e, function(t) fgp.direct(t, t^2)), TRUE))/n
fs.direct[0:2 + 1]


###
### Question 2
###

rm(list = ls())

al <- 0.5
copule <- function(u1, u2) al * min(u1, u2) + (1 - al) * max(u1 + u2 - 1, 0)
Fs <- function(x) ifelse(x < 0, al * pnorm(x/2), 1 - al + al * pnorm(x/2))

# iii)
VaRS <- function(u) optimize(function(x) abs(Fs(x) - u), c(-10, 10))$min
sapply(c(0.01, 0.99), VaRS)

# iv)
Fx12 <- function(x1, x2) copule(pnorm(x1), pnorm(x2))
c(Fx12(1, 1), pnorm(1)^2)


###
### Question 5      À TERMINER!!
###

## a)
rm(list = ls())
be <- 1/100; sig <- 1; mu <- log(100) - sig^2/2

VaRS <- function(u) qexp(u, be) + qlnorm(u, mu, sig)
VaRS(0.95)

## b)
rm(list = ls())

## c)
rm(list = ls())
n <- 2^10; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
fgp <- function(t) (0.7 + 0.2 * t + 0.1 * t^2)^10

fn <- Re(fft(sapply(e, fgp), TRUE))/n
fn[0:1 + 1]


###
### Question 8      À FAIRE!!
###

rm(list = ls())


###
### Question 9
###

rm(list = ls())

mu <- c(3, 2); sig <- 1

## b)
Fs <- function(x)
{
    roots <- Re(polyroot(c(exp(mu[2] - mu[1]), -x * exp(-mu[1]), 1)))
    plnorm(max(roots), 0, sig) - plnorm(min(roots), 0, sig)
}
Fs(40)


###
### Question 10     À FAIRE!!
###

rm(list = ls())
