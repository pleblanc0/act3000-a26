################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen partiel informatique 2012 ################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); h <- 1
lam <- 0.64; mu <- 2.08; sig <- 0.83

## a)
fb <- c(0, diff(plnorm(k, mu, sig)))
fs <- Re(fft(exp(lam * (fft(fb) - 1)), TRUE))/n
fs[c(0, 50) + 1]

## b)
Fs <- cumsum(fs)
1 - Fs[50 + 1]

## c)
VaR <- function(u) k[min(which(Fs >= u))]
VaR(0.99)


###
### Question 2
###

rm(list = ls())

# Utilisation de l'approximation Poisson composée


###
### Question 3
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); i <- 1:2
r <- 0.75 * i; q <- 0.2 * i; be <- 1/100
e <- exp(-2i * pi * k/n)

## b)
fgp <- function(t) prod((q/(1 - (1 - q) * t^i))^r)
nu <- Re(fft(sapply(e, fgp), TRUE))/n
nu[0:6 + 1]

## c)
Fs <- function(x) nu[1] + sum(nu[-1] * pgamma(x, k[-1], be))
Fs(1000)


###
### Question 4
###

rm(list = ls())

m <- 5; q <- c(0.4, 0.2); s <- 0:(2 * m); t <- -m:m

## a)
Fx12 <- function(k1, k2) min(pbinom(k1, m, q[1]), pbinom(k2, m, q[2]))
fx12 <- function(k1, k2) (
    Fx12(k1, k2) - Fx12(k1 - 1, k2) - Fx12(k1, k2 - 1) + Fx12(k1 - 1, k2 - 1)
)

fs <- sapply(s, function(k) sum(sapply(0:k, function(j) fx12(j, k - j))))
ft <- sapply(t, function(k) sum(sapply(0:m, function(j) fx12(j + k, j))))

## b)
Fx12 <- function(k1, k2) max(pbinom(k1, m, q[1]) + pbinom(k2, m, q[2]) - 1, 0)
fx12 <- function(k1, k2) (
    Fx12(k1, k2) - Fx12(k1 - 1, k2) - Fx12(k1, k2 - 1) + Fx12(k1 - 1, k2 - 1)
)

fs <- sapply(s, function(k) sum(sapply(0:k, function(j) fx12(j, k - j))))
ft <- sapply(t, function(k) sum(sapply(0:m, function(j) fx12(j + k, j))))
