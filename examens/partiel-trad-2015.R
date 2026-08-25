################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen partiel traditionnel 2015 ################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())

n <- 2^6; k <- 0:(n - 1); e <- exp(2i * pi * k/n)

fgp <- function(t1, t2) (
    0.6 + 0.1 * t2^5 + 0.05 * t1 * t2^4 + 0.05 * t1^2 * t2^3 +
    0.05 * t1^3 * t2^2 + 0.05 * t1^4 * t2 + 0.1 * t1^5
)

## a)
fm1m2 <- Re(fft(outer(e, e, fgp)))/n^2

## b)
fm1 <- Re(fft(sapply(e, function(t) fgp(t, 1))))/n
fm2 <- Re(fft(sapply(e, function(t) fgp(1, t))))/n

## c)
CovM1M2 <- sum(outer(k, k) * fm1m2) - sum(k * fm1) * sum(k * fm2)

## d)
fn <- Re(fft(sapply(e, function(t) fgp(t, t))))/n


###
### Question 2 - il faudrait retravailler cette question
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); m <- 4

fj <- c(0.2, 0.3, 0.5, rep(0, n - 3))
ft <- Re(fft(fft(fj)^m, TRUE))/n

fs <- function(s) ft[s + m + 1]
fs(0)


###
### Question 3
###

rm(list = ls())

al <- 3; lam <- 2000; tau <- 0.5; be <- 1/2000
U <- pexp(3728, 1/1400)

X1 <- lam * ((1 - U)^(-1/al) - 1)
X2como <- (-log((1 - U)^(1/tau)))/be
X2anti <- (-log(U^(1/tau)))/be

c(X1 + X2como, X1 + X2anti)


###
### Question 4      À FAIRE!!
###

rm(list = ls())


###
### Question 5
###

rm(list = ls())

mu <- c(3, 2); sig <- 1

## a)
Fs.como <- function(x) plnorm(x, log(sum(exp(mu))), sig)
Fs.como(40)

## b)
Fs.anti <- function(x)
{
    roots <- Re(polyroot(c(exp(mu[2] - mu[1]), -x * exp(-mu[1]), 1)))
    plnorm(max(roots), 0, sig) - plnorm(min(roots), 0, sig)
}
Fs.anti(40)


###
### Question 6
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
a0 <- 0.2; lam <- c(0.3, 0.5)

fgpM12 <- function(t1, t2) exp(
    a0 * (t1 * t2 - 1) + (lam[1] - a0) * (t1 - 1) + (lam[2] - a0) * (t2 - 1)
)
fgpX12 <- function(t1, t2) fgpM12(0.6 * t1 + 0.4 * t1^2, 0.7 * t2 + 0.3 * t2^2)

## a)
fx12 <- Re(fft(outer(e, e, fgpX12)))/n^2
fx12[1 + 1, 1 + 1]

## c)
fs <- Re(fft(sapply(e, function(t) fgpX12(t, t))))/n
fs[0:5 + 1]

## d)
sum(pmin(k, 2) * fs) * 100


###
### Question 8
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
m <- 2; q <- 3/7; i <- 1:2; w <- c(0.6, 0.4); be <- 1/10

## b)
fgp <- function(t) (1 - q + q * sum(w * t^i))^m
ck <- Re(fft(sapply(e, fgp)))/n
ck[0:4 + 1]

## c)
Fs <- function(x) ck[1] + sum(ck[-1] * pgamma(x, k[-1], be))
Fs(20)

EspTr <- function(d) sum(ck * k/be * pgamma(d, k + 1, be))
EspTr(20)


###
### Question 9
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
fgp <- function(t) (1 - (t - 1))^(-2) * (1 - 2 * (t - 1))^(-1)

## d)
fn <- Re(fft(sapply(e, fgp)))/n
fn[0:2 + 1]
