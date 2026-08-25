################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen partiel informatique 2015 ################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); h <- 1
lamPoi <- 1; r <- 2; q <- 2/3
be <- 1/10; al <- 6; lamPar <- 50

## a)
fb <- c(0, diff(pexp(k, be)))
fb[10 + 1]

fc <- c(0, diff(1 - (lamPar/(lamPar + k * h))^al))
fc[10 + 1]

## b)
fx <- Re(fft(exp(lamPoi * (fft(fb) - 1)), TRUE))/n
Fx <- cumsum(fx)
Fx[c(0, 10, 20, 80) + 1]

## c)
fy <- Re(fft((q/(1 - (1 - q) * fft(fc)))^r, TRUE))/n
Fy <- cumsum(fy)
Fy[c(0, 10, 20, 80) + 1]

## d)
fs <- Re(fft(fft(fx) * fft(fy), TRUE))/n
Fs <- cumsum(fs)
Fs[c(0, 10, 20, 80) + 1]

## e)
VaR <- function(u, F) k[min(which(F >= u))]
sapply(list(Fx, Fy, Fs), VaR, u = 0.95)


###
### Question 2 (Q4 dans le corrigé)
###

## a)
rm(list = ls())
n <- 10; i <- 1:n; al <- 1.5 + 0.1 * i; lam <- 100 * (al - 1)

# ii)
VaR <- function(u) sum(lam * ((1 - u)^(-1/al) - 1))
VaR(0.99)

# iii)
Fs <- function(x) optimize(function(u) abs(VaR(u) - x), c(0, 1))$minimum
Fs(2000)

## b)
rm(list = ls())
al <- 2.5; lam <- 150

# i)
2 * qpareto(0.5, al, lam)

# iii)
VaR <- function(u) qpareto((1 - u)/2, al, lam) + qpareto((1 + u)/2, al, lam)
VaR(0.99)

# iv)
Fs <- function(x) optimize(function(u) abs(VaR(u) - x), c(0, 1))$minimum
Fs(1000)

## c)
rm(list = ls())
r0 <- 0.05; rn <- 0.1; be <- 1/10; n <- 10; k <- 0:n

Fs <- function(x) (
    (1 - r0) * sum(dbinom(k, n, rn) * pgamma(x, k, be)) + r0 * pgamma(x, n, be)
)
Fs(100)


###
### Question 3
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
a1 <- 1.2; a2 <- 1.9; a3 <- 2.5; a12 <- 0.6; a123 <- 0.2

## b)
lamN <- a1 + a2 + a3 + a12 + a123
fgpC <- function(t) ((a1 + a2 + a3)/lamN * t + a12/lamN * t^2 + a123/lamN * t^3)

fc <- Re(fft(sapply(e, fgpC), TRUE))/n
fc[1:3 + 1]

## c)
fgpN <- function(t1, t2, t3) exp(
    a1 * (t1 - 1) + a2 * (t2 - 1) + a3 * (t3 - 1) +
    a12 * (t1 * t2 - 1) + a123 * (t1 * t2 * t3 - 1)
)

fn <- Re(fft(sapply(e, function(t) fgpN(t, t, t)), TRUE))/n
fn[0:10 + 1]

## d)
VaR <- function(u) k[min(which(cumsum(fn) >= u))]
sapply(c(0.1 * 5:9, 0.95, 0.99, 0.995, 0.999, 0.9999), VaR)


###
### Question 4 (Q2 dans le corrigé)
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
p00 <- 0.6; p01 <- p10 <- 0.10; p11 <- 0.20; m <- 20

## a)
fgpM1M2 <- function(t1, t2) (p00 + p10 * t1 + p01 * t2 + p11 * t1 * t2)^m
fgpN <- function(t) fgpM1M2(t^2, t^3)

## c) à faire

## d)
fn <- Re(fft(fgpN(e), TRUE))/n
fn[c(0, 5, 20) + 1]