################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen partiel traditionnel 2018 ################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
r <- 1.5; q <- 0.4; v <- 1:3; p <- c(0.2, 0.5, 0.3)

fgp <- function(t) (q/(1 - (1 - q) * sum(p * t^v)))^r

## b)
fy <- Re(fft(sapply(e, fgp), TRUE))/n
fy[c(0, 1, 2, 10)  + 1]

## c)
EspY <- sum(k * fy)
EspY * 100

# en reconnaissant la fgp
EspY.alt <- r * (1 - q)/q * sum(p * v)
EspY.alt * 100


###
### Question 2
###

rm(list = ls())

## a)
lam <- 39/12

## b)
pexp(1, lam, lower = FALSE)


###
### Question 4
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
m <- 20; v <- 0:4; p <- c(0.4, 0.1, 0.25, 0.15, 0.10)

fgp <- function(t) sum(p * t^v)^m

# i)
fs <- Re(fft(sapply(e, fgp), TRUE))/n
fs[c(44, 45) + 1]

# ii)
Fs <- cumsum(fs)
Fs[43:45 + 1]

# iii)
VaR <- function(u) k[min(which(Fs >= u))]
VaR(0.99)


###
### Question 5
###

rm(list = ls())

k <- 1:10; tau <- 0.5; lam <- 1; al <- 2.5; be <- 1/100; t <- 0.5

## b)
EspSt <- al/be * sum(pgamma(t, tau * k, lam))
EspSt

## c)
fnt <- function(k) pgamma(t, tau * k, lam) - pgamma(t, tau * (k + 1), lam)
sapply(0:2, fnt)

## d)
FstCond <- function(x, k) pgamma(x, al * k, be)
sapply(0:2, function(k) FstCond(500, k))

## e)
v <- 0:100
fn3 <- sapply(
    v, function(v) pgamma(3, tau * v, lam) - pgamma(3, tau * (v + 1), lam)
)

Fn3 <- cumsum(fn3)
v[min(which(Fn3 >= 0.5))]
