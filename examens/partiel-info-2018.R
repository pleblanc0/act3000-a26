################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen partiel informatique 2018 ################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1)
sig <- 0.9; mu <- log(5) - sig^2/2
data <- c(2184, 882, 276, 65, 16, 4, 2)

## a)
res <- optim(
    c(1, 0.5), function(p) -sum(data * dnbinom(0:6, p[1], p[2], log = TRUE))
)
r <- res$par[1]; q <- res$par[2]

fn <- dnbinom(k, r, q)
fn[0:5 + 1]

## b)
fxl <- c(0, diff(plnorm(k, mu, sig)))
fxl[0:3 + 1]

fxu <- c(diff(plnorm(k, mu, sig)), 0)
fxu[0:3 + 1]

fsl <- Re(fft(sapply(fft(fxl), function(t) (q/(1 - (1 - q) * t))^r), TRUE))/n
fsl[c(10, 20, 30) + 1]

fsu <- Re(fft(sapply(fft(fxu), function(t) (q/(1 - (1 - q) * t))^r), TRUE))/n
fsu[c(10, 20, 30) + 1]

Fsl <- cumsum(fsl)
Fsl[c(10, 20, 30) + 1]

Fsu <- cumsum(fsu)
Fsu[c(10, 20, 30) + 1]

## c)
VaR <- function(u, F) k[min(which(F >= u))]
c(VaR(0.995, Fsl), VaR(0.995, Fsu))


###
### Question 2
###

rm(list = ls())

data <- c(
    2, 1, 2, 1, 3, 3, 1, 1, 5, 8, 5, 9, 8, 8, 3, 5, 5, 6,
    7, 5, 9, 4, 4, 4, 7, 8, 10, 11, 14, 10, 9, 7, 5, 14
)
a <- 1; b <- 0.2; k.data <- 0:length(data)

## a)
lamT <- function(t, a, b) a * t + b/2 * t^2

## b)
dpois(5, lamT(17, a, b) - lamT(16, a, b))

## c)
nll <- function(p) -sum(dpois(data, diff(lamT(k.data, p[1], p[2])), log = TRUE))

## d)
-nll(c(a, b))

## e)
res <- optim(c(1, 0.1), nll)
res$par; res$value

## f)
lam2020 <- lamT(38, res$par[1], res$par[2]) - lamT(37, res$par[1], res$par[2])
qpois(c(0.01, 0.99), lam2020)


###
### Question 4
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
be <- c(1/12, 1/2); lam <- c(0.03, 0.07); m <- 10; q <- be[1]/be[2]

## b)
lamN <- m * sum(lam)

## c)
p <- m * lam / lamN

## g)
fgp <- function(t) exp(lamN * (p[1] * q * t/(1 - (1 - q) * t) + p[2] * t - 1))
fl <- Re(fft(sapply(e, fgp)))/n
fl[c(0, 5, 10) + 1]

## h)
Fs <- function(x) fl[1] + sum(fl[-1] * pgamma(x, k[-1], be[2]))
sapply(c(0, 10, 20), Fs)
