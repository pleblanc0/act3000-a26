################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen partiel traditionnel 2012 ################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); h <- 1000; kh <- k * h
r <- 0.5; q <- 3/4; al <- 3; lam <- 4000

## a)
fc <- c(0, diff(1 - (lam/(lam + kh))^al))
fb <- c(fc[0:4 + 1], sum(fc[-(0:4 + 1)]))
fb[c(1, 5) + 1]

## c)
fbt <- c(fb, rep(0, n - length(fb)))
fx <- Re(fft((q/(1 - (1 - q) * fft(fbt)))^r, TRUE))/n
fx[c(0, 4, 5) + 1]

## d)
VaR <- function(u) kh[min(which(cumsum(fx) >= u))]
VaR(0.99)


###
### Question 2
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); j <- 1:3
lam <- 1; al <- c(1/3, 1/2, 1/6); theta <- 2^(j - 2)
fb <- c(0, 0.269, 0.178, 0.123, 0.088, 0.065)

panjer.poisson <- function(fb, lam, kmax)
{
    fx <- exp(lam * (fb[1] - 1))
    fb <- c(fb, rep(0, kmax - length(fb) + 1))

    for (i in 1:kmax) fx <- c(
        fx, sum(fb[2:(i + 1)] * fx[i:1] * (lam * (1:i)/i))
    )
    fx
}

fxCond <- sapply(theta, function(t) panjer.poisson(fb, lam * t, 5))
fx <- fxCond %*% al
fx[c(0, 5) + 1]


###
### Question 3
###

rm(list = ls())

n <- 2^7; k <- 0:(n - 1); h <- 10000; kh <- k * h; m <- 5
d <- 0.04; rente <- 2000; prest <- 50000; be <- c(0.05, 0.02)

Fzr <- function(x) ifelse(x <= 50000, 1 - (1 - d * x/rente)^(be[1]/d), 1)
Fza <- function(x) ifelse(x <= 50000, (x/prest)^(be[2]/d), 1)

## b)
fzr <- c(0, diff(Fzr(kh)))
fzr[1:5 + 1]

fza <- c(0, diff(Fza(kh)))
fza[1:5 + 1]

## c)
fsr <- Re(fft(fft(fzr)^m, TRUE))/n
fsr[5:7 + 1]

fsa <- Re(fft(fft(fza)^m, TRUE))/n
fsa[5:7 + 1]

## d)
fs <- Re(fft(fft(fsr) * fft(fsa), TRUE))/n
fs[10:12 + 1]


###
### Question 4
###

rm(list = ls())

al <- 2.5; lam <- 1500; mu <- 5; sig <- 1; be <- 1/2000

TVaRX1 <- function(u) lam * (al/(al - 1) * (1 - u)^(-1/al) - 1)
TVaRX2 <- function(u) exp(mu + sig^2/2) * (1 - pnorm(qnorm(u) - sig))/(1 - u)
TVaRX3 <- function(u) qexp(u, be) + 1/be

## a)
TVaRS <- function(u) TVaRX1(u) + TVaRX2(u) + TVaRX3(u)
TVaRS(0.99)


###
### Question 5
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
i <- 1:2; m <- 3; q <- c(0.4, 0.2)

fgpXi <- function(t, i) (1 - q[i] + q[i] * t^i)^m

## a) On utilise la v.a. St = S + 6
fgpSt <- function(t) fgpXi(t, 1) * fgpXi(1/t, 2) * t^6
fst <- Re(fft(sapply(e, fgpSt)))/n

## b)
VaR <- function(u) (k - 6)[min(which(cumsum(fst) >= u))]
VaR(0.95) * 1000

TVaR <- function(u) sum(pmax((k - 6) - VaR(u), 0) * fst)/(1 - u) + VaR(u)
TVaR(0.95) * 1000


###
### Question 6
###

rm(list = ls())
be <- c(1/5, 1/8); th <- 1

## b)
EspX12 <- 1/prod(be) * (1 + th * (1 - 1/2 - 1/2 + 1/4))
EspX <- 1/be

CovX12 <- EspX12 - prod(EspX)
CovX12

## d)
G <- function(s, g1, g2) ifelse(
    g1 == g2, pgamma(s, 2, be[1]),
    sum(c(g2/(g2 - g1), g1/(g1 - g2)) * pexp(s, c(g1, g2)))
)
Fs <- function(s) (
    (1 + th) * G(s, be[1], be[2]) - th * G(s, 2 * be[1], be[2]) -
    th * G(s, be[1], 2 * be[2]) + th * G(s, 2 * be[1], 2 * be[2])
)
Fs(20)


###
### Question 7
###

rm(list = ls())
al <- 1.5; lam <- 1000; d <- 9000
mu <- c(0.08, 0.085); sig <- c(0.2, 0.3)

Fri <- function(x, i) pnorm(x, mu[i], sig[i])
Sbd <- (lam/(lam + d))^al

# À TERMINER
# Sxd.min
# Sxd.max 
 

###
### Question 8
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
a0 <- 0.4; al <- c(0.6, 1.1); mu <- c(2.3, 2.1); sig <- c(0.8, 0.9)

## a)
CovX1X2 <- a0 * exp(sum(mu) + sum(sig^2)/2)
CovX1X2

## b)
fgpM <- function(t1, t2) exp(
    a0 * (t1 * t2 - 1) + al[1] * (t1 - 1) + al[2] * (t2 - 1)
)
fm12 <- Re(fft(outer(e, e, fgpM)))/n^2
c(fm12[0:1 + 1, 0:1 + 1])

## c)
cond_prob <- (
    fm12[1, 1] +
    fm12[2, 1] * plnorm(10, mu[1], sig[1]) +
    fm12[1, 2] * plnorm(20, mu[2], sig[2]) +
    fm12[2, 2] * plnorm(10, mu[1], sig[1]) * plnorm(20, mu[2], sig[2])
)/sum(fm12[0:1 + 1, 0:1 + 1])
cond_prob


###
### Question 9      À FAIRE!!
###

rm(list = ls())
