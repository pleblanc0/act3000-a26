################################################################################
#### Théorie du risque - Automne 2025 ##########################################
#### Solutions exercices - Chapitre 1 ##########################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

#===============================================================================
#== Exercices traditionnels ====================================================
#===============================================================================

###
### Exercice 1
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
r <- 0.2; q <- 1/2; gam <- 1/4

fgp <- function(t) (q/(1 - (1 - q) * gam * t/(1 - (1 - gam) * t)))^r
fx <- Re(fft(fgp(e))/n)
fx[c(0, 3) + 1]


###
### Exercice 2
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n); m <- c(40, 50, 30)
lam <- 0.04 * m[1]; r <- c(2 * m[2], 3 * m[3]); q <- c(0.97, 0.99)

fgp <- function(t) (exp(lam * (t - 1)) * prod((q/(1 - (1 - q) * t))^r))
fn <- Re(fft(sapply(e, fgp)))/n
fn[0:3 + 1]


###
### Exercice 3
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
m <- 100; r <- 0.2; q <- 1/1.2; gam <- 0.4

fgp <- function(t) (q/(1 - (1 - q) * gam * t/(1 - (1 - gam) * t)))^(r * m)
fs <- Re(fft(fgp(e)))/n
fs[0:3 + 1]


###
### Exercice 4
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
m <- 100; q <- 0.005; w <- c(1/2, 1); p <- c(0.6, 0.4); c <- 2E5

fgp <- function(t) (1 - q + q * sum(p * t^(w * c/1E5)))^m
fs <- Re(fft(sapply(e, fgp)))/n
fs[0:4 + 1]


###
### Exercice 5
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
lam <- c(0.036, 0.054); q <- c(0.4, 0.5); m <- c(120, 80)
lamS <- sum(lam * m); probs <- lam * m/lamS

fgp <- function(t) exp(sum(lam * m * (q * t/(1 - (1 - q) * t) - 1)))
fs <- Re(fft(sapply(e, fgp)))/n
fs[0:3 + 1]

# Solution alternative
fgp.alt <- function(t) exp(lamS * (sum(probs * (q * t/(1 - (1 - q) * t))) - 1))
fs.alt <- Re(fft(sapply(e, fgp)))/n
fs.alt[0:3 + 1]


###
### Exercice 6
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); h <- 1000; kh <- k * h
r <- 0.4; q <- 2/3; tau <- 0.8; be <- 1/1000

## a)
fc <- c(0, diff(1 - exp(-((be * kh)^tau))))
fb <- c(fc[0:4 + 1], sum(fc[6:n]), rep(0, n - 6))
fb[0:5 + 1]

## b)
fgp <- function(t) (q/(1 - (1 - q) * t))^r
fx <- Re(fft(fgp(fft(fb)), TRUE))/n
fx[0:7 + 1]

## c)
VaR <- function(u) kh[min(which(cumsum(fx) >= u))]
VaR(0.99)


###
### Exercice 7
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
th <- 0.05; lam <- 1; q <- 0.4

fgp <- function(t) 1 - th + th * exp(lam * (q * t/(1 - (1 - q) * t) - 1))
fx <- Re(fft(sapply(e, fgp)))/n
fx[0:3 + 1]


###
### Exercice 8
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
w <- c(0.8, 0.2); lam <- c(0.1, 0.25); q <- 1/3

fgp <- function(t) sum(w * exp(lam * (q * t/(1 - (1 - q) * t) - 1)))
fx <- Re(fft(sapply(e, fgp)))/n
fx[0:3 + 1]


###
### Exercice 9
###

rm(list = ls())

n <- 2^15; k <- 0:(n - 1); v <- 0.95; h <- 1000; kh <- k * h
r <- 2; be <- 1; q <- 1/(1 + be); al <- 3; lam <- 4000

Fb <- function(x, a, b) 1 - (b/(b + x))^a
fgpW <- function(t) (q/(1 - (1 - q) * t))^r

fc1 <- c(0, diff(Fb(kh, al, lam * v)))
fc2 <- c(0, diff(Fb(kh, al, lam * v^2)))

fy1 <- Re(fft(fgpW(fft(fc1)), TRUE))/n
fy2 <- Re(fft(fgpW(fft(fc2)), TRUE))/n

fz <- Re(fft(fft(fy1) * fft(fy2), TRUE))/n
fz[c(0, 10, 20) + 1]

VaRZ <- function(u) kh[min(which(cumsum(fz) >= u))]
VaRZ(0.99)


###
### Exercice 10
###

rm(list = ls())

a <- c(0.7, 0.2, 0.1); b <- c(0.3, 0.5, 0.2)
B <- rbind(c(b, 0, 0), c(0, b, 0), c(0, 0, b))

c <- a %*% B
as.vector(c)


###
### Exercice 11
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n); i <- 1:5
lam <- 0.06 - 0.01 * i; lamN <- sum(lam); p <- lam/lamN

fgp <- function(t) exp(lamN * (sum(p * t^i) - 1))
fs <- Re(fft(sapply(e, fgp)))/n
fs[c(0, 1, 2, 3) + 1]


#===============================================================================
#== Exercices informatiques ====================================================
#===============================================================================

###
### Exercice 1
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
fgp <- function(t) (0.5 + 0.1 * t + 0.35 * t^2 + 0.05 * t^3)^10

fs <- Re(fft(fgp(e)))/n
fs[c(0, 1, 10, 20) + 1]


###
### Exercice 2
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
fgpM <- function(t) (0.6 + 0.3 * exp(0.1 * (t - 1)) + 0.1 * exp(0.2 * (t - 1)))
fgpN <- function(t) fgpM(t)^100

fn <- Re(fft(fgpN(e)))/n
fn[c(0, 5, 10, 15) + 1]


###
### Exercice 3
###

rm(list = ls())

n <- 2^17; k <- 0:(n - 1)
h <- 0.1; kh <- k * h; u <- 0.99
al <- 1.5; lam <- 5

Fb <- function(x) 1 - (lam/(lam + x))^al

## a)
fxu <- c(diff(Fb(kh)), 0)
fsu <- Re(fft(fft(fxu)^2, TRUE))/n
kh[min(which(cumsum(fsu) >= u))]

## b)
fxl <- c(0, diff(Fb(kh)))
fsl <- Re(fft(fft(fxl)^2, TRUE))/n
kh[min(which(cumsum(fsl) >= u))]


###
### Exercice 4
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); h <- 1
r <- 1.5; q <- 1/3; mu <- log(10) - 0.18; sig <- 0.6

## a)
fb <- c(0, diff(plnorm(k, mu, sig)))
fb[c(0, 10) + 1]

## b)
fgp <- function(t) (q/(1 - (1 - q) * t))^r
fx <- Re(fft(fgp(fft(fb)), TRUE))/n
fx[c(0, 20) + 1]

## c)
Fx <- cumsum(fx)
Fx[50 + 1]

stoploss <- function(d) sum(pmax(k - d, 0) * fx)
stoploss(50)

## d)
50 + stoploss(50)/(1 - Fx[50 + 1])


###
### Exercice 5
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); h <- 1000; kh <- k * h
mu <- 8; sig <- 1; al <- 4; be <- 1/1000

fx1 <- c(0, diff(plnorm(kh, mu, sig)))
fx2 <- c(0, diff(pgamma(kh, al, be)))

## a)
fs <- Re(fft(fft(fx1) * fft(fx2), TRUE))/n
fs[c(10, 20, 30) + 1]

## b)
VaR <- function(u) kh[min(which(cumsum(fs) >= u))]
sapply(c(0.95, 0.995), VaR)


###
### Exercice 6
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
m <- c(10, 20); q <- c(0.2, 0.3)

fgp <- function(t) prod((1 - q + q * t)^m)
fs <- Re(fft(sapply(e, fgp)))/n


###
### Exercice 7
###

rm(list = ls())

n <- 2^18; k <- 0:(n - 1); h <- 1000; kh <- k * h
m <- 1000; r <- 0.01; q <- 1/4; al <- 2.5; lam <- 15000
Fb <- function(x) 1 - (lam/(lam + x))^al

## a)
EspS <- m * r * (1 - q)/q * lam/(al - 1)
EspS

## b)
fb <- c(0, diff(Fb(kh)))
fgp <- function(t) (q/(1 - (1 - q) * t))^(r * m)

fs <- Re(fft(fgp(fft(fb)), TRUE))/n
sum(fs * kh)

# Attention, dans le corrigé le 'n' utilisé n'est probablement pas suffisant
sum(pmax(kh - 2000000, 0) * fs)
cumsum(fs)[c(500000, 1000000, 2000000)/h + 1]


###
### Exercice 8
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
j <- 1:5; lam <- 0.6 - 0.1 * j; be <- 1/1000

## b)
lamN <- sum(lam); pj <- lam/lamN

## c)
Fc <- function(x) sum(pj * pgamma(x, j, be))
sapply(c(2000, 8000), Fc)

## e)
fk <- Re(fft(sapply(e, function(t) exp(lamN * (sum(pj * t^j) - 1)))))/n
fk[0:3 + 1]

## f)
Fs <- function(x) fk[1] + sum(fk[-1] * pgamma(x, k[-1], be))
sapply(c(0, 2000, 8000, 20000), Fs)


###
### Exercice 9
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); h <- 100; kh <- k * h
al <- 2.9 - 0:4 * 0.2; eta <- 1900 - 0:4 * 200; lam <- 0.6 - 1:5 * 0.1

fgp <- function(t) exp(lamN * (3.8 - 2.2) * (t - 1))
fb <- sapply(
    1:5, function(i) c(0, diff(ppareto(0:199 * h, al[i], eta[i])),
    ppareto(199 * h, al[i], eta[i], lower = FALSE), rep(0, n - 201))
)

## b)
lamN <- sum(lam); pj <- lam / lamN

## c)
fc <- sapply(k, function(k) sum(pj * fb[k + 1, ]))
fc[c(20, 80) + 1]

## d)
EspN <- lamN * (3.8 - 2.2)
EspC <- sum(fc * kh)
EspS <- EspN * EspC

## e)
fs <- Re(fft(fgp(fft(fc)), TRUE))/n
fs[c(0, 20, 80) + 1]

## f)
Fs <- cumsum(fs)
Fs[c(0, 20, 80, 100) + 1]

## g)
stoploss <- function(d) sum(pmax(kh - 100 * d, 0) * fs)
sapply(c(0, 20, 80, 200), stoploss)

## h)
d <- c(0, 20, 80, 100)
sapply(d, function(d) d * h + stoploss(d)/(1 - Fs[d + 1]))


###
### Exercice 10
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n); m <- 20
be <- c(0.1, 0.5); al <- 5/8; q <- be[1]/be[2]

fgp <- function(t) al * q * t/(1 - (1 - q) * t) + (1 - al) * t
fk <- Re(fft(sapply(e, fgp)))/n
fk[1:5 + 1]

fm <- Re(fft(fft(fk)^m, TRUE))/n
fm[c(50, 60, 70) + 1]

Fs <- function(x) fm[1] + sum(fm[-1] * pgamma(x, k[-1], be[2]))
sapply(c(100, 140, 200, 300), Fs)


#===============================================================================
#== Allocation du risque =======================================================
#===============================================================================

###
### Exemple avec des distributions binomiales négatives
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); valid <- 1:100
r <- c(0.5, 8); q <- c(0.2, 0.5)

f1 <- dnbinom(k, r[1], q[1])
f2 <- dnbinom(k, r[2], q[2])

fs <- Re(fft(fft(f1) * fft(f2), TRUE))/n

c1 <- Re(fft(fft(k * f1) * fft(f2), TRUE))/n
c2 <- Re(fft(fft(k * f2) * fft(f1), TRUE))/n

phi1 <- c1/fs
phi2 <- c2/fs

all.equal((phi1 + phi2)[valid], k[valid])
c(sum(c1), sum(c2))


###
### Exemple avec des distributions de Poisson
###

rm(list = ls())

## a) On a que Xi | S = k ~ Binom(k, lam_i/lam_S)

## b) On a phi_i(S) = lam_i/lam_S * S avec S ~ Poisson(lam_S)

## c) On a E[phi_i(S)] = lam_i/lam_S * E[S] = lam_i

## d) On a Var(phi_i(S)) = (lam_i/lam_S)^2 * Var(S) = lam_i^2/lam_S


#===============================================================================
#== Distributions de fréquence et FFT ==========================================
#===============================================================================

###
### Distribution de Delaporte
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
lam <- 1.5; r <- 0.5; q <- 1/4

fgp <- function(t) q^r * exp(lam * (t - 1))/(1 - (1 - q) * t)^r
fm <- Re(fft(sapply(e, fgp)))/n


###
### Distribution de Delaporte composée
###

rm(list = ls())

n <- 2^13; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
lam <- 1.5; r <- 0.5; q <- 1/4; sig <- 0.8; mu <- log(100) - sig^2/2

fgp <- function(t) q^r * exp(lam * (t - 1))/(1 - (1 - q) * t)^r
fb <- c(0, diff(plnorm(k, mu, sig)))

fx <- Re(fft(fgp(fft(fb)), TRUE))/n
Fx <- cumsum(fx)

VaR <- function(u) k[min(which(Fx >= u))]
sapply(c(0.5, 0.9, 0.99, 0.999, 0.9999), VaR)


###
### Distribution inverse gaussienne composée
###

rm(list = ls())

n <- 2^13; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
lam <- 3; be <- 0.5; sig <- 0.8; mu <- log(100) - sig^2/2

fgp <- function(t) exp((1 - sqrt(1 - 2 * be * lam * (t - 1)))/be)
fb <- c(0, diff(plnorm(k, mu, sig)))

fx <- Re(fft(fgp(fft(fb)), TRUE))/n
Fx <- cumsum(fx)

VaR <- function(u) k[min(which(Fx >= u))]
sapply(c(0.5, 0.9, 0.99, 0.999, 0.9999), VaR)


#===============================================================================
#== Allocation du risque - volet 1 =============================================
#===============================================================================

###
### Exercice 1
###

rm(list = ls())

al <- c(10, 4, 1); be <- c(1/300, 1/500, 1/1000)
m <- 1E5; set.seed(2019); u <- c(0.9, 0.99, 0.999)
U <- matrix(runif(3 * m), nrow = 3, ncol = m)

X1 <- qgamma(U[1, ], al[1], be[1])
X2 <- qgamma(U[2, ], al[2], be[2])
X3 <- qgamma(U[3, ], al[3], be[3])
S <- X1 + X2 + X3

## b)
c(mean(X1), mean(X2), mean(X3), mean(S))

## c)
VaR <- function(u) sort(S)[m * u]
sapply(u, VaR)

CVaR <- function(x, u) x[S == VaR(u)]
sapply(list(X1, X2, X3), function(x) sapply(u, function(u) CVaR(x, u)))

## d)
TVaR <- function(u) mean(S[S > VaR(u)])
sapply(u, TVaR)

CTVaR <- function(x, u) mean(x[S > VaR(u)])
sapply(list(X1, X2, X3), function(x) sapply(u, function(u) CTVaR(x, u)))


###
### Exercice 2
###

rm(list = ls())

alP <- 2.25; lam <- 125; alG <- 1/9; be <- 1/900
sig <- sqrt(log((300^2/100^2) + 1)); mu <- log(100) - sig^2/2
m <- 1E5; set.seed(2019); u <- c(0.9, 0.99, 0.999)
U <- matrix(runif(3 * m), nrow = 3, ncol = m)

X1 <- qpareto(U[1, ], alP, lam)
X2 <- qgamma(U[2, ], alG, be)
X3 <- qlnorm(U[3, ], mu, sig)
S <- X1 + X2 + X3

## b)
c(mean(X1), mean(X2), mean(X3), mean(S))

## c)
VaR <- function(u) sort(S)[m * u]
sapply(u, VaR)

CVaR <- function(x, u) x[S == VaR(u)]
sapply(list(X1, X2, X3), function(x) sapply(u, function(u) CVaR(x, u)))

## d)
TVaR <- function(u) mean(S[S > VaR(u)])
sapply(u, TVaR)

CTVaR <- function(x, u) mean(x[S > VaR(u)])
sapply(list(X1, X2, X3), function(x) sapply(u, function(u) CTVaR(x, u)))


#===============================================================================
#== Allocation du risque - volet 2 =============================================
#===============================================================================

###
### Exercice traditionnel
###

rm(list = ls())

be <- c(0.1, 0.2); ci <- sapply(1:2, function(i) prod(be[-i]/(be[-i] - be[i])))
Fs <- function(x) sum(ci * pexp(x, be))

## a)
VaRS <- function(u) optimize(function(x) abs(Fs(x) - u), c(0, 100))$min
VaRS(0.99)

TVaRS <- function(u) sum(ci * exp(-be * VaRS(u)) * (1/be + VaRS(u)))/(1 - u)
TVaRS(0.99)

## b) Long à calculer !
