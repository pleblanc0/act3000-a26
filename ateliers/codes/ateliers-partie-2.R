################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### École d'actuariat - Université Laval ######################################
#### Ateliers - Partie 2 #######################################################
#### Auteur : Philippe Leblanc #################################################
################################################################################

###
### Utiliser l’algorithme FFT adroitement 1 (p.12)
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); e <- exp(2i * pi * k/n); lam <- 3

fgp <- function(t) exp(lam * (t - 1))
fm <- Re(fft(fgp(e)))/n

# alternative
fb <- c(0, 1, rep(0, n - 2))
dftb <- fft(fb)

dftm <- fgp(dftb)
fm.alt <- Re(fft(dftm, TRUE))/n

# comparaison
all.equal(fm, fm.alt)


###
### Utiliser l’algorithme FFT adroitement 2 (p.13)
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
lam <- 2; be <- 0.5

fgp <- function(t) exp((1 - sqrt(1 - 2 * be * lam * (t - 1)))/be)
fm <- Re(fft(fgp(e)))/n

# alternative
fb <- c(0, 1, rep(0, n - 2))
dftb <- fft(fb)

dftm <- fgp(dftb)
fm.alt <- Re(fft(dftm, TRUE))/n

# comparaison
all.equal(fm, fm.alt)


###
### Exercice de base 1 (p.5)
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
i <- 1:10; r <- 2; q <- 1 - 0.01 * i

# méthode 1: utiliser directement la fgp
fgp <- function(t) prod((q/(1 - (1 - q) * t))^r)
fs <- Re(fft(sapply(e, fgp)))/n

# méthode 2: utiliser les fmps
fx <- sapply(1:10, function(i) dnbinom(k, r, q[i]))
dfts <- apply(mvfft(fx), 1, prod)
fs.alt <- Re(fft(dfts, TRUE))/n

# comparaison des deux méthodes
all.equal(fs, fs.alt)


###
### Exercice de base 2 (p.6)
###

rm(list = ls())

n <- 2^5; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
q <- c(0.2, 0.3, 0.4, 0.1, 0.5); b <- c(5, 8, 7, 9, 2)

## 1)
fgp <- function(t) prod(1 - q + q * t^b)
fs <- Re(fft(sapply(e, fgp)))/n

# méthode alternative: utiliser les fmps
fx <- sapply(
    1:5, function(i) c(1 - q[i], rep(0, b[i] - 1), q[i], rep(0, n - b[i] - 1))
)
dfts <- apply(mvfft(fx), 1, prod)
fs.alt <- Re(fft(dfts, TRUE))/n

# comparaison des deux méthodes
all.equal(fs, fs.alt)

## 2) utiliser le Théorème 2 des notes de cours (Chapitre 1)
ogfi <- function(t, i) q[i] * b[i] * t^b[i] * prod(1 - q[-i] + q[-i] * t^b[-i])
EspAll <- sapply(1:5, function(i) zapsmall(Re(fft(sapply(e, ogfi, i)))/n))
all.equal(colSums(EspAll), q * b)

EspCond <- sapply(1:5, function(i) EspAll[, i]/fs)
all.equal(rowSums(EspCond), as.vector(k * I(zapsmall(fs) > 0)))

# méthode alternative: utiliser les fmps
EA <- sapply(1:5, function(i) Re(
    fft(fft(k * fx[, i]) * apply(mvfft(fx[, -i]), 1, prod), TRUE)/n
))
all.equal(EspAll, EA)

EC <- sapply(1:5, function(i) zapsmall(EA[, i])/fs)
all.equal(EspCond, EC)


###
### Exercice de base 3 (p.7)
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); e <- exp(2i * pi * k/n); m <- 10
p <- c(0.5, 0.1, 0.2, 0.08, 0.05, 0.04, 0.02, 0.01)
b <- c(0, 1, 2, 5, 10, 20, 50, 100)

## 1)
EspX <- sum(b * p)
EspS <- m * EspX

## 2)
VarX <- sum(b^2 * p) - EspX^2
VarS <- m * VarX

## 3)
fgp <- function(t) (sum(p * t^b))^m
fs <- Re(fft(sapply(e, fgp)))/n

## 5)
stoploss <- function(d) sum(pmax(k - d, 0) * fs)
sapply(0:10 * 10, stoploss)


###
### Exercice de base 4 (p.8)
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
lam <- 5; al <- c(0.8, 0.2); nu <- c(0.5, 0.1)

## 1)
EspM <- lam
EspB <- sum(al/nu)
EspX <- EspM * EspB

## 2) Attention de ne pas faire une somme pondérée des variances pour VarB!
VarM <- lam
VarB <- sum(al * (2 - nu)/nu^2) - sum(al/nu)^2
VarX <- EspM * VarB + VarM * EspB^2

## 3)
fgp <- function(t) exp(lam * (sum(al * nu * t/(1 - (1 - nu) * t)) - 1))
fx <- Re(fft(sapply(e, fgp)))/n

## 4)
VaR <- function(u) k[min(which(cumsum(fx) >= u))]
TVaR <- function(u) sum(pmax(k - VaR(u), 0) * fx)/(1 - u) + VaR(u)

u <- c(0.5, 0.9, 0.99, 0.999)
cbind(VaR = sapply(u, VaR), TVaR = sapply(u, TVaR))

## 5)
stoploss <- function(d) sum(pmax(k - d, 0) * fx)
sapply(2:10 * 10, stoploss)


###
### Exercice de base 5 (p.9)
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
r <- c(0.5, 2.5); q <- c(1/11, 1/3)
a <- c(0.8, 0.2); nu <- c(0.5, 0.1)
b <- c(0.7, 0.3); eta <- c(0.25, 0.0625)

fgpC1 <- function(t) sum(a * nu * t/(1 - (1 - nu) * t))
fc1 <- Re(fft(sapply(e, fgpC1)))/n

fgpC2 <- function(t) sum(b * eta * t/(1 - (1 - eta) * t))
fc2 <- Re(fft(sapply(e, fgpC2)))/n

b1 <- pmin(35, pmax(k - 15, 0))
fb1 <- c(
    sum(fc1[b1 == 0]), fc1[b1 != 0 & b1 != 35],
    sum(fc1[b1 == 35]), rep(0, n - 36)
)

b2 <- pmin(50, pmax(k - 20, 0))
fb2 <- c(
    sum(fc2[b2 == 0]), fc2[b2 != 0 & b2 != 50],
    sum(fc2[b2 == 50]), rep(0, n - 51)
)

## 3)
fgpMi <- function(t, i) (q[i]/(1 - (1 - q[i]) * t))^r[i]
fxi <- function(i) Re(fft(sapply(fft(cbind(fb1, fb2)[, i]), fgpMi, i), TRUE))/n
fs <- Re(fft(fft(fxi(1)) * fft(fxi(2)), TRUE))/n

## 4)
VaR <- function(u) k[min(which(cumsum(fs) >= u))]
TVaR <- function(u) sum(pmax(k - VaR(u), 0) * fs)/(1 - u) + VaR(u)

u <- c(0.5, 0.9, 0.99, 0.999)
cbind(VaR = sapply(u, VaR), TVaR = sapply(u, TVaR))

## 5)
stoploss <- function(d) sum(pmax(k - d, 0) * fs)
sapply(5:15 * 10, stoploss)


###
### Algorithme FFT et partage de risque (p.18)
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); valid <- 1:100
lam <- 8; r <- 0.5; q <- 1/5

f1 <- dpois(k, lam)
f2 <- dnbinom(k, r, q)

fs <- Re(fft(fft(f1) * fft(f2), TRUE))/n

c1 <- Re(fft(fft(k * f1) * fft(f2), TRUE))/n
c2 <- Re(fft(fft(k * f2) * fft(f1), TRUE))/n

phi1 <- c1/fs
phi2 <- c2/fs

all.equal((phi1 + phi2)[valid], k[valid])
c(sum(c1), sum(c2))


###
### Processus de branchement et épidémie (p.21)
###

rm(list = ls())

n <- 2^13; lam <- 1.3; r <- 2; q <- 0.7

## a)
EspK <- function(k) lam^k
sapply(c(1, 2, 3, 7), EspK)

## b)
VarK <- function(k) if (k == 1) lam else lam^k + lam^2 * VarK(k - 1)
sapply(c(1, 2, 3, 7), VarK)

## Cas de base
prob_fermeture <- function(periode)
{
    fc <- fft(dpois(0:(n - 1), lam))
    for (k in 2:periode) fc <- exp(lam * (fc - 1))

    fn <- Re(fft(fc, TRUE))/n
    1 - cumsum(fn)[10 + 1]
}
sapply(c(2, 3, 7), prob_fermeture)

## Variante 1
prob_variante_1 <- function(periode)
{
    fc <- fft(dpois(0:(n - 1), lam))
    if (periode > 2)
        for (k in 2:(periode - 1)) fc <- exp(lam * (fc - 1))

    fc <- (q/(1 - (1 - q) * fc))^r
    fn <- Re(fft(fc, TRUE))/n
    1 - cumsum(fn)[10 + 1]
}
sapply(c(2, 3, 7), prob_variante_1)

## Variante 2
prob_variante_2 <- function(periode)
{
    fc <- fft(dpois(0:(n - 1), 1 + 1/periode))
    for (k in (periode - 1):1) fc <- exp((1 + 1/k) * (fc - 1))

    fn <- Re(fft(fc, TRUE))/n
    1 - cumsum(fn)[10 + 1]
}
sapply(c(2, 3, 7), prob_variante_2)


###
### Processus de Poisson (p.27)
###

rm(list = ls())

## 1)
lam <- -log(0.4)/2

## 2)
dpois(1, lam * (3.6 - 2))

## 3)
pexp(2, lam, lower = FALSE)

## 4)
EspW <- 1/lam
EspW

## 5)
lam * 1.5

## 6)
t1 <- 2; t2 <- 5;N2 <- 3
c(N2 + lam * (t2 - t1), lam * (t2 - t1))


###
### Processus de Poisson composé avec sinistres gamma (p.29)
###

rm(list = ls())

al <- 2; be <- 1/100; t <- 4 - 2.5
k0 <- 1:1000; u <- c(0.1, 0.5, 0.9, 0.99, 0.999)

## 1)
lam <- -log(0.1053992245)/(15/12)

## 2)
EspSt <- lam * t * al/be
EspSt

VarSt <- lam * t * al * (al + 1)/be^2
VarSt

## 3)
fn0 <- dpois(0, lam * t)
fnk <- dpois(k0, lam * t)

Fst <- function(x) fn0 + sum(fnk * pgamma(x, k0 * al, be))
1 - sapply(0:10 * 100, Fst)

## 4)
VaRSt <- function(u) optimize(function(x) abs(Fst(x) - u), c(0, 5000))$min
sapply(u, VaRSt)

TVaRSt <- function(u) 1/(1 - u) * sum(
    fnk * k0 * al/be * pgamma(VaRSt(u), k0 * al + 1, be, lower = FALSE)
)
sapply(u, TVaRSt)

## 5)
t1 <- 2; t2 <- 5; S2 <- 200

Fs5c2 <- function(x) dpois(0, lam * (t2 - t1)) + sum(
    dpois(k0, lam * (t2 - t1)) * pgamma(x - S2, k0 * al, be)
)
Fs5c2(500)

## 6)
EspS5c2 <- S2 + lam * (t2 - t1) * al/be
EspS5c2

VarS5c2 <- lam * (t2 - t1) * al * (al + 1)/be^2
VarS5c2


###
### Processus de Poisson composé avec sinistres discrets
###

rm(list = ls())

n <- 2^20; k <- 0:(n - 1); al <- 2.5; theta <- 300
t <- 4 - 2.5; u <- c(0.1, 0.5, 0.9, 0.99, 0.999)

Fpareto <- function(x, a, b) 1 - (b/(b + x))^a
fb <- c(0, diff(Fpareto(k, al, theta))/(Fpareto(n - 1, al, theta)))

## 1)
lam <- -log(0.1053992245)/(15/12)

## 2)
EspSt <- lam * t * sum(k * fb)
EspSt

VarSt <- lam * t * sum(k^2 * fb)
VarSt

## 3)
fs <- Re(fft(exp(lam * t * (fft(fb) - 1)), TRUE))/n
Fs <- cumsum(fs)
1 - Fs[0:10 * 100 + 1]

## 4)
VaRSt <- function(u) k[min(which(Fs >= u))]
sapply(u, VaRSt)

TVaRSt <- function(u) sum(pmax(k - VaRSt(u), 0) * fs)/(1 - u) + VaRSt(u)
sapply(u, TVaRSt)

## 5)
t1 <- 2; t2 <- 5; S2 <- 200

fs5c2 <- Re(fft(exp(lam * (t2 - t1) * (fft(fb) - 1)), TRUE))/n
cumsum(fs5c2)[500 - S2 + 1]

## 6)
EspS5c2 <- S2 + lam * (t2 - t1) * sum(k * fb)
EspS5c2

VarS5c2 <- lam * (t2 - t1) * sum(k^2 * fb)
VarS5c2
