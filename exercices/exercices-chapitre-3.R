################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions exercices - Chapitre 3 ##########################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

#===============================================================================
#== Exercices traditionnels ====================================================
#===============================================================================

###
### Exercice 1
###

rm(list = ls())

be <- 1/100; mu <- log(100) - 1/2; sig <- 1
VaR <- function(u) qexp(u, be) + qlnorm(u, mu, sig)
VaR(0.95)


###
### Exercice 2
###

rm(list = ls())

al <- 3; lam <- 2000; tau <- 0.5; be <- 1/2000
U <- pexp(3728, 1/1400)

X1 <- lam * ((1 - U)^(-1/al) - 1)
X2como <- ((-log(1 - U))^(1/tau))/be
X2anti <- ((-log(U))^(1/tau))/be

c(X1 + X2como, X1 + X2anti)


###
### Exercice 3
###

rm(list = ls())

mu <- c(3, 2); sig <- 1
muComo <- log(sum(exp(mu)))

## a)
FsComo <- function(x) plnorm(x, muComo, sig)
FsComo(40)

## b)
FsAnti <- function(x)
{
    c <- 1/sig * log((x - sqrt(x^2 - 4 * exp(sum(mu))))/(2 * exp(mu[1])))
    d <- 1/sig * log((x + sqrt(x^2 - 4 * exp(sum(mu))))/(2 * exp(mu[1])))
    pnorm(d) - pnorm(c)
}
FsAnti(40)


###
### Exercice 4
###

rm(list = ls())

mu <- 0.05; sig <- 0.015; u <- c(0.005, 0.995)
a <- c(-3.2, -4.3); b <- c(-0.12, -0.35); c <- c(1000, 2000)

## a)
EspX <- c * exp(b) * exp(mu * a + sig^2 * a^2/2)
EspS <- sum(EspX)

EspX2 <- c^2 * exp(2 * b) * exp(mu * 2 * a + sig^2 * (2 * a)^2/2)
EspX1X2 <- prod(c) * exp(sum(b)) * exp(mu * sum(a) + sig^2 * sum(a)^2/2)

VarX <- EspX2 - EspX^2
VarS <- sum(VarX) + 2 * (EspX1X2 - prod(EspX))

## b)
sapply(u, function(u) qnorm(u, mu, sig))

## c)
sapply(u, function(u) sum(c * exp(b + a * qnorm(1 - u, mu, sig))))


###
### Exercice 10
###

rm(list = ls())

q <- c(0.1, 0.2); be <- c(0.5, 1); u <- c(0.5, 0.85, 0.999)

Fxi <- function(x, i) 1 - q[i] + q[i] * pexp(x, be[i])
EspTrB <- function(d, i) exp(-be[i] * d)/be[i] + d * exp(-be[i] * d)

Fx12como <- function(x1, x2) min(Fxi(x1, 1), Fxi(x2, 2))
Fx12anti <- function(x1, x2) max(Fxi(x1, 1) + Fxi(x2, 2) - 1, 0)
cdfs <- list(Fx12como, Fx12anti)

## a)
VaRX <- function(u, i) ifelse(
    u < 1 - q[i], 0, qexp((u - (1 - q[i]))/q[i], be[i])
)

VaRS <- function(u) VaRX(u, 1) + VaRX(u, 2)
sapply(u, VaRS)

## b)
TVaRX <- function(u, i) 1/(1 - u) * (
    q[i] * EspTrB(VaRX(u, i), i) + VaRX(u, i) * (Fxi(VaRX(u, i), i) - u)
)

TVaRS <- function(u) TVaRX(u, 1) + TVaRX(u, 2)
sapply(u, TVaRS)

## c)
c(Fx12como(0, 0), Fx12anti(0, 0))

## d)
sapply(cdfs, function(Fx12) Fxi(0, 1) - Fx12(0, 0))

## e)
sapply(cdfs, function(Fx12) Fxi(0, 2) - Fx12(0, 0))

## f)
zapsmall(sapply(cdfs, function(Fx12) 1 - Fxi(0, 1) - Fxi(0, 2) + Fx12(0, 0)))


###
### Exercice 14
###

rm(list = ls())

lam <- 1; be <- 1
fgmB <- function(t) be/(be - t)

## a)
phi <- function(k, t) (lam * fgmB(t) - 1 - log(1 - k))/t

## b)
tk  <- function(k) 1 - (
    (lam - sqrt(lam^2 - lam * (1 + log(1 - k))))/(2 * (1 + log(1 - k)))
)

## c)
eVaR <- function(k) phi(k, tk(k))
eVaR(0.9)


###
### Exercice 16
###

rm(list = ls())

n <- 8; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
fgpX <- function(t1, t2) (0.18 + 0.05 * t1 + 0.05 * t2 + 0.72 * t1 * t2)^3

## a)
fx1 <- Re(fft(fgpX(e, 1)))/n
fx1[0:3 + 1]

## b)
fx2 <- Re(fft(fgpX(1, e)))/n
fx2[0:3 + 1]

## d)
s <- 0:7
fs <- Re(fft(fgpX(e, e)))/n
fs[0:6 + 1]

## e)   
fx1_k_et_x2_0 <- c(0.005832, 0.00486, 0.00135, 0.000125) 
sum(0:3 * fx1_k_et_x2_0)

## f)
sum(s * fs * exp(0.02 * s)) # E[Se^(0.02S)]
sum(fs * exp(0.02 * s)) # E[e^(0.02S)]


###
### Exercice 19
###

rm(list = ls())

n <- 8; k <- 0:(n - 1); e <- exp(2i * pi * k/n)

fgpM <- function(t1, t2) (
    0.6 + 0.1 * t2^5 + 0.05 * t1 * t2^4 + 0.05 * t1^2 * t2^3 +
    0.05 * t1^3 * t2^2 + 0.05 * t1^4 * t2 + 0.1 * t1^5
)

## a)
fm1m2 <- zapsmall(Re(fft(outer(e, e, fgpM)))/n^2)

## b)
fm1 <- zapsmall(Re(fft(sapply(e, function(t) fgpM(t, 1))))/n)
fm2 <- zapsmall(Re(fft(sapply(e, function(t) fgpM(1, t))))/n)

## c)
EspM1 <- sum(k * fm1)
EspM2 <- sum(k * fm2)

EspM1M2 <- sum(outer(k, k) * fm1m2)
CovM1M2 <- EspM1M2 - EspM1 * EspM2


###
### Exercice 20
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1)

fk1 <- dnbinom(k, 2, 1/2)
fk2 <- dnbinom(k, 1, 1/4)

fn <- Re(fft(fft(fk1) * fft(fk2), TRUE))/n
fn[0:2 + 1]


###
### Exercice 21
###

rm(list = ls())

i <- 1:5; be <- 1/(10 * i); beS <- 1/sum(1/be)

## b)
Fs <- function(x) pexp(x, beS)
Fs(200)

## c)
VaRS <- function(u) qexp(u, beS)
VaRS(0.95)


###
### Exercice 22
###

rm(list = ls())

i <- 1:5; al <- 3; lam <- 20 * i; lamS <- sum(lam)

## b)
Fs <- function(x) 1 - (lamS/(lamS + x))^al
Fs(200)

## c)
VaRS <- function(u) lamS * ((1 - u)^(-1/al) - 1)
VaRS(0.95)


###
### Exercice 23
###

rm(list = ls())

mu <- 4; sig <- 2; al <- 2.1; lam <- 440

VaRS <- function(u) qlnorm(u, mu, sig) + lam * ((1 - u)^(-1/al) - 1)
VaRS(0.95)

TVaRS <- function(u) (
    exp(mu + sig^2/2) * pnorm(qnorm(u) - sig, lower = FALSE) +
    lam * (al/(al - 1) * (1 - u)^(-1/al) - 1)
)
TVaRS(0.95)


###
### Exercice 26 À TERMINER
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
a0 <- 0.5; lam <- c(2, 1); i <- 1:2; C <- 100

fu <- sapply(1:2, function(i) c(0, sapply(1:5, function(j) choose(4, j - 1) *
    (0.6 * i - 0.4)^(j - 1) * (1.4 - 0.6 * i)^(5 - j)), rep(0, n - 6))
)

fx <- sapply(1:2, function(i) Re(
    fft(exp(lam[i] * (fft(fu[, i]) - 1)), TRUE))/n
)

fs <- Re(fft(exp(
    (lam[1] - a0) * (fft(fu[, 1]) - 1) + (lam[2] - a0) * (fft(fu[, 2]) - 1) +
    a0 * (fft(fu[, 1]) * fft(fu[, 2]) - 1)), TRUE)/n
)

## a)
EspX1 <- C * sum(k * fx[, 1])
EspX2 <- C * sum(k * fx[, 2])
EspS <- EspX1 + EspX2

## b)
CovX1X2 <- a0 * C^2 * sum(k * fu[, 1]) * sum(k * fu[, 2])

## c)
VarS <- sum(fs * (k * C - EspS)^2)

## d)
fs[0:5 + 1]

## e) à faire


###
### Exercice 29
###

rm(list = ls())

p <- rbind(c(0.42, 0.18), c(0.28, 0.12)); be <- 1/10; i <- 1:2

## a)
Fx12 <- function(x1, x2) sum(p * outer(pgamma(x1, i, be), pgamma(x2, i, be)))
Fx12(30, 20)

## b)
EspX1TrX2 <- function(x2) sum(p * outer(i/be, pgamma(x2, i, be, lower = FALSE)))
EspX1TrX2(20)

## c) réponse ne concorde pas avec le corrigé
SLgamma <- function(x, j) (
    j/be * pgamma(x, j + 1, be, low = FALSE) - x * pgamma(x, j, be, low = FALSE)
)

SLX1X2 <- function(x1, x2) sum(p * outer(SLgamma(x1, i), SLgamma(x2, i)))
SLX1X2(30, 20)

## d)
EspX1 <- sum(rowSums(p) * i/be)
EspX2 <- sum(colSums(p) * i/be)

EspX1X2 <- sum(p * outer(i/be, i/be))
CovX1X2 <- EspX1X2 - EspX1 * EspX2

VarX1 <- sum(rowSums(p) * (i * (i + 1))/be^2) - EspX1^2
VarX2 <- sum(colSums(p) * (i * (i + 1))/be^2) - EspX2^2

EspS <- EspX1 + EspX2
VarS <- VarX1 + VarX2 + 2 * CovX1X2

## e)
Fs <- function(x) (
    p[1, 1] * pgamma(x, 2, be) +
    sum(p[1, 2], p[2, 1]) * pgamma(x, 3, be) +
    p[2, 2] * pgamma(x, 4, be)
)
Fs(50)


###
### Exercice 32
###

rm(list = ls())

q <- 0.1; al <- 3; lam <- c(20, 40)

## a)

## b)

## d)
Fpareto <- function(x, a, b) 1 - (b/(b + x))^a
Fb <- function(x) (
    lam[2]/(lam[2] - lam[1]) * Fpareto(x, al, lam[2]) -
    lam[1]/(lam[2] - lam[1]) * Fpareto(x, al, lam[1])
)

Fx <- function(x) 1 - q + q * Fb(x)
1 - Fx(70)


#===============================================================================
#== Exercices informatiques ====================================================
#===============================================================================

###
### Exercice 1
###

rm(list = ls())

al <- 1.5; lam <- 50; be <- 1/100; mu <- log(100) - 1/2; sig <- 1

## b)
VaRS <- function(u) qpareto(u, al, lam) + qexp(u, be) + qlnorm(u, mu, sig)
VaRS(0.95)

## c)
Fs <- function(x) optimize(function(u) abs(VaRS(u) - x), c(0, 1))$minimum
Fs(1000)


###
### Exercice 3
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
r <- 2; lam <- c(1.5, 2); gam <- 0.8; m <- 10; q <- c(0.2, 0.3)

fgpB1 <- function(t) (1 - q[1] + q[1] * t)^m
fgpB2 <- function(t) (1 - q[2] + q[2] * t)^m

fgmT <- function(t1, t2) (
    ((1 - gam)/((1 - (1 - gam)/r * t1) * (1 - (1 - gam)/r * t2) - gam))^r
)
fgpX <- function(t1, t2) fgmT(
    lam[1] * (fgpB1(t1) - 1), lam[2] * (fgpB2(t2) - 1)
)

fgpS <- function(t) fgpX(t, t)
fs <- Re(fft(fgpS(e)))/n
cumsum(fs)[c(30, 40) + 1]


#===============================================================================
#== Exercices supplémentaires ==================================================
#===============================================================================

###
### 3. Chaîne Markov-Bernoulli et distribution Markov-binomiale
###

rm(list = ls())

q <- 1/3; al <- 1/5
p00 <- 1 - (1 - al) * q
p01 <- (1 - al) * q
p10 <- (1 - al) * (1 - q)
p11 <- al + (1 - al) * q

c(p00 + p01, p10 + p11)

# i)
EspM30 <- 30 * q
EspM30

# ii)
prob.pluit.5.consecutifs <- q * p11^4
prob.pluit.5.consecutifs

# iii)
fm2 <- c((1 - q) * p00, (1 - q) * p01 + q * p10, q * p11)
fm2


###
### 4. Distribution Bernoulli multivariée échangeable particulière
###

rm(list = ls())

## a)
p1100 <- p1010 <- p1001 <- p0110 <- p0101 <- p0011 <- 1/6

## b)
f123.100 <- f123.010 <- f123.001 <- 1/6
f123.110 <- f123.101 <- f123.011 <- 1/6

## c)
f12.00 <- 1/6
f12.10 <- f12.01 <- 1/3
f12.11 <- 1/6

## d)
f1.0 <- f1.1 <- q <- 1/2

## e)
EspN <- 2
VarN <- 0

## f)
Cov12 <- -q * (1 - q)/3         # -1/12
Cov12

## g)
rho12 <- Cov12/(q * (1 - q))    # -1/3
rho12

## i) À FAIRE!
