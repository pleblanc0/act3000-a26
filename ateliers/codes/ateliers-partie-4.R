################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### École d'actuariat - Université Laval ######################################
#### Ateliers - Partie 4 #######################################################
#### Auteur : Philippe Leblanc #################################################
################################################################################

###
### DCM et réchauffement
###

rm(list = ls())

al <- c(2.1, 3.4); be <- 1/10
m <- 0:3

fm1m2 <- matrix(
    c(
        0.30, 0.10, 0.02, 0.03,
        0.20, 0.10, 0.01, 0.01,
        0.04, 0.02, 0.05, 0.03,
        0.01, 0.02, 0.01, 0.05
    ), nrow = 4, byrow = TRUE
)

# marginales de fréquence
fm1 <- rowSums(fm1m2)
fm2 <- colSums(fm1m2)

# répartition des lois composées marginales
Fx1 <- function(x) sum(fm1 * pgamma(x, al[1], be))
Fx2 <- function(x) sum(fm2 * pgamma(x, al[2], be))

# espérances
EspM <- c(sum(fm1 * m), sum(fm2 * m))
EspB <- al/be
EspX <- EspM * EspB
EspS <- sum(EspX)

# covariances
EspM12 <- sum(outer(m, m) * fm1m2)
CovM1M2 <- EspM12 - prod(EspM)
CovX12 <- CovM1M2 * prod(EspB)

# variances
VarM <- c(sum(fm1 * m^2), sum(fm2 * m^2)) - EspM^2
VarB <- al/be^2
VarX <- EspM * VarB + VarM * EspB^2
VarS <- sum(VarX) + 2 * CovX12

# répartition conjointe
Fx1x2 <- function(x1, x2) sum(
    fm1m2 * outer(pgamma(x1, al[1] * m, be), pgamma(x2, al[2] * m, be))
)

# répartition de la somme
Fs <- function(s) sum(
    fm1m2 * pgamma(s, al[1] * m + al[2] * m, be)
)

# à continuer...


###
### DCM et Poisson bivariée Teicher
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1)
lam <- c(3, 2); y0 <- 0.5
sig <- c(0.8, 0.6); mu <- log(c(20, 30)) - sig^2/2

fgp <- function(t1, t2) exp(
    (lam[1] - y0) * (t1 - 1) + (lam[2] - y0) * (t2 - 1) + y0 * (t1 * t2 - 1)
)

# espérances
EspM <- lam
EspB <- exp(mu + sig^2/2)

EspX <- EspM * EspB
EspS <- sum(EspX)

# covariances
CovM12 <- y0
CovX12 <- CovM12 * prod(EspB)

# variances
VarM <- lam
VarB <- (exp(sig^2) - 1) * exp(2 * mu + sig^2)

VarX <- EspM * VarB + VarM * EspB^2
VarS <- sum(VarX) + 2 * CovX12

# simulation : voir le code d'Omar

# calculs avec FFT + discrétisation
fb1 <- c(0, diff(plnorm(k, mu[1], sig[1])))
fb2 <- c(0, diff(plnorm(k, mu[2], sig[2])))

# transformée de Fourier de S 
dfts <- fgp(fft(fb1), fft(fb2))

# fonction de masse de probabilité de S
fs <- Re(fft(dfts, TRUE)/n)

# fonction de répartition de S
Fs <- cumsum(fs)

# fonction stoploss
SL <- function(d) sum(pmax(k - d, 0) * fs)

# Value-at-Risk
VaR <- function(u) k[min(which(Fs >= u))]

# Tail Value-at-Risk
TVaR <- function(u) SL(VaR(u))/(1 - u) + VaR(u)

# calculs pour Zd = min(S; 120)
z <- pmin(k, 120)

EspZ <- sum(z * fs)
VarZ <- sum((z - EspZ)^2 * fs)


###
### DCM et Poisson trivariée Teicher
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1)
lam <- c(0.6, 0.9, 0.5)
y12 <- 0.1; y123 <- 0.2
p <- c(1/5, 1/8, 1/4)

# fgp conjointe de (M1, M2, M3)
fgp <- function(t1, t2, t3) exp(
    (lam[1] - y12 - y123) * (t1 - 1) + (lam[2] - y12 - y123) * (t2 - 1) +
    (lam[3] - y123) * (t3 - 1) + y12 * (t1 * t2 - 1) + y123 * (t1 * t2 * t3 - 1)
)

# espérances
EspM <- lam
EspB <- 1/p
EspX <- EspM * EspB
EspS <- sum(EspX)

# covariances
CovM12 <- y12 + y123
CovM13 <- y123
CovM23 <- y123

CovX12 <- CovM12 * prod(EspB[c(1, 2)])
CovX13 <- CovM13 * prod(EspB[c(1, 3)])
CovX23 <- CovM23 * prod(EspB[c(2, 3)])

# variances
VarM <- lam
VarB <- (1 - p)/p^2
VarX <- EspM * VarB + VarM * EspB^2
VarS <- sum(VarX) + 2 * (CovX12 + CovX13 + CovX23)

# fonctions de masse de probabilité de B
fb <- sapply(1:3, function(j) c(0, dgeom(k, p[j]))[-n])

# transformées de Fourier de B
dftb <- mvfft(fb)

# transformée de Fourier de S
dfts <- fgp(dftb[, 1], dftb[, 2], dftb[, 3])

# fonction de masse de probabilité de S
fs <- Re(fft(dfts, TRUE)/n)

# fonction de répartition de S
Fs <- cumsum(fs)

# Value-at-Risk
VaR <- function(u) k[min(which(Fs >= u))]
VaR(0.95)

# Tail Value-at-Risk
TVaR <- function(u) sum(pmax(k - VaR(u), 0) * fs)/(1 - u) + VaR(u)
TVaR(0.95)

# contributions à l'écart-type de S
# Cov(Xi, S)/Var(S) = (Var(Xi) + sum_{j != i} Cov(Xi, Xj))/Var(S)
ContrXiSd <- 1/sqrt(VarS) * c(
    VarX[1] + CovX12 + CovX13,
    VarX[2] + CovX12 + CovX23,
    VarX[3] + CovX13 + CovX23
)

# partage selon la règle de régression
# EspXi + Cov(Xi, S)/Var(S) * (S - EspS)
# EspXi + ContrXiSd * (S - EspS)/sqrt(VarS)
ContriXiReg <- sapply(
    1:3, function(j) EspX[j] + ContrXiSd[j] * (k - EspS)/sqrt(VarS)
)
zapsmall(rowSums(ContriXiReg))[1:10]


###
### Dépendance négative vs dépendance positive
###

rm(list = ls())

q <- 0.4; eta <- 0.1
gam <- 0.05; b <- 5 + (1:5)

# les 32 configurations possibles
ij <- as.matrix(expand.grid(rep(list(0:1), 5)))
nb <- rowSums(ij)

# probabilités conjointes
pi <- c(
    0.03028222, 0.04729603, 0.04679070,
    0.02345753, 0.00595927, 0.00095896
)

pj <- c(
    0.13387963, 0.04828479, 0.02807253,
    0.02104835, 0.02035259, 0.03172465
)

fi <- pi[nb + 1]
fj <- pj[nb + 1]

# vérification des masses totales
c(sum(fi), sum(fj))

# probabilités marginales
qi <- colSums(ij * fi)
qj <- colSums(ij * fj)

# matrices de covariance
CovI <- crossprod(ij, ij * fi) - outer(qi, qi)
CovJ <- crossprod(ij, ij * fj) - outer(qj, qj)

# matrices de corrélation
CorI <- cov2cor(CovI)
CorJ <- cov2cor(CovJ)

round(CorI, 1)
round(CorJ, 1)

# fmps de la somme
fm <- sapply(0:5, function(k) sum(fi[nb == k]))
fn <- sapply(0:5, function(k) sum(fj[nb == k]))

# cdfs de la somme
Fm <- cumsum(fm)
Fn <- cumsum(fn)

# stoploss de la somme
SLM <- function(d) sum(pmax(0:5 - d, 0) * fm)
SLN <- function(d) sum(pmax(0:5 - d, 0) * fn)

# fmp de la somme, où S = X1 + ... + X5, Xv = bv Iv
#                     T = Y1 + ... + Y5, Yv = bv Jv
montants <- as.vector(ij %*% b)
k <- 0:sum(b)

fs <- sapply(k, function(s) sum(fi[montants == s]))
all.equal(sum(fs * k), sum(b * qi))

ft <- sapply(k, function(t) sum(fj[montants == t]))
all.equal(sum(ft * k), sum(b * qj))

# cdfs de la somme
Fs <- cumsum(fs)
Fn <- cumsum(ft)

# stoploss de la somme
SLs <- function(d) sum(pmax(k - d, 0) * fs)
SLt <- function(d) sum(pmax(k - d, 0) * ft)

# Value-at-Risk de la somme
VaRs <- function(u) k[min(which(Fs >= u))]
VaRt <- function(u) k[min(which(Fn >= u))]
c(VaRs(0.9), VaRt(0.9))

# Tail Value-at-Risk de la somme
TVaRs <- function(u) SLs(VaRs(u))/(1 - u) + VaRs(u)
TVaRt <- function(u) SLt(VaRt(u))/(1 - u) + VaRt(u)
c(TVaRs(0.9), TVaRt(0.9))

##
## contributions
##

# espérances
EspX <- EspY <- b * q

# matrices de corrélation
CorI <- matrix(-eta, 5, 5)
CorJ <- matrix( eta, 5, 5)
diag(CorI) <- diag(CorJ) <- 1

# matrices de covariance
CovI <- q * (1 - q) * CorI
CovJ <- q * (1 - q) * CorJ

CovX <- CovI * outer(b, b) 
CovY <- CovJ * outer(b, b)

# alternative : calculer directement les covariances
#           CovX <- q * (1 - q) * outer(b, b) * CorI
#           CovY <- q * (1 - q) * outer(b, b) * CorJ

# contributions 
Allocation <- function(n, Esp, Cov)
{
    V <- Cov[1:n, 1:n, drop = FALSE]

    EspS <- sum(Esp[1:n])
    VarS <- sum(V)

    phi <- EspS + gam/2 * VarS
    contributions <- Esp[1:n] + gam/2 * rowSums(V)
    
    list(
        phi = phi, contributions = contributions
    )
}

# dépendance négative
Allocation(4, EspX, CovX)
Allocation(5, EspX, CovX)

# dépendance positive
Allocation(4, EspY, CovY)
Allocation(5, EspY, CovY)

# phi(Xv, S5) - phi(Xv, S4) = gam/2 Cov(Xv, X5
#                           = - 10 * gam/2 * nq(1 - q) * (5 + v)
#                           = -0.006 (5 + v) < 0
Allocation(5, EspX, CovX)$contributions[1:4] -
    Allocation(4, EspX, CovX)$contributions
-0.006 * (5 + 1:4)

# phi(Yv, S5) - phi(Yv, S4) = gam/2 Cov(Yv, Y5)
#                           = 10 * gam/2 * nq(1 - q) * (5 + v)
#                           = 0.006 (5 + v) > 0
Allocation(5, EspY, CovY)$contributions[1:4] -
    Allocation(4, EspY, CovY)$contributions
0.006 * (5 + 1:4)

# marginales pour I
MarginaleI <- function(d)
{
    configurations <- expand.grid(rep(list(0:1), d))
    names(configurations) <- paste0('i', 1:d)

    k <- rowSums(configurations)
    j <- 0:(5 - d)

    probabilites <- sapply(
        k, function(r) sum(choose(5 - d, j) * pi[r + j + 1])
    )

    data.frame(configurations, probabilite = probabilites)
}

# les 4 probabilités de (Iv, Iw)
MarginaleI(2)

# les 8 probabilités de (Iu, Iv, Iw)
MarginaleI(3)

# les 16 probabilités de (It, Iu, Iv, Iw)
MarginaleI(4)

# marginales pour J
MarginaleJ <- function(d)
{
    configurations <- expand.grid(rep(list(0:1), d))
    names(configurations) <- paste0('j', 1:d)

    k <- rowSums(configurations)
    j <- 0:(5 - d)

    probabilites <- sapply(
        k, function(r) sum(choose(5 - d, j) * pj[r + j + 1])
    )

    data.frame(configurations, probabilite = probabilites)
}

# les 4 probabilités de (Jv, Jw)
MarginaleJ(2)

# les 8 probabilités de (Ju, Jv, Jw)
MarginaleJ(3)

# les 16 probabilités de (Jt, Ju, Jv, Jw)
MarginaleJ(4)

# dernières vérifications
all.equal(round(sum((k - sum(fs * k))^2 * fs), 3), sum(CovX))
all.equal(round(sum((k - sum(ft * k))^2 * ft), 3), sum(CovY))
