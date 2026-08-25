################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Guide de survie - Partage de risques ######################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

#===============================================================================
#== Somme de lois binomiales négatives =========================================
#===============================================================================

rm(list = ls())

# Paramètres et préparation de l'exemple
n <- 2^16; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
r <- c(0.5, 2.5, 5); q <- c(1/11, 1/3, 1/2)

# Fonctions de masse de probabilités des Xi
fx <- sapply(1:3, function(i) dnbinom(k, r[i], q[i]))

# Fonction de masse de probabilité de S
fs <- Re(fft(apply(mvfft(fx), 1, prod), TRUE))/n

# Espérances et variances des composantes et de la somme
EspX <- r * (1 - q)/q
VarX <- r * (1 - q)/q^2
EspS <- sum(EspX)
VarS <- sum(VarX)

# Allocations espérées (expected allocations)
EA <- sapply(1:3, function(i) Re(
    fft(fft(k * fx[, i]) * apply(mvfft(fx[, -i]), 1, prod), TRUE)/n
))

# Allocations conditionelles (conditional allocations)
EC <- sapply(1:3, function(i) zapsmall(EA[, i])/fs)

# Allocations proportionelles
Yprop <- k/3

# Allocations par régression
Yreg <- sapply(1:3, function(i) EspX[i] + VarX[i]/VarS * (k - EspS))

#===============================================================================
#== Somme de lois Bernoulli composées ==========================================
#===============================================================================

rm(list = ls())

# Paramètres et préparation de l'exemple
m <- 5; n <- 2^15; h <-0.1; k <- 0:(n - 1) * h; u <- 0.99
q <- c(0.3, 0.2, 0.4, 0.15, 0.25); al <- 1:5/h; be <- 1/10
s <- c(1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 150, 200, 250)

# Espérances et variances des composantes et de la somme
EspX <- q * al/be
VarX <- q * al/be^2 + q * (1 - q) * (al/be)^2
EspS <- sum(EspX)
VarS <- sum(VarX)

# Fonction de masse de probabilité de N
fj <- matrix(0, n, m)
fj[1, ] <- 1 - q
for (i in 1:5) fj[al[i] + 1, i] <- q[i]

fn <- Re(fft(apply(mvfft(fj), 1, prod), TRUE))/n

# Fonction de densité et de répartition de S
fs <- function(x) sum(fn * dgamma(x, k, be))
Fs <- function(x) fn[1] + sum(fn[-1] * pgamma(x, k[-1], be))
Fs(50)
sapply(s, Fs)

# Allocations proportionelles
Yprop <- sapply(1:5, function(i) EspX[i]/EspS * s)

# Allocations par régression
Yreg <- sapply(1:5, function(i) EspX[i] + VarX[i]/VarS * (s - EspS))

# Allocations conditionelles
wi <- sapply(1:5, function(i) Re(
    fft(fft(k * fj[, i]) * apply(mvfft(fj[, -i]), 1, prod), TRUE)/n
))
Ycond <- sapply(1:5, function(i)
    sapply(s, function(x) sum(wi[, i]/be * dgamma(x, k + 1, be))/fs(x))
)

# Value-at-Risk
VaRS <- function(u) optimize(function(x) abs(Fs(x) - u), c(0, 400))$min
VaRS(0.99)

# Tail-Value-at-Risk
TVaRS <- function(u) sum(
    fn * k/be * pgamma(VaRS(u), k + 1, be, lower = FALSE)/(1 - u)
)
TVaRS(0.99)

# Contributions à la VaR
CVaR <- sapply(1:m, function(i)
    sum(wi[, i]/be * dgamma(VaRS(u), k + 1, be))/fs(VaRS(u))
)
c(sum(CVaR), VaRS(u))

# Contributions à la TVaR
CTVaR <- sapply(1:m, function(i)
    sum(wi[, i]/be * pgamma(VaRS(u), k + 1, be, lower = FALSE))/(1 - u)
)
c(sum(CTVaR), TVaRS(u))
