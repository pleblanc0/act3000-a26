################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### École d'actuariat - Université Laval ######################################
#### Ateliers - Partie 3 #######################################################
#### Auteur : Philippe Leblanc #################################################
################################################################################

###
### Bornes inférieure et supérieure de Fréchet (lognormales)
###

rm(list = ls())

# truc pour les deux calculs: obtenir la covariance en utilisant intégrant
# sur [0, 1] l'expression des Xi en fonction de leur fonction quantile
# avant de diviser par les écarts-types pour obtenir la corrélation

# cas comonotone
rho_como <- function(sig) (exp(prod(sig)) - 1)/sqrt(prod((exp(sig^2) - 1)))
rho_como(c(1, 2))

# cas antimonotone
rho_anti <- function(sig) (exp(-prod(sig)) - 1)/sqrt(prod((exp(sig^2) - 1)))
rho_anti(c(1, 2))


###
### Somme de 3 v.a. comonotones
###

rm(list = ls())

# Paramètres des distributions
sig <- 0.8; mu <- log(100) - sig^2/2
al.gam <- 2.5; be <- al.gam/100
al.par <- 3; lam <- 200

# Préparation de la simulation
m <- 1E6; set.seed(20260830)

U.como <- runif(m)
U.indep <- matrix(runif(m * 3), ncol = 3, byrow = TRUE)

X.indep <- matrix(numeric(m * 3), ncol = 3)
X.como <- matrix(numeric(m * 3), ncol = 3)

# Simulation des v.a. indépendantes
X.indep[, 1] <- qlnorm(U.indep[, 1], mu, sig)
X.indep[, 2] <- qgamma(U.indep[, 2], al.gam, be)
X.indep[, 3] <- lam * ((1 - U.indep[, 3])^(-1/al.par) - 1)

# Simulation des v.a. comonotones
X.como[, 1] <- qlnorm(U.como, mu, sig)
X.como[, 2] <- qgamma(U.como, al.gam, be)
X.como[, 3] <- lam * ((1 - U.como)^(-1/al.par) - 1)

# Agrégation des composantes simulées
S.indep <- rowSums(X.indep)
head(S.indep)

S.como <- rowSums(X.como)
head(S.como)

# Espérance et variance de S
mean(S.indep); var(S.indep)
mean(S.como); var(S.como)

# Valeur exacte de l'espérance de S
EspX.Exact <- c(exp(mu + sig^2/2), al.gam/be, lam/(al.par - 1))
EspS.Exact <- sum(EspX.Exact)


###
### Loi Poisson bivariée Teicher
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
lam <- c(2, 3); a0 <- 1

# comparaison
e.alt <- fft(c(0, 1, rep(0, n - 2)))
all.equal(e, e.alt)

fgpM1M2 <- function(t1, t2) exp(
    a0 * (t1 * t2 - 1) + (lam[1] - a0) * (t1 - 1) + (lam[2] - a0) * (t2 - 1)
)
fgpN <- function(t) fgpM1M2(t, t)

## a)
rho <- a0

## b)
EspM1cM2 <- lam[1] - a0 + c(2, 4, 5) * a0/lam[2]
EspM2cM1 <- lam[2] - a0 + c(2, 4, 5) * a0/lam[1]

## d)
fm1m2 <- Re(fft(outer(e, e, fgpM1M2), TRUE))/n^2

EspM1 <- sum(k * rowSums(fm1m2))
EspM2 <- sum(k * colSums(fm1m2))

## e)
Fm1m2 <- t(apply(apply(fm1m2, 2, cumsum), 1, cumsum))
Fm1m2[6, 5] - Fm1m2[6, 2] - Fm1m2[3, 5] + Fm1m2[3, 2]

## f)
sum(outer(k, k) * fm1m2)
sum(outer(pmin(k, 7), I(k <= 6)) * fm1m2)

## g)
fn <- Re(fft(sapply(e, fgpN), TRUE))/n
