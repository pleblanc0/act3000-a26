################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### École d'actuariat - Université Laval ######################################
#### Ateliers - Partie 5 #######################################################
#### Auteur : Philippe Leblanc #################################################
################################################################################

###
### Fonctions de répartition et quantile conditionnelles
###

rm(list = ls())

# paramètres de la loi de Pareto
al <- 2.5; lam <- 15

# paramètres de la copule FGM
theta <- c(-0.9, 0.9)

# fonction de répartition marginale
Fx <- function(x) 1 - (lam/(lam + x))^al

# fonction quantile marginale
Qx <- function(p) lam * ((1 - p)^(-1/al) - 1)

# copule FGM bivariée
Cfgm <- function(u, v, theta) u * v * (1 + theta * (1 - u) * (1 - v))

# fonction de répartition conjointe
Fx12 <- function(x1, x2, theta) Cfgm(Fx(x1), Fx(x2), theta)

# copule FGM conditionnelle: F_{U2|U1=u}(v) = C_{2|1}(v|u)
Cfgm2c1 <- function(u, v, theta) (
    v * (1 + theta * (1 - 2 * u)) - theta * (1 - 2 * u) * v^2
)

# fonction de répartition conditionnelle: F_{X2|X1=x1}(x2)
Fx2c1 <- function(x2, x1, theta) Cfgm2c1(Fx(x1), Fx(x2), theta)

# copule FGM conditionnelle inverse: F_{U2|U1=u}^{-1}(v) = C_{2|1}^{-1}(v|u)
Cfgm2c1.inv <- function(u, v, theta)
{
    a <- theta * (1 - 2 * u)
    delta <- (1 + a)^2 * (1 - v) + (1 - a)^2 * v

    ifelse(v == 0, 0, 2 * v/(1 + a + sqrt(delta)))
}

# Quantile conditionnel de U2 sachant U1 = u pour la copule FGM
Cfgm2c1.inv.details <- function(u, v, theta)
{
    # 1. on cherche le quantile x, c'est-à-dire la valeur qui vérifie
    #    v = x + theta * (1 - 2*u) * x * (1 - x). On pose a = theta * (1 - 2*u)
    a <- theta * (1 - 2 * u)

    # 2. on transforme l'équation en une équation quadratique.
    #    v = x + a * x * (1 - x) = x + a * x - a * x^2 = (1 + a) * x - a * x^2
    #    donc a * x^2 - (1 + a) * x + v = 0.
    #    on a que  A = a, B = -(1 + a) et C = v et delta = B^2 - 4 * A * C.
    delta <- (1 + a)^2 - 4 * a * v

    # 3. pour a =! 0, la formule quadratique donne deux racines.
    #    la racine dans [0, 1] est x = (1 + a - sqrt(delta)) / (2 * a).
    #    pour voir pourquoi on prend le signe moins, on remarque que
    #    delta = (1 + a - 2*a*x)^2. Comme |a| <= 1 et 0 <= x <= 1,
    #    1 + a - 2*a*x = 1 + a*(1 - 2*x) >= 0.
    #    donc sqrt(delta) = 1 + a - 2 * a * x.

    # 4. on simplifie la racine en multipliant le numérateur et
    #    le dénominateur par le conjugué : 1 + a + sqrt(delta).
    #    au numérateur, on obtient une différence de carrés :
    #    (1 + a - sqrt(delta)) * (1 + a + sqrt(delta))
    #    ce qui devient (1 + a)^2 - delta = 4 * a * v.
    #    en simplifiant par 2 * a, il reste x = 2 * v/(1 + a + sqrt(delta)).
    denom <- 1 + a + sqrt(delta)

    # 5. Si a = 0, la formule donne déjà x = 2*v/2 = v.
    #    Si v = 0, le quantile vaut 0. On le retourne directement
    #    pour éviter la division 0/0 lorsque a = -1 et v = 0.
    ifelse(v == 0, 0, 2 * v / denom)
}

# expression équivalente de la fonction précédente, mais plus longue
Cfgm2c1.inv.long <- function(u, v, theta)
{
    a <- theta * (1 - 2 * u)

    if (a < 0) return(sqrt(((1 + a)^2/(4 * a) - v)/a) + (1 + a)/(2 * a))

    if (a == 0) return(v)

    -sqrt(((1 + a)^2/(4 * a) - v)/a) + (1 + a)/(2 * a)
}

all.equal(
    Cfgm2c1.inv(Fx(10), Fx(20), theta[1]),
    Cfgm2c1.inv.long(Fx(10), Fx(20), theta[1])
)

# fonction de répartition conditionnelle inverse: F_{X2|X1=x1}^{-1}(x2)
Fx2c1.inv <- function(x2, x1, theta) Qx(Cfgm2c1.inv(Fx(x1), Fx(x2), theta))


###
### Simulation avec copule conditionnelle
###

rm(list = ls()) 

# paramètres de la copule FGM
theta <- c(-0.9, 0, 0.9)

# paramètres pour la simulation
n <- 5000; j <- c(1, 2, n)

# inverse conditionnelle de la copule FGM
Cfgm2c1.inv <- function(u, v, theta)
{
    a <- theta * (1 - 2 * u)
    delta <- (1 + a)^2 * (1 - v) + (1 - a)^2 * v

    ifelse(v == 0, 0, 2 * v/(1 + a + sqrt(delta)))
}

# initialisation pour reproductibilité
set.seed(20261118)

# on produit un seul échantillon de U1
U1 <- numeric(n)

# on produit trois échantillons de U2, un pour chaque valeur de theta
U2 <- matrix(0, nrow = n, ncol = length(theta))

# boucle de simulation pour produire les échantillons de U1 et U2
for (i in seq_len(n))
{   
    # on produit un échantillon de U1 et on le conserve
    U1[i] <- runif(1)

    # on produit un échantillon de V2 sans le conserver
    V2 <- runif(1)

    # on itère sur les trois valeurs de theta pour produire les échantillons
    for (k in seq_along(theta))
        U2[i, k] <- Cfgm2c1.inv(U1[i], V2, theta[k])
}

# vérification de la simulation pour les trois valeurs de theta
data.frame(
    j = j, U1 = U1[j], U2neg = U2[j, 1], U2nul = U2[j, 2], U2pos = U2[j, 3]
)

# graphiques pour visualiser l'impact de la dépendance sur les échantillons
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
plot(U1, U2[, 1], main = expression(theta == -0.9), xlab = 'U1', ylab = 'U2')
plot(U1, U2[, 2], main = expression(theta == 0), xlab = 'U1', ylab = 'U2')
plot(U1, U2[, 3], main = expression(theta == 0.9), xlab = 'U1', ylab = 'U2')


###
### Simulation avec copule conditionnelle et marginales Pareto + Exponentielle
###

# on continue avec les mêmes échantillons de U1 et U2

# paramètres des lois Pareto et Exponentielle
al <- 2.5; lam <- 15; be <- 0.1

# simulation des variables aléatoires X1 et X2 à partir de U1 et U2
X1 <- lam * ((1 - U1)^(-1/al) - 1)
X2 <- -log(1 - U2)/be
# X2 <- lam * ((1 - U2)^(-1/al) - 1)

# vérification de la simulation pour les trois valeurs de theta
data.frame(
    j = j,  X1 = X1[j], X2neg = X2[j, 1], X2nul = X2[j, 2], X2pos = X2[j, 3]
)

# graphiques pour visualiser l'impact de la dépendance sur les échantillons
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
plot(X1, X2[, 1], main = expression(theta == -0.9), xlab = 'X1', ylab = 'X2')
plot(X1, X2[, 2], main = expression(theta == 0), xlab = 'X1', ylab = 'X2')
plot(X1, X2[, 3], main = expression(theta == 0.9), xlab = 'X1', ylab = 'X2')


###
### Simulation des copules de base
###

rm(list = ls())

# nombre d'échantillons
n <- 5000

# initialisation pour reproductibilité
set.seed(42)

# échantillons de U pour les trois copules et V pour celle d'indépendance
U <- runif(n); V <- runif(n)

# graphiques pour visualiser les copules de base
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
plot(U, U, main = 'Copule borne sup. de Fréchet', xlab = 'U', ylab = 'V')
plot(U, 1 - U, main = 'Copule borne inf. de Fréchet', xlab = 'U', ylab = 'V')
plot(U, V, main = 'Copule indépendance', xlab = 'U', ylab = 'V')

###
### Copule FGM avec marginales Poisson
###

lam <- c(2, 3); n <- 2^10; k <- 0:(n - 1)
Fm1 <- ppois(k, lam[1])
Fm2 <- ppois(k, lam[2])

theta <- -0.9
C.fgm <- function(u, v) u * v + theta * u * v * (1 - u) * (1 - v)

# F(10, 12)
C.fgm(Fm1[1 + 1], Fm2[2 + 1])
C.fgm(Fm1[0 + 1], Fm2[0 + 1])

###
### Question 3 - Final informatique 2015
###

rm(list = ls())

m <- 100; k <- 0:m; s <- 0:(2 * m); al <- 10

C.clayton <- function(u, v) (u^(-al) + v^(-al) - 1)^(-1/al)

# (a) Paramètres de log(X1) et de log(X2)
mu <- log(c(30, 20)) + c(0.06, 0.05)
sig <- c(0.2, 0.25)

# (b) Fonctions de répartition (limite du support à x = 100)
Fx1 <- function(x) ifelse(x >= m, 1, plnorm(x, mu[1], sig[1]))
Fx2 <- function(x) ifelse(x >= m, 1, plnorm(x, mu[2], sig[2]))

Fx1x2 <- function(x1, x2) C.clayton(Fx1(x1), Fx2(x2))

fx1x2 <- function(x1, x2) (
    Fx1x2(x1, x2) - Fx1x2(x1 - 1, x2) -
    Fx1x2(x1, x2 - 1) + Fx1x2(x1 - 1, x2 - 1)
)
fx1x2.valeurs <- outer(k, k, fx1x2)

fx1x2(30, 20)
fx1x2(35, 25)
fx1x2(25, 15)

# (c) Loi de la somme
fs <- sapply(s, function(j) sum(fx1x2.valeurs[outer(k, k, '+') == j]))
Fs <- cumsum(fs)

# (d) Fonction de répartition de la somme
Fs[c(40, 60, 70) + 1]
Fs[50 + 1]

# (e) Espérances continue et discrétisée
EspS <- sum(exp(mu + sig^2 / 2))
EspS.disc <- sum(s * fs)

EspS
EspS.disc
EspS.disc - EspS

# (f) VaR de la somme discrétisée
VaR <- function(p) s[min(which(Fs >= p))]
sapply(c(0.01, 0.05, 0.95, 0.99), function(p) VaR(p))

# (g) TVaR de la somme discrétisée
TVaR <- function(p) VaR(p) + sum(pmax(s - VaR(p), 0) * fs) / (1 - p)
TVaR(0.95)
