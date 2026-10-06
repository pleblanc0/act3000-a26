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

# Simulation des v.a. antimonotones
U.anti <- U.como
X.anti <- matrix(numeric(m * 2), ncol = 2)

X.anti[, 1] <- qlnorm(U.anti, mu, sig)
X.anti[, 2] <- qgamma(1 - U.anti, al.gam, be)

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


###
### Paires de v.a. Bernoulli
###

rm(list = ls())

q <- c(0.10, 0.25)
b <- c(1000, 400)

# Paramètres de corrélation
al.Min <- -min(
    sqrt((q[1] * q[2])/((1 - q)[1] * (1 - q)[2])),
    sqrt(((1 - q)[1] * (1 - q)[2])/(q[1] * q[2]))
)

# Cas comonotone
al.Max <- min(
  sqrt(((1 - q)[1] * q[2])/(q[1] * (1 - q)[2])),
  sqrt((q[1] * (1 - q)[2])/((1 - q)[1] * q[2]))
)

al <- c(
    al.Min, (7 * al.Min + 2 * al.Max)/10, al.Max
)

# Espérances et variances
EspX <- b * q
VarX <- b^2 * q * (1 - q)

# Proportional risk sharing
Yprop <- EspX/sum(EspX)
Yprop

# Variance de la somme
VarS <- sum(VarX) + 2 * prod(b) * al * sqrt(prod(q * (1 - q)))

# Espérance des contributions
EspY <- Yprop * sum(EspX)

# Variances de Y1 et Y2
VarY1 <- Yprop[1]^2 * VarS
VarY2 <- Yprop[2]^2 * VarS

# Bénéfices de mutualisation
BD1 <- sqrt(VarX[1]) - sqrt(VarY1)
BD2 <- sqrt(VarX[2]) - sqrt(VarY2)


###
### Bernoulli indépendantes et comonotones
###

rm(list = ls())

q <- readRDS(file.choose())$p

d <- length(q); n <- nextn(d, 2)
k <- 0:(n - 1); e <- exp(-2i * pi * k/n)

# Cas indépendance
fgpIndep <- function(t) prod(1 - q + q * t)
fn.indep <- Re(fft(sapply(e, fgpIndep), TRUE))/n

# Cas comonotone
q.como <- sort(q)
fn.Comon <- c(
    1 - q.como[d], rev(diff(q.como)), q.como[1]
)

# alternativement
fn.Comon.alt <- -diff(c(1, sort(q, TRUE), 0))

# comparaison
all.equal(fn.Comon, fn.Comon.alt)


###
### Distributions Bernoulli multivariées, partie 1
###

probabilities <- readRDS(file.choose())

# Préparation
i <- probabilities$i
f <- as.matrix(probabilities[, -1])
colSums(f)

indices <- t(sapply(strsplit(gsub('[()]', '', i), ','), as.numeric))

# Probabilités marginales P(Ij = 1)
p <- t(indices) %*% f

# 2. Probabilités conjointes P(Ia = 1, Ib = 1)
paires <- list(c(1, 2), c(1, 3), c(2, 3))

p11 <- t(
    sapply(
        paires, function(ab)
        {
            colSums(
                f[indices[, ab[1]] == 1 & indices[, ab[2]] == 1, , drop = FALSE]
            )
        }
    )
)

# Cov(Ia, Ib) = P(Ia = 1, Ib = 1) - P(Ia = 1) P(Ib = 1)
covariances <- p11 - rbind(
    p[1, ] * p[2, ],
    p[1, ] * p[3, ],
    p[2, ] * p[3, ]
)

rownames(covariances) <- c('Cov(I1,I2)', 'Cov(I1,I3)', 'Cov(I2,I3)')
round(covariances, 8)

# 3. Trois covariances >= 0
which(colSums(covariances <= 0) == 0)

# 4. Trois covariances <= 0
which(colSums(covariances >= 0) == 0)

# 5. Fmp de N = I1 + I2 + I3
n <- rowSums(indices)

fn <- t(sapply(0:3, function(k) colSums(f[n == k, , drop = FALSE])))
rownames(fn) <- paste0('', 0:3)
colSums(fn)

# 6. Espérance et variance
k <- 0:3
EspN <- drop(k %*% fn)
VarN <- drop(k^2 %*% fn) - EspN^2

data.frame(m = seq_len(ncol(f)), EspN, VarN)

# 7. Trouver le fmp indépendance
q <- c(1/4, 1/7, 1/3)
f.indep <- apply(indices, 1, function(x) prod(q^x * (1 - q)^(1 - x)))

EspN.indep <- sum(n * f.indep)
all.equal(EspN.indep, sum(q))

VarN.indep <- sum(n^2 * f.indep) - EspN.indep^2
all.equal(VarN.indep, sum(q * (1 - q)))

# 8–9. Préparation pour les bornes inférieure et supérieure de Fréchet
borne_sup <- t(
    sapply(paires, function(ab) pmin(p[ab[1], ], p[ab[2], ]))
)

borne_inf <- t(
    sapply(paires, function(ab) pmax(0, p[ab[1], ] + p[ab[2], ] - 1))
)

comonotones <- p11 == borne_sup
antimonotones <- p11 == borne_inf

# 8. Trois composantes comonotones
which(colSums(comonotones) == 3)

# 9. Deux paires antimonotones et une paire comonotone
which(colSums(antimonotones) == 2 & colSums(comonotones) == 1)

# 10. Mutuellement exclusives : P(N >= 2) = 0
which(colSums(fn[3:4, ]) == 0)

# 11. Variance maximale de N
which(max(VarN) == VarN)

# 12. Variance minimale de N
which(min(VarN) == VarN)

###
### Distributions Bernoulli multivariées, partie 2
###

probabilities_B <- data.frame(
  i = c('(0,0,0)', '(1,0,0)', '(0,1,0)', '(1,1,0)',
        '(0,0,1)', '(1,0,1)', '(0,1,1)', '(1,1,1)'),

  f10 = c(0.4047619, 0.1666666, 0.0952381, 0,
          0.2380952, 0.0476191, 0.0119048, 0.0357143),

  f11 = c(0.4404762, 0.1309524, 0.0595239, 0.0357143,
          0.2023809, 0.0833333, 0.0476190, 0),

  f12 = c(0.4285714, 0.1428571, 0.0714286, 0.0238095,
          0.2142857, 0.0714286, 0.0357143, 0.0119048)
)

# Préparation
B <- 10:12
f <- as.matrix(probabilities_B[, -1])

indices <- t(
    sapply(strsplit(gsub('[()]', '', probabilities_B$i), ','), as.numeric)
)

# Tolérance tenant compte des arrondis du tableau
tol <- 1e-6

# Validation et marginales P(Ij = 1)
colSums(f)
p <- t(indices) %*% f

# 1. Covariances
covariances <- rbind(
    'Cov(I1,I2)' = drop(
        (indices[, 1] * indices[, 2]) %*% f
    ) - p[1, ] * p[2, ],

    'Cov(I1,I3)' = drop(
        (indices[, 1] * indices[, 3]) %*% f
    ) - p[1, ] * p[3, ],

    'Cov(I2,I3)' = drop(
        (indices[, 2] * indices[, 3]) %*% f
    ) - p[2, ] * p[3, ]
)

round(covariances, 6)

# 2. Indépendance
independants <- sapply(
    seq_along(B), function(m) {
        produit <- apply(
            indices, 1, function(x) prod(
                p[, m]^x * (1 - p[, m])^(1 - x)
            )
        )
        all(abs(f[, m] - produit) < tol)
    }
)

B[independants]  # 12

# 3. Indépendance par paire
independants_par_paire <- colSums(abs(covariances) >= tol) == 0

B[independants_par_paire]  # 10, 11, 12

# Indépendance par paire SANS indépendance mutuelle
B[independants_par_paire & !independants]  # 10, 11

# 4. Fonction de masse de N
n <- rowSums(indices)

fn <- t(sapply(0:3, function(k) colSums(f[n == k, , drop = FALSE])))

rownames(fn) <- paste0('k = ', 0:3)
fn
