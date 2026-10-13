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

###
### Distribution discrète et simulation
###

m <- seq(0, 5)

Fm1m2 <- matrix(
    c(0.12, 0.18, 0.18, 0.18, 0.24, 0.30,
      0.12, 0.24, 0.24, 0.27, 0.34, 0.40,
      0.12, 0.42, 0.52, 0.60, 0.74, 0.80,
      0.16, 0.46, 0.56, 0.64, 0.81, 0.90,
      0.18, 0.48, 0.58, 0.66, 0.83, 0.96,
      0.20, 0.50, 0.60, 0.68, 0.85, 1.00),
    byrow = TRUE, nrow = length(m)
)

## 1)
Fm1 <- Fm1m2[, length(m)]
Fm2 <- Fm1m2[length(m), ]

## 2)
fm1 <- diff(c(0, Fm1))
fm2 <- diff(c(0, Fm2))

EspM1 <- sum(m * fm1)
EspM2 <- sum(m * fm2)

VarM1 <- sum(m^2 * fm1) - EspM1^2 
VarM2 <- sum(m^2 * fm2) - EspM2^2

## 3)
Fm1m2.f <- function(m1, m2) ifelse(
    m1 < 0 | m2 < 0, 0, Fm1m2[m1 + 1, m2 + 1]
)

fm1m2.f <- function(m1, m2)
{
    Fm1m2.f(m1, m2) - Fm1m2.f(m1 - 1, m2) - 
    Fm1m2.f(m1, m2 - 1) + Fm1m2.f(m1 - 1, m2 - 1) 
}

grid.m1m2 <- expand.grid(m1 = m, m2 = m)

fm1m2 <- matrix(
    mapply(fm1m2.f, grid.m1m2$m1, grid.m1m2$m2), ncol = length(m), byrow = FALSE
)

sum(fm1m2)
rowSums(fm1m2)
colSums(fm1m2)

## 4)
EspM1M2 <- sum(outer(m, m, function(x, y) x * y * fm1m2))
CovM1M2 <- EspM1M2 - EspM1 * EspM2

## 5)
fM1cM2 <- apply(fm1m2, 2, function(col) col/sum(col))
colSums(fM1cM2)

## 6)
Fm1cm2 <- apply(fM1cM2, 2, cumsum)

## 7)
fm2cm1 <- t(apply(fm1m2, 1, function(row) row/sum(row)))
rowSums(fm2cm1)

## 8)
Fm2cm1 <- t(apply(fm2cm1, 1, cumsum))

## 9)
EspM1cM2 <- function(m2) sum(mm * fm1cm2[, m2 + 1])

## 10)
EspM2cM1 <- function(m1) sum(mm * fm2cm1[m1 + 1, ])

## 11)
V1 <- c(0.334, 0.952, 0.016, 0.781)
V2 <- c(0.412, 0.859, 0.535, 0.186)

M1 <- M2 <- numeric(length(V1))

for (j in seq(length(M1)))
{  
    M1[j] <- m[min(which(cumsum(fm1) >= V1[j]))]
    M2[j] <- m[min(which(cumsum(fm2cm1[M1[j] + 1, ]) >= V2[j]))]
}

cbind(M1, M2)


###
### Chaîne de Markov homogène en temps discret
###

rm(list = ls())

tran <- matrix(
    c(0.8, 0.15, 0.05, 0.3, 0.65, 0.05, 0.3, 0.15, 0.55),
    byrow = TRUE, nrow = 3
)
rowSums(tran)

## 1)
stat <- eigen(t(tran))
stat <- stat$vectors[, 1]
stat <- stat/sum(stat)

## 2)
m <- 0:2
EspMk <- sum(m * stat)
VarMk <- sum(m^2 * stat) - EspMk^2

## 3)
fm1m2 <- function(m1, m2) stat[m1 + 1] * tran[m1 + 1, m2 + 1]

grid2 <- expand.grid(m1 = m, m2 = m)
cbind(grid2, prob = mapply(fm1m2, grid2$m1, grid2$m2))

## 4)
fm1m2m3 <- function(m1, m2, m3) (
    stat[m1 + 1] * tran[m1 + 1, m2 + 1] *  tran[m2 + 1, m3 + 1]
)

grid3 <- expand.grid(m1 = m, m2 = m, m3 = m)
cbind(grid3, prob = mapply(fm1m2m3, grid3$m1, grid3$m2, grid3$m3))

## 5)
fm1m2m3m4 <- function(m1, m2, m3, m4) (
    stat[m1 + 1] * tran[m1 + 1, m2 + 1] *
    tran[m2 + 1, m3 + 1] * tran[m3 + 1, m4 + 1]
)
  
grid4 <- expand.grid(m1 = m, m2 = m, m3 = m, m4 = m)
cbind(grid4, prob = mapply(fm1m2m3m4, grid4$m1, grid4$m2, grid4$m3, grid4$m4))

## 6)
EspM1M2 <- sum(outer(m, m, function(x, y) x * y * mapply(fm1m2, x, y))) 
CovM1M2 <- EspM1M2 - EspMk^2

# Covariance avec distance de 2
EspM1M3 <- sum(
    mapply(function(a, b, c) a * c * fm1m2m3(a, b, c),
    grid3$m1, grid3$m2, grid3$m3)
)
CovM1M3 <- EspM1M3 - EspMk^2

# Covariance avec distance de 3
EspM1M4 <- sum(
    mapply(function(a, b, c, d) a * d * fm1m2m3m4(a, b, c, d),
    grid4$m1, grid4$m2, grid4$m3, grid4$m4)
)
CovM1M4 <- EspM1M4 - EspMk^2

## 7)
EspN4 <- 4 * EspMk
VarN4 <- 4 * VarMk + 6 * CovM1M2 + 4 * CovM1M3 + 2 * CovM1M4

fgpL <- function(t) stat[1] + stat[2] * t + stat[3] * t^2

fgp.Ni <- function(i, t1, t2, theta)
{
    if (i == 1) fgpL(t2)
    else (
        (1 - theta) * fgp.Ni(i - 1, t1, t1, theta) * fgpL(t2) +
        theta * fgp.Ni(i - 1, t1, t1 * t2, theta)
    )
}

theta <- 0.5

# Fonction de masse de N20
n <- 2^8; k <- 0:(n - 1); r <- 20
e <- exp(-2i * pi * k/n)

fn20 <- Re(fft(fgp.Ni(r, e, e, theta), TRUE))/n
sum(fn20)

EspN20 <- sum(k * fn20) # 20 * E_Mk
VarN20 <- sum(k^2 * fn20) - EspN20^2

# Simulation d'un parcours
set.seed(42)
U <- runif(r)
parcours <- numeric(r)

parcours[1] <- min(which(cumsum(stat) >= U[1])) - 1

for(i in seq(2, r)) (
    parcours[i] <- min(which(cumsum(tran[parcours[i - 1] + 1, ]) >= U[i])) - 1  
)
parcours

# simulation de 100000 parcours pour approximer la fmp de N20
m <- 1E5; r <- 20; set.seed(42)
U <- matrix(runif(m * r), ncol = r, byrow = TRUE)

parcours <- matrix(0, nrow = m, ncol = r)
parcours[, 1] <- apply(
    U[, 1, drop = FALSE], 1, function(u) min(which(cumsum(stat) >= u)) - 1
)
for (i in 2:r) parcours[, i] <- mapply(
    function(prev, u) min(
        which(cumsum(tran[prev + 1, ]) >= u)
    ) - 1, parcours[, i - 1], U[, i]
)

# Cov(M1, M2)
c(cov(parcours[, 1], parcours[, 2]), CovM1M2)

# Cov(M1, M4)
c(cov(parcours[, 1], parcours[, 4]), CovM1M4)

