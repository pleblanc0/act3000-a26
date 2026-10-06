################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### École d'actuariat - Université Laval ######################################
#### Ateliers - Distributions Bernoulli multivariées ###########################
#### Auteur : Philippe Leblanc #################################################
################################################################################

###
### Première partie
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
### Deuxième partie
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
