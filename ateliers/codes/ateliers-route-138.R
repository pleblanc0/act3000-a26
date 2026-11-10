################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### École d'actuariat - Université Laval ######################################
#### J. Bernoulli et A. Markov font un road trip sur la route 138 ##############
#### Auteur : Philippe Leblanc #################################################
################################################################################

###
### J. Bernoulli et A. Markov font un road trip sur la route 138
###

rm(list = ls())

# rappels sur les matrices
q <- 0.05; al <- 0.3
A = matrix(c(1 - q, q, 1 - q, q), 2, 2)
B = matrix(c(q, -q, -(1 - q), 1 - q), 2, 2)

# indempotence 
all.equal(A %*% A, A)
all.equal(B %*% B, B)

# commutativité
all.equal(A %*% B, B %*% A)

# matrice nulle
all.equal(A %*% B, matrix(0, 2, 2))

# matrice de transition
P = matrix(
    c(
        1 - (1 - al) * q, (1 - al) * q,
        (1 - al) * (1 - q), 1 - (1 - al) * (1 - q)
    ), 2, 2
)
all.equal(P, A + al * B)

# paramètres pour l'exercice
m <- 100; q <- 0.05; al <- c(0, 0.3, 0.5, 0.8)
n <- nextn(m, 2); k <- 0:(n - 1); e <- exp(-2i * pi * k/n)

# fgp conjointe récursive (PAS optimisée, simplfiée pour mieux comprendre)
fgpI <- function(s, q, alpha)
{
    k <- length(s)

    xi <- zeta <- numeric(k)
    xi[k] <- zeta[k] <- 1

    if (k >= 2)
    {
        for (v in (k - 1):1)
        {
            xi[v] <- (
                s[v + 1] * (q + alpha * (1 - q)) * xi[v + 1] +
                (1 - q) * (1 - alpha) * zeta[v + 1]
            )
            zeta[v] <- (
                s[v + 1] * q * (1 - alpha) * xi[v + 1] +
                (1 - q + alpha * q) * zeta[v + 1]
            )
        }
    }

    q * s[1] * xi[1] + (1 - q) * zeta[1]
}

###
### Calculs pour N100 = I1 + ... + I100
###

fgpNm <- function(t, q, al, m) fgpI(rep(t, m), q, al)

fN100 <- sapply(
    al, function(a) Re(fft(sapply(e, fgpNm, q = q, al = a, m = m), TRUE))/n
)
FN100 <- apply(fN100, 2, cumsum)

# Fonctions de répartition
v <- c(0, 1, 2, 3, 4, 5, 10, 20, 30, 50)
FN100[v + 1, ]

# Espérances
apply(fN100, 2, function(f) sum(k * f))

# Variances
apply(fN100, 2, function(f) sum((k - sum(k * f))^2 * f))

# Value at Risk à 90%
VaR <- function(u, f) k[min(which(cumsum(f) >= u))]
apply(fN100, 2, function(f) VaR(0.9, f))

# Tail Value at Risk à 90%
TVaR <- function(u, f) VaR(u, f) + 1/(1 - u) * sum(
    pmax(k - VaR(u, f), 0) * f
)
apply(fN100, 2, function(f) TVaR(0.9, f))

###
### Calculs pour W100 = S100/100, où S100 = X1 + ... + X100 et Xi = 100 Ii
###

b <- 100; m <- 10; al <- c(0.3, 0.8)
n <- nextn(m * b, 2); k <- 0:(n - 1); e <- exp(-2i * pi * k/n)

fgpWm <- function(t, q, al, m, b) fgpNm(t^(b/m), q, al, m)
# fgpWm <- function(t, q, al, m, b) fgpI(rep(t^(b/m), m), q, al)

fW10 <- sapply(
    al, function(a) Re(
        fft(sapply(e, fgpWm, q = q, al = a, m = m, b = b), TRUE)/n
    )
)

# Espérances
colSums(k * fW10)

# Variances
colSums((k - colSums(k * fW10))^2 * fW10)

# Value at Risk à 90%
apply(fW10, 2, function(f) VaR(0.9, f))

# Tail Value at Risk à 90%
apply(fW10, 2, function(f) TVaR(0.9, f))

###
### Contributions de Xi à l'écart-type de Sn
###

rho <- 0.3
VarXi <- 100^2 * q * (1 - q)

VarSn <- m * VarXi + 2 * VarXi * rho/(1 - rho) * (
    (m - 1) - rho * (1 - rho^(m - 1))/(1 - rho)
)

a <- VarXi
b <- VarXi * rho/(1 - rho)
cn <- sqrt(VarSn)

ContrXi <- sapply(
    1:m, function(i) (a + b * (2 - rho^(i - 1) - rho^(m - i)))/cn
)
cn
all.equal(sum(ContrXi), cn)

# observation: symétrie des contributions!
all.equal(ContrXi[1:(m/2)], rev(ContrXi[(m/2 + 1):m]))


###
### Calcul de la fmp de I par transformée de Walsh-Hadamard
###

# matrices de Hadamard construites récursivement
H2 <- matrix(c(1, 1, 1, -1), 2, 2)
H2k <- function(k) if (k == 1) H2 else kronecker(H2, H2k(k - 1))

# calcul de la transformée de I
evals <- expand.grid(rep(list(c(1, -1)), m))
dftI <- sapply(
    1:nrow(evals), function(i) fgpI(unlist(evals[i, ]), q = q, alpha = 0.3)
)

# calcul de la fmp complète de I
fI <- drop(H2k(m) %*% dftI)/2^m
sum(fI)

# vérification des espérances et variances
supp <- expand.grid(rep(list(0:1), m))
colSums(as.matrix(supp) * fI)
colSums((as.matrix(supp) - colSums(as.matrix(supp) * fI))^2 * fI)
