################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### École d'actuariat - Université Laval ######################################
#### Ateliers - Partie 1 #######################################################
#### Auteur : Philippe Leblanc #################################################
################################################################################

###
### Somme de v.a. indépendantes et méthode de force brute
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1)
m <- 10; i <- 1:m
r <- 0.25 * i; q <- r/(r + 2)

## 1)
EspM <- r * (1 - q)/q
EspN <- sum(EspM)

## 2)
VarM <- r * (1 - q)/q^2
VarN <- sum(VarM)

## 3)
fm <- matrix(0, nrow = m, ncol = n)
for (j in i)
    fm[j, ] <- dnbinom(k, r[j], q[j])

directconvo <- function(f1, f2)
{
    fs <- numeric(length(f1) + length(f2) - 1)

    for (i in seq_along(f1))
    {
        j <- i + seq_along(f2) - 1
        fs[j] <- fs[j] + f1[i] * f2
    }
    fs
}

directconvo.nrisks <- function(matff)
{
    rows <- lapply(seq_len(nrow(matff)), function(i) matff[i, ])
    Reduce(directconvo, rows)
}

fn <- directconvo.nrisks(fm)
s <- seq_along(fn) - 1

## 5)
VaR <- function(u) s[min(which(cumsum(fn) >= u))]
TVaR <- function(u) sum(pmax(s - VaR(u), 0) * fn)/(1 - u) + VaR(u)

u <- c(0.5, 0.9, 0.99, 0.999)
cbind(VaR = sapply(u, VaR), TVaR = sapply(u, TVaR))


###
### Somme de v.a. indépendantes, approche symbolique et fgp
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1)
q <- c(0.2, 0.3, 0.4, 0.1, 0.5); b <- c(5, 8, 7, 9, 2)

## 1)
EspX <- b * q
EspS <- sum(EspX)

## 2)
VarX <- b^2 * q * (1 - q)
VarS <- sum(VarX)

## 3)
fx <- matrix(0, nrow = 5, ncol = n)
fx[, 1] <- 1 - q
for (i in 1:5) fx[i, b[i] + 1] <- q[i]

directconvo <- function(f1, f2)
{
    fs <- numeric(length(f1) + length(f2) - 1)

    for (i in seq_along(f1))
    {
        j <- i + seq_along(f2) - 1
        fs[j] <- fs[j] + f1[i] * f2
    }
    fs
}

directconvo.nrisks <- function(matff)
{
    rows <- lapply(seq_len(nrow(matff)), function(i) matff[i, ])
    Reduce(directconvo, rows)
}

fs <- directconvo.nrisks(fx)
s <- seq_along(fs) - 1

## 5)
VaR <- function(u) s[min(which(cumsum(fs) >= u))]
TVaR <- function(u) sum(pmax(s - VaR(u), 0) * fs)/(1 - u) + VaR(u)

u <- c(0.5, 0.9, 0.99, 0.999)
cbind(VaR = sapply(u, VaR), TVaR = sapply(u, TVaR))


###
### Algorithme de De Pril - Exemple 1
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); m <- 10; h <- 1000
p <- c(0.5, 0.1, 0.2, 0.08, 0.05, 0.04, 0.02, 0.01)
b <- c(0, 1, 2, 5, 10, 20, 50, 100)

## 1)
EspX <- sum(b * p)
EspS <- m * EspX

## 2)
VarX <- sum(b^2 * p) - EspX^2
VarS <- m * VarX

## 3)
fx <- numeric(n)
fx[b + 1] <- p

depril <- function(fx, n, kmax)
{
    fs <- fx[1]^n
    fx <- c(fx, numeric(kmax + 1 - length(fx)))
    
    for(i in 1:kmax) fs[i + 1] <- 1/fx[1] * sum(
      fx[2:(i + 1)] * fs[i:1] * ((n + 1) * (1:i)/i - 1)
    )
    fs
}
# attention, les valeurs sont pour 1000k
fs <- depril(fx, m, n - 1)

## 5)
stoploss <- function(d) sum(pmax(k * h - d, 0) * fs)
sapply(0:10 * 10, stoploss)


###
### Algorithme de De Pril - Exemple 2
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); m <- 10
al <- c(0.8, 0.2); nu <- c(0.5, 8)

## 1)
EspX <- sum(al * nu)
EspS <- m * EspX

## 2)
VarX <- sum(al * (nu^2 + nu)) - EspX^2
VarS <- m * VarX

## 3)
depril <- function(fx, n, kmax)
{
    fs <- fx[1]^n
    fx <- c(fx, numeric(kmax + 1 - length(fx)))
    
    for(i in 1:kmax) fs[i + 1] <- 1/fx[1] * sum(
      fx[2:(i + 1)] * fs[i:1] * ((n + 1) * (1:i)/i - 1)
    )
    fs
}
fx <- al[1] * dpois(k, nu[1]) + al[2] * dpois(k, nu[2])
fs <- depril(fx, m, n - 1)

## 5)
stoploss <- function(d) sum(pmax(k - d, 0) * fs)
sapply(0:10 * 10, stoploss)


###
### Approche symbolique et fgp
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); m <- 10; l <- 0:m
al <- c(0.8, 0.2); nu <- c(0.5, 8); gam <- l * nu[1] + (m - l) * nu[2]

fmp <- function(k) sum(
    choose(m, l) * (al[1]^l * al[2]^(m - l)) * dpois(k, gam)
)

## 3)
fs <- sapply(k, fmp)
fs[0:10 + 1]

## 5)
stoploss <- function(d) sum(pmax(k - d, 0) * fs)
sapply(0:10 * 10, stoploss)


###
### Algorithme de De Pril et translation
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
r <- 4; nu <- 0.4; m <- 10

## 1)
EspX <- r/nu
EspS <- m * EspX

## 2)
VarX <- r * (1 - nu)/nu^2
VarS <- m * VarX

## 3)
depril <- function(fx, n, kmax)
{
    fs <- fx[1]^n
    fx <- c(fx, numeric(kmax + 1 - length(fx)))
    
    for(i in 1:kmax) fs[i + 1] <- 1/fx[1] * sum(
      fx[2:(i + 1)] * fs[i:1] * ((n + 1) * (1:i)/i - 1)
    )
    fs
}
fx <- c(rep(0, r), dnbinom(k, r, nu))[1:n]

fs <- depril(fx[(r + 1):n], m, n)
fs <- c(rep(0, r * m), fs)[1:n]

## 5)
stoploss <- function(d) sum(pmax(k - d, 0) * fs)
sapply(0:10 * 10, stoploss)


###
### Algorithme de Panjer - Exemple 1
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1)
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
panjer.poisson <- function(lam, fb, smax)
{
    fs <- exp(lam * (fb[1] - 1))
    fb <- c(fb, numeric(smax + 1 - length(fb)))

    for(i in 1:smax) fs[i + 1] <- sum(
        fb[2:(i + 1)] * fs[i:1] * lam * (1:i)/i
    )
    fs
}

fb <- c(0, sapply(k[-1], function(k) sum(al * nu * (1 - nu)^(k - 1))))
fx <- panjer.poisson(lam, fb, n - 1)

## 4)
VaR <- function(u) k[min(which(cumsum(fx) >= u))]
TVaR <- function(u) sum(pmax(k - VaR(u), 0) * fx)/(1 - u) + VaR(u)

u <- c(0.5, 0.9, 0.99, 0.999)
cbind(VaR = sapply(u, VaR), TVaR = sapply(u, TVaR))

## 5)
stoploss <- function(d) sum(pmax(k - d, 0) * fx)
sapply(2:10 * 10, stoploss)


###
### Algorithme de Panjer - Exemple 2
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1)
r <- 2.5; q <- 1/3; al <- c(0.8, 0.2); nu <- c(0.5, 0.1)

## 1)
EspM <- r * (1 - q)/q
EspB <- sum(al/nu)
EspX <- EspM * EspB

## 2) Attention de ne pas faire une somme pondérée des variances pour VarB!
VarM <- r * (1 - q)/q^2
VarB <- sum(al * (2 - nu)/nu^2) - sum(al/nu)^2
VarX <- EspM * VarB + VarM * EspB^2

## 3)
panjer.nbinom <- function(r, q, fb, smax)
{
    a <- 1 - q; b <- a * (r - 1)
    fs <- (q/(1 - a * fb[1]))^r
    fb <- c(fb, numeric(smax + 1 - length(fb)))
    
    for(i in 1:smax) fs[i + 1] <- 1/(1 - a * fb[1]) * sum(
        fb[2:(i + 1)] * fs[i:1] * (b * (1:i)/i + a)
    )
    fs
}

fb <- c(0, sapply(k[-1], function(k) sum(al * nu * (1 - nu)^(k - 1))))
fx <- panjer.nbinom(r, q, fb, n - 1)

## 4)
VaR <- function(u) k[min(which(cumsum(fx) >= u))]
TVaR <- function(u) sum(pmax(k - VaR(u), 0) * fx)/(1 - u) + VaR(u)

u <- c(0.5, 0.9, 0.99, 0.999)
cbind(VaR = sapply(u, VaR), TVaR = sapply(u, TVaR))

## 5)
stoploss <- function(d) sum(pmax(k - d, 0) * fx)
sapply(2:10 * 10, stoploss)


###
### Algorithme de Panjer - Exemple 3
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1)
r <- 0.5; q <- 1/11; al <- c(0.8, 0.2); nu <- c(0.5, 0.1)

## 1)
EspM <- r * (1 - q)/q
EspB <- sum(al/nu)
EspX <- EspM * EspB

## 2) Attention de ne pas faire une somme pondérée des variances pour VarB!
VarM <- r * (1 - q)/q^2
VarB <- sum(al * (2 - nu)/nu^2) - sum(al/nu)^2
VarX <- EspM * VarB + VarM * EspB^2

## 3)
panjer.nbinom <- function(r, q, fb, smax)
{
    a <- 1 - q; b <- a * (r - 1)
    fs <- (q/(1 - a * fb[1]))^r
    fb <- c(fb, numeric(smax + 1 - length(fb)))
    
    for(i in 1:smax) fs[i + 1] <- 1/(1 - a * fb[1]) * sum(
        fb[2:(i + 1)] * fs[i:1] * (b * (1:i)/i + a)
    )
    fs
}

fb <- c(0, sapply(k[-1], function(k) sum(al * nu * (1 - nu)^(k - 1))))
fx <- panjer.nbinom(r, q, fb, n - 1)

## 4)
VaR <- function(u) k[min(which(cumsum(fx) >= u))]
TVaR <- function(u) sum(pmax(k - VaR(u), 0) * fx)/(1 - u) + VaR(u)

u <- c(0.5, 0.9, 0.99, 0.999)
cbind(VaR = sapply(u, VaR), TVaR = sapply(u, TVaR))

## 5)
stoploss <- function(d) sum(pmax(k - d, 0) * fx)
sapply(2:10 * 10, stoploss)


###
### Méli-mélo - Exemple 1
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1)
r <- c(0.5, 2.5); q <- c(1/11, 1/3)
a <- c(0.8, 0.2); nu <- c(0.5, 0.1)
b <- c(0.7, 0.3); eta <- c(0.25, 0.0625)

## 1)
EspM <- r * (1 - q)/q
EspB <- c(sum(a/nu), sum(b/eta))
EspX <- EspM * EspB
EspS <- sum(EspX)

## 2) Attention de ne pas faire une somme pondérée des variances pour VarB!
VarM <- r * (1 - q)/q^2
VarB <- c(
    sum(a * (2 - nu)/nu^2) - sum(a/nu)^2,
    sum(b * (2 - eta)/eta^2) - sum(b/eta)^2
)
VarX <- EspM * VarB + VarM * EspB^2
VarS <- sum(VarX)

## 3)
panjer.nbinom <- function(r, q, fb, smax)
{
    a <- 1 - q; b <- a * (r - 1)
    fs <- (q/(1 - a * fb[1]))^r
    fb <- c(fb, numeric(smax + 1 - length(fb)))
    
    for(i in 1:smax) fs[i + 1] <- 1/(1 - a * fb[1]) * sum(
        fb[2:(i + 1)] * fs[i:1] * (b * (1:i)/i + a)
    )
    fs
}

fb1 <- c(0, sapply(k[-1], function(k) sum(a * nu * (1 - nu)^(k - 1))))
fb2 <- c(0, sapply(k[-1], function(k) sum(b * eta * (1 - eta)^(k - 1))))

fx1 <- panjer.nbinom(r[1], q[1], fb1, n - 1)
fx2 <- panjer.nbinom(r[2], q[2], fb2, n - 1)

directconvo <- function(f1, f2)
{
    fs <- numeric(length(f1) + length(f2) - 1)

    for (i in seq_along(f1))
    {
        j <- i + seq_along(f2) - 1
        fs[j] <- fs[j] + f1[i] * f2
    }
    fs
}

fs <- directconvo(fx1, fx2)[1:n]

## 4)
VaR <- function(u) k[min(which(cumsum(fs) >= u))]
TVaR <- function(u) sum(pmax(k - VaR(u), 0) * fs)/(1 - u) + VaR(u)

u <- c(0.5, 0.9, 0.99, 0.999)
cbind(VaR = sapply(u, VaR), TVaR = sapply(u, TVaR))

## 5)
stoploss <- function(d) sum(pmax(k - d, 0) * fs)
sapply(5:15 * 10, stoploss)


###
### Méli-mélo - Exemple 2
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
r <- c(0.5, 2.5); q <- c(1/11, 1/3)
a <- c(0.8, 0.2); nu <- c(0.5, 0.1)
b <- c(0.7, 0.3); eta <- c(0.25, 0.0625)

fc1 <- c(0, sapply(k[-1], function(k) sum(a * nu * (1 - nu)^(k - 1))))
fc2 <- c(0, sapply(k[-1], function(k) sum(b * eta * (1 - eta)^(k - 1))))


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
panjer.nbinom <- function(r, q, fb, smax)
{
    a <- 1 - q; b <- a * (r - 1)
    fs <- (q/(1 - a * fb[1]))^r
    fb <- c(fb, numeric(smax + 1 - length(fb)))
    
    for(i in 1:smax) fs[i + 1] <- 1/(1 - a * fb[1]) * sum(
        fb[2:(i + 1)] * fs[i:1] * (b * (1:i)/i + a)
    )
    fs
}
fx1 <- panjer.nbinom(r[1], q[1], fb1, n - 1)
fx2 <- panjer.nbinom(r[2], q[2], fb2, n - 1)

directconvo <- function(f1, f2)
{
    fs <- numeric(length(f1) + length(f2) - 1)

    for (i in seq_along(f1))
    {
        j <- i + seq_along(f2) - 1
        fs[j] <- fs[j] + f1[i] * f2
    }
    fs
}
fs <- directconvo(fx1, fx2)[1:n]

## 4)
VaR <- function(u) k[min(which(cumsum(fs) >= u))]
TVaR <- function(u) sum(pmax(k - VaR(u), 0) * fs)/(1 - u) + VaR(u)

u <- c(0.5, 0.9, 0.99, 0.999)
cbind(VaR = sapply(u, VaR), TVaR = sapply(u, TVaR))

## 5)
stoploss <- function(d) sum(pmax(k - d, 0) * fs)
sapply(5:15 * 10, stoploss)


###
### Immeuble avec 4 unités en rangée et péril incendie
###

rm(list = ls())

n <- 2^8; k <- 0:(n - 1); e <- exp(2i * pi * k/n); b <- 1:4
lam <- 0.4; a12 <- 0.5; a23 <- 0.7; a34 <- 0.2

fbi1 <- c(1 - a12, a12 * (1 - a23), a12 * a23 * (1 - a34), a12 * a23 * a34)
fbi2 <- c(
    (1 - a12) * (1 - a23), a12 * (1 - a23) + (1 - a12) * a23 * (1 - a34),
    a12 * a23 * (1 - a34) + (1 - a12) * a23 * a34, a12 * a23 * a34
)
fbi3 <- c(
    (1 - a23) * (1 - a34), a34 * (1 - a23) + (1 - a34) * a23 * (1 - a12),
    a23 * a34 * (1 - a12) + (1 - a34) * a23 * a12, a23 * a34 * a12
)
fbi4 <- c(1 - a34, a34 * (1 - a23), a34 * a23 * (1 - a12), a34 * a23 * a12)
fb <- (fbi1 + fbi2 + fbi3 + fbi4)/4

EspB <- sum(b * fb)
VarB <- sum(b^2 * fb) - EspB^2

EspX <- lam * EspB
VarX <- lam * (VarB + EspB^2)

panjer.poisson <- function(lam, fb, smax)
{
    fs <- exp(lam * (fb[1] - 1))
    fb <- c(fb, numeric(smax + 1 - length(fb)))

    for(i in 1:smax) fs[i + 1] <- sum(
        fb[2:(i + 1)] * fs[i:1] * lam * (1:i)/i
    )
    fs
}

fb <- c(0, fb, rep(0, n - length(fb) - 1))
fx <- panjer.poisson(lam, fb, n - 1)
