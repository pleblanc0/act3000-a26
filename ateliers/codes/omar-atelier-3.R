rm(list = ls())

n <- 2^10; k <- 0:(n - 1)
lam <- c(2, 3); a0 <- 1
m1 <- 0:50
m2 <- 0:50

#### Exemple no 1 ####

## Méthode 1

dPoTeibiv <- function(k1, k2, lam, a0){
  f <- function(i, j){
    k <- 0:min(i, j)
    sum(dpois(k, a0) * dpois(i - k, lam[1] - a0) * dpois(j - k, lam[2] - a0))
  }
  outer(k1, k2, Vectorize(f))
}
fm1m2 <- dPoTeibiv(m1, m2, lam, a0)
round(fm1m2[1:6,1:6],4)

## Méthode 2

# comparaison
e <- exp(-2i * pi * k/n)
e.alt <- fft(c(0, 1, rep(0, n - 2)))
all.equal(e, e.alt)

fgpM1M2 <- function(t1, t2) exp(
  a0 * (t1 * t2 - 1) + (lam[1] - a0) * (t1 - 1) + (lam[2] - a0) * (t2 - 1)
)
fm1m2.alt <- Re(fft(outer(e, e, fgpM1M2), TRUE))/n^2

all.equal(fm1m2,fm1m2.alt[1:51,1:51])

round(fm1m2[1:6,1:6],4)

#### Fonction de repartition ####
Fm1m2 <- t(apply(apply(fm1m2.alt, 2, cumsum), 1, cumsum))
round(Fm1m2[1:6,1:6],4)

# Alternative
f <- function(i, j) sum(fm1m2.alt[0:i + 1, 0:j + 1])
Fm1m2.alt <- outer(m1, m2, Vectorize(f))

all.equal(Fm1m2[m1+1,m2+1],Fm1m2.alt)


####  Loi de la somme  ####

#  Méthode 1 : La loi de N est connue

## Panjer
panjer.poisson <- function(lam, fb, smax)
{
  fs <- exp(lam * (fb[1] - 1))
  fb <- c(fb, numeric(smax + 1 - length(fb)))
  
  for(i in 1:smax) fs[i + 1] <- sum(
    fb[2:(i + 1)] * fs[i:1] * lam * (1:i)/i
  )
  fs
}

fn_1 <- function(a0) {
  lambdaN <- sum(lam) - a0
  
  fb <- c(0, (sum(lam) - 2*a0)/lambdaN,
          a0/lambdaN, rep(0, n - 3))
  
  panjer.poisson(lambdaN, fb, n - 1)
}
fn1 <- fn_1(a0)
sum (fn1)


## fft

fn_2 <- function(a0) {
  lambdaN <- sum(lam) - a0
  
  fb <- c(0, (sum(lam) - 2*a0)/lambdaN,
          a0/lambdaN, rep(0, n - 3))
  
  fnt <- exp(lambdaN * (fft(fb) - 1))
  
  Re(fft(fnt, TRUE))/n
}
fn2 <- fn_2(a0)
sum(fn2)

#  Méthode 2 : Directement par la fgp
fn_3 <- function(a0) {
  
  fgpN <- function(t) fgpM1M2(t, t)
  
  Re(fft(sapply(e, fgpN), TRUE))/n
}

fn3 <- fn_3(a0)
sum(fn3)

#  Méthode 3 : N = K1 + K2 + 2K0
fn_4 <- function(a0) {
  fk12 <- dpois(k, sum(lam - a0))
  
  f2k0 <- numeric(n)
  f2k0[k %% 2 == 0] <- dpois(k[k %% 2 == 0]/2, a0)
  
  fnt <- fft(fk12) * fft(f2k0)
  
  Re(fft(fnt, TRUE))/n
}

fn4 <- fn_4(a0)
sum(fn4)

all.equal(fn1, fn2)
all.equal(fn1, fn3)
all.equal(fn1, fn4)

# Cas indépendant : N ~ Poisson(sum(lam))
fn_ind <- dpois(k, sum(lam))

## Fonction de répartition
Fn_ind <- cumsum(fn_ind)
Fn <- cumsum(fn1)

## Prime stop-loss
pi_ind <- sapply(k, function(s)
  sum(pmax(k - s, 0) * fn_ind)
)

pi <- sapply(k, function(s)
  sum(pmax(k - s, 0) * fn1)
)

round(
  cbind(
    k,
    "fN_ind" = fn_ind,
    "FN_ind" = Fn_ind,
    "piN_ind" = pi_ind,
    "fN" = fn1,
    "FN" = Fn,
    "piN" = pi
  )[1:11, ],
  4
)


#### Exemple no 2 ####
lam <- c(5, 10); a <- 0:5

#### Comonotonne et Antimonotonne ####

### Fonction de répartition marginales
Fm1 <- c(0,ppois(m1, lam[1]))
Fm2 <- c(0,ppois(m2, lam[2]))

### Comonotonne

## Fonction de masse de probabilité conjointe
Fm1m2c <- outer(Fm1, Fm2, pmin)

fm1m2c <- t(apply(apply(Fm1m2c, 2, diff), 1, diff))

#Alternative

fm1m2c.alt1 <-  Fm1m2c[-1, -1] -
  Fm1m2c[-nrow(Fm1m2c), -1] -
  Fm1m2c[-1, -ncol(Fm1m2c)] +
  Fm1m2c[-nrow(Fm1m2c), -ncol(Fm1m2c)]

#Alternative
pcomo <- function(x1, x2) pmin(Fm1[x1 + 2], Fm2[x2 + 2])

f12c <- function(x1, x2) pcomo(x1, x2) - pcomo(x1-1, x2) - pcomo(x1, x2-1) + pcomo(x1-1, x2-1)

fm1m2c.alt2 <- outer(m1, m2, Vectorize(f12c))

all.equal(fm1m2c,fm1m2c.alt1)
all.equal(fm1m2c,fm1m2c.alt2)

## Espérance du produit
Em12c <- sum(outer(m1, m2, "*") * fm1m2c)

## Covariance et corréliation de Pearson
covc <- Em12c - prod(lam)
ppc <- covc/sqrt(prod(lam))
cbind(Em12c, covc, ppc)

## V.a. N

kmax <- max(m1)+max(m1)
kn <- 0:kmax
fnc <- sapply(kn,
              function(s) sum(fm1m2c[outer(m1, m2, "+") == s]))
sum(fnc)

### Antimonotone

## Fonction de masse de probabilité conjointe
Fm1m2a <- pmax(outer(Fm1, Fm2, "+") - 1, 0)

fm1m2a <- t(apply(apply(Fm1m2a, 2, diff), 1, diff))

#Alternative
panti <- function(x1, x2) pmax(Fm1[x1 + 2] + Fm2[x2 + 2] - 1, 0)

f12a <- function(x1, x2) panti(x1, x2) - panti(x1-1, x2) - panti(x1, x2-1) + panti(x1-1, x2-1)

fm1m2a <- outer(m1, m2, Vectorize(f12a))

## Espérance du produit
Em12a <- sum(outer(m1, m2, "*") * fm1m2a)

## Covariance et corréliation de Pearson
cova <- Em12a - prod(lam)
ppa <- covc/sqrt(prod(lam))
cbind(Em12a, cova, ppa)

## V.a. N

kmax <- max(m1)+max(m1)
kn <- 0:kmax
fna <- sapply(kn,
              function(s) sum(fm1m2a[outer(m1, m2, "+") == s]))
sum(fna)


## Covariance et corréliation de Pearson

CovM1M2 <- a
CorM1M2 <- CovM1M2 / sqrt(prod(lam))
cbind(a, CovM1M2, CorM1M2)

## Valeurs de k
ks <- c(0, 5, 10, 15, 20, 30, 40, 50)

## Loi de N pour a0 = 0, 1, 3, 5
fn <- sapply(c(0, 1, 3, 5), fn_1)

## Primes stop-loss
piN <- function(f, k) {
  sapply(k, function(s)
    sum(pmax(0:(length(f)-1) - s, 0) * f)
  )
}

round(cbind(
  k = ks,
  "piN-" = piN(fna, ks),
  "piN(0)" = piN(fn[, 1], ks),
  "piN(1)" = piN(fn[, 2], ks),
  "piN(3)" = piN(fn[, 3], ks),
  "piN(5)" = piN(fn[, 4], ks),
  "piN+" = piN(fnc, ks)
), 5)

