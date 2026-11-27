################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### École d'actuariat - Université Laval ######################################
#### Ateliers - Exemple copules #######################################################
#### Auteur : Omar Essakine #################################################
################################################################################

###
### Exemple (page 154)
###

rm(list = ls())


alpha <- 2
beta1 <- 1/10
beta2 <- 1/10

# fonctions de répartition marginales
F1 <- function(x) pexp(x, beta1)
F2 <- function(x) pexp(x, beta2)

# copule de Clayton
C.clayton <- function(u1, u2,alpha)
  (u1^(-alpha) + u2^(-alpha) - 1)^(-1/alpha)

# copule de Frank
C.frank <- function(u1, u2,alpha)
  -1/alpha * log(1 + ((exp(-alpha*u1) - 1) *
                        (exp(-alpha*u2) - 1)) / (exp(-alpha) - 1))

# copule de Gumbel
C.gumbel <- function(u1, u2,alpha)
  exp(-((-log(u1))^alpha + (-log(u2))^alpha)^(1/alpha))

# fonctions de répartition bivariées)
F.clayton <- function(x1, x2) C.clayton(F1(x1), F2(x2),alpha)
F.frank <- function(x1, x2) C.frank(F1(x1), F2(x2),alpha)
F.gumbel <- function(x1, x2) C.gumbel(F1(x1), F2(x2),alpha)

### Calcul du gamma

# copule de Clayton
F.clayton(10, 15) -
  F.clayton(5, 15) -
  F.clayton(10, 7) +
  F.clayton(5, 7)

# copule de Frank
F.frank(10, 15) -
  F.frank(5, 15) -
  F.frank(10, 7) +
  F.frank(5, 7)

# copule de Gumbel
F.gumbel(10, 15) -
  F.gumbel(5, 15) -
  F.gumbel(10, 7) +
  F.gumbel(5, 7)

# fonctions de densité marginales
f1 <- function(x) dexp(x, beta1)
f2 <- function(x) dexp(x, beta2)

# densité de la copule de Clayton
c.clayton <- function(u1, u2,alpha)
  (1 + alpha) * (u1*u2)^(-alpha - 1) *
  (u1^(-alpha) + u2^(-alpha) - 1)^(-2 - 1/alpha)

# densité de la copule de Frank
c.frank <- function(u1, u2,alpha)
  alpha * exp(-alpha * (u1 + u2)) * (1 - exp(-alpha)) /
  (exp(-alpha * (u1 + u2)) - exp(-alpha * u1) -
     exp(-alpha * u2) + exp(-alpha))^2

# densité de la copule de Gumbel
c.gumbel <- function(u1, u2,alpha)
  C.gumbel(u1, u2,alpha) *
  ((-log(u1))^(alpha - 1) * (-log(u2))^(alpha - 1)) / (u1 * u2) *
  ((-log(u1))^alpha + (-log(u2))^alpha)^(1/alpha - 2) *
  (alpha - 1 + ((-log(u1))^alpha + (-log(u2))^alpha)^(1/alpha))

# fonctions de densité bivariées pour x1, x2 > 0
f.clayton <- function(x1, x2)
  c.clayton(F1(x1), F2(x2),alpha)*f1(x1)*f2(x2)

f.frank <- function(x1, x2)
  c.frank(F1(x1), F2(x2),alpha)*f1(x1)*f2(x2)

f.gumbel <- function(x1, x2)
  c.gumbel(F1(x1), F2(x2),alpha)*f1(x1)*f2(x2)

c(f.clayton(10, 15), f.frank(10, 15), f.gumbel(10, 15))

# copule conditionnelle de Clayton
C.clayton2c1 <- function(u1, u2,alpha)
  u1^(-alpha - 1) *
  (u1^(-alpha) + u2^(-alpha) - 1)^(-1/alpha - 1)

# copule conditionnelle de Frank
C.frank2c1 <- function(u1, u2,alpha)
  exp(-alpha*u1)*(exp(-alpha*u2) - 1) /
  (exp(-alpha) - 1 +
     (exp(-alpha*u1) - 1)*(exp(-alpha*u2) - 1))

# copule conditionnelle de Gumbel
C.gumbel2c1 <- function(u1, u2,alpha)
  C.gumbel(u1, u2,alpha) * (-log(u1))^(alpha - 1) / u1 *
  ((-log(u1))^alpha + (-log(u2))^alpha)^(1/alpha - 1)

# fonctions de répartition conditionnelles
F.clayton2c1 <- function(x1, x2) C.clayton2c1(F1(x1), F2(x2),alpha)
F.frank2c1 <- function(x1, x2) C.frank2c1(F1(x1), F2(x2),alpha)
F.gumbel2c1 <- function(x1, x2) C.gumbel2c1(F1(x1), F2(x2),alpha)

c(F.clayton2c1(10, 15),
  F.frank2c1(10, 15),
  F.gumbel2c1(10, 15))


###
### Simulation de distributions exponentielles bivariées
###
# initialisation pour reproductibilité
set.seed(20261127)
nsim <- 1000

### Copule Clayton (methode conditionelle)

# on produit un seul échantillon de U1
U1 <- runif(nsim)

# on produit un seul échantillon de V2
V2 <- runif(nsim)

# on produit deux échantillons de U2, un pour chaque valeur de alpha
U2.clayton <- ((V2 * U1^(alpha + 1))^(-alpha/(alpha + 1)) -
                 U1^(-alpha) + 1)^(-1/alpha)

# on utilise la méthode inverse pour les marginales exponentielles
X1.clayton <- qexp(U1, beta1)
X2.clayton <- qexp(U2.clayton, beta2)

###
### Copule de Frank (méthode copule conditionnelle)
###


# on simule U2 à partir des mêmes observations de U1 et V2
U2.frank <- -1/alpha * log(1 +
                             (V2 * (exp(-alpha) - 1)) /
                             (exp(-alpha * U1) * (1 - V2) + V2))

# on utilise la méthode inverse pour les marginales exponentielles
X1.frank <- qexp(U1, beta1)
X2.frank <- qexp(U2.frank, beta2)

# Verification des valeurs
j <- c(1, 2, nsim)

cbind(
  j = j,
  X1.clayton = X1.clayton[j],
  X2.clayton = X2.clayton[j],
  X1.frank = X1.frank[j],
  X2.frank = X2.frank[j]
)

# simulation de S = X1 + X2
S.clayton <- X1.clayton + X2.clayton
S.frank <- X1.frank + X2.frank

# on affiche les réalisations de S
cbind(
  j = j,
  S.clayton = S.clayton[j],
  S.frank = S.frank[j]
)

# approximation de F_S(50)
FS.clayton <- mean(S.clayton <= 50)
FS.frank <- mean(S.frank <= 50)

c(Clayton = FS.clayton, Frank = FS.frank)

kappa <- 0.95

S.clayton.sort <- sort(S.clayton)
S.frank.sort <- sort(S.frank)

k <- ceiling(nsim * kappa)

# approximation de la VaR
VaR.clayton <- S.clayton.sort[k]
VaR.frank <- S.frank.sort[k]

# approximation de la TVaR
TVaR.clayton <- mean(S.clayton.sort[(k + 1):nsim])
TVaR.frank <- mean(S.frank.sort[(k + 1):nsim])

rbind(
  Clayton = c(VaR = VaR.clayton, TVaR = TVaR.clayton),
  Frank = c(VaR = VaR.frank, TVaR = TVaR.frank)
)

###
### Methode des rectangles
###

Fy1 <- function(y1) {
  pexp(y1, rate = 0.1)
}

Fy2 <- function(y2) {
  pexp(y2, rate = 0.1)
}

Fy1y2 <- function(y1, y2) {
  C.clayton(Fy1(y1), Fy2(y2),alpha=4)
}

CalculerRepartition <- function(s, m=20) {
  n <- 2^m
  h <- s/n
  a <- (0:(n-1))*h
  b <- (1:n)*h
  lower <- sum(Fy1y2(b, s-b) - Fy1y2(a, s-b))
  upper <- sum(Fy1y2(b, s-a) - Fy1y2(a, s-a))
  c(Lower = lower, Upper = upper)
}


s <- c(5, 10, 20, 30, 50, 100)
resultats <- data.frame( s = s,
                         t(sapply(s, CalculerRepartition))
)

print(resultats, digits = 7)

# Copule de Gumbel

Fy1y2 <- function(y1, y2) {
  C.gumbel(Fy1(y1), Fy2(y2),alpha=2)
}
s <- c(5, 10, 20, 30, 50, 100)
resultats <- data.frame(
  s = s,
  t(sapply(s, CalculerRepartition))
)

print(resultats, digits = 7)



# Copule de Clayton survie

C.ClaytonSurvie <- function(u, v, alpha=2) {
  u + v - 1 + ((1-u)^(-alpha) +
                 (1-v)^(-alpha) - 1)^(-1/alpha)
}

Fy1y2 <- function(y1, y2) {
  C.ClaytonSurvie(Fy1(y1), Fy2(y2))
}
s <- c(5, 10, 20, 30, 50, 100)
resultats <- data.frame(
  s = s,
  t(sapply(s, CalculerRepartition))
)
print(resultats, digits = 7)

# approximation de la VaR

s <- seq(0, 100, by = 0.1)
kappa <- 0.95
vm <- c(2, 5, 8, 10, 12)
VaR <- function(u, Fs) s[min(which(Fs >= u))]

Fy1y2 <- function(y1, y2) {
  C.clayton(Fy1(y1), Fy2(y2),alpha=4)
}

# calcul des bornes de la VaR pour chaque valeur de m
resultats.VaR <- t(sapply(vm, function(m) {
  
  FS <- t(sapply(s, CalculerRepartition, m = m))
  
  c(Lower = VaR(kappa, FS[, "Lower"]),
    Upper = VaR(kappa, FS[, "Upper"]))
}))

# comparaison des résultats
cbind(
  m = vm,
  resultats.VaR)


###
### Copule Clayton et marginales binomiales
###

n <- c(5, 5)
q <- c(0.2, 0.3)
alpha <- 5

m1 <- 0:n[1]
m2 <- 0:n[2]

### Fonctions de répartition marginales
Fm1 <- c(0, pbinom(m1, n[1], q[1]))
Fm2 <- c(0, pbinom(m2, n[2], q[2]))


### Fonction de répartition conjointe
Fm1m2 <- outer(Fm1, Fm2, C.clayton, alpha = alpha)

### Fonction de masse de probabilité conjointe
fm1m2 <- t(apply(apply(Fm1m2, 2, diff), 1, diff))

# Alternative 
fm1m2.alt1 <- Fm1m2[-1, -1] -
  Fm1m2[-nrow(Fm1m2), -1] -
  Fm1m2[-1, -ncol(Fm1m2)] +
  Fm1m2[-nrow(Fm1m2), -ncol(Fm1m2)]

### V.a. N

kmax <- max(m1) + max(m2)
kn <- 0:kmax

fn <- sapply(kn,
             function(s) sum(fm1m2[outer(m1, m2, "+") == s]))

# vérification
sum(fn)

# fonction de masse de probabilité de N
round(cbind(k = kn, fn = fn), 7)

### Espérances et variances marginales
EM <- n*q
VarM <- n*q*(1-q)

### Espérance du produit
Em12 <- sum(outer(m1, m2, "*") * fm1m2)

### Covariance et corrélation de Pearson
cov12 <- Em12 - prod(EM)
pp <- cov12/sqrt(prod(VarM))

cbind(Em12, cov12, pp)

# corrélation de Pearson en fonction de alpha
CorM1M2 <- function(alpha,copule) {
  
  # fonction de répartition conjointe
  Fm1m2 <- outer(Fm1, Fm2, copule, alpha = alpha)
  
  # fonction de masse de probabilité conjointe
  fm1m2 <- t(apply(apply(Fm1m2, 2, diff), 1, diff))
  
  # espérance du produit
  Em12 <- sum(outer(m1, m2, "*") * fm1m2)
  
  # covariance
  cov12 <- Em12 - prod(EM)
  
  # corrélation de Pearson
  cov12/sqrt(prod(VarM))
}

# vérification pour alpha = 5
CorM1M2(5,C.clayton)

# valeurs souhaitées de la corrélation de Pearson
rho <- c(0.02, 0.5, 0.8)

# on trouve alpha pour chaque valeur de rho
alpha <- sapply(rho, function(r) {
  uniroot(function(a) CorM1M2(a,C.clayton) - r,
          interval = c(0.0001, 20),
          tol = 1e-10)$root
})

# vérification des résultats
cbind(
  rho = rho,
  alpha = alpha,
  rho.calcul = sapply(alpha, CorM1M2, copule = C.clayton)
)

###
### Copules archimédiennes et marginales binomiales négatives
###

alpha <- c(3, 2)
q <- c(1/6, 1/9)
rho <- 0.5

# Espérances et variances marginales
EM <- alpha*(1-q)/q
VarM <- alpha*(1-q)/q^2

# Covariance 
CovM1M2 <- rho*sqrt(prod(VarM))

# Espérance et variance de N
EN <- sum(EM)
VarN <- sum(VarM) + 2*CovM1M2

c(EN = EN, VarN = VarN)

m1 <- 0:200
m2 <- 0:200

Fm1 <- c(0, pnbinom(m1, alpha[1], q[1]))
Fm2 <- c(0, pnbinom(m2, alpha[2], q[2]))

rho <- 0.5

# Copule de Frank
theta.frank <- uniroot(
  function(theta) CorM1M2(theta, C.frank) - rho,
  interval = c(0.001, 15), tol = 1e-10
)$root

# Copule de Clayton
theta.clayton <- uniroot(
  function(theta) CorM1M2(theta, C.clayton) - rho,
  interval = c(0.001, 15), tol = 1e-10
)$root

# Copule de Gumbel
theta.gumbel <- uniroot(
  function(theta) CorM1M2(theta, C.gumbel) - rho,
  interval = c(1, 5), tol = 1e-10
)$root

# comparaison des paramètres de dépendance
c(Frank = theta.frank,
  Clayton = theta.clayton,
  Gumbel = theta.gumbel)

# vérification des corrélations de Pearson
c(Frank = CorM1M2(theta.frank, C.frank),
  Clayton = CorM1M2(theta.clayton, Clayton),
  Gumbel = CorM1M2(theta.gumbel, Gumbel))

### Fonctions de masse de probabilité conjointes

Fm1m2.frank <- outer(Fm1, Fm2, C.frank, alpha = theta.frank)
Fm1m2.clayton <- outer(Fm1, Fm2, C.clayton, alpha = theta.clayton)
Fm1m2.gumbel <- outer(Fm1, Fm2, C.gumbel, alpha = theta.gumbel)

fm1m2.frank <- t(apply(apply(Fm1m2.frank, 2, diff), 1, diff))
fm1m2.clayton <- t(apply(apply(Fm1m2.clayton, 2, diff), 1, diff))
fm1m2.gumbel <- t(apply(apply(Fm1m2.gumbel, 2, diff), 1, diff))

### V.a. N

kmax <- max(m1) + max(m2)
kn <- 0:kmax

# Fonction de masse de probabilité de N
fn <- function(fm1m2) {
  sapply(kn, function(s)
    sum(fm1m2[outer(m1, m2, "+") == s]))
}

fn.frank <- fn(fm1m2.frank)
fn.clayton <- fn(fm1m2.clayton)
fn.gumbel <- fn(fm1m2.gumbel)

# Vérification
c(sum(fn.frank), sum(fn.clayton), sum(fn.gumbel))

kappa <- c(0.95, 0.995)

# Fonction VaR
VaR <- function(u, f)
  kn[min(which(cumsum(f) >= u))]

# Calcul de la VaR pour les trois copules
VaR.frank <- sapply(kappa, VaR, f = fn.frank)
VaR.clayton <- sapply(kappa, VaR, f = fn.clayton)
VaR.gumbel <- sapply(kappa, VaR, f = fn.gumbel)

# Fonction TVaR
TVaR <- function(u, f) {
  v <- VaR(u, f)
  v + sum(pmax(kn - v, 0)*f)/(1-u)
}

# Calcul de la TVaR pour les trois copules
TVaR.frank <- sapply(kappa, TVaR, f = fn.frank)
TVaR.clayton <- sapply(kappa, TVaR, f = fn.clayton)
TVaR.gumbel <- sapply(kappa, TVaR, f = fn.gumbel)

resultats <- rbind(
  theta = c(theta.frank, theta.clayton, theta.gumbel),
  "VaR(0.95)" = c(VaR.frank[1], VaR.clayton[1], VaR.gumbel[1]),
  "VaR(0.995)" = c(VaR.frank[2], VaR.clayton[2], VaR.gumbel[2]),
  "TVaR(0.95)" = c(TVaR.frank[1], TVaR.clayton[1], TVaR.gumbel[1]),
  "TVaR(0.995)" = c(TVaR.frank[2], TVaR.clayton[2], TVaR.gumbel[2])
)

colnames(resultats) <- c("Frank", "Clayton", "Gumbel")

round(resultats, 5)
