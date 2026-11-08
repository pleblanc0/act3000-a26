rm(list = ls())


xx <- seq(0, 100, 20)
ss <- seq(0, 200, 20)

kappa <- c(0.01, 0.1, 0.5, 0.9,
           0.99, 0.999, 0.9999)

#### Exemple no 1 ####
#### Serie 1 ####

alpha <- c(2.1, 3.4)
beta <- 1/10

mm <- seq(0, 3)
fm1m2 <- matrix(c(
  0.30, 0.10, 0.02, 0.03,
  0.20, 0.10, 0.01, 0.01,
  0.04, 0.02, 0.05, 0.03,
  0.01, 0.02, 0.01, 0.05
), nrow = 4, byrow = TRUE)

fm1 <- rowSums(fm1m2)
fm2 <- colSums(fm1m2)

#### Fonctions de répartition de X1 et X2 ####

FX1 <- function(x) fm1[1] + sum(fm1[-1] * pgamma(x, mm[-1] * alpha[1], beta))

FX2 <- function(x) fm2[1] + sum(fm2[-1] * pgamma(x, mm[-1] * alpha[2], beta))

round(cbind(xx,
            FX1 = sapply(xx, FX1),
            FX2 = sapply(xx, FX2)
), 4)

#### Espérances et variances ####

EspM <- c(sum(mm * fm1), sum(mm * fm2))
VarM <- c(sum(mm^2 * fm1), sum(mm^2 * fm2)) - EspM^2
EspB <- alpha/beta
VarB <- alpha/beta^2

EspX <- EspM * EspB
VarX <- EspM * VarB + VarM * EspB^2

round(cbind(EspM, VarM, EspB, VarB, EspX, VarX), 4)


#### Covariances ####

EspM12 <- sum(outer(mm, mm, "*") * fm1m2)

CovM1M2 <- EspM12 - prod(EspM)
CovX1X2 <- prod(EspB) * CovM1M2

CorM1M2 <- CovM1M2/sqrt(prod(VarM))
CorX1X2 <- CovX1X2/sqrt(prod(VarX))

round(cbind(EspM12, CovM1M2, CovX1X2, CorM1M2, CorX1X2), 4)

#### Espérance et variance de S ####

EspS <- sum(EspX)
VarS <- sum(VarX) + 2 * CovX1X2
round(cbind(EspS, VarS), 4)

## Fonction de repartition conjointe

FX1X2 <- function(x1, x2){
  F1 <- c(1, pgamma(x1, alpha[1] * mm[-1], beta))
  F2 <- c(1, pgamma(x2, alpha[2] * mm[-1], beta))
  
  sum(fm1m2 * outer(F1, F2,'*'))
}

Fx1x2 <- outer(xx, xx, Vectorize(FX1X2))
round(Fx1x2, 4)

## Fonction de repartition de S

aS <- outer(mm, mm, function(i, j)
  alpha[1] * i + alpha[2] * j)[-1]

fm1m2_0 <- fm1m2[1, 1]
fm1m2_k <- fm1m2[-1]

FS <- function(x) {
  fm1m2_0 + sum(fm1m2_k * pgamma(x, aS, beta))
}

round(cbind(ss, FS = sapply(ss, FS)), 4)

## prime stop loss de S

piS <- function(x) sum(
  fm1m2_k * (aS/beta * pgamma(x, aS + 1, beta, lower = FALSE) -
           x * pgamma(x, aS, beta, lower = FALSE))
)

sapply(ss, piS)

## VaR et TVaR

VaRS <- function(u) {
  if(u <= fm1m2[1, 1]) return(0)
  optimize(function(x) abs(FS(x) - u), c(0, 1000))$minimum
}

sapply(kappa, VaRS)

TVaRS <- function(u) VaRS(u) + piS(VaRS(u))/(1 - u)

TVaRS <- function(u) 1/(1 - u) * sum(
  fm1m2_k * aS/beta * pgamma(VaRS(u), aS + 1, beta, lower = FALSE)
)

sapply(kappa, TVaRS)

#### Exemple no 1 ####
#### Serie 2 ####

## Loi conditionnelle de S sachant M2 =2

pM2_2 <- fm2[3]

FS_M2_2 <- function(x) sum(
  fm1m2[, 3]/pM2_2 *
    pgamma(x, alpha[1] * mm + 2 * alpha[2], beta)
)

round(cbind(ss, FS_M2_2 = sapply(ss, FS_M2_2)), 4)

## Loi conditionnelle de S sachant M1+M2 = 3 

pN3 <- sapply(mm, function(i) fm1m2[i + 1, 4 - i])

FS_N3 <- function(x) sum(
  pN3 * pgamma(x, alpha[1] * mm + alpha[2] * (3 - mm), beta)
)/sum(pN3)

round(cbind(ss, FS_N3 = sapply(ss, FS_N3)), 4)

## Mesures entropiques de X1 et X2

rho <- c(0.001, 0.005, 0.01, 0.05)
phiX1 <- function(r) (1/r)*log(sum(fm1 *
                               (beta/(beta - r))^(alpha[1] * mm)))
phiX2 <- function(r) (1/r)*log(sum(fm2 *
                               (beta/(beta - r))^(alpha[2] * mm)))

round(cbind(rho,
            phiX1 = sapply(rho, phiX1),
            phiX2 = sapply(rho, phiX2)
), 4)

#### Mesure entropique de cX1 + (1-c)X2 ####

phiComb <- function(c, r) {
  M_comb <-  sum(fm1m2 * outer(
    (beta/(beta - r*c))^(alpha[1] * mm),
    (beta/(beta - r*(1-c)))^(alpha[2] * mm)
  ))
  (1/r)*log(M_comb)
}
cr <- expand.grid(c = c(1/4, 1/2, 3/4), rho = rho)

phiLin <- cr$c * sapply(cr$rho, phiX1) +
  (1 - cr$c) * sapply(cr$rho, phiX2)

phiMix <- mapply(phiComb, cr$c, cr$rho)

round(cbind(cr, phiLin, phiMix,
            ecart = phiLin - phiMix), 4)


#### Exemple no 2 ####

lam <- c(3, 2); a0 <- c(0, 0.5, 1, 1.5)
mu <- c(log(20) - 0.32, log(30) - 0.18)
sigma <- c(0.8, 0.6)
x <- seq(0, 100, 20); s <- seq(0, 200, 20)
kappa <- c(0.01, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999)
rho <- 0.01; d <- 120


## Esperances et variances

Espm <- lam
Varm <- lam
Espb <- exp(mu + sigma^2/2)
Varb <- (exp(sigma^2) - 1) * Espb^2
Espx <- Espm * Espb
Varx <- Espm * (Varb + Espb^2)
round(cbind(Espm, Varm, Espb, Varb, Espx, Varx), 4)

## Covariances

CovM1M2 <- a0
CovX1X2 <- prod(Espb) * CovM1M2
CorM1M2 <- CovM1M2/sqrt(prod(Varm))
CorX1X2 <- CovX1X2/sqrt(prod(Varx))
round(cbind(a0, CovM1M2, CovX1X2, CorM1M2, CorX1X2), 4)

## Esperance et variance de S

EspS <- sum(Espx)
VarS <- sum(Varx) + 2 * CovX1X2
round(cbind(a0, EspS, VarS), 4)

# Préparation de la simulation
m <- 1E6; set.seed(20260830)

X1.sim <- matrix(numeric(m * length(a0)), ncol = length(a0))
X2.sim <- matrix(numeric(m * length(a0)), ncol = length(a0))
S.sim <- matrix(numeric(m * length(a0)), ncol = length(a0))


for(k in 1:length(a0)){
  
  # Simulation des fréquences de Poisson
  U <- matrix(runif(m * 3), ncol = 3)
  
  K0 <- qpois(U[, 1], a0[k])
  K1 <- qpois(U[, 2], lam[1] - a0[k])
  K2 <- qpois(U[, 3], lam[2] - a0[k])
  
  M <- cbind(K1 + K0, K2 + K0)
  
  # Simulation des montants composés
  X <- matrix(numeric(m * 2), ncol = 2)
  
  for(j in 1:m){
    X[j, 1] <- sum(qlnorm(runif(M[j, 1]),
                          mu[1], sigma[1]))
    
    X[j, 2] <- sum(qlnorm(runif(M[j, 2]),
                          mu[2], sigma[2]))
  }
  
  X1.sim[, k] <- X[, 1]
  X2.sim[, k] <- X[, 2]
  S.sim[, k] <- rowSums(X)
}

colnames(S.sim) <- a0

head(S.sim)

i <- 3 # gamma0 = 1

#### Fonction de répartition conjointe

FX1X2 <- function(x1, x2) mean(
  X1.sim[, i] <= x1 & X2.sim[, i] <= x2
)

Fx1x2 <- outer(xx, xx, Vectorize(FX1X2))
round(Fx1x2, 4)

#### Fonction de répartition de S 

FS <- function(x) mean(S.sim[, i] <= x)

round(cbind(ss, FS = sapply(ss, FS)), 4)

### Prime stop-loss

piS <- function(x) mean(pmax(S.sim[, i] - x, 0))

round(cbind(ss, piS = sapply(ss, piS)), 4)


#### vaR et TVaR

VaRS <- function(u) sort(S.sim[, i])[u * nrow(S.sim)]

round(cbind(kappa, VaR = sapply(kappa, VaRS)), 4)


TVaRS <- function(u) VaRS(u) + piS(VaRS(u))/(1 - u)

round(cbind(
  kappa,
  VaR = sapply(kappa, VaRS),
  TVaR = sapply(kappa, TVaRS)
), 4)


####  Z = min(S, d) 

d <- 120

Z <- pmin(S.sim[, i], d)

## Espérance et variance

EspZ <- mean(Z)
VarZ <- var(Z)

round(cbind(EZ, VarZ), 4)


## Mesure entropique

rho <- c(0.001, 0.005, 0.01, 0.05)

phiZ <- function(r) log(mean(exp(r * Z)))/r

round(cbind(rho, phiZ = sapply(rho, phiZ)), 4)

#### Exemple no 3 ####

#### Discrétisation lower

fb1 <- plnorm((k + 1)*h, mu[1], sig[1]) -
  plnorm(k*h, mu[1], sig[1])

fb2 <- plnorm((k + 1)*h, mu[2], sig[2]) -
  plnorm(k*h, mu[2], sig[2])

sum(fb1)
sum(fb2)

#### Loi de B1 + B2 

t1 <- fft(fb1)
t2 <- fft(fb2)

fb12 <- Re(fft(t1 * t2, TRUE))/n

sum(fb12)

#### Loi de la somme S 

fS <- function(a0){
  
  lambdaN <- sum(lam) - a0
  
  fb <- ((lam[1] - a0)*fb1 +
           (lam[2] - a0)*fb2 +
           a0*fb12)/lambdaN
  
  fSt <- exp(lambdaN * (fft(fb) - 1))
  
  Re(fft(fSt, TRUE))/n
}

fs <- sapply(a0, fS_1)

round(colSums(fs1), 6)

#### Pour gamma0 = 1 ####

fs <- fs1[, 3]

## Prime stop-loss 

stoploss <- function(g) sum(pmax(k * h - g, 0) * fs)

sapply(ss, stoploss)


## VaR 

VaR <- function(u) k[min(which(cumsum(fs) >= u))] * h

sapply(u, VaR)


## TVaR 

TVaR <- function(u) sum(
  pmax(k * h - VaR(u), 0) * fs
)/(1 - u) + VaR(u)

sapply(kappa, TVaR)

#### Z = min(S, d)

d <- 120
z <- pmin(k * h, d)

## Espérance

EspZ <- sum(z * fs)

## Variance

VarZ <- sum(z^2 * fs) - EspZ^2

round(cbind(EspZ, VarZ), 4)


## Mesure entropique 

rho <- c(0.001, 0.005, 0.01, 0.05)

phiZ <- function(r) log(sum(exp(r * z) * fs))/r

round(cbind(rho, phiZ = sapply(rho, phiZ)), 4)

