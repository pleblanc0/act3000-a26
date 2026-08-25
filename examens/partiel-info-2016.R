################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen partiel informatique 2016 ################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Question 1
###

rm(list = ls())

a <- 41358; m <- 2^31 - 1; x0 <- 20150418; n <- 10000
gam <- c(2, 0.5, 1); be <- c(1/100, 1/200)

x <- numeric(3 * n); x[1] <- (a * x0) %% m
for (i in 2:(3 * n)) x[i] <- (a * x[i - 1]) %% m

U <- matrix(x/m, ncol = 3, byrow = TRUE)
Y1 <- qgamma(U[, 1], gam[1], 1)
Y2 <- qgamma(U[, 2], gam[2], 1)
Y3 <- qgamma(U[, 3], gam[3], 1)

## c)
X1 <- (Y1 + Y3)/be[1]
X2 <- (Y2 + Y3)/be[2]

S <- X1 + X2
c(X1[3], X2[3]); c(X1[4], X2[4])

## d)
mean(I(X1 <= 500) * I(X2 <= 500))
mean(pmax(X1 - 500, 0) * I(X2 > 500))

VaR <- function(u) sort(S)[u * n]
VaR(0.99)

TVaR <- function(u) mean(S[S > VaR(u)])
TVaR(0.99)


###
### Question 2
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
phi <- function(t, a) exp(a * (t - 1))

fgp <- function(t1, t2, t3) (
    phi(t1, 0.5) * phi(t2, 1.2) * phi(t3, 0.7) *
    phi(t1 * t2, 0.3) * phi(t2 * t3, 0.4) * phi(t1 * t2 * t3, 0.2)
)

## c)
fn <- Re(fft(sapply(e, function(t) fgp(t, t, t))))/n
fn[0:2 + 1]

## d)
fm1 <- Re(fft(sapply(e, function(t) fgp(t, 1, 1))))/n
fm2 <- Re(fft(sapply(e, function(t) fgp(1, t, 1))))/n
fm3 <- Re(fft(sapply(e, function(t) fgp(1, 1, t))))/n
c(fm1[5 + 1], fm2[5 + 1], fm3[5 + 1])

## e)
VaR <- function(u, f) k[min(which(cumsum(f) >= u))]
TVaR <- function(u, f) VaR(u, f) + sum(pmax(k - VaR(u, f), 0) * f)/(1 - u)

BM <- function(u) sum(sapply(list(fm1, fm2, fm3), TVaR, u = u)) - TVaR(u, fn)
BM(0.99)


###
### Question 3
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); h <- 1
r <- 1.5; q <- 1/3; sig <- 0.6; mu <- log(10) - sig^2/2

## a)
fb <- c(0, diff(plnorm(k, mu, sig)))
fb[c(0, 10) + 1]

## b)
fx <- Re(fft((q/(1 - (1 - q) * fft(fb)))^r, TRUE))/n
fx[c(0, 20) + 1]

## c)
Fx <- cumsum(fx)
Fx[50 + 1]

SL <- function(d) sum(pmax(k - d, 0) * fx)
SL(50)

## e)
50 + SL(50)/(1 - Fx[50 + 1])


###
### Question 4
###

rm(list = ls())

m <- 1E5; set.seed(19661122); U <- runif(m)
al <- c(3, 2.5); lam <- c(200, 150)
qpareto <- function(u, a, l) l * ((1 - u)^(-1/a) - 1)

## a)
X1anti <- qpareto(U, al[1], lam[1])
X2anti <- qpareto(1 - U, al[2], lam[2])

Santi <- X1anti + X2anti
Santi[2]

Eanti <- 2 * mean(Santi * I(Santi < 100) + Santi * I(Santi > 1000))
Eanti

## b)
X1como <- qpareto(U, al[1], lam[1])
X2como <- qpareto(U, al[2], lam[2])

Scomo <- X1como + X2como
Scomo[2]

Ecomo <- 2 * mean(Scomo * I(Scomo < 100) + Scomo * I(Scomo > 1000))
Ecomo

## d)
c(var(Santi), var(Scomo))
c(cor(X1anti, X2anti), cor(X1como, X2como))


###
### Question 5
###

rm(list = ls())

n <- 2^12; k <- 0:(n - 1); e <- exp(2i * pi * k/n)
r <- 2; lam <- c(1.5, 2); gam <- 0.8; m <- 10; q <- c(0.2, 0.3)

fgpB1 <- function(t) (1 - q[1] + q[1] * t)^m
fgpB2 <- function(t) (1 - q[2] + q[2] * t)^m

fgmT <- function(t1, t2) (
    ((1 - gam)/((1 - (1 - gam)/r * t1) * (1 - (1 - gam)/r * t2) - gam))^r
)
fgpX <- function(t1, t2) fgmT(
    lam[1] * (fgpB1(t1) - 1), lam[2] * (fgpB2(t2) - 1)
)

# i)
fgpS <- function(t) fgpX(t, t)
fs <- Re(fft(fgpS(e)))/n

# ii)
Fs <- cumsum(fs)
Fs[c(30, 40) + 1]
