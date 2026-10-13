################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### École d'actuariat - Université Laval ######################################
#### Ateliers - Préparation aux examens ########################################
#### Auteur : Philippe Leblanc #################################################
################################################################################

###
### Partiel informatique 2012 Q1
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); h <- 1
lam <- 0.64; mu <- 2.08; sig <- 0.83

## a)
fb <- c(0, diff(plnorm(k, mu, sig)))
fs <- Re(fft(exp(lam * (fft(fb) - 1)), TRUE))/n
fs[c(0, 50) + 1]

## b)
Fs <- cumsum(fs)
1 - Fs[50 + 1]

## c)
VaR <- function(u) k[min(which(Fs >= u))]
VaR(0.99)


###
### Partiel informatique 2012 Q4
###

rm(list = ls())

m <- 5; q <- c(0.4, 0.2); s <- 0:(2 * m); t <- -m:m

## a)
Fx12 <- function(k1, k2) min(pbinom(k1, m, q[1]), pbinom(k2, m, q[2]))
fx12 <- function(k1, k2) (
    Fx12(k1, k2) - Fx12(k1 - 1, k2) - Fx12(k1, k2 - 1) + Fx12(k1 - 1, k2 - 1)
)

# S = X1 + X2 -> la condition est {X1 = j, X2 = k - j}
fs <- sapply(s, function(k) sum(sapply(0:k, function(j) fx12(j, k - j))))

# T = X1 - X2 -> la condition est {X1 = j + k, X2 = j}
ft <- sapply(t, function(k) sum(sapply(0:m, function(j) fx12(j + k, j))))

## b)
Fx12 <- function(k1, k2) max(pbinom(k1, m, q[1]) + pbinom(k2, m, q[2]) - 1, 0)
fx12 <- function(k1, k2) (
    Fx12(k1, k2) - Fx12(k1 - 1, k2) - Fx12(k1, k2 - 1) + Fx12(k1 - 1, k2 - 1)
)

fs <- sapply(s, function(k) sum(sapply(0:k, function(j) fx12(j, k - j))))
ft <- sapply(t, function(k) sum(sapply(0:m, function(j) fx12(j + k, j))))


###
### Partiel traditionnel 2012 Q1
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); h <- 1000; kh <- k * h
r <- 0.5; q <- 3/4; al <- 3; lam <- 4000

## a)
Fc <- function(x) 1 - (lam/(lam + kh))^al
fc <- c(0, diff(Fc(kh)))

fb <- c(fc[0:4 + 1], sum(fc[-(0:4 + 1)]))
fb[c(1, 5) + 1]

## c)
fbt <- c(fb, rep(0, n - length(fb)))
fx <- Re(fft((q/(1 - (1 - q) * fft(fbt)))^r, TRUE))/n
fx[c(0, 4, 5) + 1]

## d)
VaR <- function(u) kh[min(which(cumsum(fx) >= u))]
VaR(0.99)


###
### Partiel traditionnel 2012 Q2
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); j <- 1:3
lam <- 1; al <- c(1/3, 1/2, 1/6); theta <- 2^(j - 2)
fb <- c(0, 0.269, 0.178, 0.123, 0.088, 0.065) # somme != 1

panjer.poisson <- function(fb, lam, kmax)
{
    fx <- exp(lam * (fb[1] - 1))
    fb <- c(fb, rep(0, kmax - length(fb) + 1))

    for (i in 1:kmax) fx <- c(
        fx, sum(fb[2:(i + 1)] * fx[i:1] * (lam * (1:i)/i))
    )
    fx
}

# a) X|Theta = theta_j  ~ PoisComp (lambda * theta_j, Fb)

# b)
fxCond <- sapply(theta, function(t) panjer.poisson(fb, lam * t, 5))
fx <- fxCond %*% al
fx[c(0, 5) + 1]


###
### Partiel traditionnel 2012 Q4
###

rm(list = ls())

al <- 2.5; lam <- 1500; mu <- 5; sig <- 1; be <- 1/2000

TVaRX1 <- function(u) lam * (al/(al - 1) * (1 - u)^(-1/al) - 1)
TVaRX2 <- function(u) exp(mu + sig^2/2) * (1 - pnorm(qnorm(u) - sig))/(1 - u)
TVaRX3 <- function(u) qexp(u, be) + 1/be

## a)
TVaRS <- function(u) TVaRX1(u) + TVaRX2(u) + TVaRX3(u)
TVaRS(0.99)


###
### Final traditionnel 2012 Q1
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); i <- 1:5
lam <- 0.06 - 0.01 * i; lamN <- sum(lam)

## b)
fb <- c(0, lam/lamN, rep(0, n - length(lam) - 1))

## c)
fs <- Re(fft(exp(lamN * (fft(fb) - 1)), TRUE))/n


###
### Partiel informatique 2015 Q1
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); h <- 1
lamPoi <- 1; r <- 2; q <- 2/3
be <- 1/10; al <- 6; lamPar <- 50

## a)
fb <- c(0, diff(pexp(k, be)))
fb[10 + 1]

fc <- c(0, diff(1 - (lamPar/(lamPar + k * h))^al))
fc[10 + 1]

## b)
fx <- Re(fft(exp(lamPoi * (fft(fb) - 1)), TRUE))/n
Fx <- cumsum(fx)
Fx[c(0, 10, 20, 80) + 1]

## c)
fy <- Re(fft((q/(1 - (1 - q) * fft(fc)))^r, TRUE))/n
Fy <- cumsum(fy)
Fy[c(0, 10, 20, 80) + 1]

## d)
fs <- Re(fft(fft(fx) * fft(fy), TRUE))/n
Fs <- cumsum(fs)
Fs[c(0, 10, 20, 80) + 1]

## e)
VaR <- function(u, F) k[min(which(F >= u))]
sapply(list(Fx, Fy, Fs), VaR, u = 0.95)


###
### Partiel informatique 2015 Q4
###

rm(list = ls())

n <- 2^10; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
p00 <- 0.6; p01 <- p10 <- 0.10; p11 <- 0.20; m <- 20

## a)
fgpM1M2 <- function(t1, t2) (p00 + p10 * t1 + p01 * t2 + p11 * t1 * t2)^m
fgpN <- function(t) fgpM1M2(t^2, t^3)

## c)
fk <- c(0.6, 0, 0.1, 0.1, 0, 0.2)
fk <- c(fk, rep(0, n - length(fk)))

## d)
fn <- Re(fft(fgpN(e), TRUE))/n
fn[c(0, 5, 20) + 1]

# en utilisant fk
fn.alt <- Re(fft(fft(fk)^m, TRUE))/n
all.equal(fn, fn.alt)


###
### Partiel traditionnel 2015 Q1
###

rm(list = ls())

n <- 2^6; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)

fgp <- function(t1, t2) (
    0.6 + 0.1 * t2^5 + 0.05 * t1 * t2^4 + 0.05 * t1^2 * t2^3 +
    0.05 * t1^3 * t2^2 + 0.05 * t1^4 * t2 + 0.1 * t1^5
)

## a)
fm1m2 <- Re(fft(outer(e, e, fgp), TRUE))/n^2

## b)
fm1 <- Re(fft(sapply(e, function(t) fgp(t, 1)), TRUE))/n
fm2 <- Re(fft(sapply(e, function(t) fgp(1, t)), TRUE))/n

## c)
CovM1M2 <- sum(outer(k, k) * fm1m2) - sum(k * fm1) * sum(k * fm2)
CovM1M2

## d)
fn <- Re(fft(sapply(e, function(t) fgp(t, t)), TRUE))/n


###
### Partiel informatique 2016 Q3
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
### Partiel informatique 2016 Q4
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
### Final informatique 2016 Q5
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
lam <- c(-log(0.01)/3, -log(0.01)/2); lamN <- sum(lam); q <- c(0.2, 0.4)

## d)
fgpC <- function(t) sum(lam/lamN * q * t/(1 - (1 - q) * t))
fc <- Re(fft(sapply(e, fgpC), TRUE))/n

## f)
fgpS <- function(t) exp(4 * lamN * (fgpC(t) - 1))
fs <- Re(fft(sapply(e, fgpS), TRUE))/n

Fs <- cumsum(fs)
1 - Fs[105 + 1]

## g)
VaR <- function(u) k[min(which(Fs >= u))]
VaR(0.9)


###
### Final traditionnel 2016 Q1
###

rm(list = ls())

q <- c(0.2, 0.1); be <- c(1, 1/2); U <- 0.85

Fxi <- function(x, i) 1 - q[i] + q[i] * pexp(be[i] * x)
VaRi <- function(u, i) ifelse(
    u < 1 - q[i], 0, qexp((u - (1 - q[i]))/q[i], be[i])
)

## a)
VaRi(U, 1)

## b)
VaRi(U, 2)

## c)
VaRi(U, 1) + VaRi(U, 2)

## d)
VaRi(U, 1) + VaRi(1 - U, 2)
