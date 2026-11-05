################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### École d'actuariat - Université Laval ######################################
#### Ateliers - Solutions aux questions manquées ###############################
#### Auteur : Philippe Leblanc #################################################
################################################################################

###
### Examen partiel informatique 2026 - Question 8
###

rm(list = ls())

n <- 2^4; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
q <- c(0.19, 0.26, 0.22, 0.37, 0.27, 0.39, 0.25, 0.35)

##
## a)
##

# i)

pI1I2_00 <- 0.57
pI1I2_10 <- (1 - pI1I2_00 - q[2])
pI1I2_01 <- (1 - pI1I2_00 - q[1])
pI1I2_11 <- (1 - pI1I2_00 - pI1I2_10 - pI1I2_01)

fgpI1I2 <- function(t1, t2) (
    pI1I2_00 + pI1I2_10 * t1 + pI1I2_01 * t2 + pI1I2_11 * t1 * t2
)

pI3I4_11 <- 0.15
pI3I4_01 <- (q[3] - pI3I4_11)
pI3I4_10 <- (q[4] - pI3I4_11)
pI3I4_00 <- (1 - q[3] - q[4] + pI3I4_11)

fgpI3I4 <- function(t3, t4) (
    pI3I4_00 + pI3I4_10 * t3 + pI3I4_01 * t4 + pI3I4_11 * t3 * t4
)

pI5I6_11 <- 0.05
pI5I6_01 <- (q[5] - pI5I6_11)
pI5I6_10 <- (q[6] - pI5I6_11)
pI5I6_00 <- (1 - q[5] - q[6] + pI5I6_11)

fgpI5I6 <- function(t5, t6) (
    pI5I6_00 + pI5I6_10 * t5 + pI5I6_01 * t6 + pI5I6_11 * t5 * t6
)

pI7I8_00 <- 0.55
pI7I8_10 <- (1 - pI7I8_00 - q[8])
pI7I8_01 <- (1 - pI7I8_00 - q[7])
pI7I8_11 <- (1 - pI7I8_00 - pI7I8_10 - pI7I8_01)

fgpI7I8 <- function(t7, t8) (
    pI7I8_00 + pI7I8_10 * t7 + pI7I8_01 * t8 + pI7I8_11 * t7 * t8
)

# ii)
fgpN <- function(t) (
    fgpI1I2(t, t) * fgpI3I4(t, t) * fgpI5I6(t, t) * fgpI7I8(t, t)
)

##
## b)
##

# i)
fn <- Re(fft(fgpN(e), TRUE))/n
fn[c(3, 5) + 1]

# ii)
SL <- function(d) sum(pmax(k - d, 0) * fn)
SL(5)


###
### Examen partiel traditionnel  2026 - Question 4
###

rm(list = ls())

lam <- c(12, 15); y12 <- 3

###
### d)
###

## i)
lam[1] - y12 + y12/lam[2] * 20

## ii)
p10 <- (lam[1] - y12)/(sum(lam) - y12)
p01 <- (lam[2] - y12)/(sum(lam) - y12)
p11  <- y12/(sum(lam) - y12)
p00 <- 0

# vérification
p00 + p10 + p01 + p11 # 1

# marginales
q1 <- p10 + p11
q2 <- p01 + p11

# vérifications
all.equal(p00, max((1 - q1) + (1 - q2) - 1, 0))         # F(0, 0)
all.equal(p00 + p10, max(1 + (1 - q2) - 1, 0))          # F(1, 0)
all.equal(p00 + p01, max((1 - q1) + 1 - 1, 0))          # F(0, 1)
all.equal(p00 + p10 + p01 + p11, max(1 + 1 - 1, 0))     # F(1, 1)

# vérifications équivalentes ici (F(0, 0) et F(1, 1) sont triviales)
all.equal(p00 + p10, 1 - q2) 
all.equal(p00 + p01, 1 - q1)

## iii)
(lam[1] - y12) + (lam[2] - y12) + y12


###
### Examen partiel traditionnel  2026 - Question 8
###

rm(list = ls())

w <- 0:4/4
i <- expand.grid(rep(list(0:1), 4))
km <- list(c(3, 4), c(5, 6), c(7, 8))

fI1234 <- c(
    0.3, 0, 0, 0, 0, 0.1, 0.1, 0, 0, 0.1, 0.1, 0, 0, 0, 0, 0.3
)
fI1256 <- c(
    0, 0, 0, 0.3, 0, 0.1, 0.1, 0, 0, 0.1, 0.1, 0, 0.3, 0, 0, 0
)
fI1278 <- c(
    0.09, 0.06, 0.06, 0.09, 0.06, 0.04, 0.04, 0.06,
    0.06, 0.04,0.04, 0.06, 0.09, 0.06, 0.06, 0.09
)
sum(fI1234); sum(fI1256); sum(fI1278)

###
### b)
###

fW12km <- function(f, w) sum(
    f[i[, 1] + i[, 2] + i[, 3] + i[, 4] == w * 4]
)

# Calcul des fmps
fW12 <- c(0.3, 0, 0.4, 0, 0.3)
fW1234 <- sapply(w, fW12km, f = fI1234)
fW1256 <- sapply(w, fW12km, f = fI1256)
fW1278 <- sapply(w, fW12km, f = fI1278)

# Calcul des variances
VarW12 <- sum((w - sum(w * fW12))^2 * fW12)
VarW1234 <- sum((w - sum(w * fW1234))^2 * fW1234)
VarW1256 <- sum((w - sum(w * fW1256))^2 * fW1256)
VarW1278 <- sum((w - sum(w * fW1278))^2 * fW1278)

###
### c)
###

km[which(c(VarW1234, VarW1256, VarW1278) == 0)]
VarW12 - VarW1256

###
### d)
###

km[which(VarW12 - c(VarW1234, VarW1256, VarW1278) == 0)]
