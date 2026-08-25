################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Guide de survie - FFT #####################################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

#===============================================================================
#== Utilisation de la fonction FFT =============================================
#===============================================================================

# définir une longeur de vecteur (puissance de 2 pour maximiser l'efficacité)
n <- 2^5

# créer le vecteur (en allant jusqu'à n - 1, on obtient un vecteur de longeur n)
k <- 0:(n - 1)

# calculer la transformée discrète de Fourier (retourne un vecteur complexe)
dft.k <- fft(k)

# calculer l'inverse de la DFT (normaliser et conserver la partie réelle)
inv.k <- Re(fft(dft.k, TRUE)/n)

# comparaison entre x et l'inverse de la DFT de x (doit être égal!)
all.equal(k, inv.k)


#===============================================================================
#== Convolution entre deux vecteurs indépendants ===============================
#===============================================================================

# préparation de l'exemple
n <- 2^10; k <- 0:(n - 1)

# paramètres des lois discrètes
lam <- 2; r <- 2; q <- 1/2

# vecteurs de fmp des composantes de la somme
fx <- dpois(k, lam)
fy <- dnbinom(k, r, q)

# transformée de fourier discrète de S
dft.s <- fft(fx) * fft(fy)

# inversion pour obtenir la fmp de S
fs <- Re(fft(dft.s, TRUE))/n

# nettoyage numérique (optionnel)
fs <- zapsmall(fs)

# diagnostics (pour se vérifier!)
sum(fs) # == 1
min(fs) # >= 0
max(fs) # <= 1


#===============================================================================
#== Obtenir la fmp depuis une fgp ==============================================
#===============================================================================

# préparation de l'exemple
n <- 2^10; k <- 0:(n - 1)
e <- exp(-2i * pi * k/n)

# fmp de la loi de Poisson
lam <- 2
fx.direct <- dpois(k, lam)

# fonction génératrice des probabilités
fgp <- function(t) exp(lam * (t - 1))

# test diagnostic (P_X(0) = Pr(X = 0))
fgp(0) == fx.direct[1]

# obtenir la fmp
fx.fgp <- Re(fft(fgp(e), TRUE))/n

# comparaison
all.equal(fx.direct, fx.fgp)


#===============================================================================
#== Travailler avec des lois composées =========================================
#===============================================================================

# préparation de l'exemple
n <- 2^10; k <- 0:(n - 1)
e <- exp(-2i * pi * k/n)

# paramètres des lois discrètes
lam <- 2; r <- 2; q <- 1/2

# fgp de X ~ PoissonComposée(lam, Fb), B ~ BinomNegative(r, q)
fgpM <- function(t) exp(lam * (t - 1))
fgpB <- function(t) (q/(1 - (1 - q) * t))^r
fgpX <- function(t) fgpM(fgpB(t))

# fonction de masse de probabilité de X
dft.x <- fgpX(e)
fx <- Re(fft(dft.x, TRUE))/n

# test diagnostic
all.equal(fgpX(0), fx[1])


#===============================================================================
#== Lois bivariées (avec dépendance) ===========================================
#===============================================================================

# préparation de l'exemple
n <- 2^10; k <- 0:(n - 1)

# paramètres de la loi
a0 <- 1; lam <- c(2, 3)

# fgp bivariée d'une loi Poisson Teicher
fgpM <- function(t1, t2) exp(
    (lam[1] - a0) * (t1 - 1) + (lam[2] - a0) * (t2 - 1) + a0 * (t1 * t2 - 1)
)

# question: comment trouver la fgp de la somme (afin d'obtenir sa fmp)
# option 1: trouver les fgp marginales et faire le produit (très mauvaise idée!)
# option 2: trouver la fgp de la somme et ensuite faire la procédure habituelle

# puisqu'il y a de la dépendance, on prend l'option 2!
fgpN <- function(t) fgpM(t, t)

# fmp de la somme
dft.n <- fgpN(e)
fn <- Re(fft(dft.n, TRUE))/n

# test diagnostic
all.equal(exp(a0 - sum(lam)), fn[1])
