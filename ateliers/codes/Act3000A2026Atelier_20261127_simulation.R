################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### École d'actuariat - Université Laval ######################################
#### Ateliers - Simulation Copule #######################################################
#### Auteur : Omar Essakine #################################################
################################################################################

###
### Simulation copule FGM
###

nsim <- 1000
theta <- 0.1

# on produit un échantillon de U1
U1 <- runif(nsim)

# on produit un échantillon de V2
V2 <- runif(nsim)

# on simule W1 à partir de U1
W1 <- theta * (2 * U1 - 1) - 1

# on simule W2 à partir de U1 et V2
W2 <- (1 - theta * (2 * U1 - 1))^2 +
  4 * theta * V2 * (2 * U1 - 1)

# on calcule U2 à partir de W1, W2 et V2
U2 <- 2 * V2 / (sqrt(W2) - W1)

# matrice de simulations
vU <- cbind(U1, U2)

plot(U1, U2,
     main = expression("Copule FGM (" * theta == 0.1 * ")"),
     xlab = expression(U[1]),
     ylab = expression(U[2]))

###
### Simulation copule Clayton 
###

nsim <- 1000
alpha <- c(0.5,3.5)
j <- c(1,2,nsim)

# on produit trois échantillons de (U1,U2), un pour chaque valeur de alpha
vU <- matrix(0, nrow = nsim, ncol = 2*length(alpha))

# On utilise les mêmes observations uniformes pour simuler
vV <- matrix(runif(nsim * 3), nsim, 3, byrow = TRUE)

for (k in seq_along(alpha)) {
    
  # on simule Theta par le loi gamma
    Theta <- qgamma(vV[,1],1/alpha[k],1)
    
  # On simule U1 et U2
    for (i in 1:2) {
      # observation de V_i
      V <- vV[, i + 1]
      
      # méthode inverse de simulation pour la distribution exponentielle
      Y <- qexp(V, rate = Theta)
      
      # on utilise le Corrolaire 2
      U <- (1 + Y)^(-1 / alpha[k])
      
      # On simule les observations pour la k-ième valeur de alpha.
      # Chaque valeur de alpha occupe 2 colonnes, donc 2*(k-1) colonnes
      # sont déjà remplies. On ajoute i pour obtenir l'indice de la colonne à remplir.
      vU[, 2 * (k - 1) + i] <- U
    }
}

# Vérification de la simulation pour les trois valeurs de alpha
data.frame(
  j = j,
  U1_alpha1 = vU[j, 1], U2_alpha1 = vU[j, 2],
  U1_alpha2 = vU[j, 3], U2_alpha2 = vU[j, 4])

# Graphiques pour visualiser l'impact de la dépendance sur les échantillons
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

plot(vU[, 1], vU[, 2], main = expression(alpha == 0.5),
     xlab = 'U1', ylab = 'U2', xlim = c(0, 1), ylim = c(0, 1))

plot(vU[, 3], vU[, 4], main = expression(alpha == 3.5),
     xlab = 'U1', ylab = 'U2', xlim = c(0, 1), ylim = c(0, 1))


###
### Simulation copule Clayton (méthode copule conditionnelle)
###

# on produit un seul échantillon de U1
U1 <- runif(nsim)

# on produit un seul échantillon de V2
V2 <- runif(nsim)

# on produit deux échantillons de U2, un pour chaque valeur de alpha
U2 <- sapply(alpha, function(a) {
  ((V2 * U1^(a + 1))^(-a / (a + 1)) -
     U1^(-a) + 1)^(-1 / a)
})


###
### Simulation copule Frank (methode conditionelle)
###

nsim <- 1000
alpha <- c(2,7)
j <- c(1, 2, nsim)

# on produit un seul échantillon de U1
U1 <- runif(nsim)

# on produit un seul échantillon de V2
V2 <- runif(nsim)

# on produit deux échantillons de U2, un pour chaque valeur de alpha
U2 <- matrix(0, nrow = nsim, ncol = length(alpha))

# on utilise les mêmes observations uniformes pour simuler
for (k in seq_along(alpha)) {
  
  # on calcule U2 à partir de U1 et V2
  U2[, k] <- -1 / alpha[k] * log(1 +
                                   (V2 * (exp(-alpha[k]) - 1)) /
                                   (exp(-alpha[k] * U1) * (1 - V2) + V2))
}

# Vérification de la simulation pour les deux valeurs de alpha
data.frame(
  j = j,
  U1 = U1[j],
  U2_alpha1 = U2[j, 1],
  U2_alpha2 = U2[j, 2]
)

# Graphiques pour visualiser l'impact de la dépendance sur les échantillons
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

plot(U1, U2[, 1], main = expression(alpha == 0.5),
     xlab = 'U1', ylab = 'U2', xlim = c(0, 1), ylim = c(0, 1))

plot(U1, U2[, 2], main = expression(alpha == 3.5),
     xlab = 'U1', ylab = 'U2', xlim = c(0, 1), ylim = c(0, 1))

###
### Simulation copule AMH (méthode copule archimédienne)
###

nsim <- 1000
alpha <- c(0.5, 0.9)
j <- c(1, 2, nsim)

# on produit deux échantillons de (U1,U2), un pour chaque valeur de alpha
vU <- matrix(0, nrow = nsim, ncol = 2*length(alpha))

# on utilise les mêmes observations uniformes pour simuler
vV <- matrix(runif(nsim * 3), nsim, 3, byrow = TRUE)

for (k in seq_along(alpha)) {
  
  # on simule Theta par la loi géométrique de type 2
  Theta <- 1 + qgeom(vV[, 1], prob = 1 - alpha[k])
  
  # on simule U1 et U2
  for (i in 1:2) {
    
    # observation de V_i
    V <- vV[, i + 1]
    
    # méthode inverse de simulation pour la distribution exponentielle
    Y <- qexp(V, rate = Theta)
    
    # on utilise le Corollaire 2
    U <- (1 - alpha[k]) / (exp(Y) - alpha[k])
    
    # on simule les observations pour la k-ième valeur de alpha.
    # Chaque valeur de alpha occupe 2 colonnes, donc 2*(k-1) colonnes
    # sont déjà remplies. On ajoute i pour obtenir l'indice de la colonne à remplir.
    vU[, 2 * (k - 1) + i] <- U
  }
}

# Vérification de la simulation pour les deux valeurs de alpha
data.frame(
  j = j,
  U1_alpha1 = vU[j, 1], U2_alpha1 = vU[j, 2],
  U1_alpha2 = vU[j, 3], U2_alpha2 = vU[j, 4]
)

# Graphiques pour visualiser l'impact de la dépendance sur les échantillons
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

plot(vU[, 1], vU[, 2], main = expression(alpha == 0.5),
     xlab = 'U1', ylab = 'U2', xlim = c(0, 1), ylim = c(0, 1))

plot(vU[, 3], vU[, 4], main = expression(alpha == 0.9),
     xlab = 'U1', ylab = 'U2', xlim = c(0, 1), ylim = c(0, 1))


###
### Simulation copule Gumbel (méthode copule archimédienne) (pas au cours)
###

# Simulation de la loi stable positive
rstable <- function(n, phi) { 
  
  # Document de reference obligatoire (page 329)
  R <- runif(n, min = -pi / 2, max = pi / 2)
  V <- rexp(n, 1)
  
  A <- sin(phi * (R + pi / 2)) / cos(R)^(1 / phi)
  B <- (cos(phi * pi / 2 + (phi - 1) * R) / V)^((1 - phi) / phi)
  
  A * B
}

nsim <- 1000
alpha <- c(1.25, 2.5)
j <- c(1, 2, nsim)

# on produit deux échantillons de (U1,U2), un pour chaque valeur de alpha
vU <- matrix(0, nrow = nsim, ncol = 2*length(alpha))

# on utilise les mêmes observations uniformes pour simuler
vV <- matrix(runif(nsim * 2), nsim, 2, byrow = TRUE)

for (k in seq_along(alpha)) {
  
  # on simule Theta par la loi stable positive
  Theta <- rstable(nsim, phi = 1 / alpha[k])
  
  # on simule U1 et U2
  for (i in 1:2) {
    
    # observation de V_i
    V <- vV[, i]
    
    # méthode inverse de simulation pour la distribution exponentielle
    Y <- qexp(V, rate = Theta)
    
    # on utilise le Corollaire 2
    U <- exp(-Y^(1 / alpha[k]))
    
    # on simule les observations pour la k-ième valeur de alpha.
    # Chaque valeur de alpha occupe 2 colonnes, donc 2*(k-1) colonnes
    # sont déjà remplies. On ajoute i pour obtenir l'indice de la colonne à remplir.
    vU[, 2 * (k - 1) + i] <- U
  }
}

# Graphiques pour visualiser l'impact de la dépendance sur les échantillons
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

plot(vU[, 1], vU[, 2], main = expression(alpha == 1.5),
     xlab = 'U1', ylab = 'U2', xlim = c(0, 1), ylim = c(0, 1))

plot(vU[, 3], vU[, 4], main = expression(alpha == 3.5),
     xlab = 'U1', ylab = 'U2', xlim = c(0, 1), ylim = c(0, 1))

