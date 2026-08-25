################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Solutions examen final informatique 2017 ##################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

###
### Partie 1: Question 1
###

rm(list = ls())

n <- 2^14; k <- 0:(n - 1); e <- exp(-2i * pi * k/n)
gam <- c(10, 25, 15, 6, 8)

## b)
fgp <- function(t1, t2, t3) exp(
    gam[1] * (t1 - 1) + gam[2] * (t2 - 1) + gam[3] * (t3 - 1) +
    gam[4] * (t2 * t3 - 1) + gam[5] * (t1 * t2 * t3 - 1)
)

## g) à refaire avec une poisson composée
fs <- Re(fft(sapply(e, function(t) fgp(t, t, t)), TRUE))/n

# i)
Fs <- cumsum(fs)
Fs[c(70, 80, 90) + 1]

# ii)
VaR <- function(u) k[min(which(Fs >= u))]
VaR(0.99)


###
### Partie 2: Question 3        À FAIRE!!
###

rm(list = ls())
