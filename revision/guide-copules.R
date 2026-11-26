################################################################################
#### Théorie du risque - Automne 2026 ##########################################
#### Guide de survie - copules ################################################
#### Auteur: Philippe Leblanc ##################################################
################################################################################

#===============================================================================
#== Modélisation avec des copules  =============================================
#===============================================================================

# exemple avec copule clayton
al <- 2; be <- 0.1

# marginales exponentielles
Fx <- function(x) pexp(x, be)

# expression de la copule
C.clayton <- function(u, v) (u^(-al) + v^(-al) - 1)^(-1/al)

# fonction de répartition conjointe
Fx12 <- function(x1, x2) C.clayton(Fx(x1), Fx(x2))

# méthode des rectangles
rectangles <- function(s, m, F)
{
    n <- 2^m; h <- s/n
    masse <- function(i, y) F(i * h, y) - F((i - 1) * h, y)

    lower <- sum(sapply(seq_len(n - 1), function(i) masse(i, (n - i) * h)))
    upper <- sum(sapply(seq_len(n), function(i) masse(i, (n + 1 - i) * h)))

    c(lower = lower, upper = upper)
}
rectangles(10, 10, Fx12)

###
### à poursuivre...
###