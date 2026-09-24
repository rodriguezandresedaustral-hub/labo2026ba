# ===========================================================================
# cazatalentos_simulacion.R
#
# Desafio CazaTalentos: ordenar a las candidatas segun cuantos encestes van a
# hacer en una ronda NUEVA de 100 tiros.
#
# LA IDEA
# -------
# El numero que trae cada cazatalentos (80 o 79) no alcanza, porque depende
# de COMO eligio a su candidata. Si eligio a la mejor de muchas jugadoras
# parecidas, la ganadora probablemente tuvo suerte ese dia (la maldicion del
# ganador del notebook, secciones 2.4 y 2.5) y en tiros nuevos va a meter
# menos.
#
# Para medir cuanto menos, se repite la historia de cada cazatalentos
# millones de veces:
#
#   1. Se inventa un pueblo: se sortea su nivel y el indice de enceste de
#      cada una de sus jugadoras.
#   2. Todas tiran sus 100 tiros, como en la historia.
#   3. Se mira si el resultado es EXACTAMENTE el que conto la cazatalentos
#      (por ejemplo, para C1: una con 80, dos con 79, dos con 78 y el resto
#      con menos). Si no lo es, ese pueblo se descarta.
#   4. En los pueblos que quedan, la candidata tira 100 tiros nuevos y se
#      anota cuanto mete.
#
# El promedio de lo anotado en el paso 4 es cuanto se espera que meta la
# candidata (ley de los grandes numeros: con muchas repeticiones el
# promedio se acerca al valor esperado). El filtro del paso 3 es lo que hace
# que ese promedio sea el de "una candidata elegida como la eligio esa
# cazatalentos", y no el de una jugadora cualquiera.
#
# EL SUPUESTO
# -----------
# El paso 1 necesita decidir como se inventa un pueblo, y eso el enunciado no
# lo dice. Se usa esto:
#   - cada pueblo tiene un nivel (el indice promedio de sus jugadoras),
#     sorteado entre 0,40 y 0,80
#   - dentro del pueblo, las jugadoras se reparten alrededor de ese nivel con
#     un desvio sorteado entre 0,005 (todas casi iguales) y 0,10 (muy
#     distintas)
#   - ninguna jugadora llega a 0,85, porque el enunciado dice que Diana
#     Taurasi (0,85) es superior a todas
# Otras formas razonables de inventar los pueblos mueven un poco los numeros
# pero no el orden por grupos.
#
# EL CRITERIO DEL ENUNCIADO
# -------------------------
#   Ci < Cj  si  P( encestes de Ci < encestes de Cj ) > 0,5  en 100 tiros
# Al final se estima esa probabilidad para cada par, tambien simulando.
#
# COMO CORRERLO
# -------------
#   Rscript cazatalentos_simulacion.R
# C2 es lo mas lento: que un pueblo de 200 termine EXACTAMENTE con una 80,
# seis 79 y cinco 78 es muy raro (unos 3 de cada millon), y hacen falta unos
# 50 millones de pueblos. En una maquina Linux de 8 nucleos tarda unos 5
# minutos; en Windows corre en un solo nucleo y tarda mucho mas. Para una
# prueba rapida, bajar PUEBLOS$C2.
# ===========================================================================

require("parallel")

RNGkind("L'Ecuyer-CMRG")          # generador que da resultados reproducibles aun repartiendo en varios nucleos
set.seed(100213)                  # mi semilla primigenia

TIROS <- 100                      # tiros por ronda, como en el enunciado
TOPE  <- 0.85                     # ninguna jugadora llega al indice de Taurasi

# rango del nivel de los pueblos y de que tan parejas son sus jugadoras
NIVEL_MIN  <- 0.40 ; NIVEL_MAX  <- 0.80
DESVIO_MIN <- 0.005; DESVIO_MAX <- 0.10

# cuantos pueblos inventar para cada cazatalentos
PUEBLOS <- list(C1 = 1e6, C2 = 5e7, C4 = 4e6, C5 = 2e7, C6 = 3e6)
LOTE    <- 1e5                    # se simulan de a 100.000 pueblos para no llenar la memoria

# en Linux se reparte el trabajo en todos los nucleos; en Windows, uno solo
NUCLEOS <- if (.Platform$OS.type == "windows") 1 else detectCores()


# ---------------------------------------------------------------------------
# PASO 1: inventar pueblos
#   devuelve una matriz: una fila por pueblo, una columna por jugadora, con
#   el indice de enceste de cada una
# ---------------------------------------------------------------------------
inventar_pueblos <- function(n_pueblos, n_jugadoras) {
  nivel  <- runif(n_pueblos, NIVEL_MIN, NIVEL_MAX)
  desvio <- runif(n_pueblos, DESVIO_MIN, DESVIO_MAX)

  # distribucion Beta con ese promedio y ese desvio (es la forma natural de
  # sortear un numero entre 0 y 1 con un promedio dado)
  k <- nivel * (1 - nivel) / desvio^2 - 1
  ok <- k > 0                                   # combinaciones imposibles se descartan
  nivel <- nivel[ok]; k <- k[ok]

  indices <- matrix(rbeta(length(nivel) * n_jugadoras,
                          rep(nivel * k, n_jugadoras),
                          rep((1 - nivel) * k, n_jugadoras)),
                    ncol = n_jugadoras)

  # si alguna jugadora supera a Taurasi, ese pueblo no puede existir
  indices[rowSums(indices >= TOPE) == 0, , drop = FALSE]
}

# PASO 2: cada jugadora tira sus 100 tiros (misma idea que ftirar() del notebook)
tirar <- function(indices) {
  matrix(rbinom(length(indices), TIROS, indices), ncol = ncol(indices))
}


# ---------------------------------------------------------------------------
# C1 y C2: la mejor de muchas
#   se simula el pueblo entero y se queda solo con los que terminaron
#   exactamente como lo conto la cazatalentos
# ---------------------------------------------------------------------------
simular_mejor_de_muchas <- function(n_jugadoras, con_79, con_78, n_pueblos) {
  un_lote <- function(l) {
    indices <- inventar_pueblos(LOTE, n_jugadoras)
    aciertos <- tirar(indices)

    # PASO 3: el filtro. Una sola con 80, nadie arriba de 80, y la cantidad
    # exacta de 79 y de 78 que conto la cazatalentos
    igual_a_la_historia <-
      rowSums(aciertos == 80) == 1 &
      rowSums(aciertos >  80) == 0 &
      rowSums(aciertos == 79) == con_79 &
      rowSums(aciertos == 78) == con_78

    # el indice de la ganadora (la que hizo 80) en cada pueblo que quedo
    filas <- which(igual_a_la_historia)
    if (length(filas) == 0) return(numeric(0))
    columna_ganadora <- apply(aciertos[filas, , drop = FALSE], 1, function(fila) which(fila == 80))
    indices[cbind(filas, columna_ganadora)]
  }
  unlist(mclapply(seq_len(n_pueblos / LOTE), un_lote, mc.cores = NUCLEOS))
}


# ---------------------------------------------------------------------------
# C4: la mejor de solo dos, una hizo 80 y la otra 75
# ---------------------------------------------------------------------------
simular_C4 <- function(n_pueblos) {
  indices  <- inventar_pueblos(n_pueblos, 2)
  aciertos <- tirar(indices)

  # PASO 3: una hizo 80 y la otra 75 (en cualquier orden)
  primera_80 <- aciertos[, 1] == 80 & aciertos[, 2] == 75
  segunda_80 <- aciertos[, 1] == 75 & aciertos[, 2] == 80

  c(indices[primera_80, 1], indices[segunda_80, 2])      # el indice de la que hizo 80
}


# ---------------------------------------------------------------------------
# C5: la sobrina, 10 rondas de 100
#   la entrenadora conoce las 10 rondas. Como cada tiro es independiente, lo
#   unico que importa para saber que tan buena es la sobrina es el total:
#   701 encestes en 1.000 tiros. El filtro se queda con las jugadoras que
#   hicieron exactamente 701 en 1.000.
# ---------------------------------------------------------------------------
simular_C5 <- function(n_jugadoras) {
  indices <- inventar_pueblos(n_jugadoras, 1)[, 1]      # una jugadora sola por pueblo
  total <- rbinom(length(indices), 1000, indices)
  indices[total == 701]
}


# ---------------------------------------------------------------------------
# C6 y C7: una jugadora que NO fue elegida por ganarle a otras
#   C6 la eligio por la quiniela antes de que tirara; C7 testeo sola a la
#   unica del paraje. En los dos casos el 79 no paso por ningun filtro de
#   "la mejor", asi que se simula una jugadora sola que haya hecho 79.
#   C7 elige presentarla despues de ver su 79, pero esa decision no cambia lo
#   que dice el 79 sobre ella: es la misma simulacion que C6.
# ---------------------------------------------------------------------------
simular_una_sola <- function(n_jugadoras, aciertos_observados) {
  indices <- inventar_pueblos(n_jugadoras, 1)[, 1]
  aciertos <- rbinom(length(indices), TIROS, indices)
  indices[aciertos == aciertos_observados]
}


# ---------------------------------------------------------------------------
# CORRER TODO
# ---------------------------------------------------------------------------
t0 <- Sys.time()
indice_candidata <- list(
  C1 = simular_mejor_de_muchas(100, 2, 2, PUEBLOS$C1),
  C2 = simular_mejor_de_muchas(200, 6, 5, PUEBLOS$C2),
  C4 = simular_C4(PUEBLOS$C4),
  C5 = simular_C5(PUEBLOS$C5),
  C6 = simular_una_sola(PUEBLOS$C6, 79)
)
indice_candidata$C7 <- indice_candidata$C6        # C7 es exactamente el mismo caso que C6

# si para alguna candidata no quedo ningun pueblo, hacen falta mas pueblos
vacias <- names(indice_candidata)[sapply(indice_candidata, length) == 0]
if (length(vacias) > 0) stop("No quedo ningun pueblo para ", paste(vacias, collapse = ", "),
                             ": subir su cantidad en PUEBLOS")

# PASO 4: cada candidata que quedo tira 100 tiros nuevos
tiros_nuevos <- lapply(indice_candidata, function(p) rbinom(length(p), TIROS, p))


# ---------------------------------------------------------------------------
# RESULTADOS
# ---------------------------------------------------------------------------
cat("Pueblos que quedaron despues del filtro y encestes esperados en 100 tiros nuevos\n\n")
resumen <- data.frame(
  candidata   = names(tiros_nuevos),
  quedaron    = sapply(tiros_nuevos, length),
  encestes    = round(sapply(tiros_nuevos, mean), 1),
  margen_95   = round(1.96 * sapply(tiros_nuevos, sd) / sqrt(sapply(tiros_nuevos, length)), 1)
)
resumen <- resumen[order(resumen$encestes), ]
print(resumen, row.names = FALSE)

# el criterio del enunciado: P(Ci mete menos que Cj) en 100 tiros.
# Se estima enfrentando candidatas simuladas al azar, 200.000 veces por par.
p_menos <- function(a, b, veces = 2e5) {
  mean(a[sample.int(length(a), veces, replace = TRUE)] < b[sample.int(length(b), veces, replace = TRUE)])
}
orden <- resumen$candidata
P <- sapply(orden, function(j) sapply(orden, function(i) if (i == j) NA else p_menos(tiros_nuevos[[i]], tiros_nuevos[[j]])))
dimnames(P) <- list(orden, orden)
cat("\nP(la de la fila mete menos que la de la columna) en 100 tiros nuevos\n")
cat("(el enunciado dice Ci < Cj cuando esto es mayor a 0,5)\n\n")
print(round(P, 3))

cat("\nTiempo:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutos\n")
