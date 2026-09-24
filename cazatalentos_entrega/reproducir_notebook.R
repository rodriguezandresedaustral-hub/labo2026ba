# ===========================================================================
# reproducir_notebook.R
#
# Control de que la forma de simular es la misma que usa el notebook de la
# catedra. Corre el codigo del notebook tal cual, con su semilla (102191),
# y compara con los numeros que el notebook escribe en su texto:
#
#   - seccion 2.3: con 10 tiros, la verdadera mejor (0,70) gana apenas
#     1450 de 10000 veces contra un peloton de 99 jugadoras de 0,501 a 0,599
#   - seccion 2.7: con 415 tiros, la verdadera mejor gana el 99,10% de las
#     veces
#
# Si estos dos numeros salen iguales, la simulacion de
# cazatalentos_simulacion.R esta hecha con la misma mecanica que el notebook.
#
#   Rscript reproducir_notebook.R      (tarda unos 10 minutos)
# ===========================================================================

# la funcion del notebook: cuantos encestes hace una jugadora de indice
# 'prob' en 'qty' tiros. Cada tiro es un numero al azar entre 0 y 1; si cae
# por debajo del indice, es enceste
ftirar <- function(prob, qty) sum(runif(qty) < prob)


# ---------------------------------------------------------------------------
# Seccion 2.3: 100 jugadoras, 10 tiros cada una, 10000 torneos.
# Cuantas veces gana la verdadera mejor (la jugadora 1, con 0,70)
# ---------------------------------------------------------------------------
set.seed(102191)
jugadoras <- c(0.7, (501:599) / 1000)
invisible(mapply(ftirar, jugadoras, 10))           # el notebook hace una tirada de muestra antes: se repite para que la semilla quede igual

primera_ganadora <- 0
for (i in 1:10000) {
  aciertos <- mapply(ftirar, jugadoras, 10)
  if (which.max(aciertos) == 1) primera_ganadora <- primera_ganadora + 1
}
cat("Seccion 2.3: la mejor gana", primera_ganadora, "de 10000   (el notebook dice 1450)\n")


# ---------------------------------------------------------------------------
# Seccion 2.7: la misma idea con mas tiros. La mejor va al final de la
# lista porque which.max() se queda con la primera si hay empate.
# El notebook recorre varias cantidades de tiros en orden; se recorren las
# mismas para que la semilla llegue igual a 415
# ---------------------------------------------------------------------------
set.seed(102191)
jugadoras <- c((501:599) / 1000, 0.7)
for (tiros in c(10, 20, 50, 100, 200, 300, 400, 415)) {
  gana_la_mejor <- 0
  for (i in 1:10000) {
    if (which.max(mapply(ftirar, jugadoras, tiros)) == 100) gana_la_mejor <- gana_la_mejor + 1
  }
  if (tiros == 415) cat("Seccion 2.7: con 415 tiros la mejor gana el", gana_la_mejor / 10000, "  (el notebook dice 0,9910)\n")
}
