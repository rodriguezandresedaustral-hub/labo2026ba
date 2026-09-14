# ===========================================================================
# 420_gridsearch_azarosos.R
#
# Grid search en DOS ETAPAS sobre el ensemble de arboles azarosos (z420),
# con la particion y la ganancia del 291. Pensado para correr desatendido en
# la VM Spot de GCP con:
#
#   nohup Rscript 420_gridsearch_azarosos.R >> ~/buckets/b1/exp/HT41/log.txt 2>&1 &
#
# COMO SE COMBINAN LOS 32 ARBOLES (regla B, neutral al peso):
#   a la probabilidad de cada arbol se le SACA el peso, r = p/(w - p*(w-1)),
#   se promedian las 32, se le VUELVE A PONER el peso al promedio, y se envia
#   si supera w/(w+39), igual que en el arbol solo. Asi el peso solo afecta que
#   cortes hace cada arbol, no la decision. Con w = 1 es identico al z420.
#
#   La regla anterior (promedio de probs ponderadas > w/(w+39), HT40) es exacta
#   para UN arbol pero no para el promedio: con peso 500 enviaba 1.299 clientes
#   contra 3.536 con peso 1, y la ganancia caia de 582 a 402 M. Medido en
#   test_regla_ensemble.R. Con la regla B: 582 vs 561. Las filas de peso 1 son
#   identicas bajo las dos reglas y se reutilizan.
#
# Si la VM se expropia, el startup script lo relanza y retoma donde quedo:
# cada combinacion se graba apenas termina, y cada etapa terminada se marca.
#
# Etapa 1: la grilla coarse (los valores del v1 del 291).
# Etapa 2: refinada a partir de la 1. Por cada eje, el centro es la MODA del
#          top 10 (no el puesto 1, que puede haber subido por suerte de la
#          particion), mas los puntos medios hacia sus dos vecinos; si el
#          centro quedo en un borde, se extiende para ese lado. El valor del
#          puesto 1 se agrega igual, para que el lider quede adentro.
#
# Manda un mail al arrancar, uno de avance cada 6 h con el tiempo que falta,
# y uno al terminar cada etapa (top 10 + grilla
# siguiente), por Gmail SMTP con la app password guardada en ~/.gmail_app_password.
# ===========================================================================

suppressMessages({
  require("data.table")
  require("rpart")
  require("parallel")
  if (!require("primes")) { install.packages("primes", repos= "https://cloud.r-project.org"); require("primes") }
})

# ---------------------------------------------------------------------------
# CONFIGURACION
# ---------------------------------------------------------------------------
PARAM <- list()

PARAM$bucket <- "/home/ds/buckets/b1"         # VM de la catedra.  Colab: "/content/buckets/b1"
PARAM$job    <- "HT41"                        # carpeta del trabajo: exp/HT41
PARAM$sembrar_desde <- "HT4010"               # filas de peso 1 de la corrida con regla A: valen igual
PARAM$etapas <- 2                             # cuantas etapas encadenar

PARAM$semilla_primigenia <- 100213
PARAM$training_pct <- 70L

PARAM$num_trees <- 32                         # exactamente 32 arboles por combinacion
PARAM$cores <- detectCores()

PARAM$cp <- -0.5                              # fijo: todos los negativos dan el mismo arbol

PARAM$email_desde <- "REMITENTE@gmail.com"         # la cuenta duena de la app password
PARAM$email_a     <- "DESTINO@gmail.com"                  # a donde llega
# app password de Gmail, guardada en la VM (no en el bucket): chmod 600
PARAM$gmail_pass_file <- "~/.gmail_app_password"
PARAM$mail_cada_horas <- 6                    # mail de avance con el tiempo que falta
PARAM$apagar_al_terminar <- TRUE              # la VM cobra por hora: se apaga sola al final

# etapa 1: la grilla coarse, los valores del v1 del 291
PARAM$grid_coarse <- list(
  peso_baja2       = c(1, 10, 100, 250, 500),
  feature_fraction = c(0.25, 0.5, 0.75),
  maxdepth         = c(4, 8, 12, 14),
  minsplit         = c(20, 100, 400, 1000),
  minbucket        = c(5, 20, 50, 100)
)

# limites duros de cada eje, para que la extension en el borde no se vaya
LIMITES <- list(
  peso_baja2       = c(1, 2000),
  feature_fraction = c(0.10, 1.00),
  maxdepth         = c(2, 30),
  minsplit         = c(2, 5000),
  minbucket        = c(1, 2000)
)
# ejes que se refinan en escala log (punto medio = media geometrica)
EJES_LOG <- c("peso_baja2", "minsplit", "minbucket")
EJES     <- names(PARAM$grid_coarse)

# ---------------------------------------------------------------------------
# CARPETAS Y MARCAS
# ---------------------------------------------------------------------------
job_dir <- file.path(PARAM$bucket, "exp", PARAM$job)
dir.create(job_dir, showWarnings= FALSE, recursive= TRUE)
carpeta_etapa <- function(k) file.path(PARAM$bucket, "exp", sprintf("%s%d0", PARAM$job, k))

log <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")

# marca de trabajo en curso: el startup script la usa para relanzar
writeLines(format(Sys.time()), file.path(job_dir, "CORRIENDO"))

# ---------------------------------------------------------------------------
# MAIL por Gmail SMTP con curl. Necesita ~/.gmail_app_password en la VM.
# ---------------------------------------------------------------------------
mandar_mail <- function(titulo, cuerpo) {
  pass_file <- path.expand(PARAM$gmail_pass_file)
  if (!file.exists(pass_file)) { log("SIN MAIL: no existe", pass_file); return(invisible(FALSE)) }
  pass <- trimws(readLines(pass_file, warn= FALSE)[1])

  archivo <- tempfile(fileext= ".txt")
  writeLines(c(
    paste0("From: ", PARAM$email_desde),
    paste0("To: ", PARAM$email_a),
    paste0("Subject: ", titulo),
    "Content-Type: text/plain; charset=UTF-8",
    "",
    cuerpo), archivo)

  rc <- system2("curl", c("-s", "--ssl-reqd",
    "--url", "smtps://smtp.gmail.com:465",
    "--mail-from", PARAM$email_desde,
    "--mail-rcpt", PARAM$email_a,
    "--user", shQuote(paste0(PARAM$email_desde, ":", pass)),
    "-T", archivo), stdout= FALSE, stderr= FALSE)
  unlink(archivo)
  if (rc == 0) log("mail enviado:", titulo) else log("FALLO el mail (curl rc =", rc, "):", titulo)
  invisible(rc == 0)
}

formatear_grilla <- function(grid) {
  paste(sapply(EJES, function(e) sprintf("  %-17s %s", e, paste(grid[[e]], collapse= ", "))),
        collapse= "\n")
}

formatear_top10 <- function(tb) {
  top <- head(tb[order(-ganancia_test)], 10)
  lineas <- sprintf("%2d. peso %-5g ff %-4g md %-2d ms %-4d mb %-3d -> %s",
    seq_len(nrow(top)), top$peso_baja2, top$feature_fraction, top$maxdepth,
    top$minsplit, top$minbucket,
    formatC(round(top$ganancia_test), format= "d", big.mark= ".", decimal.mark= ","))
  paste(lineas, collapse= "\n")
}

# ---------------------------------------------------------------------------
# PARTICION, igual que el 291
# ---------------------------------------------------------------------------
particionar <- function(data, division, agrupa= "", campo= "fold", start= 1, seed= NA) {
  if (!is.na(seed)) set.seed(seed)
  bloque <- unlist(mapply(function(x, y) rep(y, x),
    division, seq(from= start, length.out= length(division))))
  data[, (campo) := sample(rep(bloque, ceiling(.N / length(bloque))))[1:.N], by= agrupa]
}

primos <- generate_primes(min= 100000, max= 1000000)
set.seed(PARAM$semilla_primigenia)
PARAM$semilla_particion <- sample(primos, 1)            # 570937, igual que el 291
PARAM$semillas_arboles  <- sample(primos, PARAM$num_trees)
log("semilla de particion:", PARAM$semilla_particion)

dataset <- fread(file.path(PARAM$bucket, "datasets/dataset_pequeno.csv"))
dataset <- dataset[foto_mes == 202107]
particionar(dataset, division= c(PARAM$training_pct, 100L - PARAM$training_pct),
  agrupa= "clase_ternaria", seed= PARAM$semilla_particion)
dtrain <- dataset[fold == 1]
dtest  <- dataset[fold == 2]
campos_buenos <- copy(setdiff(colnames(dtrain), c("clase_ternaria", "fold")))
log("training:", nrow(dtrain), " testing:", nrow(dtest), " campos:", length(campos_buenos))

mandar_mail(sprintf("%s arranco", PARAM$job),
  paste0(PARAM$job, " arranco en ", format(Sys.time()), "\n",
         "training ", nrow(dtrain), " / testing ", nrow(dtest), "\n\n",
         "Grilla de la etapa 1:\n", formatear_grilla(PARAM$grid_coarse)))

# ---------------------------------------------------------------------------
# EL ENSEMBLE
# ---------------------------------------------------------------------------
UnArbol <- function(semilla, param) {
  set.seed(semilla)
  qty <- as.integer(length(campos_buenos) * param$feature_fraction)
  campos_random <- sample(campos_buenos, qty)
  formulita <- as.formula(paste0("clase_ternaria ~ ", paste(campos_random, collapse= " + ")))
  environment(formulita) <- environment()
  pesos <- dtrain[, ifelse(clase_ternaria == "BAJA+2", param$peso_baja2, 1.0)]
  modelo <- rpart(formulita, data= dtrain, xval= 0, weights= pesos,
    control= list(cp= param$cp, maxdepth= param$maxdepth,
                  minsplit= param$minsplit, minbucket= param$minbucket))
  predict(modelo, dtest, type= "prob")[, "BAJA+2"]
}

EnsembleEstimarGanancia <- function(param) {
  t0 <- Sys.time()
  probs <- mclapply(PARAM$semillas_arboles, UnArbol, param= param, mc.cores= PARAM$cores)
  w <- param$peso_baja2

  # regla B: sacar el peso a cada arbol, promediar, volver a ponderar el
  # promedio, y comparar contra w/(w+39) como en el arbol solo
  tasas <- lapply(probs, function(p) p / (w - p * (w - 1)))     # sin peso
  tasa_media <- Reduce(`+`, tasas) / length(tasas)             # promedio
  prob_media <- w * tasa_media / (w * tasa_media + 1 - tasa_media)   # con peso otra vez
  umbral <- w / (w + 39)
  ganancia_test <- dtest[, sum(ifelse(prob_media > umbral,
    ifelse(clase_ternaria == "BAJA+2", 975000, -25000), 0))]

  # regla A, solo como referencia (mismos arboles, no cuesta nada)
  prob_media_A <- Reduce(`+`, probs) / length(probs)
  umbral_A <- w / (w + 39)
  ganancia_A <- dtest[, sum(ifelse(prob_media_A > umbral_A,
    ifelse(clase_ternaria == "BAJA+2", 975000, -25000), 0))]
  data.table(
    semilla= PARAM$semilla_particion,
    peso_baja2= param$peso_baja2, feature_fraction= param$feature_fraction, cp= param$cp,
    maxdepth= param$maxdepth, minsplit= param$minsplit, minbucket= param$minbucket,
    num_trees= length(probs),
    ganancia_test= ganancia_test / ((100 - PARAM$training_pct) / 100),
    enviados= sum(prob_media > umbral),
    ganancia_reglaA= ganancia_A / ((100 - PARAM$training_pct) / 100),
    enviados_reglaA= sum(prob_media_A > umbral_A),
    segundos= as.numeric(difftime(Sys.time(), t0, units= "secs")))
}

# ---------------------------------------------------------------------------
# UNA ETAPA: recorre la grilla, graba cada combinacion, retoma si ya hay filas
# ---------------------------------------------------------------------------
correr_etapa <- function(grid, carpeta, sembrar_desde= NULL) {
  dir.create(carpeta, showWarnings= FALSE, recursive= TRUE)
  archivo <- file.path(carpeta, "gridsearch_detalle.txt")
  if (!file.exists(archivo) && !is.null(sembrar_desde)) {
    origen <- file.path(PARAM$bucket, "exp", sembrar_desde, "gridsearch_detalle.txt")
    if (file.exists(origen)) {
      previo <- fread(origen)[abs(peso_baja2 - 1) < 1e-9]     # solo peso 1: identicas bajo regla A y B
      previo[, ganancia_reglaA := ganancia_test]
      previo[, enviados_reglaA := enviados]
      fwrite(previo, file= archivo, sep= "\t")
      log("sembrado desde", sembrar_desde, "con", nrow(previo), "filas de peso 1")
    }
  }
  tb <- if (file.exists(archivo)) fread(archivo) else data.table()
  log("etapa en", basename(carpeta), "- filas al arrancar:", nrow(tb))

  # cuantas combinaciones reales tiene la grilla (sin los duplicados)
  total <- nrow(CJ(mb= grid$minbucket, ms= grid$minsplit)[ms >= 2L * mb]) *
    length(grid$peso_baja2) * length(grid$feature_fraction) * length(grid$maxdepth)

  combos <- 0; ajustados <- 0
  seg_sesion <- 0                       # segundos gastados en esta sesion
  t_ultimo_mail <- Sys.time()
  for (vpeso in grid$peso_baja2) for (vff in grid$feature_fraction)
  for (vmd in grid$maxdepth) for (vms in grid$minsplit) for (vmb in grid$minbucket) {
    if (vms < 2L * vmb) next                       # arbol identico: duplicado
    combos <- combos + 1
    if (nrow(tb) > 0 && nrow(tb[
        abs(peso_baja2 - vpeso) < 1e-9 & abs(feature_fraction - vff) < 1e-9 &
        maxdepth == vmd & minsplit == vms & minbucket == vmb &
        num_trees == PARAM$num_trees]) > 0) next  # ya calculada
    fila <- EnsembleEstimarGanancia(list(peso_baja2= vpeso, feature_fraction= vff,
      cp= PARAM$cp, maxdepth= vmd, minsplit= vms, minbucket= vmb))
    tb <- rbindlist(list(tb, fila), use.names= TRUE)
    ajustados <- ajustados + 1
    fwrite(tb, file= archivo, sep= "\t")          # graba en CADA combinacion
    log(basename(carpeta), sprintf("%d/%d", ajustados, combos),
        sprintf("peso %g ff %g md %d ms %d mb %d -> %s (%ds)", vpeso, vff, vmd, vms, vmb,
                formatC(round(fila$ganancia_test), format= "d", big.mark= ".", decimal.mark= ","),
                round(fila$segundos)))

    # mail de avance: al terminar la PRIMERA combinacion de la sesion (para
    # saber que todo anda) y despues cada PARAM$mail_cada_horas
    seg_sesion <- seg_sesion + fila$segundos
    if (ajustados == 1 ||
        difftime(Sys.time(), t_ultimo_mail, units= "hours") >= PARAM$mail_cada_horas) {
      faltan <- total - combos
      prom   <- seg_sesion / ajustados
      eta_h  <- faltan * prom / 3600
      cuerpo <- paste0(
        basename(carpeta), ": ", combos, " de ", total, " combinaciones (",
        round(100 * combos / total), "%).\n",
        "Promedio en esta sesion: ", round(prom), " s por combinacion.\n",
        "Faltan ", faltan, " -> ~", round(eta_h, 1), " horas para terminar la etapa.\n",
        "Termina aprox: ", format(Sys.time() + eta_h * 3600, "%a %d %H:%M"), "\n\n",
        "Top 3 hasta ahora:\n",
        paste(head(strsplit(formatear_top10(tb), "\n")[[1]], 3), collapse= "\n"))
      mandar_mail(sprintf("%s avance: %d/%d, faltan ~%.1f h", basename(carpeta), combos, total, eta_h), cuerpo)
      t_ultimo_mail <- Sys.time()
    }
  }
  log("etapa", basename(carpeta), "terminada:", ajustados, "ajustados de", combos)
  tb_rank <- tb[order(-ganancia_test)]; tb_rank[, ranking := .I]
  fwrite(tb_rank, file= file.path(carpeta, "ranking.txt"), sep= "\t")
  tb_rank
}

# ---------------------------------------------------------------------------
# LA GRILLA SIGUIENTE, a partir del ranking de la etapa anterior
#
# Por eje se refina alrededor de dos centros: la MODA del top 10 (el valor
# que mas aparece entre las diez mejores, desempate por ganancia promedio) y
# el LIDER (puesto 1). Cuando coinciden es un solo centro. Cuando el top 10
# esta partido entre dos zonas, cubre las dos en vez de elegir una a ciegas.
# Alrededor de cada centro: los puntos medios hacia sus vecinos, o extension
# para afuera si esta en un borde. maxdepth es la excepcion: se prueban los
# enteros a +-2 del centro (sin pasar los vecinos), porque el paso par
# escondio el optimo (5).
# ---------------------------------------------------------------------------
medio <- function(a, b, eje) if (eje %in% EJES_LOG) sqrt(a * b) else (a + b) / 2

redondear <- function(x, eje) {
  x <- pmin(pmax(x, LIMITES[[eje]][1]), LIMITES[[eje]][2])
  if (eje == "feature_fraction") round(x, 2) else as.integer(round(x))
}

derivar_grilla <- function(tb_rank, grid_anterior) {
  top <- head(tb_rank[order(-ganancia_test)], 10)
  lider <- tb_rank[order(-ganancia_test)][1]
  nueva <- list()
  for (e in EJES) {
    v <- sort(unique(grid_anterior[[e]])); n <- length(v)
    if (n == 1) { nueva[[e]] <- v; next }

    # la moda del top 10 en este eje; empate -> mayor ganancia promedio
    moda <- top[, .(veces= .N, g= mean(ganancia_test)), by= e][order(-veces, -g)][1][[e]]

    # se refina alrededor de DOS centros: la moda y el lider. Cuando coinciden
    # es un solo centro; cuando el top 10 esta partido, cubre los dos lados.
    vals <- c()
    for (w in unique(c(moda, lider[[e]]))) {
      i <- which(v == w)
      if (e == "maxdepth") {
        # entero de rango chico: TODOS los enteros entre los vecinos, no solo
        # el punto medio. Seis corridas con pares escondieron que el optimo era 5.
        lo <- if (i == 1) v[1] - (v[2] - v[1]) else v[i - 1]
        hi <- if (i == n) v[n] + (v[n] - v[n - 1]) else v[i + 1]
        vals <- c(vals, seq(max(lo, w - 2), min(hi, w + 2)))   # +-2 del centro, sin pasar los vecinos
      } else if (i == 1) {                           # borde inferior: extiende abajo
        paso <- if (e %in% EJES_LOG) v[2] / v[1] else v[2] - v[1]
        abajo <- if (e %in% EJES_LOG) v[1] / paso else v[1] - paso
        vals <- c(vals, abajo, w, medio(w, v[2], e))
      } else if (i == n) {                           # borde superior: extiende arriba
        paso <- if (e %in% EJES_LOG) v[n] / v[n - 1] else v[n] - v[n - 1]
        arriba <- if (e %in% EJES_LOG) v[n] * paso else v[n] + paso
        vals <- c(vals, medio(v[n - 1], w, e), w, arriba)
      } else {                                       # interior: los dos puntos medios
        vals <- c(vals, medio(v[i - 1], w, e), w, medio(w, v[i + 1], e))
      }
    }
    nueva[[e]] <- sort(unique(redondear(vals, e)))
  }
  nueva
}

# ---------------------------------------------------------------------------
# MAIN: las etapas encadenadas
# ---------------------------------------------------------------------------
grid <- PARAM$grid_coarse
mejor_anterior <- NA_real_

for (k in seq_len(PARAM$etapas)) {
  carpeta <- carpeta_etapa(k)
  marca   <- file.path(carpeta, "ETAPA_TERMINADA")

  log("========== ETAPA", k, "de", PARAM$etapas, "==========")
  log("grilla:\n", formatear_grilla(grid))

  if (file.exists(marca)) {
    log("etapa", k, "ya estaba terminada, la salteo")
    tb_rank <- fread(file.path(carpeta, "ranking.txt"))
  } else {
    tb_rank <- correr_etapa(grid, carpeta,
      sembrar_desde= if (k == 1) PARAM$sembrar_desde else NULL)
    writeLines(format(Sys.time()), marca)
  }

  mejor <- tb_rank[1, ganancia_test]
  mejora <- if (is.na(mejor_anterior)) NA else mejor - mejor_anterior
  mejor_anterior <- mejor

  # la grilla de la etapa siguiente, o fin
  if (k < PARAM$etapas) {
    grid_sig <- derivar_grilla(tb_rank, grid)
    texto_sig <- paste0("Grilla de la etapa ", k + 1, ":\n", formatear_grilla(grid_sig))
  } else {
    grid_sig <- NULL
    texto_sig <- "Ultima etapa. Fin."
  }

  cuerpo <- paste0(
    PARAM$job, " etapa ", k, " terminada. ", nrow(tb_rank), " combinaciones.\n",
    if (!is.na(mejora)) sprintf("Mejora sobre la etapa anterior: %s\n",
      formatC(round(mejora), format= "d", big.mark= ".", decimal.mark= ",")) else "",
    "\nTop 10:\n", formatear_top10(tb_rank), "\n\n",
    "Grilla que se uso en esta etapa:\n", formatear_grilla(grid), "\n\n",
    texto_sig)
  mandar_mail(sprintf("%s etapa %d/%d terminada", PARAM$job, k, PARAM$etapas), cuerpo)
  writeLines(cuerpo, file.path(carpeta, "resumen.txt"))

  if (!is.null(grid_sig)) grid <- grid_sig
}

writeLines(format(Sys.time()), file.path(job_dir, "TERMINADO.txt"))
unlink(file.path(job_dir, "CORRIENDO"))
log("TRABAJO TERMINADO")

# la VM no es Spot y cobra por hora: se apaga sola al terminar.
# Todo lo que importa ya esta en el bucket; el mail final ya salio.
if (isTRUE(PARAM$apagar_al_terminar)) {
  log("apagando la VM en 60 segundos")
  Sys.sleep(60)
  system("sudo shutdown -h now")
}
