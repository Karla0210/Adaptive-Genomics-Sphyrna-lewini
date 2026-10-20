# ============================================================================
# CKMRsim: integrated kinship analysis and final figures
# Neutral panel of 9,529 SNPs from Sphyrna lewini; 38 juveniles
# Comparisons: FS/U, HS/U, and FS/HS
# ============================================================================

# This script performs the following steps in a single workflow:
#   1) VCF import and quality control;
#   2) genotype and allele-frequency preparation;
#   3) missing-data assessment per individual and pair;
#   4) CKMRsim simulations calibrated for two missing-data strata;
#   5) threshold estimation with FNR = 0.01;
#   6) evaluation of the 703 observed pairs;
#   7) within- and between-locality summaries and a permutation test;
#   8) generation of three individual figures and one combined figure.
#
# IMPORTANT:
# - Both missing-data strata belong to the same panel of 9,529 SNPs.
# - Thresholds and classification preserve both missing-data calibrations.
# - The figures use a single simulation of the complete panel of 9,502
#   polymorphic SNPs, so each relationship is represented by one curve.
# - If the VCF, PopMap, or any parameter changes, set
#   reutilizar_calculos <- FALSE before rerunning the script.


# 00. INSTALLATION: run only once ---------------------------------------------

# install.packages(c(
#   "remotes", "vcfR", "dplyr", "tidyr", "stringr", "tibble",
#   "ggplot2", "patchwork"
# ))
# remotes::install_github("eriqande/CKMRsim", build_vignettes = TRUE)


# 01. SESSION SETUP, PACKAGES, AND PARAMETERS ---------------------------------

rm(list = ls())
graphics.off()
gc()

library(vcfR)
library(CKMRsim)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)

# Change these paths only if the input files are stored elsewhere.
carpeta_datos <- paste0(
  "E:/Filtrado populations y VCF p7 r8/",
  "Filtrados_new/CKMRsim/"
)

archivo_vcf <- file.path(
  carpeta_datos,
  "Resumen_SNPn2.recode.vcf"
)

archivo_popmap <- file.path(
  carpeta_datos,
  "PopMap.txt"
)

carpeta_salida <- file.path(
  carpeta_datos,
  "salidas_CKMRsim_9529"
)

dir.create(
  carpeta_salida,
  recursive = TRUE,
  showWarnings = FALSE
)

# Parameters used in the final analysis.
numero_snp_entrada <- 9529L
numero_snp_polimorficos_esperado <- 9502L
numero_individuos_esperado <- 38L
umbral_missing_individual <- 0.30
epsilon_genotipificacion <- 0.005
numero_replicas <- 10000L
FNR_objetivo <- 0.01
numero_permutaciones <- 9999L
numero_cores <- 1L

# TRUE reuses existing RDS files to speed up subsequent runs.
# Set to FALSE whenever data or analysis parameters change.
reutilizar_calculos <- TRUE

stopifnot(
  file.exists(archivo_vcf),
  file.exists(archivo_popmap)
)

message("CKMRsim version: ", packageVersion("CKMRsim"))
message("Output directory: ", carpeta_salida)


# 02. IMPORT VCF AND POPMAP ---------------------------------------------------

vcf_9529 <- read.vcfR(
  archivo_vcf,
  verbose = FALSE
)

popmap <- read.table(
  archivo_popmap,
  header = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(popmap) < 2L) {
  stop("The PopMap must contain at least two columns: individual and locality.")
}

popmap <- popmap[, 1:2, drop = FALSE]
colnames(popmap) <- c("Indiv", "Localidad")

if (anyDuplicated(popmap$Indiv)) {
  stop("The PopMap contains duplicate individual IDs.")
}

GT_9529 <- extract.gt(
  vcf_9529,
  element = "GT",
  as.numeric = FALSE
)

stopifnot(
  nrow(vcf_9529@fix) == numero_snp_entrada,
  ncol(vcf_9529@gt) - 1L == numero_individuos_esperado,
  setequal(colnames(GT_9529), popmap$Indiv)
)

message(
  "VCF loaded: ",
  nrow(vcf_9529@fix),
  " SNPs and ",
  ncol(vcf_9529@gt) - 1L,
  " individuals."
)


# 03. BIALLELIC SNP AND MISSING-DATA QUALITY CONTROL --------------------------

ref_9529 <- vcf_9529@fix[, "REF"]
alt_9529 <- vcf_9529@fix[, "ALT"]

snp_bialelico_9529 <- (
  nchar(ref_9529) == 1L &
    nchar(alt_9529) == 1L &
    !grepl(",", alt_9529)
)

if (!all(snp_bialelico_9529)) {
  stop("The VCF contains variants that are not biallelic SNPs.")
}

GT_missing_9529 <- (
  is.na(GT_9529) |
    grepl("\\.", GT_9529)
)

missing_individuos_9529 <- tibble(
  Indiv = colnames(GT_9529),
  missing = colMeans(GT_missing_9529)
) %>%
  left_join(popmap, by = "Indiv") %>%
  arrange(desc(missing))

missing_locus_9529 <- rowMeans(GT_missing_9529)

message(
  "Overall proportion of missing genotypes: ",
  round(mean(GT_missing_9529), 6)
)

print(missing_individuos_9529, n = Inf, width = Inf)
print(summary(missing_locus_9529))


# 04. CONVERT GENOTYPES TO LONG FORMAT ----------------------------------------

info_loci_9529 <- tibble(
  Chrom = as.character(vcf_9529@fix[, "CHROM"]),
  Pos = as.numeric(vcf_9529@fix[, "POS"]),
  Locus = paste(
    vcf_9529@fix[, "CHROM"],
    vcf_9529@fix[, "POS"],
    sep = ":"
  ),
  REF = ref_9529,
  ALT = alt_9529
)

if (anyDuplicated(info_loci_9529$Locus)) {
  stop("Duplicate CHROM:POS identifiers were found in the VCF.")
}

rownames(GT_9529) <- info_loci_9529$Locus
rownames(GT_missing_9529) <- info_loci_9529$Locus

long_genos_9529 <- as.data.frame(
  GT_9529,
  check.names = FALSE
) %>%
  rownames_to_column("Locus") %>%
  pivot_longer(
    cols = -Locus,
    names_to = "Indiv",
    values_to = "GT"
  ) %>%
  mutate(
    GT = str_replace_all(GT, "\\|", "/"),
    GT = if_else(
      is.na(GT) | grepl("\\.", GT),
      NA_character_,
      GT
    ),
    A1 = if_else(
      is.na(GT),
      NA_character_,
      str_extract(GT, "^[01]")
    ),
    A2 = if_else(
      is.na(GT),
      NA_character_,
      str_extract(GT, "[01]$")
    )
  ) %>%
  select(Locus, Indiv, A1, A2) %>%
  pivot_longer(
    cols = c(A1, A2),
    names_to = "gene_copy",
    values_to = "GT"
  ) %>%
  left_join(
    info_loci_9529 %>% select(Locus, REF, ALT),
    by = "Locus"
  ) %>%
  mutate(
    Allele = case_when(
      GT == "0" ~ REF,
      GT == "1" ~ ALT,
      TRUE ~ NA_character_
    ),
    # CKMRsim expects gene copies identified as "1" and "2".
    gene_copy = str_remove(gene_copy, "^A")
  ) %>%
  select(Indiv, Locus, gene_copy, Allele)

stopifnot(
  nrow(long_genos_9529) ==
    numero_snp_entrada * numero_individuos_esperado * 2L,
  n_distinct(long_genos_9529$Indiv) == numero_individuos_esperado,
  n_distinct(long_genos_9529$Locus) == numero_snp_entrada,
  all(table(long_genos_9529$Indiv, long_genos_9529$Locus) == 2L),
  setequal(unique(long_genos_9529$gene_copy), c("1", "2"))
)


# 05. RETAIN POLYMORPHIC LOCI -------------------------------------------------

loci_polimorficos_9529 <- long_genos_9529 %>%
  filter(!is.na(Allele)) %>%
  distinct(Locus, Allele) %>%
  count(Locus, name = "n_alelos_observados") %>%
  filter(n_alelos_observados == 2L) %>%
  pull(Locus)

long_genos_filtrado_9529 <- long_genos_9529 %>%
  filter(Locus %in% loci_polimorficos_9529)

stopifnot(
  length(loci_polimorficos_9529) ==
    numero_snp_polimorficos_esperado,
  nrow(long_genos_filtrado_9529) ==
    numero_snp_polimorficos_esperado * numero_individuos_esperado * 2L
)

message(
  "Polymorphic loci used by CKMRsim: ",
  length(loci_polimorficos_9529)
)


# 06. COMPUTE AND REINDEX ALLELE FREQUENCIES ----------------------------------

allele_freqs_9529 <- long_genos_filtrado_9529 %>%
  filter(!is.na(Allele)) %>%
  count(Locus, Allele, name = "n_alelos") %>%
  group_by(Locus) %>%
  mutate(Freq = n_alelos / sum(n_alelos)) %>%
  ungroup() %>%
  left_join(
    info_loci_9529 %>% select(Chrom, Pos, Locus),
    by = "Locus"
  ) %>%
  select(Chrom, Pos, Locus, Allele, Freq) %>%
  arrange(Chrom, Pos, desc(Freq)) %>%
  mutate(
    AlleIdx = NA_integer_,
    LocIdx = NA_integer_
  )

control_suma_frecuencias_9529 <- allele_freqs_9529 %>%
  group_by(Locus) %>%
  summarise(suma = sum(Freq), .groups = "drop")

stopifnot(all(abs(control_suma_frecuencias_9529$suma - 1) < 1e-10))

afreqs_ready_9529 <- reindex_markers(allele_freqs_9529)

stopifnot(
  nrow(afreqs_ready_9529) ==
    numero_snp_polimorficos_esperado * 2L,
  n_distinct(afreqs_ready_9529$Locus) ==
    numero_snp_polimorficos_esperado,
  n_distinct(afreqs_ready_9529$LocIdx) ==
    numero_snp_polimorficos_esperado,
  all(table(afreqs_ready_9529$Locus) == 2L)
)

saveRDS(
  long_genos_filtrado_9529,
  file.path(carpeta_salida, "long_genos_filtrado_9529.rds")
)

saveRDS(
  afreqs_ready_9529,
  file.path(carpeta_salida, "frecuencias_CKMRsim_9529.rds")
)


# 07. PAIRWISE COMPARABLE LOCI AND MISSING-DATA STRATA ------------------------

orden_loci_CKMRsim_9529 <- afreqs_ready_9529 %>%
  distinct(LocIdx, Locus) %>%
  arrange(LocIdx) %>%
  pull(Locus)

missing_GT_filtrado_9529 <- GT_missing_9529[
  orden_loci_CKMRsim_9529,
  ,
  drop = FALSE
]

combinaciones_9529 <- combn(
  colnames(missing_GT_filtrado_9529),
  2
)

n_loci_comparables_9529 <- vapply(
  seq_len(ncol(combinaciones_9529)),
  function(i) {
    individuo1 <- combinaciones_9529[1, i]
    individuo2 <- combinaciones_9529[2, i]
    
    sum(
      !missing_GT_filtrado_9529[, individuo1] &
        !missing_GT_filtrado_9529[, individuo2]
    )
  },
  integer(1)
)

individuos_missing_muy_alto_9529 <- missing_individuos_9529 %>%
  filter(missing > umbral_missing_individual) %>%
  pull(Indiv)

loci_por_par_9529 <- tibble(
  Indiv1 = combinaciones_9529[1, ],
  Indiv2 = combinaciones_9529[2, ],
  n_loci_comparables = n_loci_comparables_9529,
  n_loci_faltantes =
    numero_snp_polimorficos_esperado - n_loci_comparables_9529
) %>%
  mutate(
    grupo_missing = if_else(
      Indiv1 %in% individuos_missing_muy_alto_9529 |
        Indiv2 %in% individuos_missing_muy_alto_9529,
      "Missing_alto",
      "Missing_normal"
    ),
    Pareja1 = pmin(Indiv1, Indiv2),
    Pareja2 = pmax(Indiv1, Indiv2)
  ) %>%
  select(
    Indiv1 = Pareja1,
    Indiv2 = Pareja2,
    n_loci_comparables,
    n_loci_faltantes,
    grupo_missing
  )

# Verify the 38 choose 2 = 703 pairs.
stopifnot(
  nrow(loci_por_par_9529) == choose(numero_individuos_esperado, 2),
  !anyDuplicated(loci_por_par_9529[c("Indiv1", "Indiv2")])
)

resumen_grupos_missing_9529 <- loci_por_par_9529 %>%
  group_by(grupo_missing) %>%
  summarise(
    n_parejas = n(),
    minimo = min(n_loci_comparables),
    Q1 = quantile(n_loci_comparables, 0.25),
    mediana = median(n_loci_comparables),
    promedio = mean(n_loci_comparables),
    Q3 = quantile(n_loci_comparables, 0.75),
    maximo = max(n_loci_comparables),
    .groups = "drop"
  )

print(resumen_grupos_missing_9529, n = Inf, width = Inf)


# 08. EMPIRICAL PATTERN OF MISSING LOCI BY STRATUM ----------------------------

calcular_parametros_missing <- function(tabla_pares, matriz_missing) {
  if (nrow(tabla_pares) == 0L) {
    stop("One of the missing-data strata contains no pairs.")
  }
  
  conteo_missing <- numeric(nrow(matriz_missing))
  
  for (i in seq_len(nrow(tabla_pares))) {
    conteo_missing <- conteo_missing + (
      matriz_missing[, tabla_pares$Indiv1[i]] |
        matriz_missing[, tabla_pares$Indiv2[i]]
    )
  }
  
  pesos <- conteo_missing / nrow(tabla_pares)
  
  if (sum(pesos) == 0) {
    pesos <- rep(1, length(pesos))
  }
  
  list(
    rando_miss_n = as.integer(round(
      median(tabla_pares$n_loci_faltantes)
    )),
    rando_miss_wts = pesos
  )
}

pares_missing_normal_9529 <- loci_por_par_9529 %>%
  filter(grupo_missing == "Missing_normal")

pares_missing_alto_9529 <- loci_por_par_9529 %>%
  filter(grupo_missing == "Missing_alto")

param_missing_normal_9529 <- calcular_parametros_missing(
  pares_missing_normal_9529,
  missing_GT_filtrado_9529
)

param_missing_alto_9529 <- calcular_parametros_missing(
  pares_missing_alto_9529,
  missing_GT_filtrado_9529
)

message(
  "Loci removed per replicate, normal missingness: ",
  param_missing_normal_9529$rando_miss_n
)

message(
  "Loci removed per replicate, high missingness: ",
  param_missing_alto_9529$rando_miss_n
)


# 09. CREATE CKMRsim OBJECT AND CHECK FOR DUPLICATES --------------------------

relaciones_CKMRsim <- c("FS", "HS", "U")

ckmr_slew_9529 <- create_ckmr(
  D = afreqs_ready_9529,
  kappa_matrix = kappas[c("MZ", relaciones_CKMRsim), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(
    epsilon = epsilon_genotipificacion
  ),
  ge_mod_true_pars_list = list(
    epsilon = epsilon_genotipificacion
  )
)

print(ckmr_slew_9529)

duplicados_9529 <- find_close_matching_genotypes(
  LG = long_genos_filtrado_9529,
  CK = ckmr_slew_9529,
  max_mismatch = 6
)

print(duplicados_9529, n = Inf, width = Inf)

saveRDS(
  ckmr_slew_9529,
  file.path(carpeta_salida, "objeto_CKMRsim_9529.rds")
)


# 10. FINAL SIMULATIONS BY MISSING-DATA STRATUM -------------------------------

archivo_sim_normal_9529 <- file.path(
  carpeta_salida,
  "Simulacion_CKMRsim_9529_missing_normal.rds"
)

archivo_sim_alto_9529 <- file.path(
  carpeta_salida,
  "Simulacion_CKMRsim_9529_missing_alto.rds"
)

if (
  reutilizar_calculos &&
  file.exists(archivo_sim_normal_9529) &&
  file.exists(archivo_sim_alto_9529)
) {
  message("Loading saved CKMRsim simulations.")
  sim_normal_9529 <- readRDS(archivo_sim_normal_9529)
  sim_alto_9529 <- readRDS(archivo_sim_alto_9529)
} else {
  message("Running simulation for the normal missingness stratum.")
  set.seed(123)
  
  sim_normal_9529 <- simulate_Qij(
    C = ckmr_slew_9529,
    calc_relats = relaciones_CKMRsim,
    sim_relats = relaciones_CKMRsim,
    reps = numero_replicas,
    rando_miss_n = param_missing_normal_9529$rando_miss_n,
    rando_miss_wts = param_missing_normal_9529$rando_miss_wts
  )
  
  saveRDS(sim_normal_9529, archivo_sim_normal_9529)
  
  message("Running simulation for the high missingness stratum.")
  set.seed(456)
  
  sim_alto_9529 <- simulate_Qij(
    C = ckmr_slew_9529,
    calc_relats = relaciones_CKMRsim,
    sim_relats = relaciones_CKMRsim,
    reps = numero_replicas,
    rando_miss_n = param_missing_alto_9529$rando_miss_n,
    rando_miss_wts = param_missing_alto_9529$rando_miss_wts
  )
  
  saveRDS(sim_alto_9529, archivo_sim_alto_9529)
}

print(sim_normal_9529)
print(sim_alto_9529)


# 11. CLASSIFICATION POWER AND THRESHOLDS AT FNR = 0.01 -----------------------

obtener_umbrales <- function(simulacion, grupo, FNR_seleccionado) {
  resultados <- bind_rows(
    mc_sample_simple(simulacion, nu = "FS", de = "U") %>%
      mutate(comparacion = "FS/U"),
    mc_sample_simple(simulacion, nu = "HS", de = "U") %>%
      mutate(comparacion = "HS/U"),
    mc_sample_simple(simulacion, nu = "FS", de = "HS") %>%
      mutate(comparacion = "FS/HS")
  ) %>%
    filter(near(FNR, FNR_seleccionado)) %>%
    mutate(grupo_missing = grupo) %>%
    select(
      grupo_missing,
      comparacion,
      Lambda_star,
      FNR,
      FPR,
      everything()
    )
  
  if (nrow(resultados) != 3L) {
    stop(
      "Exactly one threshold per comparison was not obtained for ",
      grupo,
      "."
    )
  }
  
  resultados
}

umbrales_normal_9529 <- obtener_umbrales(
  sim_normal_9529,
  "Missing_normal",
  FNR_objetivo
)

umbrales_alto_9529 <- obtener_umbrales(
  sim_alto_9529,
  "Missing_alto",
  FNR_objetivo
)

umbrales_9529 <- bind_rows(
  umbrales_normal_9529,
  umbrales_alto_9529
)

stopifnot(
  nrow(umbrales_9529) == 6L,
  all(table(
    umbrales_9529$grupo_missing,
    umbrales_9529$comparacion
  ) == 1L)
)

print(
  umbrales_9529 %>%
    select(
      grupo_missing,
      comparacion,
      Lambda_star,
      FNR,
      FPR
    ),
  n = Inf,
  width = Inf
)


# 12. LIKELIHOOD RATIOS FOR THE 703 OBSERVED PAIRS ----------------------------

archivo_LLR_reales_9529 <- file.path(
  carpeta_salida,
  "LLR_parejas_reales_9529.rds"
)

calcular_LLR_parejas <- function(
    numerador,
    denominador,
    nombre_LLR,
    nombre_num_loci
) {
  salida <- pairwise_kin_logl_ratios(
    D1 = long_genos_filtrado_9529,
    D2 = long_genos_filtrado_9529,
    CK = ckmr_slew_9529,
    numer = numerador,
    denom = denominador,
    num_cores = numero_cores
  ) %>%
    transmute(
      Indiv1 = pmin(D1_indiv, D2_indiv),
      Indiv2 = pmax(D1_indiv, D2_indiv),
      LLR = logl_ratio,
      n_loci = num_loc
    )
  
  names(salida)[3:4] <- c(nombre_LLR, nombre_num_loci)
  salida
}

if (reutilizar_calculos && file.exists(archivo_LLR_reales_9529)) {
  message("Loading saved LLR values for observed pairs.")
  resultados_LLR_9529 <- readRDS(archivo_LLR_reales_9529)
} else {
  message("Calculating FS/U for observed pairs.")
  FS_U_real_9529 <- calcular_LLR_parejas(
    "FS", "U", "FS_U", "num_loc_FS_U"
  )
  
  message("Calculating HS/U for observed pairs.")
  HS_U_real_9529 <- calcular_LLR_parejas(
    "HS", "U", "HS_U", "num_loc_HS_U"
  )
  
  message("Calculating FS/HS for observed pairs.")
  FS_HS_real_9529 <- calcular_LLR_parejas(
    "FS", "HS", "FS_HS", "num_loc_FS_HS"
  )
  
  resultados_LLR_9529 <- FS_U_real_9529 %>%
    full_join(
      HS_U_real_9529,
      by = c("Indiv1", "Indiv2")
    ) %>%
    full_join(
      FS_HS_real_9529,
      by = c("Indiv1", "Indiv2")
    )
  
  saveRDS(resultados_LLR_9529, archivo_LLR_reales_9529)
}

stopifnot(
  nrow(resultados_LLR_9529) ==
    choose(numero_individuos_esperado, 2),
  !anyDuplicated(resultados_LLR_9529[c("Indiv1", "Indiv2")])
)

resultados_pares_9529 <- resultados_LLR_9529 %>%
  left_join(
    loci_por_par_9529,
    by = c("Indiv1", "Indiv2")
  )

stopifnot(
  nrow(resultados_pares_9529) ==
    choose(numero_individuos_esperado, 2),
  !any(is.na(resultados_pares_9529$grupo_missing))
)


# 13. APPLY STRATUM-SPECIFIC THRESHOLDS AND CLASSIFY PAIRS --------------------

umbrales_wide_9529 <- umbrales_9529 %>%
  select(grupo_missing, comparacion, Lambda_star) %>%
  pivot_wider(
    names_from = comparacion,
    values_from = Lambda_star,
    names_prefix = "umbral_"
  ) %>%
  rename(
    umbral_FS_U = `umbral_FS/U`,
    umbral_HS_U = `umbral_HS/U`,
    umbral_FS_HS = `umbral_FS/HS`
  )

resultados_pares_9529 <- resultados_pares_9529 %>%
  left_join(
    umbrales_wide_9529,
    by = "grupo_missing"
  ) %>%
  mutate(
    relacion_asignada = case_when(
      FS_U >= umbral_FS_U &
        FS_HS >= umbral_FS_HS ~ "Full siblings",
      HS_U >= umbral_HS_U &
        FS_HS < umbral_FS_HS ~ "Half siblings",
      FS_U < umbral_FS_U &
        HS_U < umbral_HS_U ~ "Unrelated",
      TRUE ~ "Uncertain"
    )
  )

mapa_localidades <- setNames(
  popmap$Localidad,
  popmap$Indiv
)

resultados_pares_9529 <- resultados_pares_9529 %>%
  mutate(
    Localidad1 = unname(mapa_localidades[Indiv1]),
    Localidad2 = unname(mapa_localidades[Indiv2]),
    distribucion = if_else(
      Localidad1 == Localidad2,
      "Dentro de localidad",
      "Entre localidades"
    ),
    FS_U_por_locus = FS_U / num_loc_FS_U,
    HS_U_por_locus = HS_U / num_loc_HS_U,
    FS_HS_por_locus = FS_HS / num_loc_FS_HS,
    margen_FS_U = FS_U - umbral_FS_U,
    margen_HS_U = HS_U - umbral_HS_U,
    margen_FS_HS = FS_HS - umbral_FS_HS,
    margen_HS_U_por_locus =
      (HS_U - umbral_HS_U) / num_loc_HS_U
  )

pares_emparentados_9529 <- resultados_pares_9529 %>%
  filter(
    relacion_asignada %in%
      c("Full siblings", "Half siblings")
  ) %>%
  arrange(
    relacion_asignada,
    desc(FS_U),
    desc(HS_U)
  )

resumen_final_9529 <- resultados_pares_9529 %>%
  count(
    grupo_missing,
    relacion_asignada,
    name = "n_parejas"
  )

print(
  table(
    resultados_pares_9529$relacion_asignada,
    useNA = "ifany"
  )
)

print(resumen_final_9529, n = Inf, width = Inf)
print(pares_emparentados_9529, n = Inf, width = Inf)


# 14. WITHIN- AND BETWEEN-LOCALITY SUMMARIES ----------------------------------

resumen_parentesco_localidad_9529 <- resultados_pares_9529 %>%
  filter(Localidad1 == Localidad2) %>%
  group_by(Localidad = Localidad1) %>%
  summarise(
    n_individuos = n_distinct(c(Indiv1, Indiv2)),
    n_parejas = n(),
    max_FS_U = max(FS_U, na.rm = TRUE),
    promedio_FS_U = mean(FS_U, na.rm = TRUE),
    max_HS_U = max(HS_U, na.rm = TRUE),
    promedio_HS_U = mean(HS_U, na.rm = TRUE),
    max_FS_HS = max(FS_HS, na.rm = TRUE),
    n_FS = sum(relacion_asignada == "Full siblings"),
    n_HS = sum(relacion_asignada == "Half siblings"),
    .groups = "drop"
  ) %>%
  arrange(desc(max_HS_U))

parejas_mas_cercanas_localidad_9529 <- resultados_pares_9529 %>%
  filter(Localidad1 == Localidad2) %>%
  group_by(Localidad = Localidad1) %>%
  slice_max(
    order_by = HS_U,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    Localidad,
    Indiv1,
    Indiv2,
    n_loci_comparables,
    FS_U,
    HS_U,
    FS_HS,
    umbral_HS_U,
    relacion_asignada
  ) %>%
  arrange(desc(HS_U))

resumen_estandarizado_localidad_9529 <- resultados_pares_9529 %>%
  filter(Localidad1 == Localidad2) %>%
  group_by(Localidad = Localidad1) %>%
  summarise(
    n_parejas = n(),
    promedio_HS_U_por_locus =
      mean(HS_U_por_locus, na.rm = TRUE),
    max_HS_U_por_locus =
      max(HS_U_por_locus, na.rm = TRUE),
    mejor_margen_HS_U =
      max(margen_HS_U, na.rm = TRUE),
    mejor_margen_HS_U_por_locus =
      max(margen_HS_U_por_locus, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mejor_margen_HS_U_por_locus))

comparacion_espacial_9529 <- resultados_pares_9529 %>%
  group_by(distribucion) %>%
  summarise(
    n_parejas = n(),
    promedio_HS_U_por_locus =
      mean(HS_U_por_locus, na.rm = TRUE),
    mediana_HS_U_por_locus =
      median(HS_U_por_locus, na.rm = TRUE),
    maximo_HS_U_por_locus =
      max(HS_U_por_locus, na.rm = TRUE),
    .groups = "drop"
  )

print(resumen_parentesco_localidad_9529, n = Inf, width = Inf)
print(parejas_mas_cercanas_localidad_9529, n = Inf, width = Inf)
print(resumen_estandarizado_localidad_9529, n = Inf, width = Inf)
print(comparacion_espacial_9529, n = Inf, width = Inf)


# 15. PERMUTATION TEST FOR SPATIAL PATTERNS -----------------------------------

estadistico_observado_9529 <-
  mean(
    resultados_pares_9529$HS_U_por_locus[
      resultados_pares_9529$Localidad1 ==
        resultados_pares_9529$Localidad2
    ],
    na.rm = TRUE
  ) -
  mean(
    resultados_pares_9529$HS_U_por_locus[
      resultados_pares_9529$Localidad1 !=
        resultados_pares_9529$Localidad2
    ],
    na.rm = TRUE
  )

localidades_individuos <- setNames(
  popmap$Localidad,
  popmap$Indiv
)

set.seed(789)

estadisticos_permutados_9529 <- replicate(
  numero_permutaciones,
  {
    localidades_permutadas <- sample(
      localidades_individuos,
      replace = FALSE
    )
    
    names(localidades_permutadas) <-
      names(localidades_individuos)
    
    misma_localidad_permutada <- (
      localidades_permutadas[resultados_pares_9529$Indiv1] ==
        localidades_permutadas[resultados_pares_9529$Indiv2]
    )
    
    mean(
      resultados_pares_9529$HS_U_por_locus[
        misma_localidad_permutada
      ],
      na.rm = TRUE
    ) -
      mean(
        resultados_pares_9529$HS_U_por_locus[
          !misma_localidad_permutada
        ],
        na.rm = TRUE
      )
  }
)

p_permutacion_9529 <- (
  1 +
    sum(
      estadisticos_permutados_9529 >=
        estadistico_observado_9529
    )
) / (
  numero_permutaciones + 1
)

resultado_permutacion_9529 <- tibble(
  estadistico_observado = estadistico_observado_9529,
  n_permutaciones = numero_permutaciones,
  p_unilateral = p_permutacion_9529
)

print(resultado_permutacion_9529, n = Inf, width = Inf)


# 16. SINGLE FULL-PANEL SIMULATION FOR FIGURES --------------------------------

# This simulation illustrates the discriminatory power of the full
# panel of 9,502 polymorphic SNPs. It does not replace the
# missing-data-calibrated simulations used to classify the 703 pairs.

archivo_sim_figuras_9502 <- file.path(
  carpeta_salida,
  "Simulacion_CKMRsim_9502_panel_completo_figuras.rds"
)

if (reutilizar_calculos && file.exists(archivo_sim_figuras_9502)) {
  message("Loading saved full-panel simulation for figures.")
  sim_figuras_9502 <- readRDS(archivo_sim_figuras_9502)
} else {
  message(
    "Running the full-panel simulation of 9,502 SNPs for figures."
  )
  
  set.seed(321)
  
  sim_figuras_9502 <- simulate_Qij(
    C = ckmr_slew_9529,
    calc_relats = relaciones_CKMRsim,
    sim_relats = relaciones_CKMRsim,
    reps = numero_replicas
  )
  
  saveRDS(
    sim_figuras_9502,
    archivo_sim_figuras_9502
  )
}

print(sim_figuras_9502)

extraer_comparaciones_simuladas <- function(simulacion) {
  bind_rows(
    extract_logls(
      simulacion,
      numer = c(FS = 1),
      denom = c(U = 1)
    ) %>%
      mutate(comparacion = "FS/U"),
    extract_logls(
      simulacion,
      numer = c(HS = 1),
      denom = c(U = 1)
    ) %>%
      mutate(comparacion = "HS/U"),
    extract_logls(
      simulacion,
      numer = c(FS = 1),
      denom = c(HS = 1)
    ) %>%
      mutate(comparacion = "FS/HS")
  )
}

simulaciones_logls_9529 <- extraer_comparaciones_simuladas(
  sim_figuras_9502
)

# The shaded band represents the range of classification thresholds
# obtained from the two missing-data calibrations.
referencias_umbrales_9529 <- umbrales_9529 %>%
  group_by(comparacion) %>%
  summarise(
    umbral_min = min(Lambda_star),
    umbral_max = max(Lambda_star),
    .groups = "drop"
  )

maximos_observados_9529 <- tibble(
  comparacion = c("FS/U", "HS/U", "FS/HS"),
  max_observado = c(
    max(resultados_pares_9529$FS_U, na.rm = TRUE),
    max(resultados_pares_9529$HS_U, na.rm = TRUE),
    max(resultados_pares_9529$FS_HS, na.rm = TRUE)
  )
)

referencias_figuras_9529 <- referencias_umbrales_9529 %>%
  left_join(
    maximos_observados_9529,
    by = "comparacion"
  )

print(referencias_figuras_9529, n = Inf, width = Inf)


# 17. FUNCTION TO GENERATE EACH FIGURE PANEL ----------------------------------

colores_relaciones <- c(
  "FS" = "#F8766D",
  "HS" = "#7FBF7B",
  "U" = "#A9D7E5"
)

etiquetas_relaciones <- c(
  "FS" = "Full siblings",
  "HS" = "Half siblings",
  "U" = "Unrelated"
)

crear_panel_CKMRsim <- function(
    comparacion_seleccionada,
    relaciones_mostrar,
    titulo_x
) {
  datos_panel <- simulaciones_logls_9529 %>%
    filter(
      comparacion == comparacion_seleccionada,
      true_relat %in% relaciones_mostrar
    ) %>%
    mutate(
      true_relat = factor(
        true_relat,
        levels = relaciones_mostrar
      )
    )
  
  referencia_panel <- referencias_figuras_9529 %>%
    filter(comparacion == comparacion_seleccionada)
  
  stopifnot(
    nrow(datos_panel) > 0L,
    nrow(referencia_panel) == 1L
  )
  
  umbral_min <- referencia_panel$umbral_min
  umbral_max <- referencia_panel$umbral_max
  max_observado <- referencia_panel$max_observado
  
  # First, estimate the density height to position the annotations.
  grafica_temporal <- ggplot(
    datos_panel,
    aes(
      x = logl_ratio,
      fill = true_relat
    )
  ) +
    geom_density(
      alpha = 0.72,
      color = "black",
      linewidth = 1.0
    )
  
  y_max <- max(
    ggplot_build(grafica_temporal)$data[[1]]$y,
    na.rm = TRUE
  )
  
  cuantiles_x <- quantile(
    datos_panel$logl_ratio,
    probs = c(0.001, 0.999),
    na.rm = TRUE
  )
  
  rango_base <- range(
    cuantiles_x,
    umbral_min,
    umbral_max,
    max_observado,
    na.rm = TRUE
  )
  
  amplitud <- diff(rango_base)
  if (!is.finite(amplitud) || amplitud <= 0) {
    amplitud <- 1
  }
  
  limites_x <- c(
    rango_base[1] - 0.06 * amplitud,
    rango_base[2] + 0.06 * amplitud
  )
  
  desplazamiento_observado <- 0.015 * diff(limites_x)
  punto_medio_umbral <- mean(c(umbral_min, umbral_max))
  
  ggplot(
    datos_panel,
    aes(
      x = logl_ratio,
      fill = true_relat
    )
  ) +
    annotate(
      "rect",
      xmin = umbral_min,
      xmax = umbral_max,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey70",
      alpha = 0.30
    ) +
    geom_density(
      alpha = 0.72,
      color = "black",
      linewidth = 1.0
    ) +
    geom_vline(
      xintercept = c(umbral_min, umbral_max),
      color = "grey40",
      linetype = "dashed",
      linewidth = 1.1
    ) +
    geom_vline(
      xintercept = max_observado,
      color = "#8A1538",
      linetype = "longdash",
      linewidth = 1.3
    ) +
    annotate(
      "text",
      x = max_observado + desplazamiento_observado,
      y = y_max * 1.085,
      label = paste0(
        "Observed max\n",
        round(max_observado, 1)
      ),
      color = "#8A1538",
      size = 6.2,
      fontface = "bold",
      hjust = 0,
      lineheight = 0.90
    ) +
    annotate(
      "text",
      x = punto_medio_umbral,
      y = y_max * 1.085,
      label = paste0(
        "Classification thresholds\n",
        round(umbral_min, 1),
        "-",
        round(umbral_max, 1)
      ),
      color = "grey25",
      size = 6.2,
      fontface = "bold",
      hjust = 0.5,
      lineheight = 0.90
    ) +
    scale_fill_manual(
      values = colores_relaciones,
      breaks = relaciones_mostrar,
      labels = etiquetas_relaciones[relaciones_mostrar],
      drop = TRUE
    ) +
    labs(
      x = titulo_x,
      y = "Density",
      fill = NULL
    ) +
    coord_cartesian(
      xlim = limites_x,
      ylim = c(0, y_max * 1.18),
      clip = "off"
    ) +
    theme_classic(base_size = 22) +
    theme(
      text = element_text(color = "black"),
      axis.title = element_text(
        size = 25,
        color = "black",
        face = "bold"
      ),
      axis.text = element_text(
        size = 20,
        color = "black"
      ),
      legend.text = element_text(
        size = 20,
        color = "black"
      ),
      legend.key.size = grid::unit(0.80, "cm"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box.spacing = grid::unit(0, "pt"),
      legend.spacing.x = grid::unit(0.25, "cm"),
      legend.margin = margin(t = 7, r = 0, b = 0, l = 0),
      plot.margin = margin(t = 18, r = 14, b = 10, l = 12),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}


# 18. GENERATE AND SAVE FIGURES -----------------------------------------------

grafica_FS_U_final <- crear_panel_CKMRsim(
  comparacion_seleccionada = "FS/U",
  relaciones_mostrar = c("FS", "U"),
  titulo_x = "Log-likelihood ratio (FS/U)"
)

grafica_HS_U_final <- crear_panel_CKMRsim(
  comparacion_seleccionada = "HS/U",
  relaciones_mostrar = c("HS", "U"),
  titulo_x = "Log-likelihood ratio (HS/U)"
)

grafica_FS_HS_final <- crear_panel_CKMRsim(
  comparacion_seleccionada = "FS/HS",
  relaciones_mostrar = c("FS", "HS"),
  titulo_x = "Log-likelihood ratio (FS/HS)"
)



print(grafica_FS_U_final)
print(grafica_HS_U_final)
print(grafica_FS_HS_final)
print(figura_CKMRsim_9529)

ggsave(
  filename = file.path(
    carpeta_salida,
    "CKMRsim_FS_vs_U_9529.png"
  ),
  plot = grafica_FS_U_final,
  width = 16,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    carpeta_salida,
    "CKMRsim_HS_vs_U_9529.png"
  ),
  plot = grafica_HS_U_final,
  width = 16,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    carpeta_salida,
    "CKMRsim_FS_vs_HS_9529.png"
  ),
  plot = grafica_FS_HS_final,
  width = 16,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)



# 19. EXPORT FINAL TABLES AND OBJECTS -----------------------------------------

write.csv(
  missing_individuos_9529,
  file.path(
    carpeta_salida,
    "QC_missing_individuos_panel_9529.csv"
  ),
  row.names = FALSE
)

write.csv(
  resumen_grupos_missing_9529,
  file.path(
    carpeta_salida,
    "QC_loci_comparables_por_grupo_panel_9529.csv"
  ),
  row.names = FALSE
)

write.csv(
  umbrales_9529,
  file.path(
    carpeta_salida,
    "Umbrales_CKMRsim_panel_9529.csv"
  ),
  row.names = FALSE
)

write.csv(
  resultados_pares_9529,
  file.path(
    carpeta_salida,
    "CKMRsim_panel_9529_todas_las_parejas.csv"
  ),
  row.names = FALSE
)

write.csv(
  pares_emparentados_9529,
  file.path(
    carpeta_salida,
    "CKMRsim_panel_9529_pares_emparentados.csv"
  ),
  row.names = FALSE
)

write.csv(
  resumen_final_9529,
  file.path(
    carpeta_salida,
    "Resumen_CKMRsim_panel_9529.csv"
  ),
  row.names = FALSE
)

write.csv(
  resumen_parentesco_localidad_9529,
  file.path(
    carpeta_salida,
    "Resumen_CKMRsim_por_localidad_9529.csv"
  ),
  row.names = FALSE
)

write.csv(
  comparacion_espacial_9529,
  file.path(
    carpeta_salida,
    "Comparacion_dentro_entre_localidades_9529.csv"
  ),
  row.names = FALSE
)

write.csv(
  resumen_estandarizado_localidad_9529,
  file.path(
    carpeta_salida,
    "Resumen_estandarizado_por_localidad_9529.csv"
  ),
  row.names = FALSE
)

write.csv(
  parejas_mas_cercanas_localidad_9529,
  file.path(
    carpeta_salida,
    "Parejas_mas_cercanas_por_localidad_9529.csv"
  ),
  row.names = FALSE
)

write.csv(
  duplicados_9529,
  file.path(
    carpeta_salida,
    "Revision_posibles_duplicados_9529.csv"
  ),
  row.names = FALSE
)

write.csv(
  resultado_permutacion_9529,
  file.path(
    carpeta_salida,
    "Prueba_permutacion_espacial_9529.csv"
  ),
  row.names = FALSE
)

saveRDS(
  estadisticos_permutados_9529,
  file.path(
    carpeta_salida,
    "Distribucion_permutada_espacial_9529.rds"
  )
)

write.csv(
  referencias_figuras_9529,
  file.path(
    carpeta_salida,
    "Referencias_figuras_CKMRsim_9529.csv"
  ),
  row.names = FALSE
)

saveRDS(
  resultados_pares_9529,
  file.path(
    carpeta_salida,
    "Resultados_parejas_CKMRsim_panel_9529.rds"
  )
)

saveRDS(
  list(
    parametros = list(
      numero_snp_entrada = numero_snp_entrada,
      numero_snp_polimorficos =
        numero_snp_polimorficos_esperado,
      numero_individuos = numero_individuos_esperado,
      umbral_missing_individual =
        umbral_missing_individual,
      epsilon = epsilon_genotipificacion,
      replicas = numero_replicas,
      FNR_objetivo = FNR_objetivo,
      permutaciones = numero_permutaciones
    ),
    resumen_missing = resumen_grupos_missing_9529,
    umbrales = umbrales_9529,
    resultados_pares = resultados_pares_9529,
    resumen_final = resumen_final_9529,
    comparacion_espacial = comparacion_espacial_9529,
    prueba_permutacion = resultado_permutacion_9529,
    referencias_figuras = referencias_figuras_9529
  ),
  file.path(
    carpeta_salida,
    "Resultados_integrados_CKMRsim_9529.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(carpeta_salida, "sessionInfo_CKMRsim_9529.txt")
)

message("Analysis completed successfully.")
message("Results saved to: ", carpeta_salida)

# END OF SCRIPT
