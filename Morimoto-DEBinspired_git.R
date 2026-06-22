
## =========================================================
## 0. Packages and user settings
## =========================================================

packages_needed <- c(
  "tidyverse",
  "mgcv",
  "glmmTMB",
  "coxme",
  "survival",
  "splines",
  "performance",
  "DHARMa",
  "patchwork",
  "viridis",
  "scales",
  "broom",
  "broom.mixed",
  "ggforce"
)

packages_to_install <- packages_needed[
  !(packages_needed %in% installed.packages()[, "Package"])
]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

library(tidyverse)
library(mgcv)
library(glmmTMB)
library(coxme)
library(survival)
library(splines)
library(performance)
library(DHARMa)
library(patchwork)
library(viridis)
library(scales)
library(broom)
library(broom.mixed)
library(ggforce)


input_file <- "Dataset_Lee_Nutrigonometry.csv"
output_dir <- "DEB-GFN-Drosophila-ouput"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


## ---------------------------------------------------------
## Biological constants


eP <- 17                 # J / mg protein
eC <- 17                 # J / mg carbohydrate
biomass_energy <- 22     # J / mg dry biomass


## ---------------------------------------------------------
## Baseline biological assumptions


egg_wet_mass_ug <- 9
egg_dry_fraction <- 0.25
egg_dry_mass_ug <- egg_wet_mass_ug * egg_dry_fraction
adult_dry_mass_ug <- 300

eta_P <- 0.75
eta_C <- 0.75
kap_R <- 0.85
maint_frac_adult_energy_per_day <- 0.10

egg_protein_fraction_dry <- 0.50
eta_protein_to_egg <- 0.85



## ---------------------------------------------------------
## Nutritional target
## ---------------------------------------------------------
## Lee et al. Drosophila target: P:C = 1:4.
## Target point uses mean total intake of observed 1:4 treatments.
## ---------------------------------------------------------

target_ratio_P <- 1
target_ratio_C <- 4

target_fraction_P <- target_ratio_P / (target_ratio_P + target_ratio_C)
target_fraction_C <- target_ratio_C / (target_ratio_P + target_ratio_C)

distance_method <- "euclidean"

## TRUE if dailyeggs was calculated as lifetimeegg/lifespan.
## Then reproductive burden models are descriptive/sensitivity only.
dailyeggs_is_lifetime_over_lifespan <- TRUE


## =========================================================
## 1. Read and clean data
## =========================================================

dat_raw <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  na.strings = c(".", "NA", "")
)

dat <- dat_raw %>%
  mutate(
    carb_eaten    = as.numeric(carb_eaten),
    protein_eaten = as.numeric(protein_eaten),
    lifespan      = as.numeric(lifespan),
    lifetimeegg   = as.numeric(lifetimeegg),
    dailyeggs     = as.numeric(dailyeggs),
    treatment     = as.factor(treatment),
    Ratio         = as.factor(Ratio),
    Food          = as.factor(Food),
    P_fixed       = as.numeric(P_fixed),
    C_fixed       = as.numeric(C_fixed),
    total_eaten   = carb_eaten + protein_eaten,
    PC_ratio      = ifelse(carb_eaten > 0, protein_eaten / carb_eaten, NA_real_),
    fly_id        = row_number()
  ) %>%
  filter(
    is.finite(carb_eaten),
    is.finite(protein_eaten),
    is.finite(lifespan),
    is.finite(lifetimeegg),
    is.finite(dailyeggs),
    lifespan > 0,
    carb_eaten >= 0,
    protein_eaten >= 0,
    dailyeggs >= 0
  )

data_cleaning_summary <- tibble(
  n_rows_raw = nrow(dat_raw),
  n_rows_clean = nrow(dat),
  n_treatments = n_distinct(dat$treatment),
  n_ratio_levels = n_distinct(dat$Ratio),
  mean_lifespan = mean(dat$lifespan, na.rm = TRUE),
  mean_dailyeggs = mean(dat$dailyeggs, na.rm = TRUE),
  mean_protein_ug_day = mean(dat$protein_eaten, na.rm = TRUE),
  mean_carb_ug_day = mean(dat$carb_eaten, na.rm = TRUE)
)


## =========================================================
## 2. Treatment summary and target point
## =========================================================

diet_summary <- dat %>%
  group_by(treatment, Ratio, Food, P_fixed, C_fixed) %>%
  summarise(
    n = n(),
    P_daily_ug = mean(protein_eaten, na.rm = TRUE),
    C_daily_ug = mean(carb_eaten, na.rm = TRUE),
    total_daily_ug = mean(total_eaten, na.rm = TRUE),
    PC_daily_ratio = ifelse(C_daily_ug > 0, P_daily_ug / C_daily_ug, NA_real_),
    lifespan_obs = mean(lifespan, na.rm = TRUE),
    lifespan_sd  = sd(lifespan, na.rm = TRUE),
    lifespan_se  = lifespan_sd / sqrt(n),
    daily_eggs_obs = mean(dailyeggs, na.rm = TRUE),
    daily_eggs_sd  = sd(dailyeggs, na.rm = TRUE),
    daily_eggs_se  = daily_eggs_sd / sqrt(n),
    lifetime_eggs_obs = mean(lifetimeegg, na.rm = TRUE),
    lifetime_eggs_sd  = sd(lifetimeegg, na.rm = TRUE),
    lifetime_eggs_se  = lifetime_eggs_sd / sqrt(n),
    observed_hazard = 1 / lifespan_obs,
    .groups = "drop"
  )

target_total_intake <- diet_summary %>%
  filter(Ratio == "(1:4)") %>%
  summarise(target_total = mean(total_daily_ug, na.rm = TRUE)) %>%
  pull(target_total)

if (length(target_total_intake) == 0 || is.na(target_total_intake)) {
  stop("Could not identify Ratio == '(1:4)'. Check Ratio coding.")
}

target_P_ug <- target_fraction_P * target_total_intake
target_C_ug <- target_fraction_C * target_total_intake

target_point <- tibble(
  target_definition = "mean total intake of observed 1:4 treatments",
  target_total_intake = target_total_intake,
  target_P_ug = target_P_ug,
  target_C_ug = target_C_ug,
  target_PC_ratio = target_P_ug / target_C_ug
)


## =========================================================
## 3. DEB-inspired variables
## =========================================================

compute_nutritional_distance <- function(
    P,
    C,
    target_P,
    target_C,
    method = "euclidean",
    scale_by = NULL
) {
  if (method == "euclidean") {
    d <- sqrt((P - target_P)^2 + (C - target_C)^2)
  } else {
    stop("Only Euclidean distance is implemented in this clean script.")
  }
  
  if (!is.null(scale_by)) {
    d <- d / scale_by
  }
  
  d
}


calculate_individual_deb <- function(
    data,
    target_P,
    target_C,
    target_total,
    distance_method = "euclidean",
    eP = 17,
    eC = 17,
    biomass_energy = 22,
    egg_dry_mass_ug = 0.50,
    adult_dry_mass_ug = 50,
    eta_P = 0.75,
    eta_C = 0.75,
    kap_R = 0.85,
    maint_frac_adult_energy_per_day = 0.10,
    egg_protein_fraction_dry = 0.50,
    eta_protein_to_egg = 0.85
) {
  
  data %>%
    mutate(
      P_daily_ug = protein_eaten,
      C_daily_ug = carb_eaten,
      total_daily_ug = P_daily_ug + C_daily_ug,
      
      distance_from_target =
        compute_nutritional_distance(
          P = P_daily_ug,
          C = C_daily_ug,
          target_P = target_P,
          target_C = target_C,
          method = distance_method,
          scale_by = NULL
        ),
      
      distance_from_target_scaled =
        compute_nutritional_distance(
          P = P_daily_ug,
          C = C_daily_ug,
          target_P = target_P,
          target_C = target_C,
          method = distance_method,
          scale_by = target_total
        ),
      
      distance_from_target_squared = distance_from_target_scaled^2,
      
      gross_weighted_intake_J_day =
        0.001 * (
          eta_P * eP * P_daily_ug +
            eta_C * eC * C_daily_ug
        ),
      
      assimilated_energy_J_day = gross_weighted_intake_J_day,
      
      egg_energy_J = egg_dry_mass_ug * 0.001 * biomass_energy,
      adult_body_energy_J = adult_dry_mass_ug * 0.001 * biomass_energy,
      
      maintenance_J_day =
        maint_frac_adult_energy_per_day * adult_body_energy_J,
      
      reproductive_cost_J_day =
        dailyeggs * egg_energy_J / kap_R,
      
      egg_protein_ug =
        egg_dry_mass_ug * egg_protein_fraction_dry,
      
      assimilated_protein_ug_day =
        eta_P * P_daily_ug,
      
      reproductive_protein_required_ug_day =
        dailyeggs * egg_protein_ug / eta_protein_to_egg,
      
      reproductive_protein_burden =
        ifelse(
          assimilated_protein_ug_day > 0,
          reproductive_protein_required_ug_day / assimilated_protein_ug_day,
          NA_real_
        ),
      
      reproductive_protein_safety_margin_ug_day =
        assimilated_protein_ug_day - reproductive_protein_required_ug_day,
      
      reproductive_protein_safety_margin_fraction =
        ifelse(
          assimilated_protein_ug_day > 0,
          reproductive_protein_safety_margin_ug_day / assimilated_protein_ug_day,
          NA_real_
        ),
      
      protein_feasible_for_reproduction =
        ifelse(
          assimilated_protein_ug_day > 0,
          reproductive_protein_safety_margin_ug_day >= 0,
          NA
        ),
      
      reproductive_energetic_burden =
        reproductive_cost_J_day / assimilated_energy_J_day,
      
      non_reproductive_fraction_of_assimilated =
        1 - reproductive_energetic_burden,
      
      reproductive_share_of_explicit_costs =
        reproductive_cost_J_day /
        (reproductive_cost_J_day + maintenance_J_day),
      
      somatic_share_of_explicit_costs =
        maintenance_J_day /
        (reproductive_cost_J_day + maintenance_J_day),
      
      energetic_safety_margin_J_day =
        assimilated_energy_J_day -
        reproductive_cost_J_day -
        maintenance_J_day,
      
      energetic_safety_margin_fraction =
        energetic_safety_margin_J_day / assimilated_energy_J_day,
      
      energetically_feasible =
        energetic_safety_margin_J_day >= 0,
      
      explicit_cost_fraction_of_assimilation =
        (reproductive_cost_J_day + maintenance_J_day) /
        assimilated_energy_J_day,
      
      death_event = 1,
      lifespan_day = pmax(ceiling(lifespan), 1),
      observed_hazard_proxy = 1 / lifespan,
      log_lifespan = log(lifespan),
      log_hazard_proxy = log(observed_hazard_proxy),
      
      treatment = factor(treatment),
      Ratio = factor(Ratio),
      Food = factor(Food)
    ) %>%
    filter(
      is.finite(lifespan),
      is.finite(lifespan_day),
      lifespan_day > 0,
      is.finite(assimilated_energy_J_day),
      assimilated_energy_J_day > 0,
      is.finite(reproductive_energetic_burden),
      is.finite(energetic_safety_margin_fraction),
      is.finite(distance_from_target_scaled),
      is.finite(reproductive_share_of_explicit_costs)
    )
}


deb_individual <- calculate_individual_deb(
  data = dat,
  target_P = target_P_ug,
  target_C = target_C_ug,
  target_total = target_total_intake,
  distance_method = distance_method,
  eP = eP,
  eC = eC,
  biomass_energy = biomass_energy,
  egg_dry_mass_ug = egg_dry_mass_ug,
  adult_dry_mass_ug = adult_dry_mass_ug,
  eta_P = eta_P,
  eta_C = eta_C,
  kap_R = kap_R,
  maint_frac_adult_energy_per_day = maint_frac_adult_energy_per_day,
  egg_protein_fraction_dry = egg_protein_fraction_dry,
  eta_protein_to_egg = eta_protein_to_egg
)


parameter_energy_check <- tibble(
  egg_wet_mass_ug = egg_wet_mass_ug,
  egg_dry_fraction = egg_dry_fraction,
  egg_dry_mass_ug = egg_dry_mass_ug,
  adult_dry_mass_ug = adult_dry_mass_ug,
  egg_energy_J = egg_dry_mass_ug * 0.001 * biomass_energy,
  adult_body_energy_J = adult_dry_mass_ug * 0.001 * biomass_energy,
  maintenance_J_day =
    maint_frac_adult_energy_per_day *
    adult_dry_mass_ug * 0.001 * biomass_energy,
  eggs_equivalent_to_one_adult_body =
    (adult_dry_mass_ug * 0.001 * biomass_energy) /
    (egg_dry_mass_ug * 0.001 * biomass_energy),
  maintenance_in_egg_equivalents_per_day =
    (
      maint_frac_adult_energy_per_day *
        adult_dry_mass_ug * 0.001 * biomass_energy
    ) /
    (egg_dry_mass_ug * 0.001 * biomass_energy),
  egg_protein_fraction_dry = egg_protein_fraction_dry,
  eta_protein_to_egg = eta_protein_to_egg
)


energetic_plausibility_check <- deb_individual %>%
  summarise(
    n = n(),
    mean_assimilated_energy_J_day = mean(assimilated_energy_J_day, na.rm = TRUE),
    mean_reproductive_cost_J_day = mean(reproductive_cost_J_day, na.rm = TRUE),
    mean_maintenance_J_day = mean(maintenance_J_day, na.rm = TRUE),
    median_phi_R = median(reproductive_energetic_burden, na.rm = TRUE),
    median_rho_R = median(reproductive_share_of_explicit_costs, na.rm = TRUE),
    median_safety_margin_fraction = median(energetic_safety_margin_fraction, na.rm = TRUE),
    feasible_fraction = mean(energetically_feasible, na.rm = TRUE),
    fraction_reproduction_exceeds_assimilation =
      mean(reproductive_energetic_burden > 1, na.rm = TRUE),
    fraction_explicit_costs_exceed_assimilation =
      mean(explicit_cost_fraction_of_assimilation > 1, na.rm = TRUE)
  )


protein_reproductive_plausibility_check <- deb_individual %>%
  summarise(
    n = n(),
    mean_assimilated_protein_ug_day =
      mean(assimilated_protein_ug_day, na.rm = TRUE),
    mean_reproductive_protein_required_ug_day =
      mean(reproductive_protein_required_ug_day, na.rm = TRUE),
    median_reproductive_protein_burden =
      median(reproductive_protein_burden, na.rm = TRUE),
    median_reproductive_protein_safety_margin_fraction =
      median(reproductive_protein_safety_margin_fraction, na.rm = TRUE),
    protein_feasible_fraction =
      mean(protein_feasible_for_reproduction, na.rm = TRUE),
    fraction_reproductive_protein_requirement_exceeds_assimilated_protein =
      mean(reproductive_protein_burden > 1, na.rm = TRUE)
  )


deb_summary_table <- deb_individual %>%
  summarise(
    n = n(),
    mean_lifespan = mean(lifespan, na.rm = TRUE),
    sd_lifespan = sd(lifespan, na.rm = TRUE),
    mean_P_daily_ug = mean(P_daily_ug, na.rm = TRUE),
    mean_C_daily_ug = mean(C_daily_ug, na.rm = TRUE),
    mean_daily_eggs = mean(dailyeggs, na.rm = TRUE),
    mean_lifetime_eggs = mean(lifetimeegg, na.rm = TRUE),
    mean_assimilated_energy_J_day = mean(assimilated_energy_J_day, na.rm = TRUE),
    mean_reproductive_cost_J_day = mean(reproductive_cost_J_day, na.rm = TRUE),
    mean_maintenance_J_day = mean(maintenance_J_day, na.rm = TRUE),
    mean_reproductive_energetic_burden =
      mean(reproductive_energetic_burden, na.rm = TRUE),
    mean_reproductive_share_explicit_costs =
      mean(reproductive_share_of_explicit_costs, na.rm = TRUE),
    mean_safety_margin_fraction =
      mean(energetic_safety_margin_fraction, na.rm = TRUE),
    mean_distance_scaled =
      mean(distance_from_target_scaled, na.rm = TRUE),
    mean_assimilated_protein_ug_day =
      mean(assimilated_protein_ug_day, na.rm = TRUE),
    mean_reproductive_protein_required_ug_day =
      mean(reproductive_protein_required_ug_day, na.rm = TRUE),
    mean_reproductive_protein_burden =
      mean(reproductive_protein_burden, na.rm = TRUE),
    mean_reproductive_protein_safety_margin_fraction =
      mean(reproductive_protein_safety_margin_fraction, na.rm = TRUE),
    protein_feasible_fraction =
      mean(protein_feasible_for_reproduction, na.rm = TRUE)
  )


## =========================================================
## 4. Core statistical models
## =========================================================

deb_treatment <- deb_individual %>%
  group_by(treatment, Ratio, Food, P_fixed, C_fixed) %>%
  summarise(
    n = n(),
    P_daily_ug = mean(P_daily_ug, na.rm = TRUE),
    C_daily_ug = mean(C_daily_ug, na.rm = TRUE),
    distance_from_target_scaled = mean(distance_from_target_scaled, na.rm = TRUE),
    reproductive_energetic_burden = mean(reproductive_energetic_burden, na.rm = TRUE),
    reproductive_protein_burden = median(reproductive_protein_burden, na.rm = TRUE),
    reproductive_protein_safety_margin_fraction =
      median(reproductive_protein_safety_margin_fraction, na.rm = TRUE),
    protein_feasible_fraction = mean(protein_feasible_for_reproduction, na.rm = TRUE),
    reproductive_share_of_explicit_costs =
      mean(reproductive_share_of_explicit_costs, na.rm = TRUE),
    energetic_safety_margin_fraction =
      mean(energetic_safety_margin_fraction, na.rm = TRUE),
    lifespan_obs = mean(lifespan, na.rm = TRUE),
    lifetime_eggs_obs = mean(lifetimeegg, na.rm = TRUE),
    log_hazard = log(1 / lifespan_obs),
    .groups = "drop"
  )



## ---------------------------------------------------------
## Discrete-time survival dataset
## ---------------------------------------------------------

surv_daily <- deb_individual %>%
  select(
    fly_id,
    treatment,
    Ratio,
    Food,
    lifespan,
    lifespan_day,
    P_daily_ug,
    C_daily_ug,
    total_daily_ug,
    energetic_safety_margin_fraction,
    distance_from_target_scaled,
    reproductive_share_of_explicit_costs,
    reproductive_protein_burden
  ) %>%
  rowwise() %>%
  mutate(day = list(seq_len(lifespan_day))) %>%
  unnest(day) %>%
  ungroup() %>%
  mutate(
    death = ifelse(day == lifespan_day, 1, 0),
    treatment = droplevels(factor(treatment)),
    Ratio = droplevels(factor(Ratio)),
    Food = droplevels(factor(Food))
  )


safe_z <- function(x, center, scale) {
  if (!is.finite(scale) || scale == 0) {
    return(rep(0, length(x)))
  }
  (x - center) / scale
}


make_scale_ref <- function(data) {
  data %>%
    summarise(
      P_mean = mean(P_daily_ug, na.rm = TRUE),
      P_sd   = sd(P_daily_ug, na.rm = TRUE),
      C_mean = mean(C_daily_ug, na.rm = TRUE),
      C_sd   = sd(C_daily_ug, na.rm = TRUE),
      margin_mean = mean(energetic_safety_margin_fraction, na.rm = TRUE),
      margin_sd   = sd(energetic_safety_margin_fraction, na.rm = TRUE),
      dist_mean = mean(distance_from_target_scaled, na.rm = TRUE),
      dist_sd   = sd(distance_from_target_scaled, na.rm = TRUE),
      rho_mean = mean(reproductive_share_of_explicit_costs, na.rm = TRUE),
      rho_sd   = sd(reproductive_share_of_explicit_costs, na.rm = TRUE),
      protein_burden_mean = mean(reproductive_protein_burden, na.rm = TRUE),
      protein_burden_sd   = sd(reproductive_protein_burden, na.rm = TRUE)
    )
}


add_scaled_predictors <- function(data, scale_ref) {
  data %>%
    mutate(
      P_daily_ug_z =
        safe_z(P_daily_ug, scale_ref$P_mean, scale_ref$P_sd),
      C_daily_ug_z =
        safe_z(C_daily_ug, scale_ref$C_mean, scale_ref$C_sd),
      energetic_safety_margin_fraction_z =
        safe_z(
          energetic_safety_margin_fraction,
          scale_ref$margin_mean,
          scale_ref$margin_sd
        ),
      distance_from_target_scaled_z =
        safe_z(
          distance_from_target_scaled,
          scale_ref$dist_mean,
          scale_ref$dist_sd
        ),
      reproductive_share_of_explicit_costs_z =
        safe_z(
          reproductive_share_of_explicit_costs,
          scale_ref$rho_mean,
          scale_ref$rho_sd
        ),
      reproductive_protein_burden_z =
        safe_z(
          reproductive_protein_burden,
          scale_ref$protein_burden_mean,
          scale_ref$protein_burden_sd
        )
    )
}

scale_ref <- make_scale_ref(surv_daily)

surv_daily <- surv_daily %>%
  add_scaled_predictors(scale_ref)


## ---------------------------------------------------------
## Final time-varying model
## ---------------------------------------------------------
## This is the primary inferential model if dailyeggs is derived
## from lifetimeegg/lifespan.
## ---------------------------------------------------------

dt7_timevarying <- glmmTMB(
  death ~
    Ratio * ns(day, df = 5) +
    Food * ns(day, df = 5) +
    P_daily_ug_z +
    C_daily_ug_z +
    energetic_safety_margin_fraction_z +
    distance_from_target_scaled_z +
    P_daily_ug_z:ns(day, df = 3) +
    C_daily_ug_z:ns(day, df = 3) +
    distance_from_target_scaled_z:ns(day, df = 3) +
    (1 | treatment),
  data = surv_daily,
  family = binomial(link = "cloglog")
)
Anova(dt7_timevarying)


## Reproductive-burden sensitivity model.
dt9_timevarying <- glmmTMB(
  death ~
    Ratio * ns(day, df = 5) +
    Food * ns(day, df = 5) +
    P_daily_ug_z +
    C_daily_ug_z +
    energetic_safety_margin_fraction_z +
    distance_from_target_scaled_z +
    reproductive_share_of_explicit_costs_z +
    P_daily_ug_z:ns(day, df = 3) +
    C_daily_ug_z:ns(day, df = 3) +
    distance_from_target_scaled_z:ns(day, df = 3) +
    reproductive_share_of_explicit_costs_z:ns(day, df = 3) +
    (1 | treatment),
  data = surv_daily,
  family = binomial(link = "cloglog")
)
summary(dt9_timevarying)
performance::check_collinearity(dt9_timevarying)



## Reproductive-burden sensitivity model.
dt8_timevarying <- glmmTMB(
  death ~
    Ratio * ns(day, df = 5) +
    Food * ns(day, df = 5) +
    P_daily_ug_z +
    C_daily_ug_z +
    distance_from_target_scaled_z +
    P_daily_ug_z:ns(day, df = 3) +
    C_daily_ug_z:ns(day, df = 3) +
    distance_from_target_scaled_z:ns(day, df = 3) +
    (1 | treatment),
  data = surv_daily,
  family = binomial(link = "cloglog")
)
summary(dt8_timevarying)

dt_model_table <- AIC(
  dt7_timevarying,
  dt8_timevarying,
  dt9_timevarying
) 
anova(dt7_timevarying, dt8_timevarying, dt9_timevarying)

## Reproductive-burden sensitivity model WITHOUT interaction terms for VIF assessment
dt9_no_interactions <- glmmTMB(
  death ~
    Ratio +
    Food +
    ns(day, df = 5) +
    P_daily_ug_z +
    C_daily_ug_z +
    energetic_safety_margin_fraction_z +
    distance_from_target_scaled_z +
    reproductive_share_of_explicit_costs_z +
    (1 | treatment),
  data = surv_daily,
  family = binomial(link = "cloglog")
)

performance::check_collinearity(dt9_no_interactions)



final_dt_model <- dt9_timevarying
final_dt_model_name <- "dt9_timevarying"


final_model_tidy <- broom.mixed::tidy(
  final_dt_model,
  effects = "fixed"
)


## =========================================================
## 5. Prediction and calibration objects
## =========================================================

surv_daily <- surv_daily %>%
  mutate(
    pred_death_prob_final =
      predict(
        final_dt_model,
        newdata = surv_daily,
        type = "response",
        re.form = NULL
      )
  )

calibration_by_day <- surv_daily %>%
  group_by(day) %>%
  summarise(
    n_at_risk = n(),
    observed_death_prob = mean(death, na.rm = TRUE),
    predicted_death_prob = mean(pred_death_prob_final, na.rm = TRUE),
    observed_deaths = sum(death, na.rm = TRUE),
    predicted_deaths = sum(pred_death_prob_final, na.rm = TRUE),
    .groups = "drop"
  )

calibration_by_risk <- surv_daily %>%
  mutate(risk_decile = ntile(pred_death_prob_final, 10)) %>%
  group_by(risk_decile) %>%
  summarise(
    n = n(),
    observed_deaths = sum(death, na.rm = TRUE),
    expected_deaths = sum(pred_death_prob_final, na.rm = TRUE),
    observed_rate = mean(death, na.rm = TRUE),
    predicted_rate = mean(pred_death_prob_final, na.rm = TRUE),
    .groups = "drop"
  )


max_day <- max(deb_individual$lifespan_day, na.rm = TRUE)

pred_ind_daily <- deb_individual %>%
  select(
    fly_id,
    treatment,
    Ratio,
    Food,
    lifespan,
    lifespan_day,
    P_daily_ug,
    C_daily_ug,
    total_daily_ug,
    energetic_safety_margin_fraction,
    distance_from_target_scaled,
    reproductive_share_of_explicit_costs,
    reproductive_protein_burden
  ) %>%
  tidyr::crossing(day = 1:max_day) %>%
  mutate(
    treatment = factor(treatment, levels = levels(surv_daily$treatment)),
    Ratio = factor(Ratio, levels = levels(surv_daily$Ratio)),
    Food = factor(Food, levels = levels(surv_daily$Food))
  ) %>%
  add_scaled_predictors(scale_ref) %>%
  mutate(
    pred_death_prob =
      predict(
        final_dt_model,
        newdata = .,
        type = "response",
        re.form = NULL
      )
  ) %>%
  group_by(fly_id) %>%
  arrange(day, .by_group = TRUE) %>%
  mutate(
    survival_end_day = cumprod(1 - pred_death_prob),
    survival_start_day = lag(survival_end_day, default = 1)
  ) %>%
  ungroup()

pred_ind_lifespan <- pred_ind_daily %>%
  group_by(fly_id) %>%
  summarise(
    predicted_lifespan = sum(survival_start_day, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    deb_individual %>% select(fly_id, lifespan, treatment, Ratio, Food),
    by = "fly_id"
  ) %>%
  mutate(residual_lifespan = lifespan - predicted_lifespan)


pred_treatment_lifespan <- pred_ind_lifespan %>%
  group_by(treatment, Ratio, Food) %>%
  summarise(
    n = n(),
    observed_mean_lifespan = mean(lifespan, na.rm = TRUE),
    predicted_mean_lifespan = mean(predicted_lifespan, na.rm = TRUE),
    residual_mean_lifespan = observed_mean_lifespan - predicted_mean_lifespan,
    .groups = "drop"
  )


profile_data <- deb_individual %>%
  group_by(treatment, Ratio, Food) %>%
  summarise(
    P_daily_ug = mean(P_daily_ug, na.rm = TRUE),
    C_daily_ug = mean(C_daily_ug, na.rm = TRUE),
    total_daily_ug = mean(total_daily_ug, na.rm = TRUE),
    energetic_safety_margin_fraction =
      median(energetic_safety_margin_fraction, na.rm = TRUE),
    distance_from_target_scaled =
      median(distance_from_target_scaled, na.rm = TRUE),
    reproductive_share_of_explicit_costs =
      median(reproductive_share_of_explicit_costs, na.rm = TRUE),
    reproductive_protein_burden =
      median(reproductive_protein_burden, na.rm = TRUE),
    lifespan = mean(lifespan, na.rm = TRUE),
    .groups = "drop"
  )

profiles <- bind_rows(
  profile_data %>%
    slice_min(distance_from_target_scaled, n = 1, with_ties = FALSE) %>%
    mutate(profile = "Near P:C 1:4 target"),
  profile_data %>%
    slice_max(C_daily_ug, n = 1, with_ties = FALSE) %>%
    mutate(profile = "High carbohydrate"),
  profile_data %>%
    slice_max(P_daily_ug, n = 1, with_ties = FALSE) %>%
    mutate(profile = "High protein"),
  profile_data %>%
    slice_min(total_daily_ug, n = 1, with_ties = FALSE) %>%
    mutate(profile = "Low intake"),
  profile_data %>%
    slice_max(distance_from_target_scaled, n = 1, with_ties = FALSE) %>%
    mutate(profile = "Farthest from target"),
  profile_data %>%
    slice_max(reproductive_share_of_explicit_costs, n = 1, with_ties = FALSE) %>%
    mutate(profile = "High average reproductive burden")
) %>%
  select(
    profile,
    treatment,
    P_daily_ug,
    C_daily_ug,
    total_daily_ug,
    energetic_safety_margin_fraction,
    distance_from_target_scaled,
    reproductive_share_of_explicit_costs,
    reproductive_protein_burden
  ) %>%
  add_scaled_predictors(scale_ref)

pred_days <- tibble(day = 1:max_day)

pred_surv_profiles <- tidyr::crossing(profiles, pred_days) %>%
  mutate(treatment = factor(treatment, levels = levels(surv_daily$treatment))) %>%
  mutate(
    death_prob = predict(
      final_dt_model,
      newdata = .,
      type = "response",
      re.form = NA
    )
  ) %>%
  group_by(profile) %>%
  arrange(day, .by_group = TRUE) %>%
  mutate(survival_prob = cumprod(1 - death_prob)) %>%
  ungroup()







## ---------------------------------------------------------
## Supplementary exploratory lipid-adjusted energetic scenario
## Lee et al. lipid data available only for:
## P:C = 1:16, 1:4, 1:2 at Food = 180
## ---------------------------------------------------------

## Lipid data from Lee et al.
## Values are mean lipid fraction of dry mass ± SE.
lee_lipid_data <- tibble::tribble(
  ~Ratio,   ~Food_num, ~lipid_fraction_dry_mass, ~lipid_se,
  "(1:16)", 180,       0.300,                     0.011,
  "(1:4)",  180,       0.151,                     0.007,
  "(1:2)",  180,       0.140,                     0.009
)

## Energy-density assumptions for lipid-adjusted body energy.
## Lipid energy density is higher than lean biomass.
e_lipid <- 39  # J / mg lipid
e_lean  <- 17  # J / mg lean dry mass

## Convert adult dry mass to mg, matching the rest of the energetic model.
adult_dry_mass_mg <- adult_dry_mass_ug * 0.001

## ---------------------------------------------------------
## 1. Create individual-level lipid-adjusted scenario dataset
## ---------------------------------------------------------

lipid_scenario_individual <- deb_individual %>%
  mutate(
    Ratio_chr = as.character(Ratio),
    Food_num = readr::parse_number(as.character(Food))
  ) %>%
  inner_join(
    lee_lipid_data,
    by = c("Ratio_chr" = "Ratio", "Food_num" = "Food_num")
  ) %>%
  mutate(
    ## Lipid and lean dry mass
    lipid_mass_mg = adult_dry_mass_mg * lipid_fraction_dry_mass,
    lean_mass_mg  = adult_dry_mass_mg * (1 - lipid_fraction_dry_mass),
    
    ## Lipid-adjusted dry-body energy density
    body_energy_density_lipid_adjusted =
      lipid_fraction_dry_mass * e_lipid +
      (1 - lipid_fraction_dry_mass) * e_lean,
    
    ## Baseline and lipid-adjusted adult body energy
    adult_body_energy_baseline_J =
      adult_dry_mass_mg * biomass_energy,
    
    adult_body_energy_lipid_adjusted_J =
      adult_dry_mass_mg * body_energy_density_lipid_adjusted,
    
    ## Baseline and lipid-adjusted maintenance
    maintenance_baseline_J_day =
      maintenance_J_day,
    
    maintenance_lipid_adjusted_J_day =
      maint_frac_adult_energy_per_day *
      adult_body_energy_lipid_adjusted_J,
    
    ## Recalculate safety margin under lipid-adjusted maintenance
    energetic_safety_margin_lipid_adjusted_J_day =
      assimilated_energy_J_day -
      reproductive_cost_J_day -
      maintenance_lipid_adjusted_J_day,
    
    energetic_safety_margin_fraction_lipid_adjusted =
      energetic_safety_margin_lipid_adjusted_J_day /
      assimilated_energy_J_day,
    
    energetically_feasible_lipid_adjusted =
      energetic_safety_margin_lipid_adjusted_J_day >= 0,
    
    ## Differences relative to baseline model
    delta_body_energy_density =
      body_energy_density_lipid_adjusted - biomass_energy,
    
    delta_maintenance_J_day =
      maintenance_lipid_adjusted_J_day - maintenance_baseline_J_day,
    
    delta_safety_margin_fraction =
      energetic_safety_margin_fraction_lipid_adjusted -
      energetic_safety_margin_fraction
  )

if (nrow(lipid_scenario_individual) == 0) {
  stop("No matching lipid-scenario rows found. Check Ratio and Food coding.")
}

## ---------------------------------------------------------
## 2. Treatment-level summary
## ---------------------------------------------------------

lipid_scenario_summary <- lipid_scenario_individual %>%
  group_by(Ratio_chr, Food_num) %>%
  summarise(
    n = n(),
    
    lipid_fraction_dry_mass = first(lipid_fraction_dry_mass),
    lipid_se = first(lipid_se),
    
    body_energy_density_baseline = biomass_energy,
    body_energy_density_lipid_adjusted =
      first(body_energy_density_lipid_adjusted),
    
    maintenance_baseline_J_day =
      median(maintenance_baseline_J_day, na.rm = TRUE),
    
    maintenance_lipid_adjusted_J_day =
      median(maintenance_lipid_adjusted_J_day, na.rm = TRUE),
    
    energetic_safety_margin_baseline =
      median(energetic_safety_margin_fraction, na.rm = TRUE),
    
    energetic_safety_margin_lipid_adjusted =
      median(energetic_safety_margin_fraction_lipid_adjusted, na.rm = TRUE),
    
    feasible_fraction_baseline =
      mean(energetically_feasible, na.rm = TRUE),
    
    feasible_fraction_lipid_adjusted =
      mean(energetically_feasible_lipid_adjusted, na.rm = TRUE),
    
    observed_mean_lifespan =
      mean(lifespan, na.rm = TRUE),
    
    observed_median_lifespan =
      median(lifespan, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    Ratio_chr = factor(
      Ratio_chr,
      levels = c("(1:16)", "(1:4)", "(1:2)")
    ),
    Ratio_label = dplyr::recode(
      as.character(Ratio_chr),
      "(1:16)" = "1:16",
      "(1:4)"  = "1:4",
      "(1:2)"  = "1:2"
    ),
    Ratio_label = factor(Ratio_label, levels = c("1:16", "1:4", "1:2"))
  )

lipid_scenario_summary

## ---------------------------------------------------------
## 3. Prepare plotting data
## ---------------------------------------------------------

lipid_body_energy_plot <- lipid_scenario_summary %>%
  select(
    Ratio_label,
    body_energy_density_baseline,
    body_energy_density_lipid_adjusted
  ) %>%
  pivot_longer(
    cols = c(
      body_energy_density_baseline,
      body_energy_density_lipid_adjusted
    ),
    names_to = "scenario",
    values_to = "body_energy_density"
  ) %>%
  mutate(
    scenario = dplyr::recode(
      scenario,
      body_energy_density_baseline = "Baseline",
      body_energy_density_lipid_adjusted = "Lipid-adjusted"
    ),
    scenario = factor(scenario, levels = c("Baseline", "Lipid-adjusted"))
  )

lipid_maintenance_plot <- lipid_scenario_summary %>%
  select(
    Ratio_label,
    maintenance_baseline_J_day,
    maintenance_lipid_adjusted_J_day
  ) %>%
  pivot_longer(
    cols = c(
      maintenance_baseline_J_day,
      maintenance_lipid_adjusted_J_day
    ),
    names_to = "scenario",
    values_to = "maintenance_J_day"
  ) %>%
  mutate(
    scenario = dplyr::recode(
      scenario,
      maintenance_baseline_J_day = "Baseline",
      maintenance_lipid_adjusted_J_day = "Lipid-adjusted"
    ),
    scenario = factor(scenario, levels = c("Baseline", "Lipid-adjusted"))
  )

lipid_margin_plot <- lipid_scenario_summary %>%
  select(
    Ratio_label,
    energetic_safety_margin_baseline,
    energetic_safety_margin_lipid_adjusted
  ) %>%
  pivot_longer(
    cols = c(
      energetic_safety_margin_baseline,
      energetic_safety_margin_lipid_adjusted
    ),
    names_to = "scenario",
    values_to = "energetic_safety_margin_fraction"
  ) %>%
  mutate(
    scenario = dplyr::recode(
      scenario,
      energetic_safety_margin_baseline = "Baseline",
      energetic_safety_margin_lipid_adjusted = "Lipid-adjusted"
    ),
    scenario = factor(scenario, levels = c("Baseline", "Lipid-adjusted"))
  )

lipid_feasible_plot <- lipid_scenario_summary %>%
  select(
    Ratio_label,
    feasible_fraction_baseline,
    feasible_fraction_lipid_adjusted
  ) %>%
  pivot_longer(
    cols = c(
      feasible_fraction_baseline,
      feasible_fraction_lipid_adjusted
    ),
    names_to = "scenario",
    values_to = "feasible_fraction"
  ) %>%
  mutate(
    scenario = dplyr::recode(
      scenario,
      feasible_fraction_baseline = "Baseline",
      feasible_fraction_lipid_adjusted = "Lipid-adjusted"
    ),
    scenario = factor(scenario, levels = c("Baseline", "Lipid-adjusted"))
  )

## ---------------------------------------------------------
## 4. Figure
## ---------------------------------------------------------

p_lipid_body_energy <- ggplot(
  lipid_body_energy_plot,
  aes(
    x = Ratio_label,
    y = body_energy_density,
    colour = scenario,
    group = scenario
  )
) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  scale_colour_viridis_d(option = "D", end = 0.85) +
  labs(
    x = "P:C rail",
    y = expression("Body energy density (J mg"^{-1}*")"),
    colour = "Scenario",
    title = "a"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p_lipid_maintenance <- ggplot(
  lipid_maintenance_plot,
  aes(
    x = Ratio_label,
    y = maintenance_J_day,
    colour = scenario,
    group = scenario
  )
) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  scale_colour_viridis_d(option = "D", end = 0.85) +
  labs(
    x = "P:C rail",
    y = expression("Maintenance (J day"^{-1}*")"),
    colour = "Scenario",
    title = "b"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p_lipid_margin <- ggplot(
  lipid_margin_plot,
  aes(
    x = Ratio_label,
    y = energetic_safety_margin_fraction,
    colour = scenario,
    group = scenario
  )
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  scale_colour_viridis_d(option = "D", end = 0.85) +
  labs(
    x = "P:C rail",
    y = expression("Median safety margin " * Delta[i] / A[i]),
    colour = "Scenario",
    title = "c"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p_lipid_feasible <- ggplot(
  lipid_feasible_plot,
  aes(
    x = Ratio_label,
    y = feasible_fraction,
    colour = scenario,
    group = scenario
  )
) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  scale_colour_viridis_d(option = "D", end = 0.85) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25)
  ) +
  labs(
    x = "P:C rail",
    y = "Energetically feasible fraction",
    colour = "Scenario",
    title = "d"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

figS_lipid_adjusted_scenario <- (
  p_lipid_body_energy + p_lipid_maintenance
) / (
  p_lipid_margin + p_lipid_feasible
) +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

figS_lipid_adjusted_scenario



## =========================================================
## 6. Figure objects
## =========================================================

theme_pub <- theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(size = 0, face = "bold", hjust = 0),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(0.45, "cm"),
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 15),
    strip.background = element_rect(fill = "grey90"),
    plot.margin = margin(8, 8, 8, 8)
  )


## ---------------------------------------------------------
## Main Figure 1. Conceptual framework (not used in the text but useful to replicate)
## ---------------------------------------------------------

fig1_framework <- ggplot() +
  
  ## -----------------------------
## Input / observed boxes
## -----------------------------
annotate("rect", xmin = 0.5, xmax = 2.7, ymin = 6.2, ymax = 7.2,
         fill = "grey95", colour = "black") +
  annotate("text", x = 1.6, y = 6.7,
           label = "Adult intake\nP[d], C[d]", size = 5) +
  
  annotate("rect", xmin = 0.5, xmax = 2.7, ymin = 4.6, ymax = 5.6,
           fill = "grey95", colour = "black") +
  annotate("text", x = 1.6, y = 5.1,
           label = "Observed egg\nproduction E[d]", size = 5) +
  
  annotate("rect", xmin = 0.5, xmax = 2.7, ymin = 3.0, ymax = 4.0,
           fill = "grey95", colour = "black") +
  annotate("text", x = 1.6, y = 3.5,
           label = "Adult dry mass\nW", size = 5) +
  
  annotate("rect", xmin = 0.5, xmax = 2.7, ymin = 1.4, ymax = 2.4,
           fill = "grey95", colour = "black") +
  annotate("text", x = 1.6, y = 1.9,
           label = "Nutritional target\nP*:C* = 1:4", size = 5) +
  
  ## -----------------------------
## Intermediate / derived boxes
## -----------------------------
annotate("rect", xmin = 4.0, xmax = 6.4, ymin = 6.2, ymax = 7.2,
         fill = "grey90", colour = "black") +
  annotate("text", x = 5.2, y = 6.7,
           label = "Assimilated energy\nA[d]", size = 5) +
  
  annotate("rect", xmin = 4.0, xmax = 6.4, ymin = 4.8, ymax = 5.8,
           fill = "grey90", colour = "black") +
  annotate("text", x = 5.2, y = 5.3,
           label = "Reproductive energetic\ncost R[d]", size = 5) +
  
  annotate("rect", xmin = 4.0, xmax = 6.4, ymin = 3.2, ymax = 4.2,
           fill = "grey90", colour = "black") +
  annotate("text", x = 5.2, y = 3.7,
           label = "Reproductive protein\ndemand P[R]", size = 5) +
  
  annotate("rect", xmin = 4.0, xmax = 6.4, ymin = 1.6, ymax = 2.6,
           fill = "grey90", colour = "black") +
  annotate("text", x = 5.2, y = 2.1,
           label = "Somatic maintenance\nM[d]", size = 5) +
  
  annotate("rect", xmin = 8.0, xmax = 10.8, ymin = 5.8, ymax = 7.0,
           fill = "grey90", colour = "black") +
  annotate("text", x = 9.4, y = 6.4,
           label = "Energetic safety margin\nDelta/A[d] = (A[d]-R[d]-M[d])/A[d]",
           size = 4.5) +
  
  annotate("rect", xmin = 8.0, xmax = 10.8, ymin = 4.1, ymax = 5.3,
           fill = "grey90", colour = "black") +
  annotate("text", x = 9.4, y = 4.7,
           label = "Reproductive energetic\nburden phi[R] = R[d]/A[d]",
           size = 4.5) +
  
  annotate("rect", xmin = 8.0, xmax = 10.8, ymin = 2.4, ymax = 3.6,
           fill = "grey90", colour = "black") +
  annotate("text", x = 9.4, y = 3.0,
           label = "Reproductive protein\nburden phi[P] = P[R]/(eta[P] P[d])",
           size = 4.3) +
  
  annotate("rect", xmin = 8.0, xmax = 10.8, ymin = 0.8, ymax = 2.0,
           fill = "grey90", colour = "black") +
  annotate("text", x = 9.4, y = 1.4,
           label = "Nutritional displacement\nD = distance from P:C = 1:4 target",
           size = 4.5) +
  
  ## -----------------------------
## Response box
## -----------------------------
annotate("rect", xmin = 12.2, xmax = 14.7, ymin = 3.5, ymax = 5.1,
         fill = "grey90", colour = "black") +
  annotate("text", x = 13.45, y = 4.3,
           label = "Age-specific\nmortality h(t)", size = 5) +
  
  ## -----------------------------
## Arrows: inputs to derived quantities
## -----------------------------
geom_segment(aes(x = 2.7, xend = 4.0, y = 6.7, yend = 6.7),
             arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  annotate("text", x = 3.35, y = 7.05,
           label = "eta[P], eta[C], e[P], e[C]", size = 3.6) +
  
  geom_segment(aes(x = 2.7, xend = 4.0, y = 5.1, yend = 5.3),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  annotate("text", x = 3.35, y = 5.55,
           label = "m[egg], e[B], kappa[R]", size = 3.6) +
  
  geom_segment(aes(x = 2.7, xend = 4.0, y = 5.1, yend = 3.7),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  annotate("text", x = 3.35, y = 4.15,
           label = "m[egg], q[egg], eta[P->R]", size = 3.5) +
  
  geom_segment(aes(x = 2.7, xend = 4.0, y = 3.5, yend = 2.1),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  annotate("text", x = 3.35, y = 2.55,
           label = "f[M], e[B]", size = 3.6) +
  
  ## intake + target -> distance
  geom_segment(aes(x = 2.7, xend = 8.0, y = 1.9, yend = 1.4),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  geom_segment(aes(x = 2.7, xend = 8.0, y = 6.3, yend = 1.6),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.5, linetype = 2) +
  
  ## -----------------------------
## Arrows: derived quantities to state variables
## -----------------------------
geom_segment(aes(x = 6.4, xend = 8.0, y = 6.7, yend = 6.4),
             arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  geom_segment(aes(x = 6.4, xend = 8.0, y = 5.3, yend = 6.3),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  geom_segment(aes(x = 6.4, xend = 8.0, y = 2.1, yend = 6.0),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  
  geom_segment(aes(x = 6.4, xend = 8.0, y = 5.3, yend = 4.7),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  
  geom_segment(aes(x = 6.4, xend = 8.0, y = 3.7, yend = 3.0),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  geom_segment(aes(x = 6.4, xend = 8.0, y = 6.7, yend = 3.2),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7, linetype = 2) +
  
  ## -----------------------------
## Arrows: state variables to mortality
## -----------------------------
geom_segment(aes(x = 10.8, xend = 12.2, y = 6.4, yend = 4.6),
             arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  geom_segment(aes(x = 10.8, xend = 12.2, y = 4.7, yend = 4.4),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  geom_segment(aes(x = 10.8, xend = 12.2, y = 3.0, yend = 4.2),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  geom_segment(aes(x = 10.8, xend = 12.2, y = 1.4, yend = 4.0),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  
  xlim(0, 15.2) +
  ylim(0.4, 7.6) +
  theme_void(base_size = 16) +
  labs(
    title = "Conceptual framework for DEB-inspired nutritional energetics and age-specific mortality"
  )

fig1_framework


### 

## Figure 2

treatment_points <- deb_individual %>%
  group_by(treatment, Ratio, Food) %>%
  summarise(
    P_daily_ug = mean(P_daily_ug, na.rm = TRUE),
    C_daily_ug = mean(C_daily_ug, na.rm = TRUE),
    distance_from_target_scaled = mean(distance_from_target_scaled, na.rm = TRUE),
    reproductive_energetic_burden = median(reproductive_energetic_burden, na.rm = TRUE),
    reproductive_protein_burden = median(reproductive_protein_burden, na.rm = TRUE),
    energetic_safety_margin_fraction = median(energetic_safety_margin_fraction, na.rm = TRUE),
    lifespan = mean(lifespan, na.rm = TRUE),
    .groups = "drop"
  )

x_lim <- range(deb_individual$P_daily_ug, na.rm = TRUE)
y_lim <- range(deb_individual$C_daily_ug, na.rm = TRUE)
x_pad <- diff(x_lim) * 0.04
y_pad <- diff(y_lim) * 0.04
x_lim <- c(max(0, x_lim[1] - x_pad), x_lim[2] + x_pad)
y_lim <- c(max(0, y_lim[1] - y_pad), y_lim[2] + y_pad)

grid_n <- 220

plot_grid <- expand.grid(
  P_daily_ug = seq(x_lim[1], x_lim[2], length.out = grid_n),
  C_daily_ug = seq(y_lim[1], y_lim[2], length.out = grid_n)
) %>%
  as_tibble()

fit_surface_gam <- function(data, response, k = 10) {
  response <- rlang::ensym(response)
  response_name <- rlang::as_string(response)
  
  data_fit <- data %>%
    filter(
      is.finite(.data[[response_name]]),
      is.finite(P_daily_ug),
      is.finite(C_daily_ug)
    )
  
  if (nrow(data_fit) < 10) {
    stop("Not enough finite treatment-level points to fit surface for: ", response_name)
  }
  
  if (sd(data_fit[[response_name]], na.rm = TRUE) == 0) {
    stop("Response has zero variation for surface: ", response_name)
  }
  
  gam(
    as.formula(
      paste0(response_name, " ~ s(P_daily_ug, C_daily_ug, k = ", k, ")")
    ),
    data = data_fit,
    method = "REML"
  )
}

gam_phiR <- fit_surface_gam(treatment_points, reproductive_energetic_burden, k = 10)
gam_phiPR <- fit_surface_gam(treatment_points, reproductive_protein_burden, k = 10)
gam_margin <- fit_surface_gam(treatment_points, energetic_safety_margin_fraction, k = 10)
gam_lifespan <- fit_surface_gam(treatment_points, lifespan, k = 10)

surface_phiR <- plot_grid %>% mutate(z = predict(gam_phiR, newdata = plot_grid))
surface_phiPR <- plot_grid %>% mutate(z = predict(gam_phiPR, newdata = plot_grid))
surface_margin <- plot_grid %>% mutate(z = predict(gam_margin, newdata = plot_grid))
surface_lifespan <- plot_grid %>% mutate(z = predict(gam_lifespan, newdata = plot_grid))

individual_layer <- geom_point(
  data = deb_individual,
  aes(x = P_daily_ug, y = C_daily_ug),
  inherit.aes = FALSE,
  colour = "black",
  alpha = 0.2,
  size = 1
)

treatment_layer <- geom_point(
  data = treatment_points,
  aes(x = P_daily_ug, y = C_daily_ug),
  inherit.aes = FALSE,
  shape = 21,
  fill = "white",
  colour = "black",
  stroke = 1.1,
  size = 2.5
)

target_layer <- geom_point(
  data = target_point,
  aes(x = target_P_ug, y = target_C_ug),
  inherit.aes = FALSE,
  shape = 18,
  colour = "red",
  size = 4,
  stroke = 1.4
) 

make_surface_panel <- function(
    surface_data,
    fill_label,
    panel_title,
    fill_limits = NULL,
    fill_oob = scales::squish,
    fill_trans = "identity"
) {
  p <- ggplot() +
    geom_raster(
      data = surface_data,
      aes(x = P_daily_ug, y = C_daily_ug, fill = z),
      interpolate = TRUE,
      alpha = 0.96
    ) +
    geom_contour(
      data = surface_data,
      aes(x = P_daily_ug, y = C_daily_ug, z = z),
      colour = "black",
      alpha = 0.85,
      linewidth = 0.5,
      bins = 8
    ) + 
    individual_layer +
    treatment_layer +
    target_layer +
    geom_abline(aes(intercept = 0, slope = 4), linewidth = 0.5, colour = "red") +
    coord_cartesian(xlim = x_lim, ylim = y_lim, expand = FALSE) +
    labs(
      x = "Daily protein intake (ug/fly/day)",
      y = "Daily carbohydrate intake (ug/fly/day)",
      fill = fill_label,
      title = panel_title
    ) +
    theme_pub
  
  if (is.null(fill_limits)) {
    p <- p +
      scale_fill_viridis_c(
        option = "D",
        trans = fill_trans,
        oob = fill_oob,
        na.value = "grey80"
      )
  } else {
    p <- p +
      scale_fill_viridis_c(
        option = "D",
        trans = fill_trans,
        limits = fill_limits,
        oob = fill_oob,
        na.value = "grey80"
      )
  }
  
  p
}


p_fig2C <- make_surface_panel(
  surface_data = surface_phiR,
  fill_label = expression(phi[R] == R[i] / A[i]),
  panel_title = "",
  fill_limits = range(surface_phiR$z, na.rm = TRUE)
)

p_fig2D <- make_surface_panel(
  surface_data = surface_phiPR,
  fill_label = expression(phi[P]),
  panel_title = "",
  fill_limits = range(surface_phiPR$z, na.rm = TRUE)
)

p_fig2B <- make_surface_panel(
  surface_data = surface_margin,
  fill_label = expression(Delta / A[d]),
  panel_title = "",
  fill_limits = range(surface_margin$z, na.rm = TRUE)
)

p_fig2A <- make_surface_panel(
  surface_data = surface_lifespan,
  fill_label = "Lifespan\n(days)",
  panel_title = "",
  fill_limits = range(surface_lifespan$z, na.rm = TRUE)
)

fig2_surface_overlay_polished <- (p_fig2A + p_fig2B ) /
  (p_fig2C + p_fig2D) +
  plot_annotation(theme = theme(
    plot.title = element_text(size = 24, face = "bold"),
    plot.subtitle = element_text(size = 17),
    plot.margin = margin(10, 10, 10, 10)
  )
  )

fig2_surface_overlay_polished

## Remove repeated axis titles from individual panels
p_fig2A <- p_fig2A + labs(x = NULL, y = NULL)
p_fig2B <- p_fig2B + labs(x = NULL, y = NULL)
p_fig2C <- p_fig2C + labs(x = NULL, y = NULL)
p_fig2D <- p_fig2D + labs(x = NULL, y = NULL)

## Combine panels
fig2_inner <- (p_fig2A + p_fig2B) /
  (p_fig2C + p_fig2D) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 24, face = "bold"),
      plot.subtitle = element_text(size = 17),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

## Add shared x and y labels
fig2_surface_overlay_polished <- cowplot::ggdraw(fig2_inner) +
  cowplot::draw_label(
    "Daily protein intake (\u00b5g/fly/day)",
    x = 0.5,
    y = 0.02,
    hjust = 0.5,
    vjust = 0.5,
    size = 18
  ) +
  cowplot::draw_label(
    "Daily carbohydrate intake (\u00b5g/fly/day)",
    x = 0.02,
    y = 0.5,
    angle = 90,
    hjust = 0.5,
    vjust = -0.2,
    size = 18
  )

fig2_surface_overlay_polished 
ggsave(
  file.path(output_dir, "Figure2.png"),
  fig2_surface_overlay_polished,
  width = 12,
  height = 10,
  dpi = 350
)
fig2_surface_overlay_polished




## 
## Figure 3

profile_data_rail_concentration <- deb_individual %>%
  group_by(treatment, Ratio, Food) %>%
  summarise(
    P_daily_ug = mean(P_daily_ug, na.rm = TRUE),
    C_daily_ug = mean(C_daily_ug, na.rm = TRUE),
    total_daily_ug = mean(total_daily_ug, na.rm = TRUE),
    energetic_safety_margin_fraction =
      median(energetic_safety_margin_fraction, na.rm = TRUE),
    distance_from_target_scaled =
      median(distance_from_target_scaled, na.rm = TRUE),
    reproductive_share_of_explicit_costs =
      median(reproductive_share_of_explicit_costs, na.rm = TRUE),
    reproductive_protein_burden =
      median(reproductive_protein_burden, na.rm = TRUE),
    lifespan = mean(lifespan, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    treatment = factor(treatment, levels = levels(surv_daily$treatment)),
    Ratio = factor(Ratio, levels = levels(surv_daily$Ratio)),
    Food = factor(Food, levels = levels(surv_daily$Food))
  )

pred_surv_rail_concentration <- profile_data_rail_concentration %>%
  tidyr::crossing(day = 1:max_day) %>%
  add_scaled_predictors(scale_ref) %>%
  mutate(
    death_prob = predict(
      final_dt_model,
      newdata = .,
      type = "response",
      re.form = NA
    )
  ) %>%
  group_by(treatment, Ratio, Food) %>%
  arrange(day, .by_group = TRUE) %>%
  mutate(
    survival_prob = cumprod(1 - death_prob)
  ) %>%
  ungroup()



ratio_order <- c("(0:1)", "(1:16)", "(1:8)", "(1:4)", "(1:2)", "(1:1)", "(2:1)")

## Use only levels actually present, in the desired order
ratio_order_use <- c(
  intersect(ratio_order, unique(as.character(pred_surv_rail_concentration$Ratio))),
  setdiff(unique(as.character(pred_surv_rail_concentration$Ratio)), ratio_order)
)

food_order_use <- pred_surv_rail_concentration %>%
  mutate(Food_chr = as.character(Food)) %>%
  distinct(Food_chr) %>%
  arrange(suppressWarnings(as.numeric(Food_chr))) %>%
  pull(Food_chr)

pred_surv_fig3 <- pred_surv_rail_concentration %>%
  mutate(
    Ratio_plot = factor(as.character(Ratio), levels = ratio_order_use),
    Food_plot  = factor(as.character(Food), levels = food_order_use)
  )


## ---------------------------------------------------------
## Panel A.
## ---------------------------------------------------------

p_death_prob <- ggplot(
  pred_surv_fig3,
  aes(
    x = day,
    y = death_prob,
    colour = Ratio_plot,
    group = Ratio_plot
  )
) +
  geom_line(linewidth = 1.05) +
  facet_wrap(~ Food_plot, ncol = length(food_order_use)) +
  labs(
    x = "Adult age (days)",
    y = "Predicted daily death probability",
    colour = "P:C rail",
    title = "a"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  scale_colour_viridis_d(option = "viridis")

p_death_prob


## ---------------------------------------------------------
## Panel B. Predicted survival probability
## ---------------------------------------------------------

p_survival_prob <- ggplot(
  pred_surv_fig3,
  aes(
    x = day,
    y = survival_prob,
    colour = Ratio_plot,
    group = Ratio_plot
  )
) +
  geom_line(linewidth = 1.05) +
  facet_wrap(~ Food_plot, ncol = length(food_order_use)) +
  labs(
    x = "Adult age (days)",
    y = "Predicted survival probability",
    colour = "P:C rail",
    title = "b"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  scale_colour_viridis_d()

p_survival_prob


## ---------------------------------------------------------
## Combined Figure 3
## ---------------------------------------------------------

fig3_survival_profiles <- p_death_prob / p_survival_prob +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig3_survival_profiles



ggsave(
  file.path(output_dir, "Figure3.png"),
  fig3_survival_profiles,
  width = 16,
  height = 9,
  dpi = 350
)






### Figure 4

ratio_order <- c("(0:1)", "(1:16)", "(1:8)", "(1:4)", "(1:2)", "(1:1)", "(2:1)")
food_order <- c("45", "90", "180", "360")


calculate_ccc <- function(obs, pred) {
  
  ccc <- DescTools::CCC(
    x = obs,
    y = pred,
    ci = "z-transform"
  )
  
  tibble(
    CCC = ccc$rho.c[1],
    lower_CI = ccc$rho.c[2],
    upper_CI = ccc$rho.c[3]
  )
}




observed_km_rail_food <- deb_individual %>%
  group_by(Ratio, Food) %>%
  group_modify(~{
    km <- survival::survfit(
      survival::Surv(lifespan, death_event) ~ 1,
      data = .x
    )
    
    tibble(
      day = km$time,
      observed_survival = km$surv,
      n_risk = km$n.risk,
      n_event = km$n.event
    )
  }) %>%
  ungroup()

pred_ind_daily_rail_food <- deb_individual %>%
  select(
    fly_id,
    treatment,
    Ratio,
    Food,
    lifespan,
    lifespan_day,
    P_daily_ug,
    C_daily_ug,
    total_daily_ug,
    energetic_safety_margin_fraction,
    distance_from_target_scaled,
    reproductive_share_of_explicit_costs,
    reproductive_protein_burden
  ) %>%
  tidyr::crossing(day = 1:max_day) %>%
  mutate(
    treatment = factor(treatment, levels = levels(surv_daily$treatment)),
    Ratio = factor(Ratio, levels = levels(surv_daily$Ratio)),
    Food = factor(Food, levels = levels(surv_daily$Food))
  ) %>%
  add_scaled_predictors(scale_ref) %>%
  mutate(
    pred_death_prob = predict(
      final_dt_model,
      newdata = .,
      type = "response",
      re.form = NA
    )
  ) %>%
  group_by(fly_id) %>%
  arrange(day, .by_group = TRUE) %>%
  mutate(
    predicted_survival = cumprod(1 - pred_death_prob)
  ) %>%
  ungroup()

predicted_survival_rail_food <- pred_ind_daily_rail_food %>%
  group_by(Ratio, Food, day) %>%
  summarise(
    predicted_survival = mean(predicted_survival, na.rm = TRUE),
    predicted_survival_se =
      sd(predicted_survival, na.rm = TRUE) / sqrt(n()),
    predicted_lower =
      pmax(0, predicted_survival - 1.96 * predicted_survival_se),
    predicted_upper =
      pmin(1, predicted_survival + 1.96 * predicted_survival_se),
    .groups = "drop"
  )



predicted_lifespan_by_individual <- pred_ind_daily_rail_food %>%
  group_by(fly_id) %>%
  summarise(
    predicted_lifespan = sum(predicted_survival, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    deb_individual %>%
      select(
        fly_id,
        lifespan,
        treatment,
        Ratio,
        Food,
        P_daily_ug,
        C_daily_ug
      ),
    by = "fly_id"
  ) %>%
  mutate(
    Ratio = factor(as.character(Ratio), levels = ratio_order),
    Food  = factor(as.character(Food), levels = food_order)
  )


food_order <- predicted_lifespan_by_individual %>%
  mutate(Food_chr = as.character(Food)) %>%
  distinct(Food_chr) %>%
  arrange(suppressWarnings(as.numeric(Food_chr))) %>%
  pull(Food_chr)

## ---------------------------------------------------------
## 4. Treatment-level observed vs predicted mean lifespan
## ---------------------------------------------------------

treatment_observed_predicted_lifespan <- predicted_lifespan_by_individual %>%
  group_by(treatment, Ratio, Food) %>%
  summarise(
    n = n(),
    observed_lifespan = mean(lifespan, na.rm = TRUE),
    predicted_lifespan = mean(predicted_lifespan, na.rm = TRUE),
    residual_lifespan = observed_lifespan - predicted_lifespan,
    P_daily_ug = mean(P_daily_ug, na.rm = TRUE),
    C_daily_ug = mean(C_daily_ug, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Ratio = factor(as.character(Ratio), levels = ratio_order),
    Food  = factor(as.character(Food), levels = food_order)
  )

ccc_mean <- calculate_ccc(
  treatment_observed_predicted_lifespan$observed_lifespan,
  treatment_observed_predicted_lifespan$predicted_lifespan
)

## ---------------------------------------------------------
## 5. Treatment-level observed vs predicted median lifespan
## ---------------------------------------------------------

treatment_observed_predicted_median <- predicted_lifespan_by_individual %>%
  group_by(treatment, Ratio, Food) %>%
  summarise(
    observed_median_lifespan =
      median(lifespan, na.rm = TRUE),
    
    predicted_median_lifespan =
      median(predicted_lifespan, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    Ratio = factor(as.character(Ratio), levels = ratio_order),
    Food  = factor(as.character(Food), levels = food_order)
  )

ccc_median <- calculate_ccc(
  treatment_observed_predicted_median$observed_median_lifespan,
  treatment_observed_predicted_median$predicted_median_lifespan
)

#### Figure

fig4a_mean_fit <- ggplot(
  treatment_observed_predicted_lifespan,
  aes(
    x = observed_lifespan,
    y = predicted_lifespan,
    colour = Ratio,
    shape = Food
  )
) +
  geom_point(size = 3.5) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.1,
    vjust = -0.5,
    size = 4.5,
    label = paste0(
      "CCC = ",
      round(ccc_mean$CCC, 3),
      "\n95% CI: ",
      round(ccc_mean$lower_CI, 3),
      "–",
      round(ccc_mean$upper_CI, 3)
    )
  ) +
  labs(
    x = "Observed mean lifespan (days)",
    y = "Predicted mean lifespan (days)",
    colour = "P:C rail",
    shape = "Food concentration"
  ) +
  scale_color_viridis_d(
    option = "D",
    drop = FALSE
  ) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "bottom"
  )




fig4b_median_fit <- ggplot(
  treatment_observed_predicted_median,
  aes(
    x = observed_median_lifespan,
    y = predicted_median_lifespan,
    colour = Ratio,
    shape = Food
  )
) +
  geom_point(size = 3.5) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.1,
    vjust = -0.5,
    size = 4.5,
    label = paste0(
      "CCC = ",
      round(ccc_median$CCC, 3),
      "\n95% CI: ",
      round(ccc_median$lower_CI, 3),
      "–",
      round(ccc_median$upper_CI, 3)
    )
  ) +
  labs(
    x = "Observed median lifespan (days)",
    y = "Predicted median lifespan (days)",
    colour = "P:C rail",
    shape = "Food concentration"
  ) +
  scale_color_viridis_d(
    option = "D",
    drop = FALSE
  ) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "bottom"
  )


fig4_lifespan_fit <- (fig4a_mean_fit + fig4b_median_fit) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom"
  )

fig4_lifespan_fit

ggsave(
  file.path(output_dir, "Figure4_observed_vs_predicted_lifespan.png"),
  fig4_lifespan_fit,
  width = 11,
  height = 5.5,
  dpi = 350
)



## ---------------------------------------------------------
## Supplementary Figure S1. Observed vs predicted survival by rail x concentration
## ---------------------------------------------------------
observed_km_rail_food <- deb_individual %>%
  group_by(Ratio, Food) %>%
  group_modify(~{
    km <- survival::survfit(
      survival::Surv(lifespan, death_event) ~ 1,
      data = .x
    )
    
    tibble(
      day = km$time,
      observed_survival = km$surv,
      n_risk = km$n.risk,
      n_event = km$n.event
    )
  }) %>%
  ungroup()

pred_ind_daily_rail_food <- deb_individual %>%
  select(
    fly_id,
    treatment,
    Ratio,
    Food,
    lifespan,
    lifespan_day,
    P_daily_ug,
    C_daily_ug,
    total_daily_ug,
    energetic_safety_margin_fraction,
    distance_from_target_scaled,
    reproductive_share_of_explicit_costs,
    reproductive_protein_burden
  ) %>%
  tidyr::crossing(day = 1:max_day) %>%
  mutate(
    treatment = factor(treatment, levels = levels(surv_daily$treatment)),
    Ratio = factor(Ratio, levels = levels(surv_daily$Ratio)),
    Food = factor(Food, levels = levels(surv_daily$Food))
  ) %>%
  add_scaled_predictors(scale_ref) %>%
  mutate(
    pred_death_prob = predict(
      final_dt_model,
      newdata = .,
      type = "response",
      re.form = NA
    )
  ) %>%
  group_by(fly_id) %>%
  arrange(day, .by_group = TRUE) %>%
  mutate(
    predicted_survival = cumprod(1 - pred_death_prob)
  ) %>%
  ungroup()

predicted_survival_rail_food <- pred_ind_daily_rail_food %>%
  group_by(Ratio, Food, day) %>%
  summarise(
    predicted_survival = mean(predicted_survival, na.rm = TRUE),
    predicted_survival_se =
      sd(predicted_survival, na.rm = TRUE) / sqrt(n()),
    predicted_lower =
      pmax(0, predicted_survival - 1.96 * predicted_survival_se),
    predicted_upper =
      pmin(1, predicted_survival + 1.96 * predicted_survival_se),
    .groups = "drop"
  )


ratio_order <- c("(0:1)", "(1:16)", "(1:8)", "(1:4)", "(1:2)", "(1:1)", "(2:1)")
predicted_survival_rail_food$Ratio <- factor(as.character(predicted_survival_rail_food$Ratio), levels = ratio_order)

figS1_observed_predicted_survival_rail_food <- ggplot() +
  geom_ribbon(
    data = predicted_survival_rail_food,
    aes(
      x = day,
      ymin = predicted_lower,
      ymax = predicted_upper
    ),
    fill = "steelblue3",
    alpha = 0.12
  ) +
  geom_line(
    data = predicted_survival_rail_food,
    aes(
      x = day,
      y = predicted_survival
    ),
    colour = "steelblue3",
    linewidth = 0.9
  ) +
  geom_step(
    data = observed_km_rail_food,
    aes(
      x = day,
      y = observed_survival
    ),
    colour = "grey2",
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  facet_grid(Ratio ~ Food) +
  labs(
    x = "Adult age (days)",
    y = "Survival probability",
    title = "",
    subtitle = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )


figS1_observed_predicted_survival_rail_food

ggsave(
  file.path(output_dir, "FigureS1.png"),
  figS1_observed_predicted_survival_rail_food,
  width = 12,
  height = 10,
  dpi = 350
)


## =========================================================
## 7. Sensitivity analyses
## =========================================================


## ---------------------------------------------------------
## 7A. One-at-a-time wide biological sensitivity
## ---------------------------------------------------------

baseline_parameters <- tibble(
  scenario_name = "baseline",
  egg_dry_mass_ug = egg_dry_mass_ug,
  adult_dry_mass_ug = adult_dry_mass_ug,
  eta_P = eta_P,
  eta_C = eta_C,
  kap_R = kap_R,
  maint_frac_adult_energy_per_day = maint_frac_adult_energy_per_day,
  egg_protein_fraction_dry = egg_protein_fraction_dry,
  eta_protein_to_egg = eta_protein_to_egg
)

oat_sensitivity_grid <- bind_rows(
  tibble(parameter = "egg_dry_mass_ug", value = c(0.25, 0.50, 1.00, 2.25, 3.00, 5.00, 8.00)),
  tibble(parameter = "adult_dry_mass_ug", value = c(100, 150, 200, 300, 400, 500, 750)),
  tibble(parameter = "eta_P", value = c(0.10, 0.25, 0.50, 0.75, 0.95, 1.00)),
  tibble(parameter = "eta_C", value = c(0.10, 0.25, 0.50, 0.75, 0.95, 1.00)),
  tibble(parameter = "kap_R", value = c(0.40, 0.50, 0.70, 0.85, 0.95, 1.00)),
  tibble(parameter = "maint_frac_adult_energy_per_day", value = c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50)),
  tibble(parameter = "egg_protein_fraction_dry", value = c(0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80)),
  tibble(parameter = "eta_protein_to_egg", value = c(0.30, 0.50, 0.70, 0.85, 0.95, 1.00))
) %>%
  mutate(
    scenario_id = row_number(),
    egg_dry_mass_ug_sens =
      ifelse(parameter == "egg_dry_mass_ug", value, baseline_parameters$egg_dry_mass_ug),
    adult_dry_mass_ug_sens =
      ifelse(parameter == "adult_dry_mass_ug", value, baseline_parameters$adult_dry_mass_ug),
    eta_P_sens =
      ifelse(parameter == "eta_P", value, baseline_parameters$eta_P),
    eta_C_sens =
      ifelse(parameter == "eta_C", value, baseline_parameters$eta_C),
    kap_R_sens =
      ifelse(parameter == "kap_R", value, baseline_parameters$kap_R),
    maint_frac_sens =
      ifelse(parameter == "maint_frac_adult_energy_per_day", value, baseline_parameters$maint_frac_adult_energy_per_day),
    egg_protein_fraction_sens =
      ifelse(parameter == "egg_protein_fraction_dry", value, baseline_parameters$egg_protein_fraction_dry),
    eta_protein_to_egg_sens =
      ifelse(parameter == "eta_protein_to_egg", value, baseline_parameters$eta_protein_to_egg),
    add_my_pet_reference =
      case_when(
        parameter == "kap_R" & value == amp_kap_R ~ TRUE,
        TRUE ~ FALSE
      )
  )


summarise_sensitivity_deb <- function(tmp) {
  tmp %>%
    summarise(
      n = n(),
      phi_R_median = median(reproductive_energetic_burden, na.rm = TRUE),
      phi_R_low = quantile(reproductive_energetic_burden, 0.025, na.rm = TRUE),
      phi_R_high = quantile(reproductive_energetic_burden, 0.975, na.rm = TRUE),
      phi_PR_median = median(reproductive_protein_burden, na.rm = TRUE),
      phi_PR_low = quantile(reproductive_protein_burden, 0.025, na.rm = TRUE),
      phi_PR_high = quantile(reproductive_protein_burden, 0.975, na.rm = TRUE),
      rho_R_median = median(reproductive_share_of_explicit_costs, na.rm = TRUE),
      rho_R_low = quantile(reproductive_share_of_explicit_costs, 0.025, na.rm = TRUE),
      rho_R_high = quantile(reproductive_share_of_explicit_costs, 0.975, na.rm = TRUE),
      safety_margin_fraction_median = median(energetic_safety_margin_fraction, na.rm = TRUE),
      safety_margin_fraction_low = quantile(energetic_safety_margin_fraction, 0.025, na.rm = TRUE),
      safety_margin_fraction_high = quantile(energetic_safety_margin_fraction, 0.975, na.rm = TRUE),
      feasible_fraction = mean(energetically_feasible, na.rm = TRUE),
      protein_feasible_fraction = mean(protein_feasible_for_reproduction, na.rm = TRUE),
      explicit_costs_exceed_assimilation_fraction =
        mean(explicit_cost_fraction_of_assimilation > 1, na.rm = TRUE)
    )
}


run_oat_scenario <- function(row_i) {
  
  scenario <- oat_sensitivity_grid[row_i, ]
  
  tmp <- calculate_individual_deb(
    data = dat,
    target_P = target_P_ug,
    target_C = target_C_ug,
    target_total = target_total_intake,
    distance_method = distance_method,
    eP = eP,
    eC = eC,
    biomass_energy = biomass_energy,
    egg_dry_mass_ug = scenario$egg_dry_mass_ug_sens,
    adult_dry_mass_ug = scenario$adult_dry_mass_ug_sens,
    eta_P = scenario$eta_P_sens,
    eta_C = scenario$eta_C_sens,
    kap_R = scenario$kap_R_sens,
    maint_frac_adult_energy_per_day = scenario$maint_frac_sens,
    egg_protein_fraction_dry = scenario$egg_protein_fraction_sens,
    eta_protein_to_egg = scenario$eta_protein_to_egg_sens
  )
  
  summarise_sensitivity_deb(tmp) %>%
    mutate(
      scenario_id = scenario$scenario_id,
      parameter = scenario$parameter,
      value = scenario$value,
      egg_dry_mass_ug = scenario$egg_dry_mass_ug_sens,
      adult_dry_mass_ug = scenario$adult_dry_mass_ug_sens,
      eta_P = scenario$eta_P_sens,
      eta_C = scenario$eta_C_sens,
      kap_R = scenario$kap_R_sens,
      maint_frac_adult_energy_per_day = scenario$maint_frac_sens,
      egg_protein_fraction_dry = scenario$egg_protein_fraction_sens,
      eta_protein_to_egg = scenario$eta_protein_to_egg_sens,
      add_my_pet_reference = scenario$add_my_pet_reference
    )
}

oat_sensitivity_summary <- map_dfr(
  seq_len(nrow(oat_sensitivity_grid)),
  run_oat_scenario
)



convert_pM_to_maintenance_fraction <- function(
    p_M_J_d_cm3,
    adult_dry_mass_ug,
    biomass_energy,
    dry_mass_density_mg_cm3 = 200
) {
  
  adult_dry_mass_mg <- adult_dry_mass_ug / 1000
  
  structural_volume_cm3 <-
    adult_dry_mass_mg / dry_mass_density_mg_cm3
  
  maintenance_J_day <-
    p_M_J_d_cm3 * structural_volume_cm3
  
  adult_body_energy_J <-
    adult_dry_mass_ug * 0.001 * biomass_energy
  
  maintenance_J_day / adult_body_energy_J
}


amp_sensitivity_grid <- tidyr::crossing(
  p_M_multiplier = c(0.05, 0.10, 0.25, 0.50, 1, 2, 4, 8),
  dry_mass_density_mg_cm3 = c(100, 200, 400),
  kap_R_sens = c(0.70, 0.85, amp_kap_R, 1.00),
  adult_dry_mass_ug_sens = c(150, 300, 500),
  egg_dry_mass_ug_sens = c(1.00, egg_dry_mass_ug, 5.00)
) %>%
  mutate(
    scenario_id = row_number(),
    p_M_sens = amp_p_M_J_d_cm3 * p_M_multiplier,
    maint_frac_from_pM =
      convert_pM_to_maintenance_fraction(
        p_M_J_d_cm3 = p_M_sens,
        adult_dry_mass_ug = adult_dry_mass_ug_sens,
        biomass_energy = biomass_energy,
        dry_mass_density_mg_cm3 = dry_mass_density_mg_cm3
      )
  )


run_amp_scenario <- function(row_i) {
  
  scenario <- amp_sensitivity_grid[row_i, ]
  
  tmp <- calculate_individual_deb(
    data = dat,
    target_P = target_P_ug,
    target_C = target_C_ug,
    target_total = target_total_intake,
    distance_method = distance_method,
    eP = eP,
    eC = eC,
    biomass_energy = biomass_energy,
    egg_dry_mass_ug = scenario$egg_dry_mass_ug_sens,
    adult_dry_mass_ug = scenario$adult_dry_mass_ug_sens,
    eta_P = eta_P,
    eta_C = eta_C,
    kap_R = scenario$kap_R_sens,
    maint_frac_adult_energy_per_day = scenario$maint_frac_from_pM,
    egg_protein_fraction_dry = egg_protein_fraction_dry,
    eta_protein_to_egg = eta_protein_to_egg
  )
  
  summarise_sensitivity_deb(tmp) %>%
    mutate(
      scenario_id = scenario$scenario_id,
      p_M_multiplier = scenario$p_M_multiplier,
      p_M_sens = scenario$p_M_sens,
      dry_mass_density_mg_cm3 = scenario$dry_mass_density_mg_cm3,
      kap_R = scenario$kap_R_sens,
      adult_dry_mass_ug = scenario$adult_dry_mass_ug_sens,
      egg_dry_mass_ug = scenario$egg_dry_mass_ug_sens,
      maint_frac_from_pM = scenario$maint_frac_from_pM
    )
}

amp_sensitivity_summary <- map_dfr(
  seq_len(nrow(amp_sensitivity_grid)),
  run_amp_scenario
)




set.seed(123)

n_global_scenarios <- 1000

global_sensitivity_grid <- tibble(
  scenario_id = seq_len(n_global_scenarios),
  egg_dry_mass_ug = runif(n_global_scenarios, 0.5, 6.0),
  adult_dry_mass_ug = runif(n_global_scenarios, 100, 600),
  eta_P = runif(n_global_scenarios, 0.25, 1.00),
  eta_C = runif(n_global_scenarios, 0.25, 1.00),
  kap_R = runif(n_global_scenarios, 0.50, 1.00),
  egg_protein_fraction_dry = runif(n_global_scenarios, 0.25, 0.75),
  eta_protein_to_egg = runif(n_global_scenarios, 0.50, 1.00),
  use_amp_pM = sample(c(TRUE, FALSE), n_global_scenarios, replace = TRUE),
  p_M_multiplier = exp(runif(n_global_scenarios, log(0.10), log(4.00))),
  dry_mass_density_mg_cm3 = runif(n_global_scenarios, 100, 400)
) %>%
  mutate(
    p_M_sens = amp_p_M_J_d_cm3 * p_M_multiplier,
    maint_frac_from_pM =
      convert_pM_to_maintenance_fraction(
        p_M_J_d_cm3 = p_M_sens,
        adult_dry_mass_ug = adult_dry_mass_ug,
        biomass_energy = biomass_energy,
        dry_mass_density_mg_cm3 = dry_mass_density_mg_cm3
      ),
    maint_frac_direct = exp(runif(n_global_scenarios, log(0.005), log(0.50))),
    maint_frac_adult_energy_per_day =
      ifelse(use_amp_pM, maint_frac_from_pM, maint_frac_direct)
  )


run_global_scenario <- function(row_i) {
  
  scenario <- global_sensitivity_grid[row_i, ]
  
  tmp <- calculate_individual_deb(
    data = dat,
    target_P = target_P_ug,
    target_C = target_C_ug,
    target_total = target_total_intake,
    distance_method = distance_method,
    eP = eP,
    eC = eC,
    biomass_energy = biomass_energy,
    egg_dry_mass_ug = scenario$egg_dry_mass_ug,
    adult_dry_mass_ug = scenario$adult_dry_mass_ug,
    eta_P = scenario$eta_P,
    eta_C = scenario$eta_C,
    kap_R = scenario$kap_R,
    maint_frac_adult_energy_per_day = scenario$maint_frac_adult_energy_per_day,
    egg_protein_fraction_dry = scenario$egg_protein_fraction_dry,
    eta_protein_to_egg = scenario$eta_protein_to_egg
  )
  
  summarise_sensitivity_deb(tmp) %>%
    mutate(
      scenario_id = scenario$scenario_id,
      egg_dry_mass_ug = scenario$egg_dry_mass_ug,
      adult_dry_mass_ug = scenario$adult_dry_mass_ug,
      eta_P = scenario$eta_P,
      eta_C = scenario$eta_C,
      kap_R = scenario$kap_R,
      egg_protein_fraction_dry = scenario$egg_protein_fraction_dry,
      eta_protein_to_egg = scenario$eta_protein_to_egg,
      use_amp_pM = scenario$use_amp_pM,
      p_M_multiplier = scenario$p_M_multiplier,
      p_M_sens = scenario$p_M_sens,
      dry_mass_density_mg_cm3 = scenario$dry_mass_density_mg_cm3,
      maint_frac_adult_energy_per_day =
        scenario$maint_frac_adult_energy_per_day
    )
}

global_sensitivity_summary <- map_dfr(
  seq_len(nrow(global_sensitivity_grid)),
  run_global_scenario
)











## ---------------------------------------------------------
## Fig S
## ---------------------------------------------------------

parameter_labels <- c(
  egg_dry_mass_ug = "Egg dry mass (ug)",
  adult_dry_mass_ug = "Adult dry mass (ug)",
  eta_P = "Protein assimilation efficiency",
  eta_C = "Carbohydrate assimilation efficiency",
  kap_R = "Reproductive efficiency (kap_R)",
  maint_frac_adult_energy_per_day = "Maintenance fraction per day",
  egg_protein_fraction_dry = "Egg protein fraction",
  eta_protein_to_egg = "Protein-to-egg efficiency"
)

figS5_oat_phiR <- ggplot(oat_sensitivity_summary, aes(x = value, y = phi_R_median)) +
  geom_ribbon(aes(ymin = phi_R_low, ymax = phi_R_high), alpha = 0.25) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ parameter, scales = "free_x", labeller = as_labeller(parameter_labels)) +
  labs(
    x = "Parameter value",
    y = expression(phi[R] == R[d] / A[d]),
    title = "a"
  ) +
  theme_bw(base_size = 13)

figS5_oat_phiPR <- ggplot(oat_sensitivity_summary, aes(x = value, y = phi_PR_median)) +
  geom_ribbon(aes(ymin = phi_PR_low, ymax = phi_PR_high), alpha = 0.25) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ parameter, scales = "free_x", labeller = as_labeller(parameter_labels)) +
  labs(
    x = "Parameter value",
    y = expression(phi[P*","*R]),
    title = "b"
  ) +
  theme_bw(base_size = 13)

figS5_oat_margin <- ggplot(oat_sensitivity_summary, aes(x = value, y = safety_margin_fraction_median)) +
  geom_ribbon(aes(ymin = safety_margin_fraction_low, ymax = safety_margin_fraction_high), alpha = 0.25) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ parameter, scales = "free_x", labeller = as_labeller(parameter_labels)) +
  labs(
    x = "Parameter value",
    y = expression(Delta / A[d]),
    title = "c"
  ) +
  theme_bw(base_size = 13)

figS5_oat_feasible <- ggplot(oat_sensitivity_summary, aes(x = value, y = feasible_fraction)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ parameter, scales = "free_x", labeller = as_labeller(parameter_labels)) +
  ylim(0, 1) +
  labs(
    x = "Parameter value",
    y = "Fraction energetically feasible",
    title = "d"
  ) +
  theme_bw(base_size = 13)

figS5_oat_sensitivity <- (figS5_oat_phiR + figS5_oat_phiPR) /
  (figS5_oat_margin + figS5_oat_feasible)



ggsave(
  file.path(output_dir, "FigureS2.png"),
  figS5_oat_sensitivity,
  width = 17,
  height = 15,
  dpi = 350
)






figS7_global_phiR <- ggplot(
  global_sensitivity_summary,
  aes(x = phi_R_median)
) +
  geom_histogram(bins = 50) +
  labs(
    x = expression("Median " * phi[R]),
    y = "Number of scenarios",
    title = ""
  ) +
  theme_bw(base_size = 14)

figS7_global_margin <- ggplot(
  global_sensitivity_summary,
  aes(x = safety_margin_fraction_median)
) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(
    x = expression("Median " * Delta / A[d]),
    y = "Number of scenarios",
    title = ""
  ) +
  theme_bw(base_size = 14)

figS7_global_feasible <- ggplot(
  global_sensitivity_summary,
  aes(x = feasible_fraction)
) +
  geom_histogram(bins = 50) +
  labs(
    x = "Energetically feasible fraction",
    y = "Number of scenarios",
    title = ""
  ) +
  theme_bw(base_size = 14)

figS7_global_amp <- ggplot(
  global_sensitivity_summary,
  aes(
    x = maint_frac_adult_energy_per_day,
    y = safety_margin_fraction_median
  )
) +
  geom_point(alpha = 0.45, size = 1.4) +
  scale_x_log10() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(
    x = "Maintenance fraction per day",
    y = expression("Median " * Delta / A[d]),
    colour = "Maintenance\nfrom p_M?",
    title = ""
  ) +
  theme_bw(base_size = 14)

figS7_global_sensitivity <- (figS7_global_phiR + figS7_global_margin) /
  (figS7_global_feasible + figS7_global_amp)


ggsave(
  file.path(output_dir, "FigureS3.png"),
  figS7_global_sensitivity,
  width = 15,
  height = 9,
  dpi = 350
)





## ---------------------------------------------------------
## Save lipid-adjusted scenario figure and table
## ---------------------------------------------------------

ggsave(
  filename = file.path(output_dir, "FigureS4.png"),
  plot = figS_lipid_adjusted_scenario,
  width = 12,
  height = 7,
  dpi = 600
)
