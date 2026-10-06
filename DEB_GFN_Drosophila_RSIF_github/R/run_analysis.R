# Reproducible analysis for the DEB-GFN Drosophila manuscript
# Run from the repository root with: Rscript R/run_analysis.R

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
  "ggforce",
  "car",
  "DescTools",
  "cowplot"
)

missing_packages <- packages_needed[
  !vapply(packages_needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    "Install required packages before running the analysis: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

invisible(lapply(packages_needed, library, character.only = TRUE))

input_file <- file.path("data", "Dataset_Lee_Nutrigonometry.csv")
output_dir <- "outputs"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

eP <- 17                 # J / mg protein
eC <- 17                 # J / mg carbohydrate
biomass_energy <- 22     # J / mg dry biomass

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

target_ratio_P <- 1
target_ratio_C <- 4

target_fraction_P <- target_ratio_P / (target_ratio_P + target_ratio_C)
target_fraction_C <- target_ratio_C / (target_ratio_P + target_ratio_C)

distance_method <- "euclidean"

dailyeggs_is_lifetime_over_lifespan <- TRUE

parameter_provenance <- tibble::tribble(
  ~parameter, ~symbol, ~baseline_value, ~used_in, ~source_type, ~literature_basis, ~revision_role,

  "Protein energy density", "eP", "17 J mg^-1",
  "Baseline energetic accounting",
  "Biochemical / energetic convention",
  "Retained because protein and carbohydrate have similar gross energy density; not estimated from Lee et al.",
  "Defines gross energy contribution of protein intake in the baseline model.",

  "Carbohydrate energy density", "eC", "17 J mg^-1",
  "Baseline energetic accounting",
  "Biochemical / energetic convention",
  "Retained because simple carbohydrates and protein have similar gross energy density; not estimated from Lee et al.",
  "Defines gross energy contribution of carbohydrate intake in the baseline model.",

  "Protein assimilation efficiency", "eta_P", "0.75 baseline; varied in revision scenarios",
  "Baseline model and unequal-assimilation scenarios",
  "Literature-informed sensitivity/scenario parameter",
  "The literature audit supports nutrient-specific physiology but did not identify directly transferable adult female, diet-specific eta_P values matching the Lee et al. design.",
  "Reviewer-requested parameter; varied independently from eta_C to test dependence on the equal-assimilation assumption.",

  "Carbohydrate assimilation efficiency", "eta_C", "0.75 baseline; varied in revision scenarios",
  "Baseline model and unequal-assimilation scenarios",
  "Literature-informed sensitivity/scenario parameter",
  "The literature audit supports nutrient-specific physiology but did not identify directly transferable adult female, diet-specific eta_C values matching the Lee et al. design.",
  "Reviewer-requested parameter; varied independently from eta_P to test dependence on the equal-assimilation assumption.",

  "Adult dry mass", "W", "300 ug baseline; varied in sensitivity analysis",
  "Maintenance calculation and sensitivity analysis",
  "Literature-informed sensitivity parameter",
  "Adult dry-mass and metabolic-physiology studies were used to justify broad sensitivity analyses rather than to replace the Lee et al. baseline.",
  "Evaluates how uncertainty in adult body mass affects maintenance and energetic feasibility.",

  "Daily maintenance fraction", "f_M", "0.10 day^-1 baseline; varied in sensitivity and penalty scenarios",
  "Maintenance calculation and reviewer-response penalty scenarios",
  "Exploratory / literature-informed sensitivity parameter",
  "No directly transferable adult female, diet-specific DEB-style maintenance flux was identified for the Lee et al. design.",
  "Retained as baseline, then extended with diet-imbalance and protein-excess penalty scenarios.",

  "Diet-imbalance maintenance penalty", "lambda_D", "0 in baseline; varied in revised scenario analysis",
  "Diet-imbalance maintenance-penalty scenario",
  "Exploratory scenario parameter",
  "Motivated by literature showing diet imbalance affects feeding, body composition, metabolism and lifespan, but not estimated from a single transferable study.",
  "Tests whether conclusions depend on diet-independent maintenance.",

  "Protein-excess maintenance penalty", "lambda_P", "0 in baseline; varied in revised scenario analysis",
  "Protein-excess maintenance-penalty scenario",
  "Exploratory scenario parameter",
  "Motivated by Drosophila high-protein diet literature and geometric-stoichiometry arguments about costs of protein excess.",
  "Tests whether high-protein rails impose additional explicit energetic costs.",

  "Reproductive efficiency", "kappa_R", "0.85 baseline; varied in sensitivity analysis",
  "Reproductive energetic cost",
  "DEB-inspired assumption / sensitivity parameter",
  "No Lee-specific reproductive efficiency estimate was identified.",
  "Retained as sensitivity parameter rather than claimed as directly measured.",

  "Egg dry mass", "m_egg", "2.25 ug baseline; varied in sensitivity analysis",
  "Reproductive energetic cost and protein burden",
  "Literature-informed assumption / sensitivity parameter",
  "Egg-size literature motivates uncertainty in reproductive-cost calculations, but exact Lee-specific egg dry mass was not available.",
  "Tests robustness of reproductive energetic burden and protein burden.",

  "Egg protein fraction", "q_egg", "0.50 baseline; varied in sensitivity analysis",
  "Protein-specific reproductive burden",
  "Sensitivity parameter",
  "No directly matched Lee-specific egg composition estimate was identified.",
  "Tests robustness of protein-specific reproductive burden.",

  "Protein-to-egg conversion efficiency", "eta_P_to_R", "0.85 baseline; varied in sensitivity analysis",
  "Protein-specific reproductive burden",
  "Sensitivity parameter",
  "No directly matched Lee-specific conversion efficiency estimate was identified.",
  "Tests robustness of protein-specific reproductive burden.",

  "Lipid fraction of adult dry mass", "q_L", "Diet-specific for three Lee et al. treatments",
  "Lipid-adjusted body-composition scenario",
  "Directly extracted from Lee et al.",
  "Lee et al. reported lipid content for P:C = 1:16, 1:4 and 1:2 at food concentration 180 g L^-1.",
  "Moved to main-text analysis as requested by Reviewer 1.",

  "Lipid energy density", "e_L", "39 J mg^-1",
  "Lipid-adjusted scenario",
  "Biochemical / energetic convention",
  "Used to contrast lipid-rich body composition with lean biomass assumptions.",
  "Recalculates body energy density and maintenance in the lipid-adjusted scenario.",

  "Lean biomass energy density", "e_lean", "17 J mg^-1",
  "Lipid-adjusted scenario",
  "Biochemical / energetic convention",
  "Used as simplified lean dry-mass energy density in the lipid-adjusted scenario.",
  "Recalculates body energy density and maintenance in the lipid-adjusted scenario."
)

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

lee_lipid_data <- tibble::tribble(
  ~Ratio,   ~Food_num, ~lipid_fraction_dry_mass, ~lipid_se,
  "(1:16)", 180,       0.300,                     0.011,
  "(1:4)",  180,       0.151,                     0.007,
  "(1:2)",  180,       0.140,                     0.009
)

e_lipid <- 39  # J / mg lipid
e_lean  <- 17  # J / mg lean dry mass

adult_dry_mass_mg <- adult_dry_mass_ug * 0.001

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
    lipid_mass_mg = adult_dry_mass_mg * lipid_fraction_dry_mass,
    lean_mass_mg  = adult_dry_mass_mg * (1 - lipid_fraction_dry_mass),

    body_energy_density_lipid_adjusted =
      lipid_fraction_dry_mass * e_lipid +
      (1 - lipid_fraction_dry_mass) * e_lean,

    adult_body_energy_baseline_J =
      adult_dry_mass_mg * biomass_energy,

    adult_body_energy_lipid_adjusted_J =
      adult_dry_mass_mg * body_energy_density_lipid_adjusted,

    maintenance_baseline_J_day =
      maintenance_J_day,

    maintenance_lipid_adjusted_J_day =
      maint_frac_adult_energy_per_day *
      adult_body_energy_lipid_adjusted_J,

    energetic_safety_margin_lipid_adjusted_J_day =
      assimilated_energy_J_day -
      reproductive_cost_J_day -
      maintenance_lipid_adjusted_J_day,

    energetic_safety_margin_fraction_lipid_adjusted =
      energetic_safety_margin_lipid_adjusted_J_day /
      assimilated_energy_J_day,

    energetically_feasible_lipid_adjusted =
      energetic_safety_margin_lipid_adjusted_J_day >= 0,

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

fig1_framework <- ggplot() +

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

annotate("rect", xmin = 12.2, xmax = 14.7, ymin = 3.5, ymax = 5.1,
         fill = "grey90", colour = "black") +
  annotate("text", x = 13.45, y = 4.3,
           label = "Age-specific\nmortality h(t)", size = 5) +

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

  geom_segment(aes(x = 2.7, xend = 8.0, y = 1.9, yend = 1.4),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.7) +
  geom_segment(aes(x = 2.7, xend = 8.0, y = 6.3, yend = 1.6),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.5, linetype = 2) +

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

p_fig2A <- p_fig2A + labs(x = NULL, y = NULL)
p_fig2B <- p_fig2B + labs(x = NULL, y = NULL)
p_fig2C <- p_fig2C + labs(x = NULL, y = NULL)
p_fig2D <- p_fig2D + labs(x = NULL, y = NULL)

fig2_inner <- (p_fig2A + p_fig2B) /
  (p_fig2C + p_fig2D) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 24, face = "bold"),
      plot.subtitle = element_text(size = 17),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

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
  file.path(output_dir, "Figure4.png"),
  fig4_lifespan_fit,
  width = 11,
  height = 5.5,
  dpi = 350
)

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
  tibble(
    parameter = "egg_dry_mass_ug",
    value = c(0.25, 0.50, 1.00, 2.25, 3.00, 5.00, 8.00)
  ),
  tibble(
    parameter = "adult_dry_mass_ug",
    value = c(100, 150, 200, 300, 400, 500, 750)
  ),
  tibble(
    parameter = "eta_P",
    value = c(0.10, 0.25, 0.50, 0.75, 0.95, 1.00)
  ),
  tibble(
    parameter = "eta_C",
    value = c(0.10, 0.25, 0.50, 0.75, 0.95, 1.00)
  ),
  tibble(
    parameter = "kap_R",
    value = c(0.40, 0.50, 0.70, 0.85, 0.95, 1.00)
  ),
  tibble(
    parameter = "maint_frac_adult_energy_per_day",
    value = c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50)
  ),
  tibble(
    parameter = "egg_protein_fraction_dry",
    value = c(0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80)
  ),
  tibble(
    parameter = "eta_protein_to_egg",
    value = c(0.30, 0.50, 0.70, 0.85, 0.95, 1.00)
  )
) %>%
  mutate(
    scenario_id = row_number(),

    egg_dry_mass_ug_sens =
      if_else(
        parameter == "egg_dry_mass_ug",
        value,
        baseline_parameters$egg_dry_mass_ug
      ),

    adult_dry_mass_ug_sens =
      if_else(
        parameter == "adult_dry_mass_ug",
        value,
        baseline_parameters$adult_dry_mass_ug
      ),

    eta_P_sens =
      if_else(
        parameter == "eta_P",
        value,
        baseline_parameters$eta_P
      ),

    eta_C_sens =
      if_else(
        parameter == "eta_C",
        value,
        baseline_parameters$eta_C
      ),

    kap_R_sens =
      if_else(
        parameter == "kap_R",
        value,
        baseline_parameters$kap_R
      ),

    maint_frac_sens =
      if_else(
        parameter == "maint_frac_adult_energy_per_day",
        value,
        baseline_parameters$maint_frac_adult_energy_per_day
      ),

    egg_protein_fraction_sens =
      if_else(
        parameter == "egg_protein_fraction_dry",
        value,
        baseline_parameters$egg_protein_fraction_dry
      ),

    eta_protein_to_egg_sens =
      if_else(
        parameter == "eta_protein_to_egg",
        value,
        baseline_parameters$eta_protein_to_egg
      )
  )

summarise_sensitivity_deb <- function(tmp) {

  tmp %>%
    summarise(
      n = n(),

      phi_R_median =
        median(reproductive_energetic_burden, na.rm = TRUE),
      phi_R_low =
        quantile(reproductive_energetic_burden, 0.025, na.rm = TRUE),
      phi_R_high =
        quantile(reproductive_energetic_burden, 0.975, na.rm = TRUE),

      phi_PR_median =
        median(reproductive_protein_burden, na.rm = TRUE),
      phi_PR_low =
        quantile(reproductive_protein_burden, 0.025, na.rm = TRUE),
      phi_PR_high =
        quantile(reproductive_protein_burden, 0.975, na.rm = TRUE),

      rho_R_median =
        median(reproductive_share_of_explicit_costs, na.rm = TRUE),
      rho_R_low =
        quantile(reproductive_share_of_explicit_costs, 0.025, na.rm = TRUE),
      rho_R_high =
        quantile(reproductive_share_of_explicit_costs, 0.975, na.rm = TRUE),

      safety_margin_fraction_median =
        median(energetic_safety_margin_fraction, na.rm = TRUE),
      safety_margin_fraction_low =
        quantile(energetic_safety_margin_fraction, 0.025, na.rm = TRUE),
      safety_margin_fraction_high =
        quantile(energetic_safety_margin_fraction, 0.975, na.rm = TRUE),

      feasible_fraction =
        mean(energetically_feasible, na.rm = TRUE),

      protein_feasible_fraction =
        mean(protein_feasible_for_reproduction, na.rm = TRUE),

      explicit_costs_exceed_assimilation_fraction =
        mean(explicit_cost_fraction_of_assimilation > 1, na.rm = TRUE),

      .groups = "drop"
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
      eta_protein_to_egg = scenario$eta_protein_to_egg_sens
    )
}

oat_sensitivity_summary <- map_dfr(
  seq_len(nrow(oat_sensitivity_grid)),
  run_oat_scenario
)

set.seed(123)

n_global_scenarios <- 1000

global_sensitivity_grid <- tibble(
  scenario_id = seq_len(n_global_scenarios),

  egg_dry_mass_ug =
    runif(n_global_scenarios, 0.5, 6.0),

  adult_dry_mass_ug =
    runif(n_global_scenarios, 100, 600),

  eta_P =
    runif(n_global_scenarios, 0.25, 1.00),

  eta_C =
    runif(n_global_scenarios, 0.25, 1.00),

  kap_R =
    runif(n_global_scenarios, 0.50, 1.00),

  egg_protein_fraction_dry =
    runif(n_global_scenarios, 0.25, 0.75),

  eta_protein_to_egg =
    runif(n_global_scenarios, 0.50, 1.00),

  maint_frac_adult_energy_per_day =
    exp(runif(n_global_scenarios, log(0.005), log(0.50)))
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
    maint_frac_adult_energy_per_day =
      scenario$maint_frac_adult_energy_per_day,
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
      maint_frac_adult_energy_per_day =
        scenario$maint_frac_adult_energy_per_day
    )
}

global_sensitivity_summary <- map_dfr(
  seq_len(nrow(global_sensitivity_grid)),
  run_global_scenario
)

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

figS5_oat_phiR <- ggplot(
  oat_sensitivity_summary,
  aes(x = value, y = phi_R_median)
) +
  geom_ribbon(aes(ymin = phi_R_low, ymax = phi_R_high), alpha = 0.25) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(
    ~ parameter,
    scales = "free_x",
    labeller = as_labeller(parameter_labels)
  ) +
  labs(
    x = "Parameter value",
    y = expression(phi[R] == R[d] / A[d]),
    title = "a"
  ) +
  theme_bw(base_size = 13)

figS5_oat_phiPR <- ggplot(
  oat_sensitivity_summary,
  aes(x = value, y = phi_PR_median)
) +
  geom_ribbon(aes(ymin = phi_PR_low, ymax = phi_PR_high), alpha = 0.25) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(
    ~ parameter,
    scales = "free_x",
    labeller = as_labeller(parameter_labels)
  ) +
  labs(
    x = "Parameter value",
    y = expression(phi[P*","*R]),
    title = "b"
  ) +
  theme_bw(base_size = 13)

figS5_oat_margin <- ggplot(
  oat_sensitivity_summary,
  aes(x = value, y = safety_margin_fraction_median)
) +
  geom_ribbon(
    aes(
      ymin = safety_margin_fraction_low,
      ymax = safety_margin_fraction_high
    ),
    alpha = 0.25
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(
    ~ parameter,
    scales = "free_x",
    labeller = as_labeller(parameter_labels)
  ) +
  labs(
    x = "Parameter value",
    y = expression(Delta / A[d]),
    title = "c"
  ) +
  theme_bw(base_size = 13)

figS5_oat_feasible <- ggplot(
  oat_sensitivity_summary,
  aes(x = value, y = feasible_fraction)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(
    ~ parameter,
    scales = "free_x",
    labeller = as_labeller(parameter_labels)
  ) +
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

figS7_global_maintenance <- ggplot(
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
    title = ""
  ) +
  theme_bw(base_size = 14)

figS7_global_sensitivity <- (figS7_global_phiR + figS7_global_margin) /
  (figS7_global_feasible + figS7_global_maintenance)

ggsave(
  file.path(output_dir, "FigureS3.png"),
  figS7_global_sensitivity,
  width = 15,
  height = 9,
  dpi = 350
)

ggsave(
  filename = file.path(output_dir, "Figure5.png"),
  plot = figS_lipid_adjusted_scenario,
  width = 12,
  height = 7,
  dpi = 600
)

eta_scenario_values <- c(0.50, 0.65, 0.75, 0.90)

unequal_assimilation_grid <- tidyr::crossing(
  eta_P_scenario = eta_scenario_values,
  eta_C_scenario = eta_scenario_values
) %>%
  mutate(
    scenario_id = row_number(),
    scenario_name = paste0(
      "etaP_", eta_P_scenario,
      "_etaC_", eta_C_scenario
    ),
    scenario_class = case_when(
      eta_P_scenario == eta_P & eta_C_scenario == eta_C ~ "baseline_equal",
      eta_P_scenario == eta_C_scenario ~ "equal_nonbaseline",
      TRUE ~ "unequal"
    )
  )

run_unequal_assimilation_scenario <- function(row_i) {

  scenario <- unequal_assimilation_grid[row_i, ]

  tmp <- calculate_individual_deb(
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
    eta_P = scenario$eta_P_scenario,
    eta_C = scenario$eta_C_scenario,
    kap_R = kap_R,
    maint_frac_adult_energy_per_day = maint_frac_adult_energy_per_day,
    egg_protein_fraction_dry = egg_protein_fraction_dry,
    eta_protein_to_egg = eta_protein_to_egg
  )

  tmp %>%
    summarise_sensitivity_deb() %>%
    mutate(
      scenario_id = scenario$scenario_id,
      scenario_name = scenario$scenario_name,
      scenario_class = scenario$scenario_class,
      eta_P = scenario$eta_P_scenario,
      eta_C = scenario$eta_C_scenario,
      eta_difference = eta_P - eta_C,
      eta_ratio = eta_P / eta_C
    )
}

unequal_assimilation_summary <- purrr::map_dfr(
  seq_len(nrow(unequal_assimilation_grid)),
  run_unequal_assimilation_scenario
)

unequal_assimilation_treatment_summary <- purrr::map_dfr(
  seq_len(nrow(unequal_assimilation_grid)),
  function(row_i) {

    scenario <- unequal_assimilation_grid[row_i, ]

    calculate_individual_deb(
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
      eta_P = scenario$eta_P_scenario,
      eta_C = scenario$eta_C_scenario,
      kap_R = kap_R,
      maint_frac_adult_energy_per_day = maint_frac_adult_energy_per_day,
      egg_protein_fraction_dry = egg_protein_fraction_dry,
      eta_protein_to_egg = eta_protein_to_egg
    ) %>%
      group_by(treatment, Ratio, Food) %>%
      summarise(
        n = n(),
        P_daily_ug = mean(P_daily_ug, na.rm = TRUE),
        C_daily_ug = mean(C_daily_ug, na.rm = TRUE),
        lifespan = mean(lifespan, na.rm = TRUE),
        phi_R_median = median(reproductive_energetic_burden, na.rm = TRUE),
        phi_P_median = median(reproductive_protein_burden, na.rm = TRUE),
        safety_margin_median =
          median(energetic_safety_margin_fraction, na.rm = TRUE),
        feasible_fraction =
          mean(energetically_feasible, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        scenario_id = scenario$scenario_id,
        scenario_name = scenario$scenario_name,
        scenario_class = scenario$scenario_class,
        eta_P = scenario$eta_P_scenario,
        eta_C = scenario$eta_C_scenario
      )
  }
)

fig_revision_unequal_assimilation <- unequal_assimilation_summary %>%
  mutate(
    eta_P_plot = factor(eta_P, levels = eta_scenario_values),
    eta_C_plot = factor(eta_C, levels = eta_scenario_values),
    feasible_label = formatC(feasible_fraction, format = "f", digits = 2)
  ) %>%
  ggplot(aes(x = eta_P_plot, y = eta_C_plot)) +
  geom_tile(
    aes(fill = safety_margin_fraction_median),
    colour = "white",
    linewidth = 0.4
  ) +
  geom_tile(
    data = . %>% filter(scenario_class == "baseline_equal"),
    aes(x = eta_P_plot, y = eta_C_plot),
    fill = NA,
    colour = "black",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(label = feasible_label),
    size = 4
  ) +
  scale_fill_viridis_c(option = "D") +
  labs(
    x = expression(eta[P]),
    y = expression(eta[C]),
    fill = expression("Median " * Delta / A[d])
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

penalty_strength_values <- c(0.25, 0.50, 1.00, 2.00)

maintenance_penalty_grid <- bind_rows(
  tibble(
    penalty_type = "none_baseline",
    penalty_strength = 0,
    lambda_D = 0,
    lambda_P = 0
  ),
  tibble(
    penalty_type = "distance_from_target",
    penalty_strength = penalty_strength_values,
    lambda_D = penalty_strength_values,
    lambda_P = 0
  ),
  tibble(
    penalty_type = "protein_excess",
    penalty_strength = penalty_strength_values,
    lambda_D = 0,
    lambda_P = penalty_strength_values
  ),
  tibble(
    penalty_type = "combined_distance_and_protein_excess",
    penalty_strength = penalty_strength_values,
    lambda_D = penalty_strength_values,
    lambda_P = penalty_strength_values
  )
) %>%
  mutate(
    scenario_id = row_number(),
    scenario_name = paste0(
      penalty_type,
      "_lambdaD_", lambda_D,
      "_lambdaP_", lambda_P
    )
  )

apply_maintenance_penalty <- function(data, lambda_D = 0, lambda_P = 0) {

  target_ratio_PC <- target_P_ug / target_C_ug

  data %>%
    mutate(
      target_P_given_C_ug = C_daily_ug * target_ratio_PC,

      protein_excess_fraction = case_when(
        target_P_given_C_ug > 0 ~
          pmax(0, (P_daily_ug - target_P_given_C_ug) / target_P_given_C_ug),
        target_P_given_C_ug <= 0 & P_daily_ug > 0 ~ 1,
        TRUE ~ 0
      ),

      imbalance_penalty_J_day =
        lambda_D * maintenance_J_day * distance_from_target_scaled,

      protein_excess_penalty_J_day =
        lambda_P * maintenance_J_day * protein_excess_fraction,

      maintenance_penalty_J_day =
        imbalance_penalty_J_day + protein_excess_penalty_J_day,

      maintenance_penalised_J_day =
        maintenance_J_day + maintenance_penalty_J_day,

      energetic_safety_margin_penalised_J_day =
        assimilated_energy_J_day -
        reproductive_cost_J_day -
        maintenance_penalised_J_day,

      energetic_safety_margin_penalised_fraction =
        energetic_safety_margin_penalised_J_day / assimilated_energy_J_day,

      energetically_feasible_penalised =
        energetic_safety_margin_penalised_J_day >= 0,

      explicit_cost_fraction_penalised =
        (reproductive_cost_J_day + maintenance_penalised_J_day) /
        assimilated_energy_J_day,

      reproductive_share_penalised_explicit_costs =
        reproductive_cost_J_day /
        (reproductive_cost_J_day + maintenance_penalised_J_day)
    )
}

maintenance_penalty_summary <- purrr::map_dfr(
  seq_len(nrow(maintenance_penalty_grid)),
  function(row_i) {

    scenario <- maintenance_penalty_grid[row_i, ]

    tmp <- apply_maintenance_penalty(
      deb_individual,
      lambda_D = scenario$lambda_D,
      lambda_P = scenario$lambda_P
    )

    tmp %>%
      summarise(
        n = n(),

        median_original_margin =
          median(energetic_safety_margin_fraction, na.rm = TRUE),

        median_penalised_margin =
          median(energetic_safety_margin_penalised_fraction, na.rm = TRUE),

        feasible_fraction_original =
          mean(energetically_feasible, na.rm = TRUE),

        feasible_fraction_penalised =
          mean(energetically_feasible_penalised, na.rm = TRUE),

        median_original_maintenance =
          median(maintenance_J_day, na.rm = TRUE),

        median_penalised_maintenance =
          median(maintenance_penalised_J_day, na.rm = TRUE),

        median_maintenance_penalty =
          median(maintenance_penalty_J_day, na.rm = TRUE),

        median_reproductive_share_penalised =
          median(reproductive_share_penalised_explicit_costs, na.rm = TRUE),

        .groups = "drop"
      ) %>%
      mutate(
        scenario_id = scenario$scenario_id,
        scenario_name = scenario$scenario_name,
        penalty_type = scenario$penalty_type,
        penalty_strength = scenario$penalty_strength,
        lambda_D = scenario$lambda_D,
        lambda_P = scenario$lambda_P
      )
  }
)

maintenance_penalty_treatment_summary <- purrr::map_dfr(
  seq_len(nrow(maintenance_penalty_grid)),
  function(row_i) {

    scenario <- maintenance_penalty_grid[row_i, ]

    tmp <- apply_maintenance_penalty(
      deb_individual,
      lambda_D = scenario$lambda_D,
      lambda_P = scenario$lambda_P
    )

    tmp %>%
      group_by(treatment, Ratio, Food) %>%
      summarise(
        n = n(),
        P_daily_ug = mean(P_daily_ug, na.rm = TRUE),
        C_daily_ug = mean(C_daily_ug, na.rm = TRUE),
        lifespan = mean(lifespan, na.rm = TRUE),

        median_penalised_margin =
          median(energetic_safety_margin_penalised_fraction, na.rm = TRUE),

        feasible_fraction_penalised =
          mean(energetically_feasible_penalised, na.rm = TRUE),

        median_maintenance_penalty =
          median(maintenance_penalty_J_day, na.rm = TRUE),

        .groups = "drop"
      ) %>%
      mutate(
        scenario_id = scenario$scenario_id,
        scenario_name = scenario$scenario_name,
        penalty_type = scenario$penalty_type,
        penalty_strength = scenario$penalty_strength,
        lambda_D = scenario$lambda_D,
        lambda_P = scenario$lambda_P
      )
  }
)

maintenance_penalty_plot_data <- maintenance_penalty_summary %>%
  filter(penalty_type != "none_baseline") %>%
  mutate(
    penalty_label = dplyr::recode(
      penalty_type,
      distance_from_target = "Distance-from-target penalty",
      protein_excess = "Protein-excess penalty",
      combined_distance_and_protein_excess = "Combined penalty"
    )
  )

baseline_penalty_plot_data <- maintenance_penalty_summary %>%
  filter(penalty_type == "none_baseline") %>%
  tidyr::crossing(
    penalty_label = c(
      "Distance-from-target penalty",
      "Protein-excess penalty",
      "Combined penalty"
    )
  )

maintenance_penalty_plot_data <- bind_rows(
  baseline_penalty_plot_data,
  maintenance_penalty_plot_data
) %>%
  mutate(
    penalty_label = factor(
      penalty_label,
      levels = c(
        "Distance-from-target penalty",
        "Protein-excess penalty",
        "Combined penalty"
      )
    )
  )

fig_revision_maintenance_penalty <- ggplot(
  maintenance_penalty_plot_data,
  aes(x = penalty_strength, y = median_penalised_margin)
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(
    aes(group = penalty_label),
    alpha = 0.65,
    linewidth = 0.8
  ) +
  geom_point(
    aes(size = feasible_fraction_penalised),
    alpha = 0.85
  ) +
  facet_wrap(~ penalty_label, scales = "free_x") +
  scale_size_continuous(
    limits = c(0, 1),
    range = c(2, 6)
  ) +
  labs(
    x = "Penalty strength",
    y = expression("Median penalised safety margin " * Delta / A[d]),
    size = "Feasible fraction"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

quantity_quality_summary <- treatment_observed_predicted_lifespan %>%
  mutate(
    Food_num = readr::parse_number(as.character(Food)),
    Ratio_chr = as.character(Ratio),

    quality_distance = ifelse(
      P_daily_ug > 0 & C_daily_ug > 0,
      abs(log2((P_daily_ug / C_daily_ug) / (target_P_ug / target_C_ug))),
      NA_real_
    ),

    log2_food = log2(Food_num)
  ) %>%
  group_by(Ratio_chr) %>%
  arrange(Food_num, .by_group = TRUE) %>%
  mutate(
    observed_lifespan_change_from_lowest_food =
      observed_lifespan - first(observed_lifespan),

    predicted_lifespan_change_from_lowest_food =
      predicted_lifespan - first(predicted_lifespan),

    observed_lifespan_rank_within_rail =
      rank(-observed_lifespan, ties.method = "average"),

    predicted_lifespan_rank_within_rail =
      rank(-predicted_lifespan, ties.method = "average")
  ) %>%
  ungroup()

quantity_quality_slopes <- quantity_quality_summary %>%
  group_by(Ratio_chr) %>%
  summarise(
    n_food_levels = n_distinct(Food_num),
    min_food = min(Food_num, na.rm = TRUE),
    max_food = max(Food_num, na.rm = TRUE),

    observed_lifespan_lowest_food =
      observed_lifespan[which.min(Food_num)][1],

    observed_lifespan_highest_food =
      observed_lifespan[which.max(Food_num)][1],

    predicted_lifespan_lowest_food =
      predicted_lifespan[which.min(Food_num)][1],

    predicted_lifespan_highest_food =
      predicted_lifespan[which.max(Food_num)][1],

    observed_quantity_effect =
      observed_lifespan_highest_food - observed_lifespan_lowest_food,

    predicted_quantity_effect =
      predicted_lifespan_highest_food - predicted_lifespan_lowest_food,

    observed_slope_per_log2_food =
      coef(lm(observed_lifespan ~ log2_food))[2],

    predicted_slope_per_log2_food =
      coef(lm(predicted_lifespan ~ log2_food))[2],

    mean_quality_distance =
      mean(quality_distance, na.rm = TRUE),

    .groups = "drop"
  )

near_target_reference <- quantity_quality_summary %>%
  filter(Ratio_chr == "(1:4)") %>%
  summarise(
    max_near_target_observed = max(observed_lifespan, na.rm = TRUE),
    max_near_target_predicted = max(predicted_lifespan, na.rm = TRUE),
    min_near_target_observed = min(observed_lifespan, na.rm = TRUE),
    min_near_target_predicted = min(predicted_lifespan, na.rm = TRUE)
  )

quantity_overcomes_quality_table <- quantity_quality_summary %>%
  mutate(
    exceeds_best_near_target_observed =
      observed_lifespan >= near_target_reference$max_near_target_observed,

    exceeds_best_near_target_predicted =
      predicted_lifespan >= near_target_reference$max_near_target_predicted,

    exceeds_low_near_target_observed =
      observed_lifespan >= near_target_reference$min_near_target_observed,

    exceeds_low_near_target_predicted =
      predicted_lifespan >= near_target_reference$min_near_target_predicted
  )

quantity_quality_plot_data <- quantity_quality_summary %>%
  mutate(
    Ratio_plot = factor(Ratio_chr, levels = ratio_order),
    Food_num = as.numeric(Food_num)
  ) %>%
  select(
    Ratio_plot,
    Food_num,
    observed_lifespan,
    predicted_lifespan
  ) %>%
  pivot_longer(
    cols = c(observed_lifespan, predicted_lifespan),
    names_to = "lifespan_source",
    values_to = "lifespan_days"
  ) %>%
  mutate(
    lifespan_source = dplyr::recode(
      lifespan_source,
      observed_lifespan = "Observed",
      predicted_lifespan = "Predicted"
    ),
    lifespan_source = factor(
      lifespan_source,
      levels = c("Observed", "Predicted")
    )
  )

ratio_order_revision <- c("(1:16)", "(1:8)", "(1:4)", "(1:2)", "(1:1)", "(2:1)")

quantity_quality_plot_data <- quantity_quality_summary %>%
  mutate(
    Ratio_plot = factor(Ratio_chr, levels = ratio_order_revision),
    Food_num = as.numeric(Food_num)
  ) %>%
  filter(!is.na(Ratio_plot)) %>%
  select(
    Ratio_plot,
    Food_num,
    observed_lifespan,
    predicted_lifespan
  ) %>%
  pivot_longer(
    cols = c(observed_lifespan, predicted_lifespan),
    names_to = "lifespan_source",
    values_to = "lifespan_days"
  ) %>%
  mutate(
    lifespan_source = dplyr::recode(
      lifespan_source,
      observed_lifespan = "Observed",
      predicted_lifespan = "Predicted"
    ),
    lifespan_source = factor(
      lifespan_source,
      levels = c("Observed", "Predicted")
    )
  )

fig_revision_quantity_quality <- ggplot(
  quantity_quality_plot_data,
  aes(
    x = Food_num,
    y = lifespan_days,
    colour = Ratio_plot,
    group = interaction(Ratio_plot, lifespan_source),
    linetype = lifespan_source
  )
) +
  facet_wrap(~ Ratio_plot, ncol = 3) +
  geom_line(linewidth = 0.9) +
  geom_point(
    data = quantity_quality_plot_data %>%
      filter(lifespan_source == "Observed"),
    aes(
      x = Food_num,
      y = lifespan_days,
      colour = Ratio_plot
    ),
    inherit.aes = FALSE,
    size = 2.4
  ) +
  scale_x_continuous(
    trans = "log2",
    breaks = sort(unique(quantity_quality_summary$Food_num))
  ) +
  scale_colour_viridis_d(option = "D", drop = FALSE) +
  scale_linetype_manual(
    values = c(
      Observed = "solid",
      Predicted = "dashed"
    )
  ) +
  labs(
    x = expression("Food concentration (g L"^{-1}*")"),
    y = "Treatment-level lifespan (days)",
    colour = "P:C rail",
    linetype = "Lifespan"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

fig_revision_quantity_quality

revision_scenario_definitions <- bind_rows(
  tibble(
    analysis = "Baseline energetic accounting",
    scenario = "baseline",
    parameters_varied = "None",
    literature_status =
      "Uses Lee et al. intake, lifespan and reproduction; other quantities are baseline assumptions tested by sensitivity analysis.",
    reviewer_issue_addressed =
      "Original model reproduction and comparison baseline."
  ),
  tibble(
    analysis = "Unequal assimilation",
    scenario = "eta_P x eta_C grid",
    parameters_varied =
      "eta_P and eta_C varied independently across 0.50, 0.65, 0.75 and 0.90.",
    literature_status =
      "Literature-informed scenario; no directly transferable Lee-matched eta_P or eta_C values identified.",
    reviewer_issue_addressed =
      "Addresses criticism that eta_P = eta_C makes protein and carbohydrate energetically interchangeable."
  ),
  tibble(
    analysis = "Diet-imbalance maintenance penalty",
    scenario = "M* = M + lambda_D M D_scaled",
    parameters_varied =
      "lambda_D varied across 0.25, 0.50, 1.00 and 2.00.",
    literature_status =
      "Exploratory scenario motivated by diet-imbalance literature; not a directly fitted DEB maintenance flux.",
    reviewer_issue_addressed =
      "Addresses criticism that maintenance is diet-independent."
  ),
  tibble(
    analysis = "Protein-excess maintenance penalty",
    scenario =
      "M* = M + lambda_P M max(0, (P - P_target_given_C) / P_target_given_C)",
    parameters_varied =
      "lambda_P varied across 0.25, 0.50, 1.00 and 2.00.",
    literature_status =
      "Exploratory scenario motivated by high-protein and geometric-stoichiometry literature; protein excess is calculated relative to the self-selected optimal P:C target rail.",
    reviewer_issue_addressed =
      "Addresses criticism that protein excess and nutrient imbalance are not explicitly represented."
  ),
  tibble(
    analysis = "Combined maintenance penalty",
    scenario =
      "M* = M + lambda_D M D_scaled + lambda_P M max(0, (P - P_target_given_C) / P_target_given_C)",
    parameters_varied =
      "lambda_D and lambda_P varied together across 0.25, 0.50, 1.00 and 2.00.",
    literature_status =
      "Exploratory combined scenario testing whether simultaneous diet-imbalance and protein-excess penalties alter qualitative conclusions.",
    reviewer_issue_addressed =
      "Addresses the combined consequences of nutrient imbalance and excess protein."
  ),
  tibble(
    analysis = "Lipid-adjusted body composition",
    scenario = "Lee et al. lipid subset",
    parameters_varied =
      "q_L uses Lee et al. values for P:C = 1:16, 1:4 and 1:2 at 180 g L^-1.",
    literature_status =
      "Directly extracted from Lee et al.; moved to main text in revision.",
    reviewer_issue_addressed =
      "Addresses request to include lipid modelling in the main text."
  ),
  tibble(
    analysis = "Quantity-versus-quality compensation",
    scenario = "Treatment-level food-concentration analysis within P:C rails",
    parameters_varied =
      "Food concentration and P:C rail; no new energetic parameter.",
    literature_status =
      "Uses the original Lee et al. experimental design and model predictions. The self-selected optimal P:C target is not assumed to be the lifespan-maximising rail.",
    reviewer_issue_addressed =
      "Addresses request to investigate whether quantity can compensate for poor nutritional quality."
  )
)

write.csv(
  parameter_provenance,
  file.path(output_dir, "Table_revision_parameter_provenance.csv"),
  row.names = FALSE
)

write.csv(
  revision_scenario_definitions,
  file.path(output_dir, "Table_revision_scenario_definitions.csv"),
  row.names = FALSE
)

write.csv(
  unequal_assimilation_grid,
  file.path(output_dir, "Table_revision_unequal_assimilation_grid.csv"),
  row.names = FALSE
)

write.csv(
  unequal_assimilation_summary,
  file.path(output_dir, "Table_revision_unequal_assimilation_summary.csv"),
  row.names = FALSE
)

write.csv(
  unequal_assimilation_treatment_summary,
  file.path(output_dir, "Table_revision_unequal_assimilation_treatment_summary.csv"),
  row.names = FALSE
)

write.csv(
  maintenance_penalty_grid,
  file.path(output_dir, "Table_revision_maintenance_penalty_grid.csv"),
  row.names = FALSE
)

write.csv(
  maintenance_penalty_summary,
  file.path(output_dir, "Table_revision_maintenance_penalty_summary.csv"),
  row.names = FALSE
)

write.csv(
  maintenance_penalty_treatment_summary,
  file.path(output_dir, "Table_revision_maintenance_penalty_treatment_summary.csv"),
  row.names = FALSE
)

write.csv(
  quantity_quality_summary,
  file.path(output_dir, "Table_revision_quantity_quality_summary.csv"),
  row.names = FALSE
)

write.csv(
  quantity_quality_slopes,
  file.path(output_dir, "Table_revision_quantity_quality_slopes.csv"),
  row.names = FALSE
)

write.csv(
  quantity_overcomes_quality_table,
  file.path(output_dir, "Table_revision_quantity_overcomes_quality.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(output_dir, "Figure1.png"),
  plot = fig1_framework,
  width = 13,
  height = 7,
  dpi = 350
)

ggsave(
  filename = file.path(output_dir, "Figure5.png"),
  plot = figS_lipid_adjusted_scenario,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "FigureS4.png"),
  plot = fig_revision_unequal_assimilation,
  width = 8,
  height = 6,
  dpi = 350
)

ggsave(
  filename = file.path(output_dir, "FigureS5.png"),
  plot = fig_revision_maintenance_penalty,
  width = 12,
  height = 6,
  dpi = 350
)

ggsave(
  filename = file.path(output_dir, "FigureS6.png"),
  plot = fig_revision_quantity_quality,
  width = 11,
  height = 6.5,
  dpi = 350
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo_RSIF_revision.txt")
)
