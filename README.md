# DEB-inspired bioenergetic model linking nutritional geometry to age-specific mortality

This repository contains the final consolidated R script used to reproduce the analyses and figures for the manuscript:

**A DEB-inspired bioenergetic model linking nutritional geometry to age-specific mortality**


## Files

| File | Description |
|---|---|
| `Morimoto-DEBinspired_git.R` | Final manuscript-only R script. Runs the analysis and generates all figures/tables used in the main manuscript and supplementary material. |
| `Dataset_Lee_Nutrigonometry.csv` | Input dataset required by the script. This file must be in the working directory before running the script. |
| `DEB_GFN_Drosophila_final_outputs/` | Output folder created automatically by the script. Contains figures, model summaries, and tables. |

## What the script does

The consolidated script performs only the analyses retained in the final manuscript. Specifically, it:

1. reads and cleans the Lee et al. adult *Drosophila* nutritional geometry dataset;
2. calculates the nutritional target point on the observed `P:C = 1:4` rail;
3. derives DEB-inspired energetic variables, including:
   - assimilated energy;
   - maintenance expenditure;
   - reproductive energetic cost;
   - energetic safety margin;
   - reproductive energetic burden;
   - reproductive share of explicit energetic costs;
   - protein-specific reproductive burden;
   - nutritional displacement from the intake target;
4. fits the final discrete-time survival model with complementary log-log link;
5. compares candidate survival models using AIC and likelihood-ratio tests;
6. evaluates model terms using Type-II analysis of deviance;
7. checks collinearity in the final model;
8. generates treatment-level mortality and survival predictions;
9. calculates concordance between observed and predicted mean and median lifespan;
10. performs one-at-a-time and global sensitivity analyses;
11. performs the lipid-adjusted energetic scenario for the three treatments with available lipid data;
12. saves all manuscript figures and relevant model/tables outputs.

Exploratory analyses not retained in the final manuscript, such as Cox diagnostics, DHARMa diagnostics, Add-my-Pet anchored sensitivity analyses, target-point sensitivity analyses, representative profile predictions, and the protein-burden survival model, were removed from this final script.

## Software requirements

The analysis was developed using:

- R version 4.5.0 (2025-04-11)
- Platform: `aarch64-apple-darwin20`
- Running under: macOS 26.3.1

The script installs missing packages automatically if they are not already available.

Core R packages used by the final consolidated script are:

| Package | Purpose |
|---|---|
| `tidyverse` | Data wrangling, reshaping, and plotting infrastructure |
| `mgcv` | Smooth surface fitting for nutritional landscapes |
| `glmmTMB` | Discrete-time mixed-effects survival modelling |
| `splines` | Natural cubic splines for age-dependent effects |
| `performance` | Collinearity and model checks |
| `patchwork` | Figure assembly |
| `viridis` | Colour scales |
| `scales` | Plot scale handling |
| `broom.mixed` | Tidying mixed-model outputs |
| `car` | Type-II analysis of deviance |
| `DescTools` | Lin's concordance correlation coefficient |
| `cowplot` | Shared figure axis labels |

The package versions used during manuscript preparation were:

```text
R 4.5.0
DescTools 0.99.60
glmmTMB 1.1.12
car 3.1-3
performance 0.16.0
mgcv 1.9-3
patchwork 1.3.0
viridis 0.6.5
viridisLite 0.4.2
scales 1.4.0
broom.mixed 0.2.9.7
broom 1.0.8
ggplot2 3.5.2
tidyverse 2.0.0
dplyr 1.2.1
tidyr 1.3.1
readr 2.1.5
tibble 3.3.0
purrr 1.2.2
stringr 1.5.1
forcats 1.0.0
cowplot 1.1.3
```

## How to run

1. Place R file and `Dataset_Lee_Nutrigonometry.csv` in the same working directory.
2. Open R or RStudio.
3. Set the working directory to the folder containing both files.
4. Run:

```r
source("Morimoto-DEBinspired_git.R")
```

The script will create the output directory:

```text
DEB_GFN_Drosophila_final_outputs/
```

All tables, model summaries, and figures will be saved there.

## Main outputs

The script generates the following main manuscript figures:

| Output file | Manuscript figure |
|---|---|
| `Figure1_maintext.png` | Figure 1. Conceptual DEB-GFN framework |
| `Figure2_maintext.png` | Figure 2. DEB-inspired energetic quantities across nutrient space |
| `Figure3_maintext.png` | Figure 3. Predicted mortality and survival trajectories |
| `Figure4_maintext.png` | Figure 4. Observed versus predicted treatment-level lifespan |

It also generates the following supplementary figures:

| Output file | Supplementary figure |
|---|---|
| `FigureS1_observed_vs_predicted_survival_by_rail_food.png` | Figure S1. Observed versus predicted survival curves |
| `FigureS2.png` | Figure S2. One-at-a-time parameter sensitivity |
| `FigureS3.png` | Figure S3. Global sensitivity analysis |
| `FigureS4.png` | Figure S4. Lipid-adjusted energetic scenario |

## Main model

The final survival model is a discrete-time mixed-effects model with complementary log-log link, fitted using `glmmTMB`. Adult age is modelled using natural cubic splines. The final model includes:

- nutritional rail;
- food concentration;
- daily protein intake;
- daily carbohydrate intake;
- energetic safety margin;
- nutritional displacement from the `P:C = 1:4` target;
- reproductive share of explicit energetic costs;
- age-dependent interactions for selected nutritional predictors;
- treatment-level random intercept.

The model is intended as an adult-only, DEB-inspired energetic modelling framework. It is not a full DEB parameterisation and does not reconstruct larval development, reserves, maturity, or individual time-resolved oviposition dynamics.

## Notes on reproductive burden

In this dataset, reproductive output is available as average daily egg production rather than as a full individual daily oviposition time series. Therefore, reproductive burden variables should be interpreted as average adult-life reproductive burdens. The model tests whether mean reproductive demand is associated with age-specific mortality, but it does not reconstruct short-term changes in mortality risk following daily changes in egg production.

## Lipid-adjusted scenario

Lee et al. reported lipid content for only three diet treatments:

- `P:C = 1:16`, food concentration 180 g L^-1;
- `P:C = 1:4`, food concentration 180 g L^-1;
- `P:C = 1:2`, food concentration 180 g L^-1.

The script uses these data in a restricted exploratory scenario to recalculate body energy density, maintenance expenditure, energetic safety margin, and energetic feasibility. This analysis is not used to construct a full lipid landscape.

## Reproducibility notes

The script is designed so that all output-writing operations occur at the end of the file where possible. Some package installation may vary by system. Numerical results should be checked if major package or R versions differ from those listed above.

## Citation and acknowledgement

The analysis uses the adult *Drosophila melanogaster* nutritional geometry dataset from Lee et al. The raw dataset was kindly provided by Prof Kwang P. Lee and is acknowledged in the manuscript.

## AI assistance disclosure

A generative AI tool (ChatGPT-5.5 Thinking, OpenAI) was used to assist with grammar, wording, code readability, and code troubleshooting. Generative AI was used to create the README file.
