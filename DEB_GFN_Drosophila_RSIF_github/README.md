# DEB-GFN Drosophila analysis

This repository contains the analysis code and input data for the manuscript:

**A DEB-inspired bioenergetic model linking nutritional geometry to age-specific mortality**

The analysis links adult protein and carbohydrate intake, egg production and lifespan in a Drosophila nutritional-geometry dataset to DEB-inspired energetic state variables and age-specific mortality models. It also includes scenario and sensitivity analyses for assimilation, maintenance penalties, lipid-adjusted body composition and quantity-versus-quality compensation.

## Repository structure

```text
DEB_GFN_Drosophila_RSIF_github/
├── R/
│   └── run_analysis.R
├── data/
│   └── Dataset_Lee_Nutrigonometry.csv
├── outputs/
├── data_dictionary.csv
├── metadata.yml
├── README.md
└── .gitignore
```

## Requirements

The analysis was written for R and requires the following packages:

```r
tidyverse, mgcv, glmmTMB, coxme, survival, splines, performance, DHARMa,
patchwork, viridis, scales, broom, broom.mixed, ggforce, car, DescTools, cowplot
```

The script checks for required packages and stops with a list of missing packages. Install missing packages before running the analysis.

## How to run

From the repository root:

```bash
Rscript R/run_analysis.R
```

All figures, tables and session information are written to `outputs/`.

## Main outputs

The script generates:

- `Figure1.png` — conceptual DEB-GFN framework
- `Figure2.png` — nutritional surfaces and derived energetic state variables
- `Figure3.png` — predicted age-specific mortality and survival trajectories
- `Figure4.png` — observed versus predicted treatment-level lifespan
- `Figure5.png` — lipid-adjusted body-composition scenario
- `FigureS1.png` to `FigureS6.png` — supplementary figures
- CSV tables for parameter provenance, scenario definitions, unequal assimilation, maintenance penalties and quantity-versus-quality analyses
- `sessionInfo_RSIF_revision.txt` for reproducibility

## Input data

`data/Dataset_Lee_Nutrigonometry.csv` contains individual-level observations with protein intake, carbohydrate intake, lifespan, lifetime egg production, daily egg production, treatment identifiers, P:C rail and food concentration.

The data derive from the Drosophila nutritional-geometry study by Lee et al. (2008). Please cite the original data source as appropriate when reusing the dataset.

## Reproducibility notes

The analysis uses the observed adult-stage nutritional-geometry data and constructs DEB-inspired state variables. It is not a full dynamic DEB parameterisation because the dataset does not contain time-resolved reserve dynamics, excretion, metabolic-rate measurements or daily reproductive allocation.
