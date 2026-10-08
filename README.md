# Estimating missing CWT recoveries in unsampled strata — North of Falcon Chinook Salmon troll fishery

Code and data supporting:

> Kvistad, J. T. (in review). Estimating missing coded-wire tag recoveries from unsampled strata in a mixed-stock Chinook Salmon troll fishery under a non-random missingness structure. *North American Journal of Fisheries Management*.

The analysis compares two estimators of coded-wire tag (CWT) recoveries in unsampled time-area strata of the North of Falcon (Washington Marine Areas 1–4, 4B) Chinook Salmon troll fishery, 2005–2023:

1. **Tag density model** — a Tweedie GAM of CWT recoveries with sampling effort as an offset, applied to landed catch in unsampled strata, with stock and age composition drawn from a Dirichlet posterior.
2. **Ratio estimator** — nearest-donor borrowing of tag density and composition from the closest sampled stratum within 90 days.

Both are evaluated operationally and in a 1,000-replicate hold-out.

## Repository structure

```
.
├── analysis/
│   └── NOF_troll_analysis.qmd      # Main analysis (data prep, GAM, composition, ratio estimator, figures)
├── R/
│   ├── method2_uncertainty.R       # Posterior simulation and bootstrap intervals (Figure 7, Figure S1)
│   ├── method2_prior_sensitivity.R # Dirichlet prior sensitivity (Table S1)
│   ├── method_holdout_combined.R   # Hold-out evaluation (Figure 10)
│   ├── 00_export_public_data.R     # Author-only: builds data/derived/ from the PSC CAMP database
│   └── run_all.R                   # Runs everything in order
├── data/
│   ├── raw/                        # RMIS extracts (public)
│   ├── derived/                    # Stock lookup derived from CAMP
│   └── README.md                   # Data sources and dictionary
├── out/                            # Figures and tables written here (not tracked)
├── renv.lock                       # Package versions
├── CITATION.cff
└── LICENSE
```

## Reproducing the analysis

Requirements: R (≥ 4.3), Quarto, and the packages pinned in `renv.lock`.

```r
# 1. Clone or download this repository and open it as the working directory.
# 2. Restore the package library:
install.packages("renv")
renv::restore()

# 3. Run the full analysis:
source("R/run_all.R")
```

`run_all.R` extracts the R code from the Quarto document, runs it, and then runs the uncertainty, prior sensitivity, and hold-out scripts, which depend on objects created by the main analysis. Figures and tables are written to `out/`. The hold-out (1,000 replicates, each refitting the GAM) takes the longest; set `N_REP` lower in `R/method_holdout_combined.R` for a quick check.

Random seeds are fixed (`set.seed(42)`), so results should reproduce exactly given the package versions in `renv.lock`.

## Data

Catch, sampling, recovery, and release data are public extracts from the Regional Mark Information System (RMIS). Stock assignments for PST indicator stock tag codes come from the Chinook Analysis Modeling Platform (CAMP), an internal database of the Pacific Salmon Commission (PSC). Only a derived lookup of tag code, stock, and brood year for tag codes recovered in this fishery is included here; the CAMP database itself is not public. See [`data/README.md`](data/README.md).

## License

Code is released under the MIT License (see `LICENSE`). Data in `data/raw/` are RMIS extracts and remain subject to RMIS terms of use; cite RMIS as the original source. The derived stock lookup in `data/derived/` is provided for reproduction of this analysis.

## Contact

Jake Kvistad, Washington Department of Fish and Wildlife — jake.kvistad@wa.dfw.gov
