# run_all.R
#
# Reproduces every number, figure, and table in the manuscript.
# Run from the repository root:  source("R/run_all.R")

set.seed(42)
dir.create("out", showWarnings = FALSE)

# 1. Main analysis -------------------------------------------------------
# Extract the R code from the Quarto document and run it in this session,
# so that the objects it creates (cwt_fit, pred_dat, res2, m, ...) are
# available to the downstream scripts.
knitr::purl("analysis/NOF_troll_analysis.qmd",
            output = "analysis/NOF_troll_analysis.R",
            documentation = 0, quiet = TRUE)
source("analysis/NOF_troll_analysis.R", echo = FALSE)

# 2. Uncertainty: posterior simulation and bootstrap (Figure 7, Figure S1)
n_draws    <- 1000
kappa      <- 50
gam_method <- "both"
source("R/method2_uncertainty.R")

# 3. Dirichlet prior sensitivity (Table S1) -----------------------------
gam_method <- "posterior"
source("R/method2_prior_sensitivity.R")

# 4. Hold-out evaluation (Figure 10) ------------------------------------
source("R/method_holdout_validation.R")

# 5. Record the computing environment -----------------------------------
writeLines(capture.output(sessionInfo()), "out/sessionInfo.txt")
message("Done. Outputs written to out/")
