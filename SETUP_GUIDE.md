# Setup guide (for Jake — delete this file before publishing)

Steps to get from your working folder to a public, citable repository. Do them in order.

---

## Step 1. Confirm data permissions (before anything is public)

- [ ] **CAMP-derived lookup.** `data/derived/indicator_stock_lookup.csv` holds tag code → indicator stock → brood year, extracted from CAMP. Confirm with the PSC that this derived table can be published. If they decline, leave the file out, keep the code, and state in the Data Availability Statement that stock assignments are available from the PSC on request.
- [ ] **RMIS extracts.** RMIS data are public. Record the retrieval date and query parameters in `data/README.md` (the `[MONTH YEAR]` and `[REPORT TYPE / QUERY PARAMETERS]` placeholders).
- [ ] **Treaty fishery data.** Tribal catch and sample records are public in RMIS, but given the co-management context it may be worth a courtesy heads-up to NWIFC before posting them.

---

## Step 2. Copy your files into the repository

| Your file | Goes to |
|---|---|
| Current `NOF_troll_quarto_v3.qmd` | `analysis/NOF_troll_analysis.qmd` |
| `method2_uncertainty.R` | `R/method2_uncertainty.R` |
| `method2_prior_sensitivity.R` | `R/method2_prior_sensitivity.R` |
| `method_holdout_combined.R` (the **corrected** version: Equation 2 ratio estimator, `cwt_recovs` truth, common-strata scoring, three-row figure) | `R/method_holdout_combined.R` |
| `Data/Troll_CatchSample_2005_2023.txt` | `data/raw/` |
| `Data/NOF Troll Recoveries 2005-2025.txt` | `data/raw/` |
| `Data/NOF Troll Tag Releases.txt` | `data/raw/` |

Do **not** copy `Data/camp_tables.rds`, `camp2025.config`, `New Tag Codes 2025.xlsx`, or any manuscript `.docx`. The `.gitignore` blocks these as a safety net, but don't rely on it.

---

## Step 3. Edit the code for the public repository

### 3a. Build the stock lookup (once, on your machine)

In your **private** working folder (where `Data/camp_tables.rds` exists), copy in `R/00_export_public_data.R`, change its input path to point at your recovery file, and run it. Copy the resulting `indicator_stock_lookup.csv` into the repository's `data/derived/`.

### 3b. Replace the CAMP chunk in the .qmd

Replace the whole `camp-connect` chunk (the database connection, `cds`, `era_tag_codes`, `era_recovery`) with:

````
```{r}
#| label: stock-lookup
#| echo: false

# PST indicator stock assignments by tag code, derived from the PSC CAMP
# database (see data/README.md). The full CAMP database is not public.
pad_tag <- function(x) str_pad(as.character(x), width = 6, pad = '0')

stock_lookup <- read.csv("data/derived/indicator_stock_lookup.csv",
                         colClasses = c(tag_code = "character")) |>
  mutate(tag_code = pad_tag(tag_code))

stopifnot(!anyDuplicated(stock_lookup$tag_code))
```
````

`era_recovery` was dead code (nothing downstream uses it), and `cds` was read but never used, so removing them changes no results.

### 3c. Point the stock-composition chunk at the lookup

In the `stock-composition` chunk, replace

```r
cwt_releases <- camp_tables$camp_cwt_release
```

with

```r
cwt_releases <- stock_lookup
```

and delete the `pad_tag <- …` definition in that chunk (it now lives in `stock-lookup`).

### 3d. File paths and housekeeping

- [ ] Find and replace `"Data/` → `"data/raw/` throughout the .qmd.
- [ ] Remove `library(odbc)` and `source("Code/quarto_install_scripts.R")`. `renv` handles packages.
- [ ] In `R/method2_prior_sensitivity.R`, change `source('method2_uncertainty.R')` → `source('R/method2_uncertainty.R')`.
- [ ] In `R/method_holdout_combined.R`, delete the unused `VS_CAP` line.
- [ ] Make sure every `ggsave()` writes to `out/` (some still point to `Auxiliaries/NOF Troll/out/`).
- [ ] Remove the stale qmd text flagged earlier: the Conclusions & Recommendations block (93-day window, τ, compounding) and the old Method 1 references. Reviewers may read the code.

### 3e. Test from a clean session

Restart R with an empty environment, set the working directory to the repository root, and run:

```r
source("R/run_all.R")
```

Check that the headline numbers match the manuscript (1,383; 906; 1,355; 1,254; the hold-out medians). If anything fails, it is almost always a path or an object that only existed in your old session.

---

## Step 4. Lock package versions with renv

From the repository root:

```r
install.packages("renv")
renv::init()       # scans the code, records packages and versions
renv::snapshot()   # writes renv.lock
```

Commit `renv.lock`, `.Rprofile`, and `renv/activate.R`. The `.gitignore` already excludes the package library itself. Anyone can then run `renv::restore()` to get your exact package versions.

---

## Step 5. Put it on GitHub (GitHub Desktop — no command line)

1. Create a free account at github.com if you don't have one. Consider whether WDFW has an organization account it prefers you use.
2. Install **GitHub Desktop** (desktop.github.com) and sign in.
3. *File → Add local repository →* choose the `nof-troll-cwt-unsampled` folder. When it says it isn't a Git repository, click **create a repository**. Leave "Initialize with README" unchecked, since you already have one.
4. Look over the list of files in the left panel. This is your last check that nothing private (`.rds`, `.config`, `.docx`) is about to go up.
5. Type a summary ("Initial release for NAJFM submission") and click **Commit to main**.
6. Click **Publish repository**. Choose **private** for now (you can make it public at acceptance, or at submission if the journal wants reviewer access).
7. After this, the routine is: edit files → GitHub Desktop shows the changes → write a one-line summary → **Commit** → **Push**.

---

## Step 6. Get a permanent DOI with Zenodo

GitHub is not an archive: repositories can be deleted or renamed. Journals prefer a DOI, and Zenodo (free, run by CERN) provides one by archiving a snapshot of your repository.

1. Sign in at zenodo.org using your GitHub account.
2. Under *GitHub* in your Zenodo account settings, switch on the toggle for `nof-troll-cwt-unsampled`. The repository must be public for this step.
3. On GitHub, go to the repository → **Releases → Create a new release**, tag it `v1.0.0`, and publish.
4. Zenodo archives that release and issues a DOI within a few minutes. Add the DOI to `CITATION.cff`, `README.md`, and the Data Availability Statement.

If you want the repository private until acceptance, do Steps 5 and 6 in two stages: private repository at submission (shared with the editor if requested), then public plus Zenodo release at acceptance, with the DOI added at proofs.

---

## Step 7. Data Availability Statement (for the manuscript)

> The code and data that support the findings of this study are openly available on GitHub at https://github.com/USERNAME/nof-troll-cwt-unsampled and archived at Zenodo (https://doi.org/10.5281/zenodo.XXXXXXX). Catch, sampling, coded-wire tag recovery, and release data were obtained from the Regional Mark Information System (RMIS; https://www.rmpc.org), maintained by the Pacific States Marine Fisheries Commission. Stock assignments for Pacific Salmon Treaty indicator stocks were derived from the Chinook Analysis Modeling Platform (CAMP), maintained by the Pacific Salmon Commission; a derived table of the tag code, stock, and brood year assignments used in this analysis is included in the archive. The full CAMP database is not publicly available, and requests for access should be directed to the Pacific Salmon Commission.

If the PSC declines to allow the lookup to be published, replace the second-to-last sentence with: "Stock assignments for the tag codes used in this analysis are available from the Pacific Salmon Commission upon request."
