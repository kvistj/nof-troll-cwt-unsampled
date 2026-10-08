---
title: "North of Falcon troll coded-wire tag auxiliary methods for the CTC ERA"
author: "Jake T. Kvistad, Ph.D."
bibliography: references.bibtex
reference-section-title: "References"
reference-location: document
csl: https://raw.githubusercontent.com/citation-style-language/styles/master/american-fisheries-society.csl

## Rendering format parameters
format: 
  html:
    embed-resources: true
  pdf:
    keep-md: true
    keep-tex: true
editor: source

## Execute parameters are used for default code chunk options ##
execute: 
  eval: true
  message: false
  echo: true
  warning: false
  error: false
---

## Background

The North of Falcon Troll fisheries were identified in CTC Report 46 [@PSC2021] as missing sampling in the catch and sample reports hosted on RMIS. Missing sampling in the catch and sample reports results in an undefined CWT expansion factor because they are calculated as $Number Caught \div Number Sampled$ and missing sampling in the denominator results in a divide by zero situation. Because there is no sampling from these reporting strata, this also means that there are no coded-wire tag returns during these missing sampling periods. This all has significance for PST management because annual calendar year exploitation rates (CYER) are evaluated against the CYER average in the 2009\textendash 2015 base period and the CTC calculates CYERs by summing fishery-specific AEQ mortalities and dividing by the sum of all AEQ mortalities plus escapement (@eq-cyer). Unaccounted mortalities in any individual fishery can result in inflated total CYERs because the number of mortalities in the denominator will be underestimated. This is particularly important in the context of PST management because as our sampling coverage expands, sampling methods and technologies improve, and our CYER estimates get more accurate, they are still being compared to a base period average that is probably too low considering the extent of missing sampling during the 2009\textendash 2015 base period. This increases the risk of erroneously exceeding CYER limits outlined in the PST for any Attachment I stock. The auxiliary method described herein fills gaps in the CWT expansion factors on RMIS and estimates the expected number of ERA tag recoveries from unsampled strata in an attempt to improve mortality accounting in SUS ISBM fisheries.

$$CYER_{CY,F}=\frac{\sideset{}{_\text{a=minage}^{\text{maxage}}}\sum_{}\sideset{}{_{f\in{F_{ISBM}}}}\sum_{}Morts_{CY,a,f}\cdot AEQ_{BY=CY-a,a,f}}{\sideset{}{_\text{a=minage}^{\text{maxage}}}\sum_{}\sideset{}{_\text{f=1}^{\text{NumFisheries}}}\sum_{}Morts_{CY,a,f}\cdot AEQ_{BY=CY-a,a,f}+Esc_{CY,a}}$$ {#eq-cyer}


::: {.cell}

:::



::: {.cell}

:::


## Data

Catch and sample reports were retrieved for all troll fisheries from RMIS for the years spanning 2005 through 2023 and filtered to the catch location names corresponding to the North of Falcon catch areas (WA Marine Area Codes 1\textendash 4 and 4B). Both treaty and non-treaty fisheries had notable amounts of sample data either missing or reported as 0 in the catch and sample reports (@fig-missing). Since the treaty fishery catch location names in RMIS did not have a consistent definition over time, most of the treaty catch location names were renamed to a generic "3,4 TT". The treaty fishery in MA 2 did not have any sampling reported to RMIS at all since the beginning of the time series; therefore, the MA 2 treaty fishery will require its own follow up analysis to properly address.


::: {.cell}

:::



::: {.cell}
::: {.cell-output-display}
![](NOF_troll_analysis_files/figure-pdf/littles-mcar-1.pdf)
:::
:::



::: {.cell}
::: {.cell-output-display}
![Sum of unsampled catch in the North of Falcon troll fisheries by month, gear type, catch location, and fishery (10 = Non-treaty, 15 = Treaty).](NOF_troll_analysis_files/figure-pdf/fig-missing-1.pdf){#fig-missing}
:::
:::


Out of 2569 records, 1013 (39%) had catch reported but no sampling, 16 (0.6%) had sampling reported but no catch reported, and 3 (0.1%) had the number of fish caught reported as 0 but had a non-zero amount sampled. The first case (catch reported but no sampling) was the most common and believed to reflect a systemic dynamic in the sampling regimes in these fisheries while the other two cases are probably clerical errors.


::: {.cell}

:::


## Monthly Aggregation and Donor Matching

The reporting strata on RMIS are statistical weeks, but not every week is sampled in the NOF troll fishery areas. Aggregating to a monthly resolution reduces the number of unsampled strata and improves the precision of expansion factors by increasing sample sizes, at the cost of smoothing within-month variability. Because the CTC runs its models on annual time steps, monthly resolution is adequate. Once expansion factors are computed monthly ($\varsigma = N^{Caught} / N^{Sampled}$), CWT recoveries from sampled weeks are expanded to the full monthly catch, so unsampled weeks within a sampled month require no separate estimate. The strata that remain unaccounted for are months in which no sampling occurred at all, and these are the target of both estimators below.

The ratio estimator (Method 1) requires a sampled donor for each unsampled month. Donors were selected by nearest-neighbor matching within the same fishery and area (@eq-nearest-neighbor; Appendix A):

$$s^* = \underset{s \in S_{f,a}}{\operatorname{argmin}} |t_u - t_s| \quad \text{subject to} \quad |t_u - t_s| \leq \delta$$ {#eq-nearest-neighbor}

where $s^*$ is the matched donor, $S_{f,a}$ is the set of sampled strata in fishery $f$ and area $a$, $t_u$ and $t_s$ are the month midpoints of the unsampled and sampled strata, and $\delta = 90$ days is the maximum search window.

## Estimating CWT Recoveries in Unsampled Strata

Once CWT expansions are recalculated by monthly strata, the expected number of ERA-tagged CWT recoveries must be imputed somehow such that the recalculated expansion factors can be applied and more accurate mortality accounting can be achieved in the ERA. The number of ERA-tagged CWT recoveries in the observed catch can be obtained from the Chinook Analysis Modeling Platform (CAMP).


::: {.cell}

:::


### Method 1: Ratio Estimator

The ratio estimator applies the stock-and-age composition of the nearest sampled stratum, including its untagged fraction, to the landed catch of each unsampled stratum. For an unsampled stratum $r$ and stock-age category $k$:

$$\widehat{CWT}_{r,k} = L_r \cdot \frac{n_{d,k}}{N^{Sampled}_d}$$ {#eq-ratio-scalar}

where $d$ is the nearest sampled donor stratum in the same fishery and area within $\delta = 90$ days (@eq-nearest-neighbor), $n_{d,k}$ is the number of CWTs of stock-age $k$ recovered in the donor, $N^{Sampled}_d$ is the number of fish sampled in the donor, and $L_r$ is the landed catch of the recipient. The ratio $n_{d,k}/N^{Sampled}_d$ is the donor's tag density for category $k$, so the estimate is on the same scale as the tag density model: the expected number of CWT-bearing fish in the landed catch. No separate expansion factor is applied. Donor recoveries are taken from the RMIS recovery data with CAMP stock labels, so non-PST stocks are included. The estimator is unbiased if tag density by stock and age is equal between donor and recipient strata.

The estimator is applied after the stock-composition chunks, where the recovery data it requires are assembled (see the `ratio-estimator` chunk).

### Method 2: Modeled CWT Encounter Rate

An alternative approach is to model the CWT encounter rate directly and apply it to the landed catch reported in unsampled strata, using time- and area-specific stock and age compositions to resolve the resulting estimates into stock- and age-specific components. Because the encounter rate is defined per unit of sampling effort and catch is reported in every stratum regardless of whether sampling occurred, the expected number of CWT recoveries in an unsampled stratum can be estimated without reconstructing the sampling process itself.

Stock and age compositions in the sampled fisheries can be derived empirically from sampling data or estimated via a model. These compositions, combined with hatchery tagging rates from release data, are then used to calculate expected CWT recoveries in unsampled fisheries. A key limitation of the empirical approach is that stock and age compositions observed in the sampled catch are assumed to hold in the unsampled catch as well. Modeling stock and age compositions in unsampled fisheries could relax this assumption, but at the cost of added analytical complexity.

##### CWT Encounter Rate

A core assumption of Method 1 (Ratio Estimator) was that the probability of CWT encounters is consistent across sampled and unsampled strata. This assumption is violated in practice due to a variety of factors, such as seasonal variations in run timing of the tagged PST indicator stocks and fluctuations in the abundance of unmarked fish present in the fishery. As a result, the average number of CWT encounters varies considerably within a calendar year (@fig-cwt-per-100).


::: {.cell}

:::



::: {.cell}

:::



::: {.cell}
::: {.cell-output-display}
![Median CWT recoveries per 100 landed catch in the North of Falcon Chinook troll fishery (2005--2023). One outlier removed.](NOF_troll_analysis_files/figure-pdf/fig-cwt-per-100-1.pdf){#fig-cwt-per-100}
:::
:::


The CWT encounter rate is the frequency at which CWTs are recovered from the sampled catch (@fig-cwt-encounter-rate). Let the encounter rate expressed per unit of landed catch be denoted by $\xi$, then:

$$\xi = \frac{N^{CWT}}{N^{Caught}}$$ {#eq-xi}

where $N^{CWT}$ represents the number of CWTs recovered and $N^{Caught}$ represents the total number of fish caught in a particular stratum. @eq-xi describes the frequency of CWT encounters relative to the landed catch and is useful for describing seasonal and spatial variation in encounter frequency; it is not, however, the quantity that generates the observed recoveries. A coded-wire tag can only be recovered from a fish that was sampled, so the count of recoveries in a stratum is generated by the number of fish sampled and the density of tagged fish among them. Let this tag density be denoted by $\theta$:

$$\theta = \frac{N^{CWT}}{N^{Sampled}}$$ {#eq-theta}

Throughout, a circumflex denotes a model-derived estimate: $\theta$ is the unobserved tag density in a stratum and $\hat{\theta}$ is the value predicted by the model described below. The expected number of CWT recoveries present in the landed catch of an unsampled stratum ($CWT_{Est}$) then follows from applying $\hat{\theta}$ to the reported catch:

$$CWT_{Est} = \hat{\theta}\cdot N^{Caught}$$ {#eq-cwt-est}

This substitution performs the expansion from sample to catch scale; no separate expansion factor is applied. Estimating $\theta$ requires only strata in which sampling occurred, where $N^{Sampled}$ is observed by definition, so no estimate of the missing sampling is required at any point. The approach does assume that tag density in unsampled strata is consistent with tag density in sampled strata of comparable year, month, and area.


::: {.cell}
::: {.cell-output-display}
![Median CWT encounter rates {@eq-xi} in the North of Falcon troll fishery (2005--2023), stratified by recovery location and month.](NOF_troll_analysis_files/figure-pdf/fig-cwt-encounter-rate-1.pdf){#fig-cwt-encounter-rate}
:::
:::


I used a generalized additive model (GAM) to model the number of CWT recoveries directly, with sampling effort entered as an offset (@eq-gam).

$$\begin{aligned}
N^{CWT}_i &\sim \text{Tweedie}(\mu_i,\phi,p) \\
\log(\mu_i) &= \log(N^{Sampled}_i) + \mathbf{X}_i\boldsymbol{\beta}+\sum_j f_j(x_{ij})
\end{aligned}$$ {#eq-gam}

The GAM included smooth terms $[f_j(\cdot)]$ on `run_year` and `month`, and a categorical variable modeling the parametric effect of `recovery_location_name`, where $\mathbf{X}_i$ is the $i$-th row of the parametric model matrix and $\boldsymbol{\beta}$ the vector of parametric coefficients. Because the treaty catch locations were consolidated to a single label, fishery is perfectly determined by recovery location and was not included as a separate term. Sampling effort entered as an offset with its coefficient fixed at unity, which constrains expected recoveries to scale proportionally with the number of fish sampled and makes the linear predictor an estimate of tag density; $\exp(\mathbf{X}_i\hat{\boldsymbol{\beta}}+\sum_j \hat{f}_j(x_{ij}))$ recovers $\hat{\theta}_i$ directly, and $\hat{\mu}_i$ is the expected number of recoveries. Specifying sampling effort as an offset rather than as a covariate avoids the circularity that arises when a variable appears simultaneously in the denominator of the response and among the predictors. Predictions for unsampled strata were generated by substituting $\log(N^{Caught})$ for the offset, which scales tag density to the full landed catch. I used a Tweedie distribution with a log-link because of the right-skewness in the recovery counts and the increase in variance with the mean (@fig-mean-variance). The Tweedie distribution includes a dispersion parameter ($\phi$) and a power parameter ($p$), both of which are estimated when fitting the model. The GAM was fit in R using the `mgcv` package (@Wood2011).


::: {.cell}

:::



::: {.cell}
::: {.cell-output-display}
![a) Sampling distribution of CWT recovery counts with probability density function overlain. b) Binned mean of the recovery counts plotted against the binned variance (log scale), with the fitted Tweedie variance function $\hat{\phi}\mu^{\hat{p}}$ overlain in blue.](NOF_troll_analysis_files/figure-pdf/fig-mean-variance-1.pdf){#fig-mean-variance}
:::
:::



::: {.cell}
::: {.cell-output-display}
![Ten-fold cross-validated predictions plotted against observed values for (a) CWT recovery counts and (b) tag density ($\theta$, recoveries per fish sampled). Axes are on a square-root scale; dashed red lines show a 1:1 relationship. Agreement on the count scale largely reflects sampling effort, which enters the model as an offset; panel b removes sampling effort and shows agreement in the quantity the model estimates. Point area in panel b is proportional to the number of fish sampled; the wider scatter among small-sample strata reflects sampling variability in observed tag density.](NOF_troll_analysis_files/figure-pdf/fig-k-fold-cv-1.pdf){#fig-k-fold-cv fig-pos='H'}
:::
:::


The model was fit to 467 sampled strata, including 53 in which no CWTs were recovered. Because sampling effort is accounted for by the offset, deviance explained is reported relative to an offset-only model: year, month, and recovery location explained 39.5% of the remaining deviance, and the fitted model was strongly favored over the offset-only alternative ($\Delta$AIC = 197.4). To assess out-of-sample prediction capability, I conducted a k-fold cross validation routine (k=10). Cross-validated predictions tracked observed counts closely ($R^2$ = 0.91; RMSE = 60.6) and reproduced aggregate recoveries to within 2.5% (@fig-k-fold-cv a). Because counts scale with sampling effort, agreement on the count scale largely reflects the offset. On the tag density scale, cross-validated predictions explained 42% of the variation in observed $\theta$ when strata were weighted by the number of fish sampled (18% unweighted; @fig-k-fold-cv b). Calibration was 1.08 in the 56 sampled winter strata and 1.02 in summer strata. Note that cross-validation evaluates prediction within sampled strata; the assumption that tag density transports to unsampled strata cannot be tested with these data.


::: {.cell}

:::


::: {.cell}

:::


##### Stock and Age Composition

Once $CWT_{Est}$ is evaluated, empirical stock and age composition data can be used to partition the estimated CWT recoveries from the unsampled strata into their stock- and age-specific components. Composition was applied using a lagging five-year window: recoveries from years $t-4$ through $t$ inform composition for year $t$. A lagging rather than centered window was used because composition in a given year reflects brood cohorts that have already recruited to the fishery; recoveries from subsequent years carry information about cohorts not yet present. The time-pooled composition for a given month and area was used as a fallback for year-month-area combinations lacking recovery data.


::: {.cell}

:::


Empirical estimates of stock and age composition can be calculated from RMIS recovery data (@fig-stock-comp). Stocks were assigned to RMIS data by joining the stock names in CAMP based on tag code. A considerable number of the recoveries in the RMIS data were unassigned, indicating considerable influence by non-PST stocks on the stock composition of the fishery. It is notable, however, that the proportion of the stock composition comprising Puget Sound stocks rises in the winter months (Nov-Apr).


::: {.cell}

:::



::: {.cell}
::: {.cell-output-display}
![Empirical monthly stock (a) and age (b) composition estimated from coded-wire tag returns in the North of Falcon troll Chinook fishery (2005--2023).](NOF_troll_analysis_files/figure-pdf/fig-stock-comp-1.pdf){#fig-stock-comp fig-pos='H'}
:::
:::


The monthly age composition estimated from CWTs shows the fishery largely dominated by Age 3 and 4 fish, with a small increase in the number of Age 5+ fish over the summer months followed by the arrival of Age 2 and more Age 3 fish in November and December. There are currently only two tag recoveries in this fishery from the month of October and both were Age 4s.


::: {.cell}
::: {.cell-output-display}
![Predicted CWT recoveries from unsampled strata in the North of Falcon Chinook troll fishery (2005--2023) apportioned into stock and age specific components.](NOF_troll_analysis_files/figure-pdf/fig-cwt-stock-age-1.pdf){#fig-cwt-stock-age fig-pos='H'}
:::
:::



::: {.cell}

:::



::: {.cell}

:::



::: {.cell}
::: {.cell-output-display}
![Estimated CWT recoveries in unsampled strata of the North of Falcon Chinook troll fishery (2005--2023) from the nearest-donor ratio estimator, by stock region and age.](NOF_troll_analysis_files/figure-pdf/fig-res-1.pdf){#fig-res fig-pos='H'}
:::
:::


The ratio estimator found a donor within 90 days for 90 of the 105 unsampled strata; the remaining 15 strata (1,181 fish) had no sampled stratum in the same fishery and area within the search window. Across matched strata, the ratio estimator estimated 1,355 CWT recoveries, of which 632 (47%) were PST indicator stocks, compared with 1,235 from the tag density model for the same strata. Donors sampling less than 5% of their catch (4 strata) contributed 159 of the total. Of 41 winter recipients, 23 borrowed from a summer donor.

Applying the stock and age compositions to the predicted CWT recoveries yielded the results displayed in @fig-cwt-stock-age. In sum, 1,362 additional CWT recoveries were expected to have occurred between 2005 and 2023 in the unsampled strata. Of those, 1,351 could be apportioned into stock and age components; the remainder fall in 8 strata for which no CWT recoveries have been reported in the corresponding month and area over the full time series. Of the apportioned total, 839 were expected to be PST indicator stocks. To quantify the statistical uncertainty in the estimates, I performed a bootstrap simulation that resampled the model training data with replacement and drew stock and age compositions from their Dirichlet posteriors. Details are provided in [Appendix B](#sec-appB).

Auxiliaries for the CTC ERA can be prepared by mapping the additional CWT recoveries to the appropriate tag codes. Accounting for the additional fishery mortalities in the unsampled strata should reduce bias in the CYERs and more accurately measure fishery exploitation on PST indicator stocks.

## Conclusions & Recommendations

The analyses contained within this brief offer a comparison of two separate methods of estimating CWT recoveries from unsampled fishery strata. The first method is a purely analytical one in which CWT recoveries are calculated with a ratio estimator assuming constant stock and age composition between sampled and unsampled fishery strata. The second method models tag density directly as a function of year, month, and recovery area, scales it to the landed catch reported in each unsampled stratum, and partitions the resulting estimates using empirical stock and age composition data.

Methods differed substantially in the number of PST indicator stock CWTs estimated for the unsampled strata (Method 1: 632 compared to Method 2: 839). These differences may be related to the ways in which CWT recoveries are estimated in the unsampled strata. Method 1 assumes that CWT recoveries in the unsampled strata are directly proportional to the landed catch ratio of the sampled and unsampled strata (@eq-ratio-scalar); violations of this assumption, such as differences in the relative abundance of ERA-tagged fish between sampled and unsampled strata, could result in inaccurate estimates. Method 1 also requires that a temporally adjacent sampled stratum exist within the search window, and inherits whatever composition that neighbor happens to carry. By contrast, Method 2 models tag density as a smooth function of year, month, and recovery area, which relaxes the assumption of a constant encounter frequency across strata. Tag density varies systematically across months and catch areas (@fig-cwt-per-100, @fig-cwt-encounter-rate), and year, month, and area together explained 39.5% of the deviance remaining after sampling effort is accounted for ($\Delta$AIC = 197.4 relative to an offset-only model). Modeling that variation directly is therefore preferable to borrowing wholesale from an adjacent sampled stratum.

Stock composition in Method 1 was restricted to Puget Sound and WA Coast stocks (@fig-res) because they were wholly informed by the stock and age compositions of nearby sampled strata (i.e., \leq93 days away from the unsampled catch event within the same recovery area) which comprised only those stocks, whereas stock composition in Method 2 was informed by the average empirical distribution of tag recoveries that have occurred in those times and areas over the entire timescale of the dataset (2005--2023). A similar pattern was observed for age composition, although age composition differed less between methods. Differences between both methods might be explained by how each handles stock and age composition. Method 1 borrows the stock and age composition of a single nearby sampled stratum; Method 2 instead uses the composition of all CWT recoveries reported in the same month and area over a lagging five-year window, and treats that composition as uncertain by drawing it from a Dirichlet posterior rather than fixing it at the observed proportion. This matters most where recoveries are sparse: a stratum informed by one or two tags yields an observed composition that is necessarily degenerate, assigning exactly zero to every stock that happened not to appear. Additionally, the empirical stock and age composition data implies that the stock and age composition in unsampled strata (predominantly occurring over winter months) is not always consistent with the stock and age compositions in the sampled strata (predominantly occurring over the summer/fall months) (@fig-stock-comp). For example, the winter months see an increase in representation by Puget Sound stocks and trends towards younger ages in the CWT recoveries. Several limitations of Method 2 warrant note. Predictions for unsampled strata rest on the assumption that tag density is consistent with that in sampled strata of comparable year, month, and area; cross-validation evaluates prediction only within sampled strata and cannot test this assumption. Sampling effort is treated as measured without error, though a small number of strata reported recoveries in excess of the fish sampled, indicating that reported sampling is occasionally unreliable. Compositional data are thin in precisely the strata of greatest interest: 23 predicted strata, carrying 21% of the unsampled catch, have fewer than ten CWT recoveries informing their composition, and 8 strata had none at all. Finally, the treaty fishery in MA 2 is excluded from the analysis entirely and requires separate treatment.

I recommend using Method 2 for the following reasons:

1)  Modeling tag density directly (i.e., the GAM) produced substantially different estimates of total CWT recoveries of PST indicator stocks than Method 1. I hypothesize that the 93 day search window used by the nearest-neighbor approach in Method 1 to infill missing sampling data is too large, resulting in inappropriate CWT expansion factors; the landed catch ratios ($\tau$) used to scale CWT recovery estimates in the unsampled strata are probably too large for the same reasons. By contrast, Method 2 estimates recoveries from a tag density defined per unit of sampling effort and applies it to the catch actually reported in each unsampled stratum. It therefore does not depend on the quality of a temporal match to a neighboring sampled stratum, and it requires no estimate of how many fish would have been sampled had sampling occurred.

2)  Method 2 relaxes the equal stock and age composition assumption of Method 1 by using month- and area-specific empirical compositions, and propagates the uncertainty in those compositions into the reported confidence intervals rather than treating them as known.

## References

::: {#refs}
:::

\newpage

## Appendix A {#sec-appA}

::: {.callout-note appearance="minimal"}
**Method 1 Pseudocode**

For each unsampled record $u$:

1.  Identify all sampled records $S$ sharing the same fishery and catch area
2.  Compute temporal distance $d = |t_u - t_s|$ for each $s \in S$
3.  Select $s^* = \underset{s \in S}{\operatorname{argmin}}(d)$, subject to $d \leq \delta$
4.  If no $s^*$ satisfies $d \leq \delta$, leave $u$ unestimated
5.  For each stock-age category $k$, estimate $\widehat{CWT}_{u,k} = L_u \cdot n_{s^*,k} / N^{Sampled}_{s^*}$
6.  Ties in $d$ are resolved in favor of the earlier donor; $\delta = 90$ days between month midpoints
:::

\newpage

## Appendix B {#sec-appB}

**Method 2 Bootstrap Uncertainty Analysis**

I quantified uncertainty using a nonparametric bootstrap that propagated two independent sources of error to the final CWT estimates: error in the fitted tag density model, and error in the empirical stock and age composition used to apportion the resulting counts. I ran 1,000 bootstrap iterations. In each iteration the GAM training data were resampled with replacement and the model refit, predictions were generated for the unsampled strata by substituting landed catch for the sampling offset, and a stock and age composition vector was drawn for each stratum from its Dirichlet posterior.

Composition uncertainty was treated as a multinomial sampling problem. For a stratum in which the windowed recoveries yield counts $n = (n_1,\dots,n_K)$ across $K$ stock-age categories, and with a Dirichlet prior $\text{Dir}(\alpha)$ on the composition $p$, the posterior is conjugate:

$$p \mid n \sim \text{Dir}(\alpha + n)$$

Counts were assembled over a fixed category grid spanning all observed stock-age combinations, so that categories unobserved in a given stratum retain positive posterior probability rather than being treated as structurally absent. I compared three priors: Jeffreys ($\alpha_k = 1/2$), uniform ($\alpha_k = 1$), and an empirical-Bayes prior centered on the time-pooled composition for that month and area ($\alpha_k = \kappa\bar{p}_k$). The first two center the prior on equal proportions across all $K$ categories, which is inconsistent with the observed composition of the fishery; the empirical-Bayes prior instead shrinks sparse strata toward the pooled composition, formalizing the borrowing that the rolling window performs implicitly. Results are reported under the empirical-Bayes prior. The 95% confidence intervals were calculated from the 2.5th and 97.5th percentiles of the bootstrap distribution and the medians (50th percentile) are reported as the central tendencies.

The bootstrap median for the total estimated CWT recoveries across all stocks (including non-PST stocks) was 1,772 (95% CI: 1,532--2,008), closely matching the point estimate of 1,362. For PST indicator stocks only, the median estimate under the empirical-Bayes prior was 1,153 (95% CI: 976--1,325).

The PST estimate was insensitive to the strength of the empirical-Bayes prior, varying by less than 2% across $\kappa$ ranging from 5 to 200 (Table B1). Prior *strength* therefore has little influence on the result. Prior *centering*, by contrast, has a substantial one: priors centered on equal proportions across all $K$ categories produced PST totals roughly 35% higher than the empirical-Bayes prior. Because a single non-PST category accounts for a large share of observed recoveries, distributing prior mass evenly across categories mechanically inflates the estimated PST share. That the three priors do not converge indicates that the sparse strata carry limited independent compositional information -- the same limitation quantified in the main text, where 23 of the predicted strata, carrying 21% of the unsampled catch, have fewer than ten CWT recoveries informing their composition.

| Prior | $\kappa$ | Median | 95% CI |
|:------|------:|-------:|:-----------|
| Empirical Bayes | 5   | 1,136 | 931--1,348 |
| Empirical Bayes | 25  | 1,147 | 974--1,330 |
| Empirical Bayes | 50  | 1,153 | 976--1,325 |
| Empirical Bayes | 200 | 1,149 | 993--1,322 |
| Jeffreys ($\alpha_k = 1/2$) | -- | 1,508 | 1,293--1,727 |
| Uniform ($\alpha_k = 1$)    | -- | 1,568 | 1,344--1,795 |

: Sensitivity of the estimated PST indicator stock total to the choice of Dirichlet prior. {#tbl-prior-sensitivity}

The upper bound of the 95% credible interval for PST stocks remained well below the Method 1 estimate of 632, indicating that the discrepancy between methods is not attributable to statistical uncertainty in Method 2. A small number of bootstrap iterations (fewer than 1 in 500) failed to converge when a resample contained too few unique values of `run_year` or `month` to support the default basis dimension; these iterations were discarded.

\newpage

```{=latex}
\begin{center}
\includegraphics[width=0.8\textwidth]{out/method2_prior_comparison.png}

\textbf{Figure B1.} Bootstrap distributions of total estimated CWT recoveries from unsampled strata under the empirical-Bayes prior. Dashed blue lines indicate median estimates and dashed red lines indicate 95\% confidence intervals.
\end{center}
```
