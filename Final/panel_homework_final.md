---
title: "Panel Data — Homework 2025–2026"
subtitle: "Pooled OLS, Fixed Effects, and Random Effects"
author: "Dragos Florin Vasile & Wong Hei Wong"
date: "March 2026"
fontsize: 11pt
geometry: margin=1in
header-includes:
  - \usepackage{booktabs}
  - \usepackage{amsmath}
  - \usepackage{amssymb}
  - \usepackage{float}
  - \usepackage{xcolor}
---

# Exercise 1. Theory (20 points)

We consider

$$
y_{it} = \alpha_i + \beta x_{it} + u_{it},
\qquad
E(u_{it}\mid X)=0,
\qquad
\text{Var}(u_{it}\mid X)=\sigma^2,
\tag{1}
$$

where the panel may be unbalanced: individual $i$ is observed only in periods $t\in P_i$, and $T_i=|P_i|$.

Averaging equation (1) over the observed periods of individual $i$ gives the **between equation**:

$$
\bar y_{i.} = \alpha_i + \beta\,\bar x_{i.} + \bar u_{i.}\,,
\qquad
\bar y_{i.} = \frac{1}{T_i}\sum_{t\in P_i}y_{it}\,,
\qquad
\bar x_{i.} = \frac{1}{T_i}\sum_{t\in P_i}x_{it}\,.
\tag{2}
$$

## 1. Derivation of the within estimator (3 pts)

We treat the $\alpha_i$'s as unknown constants and solve the least-squares problem

$$
\min_{\beta,\,\alpha_1,\ldots,\alpha_N}
\sum_{i=1}^{N}\sum_{t\in P_i}
\left(y_{it}-\alpha_i-\beta x_{it}\right)^2.
$$

**FOC w.r.t. $\alpha_i$.** Setting $\partial S/\partial\alpha_i=0$:

$$
\sum_{t\in P_i}\left(y_{it}-\alpha_i-\beta x_{it}\right)=0,
$$

so, dividing by $T_i$:

$$
\hat\alpha_i = \bar y_{i.} - \beta\,\bar x_{i.}\,.
$$

Substituting back eliminates the individual effects. The concentrated (demeaned) problem becomes

$$
\min_\beta \sum_{i=1}^{N}\sum_{t\in P_i}
\left[(y_{it}-\bar y_{i.}) - \beta\,(x_{it}-\bar x_{i.})\right]^2.
$$

**FOC w.r.t. $\beta$.** Setting the derivative to zero:

$$
\sum_{i=1}^{N}\sum_{t\in P_i}
(x_{it}-\bar x_{i.})\bigl[(y_{it}-\bar y_{i.})-\beta\,(x_{it}-\bar x_{i.})\bigr]=0,
$$

which gives

$$
\hat\beta_W
=
\frac{W_{XY}}{W_{XX}},
\qquad
W_{ZQ}
=
\sum_{i=1}^{N}\sum_{t\in P_i}
(z_{it}-\bar z_{i.})(q_{it}-\bar q_{i.}).
$$

Therefore

$$
\hat\beta_W = W_{XX}^{-1}\,W_{XY}\,,
\qquad
\hat\alpha_i^W = \bar y_{i.} - \hat\beta_W\,\bar x_{i.}\,.
\;\;\blacksquare
$$

## 2. Impact of unbalancedness (1 pt)

The form of the within estimator does not change. The only difference is that the individual means and within sums are computed over the observed periods $P_i$, so each individual contributes according to the observations actually available. Unbalancedness changes the sample moments, not the logic of the estimator.

We do not need successive observations for each individual either. The within transformation only requires that an individual is observed at least twice ($T_i\ge 2$) **and** that $x_{it}$ varies within that individual, so that the individual contributes non-zero within variation to $W_{XX}$. Gaps in time are not a problem; the estimator uses deviations from individual means, which do not depend on continuity of the time series.

## 3. Unbiasedness and variance of $\hat\beta_W$ (4 pts)

**Unbiasedness.** Substituting eq. (1) into $W_{XY}$:

$$
W_{XY}
=
\sum_{i=1}^{N}\sum_{t\in P_i}
(x_{it}-\bar x_{i.})(y_{it}-\bar y_{i.})
=
\beta\,W_{XX} + W_{XU},
$$

so

$$
\hat\beta_W = \beta + W_{XX}^{-1}\,W_{XU}.
$$

Now $W_{XU}=\sum_i\sum_{t\in P_i}(x_{it}-\bar x_{i.})(u_{it}-\bar u_{i.})$. Since $\sum_{t\in P_i}(x_{it}-\bar x_{i.})=0$, the $\bar u_{i.}$ term drops out, giving

$$
W_{XU}
=
\sum_{i=1}^{N}\sum_{t\in P_i}
(x_{it}-\bar x_{i.})\,u_{it}\,.
$$

Conditional on $X$, strict exogeneity implies

$$
E(W_{XU}\mid X) = 0
\qquad\Longrightarrow\qquad
E(\hat\beta_W\mid X) = \beta.
$$

**Variance.** The $u_{it}$ are conditionally homoskedastic with variance $\sigma^2$ and mutually uncorrelated across all $(i,t)$ pairs (spherical errors). Therefore all cross-terms in the variance vanish and only the diagonal remains:

$$
\text{Var}(W_{XU}\mid X)
=
\sigma^2 \sum_{i=1}^{N}\sum_{t\in P_i}(x_{it}-\bar x_{i.})^2
=
\sigma^2\,W_{XX}\,.
$$

Hence

$$
\text{Var}(\hat\beta_W\mid X)
=
\frac{\text{Var}(W_{XU}\mid X)}{W_{XX}^2}
=
\sigma^2\,W_{XX}^{-1}.
\;\;\blacksquare
$$

## 4. Is $\bar u_{i.}$ homoskedastic? (2 pts)

From eq. (2), $\bar u_{i.}=T_i^{-1}\sum_{t\in P_i}u_{it}$, so $\text{Var}(\bar u_{i.}\mid X)=\sigma^2/T_i$. When $T_i$ varies across individuals the between-equation error is **heteroskedastic**; it is homoskedastic only in the balanced case.

## 5. Transformation to restore homoskedasticity (1 pt)

Multiplying eq. (2) by $\sqrt{T_i}$:

$$
\sqrt{T_i}\,\bar y_{i.}
=
\sqrt{T_i}\,\alpha_i
+ \beta\,\sqrt{T_i}\,\bar x_{i.}
+ \sqrt{T_i}\,\bar u_{i.}\,.
\tag{3}
$$

The transformed error has variance $T_i\cdot\sigma^2/T_i=\sigma^2$, which is constant. This is a standard WLS correction weighting each individual by $\sqrt{T_i}$.

## 6. First between estimator $\hat\beta_{B1}$ (4 pts)

Setting $\alpha_i=\alpha$ for all $i$ in eq. (2) and applying OLS across the $N$ individual means:

$$
\hat\beta_{B1}
=
\frac{\sum_{i=1}^{N}(\bar x_{i.}-\bar x_B)(\bar y_{i.}-\bar y_B)}
{\sum_{i=1}^{N}(\bar x_{i.}-\bar x_B)^2}\,,
\qquad
\bar x_B = \frac{1}{N}\sum_{i=1}^{N}\bar x_{i.}\,,
\quad
\bar y_B = \frac{1}{N}\sum_{i=1}^{N}\bar y_{i.}\,.
$$

The intercept is $\hat\alpha_{B1} = \bar y_B - \hat\beta_{B1}\,\bar x_B$. Because this OLS ignores the heteroskedasticity identified in Q4, the estimator is unbiased but inefficient when $T_i$ varies.

## 7. Second between estimator $\hat\beta_{B2}$ (4 pts)

Applying OLS to the $\sqrt{T_i}$-transformed model (eq. 3) is algebraically equivalent to WLS with weights $T_i$:

$$
\hat\beta_{B2}
=
\frac{\sum_{i=1}^{N}T_i(\bar x_{i.}-\bar x_W)(\bar y_{i.}-\bar y_W)}
{\sum_{i=1}^{N}T_i(\bar x_{i.}-\bar x_W)^2}\,,
\qquad
\bar x_W
=
\frac{\sum_i T_i\,\bar x_{i.}}{\sum_i T_i}\,,
\quad
\bar y_W
=
\frac{\sum_i T_i\,\bar y_{i.}}{\sum_i T_i}\,.
$$

The intercept is $\hat\alpha_{B2} = \bar y_W - \hat\beta_{B2}\,\bar x_W$. Because $\sum_i T_i\bar x_{i.}$ is just the sum over all observed observations, $\bar x_W$ and $\bar y_W$ coincide with the global means of the raw data.

## 8. Lesson from comparing $B1$ and $B2$ (1 pt)

In a balanced panel ($T_i=T$ for all $i$), the weights are all equal and cancel from the ratio, so $\hat\beta_{B1}=\hat\beta_{B2}$ and $\hat\alpha_{B1}=\hat\alpha_{B2}$. In an unbalanced panel, the two diverge: $\hat\beta_{B1}$ gives every individual the same weight regardless of how many periods it was observed, while $\hat\beta_{B2}$ up-weights individuals with more observations, whose averages are more precise. Since the error variance is $\sigma^2/T_i$, $\hat\beta_{B2}$ is the GLS-efficient (BLUE) estimator in the between equation. The practical takeaway is that when the panel is unbalanced, ignoring the heteroskedasticity in the between model leads to a loss of efficiency.

\newpage

# Exercise 2. Application (20 points)

All regressions below use cluster-robust standard errors at the state level, which allows for arbitrary heteroskedasticity and serial correlation within states. This is the standard choice for panel data of this kind (see Chapter 1).

## 1. Data description and cleaning (2 pts)

The raw file contains $51\times 15=765$ state-year observations: the 50 U.S. states plus D.C., for 1983–1997. The preprocessing is:

```stata
import delimited "../SeatBelts.csv", delimiter(";") clear
isid state year
xtset fips year
gen ln_income = ln(income)
gen sample_main = !missing(fatalityrate, sb_useage, speed65, speed70, ///
    ba08, drinkage21, ln_income, age)
```

The check `isid state year` confirms there are no duplicates; `xtset fips year` declares the panel. The only relevant data issue is that `sb_useage` is missing for 209 observations, concentrated in the early years of the sample (roughly 1983–1989), when many states had not yet begun systematically tracking belt usage. All other variables are complete. Because `sb_useage` is a key regressor, we restrict the estimation sample to the 556 observations with non-missing values. This makes the sample unbalanced, but as shown in Exercise 1 this does not create any bias. We also generate `ln_income = ln(income)`, a standard log transformation for a variable measured in dollars.

| Statistic             |      Value |
| --------------------- | ---------: |
| Raw observations      |        765 |
| States (incl. DC)     |         51 |
| Years                 | 1983–1997 |
| Missing `sb_useage` |        209 |
| Estimation sample     |        556 |
| Mean `fatalityrate` |     0.0198 |
| Mean `sb_useage`    |      0.529 |
| Mean `income`       |   \$19,572 |

## 2. Pooled OLS (2 pts)

```stata
reg fatalityrate sb_useage speed65 speed70 ba08 drinkage21 ///
    ln_income age if sample_main, vce(cluster fips)
```

| Variable   |   Coeff. | Cl. SE |      t |     p |
| ---------- | -------: | -----: | -----: | ----: |
| sb_useage  |   0.0041 | 0.0024 |   1.70 | 0.095 |
| speed65    |   0.0001 | 0.0006 |   0.25 | 0.807 |
| speed70    |   0.0024 | 0.0007 |   3.46 | 0.001 |
| ba08       | −0.0019 | 0.0007 | −2.60 | 0.012 |
| drinkage21 |   0.0001 | 0.0013 |   0.06 | 0.950 |
| ln_income  | −0.0181 | 0.0024 | −7.47 | 0.000 |
| age        | −0.0000 | 0.0004 | −0.02 | 0.987 |
| Constant   |   0.1965 | 0.0198 |   9.92 | 0.000 |

$N=556$, $R^2=0.549$.

The pooled OLS coefficient on `sb_useage` is positive (+0.0041) and only marginally significant. Taken at face value this would suggest higher belt use is associated with *more* fatalities, which is not credible. The likely explanation is omitted variable bias. States with historically high fatality rates (due to rural roads, driving culture, etc.) may have adopted seat belt laws precisely *because* of that problem. Pooled OLS mixes cross-state differences with over-time variation and therefore picks up this selection effect.

## 3. State fixed effects (2 pts)

```stata
xtreg fatalityrate sb_useage speed65 speed70 ba08 drinkage21 ///
    ln_income age if sample_main, fe vce(cluster fips)
```

| Variable   |   Coeff. | Cl. SE |      t |     p |
| ---------- | -------: | -----: | -----: | ----: |
| sb_useage  | −0.0058 | 0.0017 | −3.46 | 0.001 |
| speed65    | −0.0004 | 0.0005 | −0.93 | 0.355 |
| speed70    |   0.0012 | 0.0003 |   3.54 | 0.001 |
| ba08       | −0.0014 | 0.0004 | −3.67 | 0.001 |
| drinkage21 |   0.0007 | 0.0007 |   1.04 | 0.305 |
| ln_income  | −0.0135 | 0.0024 | −5.67 | 0.000 |
| age        |   0.0010 | 0.0007 |   1.31 | 0.196 |

$R^2_{\text{within}}=0.687$, $\hat\sigma_u=0.00383$, $\hat\sigma_e=0.00179$, $\hat\rho=0.821$.

Yes, the results change sharply. The sign on `sb_useage` reverses to −0.0058 and is now significant at the 0.1% level. Once we control for permanent state characteristics, identification comes from within-state changes over time: within a given state, years with higher seat belt use see fewer fatalities. The sign reversal is a textbook illustration of omitted variable bias in the pooled model. The high $\hat\rho=0.82$ confirms that state-level unobservables account for the bulk of residual variance, which is why controlling for them matters so much.

## 4. State + time fixed effects (2 pts)

```stata
xtreg fatalityrate sb_useage speed65 speed70 ba08 drinkage21 ///
    ln_income age i.year if sample_main, fe vce(cluster fips)
```

| Variable   |   Coeff. | Cl. SE |      t |     p |
| ---------- | -------: | -----: | -----: | ----: |
| sb_useage  | −0.0037 | 0.0015 | −2.56 | 0.013 |
| speed65    | −0.0008 | 0.0006 | −1.35 | 0.183 |
| speed70    |   0.0008 | 0.0005 |   1.76 | 0.085 |
| ba08       | −0.0008 | 0.0004 | −1.86 | 0.069 |
| drinkage21 | −0.0011 | 0.0006 | −1.82 | 0.074 |
| ln_income  |   0.0063 | 0.0067 |   0.94 | 0.354 |
| age        |   0.0013 | 0.0007 |   1.90 | 0.063 |

Joint test on year dummies: $F(14,50)=9.82$, $p < 0.001$. $R^2_{\text{within}}=0.751$.

Adding year FE absorbs national trends common to all states — safer vehicles year after year, nationwide campaigns, improvements in trauma care, etc. The year dummies are jointly very significant, confirming that common time trends matter. The coefficient on `sb_useage` shrinks somewhat to −0.0037, which makes sense: part of what looked like a seatbelt effect in Q3 was actually a national downward trend in fatalities. After removing it, the remaining within-state, within-year variation still identifies a negative and significant effect. Note also that `ln_income` becomes insignificant here because its over-time trend is largely absorbed by the year dummies.

## 5. Most reliable specification (2 pts)

|                      | Pooled OLS |    State FE    |   State + Year FE   |
| -------------------- | :--------: | :------------: | :------------------: |
| `sb_useage` coeff. |  +0.0041  |    −0.0058    |       −0.0037       |
| $R^2$              |   0.549   | 0.687 (within) |    0.751 (within)    |
| Year FE              |     No     |       No       | Yes ($p < 0.001$) |

The **two-way fixed effects** (state + year) specification is the most reliable:

1. Pooled OLS is clearly contaminated by omitted state-level heterogeneity — the positive sign on `sb_useage` is a red flag.
2. State FE fix the cross-sectional confounding, but they do not control for aggregate time trends shared across all states.
3. Two-way FE remove both state heterogeneity and common shocks. The strongly significant year dummies show that at least some time effects are empirically important; ignoring them may contaminate the state-FE estimate. The two-way specification leaves identification to variation that differs *across states and over time*, which is the most credible source of variation in this context.

## 6. Effect size and lives saved (2 pts)

From Q4, $\hat\beta_{sb}=-0.0037$. This means a one-unit increase in `sb_useage` (from 0 to 1) lowers the fatality rate by 0.0037 deaths per million traffic miles. In percentage terms, a 1-percentage-point increase in belt use reduces fatalities by about

$$
\frac{0.0037\times 0.01}{0.0198} \approx 0.19\%
$$

of the mean fatality rate. That seems small per unit, but traffic volumes are enormous, so even small rate changes add up.

**Scenario: belt use from 52% to 90%.** The change is $\Delta = 0.90 - 0.52 = 0.38$. The predicted reduction in the fatality rate is

$$
\Delta\,\text{fatalityrate}
= -0.0037 \times 0.38
= -0.00141 \text{ per million miles}.
$$

To translate that into lives, we multiply by aggregate U.S. traffic volume:

$$
0.00141 \times 2{,}113{,}834 \approx 2{,}987 \text{ lives saved per year (average VMT)},
$$

$$
0.00141 \times 2{,}560{,}372 \approx 3{,}618 \text{ lives saved per year (1997 VMT)}.
$$

For context, roughly 42,000 Americans died in traffic crashes per year during this period, so the estimate implies a **7–9% reduction in total road fatalities** — a large and policy-relevant effect.

## 7. Random effects (2 pts)

```stata
xtreg fatalityrate sb_useage speed65 speed70 ba08 drinkage21 ///
    ln_income age if sample_main, re vce(cluster fips)
```

| Variable   |   Coeff. | Cl. SE |      z |     p |
| ---------- | -------: | -----: | -----: | ----: |
| sb_useage  | −0.0045 | 0.0016 | −2.73 | 0.006 |
| speed65    | −0.0003 | 0.0004 | −0.76 | 0.447 |
| speed70    |   0.0013 | 0.0004 |   3.80 | 0.000 |
| ba08       | −0.0014 | 0.0004 | −3.66 | 0.000 |
| drinkage21 |   0.0008 | 0.0007 |   1.09 | 0.278 |
| ln_income  | −0.0126 | 0.0020 | −6.36 | 0.000 |
| age        |   0.0002 | 0.0005 |   0.45 | 0.653 |

$\hat\sigma_u=0.00302$, $\hat\sigma_e=0.00179$, $\hat\rho=0.740$.

This RE specification corresponds to Q3 (state effects, no year dummies). The critical assumption is that the state effect $\alpha_i$ is uncorrelated with the regressors:

$$
E(\alpha_i \mid x_{i1},\dots,x_{iT}) = 0.
$$

If this holds, RE is more efficient because it uses both within-state and between-state variation. The coefficient on `sb_useage` is −0.0045, sitting between the FE estimate (−0.0058) and the pooled OLS (+0.0041). This is expected: RE is a matrix-weighted average of within and between estimators, so it partially reintroduces the cross-state variation that we know is problematic here. Whether the orthogonality assumption actually holds is what we test next.

## 8. Robust Hausman test (2 pts)

The classical Hausman statistic requires homoskedastic, serially uncorrelated errors to form a valid difference of covariance matrices. Since we rely on cluster-robust standard errors throughout, we instead implement the **Mundlak (1978) device**, which is the standard robust alternative covered in Chapter 2.

**Procedure.** We augment the RE model (matching Q7, without year dummies) with the within-state means of all time-varying regressors $\bar x_{i.}$ and test their joint significance:

$$
y_{it} = x_{it}'\beta + \bar x_{i.}'\gamma + \alpha_i^* + \varepsilon_{it}.
$$

Under $H_0\!:\gamma=0$, the individual means add no explanatory power and the RE assumption holds.

```stata
foreach var in sb_useage speed65 speed70 ba08 drinkage21 ln_income age {
    bysort fips: egen mean_`var' = mean(`var')
}
xtreg fatalityrate sb_useage speed65 speed70 ba08 drinkage21 ///
    ln_income age mean_sb_useage mean_speed65 mean_speed70 mean_ba08 ///
    mean_drinkage21 mean_ln_income mean_age, re vce(cluster fips)
test mean_sb_useage mean_speed65 mean_speed70 mean_ba08 ///
    mean_drinkage21 mean_ln_income mean_age
```

**Result:** $\chi^2(7)=28.42$, $p=0.0002$. We reject $H_0$ decisively.

The state effects *are* correlated with the regressors. This means the RE estimator is inconsistent, and **fixed effects is the appropriate strategy**. This is intuitive: states' unobserved characteristics (terrain, urbanisation, driving culture) are very likely correlated with policy variables such as seat belt enforcement and income levels. The Mundlak test confirms what economic reasoning already suggested.

## 9. Enforcement and seat belt usage (2 pts)

```stata
xtreg sb_useage primary secondary speed65 speed70 ba08 drinkage21 ///
    ln_income age i.year if sample_main, fe vce(cluster fips)
lincom primary - secondary
```

| Variable  | Coeff. | Cl. SE |    t |      p |
| --------- | -----: | -----: | ---: | -----: |
| primary   | 0.2056 | 0.0232 | 8.87 | <0.001 |
| secondary | 0.1085 | 0.0134 | 8.09 | <0.001 |

$R^2_{\text{within}}=0.842$. Test $H_0$: primary = secondary: $F(1,50)=17.6$, $p=0.0001$.

Both types of enforcement raise belt use, and both are highly significant. Primary enforcement increases usage by about **20.6 percentage points** relative to no enforcement, while secondary enforcement raises it by about **10.9 percentage points**. The difference (9.7 pp) is also precisely estimated ($p=0.0001$). This makes intuitive sense: primary enforcement allows police to stop a car solely for a belt violation, creating a direct deterrent. Under secondary enforcement, a belt ticket can only be added on top of another stop, so its deterrent effect is weaker.

## 10. Lives saved in New Jersey (2 pts)

In 2000, New Jersey switched from secondary to primary enforcement. We combine the usage and fatality estimates in a two-step calculation.

**Step 1 — Change in belt usage:**

$$
\Delta\,\text{sb\_useage}_{NJ}
= \hat\beta_{\text{primary}}-\hat\beta_{\text{secondary}}
= 0.2056-0.1085
= 0.0971.
$$

Switching to primary enforcement raises belt use by about 9.7 percentage points.

**Step 2 — Change in fatality rate** (using two-way FE from Q4):

$$
\Delta\,\text{fatalityrate}_{NJ}
= \hat\beta_{sb} \times \Delta\,\text{sb\_useage}_{NJ}
= -0.0037 \times 0.0971
= -0.000361 \text{ per million miles}.
$$

**Step 3 — Lives saved.** NJ's 1997 VMT (closest observed year) is 63,308 million miles.

$$
\text{Lives saved}
= |\Delta\,\text{fatalityrate}_{NJ}| \times \text{VMT}_{NJ}
= 0.000361 \times 63{,}308
\approx \mathbf{22.9\;\text{lives per year}}.
$$

So the regulatory change is estimated to prevent roughly **23 fatalities per year** in New Jersey. That corresponds to about 3% of NJ's annual traffic deaths — a meaningful gain from what is essentially a change in how an existing law is enforced. The calculation relies on extrapolating 1997 traffic volumes to 2000, so it should be treated as an informed approximation rather than a precise forecast.

---

## Summary — Fatality-Rate Regressions

| Variable       | Pooled OLS | State FE | State+Year FE | Random Effects |
| -------------- | ---------: | -------: | ------------: | -------------: |
| sb_useage      |     0.0041 | −0.0058 |      −0.0037 |       −0.0045 |
|                |   (0.0024) | (0.0017) |      (0.0015) |       (0.0016) |
| speed65        |     0.0001 | −0.0004 |      −0.0008 |       −0.0003 |
| speed70        |     0.0024 |   0.0012 |        0.0008 |         0.0013 |
| ba08           |   −0.0019 | −0.0014 |      −0.0008 |       −0.0014 |
| drinkage21     |     0.0001 |   0.0007 |      −0.0011 |         0.0008 |
| ln_income      |   −0.0181 | −0.0135 |        0.0063 |       −0.0126 |
| age            |   −0.0000 |   0.0010 |        0.0013 |         0.0002 |
| State FE       |         No |      Yes |           Yes |             RE |
| Year FE        |         No |       No |           Yes |             No |
| N              |        556 |      556 |           556 |            556 |
| $R^2$ within |         — |    0.687 |         0.751 |             — |

*Cluster-robust SEs in parentheses (clustered at state level).*
