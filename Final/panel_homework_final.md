---
title: "Panel Data Homework"
author: "Dragos Florin Vasile & Wong Hei Wong"
date: "March 2026"
output:
  pdf_document: default
  html_document:
    df_print: paged
fontsize: 11pt
geometry: margin=1in
header-includes:
- \usepackage{booktabs}
- \usepackage{amsmath}
- \usepackage{amssymb}
- \usepackage{float}
- \usepackage{xcolor}
- "\\DeclareUnicodeCharacter{2212}{\\ensuremath{-}}"
---

# Exercise 1

## Question 1

We treat the $\alpha_i$'s as fixed, unknown constants and solve the least-squares problem with respect to both $\beta$ and $\alpha_i$, leading to the following objective function.

$$
\min_{\beta,\,\alpha_1,\ldots,\alpha_N}
\sum_{i=1}^{N}\sum_{t\in P_i}
\left(y_{it}-\alpha_i-\beta x_{it}\right)^2.
$$

**Part 1**. Deriving the within estimator $\hat{\alpha}_{iW}$

We solve the FOC with respect to $\alpha_i$, and setting to 0 gives us:

$$ \frac{\partial SSR}{\partial \alpha_{i}} = -2 \sum_{t\in P_{i}}(y_{it} - \alpha_{i} - \beta x_{it}) = 0 $$

We can simplify the terms into

$$\sum_{t\in P_{i}}y_{it} - \sum_{t\in P_{i}}\alpha_{i} - \beta \sum_{t\in P_{i}}x_{it} = 0$$

Since $\alpha_{i}$ does not depend on $t$, summing it $T_i$ times simply gives us $T_i \alpha_i$, dividing it by $T_i$ gives us,

$$\frac{1}{T_{i}}\sum_{t\in P_{i}}y_{it} - \alpha_{i} - \beta \frac{1}{T_{i}}\sum_{t\in P_{i}}x_{it} = 0$$

Using the definition for individual-specific mean, we rewrite and rearrange the terms, which gives us the following,

$$\frac{1}{T_{i}}\sum_{t\in P_{i}}y_{it} - \alpha_{i} - \beta \frac{1}{T_{i}}\sum_{t\in P_{i}}x_{it} = 0$$

$$\overline{y}_{i.} - \alpha_{i} - \beta \overline{x}_{i.} = 0$$

$$\hat{\alpha}_{iW} = \overline{y}_{i.} - \hat{\beta}_{W}\overline{x}_{i.} \quad \square$$

**Part 2**. Deriving the within estimator $\hat{\beta}_W$

We follow a similar procedure and solve the FOC with respect to $\beta$ this time. Using the results from the previous part we can rewrite the residual terms,

$$e_{it} = y_{it} - (\overline{y}_{i.} - \beta \overline{x}_{i.}) - \beta x_{it}$$

$$e_{it} = (y_{it} - \overline{y}_{i.}) - \beta(x_{it} - \overline{x}_{i.})$$

which gives us

$$\frac{\partial SSR^{*}}{\partial \beta} = -2 \sum_{i=1}^{N}\sum_{t\in P_{i}}(x_{it} - \overline{x}_{i.})[(y_{it} - \overline{y}_{i.}) - \beta(x_{it} - \overline{x}_{i.})] = 0$$

Using the given notation for $W_{ZQ}$, in a similar fashion, we expand the terms leading to the following,

$$\sum_{i=1}^{N}\sum_{t\in P_{i}}(x_{it} - \overline{x}_{i.})(y_{it} - \overline{y}_{i.}) - \beta \sum_{i=1}^{N}\sum_{t\in P_{i}}(x_{it} - \overline{x}_{i.})^{2} = 0$$

$$W_{XY} - \beta W_{XX} = 0$$

$$\hat{\beta}_{W} = W_{XX}^{-1}W_{XY} \quad \square$$

## Question 2

The form of the within estimator is algebraically identical for both balanced and unbalanced panels. The only adjustment is that when calculating individual-specific means ($\overline{y}_{i.}$ and $\overline{x}_{i.}$), we divide by the specific number of observations available for that individual, $T_i = \sum_{t=1}^{T} d_{it}$, rather than a fixed $T$.

Successive observations for each individual are not required. The within transformation only requires that an individual is observed at least twice ($T_i\ge 2$) and that $x_{it}$ varies within that individual, so that the individual contributes non-zero within variation to $W_{XX}$. Non-consecutive years are not an issue in our estimates since they do not depend on the continuity of the observations.

## Question 3

**Part 1**. Deriving the Variance for $\text{Var}(W_{XU}\mid X)$

Using the definition for $W_{ZQ}$, we can define the $W_{XU}$ term as

$$W_{XU} = \sum_{i=1}^{N}\sum_{t\in P_{i}}(x_{it}-\overline{x}_{i.})u_{it}$$

so we have

$$\text{Var}(W_{XU}\mid X) = \text{Var}\left( \sum_{i=1}^{N}\sum_{t\in P_{i}}(x_{it}-\overline{x}_{i.})u_{it} \bigg| X \right)$$

because we assume $u_{it}\mid X \sim IID(0,\sigma^{2})$, we have

$$\text{Var}(W_{XU}\mid X) = \sum_{i=1}^{N}\sum_{t\in P_{i}}(x_{it}-\overline{x}_{i.})^2 \text{Var}(u_{it}\mid X)$$

$$\text{Var}(W_{XU}\mid X) = \sigma^2 \sum_{i=1}^{N}\sum_{t\in P_{i}}(x_{it}-\overline{x}_{i.})^2$$

Using the definition for $W_{ZQ}$ again, we can collapse the latter term into the following

$$\text{Var}(W_{XU}\mid X) = \sigma^2 W_{XX}$$

**Part 2**. Unbiasedness

The within estimator we derived in Q1 is $\hat{\beta}_{W} = W_{XX}^{-1}W_{XY}$, so using the definition we have

$$W_{XY} = \sum_{i=1}^{N}\sum_{t\in P_{i}}(x_{it}-\overline{x}_{i.})((\alpha_i + \beta x_{it} + u_{it}) - (\alpha_i + \beta \overline{x}_{i.} + \overline{u}_{i.}))$$

The $\alpha_i$ terms cancel out, and we can further simplify the terms,


$$W_{XY} = \beta W_{XX} + W_{XU}$$

where

$$
\sum_{t\in P_i}(x_{it}-\overline{x}_{i.})\overline{u}_{i.}
=
\overline{u}_{i.}\sum_{t\in P_i}(x_{it}-\overline{x}_{i.})
=0,
$$

so the individual-mean error term drops out.

Plug it back to our original estimator and take expectation conditional on $X$,

$$\hat{\beta}_{W} = W_{XX}^{-1}(\beta W_{XX} + W_{XU}) = \beta + W_{XX}^{-1}W_{XU}$$

$$E(\hat{\beta}_{W}\mid X) = \beta + W_{XX}^{-1} E(W_{XU}\mid X)$$

Since $E(u_{it}\mid X) = 0$, $E(W_{XU}\mid X) = 0$, therefore $E(\hat{\beta}_{W}\mid X) = \beta$ and the estimator is unbiased. $\square$

**Part 3**. Deriving the variance

Using the property that $\text{Var}(AY\mid X) = A\,\text{Var}(Y\mid X)\,A'$, we have

$$\text{Var}(\hat{\beta}_{W}\mid X) = W_{XX}^{-1} \text{Var}(W_{XU}\mid X) W_{XX}^{-1}$$

Using the earlier result in Part 1, we have

$$\text{Var}(\hat{\beta}_{W}\mid X) = W_{XX}^{-1} (\sigma^2 W_{XX}) W_{XX}^{-1}$$

$$\text{Var}(\hat{\beta}_{W}\mid X) = \sigma^2 W_{XX}^{-1} W_{XX} W_{XX}^{-1}$$

$$\text{Var}(\hat{\beta}_{W}\mid X) = \sigma^2 W_{XX}^{-1} \quad \square$$ 


## Question 4

From the between equation, we have the averaged error term

$$\overline{u}_{i.} = \frac{1}{T_{i}} \sum_{t\in P_{i}} u_{it}$$

Since we assume that $(u_{it}|X) \sim IID(0, \sigma^2)$, the conditional variance becomes,

$$\text{Var}(\overline{u}_{i.}\mid X) = \text{Var}\left( \frac{1}{T_{i}} \sum_{t\in P_{i}} u_{it} \bigg| X \right) = \frac{1}{T_{i}^2} \sum_{t\in P_{i}} \text{Var}(u_{it}\mid X)$$

Given independent and homoskedastic errors we have,

$$\text{Var}(\overline{u}_{i.}\mid X) = \frac{T_{i} \sigma^2}{T_{i}^2} = \frac{\sigma^2}{T_{i}}$$

Therefore, if we have an unbalanced panel where $T_i$ varies across individuals, the error variance is not constant across individuals since it depends on $T_i$, so we will have a heteroskedastic error term. It is only homoskedastic under a balanced panel.

## Question 5

We can get rid of the heteroskedasticity by applying a WLS correction. By multiplying the equation by $\sqrt{T_i}$

$$\sqrt{T_i}\overline{y}_{i.} = \sqrt{T_i}\alpha_{i} + \beta\sqrt{T_i}\overline{x}_{i.} + \sqrt{T_i}\overline{u}_{i.}$$

The transformed error has variance

$$\text{Var}(\sqrt{T_i}\overline{u}_{i.}\mid X) = T_i \cdot \text{Var}(\overline{u}_{i.}\mid X) = T_i \cdot \frac{\sigma^2}{T_i} = \sigma^2$$

which is now homoskedastic again.

## Question 6

We minimize the Sum of Squared Residuals in the new equation with respect to $\alpha$ and $\beta$

$$\min_{\alpha, \beta} \sum_{i=1}^{N} (\bar{y}_{i.} - \alpha - \beta \bar{x}_{i.})^2$$

**Part 1**. Deriving the between estimator $\hat{\alpha}_{B1}$

Taking the partial derivative with respect to $\alpha$, set it to zero and rearrange the terms,

$$\frac{\partial SSR}{\partial \alpha} = -2 \sum_{i=1}^{N} (\bar{y}_{i.} - \alpha - \beta \bar{x}_{i.}) = 0$$

$$\sum_{i=1}^{N} \bar{y}_{i.} - \sum_{i=1}^{N} \alpha - \beta \sum_{i=1}^{N} \bar{x}_{i.} = 0$$


Dividing the whole equation by $N$ gives

$$\frac{1}{N}\sum_{i=1}^{N} \bar{y}_{i.} - \alpha - \beta \frac{1}{N}\sum_{i=1}^{N} \bar{x}_{i.} = 0$$

So we have, 

$$\hat{\alpha}_{B1} = \bar{\bar{y}} - \hat{\beta}_{B1} \bar{\bar{x}}$$

**Part 2**. Deriving the between estimator $\hat{\beta}_{B1}$

Taking the partial derivative with respect to $\beta$, setting it to zero and rearranging the terms, and using the result from the previous part, we have,

$$\frac{\partial SSR}{\partial \beta} = -2 \sum_{i=1}^{N} \bar{x}_{i.} (\bar{y}_{i.} - \alpha - \beta \bar{x}_{i.}) = 0$$

$$\sum_{i=1}^{N} \bar{x}_{i.} [(\bar{y}_{i.} - \bar{\bar{y}}) - \beta (\bar{x}_{i.} - \bar{\bar{x}})] = 0$$

$$\sum_{i=1}^{N} (\overline{x}_{i.} - \overline{\overline{x}}) \left[ (\overline{y}_{i.} - \overline{\overline{y}}) - \beta (\overline{x}_{i.} - \overline{\overline{x}}) \right] = 0$$

$$\sum_{i=1}^{N} (\overline{x}_{i.} - \overline{\overline{x}})(\overline{y}_{i.} - \overline{\overline{y}}) - \beta \sum_{i=1}^{N} (\overline{x}_{i.} - \overline{\overline{x}})^2 = 0$$

So we have,

$$\hat{\beta}_{B1} = \frac{\sum_{i=1}^{N} (\overline{x}_{i.} - \overline{\overline{x}})(\overline{y}_{i.} - \overline{\overline{y}})}{\sum_{i=1}^{N} (\overline{x}_{i.} - \overline{\overline{x}})^2}$$

## Question 7

The procedure is nearly the same, with our new objective function becoming,

$$\min_{\alpha, \beta} \sum_{i=1}^{N} (\sqrt{T_i}\bar{y}_{i.} - \sqrt{T_i}\alpha - \beta \sqrt{T_i}\bar{x}_{i.})^2$$

Or equivalently, 

$$\min_{\alpha, \beta} \sum_{i=1}^{N} T_i(\bar{y}_{i.} - \alpha - \beta \bar{x}_{i.})^2$$

**Part 1**. Deriving the between estimator $\hat{\alpha}_{B2}$

Taking the partial derivative with respect to $\alpha$, set it to zero and rearrange the terms,

$$\frac{\partial SSR}{\partial \alpha} = -2 \sum_{i=1}^{N} T_i(\overline{y}_{i.} - \alpha - \beta \overline{x}_{i.}) = 0$$

$$\sum_{i=1}^{N} T_i \overline{y}_{i.} - \alpha \sum_{i=1}^{N} T_i - \beta \sum_{i=1}^{N} T_i \overline{x}_{i.} = 0$$

Now instead of dividing it by $N$, we divide it by $\sum_{i=1}^{N} T_i$ to get the following,

$$\frac{\sum T_i \overline{y}_{i.}}{\sum T_i} - \alpha - \beta \frac{\sum T_i \overline{x}_{i.}}{\sum T_i} = 0$$

Using the notation for global mean, we have,

$$\hat{\alpha}_{B2} = \overline{y} - \hat{\beta}_{B2} \overline{x}$$

**Part 2**. Deriving the between estimator $\hat{\beta}_{B2}$

Taking the partial derivative with respect to $\beta$, setting it to zero and rearranging the terms, and using the result from the previous part, we have,

$$\frac{\partial SSR}{\partial \beta} = -2 \sum_{i=1}^{N} T_i \overline{x}_{i.} (\overline{y}_{i.} - \alpha - \beta \overline{x}_{i.}) = 0$$

$$\sum_{i=1}^{N} T_i \overline{x}_{i.} [(\overline{y}_{i.} - \overline{y}) - \beta (\overline{x}_{i.} - \overline{x})] = 0$$

$$\sum_{i=1}^{N} T_i (\overline{x}_{i.} - \overline{x}) [(\overline{y}_{i.} - \overline{y}) - \beta (\overline{x}_{i.} - \overline{x})] = 0$$

$$\sum_{i=1}^{N} T_i (\overline{x}_{i.} - \overline{x})(\overline{y}_{i.} - \overline{y}) = \beta \sum_{i=1}^{N} T_i (\overline{x}_{i.} - \overline{x})^2$$

So we have,

$$\hat{\beta}_{B2} = \frac{\sum_{i=1}^{N} T_i (\overline{x}_{i.} - \overline{x})(\overline{y}_{i.} - \overline{y})}{\sum_{i=1}^{N} T_i (\overline{x}_{i.} - \overline{x})^2}$$

## Question 8

In a balanced panel ($T_i=T$ for all $i$), the weights are all equal and cancel from the ratio, so $\hat\beta_{B1}=\hat\beta_{B2}$ and $\hat\alpha_{B1}=\hat\alpha_{B2}$. In an unbalanced panel, $\hat\beta_{B1}$ gives every individual the same weight regardless of how many periods it was observed, while $\hat\beta_{B2}$ weights individuals with more observations more. Since the error variance is $\sigma^2/T_i$, $\hat\beta_{B2}$ is the GLS-efficient (BLUE) estimator in the between equation. The practical takeaway is that when the panel is unbalanced, ignoring the heteroskedasticity in the between model leads to a loss of efficiency.

\newpage

# Exercise 2

All regressions below use cluster-robust standard errors at the state level, which allows for arbitrary heteroskedasticity and serial correlation within states. The code is written in Stata. Switch to your working directory in the do-file.

## Question 1

The raw file contains 765 state-year observations for 50 U.S. states plus D.C. in the period 1983–1997 (15 years). The preprocessing steps are:

``` stata
import delimited "../SeatBelts.csv", delimiter(";") clear
isid state year
xtset fips year
gen ln_income = ln(income)
gen sample_main = !missing(fatalityrate, sb_useage, speed65, speed70, ///
    ba08, drinkage21, ln_income, age)
```

`isid state year` checks whether the combination of state and year is unique in our dataset, and it is; `xtset fips year` declares the panel data structure. `sb_useage` is missing for 209 observations, concentrated in the early years of the sample (around 1983–1989). All other variables are complete. Because `sb_useage` is a key regressor, we restrict the estimation sample and keep 556 observations with non-missing values. This makes the sample unbalanced, but as shown in Exercise 1 this does not create any bias. A standard log transformation for income measured in dollars is generated in `ln_income = ln(income)`.

| Statistic           |     Value |
|---------------------|----------:|
| Raw observations    |       765 |
| States (incl. DC)   |        51 |
| Years               | 1983–1997 |
| Missing `sb_useage` |       209 |
| Estimation sample   |       556 |
| Mean `fatalityrate` |    0.0198 |
| Mean `sb_useage`    |     0.529 |
| Mean `income`       |  \$19,572 |

## Question 2

``` stata
reg fatalityrate sb_useage speed65 speed70 ba08 drinkage21 ///
    ln_income age if sample_main, vce(cluster fips)
```

| Variable   |  Coeff. | Cl. SE |     t |     p |
|------------|--------:|-------:|------:|------:|
| sb_useage  |  0.0041 | 0.0024 |  1.70 | 0.095 |
| speed65    |  0.0001 | 0.0006 |  0.25 | 0.807 |
| speed70    |  0.0024 | 0.0007 |  3.46 | 0.001 |
| ba08       | −0.0019 | 0.0007 | −2.60 | 0.012 |
| drinkage21 |  0.0001 | 0.0013 |  0.06 | 0.950 |
| ln_income  | −0.0181 | 0.0024 | −7.47 | 0.000 |
| age        | −0.0000 | 0.0004 | −0.02 | 0.987 |
| Constant   |  0.1965 | 0.0198 |  9.92 | 0.000 |

$N=556$, $R^2=0.549$.

The pooled OLS coefficient on sb_useage is positive (0.0041) and only marginally significant. This implies that a one-unit increase in seat belt usage (shifting from 0% to 100% usage) is associated with an increase of 0.0041 fatalities per million traffic miles, ceteris paribus. Since the mean fatality rate is 0.0198, an increase of 0.0041 represents roughly a 20.7% increase in the fatality rate relative to the mean. This estimate suggests higher belt use is associated with more fatalities, which is not credible. This is likely the result of omitted variable bias. States with historically high fatality rates (due to dangerous rural road networks, driving culture, etc.) may have adopted seat belt laws precisely because of those problems. Pooled OLS mixes cross-state differences with over-time variation and therefore picks up this selection effect.

## Question 3

``` stata
xtreg fatalityrate sb_useage speed65 speed70 ba08 drinkage21 ///
    ln_income age if sample_main, fe vce(cluster fips)
```

| Variable   |  Coeff. | Cl. SE |     t |     p |
|------------|--------:|-------:|------:|------:|
| sb_useage  | −0.0058 | 0.0017 | −3.46 | 0.001 |
| speed65    | −0.0004 | 0.0005 | −0.93 | 0.355 |
| speed70    |  0.0012 | 0.0003 |  3.54 | 0.001 |
| ba08       | −0.0014 | 0.0004 | −3.67 | 0.001 |
| drinkage21 |  0.0007 | 0.0007 |  1.04 | 0.305 |
| ln_income  | −0.0135 | 0.0024 | −5.67 | 0.000 |
| age        |  0.0010 | 0.0007 |  1.31 | 0.196 |

$R^2_{\text{within}}=0.687$, $\hat\sigma_u=0.00383$, $\hat\sigma_e=0.00179$, $\hat\rho=0.821$.

As predicted, the results completely change. Now the sign on sb_useage reverses to −0.0058 and is highly statistically significant. This implies that a one-unit increase in seat belt usage (shifting from 0% to 100% usage) is associated with a decrease of 0.0058 fatalities per million traffic miles, ceteris paribus. Relative to the mean fatality rate, this corresponds to a 29% reduction in fatalities. Once we control for state fixed effects, our identification now comes from within-state changes over time. Within a given state, years with higher seat belt use see fewer fatalities. The high $\hat\rho=0.82$ confirms that state-level unobservables account for the majority of residual variance, which makes controlling for them important.

## Question 4

``` stata
xtreg fatalityrate sb_useage speed65 speed70 ba08 drinkage21 ///
    ln_income age i.year if sample_main, fe vce(cluster fips)
```

| Variable   |  Coeff. | Cl. SE |     t |     p |
|------------|--------:|-------:|------:|------:|
| sb_useage  | −0.0037 | 0.0015 | −2.56 | 0.013 |
| speed65    | −0.0008 | 0.0006 | −1.35 | 0.183 |
| speed70    |  0.0008 | 0.0005 |  1.76 | 0.085 |
| ba08       | −0.0008 | 0.0004 | −1.86 | 0.069 |
| drinkage21 | −0.0011 | 0.0006 | −1.82 | 0.074 |
| ln_income  |  0.0063 | 0.0067 |  0.94 | 0.354 |
| age        |  0.0013 | 0.0007 |  1.90 | 0.063 |

Joint test on year dummies: $F(14,50)=9.82$, $p < 0.001$. $R^2_{\text{within}}=0.751$.

Adding year fixed effects absorbs national trends and common shocks that affect all states simultaneously. Examples include improvements in vehicle safety designs, nationwide road safety campaigns, and advancements in trauma care. The year dummies are jointly significant ($p < 0.001$), confirming that accounting for these common time trends is empirically necessary. By including them, the coefficient on `sb_useage` shrinks from −0.0058 (in Q3) to −0.0037. This suggests that the estimate in Q3 was slightly biased because it attributed the general national downward trend in fatalities to seat belt usage. After removing these trends, the remaining within-state, within-year variation still identifies a negative and statistically significant effect of seat belt use on fatalities. `ln_income` becomes insignificant because its over-time growth is now largely captured by the year fixed effects.

## Question 5

|                    | Pooled OLS |    State FE    |  State + Year FE  |
|--------------------|:----------:|:--------------:|:-----------------:|
| `sb_useage` coeff. |  +0.0041   |    −0.0058     |      −0.0037      |
| $R^2$              |   0.549    | 0.687 (within) |  0.751 (within)   |
| Year FE            |     No     |       No       | Yes ($p < 0.001$) |

The two-way fixed effects (state + year) specification is the most reliable because it addresses the core biases present in the other models. Pooled OLS clearly suffers from omitted variable bias due to unobserved state-level heterogeneity, which is evidenced by the positive and non-credible sign on the `sb_useage` coefficient. While the State FE model fixes this state-level confounding by controlling for state fixed effects, it fails to account for aggregate national time trends shared across all states. The final two-way FE specification removes both state heterogeneity and common time-trend effects, such as nationwide improvements in vehicle safety or trauma care. The strongly significant year dummies ($p < 0.001$) confirm that time effects are empirically important, and omitting them would likely bias our estimates. More importantly, the two-way specification leaves identification to variation that differs across states and over time simultaneously, which is arguably the most credible source of variation for this study.

## Question 6

In Q4, our estimate $\hat{\beta}_{sb} = -0.0037$ indicates that a full implementation of seat belt usage (shifting from a proportion of 0 to 1) reduces the fatality rate by 0.0037 deaths per million traffic miles. Specifically, a 1-percentage-point increase ($\Delta = 0.01$) reduces fatalities by approximately 0.19% relative to the mean fatality rate of 0.0198 ($\frac{0.0037 \times 0.01}{0.0198} \approx 0.00187$). In our scenario, a shift from 52% to 90% usage ($\Delta = 0.38$) predicts a reduction of 0.00141 fatalities per million miles. Given an average aggregate U.S. traffic volume of 2,113,834 million miles during 1983–1997, this translates to roughly 2,987 lives saved per year. This represents a substantial 7.15% reduction in total road deaths across the sampling period relative to the mean. This confirms that seat belt usage is a highly effective, policy-relevant intervention for public safety.

## Question 7

``` stata
xtreg fatalityrate sb_useage speed65 speed70 ba08 drinkage21 ///
    ln_income age if sample_main, re vce(cluster fips)
```

| Variable   |  Coeff. | Cl. SE |     z |     p |
|------------|--------:|-------:|------:|------:|
| sb_useage  | −0.0045 | 0.0016 | −2.73 | 0.006 |
| speed65    | −0.0003 | 0.0004 | −0.76 | 0.447 |
| speed70    |  0.0013 | 0.0004 |  3.80 | 0.000 |
| ba08       | −0.0014 | 0.0004 | −3.66 | 0.000 |
| drinkage21 |  0.0008 | 0.0007 |  1.09 | 0.278 |
| ln_income  | −0.0126 | 0.0020 | −6.36 | 0.000 |
| age        |  0.0002 | 0.0005 |  0.45 | 0.653 |

$\hat\sigma_u=0.00302$, $\hat\sigma_e=0.00179$, $\hat\rho=0.740$.

This RE specification corresponds to Q3 (state FE, no year dummies). The critical assumption is that the state effect $\alpha_i$ is uncorrelated with the regressors:

$$
E(\alpha_i \mid x_{i1},\dots,x_{iT}) = 0.
$$

If this orthogonality assumption holds, RE is more efficient than fixed effects because it utilizes both within-state and between-state variation. In this model, the coefficient on `sb_useage` is −0.0045, which sits between the FE estimate (−0.0058) and the pooled OLS estimate (+0.0041). This result is expected because the RE estimator is a matrix-weighted average of the within and between estimators. It partially reintroduces the cross-state variation that was identified as problematic in the pooled model due to selection effects.

## Question 8

The classic Hausman test assumes that the random-effects estimator is fully efficient under the null hypothesis, which requires the errors to be homoskedastic. Given the nature of state-level traffic fatality data, our model likely suffers from heteroskedasticity and within-state serial correlation. Therefore, to perform a robust version of the Hausman test, we used the Sargan-Hansen overidentification test (using the `xtoverid` command in Stata). We run this test on the same specification as Q7 (no year dummies) so the FE-RE comparison is fully consistent.

``` stata
xtreg fatalityrate sb_useage speed65 speed70 ba08 drinkage21 ///
    ln_income age if sample_main, re vce(cluster fips)
xtoverid
```

The hypotheses are:

$H_0$: The unobserved state-specific effects are strictly uncorrelated with the explanatory variables (both RE and FE are consistent, but RE is more efficient).

$H_1$: The unobserved state-specific effects are correlated with the explanatory variables (RE is inconsistent, FE is consistent).

Results: $\chi^2(7) = 28.421$, p-value = $0.0002$.

Because the p-value is less than the standard 0.05 significance level, we reject the null hypothesis and deduce that the unobserved state heterogeneity is correlated with our regressors. Consequently, the random-effects estimator is inconsistent. We should use the fixed-effects (FE) model, as it safely controls for this time-invariant unobserved heterogeneity.

## Question 9

``` stata
xtreg sb_useage primary secondary speed65 speed70 ba08 drinkage21 ///
    ln_income age i.year if sample_main, fe vce(cluster fips)
lincom primary - secondary
```

| Variable  | Coeff. | Cl. SE |    t |       p |
|-----------|-------:|-------:|-----:|--------:|
| primary   | 0.2056 | 0.0232 | 8.87 | 0.000 |
| secondary | 0.1085 | 0.0134 | 8.09 | 0.000 |
| speed65   | 0.0228 | 0.0205 | 1.11 | 0.271 |
| speed70   | 0.0120 | 0.0206 | 0.58 | 0.561 |
| ba08      | 0.0037 | 0.0176 | 0.21 | 0.832 |
| drinkage21 | 0.0107 | 0.0272 | 0.39 | 0.695 |
| ln_income | 0.0583 | 0.2564 | 0.23 | 0.821 |
| age       | 0.0138 | 0.0231 | 0.60 | 0.553 |

$R^2_{\text{within}}=0.842$. $\hat\sigma_u=0.079$, $\hat\sigma_e=0.057$, $\hat\rho=0.656$.

The regression results indicate that both primary and secondary enforcement laws significantly increase seat belt usage compared to the absence of such mandates. Primary enforcement is associated with a 20.6 percentage point increase in usage, while secondary enforcement yields a lower 10.9 percentage point increase, with both effects being highly statistically significant at the 1% level. For the separate equality test, we use $H_0: \beta_{primary} - \beta_{secondary} = 0$ (from `lincom primary - secondary`), and we reject it at conventional levels ($p < 0.001$). This shows that primary enforcement is more effective than secondary enforcement at driving compliance. This suggests that the stricter legal mandate leads to significantly higher seat belt usage.

## Question 10

To estimate the number of lives saved, we combine the behavioral impact of the law with its mechanical effect on mortality. First, the transition from secondary to primary enforcement yields a net increase in seat belt usage of 9.71 percentage points ($\hat\beta_{primary} - \hat\beta_{secondary}$). Multiplying this usage gain by the fatality rate coefficient from the two-way fixed effects model ($\hat\beta_{sb} = -0.0037$) results in a reduction of approximately 0.000361 deaths per million miles. When applied to New Jersey's last observed volume of 63,308 million miles, the change is estimated to prevent roughly 23 fatalities per year. This represents a 2.4% reduction relative to New Jersey's average fatality rate.

------------------------------------------------------------------------

## Summary of the results in our regressions

| Variable     | Pooled OLS | State FE | State + Year FE | Random Effects |
|--------------|-----------:|---------:|--------------:|---------------:|
| sb_useage    |     0.0041 |  −0.0058 |       −0.0037 |        −0.0045 |
|              |   (0.0024) | (0.0017) |      (0.0015) |       (0.0016) |
| speed65      |     0.0001 |  −0.0004 |       −0.0008 |        −0.0003 |
| speed70      |     0.0024 |   0.0012 |        0.0008 |         0.0013 |
| ba08         |    −0.0019 |  −0.0014 |       −0.0008 |        −0.0014 |
| drinkage21   |     0.0001 |   0.0007 |       −0.0011 |         0.0008 |
| ln_income    |    −0.0181 |  −0.0135 |        0.0063 |        −0.0126 |
| age          |    −0.0000 |   0.0010 |        0.0013 |         0.0002 |
| State FE     |         No |      Yes |           Yes |             RE |
| Year FE      |         No |       No |           Yes |             No |
| N            |        556 |      556 |           556 |            556 |
| $R^2$ within |          — |    0.687 |         0.751 |              — |

*Cluster-robust SEs in parentheses (clustered at state level).*
