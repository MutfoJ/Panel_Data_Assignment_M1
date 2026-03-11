*==============================================================================*
*  M1 TSE — Panel Data Homework 2025–2026                                     *
*  Pooled OLS, Fixed Effects, and Random Effects                               *
*  Authors: Dragos Florin Vasile & Wong Hei Wong                               *
*  Date:    March 2026                                                         *
*==============================================================================*

/*  INSTRUCTIONS
    ============
    This file replicates every empirical result reported in the accompanying
    answers sheet.  It is organised into numbered blocks that match the
    assignment questions.  All output goes to the log file so the reader can
    cross-check numbers easily.

    - Keep this .do inside the "Final" subfolder.
    - The data file SeatBelts.csv should sit one level up (i.e. ../).
    - The log is written to the same Final subfolder.
*/

clear all
set more off
capture log close _all

* ---------- Set working directory to Final subfolder -------------------------
cd "C:\Users\vdrag\OneDrive\Escritorio\Panel Data Project\Final"

log using "panel_homework_final.log", replace text

* ---------- Path to data -----------------------------------------------------
local datafile "../SeatBelts.csv"
capture confirm file "`datafile'"
if _rc {
    di as error "SeatBelts.csv not found -- keep this .do inside Final/."
    exit 601
}

* A local that stores all the time-varying controls we use throughout.
local xvars sb_useage speed65 speed70 ba08 drinkage21 ln_income age


*==============================================================================*
*  QUESTION 1 — Data description and cleaning                                  *
*==============================================================================*

/*  We import the semicolon-delimited CSV, confirm there are no duplicates,
    declare the panel structure, generate log-income, and define the main
    estimation sample.  The only cleaning issue is that sb_useage is missing
    for 209 state-years, mostly 1983–1989.  All regressions that use
    sb_useage are therefore restricted to the 556 complete-case observations.
*/

import delimited "`datafile'", delimiter(";") varnames(1) clear

* Verify unique state-year identifiers and set panel
isid state year
xtset fips year

* Generate log-income (natural functional form for a dollar-denominated variable)
gen ln_income = ln(income)
label var ln_income "log(income)"

* Mark the estimation sample — keeps rows where nothing relevant is missing
gen sample_main = !missing(fatalityrate, sb_useage, speed65, speed70, ///
    ba08, drinkage21, ln_income, age)
label var sample_main "Complete-case sample (556 obs)"

di _newline(2) as text "=============================================="
di as text "  QUESTION 1: DATA OVERVIEW"
di as text "=============================================="

count
di as result "Raw observations: " r(N)

qui count if sample_main
di as result "Estimation sample: " r(N)

qui count if missing(sb_useage)
di as result "Missing sb_useage: " r(N)

tab year if missing(sb_useage)
xtdescribe if sample_main
summarize fatalityrate sb_useage speed65 speed70 ba08 drinkage21 ///
    income ln_income age primary secondary vmt if sample_main

/*  We also compute national traffic volumes here because they will be
    needed later (Q6, Q10) to convert fatality-rate changes into lives.
*/
preserve
collapse (sum) us_vmt = vmt, by(year)
qui summarize us_vmt
scalar mean_us_vmt = r(mean)                   // average annual US VMT
qui summarize us_vmt if year == 1997
scalar us_vmt_1997 = r(mean)                   // 1997 US VMT
restore

qui summarize vmt if sample_main
scalar mean_state_vmt = r(mean)                // average state-year VMT


*==============================================================================*
*  QUESTION 2 — Pooled OLS                                                     *
*==============================================================================*

/*  Baseline cross-sectional regression.  We expect the coefficient on
    sb_useage to be misleading because of omitted state-level confounders.
*/

di _newline(2) as text "=============================================="
di as text "  QUESTION 2: POOLED OLS"
di as text "=============================================="

reg fatalityrate `xvars' if sample_main, vce(cluster fips)
estimates store pooled

scalar b_pool = _b[sb_useage]
scalar se_pool = _se[sb_useage]


*==============================================================================*
*  QUESTION 3 — State fixed effects                                            *
*==============================================================================*

/*  Adding state FE removes all time-invariant unobserved heterogeneity.
    Identification now comes from within-state variation over time.
*/

di _newline(2) as text "=============================================="
di as text "  QUESTION 3: STATE FIXED EFFECTS"
di as text "=============================================="

xtreg fatalityrate `xvars' if sample_main, fe vce(cluster fips)
estimates store fe_state

scalar b_fe = _b[sb_useage]
scalar se_fe = _se[sb_useage]


*==============================================================================*
*  QUESTION 4 — State + year fixed effects (two-way FE)                        *
*==============================================================================*

/*  Year dummies absorb national trends in fatality rates (safer cars,
    federal campaigns, etc.).  We test their joint significance.
*/

di _newline(2) as text "=============================================="
di as text "  QUESTION 4: TWO-WAY FIXED EFFECTS"
di as text "=============================================="

xtreg fatalityrate `xvars' i.year if sample_main, fe vce(cluster fips)
estimates store fe_tw

scalar b_tw = _b[sb_useage]
scalar se_tw = _se[sb_useage]

testparm i.year
scalar p_yearFE = r(p)
di as result "Joint test on year dummies: p = " %7.5f p_yearFE


*==============================================================================*
*  QUESTION 6 — Size of the sb_useage coefficient and lives saved              *
*==============================================================================*

/*  We translate the two-way FE estimate into an effect size: how many
    lives would be saved if seat belt use rose from 52 % to 90 %?
*/

di _newline(2) as text "=============================================="
di as text "  QUESTION 6: EFFECT SIZE"
di as text "=============================================="

scalar delta_sb = 0.90 - 0.52
scalar delta_rate = b_tw * delta_sb

qui summarize fatalityrate if sample_main
scalar mean_fatal = r(mean)
scalar pct_mean  = -100 * delta_rate / mean_fatal

scalar lives_us_mean = -delta_rate * mean_us_vmt
scalar lives_us_1997 = -delta_rate * us_vmt_1997

di as result "Two-way FE coefficient     = " %9.6f b_tw
di as result "Delta fatality rate (52→90) = " %9.6f delta_rate
di as result "% of mean fatality rate     = " %6.2f pct_mean "%"
di as result "Lives saved/year (avg VMT)  = " %9.1f lives_us_mean
di as result "Lives saved/year (1997 VMT) = " %9.1f lives_us_1997


*==============================================================================*
*  QUESTION 7 — Random effects                                                 *
*==============================================================================*

/*  The RE estimator corresponding to Q3 (no year dummies) so that the
    FE-vs-RE comparison in Q8 is clean.  Under RE, the state effects are
    assumed uncorrelated with all regressors.
*/

di _newline(2) as text "=============================================="
di as text "  QUESTION 7: RANDOM EFFECTS"
di as text "=============================================="

xtreg fatalityrate `xvars' if sample_main, re vce(cluster fips)
estimates store re_state

scalar b_re = _b[sb_useage]
scalar se_re = _se[sb_useage]


*==============================================================================*
*  QUESTION 8 — Robust Hausman test via Mundlak / CRE approach                 *
*==============================================================================*

/*  Because we use cluster-robust SEs, the classical Hausman statistic
    is unreliable.  Instead we follow the Mundlak (1978) device: augment
    the RE model with the within-state means of the time-varying regressors
    and test their joint significance.  The RE model here matches Q7
    (no year dummies) to keep the FE-vs-RE comparison internally consistent.
*/

di _newline(2) as text "=============================================="
di as text "  QUESTION 8: ROBUST HAUSMAN / MUNDLAK"
di as text "=============================================="

preserve
keep if sample_main
xtset fips year

foreach var in sb_useage speed65 speed70 ba08 drinkage21 ln_income age {
    bysort fips: egen mean_`var' = mean(`var')
}

xtreg fatalityrate `xvars' ///
    mean_sb_useage mean_speed65 mean_speed70 mean_ba08 ///
    mean_drinkage21 mean_ln_income mean_age, ///
    re vce(cluster fips)
estimates store cre

test mean_sb_useage mean_speed65 mean_speed70 mean_ba08 ///
    mean_drinkage21 mean_ln_income mean_age
scalar chi2_mundlak = r(chi2)
scalar p_mundlak    = r(p)

di as result "Mundlak joint chi2 = " %9.4f chi2_mundlak
di as result "Mundlak p-value    = " %9.6f p_mundlak
restore


*==============================================================================*
*  QUESTION 9 — Enforcement laws and seat belt usage                           *
*==============================================================================*

/*  Dependent variable is now sb_useage.  We include primary, secondary,
    the same controls, and state + year FE.
*/

di _newline(2) as text "=============================================="
di as text "  QUESTION 9: ENFORCEMENT AND USAGE"
di as text "=============================================="

xtreg sb_useage primary secondary speed65 speed70 ba08 drinkage21 ///
    ln_income age i.year if sample_main, fe vce(cluster fips)
estimates store usage_tw

scalar b_prim = _b[primary]
scalar b_sec  = _b[secondary]

lincom primary - secondary
scalar delta_prim_sec = b_prim - b_sec

di as result "Primary coefficient         = " %9.6f b_prim
di as result "Secondary coefficient       = " %9.6f b_sec
di as result "Primary − Secondary         = " %9.6f delta_prim_sec


*==============================================================================*
*  QUESTION 10 — New Jersey: lives saved by switching to primary               *
*==============================================================================*

/*  NJ switched from secondary to primary enforcement in 2000.  We use
    the 1997 NJ traffic volume as the nearest available proxy.
*/

di _newline(2) as text "=============================================="
di as text "  QUESTION 10: NEW JERSEY"
di as text "=============================================="

qui summarize vmt if state == "NJ" & year == 1997
scalar nj_vmt = r(mean)

scalar delta_rate_nj = b_tw * delta_prim_sec
scalar lives_nj      = -delta_rate_nj * nj_vmt

di as result "NJ 1997 VMT (millions)      = " %9.1f nj_vmt
di as result "Change in NJ sb_useage      = " %9.6f delta_prim_sec
di as result "Change in NJ fatality rate  = " %9.6f delta_rate_nj
di as result "Lives saved/year in NJ      = " %9.1f lives_nj


*==============================================================================*
*  COMPACT SUMMARY                                                             *
*==============================================================================*

di _newline(2) as text "=============================================="
di as text "  SUMMARY: ALL FATALITY MODELS"
di as text "=============================================="

estimates table pooled fe_state fe_tw re_state, ///
    keep(`xvars') b(%9.4f) se(%9.4f) stats(N)

di _newline as text "=============================================="
di as text "  KEY SCALARS"
di as text "=============================================="
di as result "Pooled  sb_useage  = " %9.6f b_pool
di as result "FE      sb_useage  = " %9.6f b_fe
di as result "TW-FE   sb_useage  = " %9.6f b_tw
di as result "RE      sb_useage  = " %9.6f b_re
di as result "Year FE joint p    = " %9.6f p_yearFE
di as result "Mundlak chi2       = " %9.4f chi2_mundlak
di as result "Mundlak p          = " %9.6f p_mundlak
di as result "Primary coeff      = " %9.6f b_prim
di as result "Secondary coeff    = " %9.6f b_sec
di as result "NJ lives/year      = " %9.1f lives_nj

log close
