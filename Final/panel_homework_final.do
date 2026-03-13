*==============================================================================
*  Panel Data Homework 2025–2026
*  Authors: Dragos Florin Vasile & Wong Hei Wong
*  Date: March 2026
*==============================================================================

/*  INSTRUCTIONS
    ============
    This file replicates all empirical results in the answers sheet. It is organized by assignment question. All output goes to a generated log file so the reader can cross-check numbers easily.

    HOW TO RUN
    ----------
     1. Change the path in the line below to point to the directory that
         contains SeatBelts.csv and this .do file.
    2. Run the entire file. A log file will be generated in the same folder.
*/

clear all
set more off
capture log close _all

* >>>>>>>  CHANGE THIS PATH to the folder where SeatBelts.csv is located  <<<<<<
global projdir "C:\Users\vdrag\OneDrive\Escritorio\Panel Data Project\Final"
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

cd "$projdir"
log using "$projdir/panel_homework_final.log", replace text

* ---------- Path to data (same folder) ---------------------------------------
local datafile "$projdir/SeatBelts.csv"
capture confirm file "`datafile'"
if _rc {
    di as error "SeatBelts.csv not found in $projdir — place it next to this .do."
    exit 601
}

* A local that stores all time-varying controls used in fatality regressions.
local xvars sb_useage speed65 speed70 ba08 drinkage21 ln_income age


*==============================================================================
*  QUESTION 1 — Data description and cleaning
*==============================================================================

/*  We import the semicolon-delimited CSV, confirm there are no duplicates,
    declare the panel structure, generate log-income, and define the main
    estimation sample. The only issue is that sb_useage is missing
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
label var sample_main "Complete-case sample"

*Summary for the raw data set

sum

* Counting Raw observations
count

* Counting Cleaned observations

count if sample_main

* Counting missing variables in our observations

count if missing(sb_useage)

* Missingness pattern and summary statistics
misstable summarize, all
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
scalar mean_us_vmt = r(mean)
di mean_us_vmt
qui summarize us_vmt if year == 1997
scalar us_vmt_1997 = r(mean)
restore


*==============================================================================
*  QUESTION 2 — Pooled OLS
*==============================================================================

/*  Baseline cross-sectional regression.  We expect the coefficient on
    sb_useage to be misleading because of omitted state-level confounders.
*/

reg fatalityrate `xvars' if sample_main, vce(cluster fips)
estimates store pooled

scalar b_pool = _b[sb_useage]
scalar se_pool = _se[sb_useage]


*==============================================================================
*  QUESTION 3 — State FE
*==============================================================================

/*  Adding state FE removes all time-invariant unobserved heterogeneity.
    Identification now comes from within-state variation over time.
*/

xtreg fatalityrate `xvars' if sample_main, fe vce(cluster fips)
estimates store fe_state

scalar b_fe = _b[sb_useage]
scalar se_fe = _se[sb_useage]


*==============================================================================
*  QUESTION 4 — State + year FE
*==============================================================================

/*  Year dummies absorb national trends in fatality rates (safer cars,
    federal campaigns, etc.).  We test their joint significance.
*/

xtreg fatalityrate `xvars' i.year if sample_main, fe vce(cluster fips)
estimates store fe_tw

scalar b_tw = _b[sb_useage]
scalar se_tw = _se[sb_useage]

testparm i.year
scalar p_yearFE = r(p)
di p_yearFE


*==============================================================================
*  QUESTION 6 — Effect size and lives saved
*==============================================================================

/*  We translate the two-way FE estimate into an effect size: how many
    lives would be saved if seat belt use rose from 52 % to 90 %?
*/

*Reduction in death rate
scalar delta_sb = 0.90 - 0.52
scalar delta_rate = b_tw * delta_sb

*Reduction in fatality rate
qui summarize fatalityrate if sample_main
scalar mean_fatal = r(mean)
scalar pct_mean  = -100 * delta_rate / mean_fatal

*Reduction in fatalities figure
scalar lives_us_mean = -delta_rate * mean_us_vmt
scalar lives_us_1997 = -delta_rate * us_vmt_1997

di b_tw
di delta_rate
di pct_mean
di lives_us_mean


*==============================================================================
*  QUESTION 7 — Random effects
*==============================================================================

/*  The RE estimator corresponding to Q3 (no year dummies) so that the
    FE-vs-RE comparison in Q8 is clean.  Under RE, the state effects are
    assumed uncorrelated with all regressors.
*/

xtreg fatalityrate `xvars' if sample_main, re vce(cluster fips)
estimates store re_state

scalar b_re = _b[sb_useage]
scalar se_re = _se[sb_useage]


*==============================================================================
*  QUESTION 8 — Robust Hausman test (xtoverid)
*==============================================================================

/*  Relies on the xtoverid package. If missing, try auto-install from SSC.
*/

* Check required package
capture which xtoverid
if _rc {
    di as text "xtoverid not found. Attempting auto-install from SSC..."
    capture noisily ssc install xtoverid, replace
    capture which xtoverid
    if _rc {
        di as error "Auto-install failed. Please run manually: ssc install xtoverid"
        exit 199
    }
}

* FE Model
xtreg fatalityrate `xvars' if sample_main, fe vce(cluster fips)
estimates store fe_state

* RE Model
xtreg fatalityrate `xvars' if sample_main, re vce(cluster fips)
estimates store re_state

* Robust Hausman Test
quietly xtreg fatalityrate `xvars' if sample_main, re vce(cluster fips)
xtoverid

* Store estimates returned by xtoverid (varies by version)
scalar chi2_robust = .
scalar p_robust    = .
capture scalar chi2_robust = r(j)
capture scalar p_robust    = r(jp)


*==============================================================================
*  QUESTION 9 — Enforcement laws and seat belt usage
*==============================================================================

/*  Dependent variable is now sb_useage.  We include primary, secondary,
    the same controls, and state + year FE.
*/

xtreg sb_useage primary secondary speed65 speed70 ba08 drinkage21 ///
    ln_income age i.year if sample_main, fe vce(cluster fips)
estimates store usage_tw

scalar b_prim = _b[primary]
scalar b_sec  = _b[secondary]

lincom primary - secondary
scalar delta_prim_sec = b_prim - b_sec

di b_prim
di b_sec
di delta_prim_sec


*==============================================================================
*  QUESTION 10 — New Jersey lives saved
*==============================================================================

/*  NJ switched from secondary to primary enforcement in 2000.  We use
    the 1997 NJ traffic volume as the nearest available proxy.
*/

qui summarize vmt if state == "NJ" & year == 1997
scalar nj_vmt = r(mean)

scalar delta_rate_nj = b_tw * delta_prim_sec
scalar lives_nj      = -delta_rate_nj * nj_vmt

qui summarize fatalityrate if fips == 34
scalar nj_mean_deaths = r(mean)
scalar save_rate_nj = -100 * delta_rate_nj / nj_mean_deaths

di nj_vmt
di delta_prim_sec
di delta_rate_nj
di lives_nj
di save_rate_nj



*==============================================================================
* Summary of Key Results
*==============================================================================

estimates table pooled fe_state fe_tw re_state, ///
    keep(`xvars') b(%9.4f) se(%9.4f) stats(N)

* Final scalar recap (custom display)
di _newline as text "=============================================="
di as text "  KEY SCALARS"
di as text "=============================================="
di as result "Pooled  sb_useage  = " %9.6f b_pool
di as result "FE      sb_useage  = " %9.6f b_fe
di as result "TW-FE   sb_useage  = " %9.6f b_tw
di as result "RE      sb_useage  = " %9.6f b_re
di as result "Year FE joint p    = " %9.6f p_yearFE
di as result "Sargan-Hansen chi2 = " %9.4f chi2_robust
di as result "Sargan-Hansen p    = " %9.6f p_robust
di as result "Primary coeff      = " %9.6f b_prim
di as result "Secondary coeff    = " %9.6f b_sec
di as result "Primary-Secondary  = " %9.6f delta_prim_sec
di as result "Lives US (avg VMT) = " %9.1f lives_us_mean
di as result "Lives US (1997)    = " %9.1f lives_us_1997
di as result "NJ lives/year      = " %9.1f lives_nj

log close
