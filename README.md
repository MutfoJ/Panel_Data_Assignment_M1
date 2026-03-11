# Panel Data Assignment — M1 (2025–2026)

Homework for the **Panel Data Econometrics** course (M1 TSE, 2025–2026).  
The assignment covers pooled OLS, fixed effects, random effects, and robust Hausman testing applied to U.S. traffic-fatality data.

## Repository structure

```
├── README.md
├── M1_Panel_Homework_2025_2026.pdf   # Assignment instructions
├── SeatBelts_Description.pdf         # Variable definitions and data documentation
├── SeatBelts.csv                     # Panel dataset (51 states × 15 years)
└── Final/
    ├── panel_homework_final.do       # Stata do-file (replicates all results)
    ├── panel_homework_final.log      # Stata log (full output)
    ├── panel_homework_final.md       # Answer sheet (Markdown source)
    └── panel_homework_final.pdf      # Answer sheet (compiled PDF)
```

## Data

**SeatBelts.csv** — semicolon-delimited panel of 765 observations (51 U.S. states including D.C., 1983–1997).  
Key variables include `fatalityrate`, `sb_useage`, `primary`, `secondary`, enforcement dummies, and socio-demographic controls.  
See `SeatBelts_Description.pdf` for full variable definitions.

## Replication

1. Place `SeatBelts.csv` one directory above `Final/` (or adjust the path in the do-file).
2. Open Stata and run:
   ```stata
   do "Final/panel_homework_final.do"
   ```
3. The log file `panel_homework_final.log` will be produced inside `Final/`.

## Methods

| Question | Model | Key result |
|----------|-------|------------|
| Q2 | Pooled OLS | sb_useage = +0.0041 (biased, wrong sign) |
| Q3 | State FE | sb_useage = −0.0058 (sign reversal) |
| Q4 | Two-way FE | sb_useage = −0.0037 (preferred spec) |
| Q7 | Random Effects | sb_useage = −0.0045 |
| Q8 | Mundlak test | χ²(7) = 28.42, reject RE |
| Q9 | Usage equation | Primary +20.6 pp, Secondary +10.9 pp |

## Authors

Dragos Florin Vasile & Wong Hei Wong
