# Panel Data Assignment — M1 (2025–2026)

Homework for the **Panel Data Econometrics** course (M1 TSE, 2025–2026).  
The assignment covers pooled OLS, fixed effects, random effects, and robust Hausman testing applied to U.S. traffic-fatality data.

## Repository structure

```
├── README.md
├── M1_Panel_Homework_2025_2026.pdf   # Assignment instructions
├── SeatBelts_Description.pdf         # Variable definitions and data documentation
└── Final/
    ├── SeatBelts.csv                 # Panel dataset (51 states × 15 years)
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

1. Open `Final/panel_homework_final.do` in Stata.
2. Change the single `global projdir` path at the top to point to `Final/` on your machine.
3. Run the entire file. The log will be written to the same folder.

Everything needed (data + code) lives inside `Final/`, so no relative-path issues.

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
