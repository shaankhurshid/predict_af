# predict_af
Data processing scripts for Predict-AF analysis

## Dependencies
R v4.0 (recommended)
Libraries: data.table, survival, rms, prodlim, timeROC, nricens, epiR, dca
devtools::install_github("ddsjoberg/dcurves")

## Data/Variable Dictionary
- Analysis assumes presence of two data.tables
  -   af_set: all eligible individuals (i.e., without prevalent AF)
  -   stroke_set: subset of af_set also without prevalent stroke
- Analysis assumes presence of a data.table with the following columns defined
  - Outcome variables
    - incd_af_5y (binary): Incident AF status at 5 years (1 = yes, 0 = no)
    - incd_af_5y.t (continuous): AF follow-up time in years (ends at AF or censored at earliest of death, last follow-up, or administrative censoring as applicable)
    - stroke_5y (binary): Incident stroke status at 5 years (1 = yes, 0 = no)
    - stroke_5y.t (continuous): Stroke follow-up time in years (ends at stroke or censored at earliest of death, last follow-up, or administrative censoring as applicable)
  - Exposure variables
    - charge (continuous): Individual's CHARGE-AF score (see PMID: 23537808 for definition)
    - prs_norm (continuous): Summed PRS divided by number of alleles available in dataset {cite paper when published}
    - agevisit0 (continuous): Age at start of follow-up
    - needs_oac_visit0 (binary): Whether the individuals CHADS2-VASc score is >= 2 (male) or >= 3 (female) at start of follow-up
    - needs_oac_at_af (binary): Whether the individuals CHADS2-VASc score is >= 2 (male) or >= 3 (female) at the time of an AF diagnosis

