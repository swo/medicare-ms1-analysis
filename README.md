% Scott Olesen
% 5 April 2017

# Subfolder/analyses

- `bayesian`: quasi-Bayesian hierarchical model of consumption-resistance
- `carrier-outpatient`
- `chronic-conditions`: Correlations between conditions, relationship to consumption. Probably deprecated in favor of `ms1`
- `death-year`: Do people have really different consumption in their last year of life? Not really.
- `dosage-form`: Analysis to help me figure out which dosage codes to use.
- `drug-drug-correlations`: Original tetrachoric correlations between drugs. Matthews correlation coefficient is probably better, or some sort of Poisson thing.
- `fit-distribution`: Attempts to fit the consumption histograms to know curves. Poisson lognormal is good; negative binomial is OK; both are better than Poisson. But I found in ms1 that Poisson and negative binomial regressions give the same results, so this probably doesn't matter. Deprecated.
- `hicks-summary`: Summary of consumption by drug and beneficiary characteristics. Rolled into ms1. Probably deprecated.
- `inequality`: Measurements of consumption unevenness by state and drug. Plots of unevenness against population. Should be rolled into ms2.
- `inequality-resistance`: Consumption unevenness against resistance. Should be rolled into ms2.
- `ms1`: Analyses in manuscript 1.
- `region-consumption`: Early visualizations of unevenness by region. Better answered in `inequality`. Deprecated.
- `sex-difference`: Will become ms3.
- `states`: Consumption and resistance by state and HRR. Should be rolled into ms2.
- `temporal`: Temporal patterns in consumption. Maybe rolled into GISP ms.
- `top-users`: Consumption patterns of top users. They mostly take the same drugs as everyone else, with slight variations. Unclear where these results will go.
- `verification`: Checks on the data, which helped me clear up the 2012 weirdness.
