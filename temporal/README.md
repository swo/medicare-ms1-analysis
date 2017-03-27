% Temporal Medicare antibiotic data
% Scott Olesen (olesen@hsph.harvard.edu)
% 23 Mar 2017

# Overview

These data are drawn from Medicare Part D beneficiaries. The original data
include information about the beneficiaries (e.g., state of residence) and
information about the prescriptions (e.g., date, quantity, and identity of the
drug).

The data cover 2011 to 2014. The data are a 20% sample of Part D beneficiaries
(i.e., all claims from 20% of beneficiaries). These data only cover
beneficiaries who were covered under Medicare for the entire data year, and it
excludes those people who had an HMO (a.k.a., Medicare Advantage, which
effectively censors some of their prescriptions from this database).

# Data files

There are four data files, `temporal_abx_20XX.tsv`, one for each year. Each
file has tab-separated data fields:

- `week`: The week of the year in which the prescriptions were filled. (They are integers 1 through 53; they appear as decimals only for technical reasons.) December 31 is called week 53.
- `state`
- `antibiotic`: Name of the drug
- `n_claims`: Number of claims made during that week
- `n_days`: Sum of the "days' supply" value supplied in the Part D data

I also included `abx_class.tsv`, which maps the individual antibiotics into
larger groups.
