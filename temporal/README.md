% Temporal Medicare antibiotic data
% Scott Olesen (olesen@hsph.harvard.edu)
% 17 Jan 2017

# Overview

These data are drawn from Medicare Part D beneficiaries. The original data
include information about the beneficiaries (e.g., state of residence) and
information about the prescriptions (e.g., date, quantity, and identity of the
drug). The drug's identity is listed using its 11-digit NDC code.  I had to
develop a translation between NDC codes and the names of the antibiotics.

The data cover 2011 to 2014. The data are a 20% sample of Part D beneficiaries
(i.e., all claims from 20% of beneficiaries). Beneficiaries enter and leave the
Part D program, but the 20% sample aims to remain the same. 70% of
beneficiaries appear in all four years. Beneficiaries might have claims that
are not represented among these data. These data are also not necessarily an
unbiased sample of Medicare beneficiaries: the more wealthy beneficiaries tend
to have alternatives to Part D.

# Data files

There are four data files, one for each year. The consumption of antibiotics in
2012 appears to be about 20% larger than in the other years; I don't have an
explanation for that.

Each file has tab-separated data fields:

- `week`: The week of the year in which the prescriptions were filled. (They are integers 1 through 53; they appear as decimals only for technical reasons.) December 31 is called week 53.
- `state`
- `abx`: Name of the drug
- `n_rx`: Number of prescriptions filled during that week
- `n_days`: Sum of the "days' supply" value supplied in the Part D data
