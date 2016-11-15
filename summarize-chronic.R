#!/usr/bin/env Rscript

# Summarize the chronic conditions data into two-column table

# Digits in the chronic condition data have a funny encoding:
# "3" means that the beneficiary has the condition; other numbers
# mean other things. Convert "3" to 1 and everything else to 0.
change_digit = function(x) ifelse(x == 3, 1, 0)

summarize_cc = function(in_fn, out_fn) {
  read_tsv(in_fn) %>%
    select(BENE_ID,
           AMI, ALZH, COPD, CHF, DIABETES, ISCHMCHT, STRKETIA,
           CHRNKIDN, CNCRBRST, CNCRCLRC, CNCRPRST, CNCRLUNG, CNCRENDM) %>%
    rename(bene=BENE_ID) %>%
    mutate_at(vars(-bene), change_digit) %>%
    mutate(comorbidity=AMI+ALZH+COPD+CHF+DIABETES+ISCHMCHT+STRKETIA+
           2*(CHRNKIDN+CNCRBRST+CNCRCLRC+CNCRPRST+CNCRLUNG+CNCRENDM)) %>%
    mutate(comorbidity=as.integer(comorbidity)) %>%
    select(bene, comorbidity) %>%
    write_tsv(out_fn)
}
 
summarize_cc("../data/cc2011.tsv", "cc2011.tsv")