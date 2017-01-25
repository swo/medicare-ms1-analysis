#!/usr/bin/env Rscript

# create summary beneficiary x abx file

hrr = read_tsv("../db/hrr/ZipHsaHrr11.txt") %>%
  select(zip=zipcode11, hrr=hrrnum)
census = read_tsv("../db/census-regions/census-regions.tsv") %>%
    select(state, region)

report = function(dat, group, out_fn) {
  dat %>%
    group_by_(group) %>%
    summarize(n=n()) %>%
    write_tsv(out_fn)
}

report_for_year = function(year) {
  abx_fn = sprintf('bene-abx-%s.tsv', year)

  abx = read_tsv(abx_fn) %>%
    select(bene, zip, state) %>%
    left_join(hrr, by='zip') %>%
    left_join(census, by='state')

  for (g in c('region', 'state', 'hrr')) {
    out_fn = sprintf('summary-bene-location-%s-%s.tsv', g, year)
    report(abx, g, out_fn)
  }
}

for (y in c('2011', '2012', '2013', '2014')) {
  report_for_year(y)
}
