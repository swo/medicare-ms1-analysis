library(purrr)
library(lazyeval)

drugs = c('ciprofloxacin', 'levofloxacin', 'trimethoprim-sulfamethoxazole', 'ampicillin', 'ampicillin-sulbactam')
bugs = c('Escherichia coli')

hrr = read_tsv('../../db/hrr/ZipHsaHrr11.txt') %>%
  select(zip=zipcode11, hrr=hrrnum)

abx = read_tsv('../bene-abx-2011.tsv') %>%
  # swo>
  # filter out the people who are not old
  filter(age >= 65) %>%
  left_join(hrr, by='zip')

abg = read_tsv('../../../antibiogram/data/abg.tsv', col_types=cols(zipcode=col_character())) %>%
  rename(zip=zipcode) %>%
  mutate(drug=tolower(drug))

# compute susceptibility information for each bug/drug combo
sus = abg %>%
  filter(bug %in% bugs, drug %in% drugs) %>%
  left_join(hrr, by='zip') %>%
  group_by(bug, drug, hrr) %>%
  summarize(n_abg=n(),
            n_isolates=sum(n_isolates, na.rm=TRUE),
            n_years=length(unique(year)),
            median_sus=median(percent_susceptible, na.rm=TRUE),
            mean_sus=mean(percent_susceptible, na.rm=TRUE),
            sd_sus=sd(percent_susceptible, na.rm=TRUE),
            min_sus=min(percent_susceptible, na.rm=TRUE),
            max_sus=max(percent_susceptible, na.rm=TRUE))

write_tsv(sus, 'tmp-sus')

known_hrrs = sus %$% hrr %>% unique %>% sort

# compute summary stats for the hrr
hrr_summary = abx %>%
  filter(hrr %in% known_hrrs) %>%
  group_by(hrr) %>%
  summarize(n_bene=n(),
            n_abx_claims=sum(n_claims),
            n_abx_days=sum(days),
            n_zero_abx=sum(days == 0),
            n_pdes=sum(all_pde),
            mean_comorbidity=mean(comorbidity),
            mean_age=mean(age),
            n_white=sum(race=='white'),
            n_male=sum(sex=='male'))

write_tsv(hrr_summary, 'tmp-hrr')

nz = function(x) x[x > 0]

# compute consumption for each bug/drug combo
drug_vars = map(drugs, as.name)
con = abx %>%
  filter(hrr %in% known_hrrs) %>%
  select_(.dots=c('bene', 'hrr', drug_vars)) %>%
  gather_('drug', 'days', drugs) %>%
  group_by(hrr, drug) %>%
  summarize(n_bene=n(),
            n_days=sum(days),
            n_bene_zero=sum(days==0),
            max_days=max(days),
            mean_nz_days=mean(nz(days)),
            median_nz_days=median(nz(days)),
            quantile_nz_95=quantile(nz(days), probs=c(0.95)),
            quantile_nz_99=quantile(nz(days), probs=c(0.99))) %>%
  mutate(did=1000/365*n_days/n_bene)

write_tsv(con, 'tmp-con')

models = function(df) {
  out = data_frame()
  for (resp in c('median_sus', 'mean_sus', 'min_sus', 'max_sus')) {
    for (expl in c('n_days', 'did', 'n_bene_zero', 'max_days', 'mean_nz_days', 'median_nz_days', 'quantile_nz_95', 'quantile_nz_99')) {
      new_rows = df %>%
        mblm(as.formula(paste0(resp, '~', expl)), dataframe=.) %>%
        tidy %>%
        mutate(response=resp, explanatory=expl)
      out = bind_rows(out, new_rows)
    }
  }
  out
}

res = con %>%
  ungroup %>%
  left_join(sus, by=c('drug', 'hrr')) %>%
  filter(!is.na(n_isolates)) %>%
  arrange(drug) %>%
  group_by(bug, drug) %>%
  do(models(.))

write_tsv(res, 'tmp-res')
