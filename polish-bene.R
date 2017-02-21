state_codes = read_tsv('../db/state-codes/state-codes.tsv') %>%
  mutate(is_state_or_dc=as.logical(is_state_or_dc)) %>%
  filter(is_state_or_dc)

parse_bene = function(x) {
  x %>%
    rename(bene_id=BENE_ID,
           state_cd=STATE_CD,
           age=AGE,
           plan_coverage_months=PLNCOVMO) %>%
    left_join(state_codes, by='state_cd') %>%
    mutate(zipcode=substr(BENE_ZIP, 1, 5)) %>%
    mutate(sex=case_when(.$SEX == 0 ~ 'unknown',
                        .$SEX == 1 ~ 'male',
                        .$SEX == 2 ~ 'female')) %>%
    mutate(race=case_when(.$RTI_RACE_CD == 0 ~ 'unknown',
                          .$RTI_RACE_CD == 1 ~ 'white', # non-Hispanic white
                          .$RTI_RACE_CD == 2 ~ 'black',
                          .$RTI_RACE_CD == 3 ~ 'other',
                          .$RTI_RACE_CD == 4 ~ 'Asian', # or Pacific Islander
                          .$RTI_RACE_CD == 5 ~ 'Hispanic',
                          .$RTI_RACE_CD == 6 ~ 'American Indian')) %>% # or Alaska Native
    select(bene_id, sex, age, race, state, zipcode, plan_coverage_months)
}

bene_col_types = 'cicccicciiici'

parse_year = function(year, n_max=Inf) {
  in_fn = sprintf('../data/bene_%g.tsv', year)
  out_fn = sprintf('bene_%g.tsv', year)
  read_tsv(in_fn, n_max=n_max, col_types=bene_col_types) %>%
    parse_bene %>% write_tsv(out_fn)
}

for (year in 2011:2014) parse_year(year)