#!/usr/bin/env Rscript

# Polish the beneficiary and antibiotic PDE files at the same time:
#  - Keep all the beneficiaries that match criteria (good geography, 12 months coverage)
#  - Keep PDEs that match known beneficiaries

state_codes = read_tsv('../db/state-codes/state-codes.tsv') %>%
  mutate(is_state_or_dc=as.logical(is_state_or_dc)) %>%
  filter(is_state_or_dc)

oral_injected_codes = scan('dosage-form/oral_injected_codes.txt', what=character())
antibiotic_names = read_tsv('abx_class.tsv')

polish_bene = function(x) {
  rename(x, bene_id=BENE_ID,
          state_cd=STATE_CD,
          age=AGE,
          hmo_months=HMO_MO,
          plan_coverage_months=PLNCOVMO) %>%
    mutate_at(vars(age, hmo_months, plan_coverage_months), as.integer) %>% # should be already done from read_tsv
    left_join(state_codes, by='state_cd') %>%
    filter(!is.na(state), plan_coverage_months==12) %>% # filter for geography and plan
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
    select(bene_id, sex, age, race, state, zipcode, hmo_months)
}

polish_pde = function(x) {
  select(x, bene_id=BENE_ID,
            service_date=SRVC_DT,
            fill_number=FILL_NUM,
            days_supply=DAYSSPLY,
            generic_name=GNN,
            dosage_form_code=GCDF,
            dosage_form_desc=GCDF_DESC) %>%
    mutate(service_date=dmy(service_date),
           days_supply=as.integer(days_supply)) %>%
    filter(dosage_form_code %in% oral_injected_codes) %>%
    left_join(antibiotic_names, by='generic_name') %>%
    select(bene_id, service_date, fill_number, antibiotic, antibiotic_class, days_supply)
}

# bene column types
# c BENE_ID zzzzzzLmmLLLPwP zzzzzzLmmLLLSPL
# i RFRNC_YR  2011  2011
# c STATE_CD
# c BENE_ZIP  999999999 999999999
# i AGE 74  39
# c BENE_DOB  08NOV1937 22NOV1972
# c DEATH_DT
# i SEX 1 2
# i RACE  3 1
# i RTI_RACE_CD 4 1
# c MS_CD 10  20
# i HMO_MO  0 0
# c STRCTFLG  05  15
# i PLNCOVMO

bene_column_types = 'ciccicciiicici'

# pde column types
# c BENE_ID zzzzzzzeLXPSzXf zzzzzzzeLXPSzXf
# c SRVC_DT 13MAR2011 29NOV2011
# c PRDSRVID  45802005635 45802005635
# c DAWPS_CD  0 0
# n QTYDSPNS  30.000  30.000
# i DAYSSPLY  30  30
# i FILL_NUM  0 1
# c DRCVSTCD  C C
# n TOTALCST  7.00  7.00
# c FORMULARY_ID  xxxkk888  xxxkk888
# c FRMLRY_RX_ID  00002013  00002013
# c TIER_ID 01  01
# c STEP
# i QUANTITY_LIMIT_YN 0 0
# i PRIOR_AUTHORIZATION_YN  0 0
# c BN  GENTAMICIN SULFATE  GENTAMICIN SULFATE
# c GNN GENTAMICIN SULFATE  GENTAMICIN SULFATE
# c STR 0.1 % 0.1 %
# c GCDF  KA  KA
# c GCDF_DESC CREAM (GRAM)  CREAM (GRAM)

pde_column_types = 'ccccniicncccciiccccc'

parse_year = function(year, n_max=Inf) {
  bene_in_fn = sprintf('../data/bene_%i.tsv', year)
  pde_in_fn = sprintf('../data/pde_%i.tsv', year)

  bene_out_fn = sprintf('bene_%i.tsv', year)
  pde_out_fn = sprintf('pde_%i.tsv', year)

  bene = read_tsv(bene_in_fn, n_max=n_max, col_types=bene_column_types) %>%
    polish_bene

  pde = read_tsv(pde_in_fn, n_max=n_max, col_types=pde_column_types) %>%
    polish_pde %>%
    semi_join(select(bene, bene_id), by='bene_id') # keep only matching entries

  write_tsv(bene, bene_out_fn)
  write_tsv(pde, pde_out_fn)
}

for (year in 2011:2014) parse_year(year)
