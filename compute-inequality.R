# create summary beneficiary x abx file

codes = read_tsv('../db/proc/code-map.txt', col_names=c('ndc', 'abx'))

all_pde_counts = read_tsv('../data/pde_freq_2011.txt') %>%
  rename(bene=BENE_ID, all_pde=COUNT)

load.abxpde = function(fn) {
  read_tsv(fn) %>%
    rename(bene=BENE_ID, ndc=PRDSRVID, days=DAYSSPLY) %>%
    mutate(days=as.integer(days)) %>%
    inner_join(codes, by='ndc') %>%
    select(bene, abx, days)
}

county_codes = read_tsv("../db/state-county-codes/county-codes.tsv")
race_codes = read_tsv("../db/race-codes/race.tsv")

load.bene = function(fn, n_max=Inf) {
  # get and rearrange raw data
  read_tsv(fn, n_max=n_max) %>%
    # join state and county codes to merge with SSA list
    unite(state_county_code, STATE_CD, CNTY_CD, sep='') %>%
    # keep only the 5-digit zip
    mutate(zip=substr(BENE_ZIP, 0, 5)) %>%
    select(bene=BENE_ID, age=AGE, state_county_code, zip,
           race_code=RACE,
           plan_coverage_months=PLNCOVMO) %>%
    # merge (and drop) state+county code
    left_join(county_codes, by='state_county_code') %>%
    select(-state_county_code) %>%
    # merge (and drop) race codes
    left_join(race_codes, by='race_code') %>%
    select(-race_code)
}

bene = load.bene('../data/bene_match_2011.txt')
abxpde = load.abxpde('../data/abx_pde_2011.tsv')

# convert to a wide abx pde data frame 
proto_wap = abxpde %>%
  group_by(bene, abx) %>%
  summarize(n_claims=n(), days=sum(days))

wap = proto_wap %>%
  ungroup %>%
  spread(abx, days, fill=0)

wap$total_abx = wap %>% select(-bene, -n_claims) %>% rowSums

# merge in the information from all the pdes
wap %<>% inner_join(all_pde_counts, by='bene')

# merge in the information about the beneficiaries
wap %<>% inner_join(bene, by='bene')

# save this new file
write_tsv(wap, 'bene-abx-2011.tsv')