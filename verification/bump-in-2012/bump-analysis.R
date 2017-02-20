summarize_year = function(year) {
  fn = sprintf('../bene-abx-%s.tsv', year)
  read_tsv(fn) %>%
    mutate(state=factor(state)) %>%
    select(bene, days, ciprofloxacin, state) %>%
    rename(abx_days=days, azi_days=ciprofloxacin)
}

dat11 = summarize_year('2011') %>% select(bene, days11=abx_days, state)
dat12 = summarize_year('2012') %>% select(bene, days12=abx_days)
dat13 = summarize_year('2013') %>% select(bene, days13=abx_days)

x = dat11 %>%
  filter(state=='Florida') %>%
  left_join(dat12, by='bene', suffix=c('.11', '.12')) %>%
  left_join(dat13, by='bene')