library(hash)
library(readr)

code_map = read_tsv('../db/proc/code-map.txt', col_names=c('ndc', 'abx'))

dat = read_tsv('../data/abx_pde_2011.txt', col_types='ccccnnn') %>%
  rename(bene_id=BENE_ID, service_dat=SRVC_DT, ndc=PRDSRVID, dawps=DAWPS_CD, quantity=QTYDSPNS, days_supply=DAYSSPLY, cost=TOTALCST) %>%
  select(bene_id, ndc)

dat %<>% left_join(code_map, by='ndc')

dat_sum = dat %>%
  select(bene_id, abx) %>%
  group_by(bene_id) %>%
  summarize(counts=n()) %>%
  arrange(desc(counts))

write_tsv(dat_sum, 'bene-sum-2011.tsv')
