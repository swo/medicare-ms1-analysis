# BENE_ID	SRVC_DT	PRDSRVID	DAWPS_CD	QTYDSPNS	DAYSSPLY	FILL_NUM	DRCVSTCD	TOTALCST	FORMULARY_ID	FRMLRY_RX_ID	TIER_ID	STEP	QUANTITY_LIMIT_YN	PRIOR_AUTHORIZATION_YN	BN	GNN	GCDF	GCDF_DESC
# zzzzzzzwPSfeewe	17MAR2011	00591081055	0	50.000	25	0	C	8.02	xxxkk888	00002001	01		0	0	SILVER SULFADIAZINE	SILVER SULFADIAZINE	1 %	KA	CREAM (GRAM)
# zzzzzzzmzmfmeSX	02AUG2011	00065000603	0	3.000	10	0	C	84.13	xxxxxxxk	00008689	01		0	0	MOXEZA	MOXIFLOXACIN HCL	0.5 %	RH	DROPS, VISCOUS (ML)

oral_injected_codes = scan('dosage-form/oral_injected_codes.txt', what=character())
antibiotic_names = read_tsv('abx_class.tsv')

polish_pde = function(x) {
  x %>%
    select(bene_id=BENE_ID,
           service_date=SRVC_DT,
           days_supply=DAYSSPLY,
           generic_name=GNN,
           dosage_form_code=GCDF,
           dosage_form_desc=GCDF_DESC) %>%
    mutate(service_date=dmy(service_date),
           days_supply=as.integer(days_supply)) %>%
    filter(dosage_form_code %in% oral_injected_codes) %>%
    left_join(antibiotic_names, by='generic_name') %>%
    select(bene_id, service_date, antibiotic, days_supply)
}

polish_year = function(year, n_max=Inf) {
  in_fn = sprintf('../data/pde_%i.tsv', year)
  out_fn = sprintf('pde_%i.tsv', year)
  read_tsv(in_fn, n_max=n_max) %>% polish_pde %>% write_tsv(out_fn)
}

for (year in 2011:2014) polish_year(year)