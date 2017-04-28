qplm = function(df, indep, expl) {
  f = formula(paste0(indep, ' ~ ', paste0(expl, collapse=' + ')))
  glm(f, data=df, family=quasipoisson(link='identity'), start=rep(0.1, length(labels(terms(f))) + 1))
}

# Quasi-Poisson (w/ identity link) partial correlation
qppc = function(df, y, x, controls=NULL) {
  stopifnot(is.character(y) && is.character(x))
  stopifnot(is.null(controls) || is.vector(controls, mode='character'))
  if (is.null(controls) || length(controls) == 0) {
    # this is just a normal correlation
    cor.test(df[[x]], df[[y]]) %>% tidy
  } else {
    ry = qplm(df, y, controls)$residuals
    rx = qplm(df, x, controls)$residuals
    cor.test(rx, ry) %>% tidy
  }
}

qppr = function(df, y, x, controls) {
  stopifnot(is.character(y) && is.character(x))
  stopifnot(is.vector(controls, mode='character'))
  
  ry = qplm(df, y, controls)$residuals
  rx = qplm(df, x, controls)$residuals
  lm(ry ~ rx) %>% tidy
}

bene = read_tsv('../bene_2011.tsv') %>%
  filter(between(age, 66, 96), hmo_months==0) %>%
  mutate(is_female=sex=='female', is_white=race=='white') %>%
  select(-zipcode, -hmo_months, -sex, -race)

cu = read_tsv('../cu_2011.tsv') %>%
  filter(bene_id %in% bene$bene_id)

cc = read_feather('../cc_2011.feather') %>%
  filter(bene_id %in% bene$bene_id)

pde = read_tsv('../pde_2011.tsv') %>%
  filter(bene_id %in% bene$bene_id) %>%
  count(bene_id) %>%
  rename(abx_pde=n)

dat = bene %>%
  left_join(pde, by='bene_id') %>%
  left_join(cu, by='bene_id') %>%
  left_join(cc, by='bene_id') %>%
  replace_na(list(abx_pde=0))

ncc = dat %>% select(AMI:HYPOTH) %>% rowSums
dat %<>% mutate(ncc=ncc)

f = function(df, group_name) {
  group_by(df, age) %>%
    summarize(mean_abx_pde=mean(abx_pde), sem=sd(abx_pde)/sqrt(n())) %>%
    mutate(group=group_name)
}
x = bind_rows(
  f(dat, 'all'),
  f(dat %>% filter(ncc==0), 'no_cc'),
  f(dat %>% filter(hospital_op_visits==0), 'no_op')
) %>% 
  filter(age >= 67)

p = ggplot(x, aes(x=age, y=mean_abx_pde, group=group)) +
  geom_point() + geom_line(size=0.5) +
  geom_linerange(aes(ymin=mean_abx_pde-sem, ymax=mean_abx_pde+sem)) +
  scale_x_continuous(limits=c(67, 96), breaks=c(67, 70, 80, 90, 96)) +
  ylim(0.5, 2.0) +
  ylab('mean num. of antibiotic claims') +
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_line(color='#f2f2f2'),
        text=element_text(size=14)) +
  annotate('text',
           x=c(75, 75, 85),
           y=c(1.65, 0.95, 0.60),
           label=c('All beneficiaries', 'No chronic conditions', 'No hospital OP visits'),
           size=5)

show(p)