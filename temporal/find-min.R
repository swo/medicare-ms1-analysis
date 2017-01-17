min_max = function(x) {
  m = loess(n_rx ~ week, x)
  f = function(x) predict(m, x)
  min_res = optimize(f, c(1, 53))
  max_res = optimize(f, c(1, 53), maximum=TRUE)
  data_frame(min_week=c(min_res$minimum), min_obj=c(min_res$objective), max_week=c(max_res$maximum), max_obj=c(max_res$objective))
}

mm = dat %>%
  ungroup %>%
  filter(abx=='azithromycin') %>%
  select(week, state, n_rx) %>%
  group_by(state, week) %>%
  summarize(n_rx=sum(n_rx)) %>%
  do(min_max(.))