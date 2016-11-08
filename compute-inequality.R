# summarize unnevenness

library(ineq)

olesen_index = function(counts) {
  x = data.frame(counts=counts) %>%
    arrange(counts) %>%
    mutate(x=1:dim(.)[1], x=x/max(x), y=cumsum(counts)/sum(counts)) %>%
    filter(x + y - 1 < 0) %$%
    x
  
  if (length(x) > 0) {
    tail(x, 1)
  } else {
    1.0
  }
}

printf <- function(fmt, ...) sprintf(fmt, ...) %>% print

codes = read_tsv('../db/proc/code-map.txt', col_names=c('ndc', 'abx'))

all_pde_counts = read_tsv('bene_sum_2011.tsv') %>%
  mutate(bene=bene_id, all_pde=counts)

load.abxpde = function(fn) {
  read_tsv(fn) %>%
    rename(bene=BENE_ID, ndc=PRDSRVID) %>%
    inner_join(codes, by='ndc') %>%
    select(bene, abx)
}

abxpde = load.abxpde('../data/abx_pde_2011.tsv')

# convert to a wide abx pde data frame 
wap = abxpde %>%
  group_by(bene, abx) %>%
  summarize(counts=n()) %>%
  ungroup() %>%
  spread(abx, counts, fill=0)

wap$total_abx = wap %>% select(-bene) %>% rowSums

# merge in the information from all the pdes
wap %<>% inner_join(all_pde_counts, by='bene')

# do a summary over beneficiaries for all abx
all_abx_counts = function(abxpde) {
  abxpde %>%
    group_by(bene) %>%
    summarize(counts=n()) %$%
    counts
}

# filter for a specific abx before summarizing
abx_counts = function(this_abx, df) {
  df %>%
    filter(abx==this_abx) %>%
    group_by(bene) %>%
    summarize(counts=n()) %$%
    counts
}

add_zeros = function(n_total, nonzero_counts) {
  n_zero_benes = n_total_benes - length(nonzero_counts)
  c(nonzero_counts, rep(0, times=n_zero_benes)) %>% .[order(.)]
}

n_total_benes = 10341351
#nonzero_counts = all_abx_counts(abxpde)
#sprintf("Gini: %f\nOlesen: %f\n", ineq(counts, type='Gini'), olesen.index(counts))

plot_cdf = function(counts) {
  data.frame(counts=counts) %>%
    arrange(counts) %>%
    mutate(x=1:dim(.)[1], x=x/max(x), y=cumsum(counts)/sum(counts)) %>%
    ggplot(aes(x=x, y=y)) + geom_line() + coord_fixed()
}

plot_cdf2 = function(counts) {
  data.frame(counts=counts) %>%
    arrange(counts) %>%
    mutate(x=1:dim(.)[1], y=cumsum(counts)) %>%
    ggplot(aes(x=x, y=y)) + geom_line() + coord_fixed()
}

gini = function(x) ineq(x, type='Gini')

f = function(abx, df) {
  nonzero_counts = abx_counts(abx, df)
  counts = add_zeros(n_total_benes, nonzero_counts)
  list(abx=abx, total_usage=sum(counts), gini_all=gini(counts), olesen_all=olesen_index(counts),
       gini_nonzero=gini(nonzero_counts), olesen_nonzero=olesen_index(nonzero_counts))
}

abxs = abxpde %$% abx %>% unique
res = lapply(abxs, function(x) f(x, abxpde)) %>%
  simplify2array %>% t %>% as.data.frame %>%
  mutate(abx=as.character(abx),
         total_usage=as.numeric(total_usage),
         gini_all=as.numeric(gini_all),
         olesen_all=as.numeric(olesen_all),
         gini_nonzero=as.numeric(gini_nonzero),
         olesen_nonzero=as.numeric(olesen_nonzero))