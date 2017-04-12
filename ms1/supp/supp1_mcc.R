
# load db data
condition_names = read_tsv('../../chronic-conditions/names.txt',
                           col_names=c('condition', 'condition_name'))

# load only data from 2011
usage = read_tsv('../usage_2011.tsv')

# mcc functions
mcc = function(true, pred) {
  ut = true - mean(true)
  up = pred - mean(pred)

  co = mean(ut * up)
  vt = mean(ut ** 2)
  vp = mean(up ** 2)

  if (vt * vp == 0) {
    NA
  } else {
    co / sqrt(vt * vp)
  }
}

mcc_cor = function(x) {
  x = as.matrix(x)
  labels = colnames(x)
  n = length(labels)
  out = matrix(nrow=n, ncol=n, dimnames=list(labels, labels))

  for (i in 1:n) out[i, i] = 1.0

  for (i in 1:(n-1)) {
    ui = x[,i] - mean(x[,i])
    vi = mean(ui ** 2)
    for (j in (i+1):n) {
      uj = x[,j] - mean(x[,j])
      vj = mean(uj ** 2)
      co = mean(ui * uj)
      val = co / sqrt(vi * vj)
      out[i, j] = val
      out[j, i] = val
    }
  }
  out
}

# Correlations between chronic conditions
correlations = usage %>%
  select(AMI:HYPOTH) %>%
  mcc_cor

save(correlations, file='correlations.RData')
