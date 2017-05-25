import::from(ineq, Gini)

gini_change = function(x) {
  sx = sort(x, decreasing=TRUE)
  max_plus = sx
  max_plus[1] = max_plus[1] + 1
  max_minus = sx
  max_minus[1] = max_minus[1] - 1
  min_plus = c(sx, 1)
  min_minus = sx[-length(sx)]
  base = Gini(x)
  data_frame(base=base,
             max_plus=Gini(max_plus) - base,
             max_minus=Gini(max_minus) - base,
             min_plus=Gini(min_plus) - base,
             min_minus=Gini(min_minus) - base)
}