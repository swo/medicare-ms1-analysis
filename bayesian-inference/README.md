# Background

Marc's model idea: Each state generates the probabilities in each of its
hospitals, and the number of isolates gives the estimation for the probability
in that hospital.

# Results

## Doing the full model is cumbersome

I tried working out a model like what is shown in the pdf. This was a huge
hassle.  The integrals tended to get really small, and the process of selecting
the correct function for the relationship between consumption and resistance
was pretty fraught: linear models would give values outside 0 and 1, for
example.

## Weighted mean resistances are slightly better than unweighted mean resistances

To make things simpler, I looked quinolones and *E. coli*.

1. I computed the mean resistances in each state.
1. I computed the mean resistances weighted by $\sqrt{n}$, the number of isolates, in each state
1. I computed the mean of the posterior beta-binomial distribution for each state, using a uniform prior over $[0, 1]$ for $\mu$ (the mean) and an exponential prior on $\nu$ (the strength).

When comparing these values, I found that the unweighted mean slightly
overestimates the mean predicted from the Bayesian method, and that the
weighted mean performs slightly better.

This might be something good to throw into a supplement.
