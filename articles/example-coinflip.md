# Coin Flip

## Introduction

This example is a write-up of an answer I gave to a [coin flip question
on
StackExchange](https://stats.stackexchange.com/questions/553289/could-any-equation-have-predicted-the-results-of-this-simulation).
Essentially, the question was: given a series of coin flips, how many
coin flips would it take (on average) to get the sequence
Heads-Tails-Heads.

It’s possible to model this scenario using absorbing Markov chains.
Metrics within the package can then be used to answer the posed
question, as well as explore other questions that may be of interest.

The complete code for this example is available on
[Github](https://github.com/andrewmarx/samc/blob/main/vignettes/scripts/example-coinflip.R).

## Libraries

``` r
library(samc)
```

## Setup

The first step is to create a P matrix representing the possible states.
The key is to define HTH as our absorbing state. In other words, once
the model reaches that state, it is complete.

![Matrix of transition probabilities between different sets of three
coin flips.](img/coinflipmatrix.svg)

Matrix of transition probabilities between different sets of three coin
flips.

``` r
# Coin flip probabilities
p <- 0.5    # Probability of heads
q <- 1 - p  # Probability of tails

p_mat <- matrix(c(0, 0, 0, 0, 0, 0, q, p,
                  q, p, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, q, p,
                  q, p, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, q, p, 0, 0,
                  0, 0, q, p, 0, 0, 0, 0,
                  0, 0, 0, 0, q, p, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 1),
                8, byrow = TRUE)

rownames(p_mat) <- c("HHT", "HHH", "THT", "THH", "TTT", "TTH", "HTT", "HTH")
colnames(p_mat) <- rownames(p_mat)

# A samc object is the core of the package
samc_obj <- samc(p_mat)
#> [1] "Warning: Some checks for manually created P matrices are still missing:"
#> [1] "1) Discontinuous data will not work with the cond_passage() function."
#> [1] "2) Every disconnected region of the graph must have at least one non-zero absorption value."
```

First, the coin flip probabilities are defined using variables, which
allows us to easily change the matrix to simulate a biased coin if we
want.

The last row and column of the matrix represents HTH, or the absorbing
state. Notice that the last row makes it so that HTH can only transition
to itself, and never to a different state. This construction is
explained in the [Overview
background](https://andrewmarx.github.io/samc/articles/overview.html#background-1).

The remaining rows/columns represent the remaining possible combinations
of 3 coin flips. Technically, this model assumes our starting point
includes a sequence of 3 previous flips. However, if we decide to use
TTT as our initial starting point, it will take at least 3 flips to get
to the target sequence, which is no different than starting from an
empty sequence of flips. Contrast this with the sequence of HHT, which
only requires one more H flip for the HTH sequence to occur.

Finally, notice that the rows sum to `1`, which is required for the P
matrix. This makes sense because each row represents the transition
probabilities for a state. Each state has a 100% probability of doing
*something*, whether it’s changing to a different transient state, an
absorbing state, or back to the current state.

For the original question, it’s possible to collapse some of these
states and simplify the matrix by only considering the last two coin
flips instead of three. For example, HTT would also require a minimum of
3 flips to reach (like TTT, which has the same last two flips), so it
could also be used as a starting point to answer the original question.
However, the current construction simplifies the interpretation a little
and allows different questions to be asked about the model.

## Metrics

To answer the original question of how many coin flips it will take, on
average, to see HTH, all we need to do is calculate the expected time to
absorption. With the samc package, that is accomplished using the
[`survival()`](https://andrewmarx.github.io/samc/reference/survival.md)
function.

``` r
# Given the last 3 flips, how many more flips until we hit HTH (absorption)?
survival(samc_obj)
#> [1]  6  8  6  8 10  8 10
```

The fifth element represents TTT and the seventh element represents HTT,
or the starting points that do not have progress toward the result. They
both have a value of `10`, which indicates that it will take an average
of `10` coin flips with a fair coin to see the sequence of HTH.

The remaining elements give us insight into what would be expected if we
have already made some coin flips with some progress toward the goal.
For example, the first and third elements represent, HHT and THT, which
are only one flip away from HTH. As a result, they have the lowest
expected averages to completion, which should intuitively make sense.
But they still are relatively high at `6` flips expected. Even though
there is a 50% probability of completing the sequence in a single flip,
a failure would put us back in the HTT state, which is an average of
`10` steps away.

There are other aspects of the model that can be explored using the
package. For example, the
[`visitation()`](https://andrewmarx.github.io/samc/reference/visitation.md)
metric calculates how many times each state is expected to be visited
before absorption. In the context of the coin flip model, this means we
can calculate how many times we expect a sequence to occur before
hitting HTH.

``` r
# Given a starting point (in this case, a sequence of 3 flips), how many times
# would we expect the different combinations of 3 flips to occur before absorption?
visitation(samc_obj, origin = "HHT")
#> [1] 1.5 0.5 0.5 0.5 1.0 1.0 1.0

sum(visitation(samc_obj, origin = "HHT")) # Compare to survival() result
#> [1] 6
```

In this case, we specify a starting sequence of HHT, which is one flip
away from HTH (absorption). As an example interpretation, with this
starting point, we expect the sequence HHT to occur one and a half
times, on average, before hitting the HTH sequence. What’s interesting
mathematically is that the sum of these visitation values is the same as
the expected remaining flips before hitting HTH.

Ultimately, the
[`visitation()`](https://andrewmarx.github.io/samc/reference/visitation.md)
metric produces a matrix of all the possible combinations of states and
specifying an origin point simply gives us a single row of the matrix.
It’s possible to get individual columns, as well as the entire matrix
(which requires disabling a safety measure in place for larger
problems).

``` r
# Instead of a start point, we can look at an endpoint and how often we expect
# it to occur for each of the possible starting points
visitation(samc_obj, dest = "THT")
#> [1] 0.5 0.5 1.5 0.5 1.0 1.0 1.0

# These results are just rows/cols of a larger matrix. We can get the entire matrix
# of the start/end possibilities but first, we have to disable some safety measures
# in place because this package is designed to work with extremely large P matrices
# (millions of rows/cols) where these types of results will consume too much RAM and
# crash R
samc_obj$override <- TRUE
visitation(samc_obj)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> [1,]  1.5  0.5  0.5  0.5    1    1    1
#> [2,]  1.5  2.5  0.5  0.5    1    1    1
#> [3,]  0.5  0.5  1.5  0.5    1    1    1
#> [4,]  1.5  1.5  0.5  1.5    1    1    1
#> [5,]  1.0  1.0  1.0  1.0    3    2    1
#> [6,]  1.0  1.0  1.0  1.0    1    2    1
#> [7,]  1.0  1.0  1.0  1.0    2    2    2

rowSums(visitation(samc_obj)) # equivalent to survival() above
#> [1]  6  8  6  8 10  8 10
```

A third relevant metric is
[`dispersal()`](https://andrewmarx.github.io/samc/reference/dispersal.md).
Rather than the expected number of times we would expect a coin flip
sequence to occur, it calculates the probability of a sequence occurring
before hitting our target sequence of HTH. The short-term variant can be
used to specify a limit to the number of flips we will make.

``` r
dispersal(samc_obj, dest = "TTT", time = 5)
#> [1] 0.28125 0.21875 0.28125 0.21875 0.00000 0.21875 0.59375
```

In this case, we’re trying to find the probability of TTT occurring
within 5 flips. The results show the probability for every possible
starting sequencing. For example, the seventh element, HTT is only one
flip away and has the highest probability of `0.59375`.
