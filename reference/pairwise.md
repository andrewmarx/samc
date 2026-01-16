# Pairwise analyses

Analysis for pairwise combinations locations

## Usage

``` r
pairwise(fun, samc, origin, dest)

# S4 method for class 'function,samc,location,location'
pairwise(fun, samc, origin, dest)

# S4 method for class 'function,samc,location,missing'
pairwise(fun, samc, origin)
```

## Arguments

- fun:

  A samc analytical function with signature fun(samc, origin, dest)

- samc:

  A
  [`samc-class`](https://andrewmarx.github.io/samc/reference/samc-class.md)
  object

- origin:

  A vector of locations

- dest:

  A vector of locations. Can be excluded to reuse the origin parameter

## Value

A 'long' format data.frame

## Details

When providing vector inputs for the \`origin\` and \`dest\` parameters
to analytical functions, the package assumes that users are providing
pairs of \`origin\` and \`dest\`. That is, \`origin\[1\]\` is paired
with \`dest\[1\]\`, \`origin\[2\]\` is paired \`dest\[2\]\`, etc.
Another way to think about it is that these two vector inputs can be
treated as columns in the same dataframe. The result of the analytical
function then is a vector of the same length as the input. This behavior
works for any situation, so it is the default for the package.

However, some users may wish to run an analytical function for all the
pairwise combinations of the values in the input vectors. That is,
\`origin\[1\]\` is paired with \`dest\[1\]\`,\`dest\[2\]\`,
\`dest\[3\]\`, etc, before moving on to the next elements in \`origin\`.
This approach has the advantage of potentially reducing the amount of
code needed for an analysis, and the results can be represented as a
pairwise matrix, but it is not suitable for all situations. To enable
this second approach more easily, the \`pairwise()\` function runs all
the combinations of the \`origin\` and \`dest\` parameters for an
analytical function and returns the results in a 'long' format
data.frame. This data.frame can then be reshaped into a pairwise matrix
or 'wide' format data.frame using tools like the reshape2 or tidyr
packages.

This function is not intended to be used with other inputs such as
\`init\` or \`time\`

## Examples

``` r
library(samc)

# Load example data
res_data <- samc::example_split_corridor$res
abs_data <- samc::example_split_corridor$abs


# Create samc-class object
samc_obj <- samc(res_data, abs_data,
                 model = list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE))

# pairwise() example
pw <- pairwise(cond_passage, samc_obj, origin = 1:4, dest = 5)
print(pw)
#>   origin dest   result
#> 1      1    5 5402.175
#> 2      2    5 4642.001
#> 3      3    5 3600.117
#> 4      4    5 2249.220

# pairwise() example without dest
pw <- pairwise(dispersal, samc_obj, origin = c(2, 7))
print(pw)
#>   origin dest    result
#> 1      2    2 0.7908536
#> 2      7    2 0.6243287
#> 3      2    7 0.6052441
#> 4      7    7 0.7972468

```
