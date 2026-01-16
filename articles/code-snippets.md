# Code Snippets

## Introduction

This page is home to code that isn’t suitable for inclusion in the samc
package or dedicated vignettes, but that users may find beneficial
nonetheless. Some of what is included here may eventually be integrated
into the package. Code here is not as thoroughly tested as the package
and may not work in all situations or versions of the package.

## Reshaping output from pairwise()

``` r
# df will be a "long" format data.frame with columns "origin", "dest", and "result"
df <- pairwise(...)

# Use reshape2 to convert to a pairwise matrix
reshape2::acast(df, origin ~ dest, value.var = "result")
```
