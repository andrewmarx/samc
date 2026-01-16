# Troubleshooting

## Introduction

This document lists common errors users might encounter and how to
address them.

## Error in raster::compareRaster

This will occur if the
[`check()`](https://andrewmarx.github.io/samc/reference/check.md)
function detects that different sources of landscape data have
mismatched properties. This could be for a variety of reasons, including
differing coordinate reference systems (CRS) if using RasterLayers,
differing dimensions, and differing locations of NA data. An important
thing to be aware of is that the package does not support mixing
landscape data in different types of objects; all landscape data used in
an analysis must either be in matrices or RasterLayers. This includes
occupancy data that is not part of the creation of the `samc-class`
object.

``` r
# Working example
r1 <- res_data
r2 <- res_data
check(r1, r2)
#> [1] TRUE


# Remove the NA's in r2 by overwriting all the elements with the number 1.
# check() doesn't check the actual values of the data, but it does check the
# location of NA's
r1 <- res_data
r2 <- res_data
r2[1] <- 1
check(r1, r2)
#> Error:
#> ! NA mismatch in input data


# Change the dimensions of r2 by subsetting it. check() ensures that the data
# inputs have the same number of rows and columns
r1 <- res_data
r2 <- res_data
r2 <- r2[1:5, 1:5]
check(r1, r2)
#> Error in `h()`:
#> ! error in evaluating the argument 'a' in selecting a method for function 'check': [rast] extents do not match
```

## Unable to find an inherited method

The SAMC package makes use of R’s S4 method dispatch system to enforce
how the different versions of each function can be used. Hopefully, this
will ensure that users are not unintentionally misusing functions by
producing an error rather than letting the code run and returning a
result that may not be obviously incorrect. If this error appears, it
means that the combination of arguments provided to the function are not
valid. This could be for a couple different reasons. The first is that
an optional parameter is skipped without specifying argument names for
the subsequent arguments. The second is that a user is passing the wrong
type of data to an argument (e.g., passing a `numeric` when the function
expects a `RasterLayer` or a `matrix`).

``` r
# Example: Skipping optional arguments. In this case, the `fidelity` argument is
# optional, so we skip it. The model argument, however, is always required, 
# so we pass the relevant to it. But because we don't specify which
# argument it is, R is trying to find a version of the function that expects it
# as the third argument, but this version does not exist.
samc_obj <- samc(res_data, abs_data, list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE))
#> Error:
#> ! unable to find an inherited method for function 'samc' for signature 'data = "matrix", absorption = "matrix", fidelity = "list", model = "missing"'

# Solution
samc_obj <- samc(res_data, abs_data, model = list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE))


# Example: Incorrect input types. In this case, we are attempting to pass a
# single numeric value as absorption data. However, the absorption data must 
# always be in a matrix or RasterLayer object
samc_obj <- samc(res_data, 0.01, model = list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE))
#> Error:
#> ! unable to find an inherited method for function 'samc' for signature 'data = "matrix", absorption = "numeric", fidelity = "missing", model = "list"'

# Solution
samc_obj <- samc(res_data, abs_data, model = list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE))
```

## Error: All disconnected regions must have at least one non-zero absorption value

See the [Disconnected
Data](https://andrewmarx.github.io/samc/articles/data-disconnected.md)
vignette.

## Warning: Input contains disconnected regions. This does not work with the cond_passage() metric

See the [Disconnected
Data](https://andrewmarx.github.io/samc/articles/data-disconnected.md)
vignette.

## Error in \<…\> : function ‘Rcpp_precious_remove’ not provided by package ‘Rcpp’

Update the Rcpp package. Then make sure to unload all of your packages
and reload them again. You do not need to load the Rcpp package (i.e, do
not [`library(Rcpp)`](https://www.rcpp.org)).
