# Conditional Mean First Passage Time

Calculate the mean number of steps to first passage

## Usage

``` r
cond_passage(samc, init, origin, dest)

# S4 method for class 'samc,missing,missing,location'
cond_passage(samc, dest)

# S4 method for class 'samc,missing,location,location'
cond_passage(samc, origin, dest)

# S4 method for class 'samc,ANY,missing,location'
cond_passage(samc, init, dest)
```

## Arguments

- samc:

  A
  [`samc-class`](https://andrewmarx.github.io/samc/reference/samc-class.md)
  object created using the
  [`samc`](https://andrewmarx.github.io/samc/reference/samc.md)
  function.

- init:

  Sets the initial state \\\psi\\ of the transients states. Input must
  be able to pass the
  [`check`](https://andrewmarx.github.io/samc/reference/check.md)
  function when compared against the
  [`samc-class`](https://andrewmarx.github.io/samc/reference/samc-class.md)
  object. Can only contain positive finite values.

- origin:

  A positive integer or character name representing transient state
  \\\mathit{i}\\. Corresponds to row \\\mathit{i}\\ of matrix
  \\\mathbf{P}\\ in the
  [`samc-class`](https://andrewmarx.github.io/samc/reference/samc-class.md)
  object. When paired with the `dest` parameter, multiple values may be
  provided as a vector.

- dest:

  A positive integer or character name representing transient state
  \\\mathit{j}\\. Corresponds to column \\\mathit{j}\\ of matrix
  \\\mathbf{P}\\ in the
  [`samc-class`](https://andrewmarx.github.io/samc/reference/samc-class.md)
  object. When paired with the `origin` parameter, multiple values may
  be provided as a vector.

## Value

See Details

## Details

\\\tilde{t}=\tilde{B}\_j^{-1}\tilde{F}\tilde{B}\_j{\cdot}1\\

- **cond_passage(samc, dest)**

  The result is a vector where each element corresponds to a cell in the
  landscape, and can be mapped back to the landscape using the
  [`map`](https://andrewmarx.github.io/samc/reference/map.md) function.
  Element *i* is the mean number of steps before absorption starting
  from location *i* conditional on absorption into *j*

  Note that mathematically, the formula actually does not return a value
  for when *i* is equal to *j*. This leads to a situation where the
  resultant vector is actually one element short and the index for some
  of the elements may be shifted. The **cond_passage()** function fills
  inserts a `0` value for vector indices corresponding to *i == j*. This
  corrects the final result so that vector indices work as expected, and
  allows the vector to be properly used in the
  [`map`](https://andrewmarx.github.io/samc/reference/map.md) function.

- **cond_passage(samc, origin, dest)**

  The result is a numeric value representing the mean number of steps
  before absorption starting from a given origin conditional on
  absorption into *j*.

  As described above, mathematically the formula does not return a
  result for when the `origin` and `dest` inputs are equal, so the
  function simply returns a `0` in this case.

**WARNING**: This function is not compatible when used with data where
there are states with total absorption present. When present, states
representing total absorption leads to unsolvable linear equations. The
only exception to this is when there is a single total absorption state
that corresponds to input to the `dest` parameter. In this case, the
total absorption is effectively ignored when the linear equations are
solved.

**WARNING**: This function will crash when used with data representing a
disconnected graph. This includes, for example, isolated pixels or
islands in raster data. This is a result of the transition matrix for
disconnected graphs leading to some equations being unsolvable.
Different options are being explored for how to best identify these
situations in data and handle them accordingly.

## Performance

Performance details are in the performance vignette:
[`vignette("performance", package = "samc")`](https://andrewmarx.github.io/samc/articles/performance.md).

## Examples

``` r
# "Load" the data. In this case we are using data built into the package.
# In practice, users will likely load raster data using the raster() function
# from the raster package.
res_data <- samc::example_split_corridor$res
abs_data <- samc::example_split_corridor$abs
init_data <- samc::example_split_corridor$init


# Make sure our data meets the basic input requirements of the package using
# the check() function.
check(res_data, abs_data)
#> [1] TRUE
check(res_data, init_data)
#> [1] TRUE

# Setup the details for a random-walk model
rw_model <- list(fun = function(x) 1/mean(x), # Function for calculating transition probabilities
                 dir = 8, # Directions of the transitions. Either 4 or 8.
                 sym = TRUE) # Is the function symmetric?


# Create a `samc-class` object with the resistance and absorption data using
# the samc() function. We use the recipricol of the arithmetic mean for
# calculating the transition matrix. Note, the input data here are matrices,
# not RasterLayers.
samc_obj <- samc(res_data, abs_data, model = rw_model)


# Convert the initial state data to probabilities
init_prob_data <- init_data / sum(init_data, na.rm = TRUE)


# Calculate short- and long-term metrics using the analytical functions
short_mort <- mortality(samc_obj, init_prob_data, time = 50)
short_dist <- distribution(samc_obj, origin = 3, time = 50)
long_disp <- dispersal(samc_obj, init_prob_data)
#> 
#> Cached diagonal not found.
#> Performing setup. This can take several minutes... Complete.
#> Calculating matrix inverse diagonal...
#> Computing: 100% (done)                         
#> Complete                                                      
#> Diagonal has been cached. Continuing with metric calculation...
visit <- visitation(samc_obj, dest = 4)
surv <- survival(samc_obj)


# Use the map() function to turn vector results into RasterLayer objects.
short_mort_map <- map(samc_obj, short_mort)
short_dist_map <- map(samc_obj, short_dist)
long_disp_map <- map(samc_obj, long_disp)
visit_map <- map(samc_obj, visit)
surv_map <- map(samc_obj, surv)
```
