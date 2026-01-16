# Calculate survival metrics

Calculates the expected time to absorption

## Usage

``` r
survival(samc, init, origin)

# S4 method for class 'samc,missing,missing'
survival(samc)

# S4 method for class 'samc,missing,location'
survival(samc, origin)

# S4 method for class 'samc,ANY,missing'
survival(samc, init)
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

## Value

See Details

## Details

\\z=(I-Q)^{-1}{\cdot}1=F{\cdot}1\\

- **survival(samc)**

  The result is a vector \\\mathbf{v}\\ where \\\mathbf{v}\_i\\ is the
  expected time to absorption if starting at transient state
  \\\mathit{i}\\.

  If the samc-class object was created using matrix or RasterLayer maps,
  then vector \\\mathbf{v}\\ can be mapped to a RasterLayer using the
  [`map`](https://andrewmarx.github.io/samc/reference/map.md) function.

\\\psi^Tz\\

- **survival(samc, init)**

  The result is a numeric that is the expected time to absorption given
  an initial state \\\psi\\.

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
