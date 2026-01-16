# Check landscape data

Check that landscape inputs have valid values and matching properties.

## Usage

``` r
check(a, b)

# S4 method for class 'Raster,missing'
check(a)

# S4 method for class 'SpatRaster,missing'
check(a)

# S4 method for class 'matrix,missing'
check(a)

# S4 method for class 'SpatRaster,SpatRaster'
check(a, b)

# S4 method for class 'Raster,Raster'
check(a, b)

# S4 method for class 'matrix,matrix'
check(a, b)

# S4 method for class 'samc,Raster'
check(a, b)

# S4 method for class 'samc,SpatRaster'
check(a, b)

# S4 method for class 'samc,matrix'
check(a, b)

# S4 method for class 'samc,numeric'
check(a, b)
```

## Arguments

- a:

  A
  [`samc-class`](https://andrewmarx.github.io/samc/reference/samc-class.md),
  [`matrix`](https://rdrr.io/r/base/matrix.html), or
  [`RasterLayer-class`](https://rdrr.io/pkg/raster/man/Raster-classes.html)
  object

- b:

  A [`matrix`](https://rdrr.io/r/base/matrix.html) or
  [`RasterLayer-class`](https://rdrr.io/pkg/raster/man/Raster-classes.html)
  object

## Value

See *Details* section.

## Details

This function is used to ensure that inputs (resistance, absorption,
fidelity, and occupancy) have valid values and the same properties. This
includes checking the CRS (if using raster inputs), dimensions, and
locations of cells with NA data. It can be used to directly compare two
matrices or two rasters, or it can be used to check a
[`samc-class`](https://andrewmarx.github.io/samc/reference/samc-class.md)
object against a matrix or raster.

It can also be used to check a numeric vector against a
[`samc-class`](https://andrewmarx.github.io/samc/reference/samc-class.md)
object created from a P matrix. In this case, the length of the vector
must be equal to the number of transient states. If the transient states
are named, the vector must contain the same names.

The function returns `TRUE` if the inputs have matching properties.
Otherwise, it will stop execution and print an error message with
details about the difference between the two inputs.

Note that the package assumes the different landscape inputs will be the
same type, either matrices or RasterLayers. Mixing RasterLayer data and
matrix data is not supported.

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
