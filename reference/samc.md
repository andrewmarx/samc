# Create an samc object

Create an samc object that contains the absorbing Markov chain data

## Usage

``` r
samc(data, absorption, fidelity, model, options = NULL)

# S4 method for class 'SpatRaster,SpatRaster,SpatRaster,list'
samc(data, absorption, fidelity, model, options = NULL)

# S4 method for class 'RasterLayer,RasterLayer,RasterLayer,list'
samc(data, absorption, fidelity, model, options = NULL)

# S4 method for class 'SpatRaster,SpatRaster,missing,list'
samc(data, absorption, model, options = NULL)

# S4 method for class 'RasterLayer,RasterLayer,missing,list'
samc(data, absorption, model, options = NULL)

# S4 method for class 'matrix,matrix,matrix,list'
samc(data, absorption, fidelity, model, options = NULL)

# S4 method for class 'matrix,matrix,missing,list'
samc(data, absorption, model, options = NULL)

# S4 method for class 'dgCMatrix,missing,missing,missing'
samc(data, options = NULL)

# S4 method for class 'matrix,missing,missing,missing'
samc(data, options = NULL)
```

## Arguments

- data:

  A
  [`SpatRaster-class`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or
  [`RasterLayer-class`](https://rdrr.io/pkg/raster/man/Raster-classes.html)
  or [`matrix`](https://rdrr.io/r/base/matrix.html) or Matrix package
  dgCMatrix sparse matrix.

- absorption:

  A
  [`SpatRaster-class`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or
  [`RasterLayer-class`](https://rdrr.io/pkg/raster/man/Raster-classes.html)
  or [`matrix`](https://rdrr.io/r/base/matrix.html)

- fidelity:

  A
  [`SpatRaster-class`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or
  [`RasterLayer-class`](https://rdrr.io/pkg/raster/man/Raster-classes.html)
  or [`matrix`](https://rdrr.io/r/base/matrix.html)

- model:

  A list with args for constructing a transition matrix.

- options:

  A list of options that changes how the samc behaves computationally.

## Value

A
[`samc-class`](https://andrewmarx.github.io/samc/reference/samc-class.md)
object

## Details

This function is used to create a
[`samc-class`](https://andrewmarx.github.io/samc/reference/samc-class.md)
object. There are multiple options for creating this object.

**Option 1: Raster or Matrix Maps**

*`samc(data = matrix, absorption = matrix, fidelity = matrix, model = list())`*

*`samc(data = SpatRaster, absorption = SpatRaster, fidelity = SpatRaster, model = list())`*

*`samc(data = RasterLayer, absorption = RasterLayer, fidelity = RasterLayer, model = list())`*

The
[`samc-class`](https://andrewmarx.github.io/samc/reference/samc-class.md)
object can be created from a combination of resistance (or conductance),
absorption, and fidelity data. These different landscape data inputs
must be the same type (a matrix, SpatRaster, or RasterLayer), and have
identical properties, including dimensions, location of NA cells, and
CRS (if using raster inputs).

The `data` and `absorption` inputs are always mandatory for this
approach. The `fidelity` input is optional. If the `fidelity` input is
not provided, then it is assumed that there is no site fidelity (i.e.,
individuals will always move to an adjacent cell each time step).

The `model` parameter is a mandatory list with the following template:
`` list(fun = `function`, dir = `numeric`, sym = `logical`) ``. It is
used when calculating the values for the transition matrix. `fun` must
be a mathematical function taking a single vector `x` as input. `x[1]`
contains the value for node *i* and `x[2]` contains the value for node
*j* from the `data` input to `samc()`. `dir` determines which
neighboring nodes are used andmust be either `4` or `8`. `symmetric` is
an optimization to reduce redundant work for when `fun` is communative.
It must be either `TRUE` or `FALSE`.

Special cases for `fun` exist to significantly speed up creation of the
samc model. They are selected using specific strings as the value for
`fun`. Currently, `"1/mean(x)"` and `"x[2]"` are implemented. Other
special cases can be implemented on request.

When using raster inputs, SpatRaster objects (from the terra package)
are recommended over RasterLayer objects (from the raster package).
Internally, samc is using SpatRaster objects, which means RasterLayer
objects are being converted to SpatRaster objects, which is a source of
memory inefficiency.

**Option 2: P Matrix**

*`samc(data = matrix)`*

*`samc(data = dgCMatrix)`*

The `data` parameter can be used alone to create a
[`samc-class`](https://andrewmarx.github.io/samc/reference/samc-class.md)
object directly from a preconstructed P matrix. This matrix must be
either a base R matrix, or a sparse matrix (dgCMatrix format) from the
Matrix package. It must meet the following requirements:

- The number of rows must equal the number of columns (a square matrix)

- Total absorption must be a single column on the right-hand side of the
  matrix

- At the bottom of the matrix, there must be a row filled with 0's
  except for the last element (bottom-right of the matrix diagonal),
  which must be set to 1

- Every disconnected region of the matrix must have at least one
  non-zero absorbing value

- Each row must sum to 1

- All values must be in the range of 0-1

Additionally, the columns and rows of the P matrix may be named (e.g.,
using dimnames(), rowname(), colnames(), etc). When specifying `origin`
or `dest` inputs to metrics, these names may be used instead of cell
numbers. This has the advantage of making the code for an analysis
easier to read and interpret, which may also help to eliminate
unintentional mistakes. There are two requirements for naming the
rows/cols of a P matrix. First, since the P matrix represents a pairwise
matrix, the row and column names must be the same. Second, there must be
no duplicate names. The exception to these rules is the very last column
and the very last row of the P matrix. Since these are not part of the
pairwise transition matrix, they may have whatever names the user
prefers.

**Additional Options**

Additional options can be passed to `samc(..., options)` as a list.
There are several possible options:
`` list(method = `character`, threads = `numeric`, override = `logical`, precision = `character`) ``.
The use of options varies depending on other inputs. Some options can be
changed after model creation (see samc-class documentation). It can be
omitted for default behaviors.

`method` controls the type of mathematical algorithm used for the model.
It can be either `"direct"` (the default), `"iter"`, or `"conv"`. See
the Computation Methods vignette for details.

`threads` is a positive integer (default `1`) that enables
parallelization in certain limited cases. See the Parallel Computing
vignette for details.

`override` is a logical value (default: `FALSE`) that enables certain
memory intensive calculations. It only applies to the `"direct"` and
`"iter"` methods. See the samc-class reference for details.

`precision` controls the numerical precision of calculations. It can
either be `"double"` (~15 digits precision; the default) or `"double"`
(~7 digits of precision). Choosing the lower precision level can lead to
significant speed improvements, and reduces the memory requirements by
approximately half. It currently only applies to the `"conv"` method.

**Additional Information**

Depending on the data used to construct the samc-class object, some
metrics may cause crashes. This is a result of the underlying P matrix
having specific properties that make some equations unsolvable. One
known case is a P matrix that represents a disconnected graph, which can
lead to the
[`cond_passage()`](https://andrewmarx.github.io/samc/reference/cond_passage.md)
function crashing. In terms of raster/matrix inputs, a disconnected
graph occurs when one or more pixels/cells are unreachable from other
pixels/cells due to the presence of a full barrier made up of NA values.
In a raster, these may be obvious as islands but can be as inconspicuous
as a rogue isolated pixel. There may be other currently unknown
situations that lead to unsolvable metrics.

Future work is planned towards identifying these issues during the
creation of the samc-class object and handling them appropriately to
prevent inadvertent crashes.

**Version 3 Changes**

Support for creating samc-class objects from TransitionLayer objects was
removed so that the package is not dependent on gdistance.

**Version 2 Changes**

Version 1.5.0 officially removed support for the deprecated
`resistance`, `tr_fun`, `directions`, `p_mat`, `latlon`, and `override`
arguments. Old code will have to be updated to the new samc() function
structure in order to work.

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
