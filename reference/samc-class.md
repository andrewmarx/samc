# samc class

S4 class to manage SAMC data.

## Details

The samc class is used to help ensure that the package is used correctly
and to minimize the possibility for users to accidentally produce
nonsensical results that may not be obviously incorrect. This class
contains the p matrix necessary for all the calculations in the package,
and enforces its type so that users are less likely to inadvertently
alter it in a way that will cause issues in calculations.

The [`samc()`](https://andrewmarx.github.io/samc/reference/samc.md)
function is used to create `samc-class` objects.

The samc-class slots are subject to change, so users should not be using
the `@` operator to access or change them. Doing so leads to the risk of
broken code in the future. Instead, where relevant, the `$` operator can
be used to get and set components of the class safely. This is a current
list of what can be accessed and modified in the class:

- **override**

  Some analyses are memory intensive and have the potential to make a
  user's system non-responsive or crash. By default, a samc-class object
  cannot be used in these analyses to prevent unintentional loss of
  work. In some cases, users may wish to use these particular analyses,
  in which case this behavior can be overridden. To get the current
  state of the override, use `samc_obj$override`. To enable the use of
  the analyses, the override can be set to `TRUE` using
  `samc_obj$override <- TRUE`. Before enabling the override, users
  should familiarize themselves with the Performance vignette.

- **q_matrix**

  Advanced users may wish to have direct access to the Q matrix for
  developing custom calculations/analyses. Assumptions should not be
  made about the internal storage and management of the P and Q matrices
  in the samc-class; these things are subject to change in the future.
  To safely access the Q matrix, use `samc_obj$q_matrix`. The Q matrix
  inside of the samc-class cannot be modified.

- **p_matrix**

  `samc_obj$p_matrix` can be used to get a copy of the P matrix.

- **abs_states**

  Used to attach additional absorbing states to an samc object. This
  does not cause P/Q matrices to be updated. Instead, it is intended to
  provide decomposed results from the
  [`mortality()`](https://andrewmarx.github.io/samc/reference/mortality.md)
  and
  [`absorption()`](https://andrewmarx.github.io/samc/reference/absorption.md)
  metrics for different sources of absorption that might be contributing
  to the total absorption values that were used to create the samc
  object.

  The input must be in the same form as the absorption inputs used in
  [`samc()`](https://andrewmarx.github.io/samc/reference/samc.md).
  Matrices are passed in as a `list`, and rasters are passed in as a
  `RasterStack`. Using `NA` as the input will reset it.

- **solver**

  `samc_obj$solver` can be used to change the default linear algebra
  solver used in some of the metrics. The default value of "direct"
  means a direct solver is used, and is what was used in previous
  versions of the package. The alternative value of "iter" switches the
  package to an iterative solver, which is significantly more memory
  efficient for larger datasets, but in general will be noticeably
  slower depending on patterns in the data.

- **threads**

  `samc_obj$threads` can be used to get or set the number of threads
  used for parallel computations. Details can be found in the Parallel
  Computing vignette.

## Slots

- `data`:

  Data associated with different components of the P matrix

- `conv_cache`:

  Convolution cache

- `model`:

  List containing model info used to build the samc object

- `source`:

  Information about the data source for the P matrix

- `nodes`:

  The number of nodes in the graph

- `map`:

  Used to verify landscape inputs and mapping of vector data

- `crw_map`:

  Matrix used to map location and direction to edges description

- `prob_mat`:

  Matric for CRW probabilities

- `names`:

  Names of the transient states

- `clumps`:

  Number of discontinuous regions in data

- `override`:

  Used to prevent accidental use of memory intensive functions

- `solver`:

  Controls the linear solver used for relevant metrics

- `threads`:

  Used for multi-threading

- `precision`:

  "single" or "double". Only applies to convolution currently.

- `.cache`:

  Cached data for performance boosts
