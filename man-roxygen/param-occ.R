#' @param occ The initial state \eqn{\psi} of the Markov chain. If the \code{\link{samc-class}}
#' objects was constructed using map inputs, then occ must be the same type of input
#' (either \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}),
#' and must have the same properties as initial map inputs (see the \code{\link{check}}
#' function).

# Descriptions in
# the Details section assume the values of \eqn{\psi} sum to 1 (\eqn{\sum_{i=1}^{n} \psi = 1}).
# In this case, the results of analyses represent probabilities. This is not required,
# however
