#' @param time A positive integer or a vector of positive integers representing
#' \eqn{\mathit{t}} time steps. Vectors must be ordered and contain no duplicates.
#' Vectors may not be used for metrics that return dense matrices. The maximum time
#' step value is capped at 10,000 due to numerical precision issues.
