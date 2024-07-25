# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

# This script is the source of the code for the coin flip example vignette
# The lines with `@knitr` are for processing by knitr

# TODO Incorporate more comments to make the script more useful on it's own
# without the vignette


## @knitr Part1
#
# Part 1 ----
#

## @knitr library_1
library(samc)


## @knitr setup_1
# Coin flip probabilities
p <- 0.5    # Probability of heads
q <- 1 - p  # Probability of tails

p_mat <- matrix(c(0, 0, 0, 0, 0, 0, q, p,
                  q, p, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, q, p,
                  q, p, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, q, p, 0, 0,
                  0, 0, q, p, 0, 0, 0, 0,
                  0, 0, 0, 0, q, p, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 1),
                8, byrow = TRUE)

rownames(p_mat) <- c("HHT", "HHH", "THT", "THH", "TTT", "TTH", "HTT", "HTH")
colnames(p_mat) <- rownames(p_mat)

# A samc object is the core of the package
samc_obj <- samc(p_mat)


## @knitr metric_1
# Given the last 3 flips, how many more flips until we hit HTH (absorption)?
survival(samc_obj)


## @knitr metric_2_1
# Given a starting point (in this case, a sequence of 3 flips), how many times
# would we expect the different combinations of 3 flips to occur before absorption?
visitation(samc_obj, origin = "HHT")

sum(visitation(samc_obj, origin = "HHT")) # Compare to survival() result


## @knitr metric_2_2
# Instead of a start point, we can look at an endpoint and how often we expect
# it to occur for each of the possible starting points
visitation(samc_obj, dest = "THT")

# These results are just rows/cols of a larger matrix. We can get the entire matrix
# of the start/end possibilities but first, we have to disable some safety measures
# in place because this package is designed to work with extremely large P matrices
# (millions of rows/cols) where these types of results will consume too much RAM and
# crash R
samc_obj$override <- TRUE
visitation(samc_obj)

rowSums(visitation(samc_obj)) # equivalent to survival() above


## @knitr metric_3
dispersal(samc_obj, dest = "TTT", time = 5)

