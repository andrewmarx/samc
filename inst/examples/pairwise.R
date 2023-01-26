library(samc)

# Load example data
res_data <- samc::example_split_corridor$res
abs_data <- samc::example_split_corridor$abs
occ_data <- samc::example_split_corridor$occ


# Create samc-class object
samc_obj <- samc(res_data, abs_data,
                 model = list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE))

# pairwise() example
pw <- pairwise(cond_passage, samc_obj, origin = 1:4, dest = 5)
print(pw)

# pairwise() example without dest
pw <- pairwise(dispersal, samc_obj, origin = c(2, 7))
print(pw)


