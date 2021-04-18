library(samc)

# Load example data
res_data <- samc::ex_res_data
abs_data <- samc::ex_abs_data
occ_data <- samc::ex_occ_data


# Create samc-class object
samc_obj <- samc(res_data, abs_data,
                 tr_args = list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE))

# pairwise() example
pw <- pairwise(cond_passage, samc_obj, origin = 1:4, dest = 5)
print(pw)

# pairwise() example without dest
pw <- pairwise(dispersal, samc_obj, origin = c(2, 7))
print(pw)


