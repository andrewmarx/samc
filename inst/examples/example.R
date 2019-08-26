# "Load" the data. In this case we are using data built into the package.
# In practice, users will likely load raster data using the raster() function
# from the raster package.
res_data <- samc::ex_res_data
abs_data <- samc::ex_abs_data
occ_data <- samc::ex_occ_data


# Make sure our data meets the basic input requirements of the package using
# the check() function.
check(res_data, abs_data)
check(res_data, occ_data)


# Create a `samc-class` object with the resistance and absorption data using
# the samc() function. We use the recipricol of the arithmetic mean for
# calculating the transition matrix. Note, the input data here are matrices,
# not RasterLayers. If using RasterLayers, the latlon parameter must be set.
samc_obj <- samc(res_data, abs_data, tr_fun = function(x) 1/mean(x))


# Convert the occupancy data to probability of occurrence
occ_prob_data <- occ_data / sum(occ_data, na.rm = TRUE)


# Calculate short- and long-term metrics using the analytical functions
short_mort <- mortality(samc_obj, occ_prob_data, time = 50)
short_dist <- distribution(samc_obj, origin = 3, time = 50)
long_disp <- dispersal(samc_obj, occ_prob_data)
visit <- visitation(samc_obj, dest = 4)
surv <- survival(samc_obj)


# Use the map() function to turn vector results into RasterLayer objects.
short_mort_map <- map(samc_obj, short_mort)
short_dist_map <- map(samc_obj, short_dist)
long_disp_map <- map(samc_obj, long_disp)
visit_map <- map(samc_obj, visit)
surv_map <- map(samc_obj, surv)
