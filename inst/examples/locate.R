library(samc)
library(raster)


# Load example data
res_data <- samc::ex_res_data
abs_data <- samc::ex_abs_data
occ_data <- samc::ex_occ_data


# Create samc-class object
samc_obj <- samc(res_data, abs_data, tr_fun = function(x) 1/mean(x))


# We can use locate() to return a raster with the cell numbers encoded as data
# in the cells
cell_raster <- locate(samc_obj)
plot(cell_raster)


# We can use a variety of spatial inputs to get cell numbers using locate()
# The simplest is a two-column data.frame
coords <- data.frame(x = c(50, 79, 22),
                     y = c(25, 11, 19))
print(coords)
locate(samc_obj, coords)

# You will get an error if you input a coordinate that does not correspond
# to a non-NA cell
coords <- data.frame(x = c(1),
                     y = c(1))
print(coords)
try(locate(samc_obj, coords))
