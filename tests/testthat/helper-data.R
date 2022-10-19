# testthat sources files beginning with helper- before running unit tests.
# This helper file contains data for each of the unit tests so that it is not
# repeated across test files.

baselist <- list(
  res = matrix(c( 1,  1,  2,  2,  3,  3,  2,  2,  1,
                  1,  2,  2,  3, 10, 10,  3,  2,  2,
                  1,  2,  2,  3, 10, 10,  3,  2,  2,
                  1,  1,  1,  3, 10, 10,  3,  1,  1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  1,  1, 10, 10, 10, 10, 10, 10,  1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1),
                nrow = 9),
  abs = matrix(c( 1,  2,  1,  1,  1,  1,  1,  1,  1,
                  1,  1,  1,  2,  2,  1,  1,  0,  1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  1,  5,  1,  1,  1,  1,  1,  0,  0,
                  1,  1,  1,  1,  1,  1,  1,  1,  0,
                  1,  1,  1,  1,  1,  1,  1,  1,  0,
                  1,  1,  1,  2,  1,  1,  1,  7,  7,
                  3,  1,  1,  1,  1,  1,  1,  1,  1),
                nrow = 9) / 1000,
  fid = matrix(c(10, 10, 10, 10, 10, 10, 10, 10, 10,
                 10, 10, 10, 10,  1,  1, 10, 10, 10,
                 10, 10, 10, 10,  1,  1, 10, 10, 10,
                 10, 10, 10, 10, 10, 10, 10, 10, 10,
                 10, 10, 10, 10, 10, 10, 50, 50, 50,
                 10, 10, 10, 10, 10, 10, 50, 50, 50,
                 10, 10, 10, 10, 10, 10, 10, 10, 10,
                 10, 10, 10, 10, 10, 10, 10, 10, 10,
                 10, 10, 10, 10, 10, 10, 10, 10, 10),
                nrow = 9) / 100,
  occ = matrix(c( 0,  0,  0,  0,  0,  0,  0,  0,  0,
                  0,  0,  0,  0,  0,  0,  0,  0,  0,
                  0,  0,  0,  0,  0,  0,  0,  0,  0,
                  0,  1,  1,  0,  0,  0,  0,  0,  0,
                  0,  1,  1,  0,  0,  0,  0,  0,  0,
                  0,  1,  1,  0,  0,  0,  0,  0,  0,
                  0,  1,  1,  0,  0,  0,  0,  0,  0,
                  0,  0,  0,  0,  0,  0,  0,  0,  0,
                  0,  0,  0,  0,  0,  0,  0,  0,  0),
                nrow = 9)
)


# Raster masks for testing different scenarios
masklist <- list(
  mask1 = matrix(c( 1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1),
                  nrow = 9),
  mask2 = matrix(c( 1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1, NA, NA, NA, NA,  1,  1,
                    1,  1,  1, NA, NA, NA, NA,  1,  1,
                    1,  1,  1, NA, NA, NA, NA,  1,  1,
                    1,  1,  1, NA, NA, NA, NA,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1),
                  nrow = 9),
  mask3 = matrix(c( 1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                   NA, NA, NA, NA, NA, NA, NA, NA, NA,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1),
                  nrow = 9),
  mask4 = matrix(c( 1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1, NA, NA, NA,  1,  1,
                    1,  1,  1,  1, NA,  1, NA,  1,  1,
                    1,  1,  1,  1, NA, NA, NA,  1,  1,
                    1, NA, NA, NA,  1,  1,  1,  1,  1,
                    1, NA,  1, NA,  1,  1,  1,  1,  1,
                    1, NA, NA, NA,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1),
                  nrow = 9)
)

p1 <- runif(81, max = 0.4)
p2 <- 1 - p1

testlist <- list()
for(i in 1:length(masklist)) {
  print(i)
  testlist[[i]] <- lapply(baselist, function(x) {masklist[[i]] * x})

  testlist[[i]]$length <- sum(!is.na(testlist[[i]]$res))

  testlist[[i]]$samc <- samc(testlist[[i]]$res,
                             testlist[[i]]$abs,
                             testlist[[i]]$fid,
                             tr_args = list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE))


  testlist[[i]]$samc@names = as.character(1:length(testlist[[i]]$samc@data@t_abs))

  testlist[[i]]$samc$abs_states <- list(testlist[[i]]$abs * p1, testlist[[i]]$abs * p2)

  testlist[[i]]$id <- i
}

# Asymmetric versions
n <- length(testlist)
for(i in (n + 1):(n + length(masklist))) {
  print(i)
  testlist[[i]] <- lapply(baselist, function(x) {masklist[[i - n]] * x})

  testlist[[i]]$length <- sum(!is.na(testlist[[i]]$res))

  testlist[[i]]$samc <- samc(testlist[[i]]$res,
                             testlist[[i]]$abs,
                             testlist[[i]]$fid,
                             tr_args = list(fun = function(x) 1/(mean(x) + x[1]), dir = 4, sym = FALSE))

  testlist[[i]]$samc@names = as.character(1:length(testlist[[i]]$samc@data@t_abs))

  testlist[[i]]$samc$abs_states <- list(testlist[[i]]$abs * p1, testlist[[i]]$abs * p2)
  testlist[[i]]$id <- i
}


time = 100
time_vec = c(3, 5, 7, 11, 13)
row_vec = c(7, 34, 5, 5)
col_vec = c(13, 13, 5, 19)
