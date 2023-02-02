context("Check")


# Fail cases for check()
check_fail <- list(
  check1 = matrix(c(5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5),
               nrow = 8),
  check2 = matrix(c(5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,Inf,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5),
                  nrow = 9),
  check3 = matrix(c(5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,NaN,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5),
                  nrow = 9),
  check4 = matrix(c(5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5, NA,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5,
                    5,  5,  5,  5,  5,  5,  5,  5,  5),
                  nrow = 9)
)

test_that("check() works as intended", {
  # Check that the testlist scenarios all work
  for (t in testlist) {
    expect_true(check(t$res, t$init)) # Not really needed because samc() uses it, but included anyway
    expect_true(check(t$samc, t$init))
  }

  # Check that fail cases are producing errors
  r <- testlist[[1]]$res
  s <- testlist[[1]]$samc
  for (cf in check_fail) {
    expect_error(check(r, cf)) # Not really needed because samc() uses it, but included anyway
    expect_error(check(s, cf))
  }
})
