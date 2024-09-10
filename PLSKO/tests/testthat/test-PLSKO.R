# test_that("multiplication works", {
#   expect_equal(2 * 2, 4)
# })

test_that("output dim is the same with input X",{
  test.X <- matrix(rnorm(100), 10, 10)
  output.X <- plsko(test.X)
  expect_equal(dim(test.X)[1], dim(output.X)[1])
  expect_equal(dim(test.X)[2], dim(output.X)[2])
})
