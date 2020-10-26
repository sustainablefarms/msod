# test bugsvar2array
context("Conversion to and from BUGS variable format")

# build test array
idx <- expand.grid(row = 1:4, col = 1:3)
bugsnames <- paste0("V[", idx$row,",", idx$col, "]")
draws <- matrix(bugsnames, ncol = nrow(idx), nrow = 5, byrow = TRUE)
draws <- apply(draws, 2, paste0, paste0("d=",1:5))
colnames(draws) <- bugsnames

test_that("bugsvar2array works for array", {
  out <- bugsvar2array(draws, "V", 1:4, 1:3)
  expect_true(all(grepl("d=5", out[,,5], fixed = TRUE)))
  expect_true(all(grepl("[1,", out[1,,], fixed = TRUE)))
  expect_true(all(grepl(",3]", out[,3,], fixed = TRUE)))
  expect_false(all(grepl(",3]", out[,,4], fixed = TRUE)))
  expect_equal(dim(out), c(4, 3, 5))
})

test_that("bugsvar2array works for a vector", {
  theta <- draws[1, , drop = TRUE]
  out <- bugsvar2array(theta, "V", 1:4, 1:3)
  expect_true(all(grepl("d=1", out[,,1], fixed = TRUE)))
  expect_true(all(grepl("[1,", out[1,,], fixed = TRUE)))
  expect_true(all(grepl(",3]", out[,3,], fixed = TRUE)))
  expect_false(all(grepl(",2]", out[2,,], fixed = TRUE)))
  expect_equal(dim(out), c(4, 3, 1))
})
