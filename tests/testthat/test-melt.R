test_that("can parse packer names", {
  res <- parse_names(c("S", "R[2]", "x[1,3]", "beta<a>", "x[2]<b>"))
  expect_equal(res$variable, c("S", "R", "x", "beta", "x"))
  expect_equal(res$index,
               rbind(c(NA, NA), c(2, NA), c(1, 3), c(NA, NA), c(2, NA)))
  expect_equal(res$group, c(NA, NA, NA, "a", "b"))
})

test_that("names from real packers round-trip", {
  p <- monty::monty_packer("a", array = list(b = 2, m = c(2, 3)))
  res <- parse_names(p$names())
  idx <- p$index()
  expect_equal(res$variable, rep(c("a", "b", "m"), c(1, 2, 6)))
  expect_equal(res$index[idx$m, ], unname(which(idx$m > 0, arr.ind = TRUE)))
})

test_that("melt matches array indexing", {
  x <- array(runif(3 * 4 * 5), c(3, 4, 5))
  res <- melt_state(x, c("a", "b[1]", "b[2]"),
                    list(particle = 1:4, time = seq(0, 2, by = 0.5)))
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("particle", "time", "variable", "index1", "value"))
  expect_equal(nrow(res), 60)
  r <- res[res$variable == "b" & res$index1 == 2 & res$particle == 3 &
             res$time == 1.5, ]
  expect_equal(r$value, x[3, 3, 4])
  expect_true(all(is.na(res$index1[res$variable == "a"])))
})

test_that("can select variables, name dims and label them", {
  x <- array(seq_len(5 * 2), c(5, 2))
  nms <- c("a", "m[1,1]", "m[2,1]", "m[1,2]", "m[2,2]")
  res <- melt_state(x, nms, list(time = 1:2), vars = "m",
                    dims = c("age", "vax"),
                    labels = list(age = c("child", "adult")))
  expect_named(res, c("time", "variable", "age", "vax", "value"))
  expect_equal(levels(res$age), c("child", "adult"))
  expect_equal(res$value[res$age == "adult" & res$vax == 2], c(5, 10))
  expect_error(melt_state(x, nms, list(time = 1:2), vars = "z"), "Unknown")
  expect_error(melt_state(x, nms, list(time = 1:2), dims = list(m = "age")),
               "has length 1")
  expect_error(
    melt_state(x, nms, list(time = 1:2), dims = "age"), "names 1 dimension")
  expect_error(
    melt_state(x, nms, list(time = 1:2), dims = c("age", "vax"),
               labels = list(age = "child")),
    "has 1 element")
})

test_that("melt validates dimensions", {
  expect_error(melt_state(array(0, c(2, 3)), c("a", "b"), list(time = 1:4)),
               "Unexpected array dimensions")
})

test_that("variables are factors in the order of the state", {
  x <- array(0, c(4, 2))
  res <- melt_state(x, c("S", "I[1]", "I[2]", "A"), list(time = 1:2))
  expect_s3_class(res$variable, "factor")
  expect_false(is.ordered(res$variable))
  expect_equal(levels(res$variable), c("S", "I", "A"))
  res <- melt_state(x, c("S", "I[1]", "I[2]", "A"), list(time = 1:2),
                    vars = c("A", "S"))
  expect_equal(levels(res$variable), c("S", "A"))
})

test_that("dims must be a character vector or a named list", {
  x <- array(0, c(2, 2))
  expect_error(melt_state(x, c("a[1]", "a[2]"), list(time = 1:2), dims = 1),
               "must be an unnamed character vector or a named list")
  expect_error(
    melt_state(x, c("a[1]", "a[2]"), list(time = 1:2), dims = c(a = "age")),
    "must be an unnamed character vector or a named list")
})
