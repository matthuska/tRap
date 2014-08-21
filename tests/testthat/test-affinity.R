context("affinity")

test_that("affinity runs example from its man page", {
  pwm <- matrix(c(5, 4, 3, 1, 10, 12, 5, 3, 3, 5, 3, 10), nrow=4)
  seq <- "ACTGACGTGTGCACACGATGCTAGCTG"
  affinity <- affinity(pwm, seq)
  expect_true(affinity > 0)
})
