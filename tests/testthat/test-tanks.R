test_that("dmsta_build_tanks partitions area and snaps cumulative fraction", {
  tanks <- dmsta_build_tanks(A_cell = 10, ttankS = 3.5, snap_last = TRUE)

  expect_type(tanks, "list")
  expect_true(tanks$Ntanks == 4L)

  # Areas add up
  expect_equal(sum(tanks$A_Tank), 10, tolerance = 1e-12)

  # Fractions sum to 1; cumulative ends at 1
  expect_equal(sum(tanks$F_Tank), 1, tolerance = 1e-12)
  expect_equal(tanks$Fcum[tanks$Ntanks], 1, tolerance = 1e-12)

  # Non-negative
  expect_true(all(tanks$A_Tank >= 0))
  expect_true(all(tanks$F_Tank >= 0))
})

test_that("dmsta_p_init_state returns vectors aligned with Ntanks", {
  tanks <- dmsta_build_tanks(A_cell = 10, ttankS = 3.5)
  st <- dmsta_p_init_state(tanks, Z_init_m = 1, C_init_ppb = 100, Y_init_mgm2 = 10)

  expect_named(st, c("M", "S"))
  expect_length(st$M, tanks$Ntanks)
  expect_length(st$S, tanks$Ntanks)

  # Sanity: positive masses when inputs positive
  expect_true(all(st$M > 0))
  expect_true(all(st$S > 0))
})
