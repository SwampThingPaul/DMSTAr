test_that("dmsta_node_route implements node routing priority cases", {

  # Common inputs
  Qi <- 10
  RecycleQ <- 2
  qtot <- Qi + RecycleQ

  # ---- 1) Seep-all (Seepout_Rate > 0 takes top priority) ----
  out1 <- dmsta_node_route(
    Qi = Qi,
    RecycleQ = RecycleQ,
    Seepout_Rate = 0.1,
    Qimax = 999,
    Qomax = 999
  )

  expect_equal(out1$Qout, 0)
  expect_equal(out1$SeepOut, qtot)
  expect_equal(out1$Bypass, 0)
  expect_equal(out1$Etest, 0)
  expect_equal(out1$SeepIn, 0)

  # ---- 2) Bypass-all via Qimax (>0) when Seepout_Rate <= 0 ----
  out2 <- dmsta_node_route(
    Qi = Qi,
    RecycleQ = RecycleQ,
    Seepout_Rate = 0,
    Qimax = 1,
    Qomax = 999
  )

  expect_equal(out2$Qout, 0)
  expect_equal(out2$SeepOut, 0)
  expect_equal(out2$Bypass, qtot)
  expect_equal(out2$Etest, 0)
  expect_equal(out2$SeepIn, 0)

  # ---- 3) Low-flow bypass via negative Qomax (uses min(Qi, -Qomax)) ----
  # Note: bypass uses Qi ONLY (not Qi + RecycleQ).
  # With Qomax = -3: bypass = min(10, 3) = 3; outflow = 12 - 3 = 9
  out3 <- dmsta_node_route(
    Qi = Qi,
    RecycleQ = RecycleQ,
    Seepout_Rate = 0,
    Qimax = 0,
    Qomax = -3
  )

  expect_equal(out3$Bypass, 3)
  expect_equal(out3$Qout, qtot - 3)
  expect_equal(out3$SeepOut, 0)
  expect_equal(out3$Etest, 0)
  expect_equal(out3$SeepIn, 0)

  # ---- 4) Outflow cap via positive Qomax (bypass remainder) ----
  # With Qomax = 7: bypass = 12 - 7 = 5; outflow = 7
  out4 <- dmsta_node_route(
    Qi = Qi,
    RecycleQ = RecycleQ,
    Seepout_Rate = 0,
    Qimax = 0,
    Qomax = 7
  )

  expect_equal(out4$Qout, 7)
  expect_equal(out4$Bypass, qtot - 7)
  expect_equal(out4$SeepOut, 0)
  expect_equal(out4$Etest, 0)
  expect_equal(out4$SeepIn, 0)

  # ---- sanity: outputs should not be negative ----
  for (x in list(out1, out2, out3, out4)) {
    expect_gte(x$Qout, 0)
    expect_gte(x$Bypass, 0)
    expect_gte(x$SeepOut, 0)
    expect_equal(x$SeepIn, 0)
    expect_equal(x$Etest, 0)
  }
})
