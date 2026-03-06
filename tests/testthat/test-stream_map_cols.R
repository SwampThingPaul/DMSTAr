test_that("stream_map_cols returns correct Q/L mapping", {
  expect_equal(stream_map_cols("bypass"),   c(Q="Q_out_bypass",         L="L_out_bypass"))
  expect_equal(stream_map_cols("release1"), c(Q="Q_out_release1",       L="L_out_release1"))
  expect_equal(stream_map_cols("release2"), c(Q="Q_out_release2",       L="L_out_release2"))
  expect_equal(stream_map_cols("seepage"),  c(Q="Q_out_seep_discharge", L="L_out_seep_discharge"))
})

test_that("stream_map_cols outflow treated vs total", {
  expect_equal(stream_map_cols("outflow","treated"), c(Q="Q_out_treated", L="L_out_treated"))
  expect_equal(stream_map_cols("outflow","total"),   c(Q="Q_out_total",   L="L_out_total"))
})

test_that("stream_map_cols errors on unknown stream", {
  expect_error(stream_map_cols("bogus"), "Unknown stream")
})
