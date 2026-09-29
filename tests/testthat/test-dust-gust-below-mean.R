# dust: gust below the mean wind -- error naming rows/values by default, or
# clamp with gust_below_mean = "clamp". Fixtures in helper-na-policy.R.

describe("dust gust below mean wind", {
  it("errors by default, naming the rows and values", {
    met <- .na_met()
    met$wind_gusts_10m[c(2, 5)] <- met$wind_speed_10m[c(2, 5)] - 0.5
    err <- tryCatch(dust_hazard(met), error = function(e) e)
    expect_s3_class(err, "meteoHazard_input_error")
    msg <- conditionMessage(err)
    expect_match(msg, "row 2 (gust 8.5 < mean 9)", fixed = TRUE)
    expect_match(msg, "row 5 (gust 3.5 < mean 4)", fixed = TRUE)
  })

  it("gust_below_mean = 'clamp' raises the gust to the mean with a warning", {
    met <- .na_met()
    met$wind_gusts_10m[2] <- met$wind_speed_10m[2] - 0.5
    fixed <- met
    fixed$wind_gusts_10m[2] <- met$wind_speed_10m[2]
    expect_warning(
      out <- dust_hazard(met, gust_below_mean = "clamp"),
      class = "meteoHazard_gust_clamped"
    )
    expect_identical(out, dust_hazard(fixed))
    expect_warning(
      out_f <- dust_flux(20L, 10, c(9, 5), c(8.5, 6), c(0.02, 0.02),
                         gust_below_mean = "clamp"),
      class = "meteoHazard_gust_clamped"
    )
    expect_identical(out_f, dust_flux(20L, 10, c(9, 5), c(9, 6), c(0.02, 0.02)))
  })

  it("clamping also feeds the saltation crust gate", {
    met <- .na_met()
    met$wind_gusts_10m[4] <- met$wind_speed_10m[4] - 1
    fixed <- met
    fixed$wind_gusts_10m[4] <- met$wind_speed_10m[4]
    out <- suppressWarnings(dust_hazard(met, crust = TRUE, crust_decay = "saltation",
                                        gust_below_mean = "clamp"))
    expect_identical(out, dust_hazard(fixed, crust = TRUE, crust_decay = "saltation"))
  })
})
