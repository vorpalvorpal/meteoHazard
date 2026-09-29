# Odour input range validation: invalid (non-missing) values error, as they do
# for dust, litter and TWL; NA stays a missing input (NA row + one warning).
# Fixtures in helper-na-policy.R.

odour_bad_values <- list(
  wind_speed_10m         = -1,
  wind_direction_10m     = 400,
  direct_radiation       = -50,
  cloud_cover            = 120,
  boundary_layer_height  = -10,
  relative_humidity_2m   = 150,
  precipitation          = -0.5,
  soil_moisture_0_to_1cm = 1.5,
  soil_moisture_1_to_3cm = -0.1,
  pressure_msl           = -1
)

describe("odour input ranges", {
  for (col in names(odour_bad_values)) {
    local({
      col <- col
      it(paste("odour_risk() and odour_hazard() reject an invalid", col), {
        met <- .na_met()
        met[[col]][2] <- odour_bad_values[[col]]
        expect_error(odour_risk(met, .na_odour_site()), col)
        if (col != "wind_direction_10m") {
          expect_error(odour_hazard(met), col)
        }
      })
    })
  }

  it("rejects an invalid optional multi-level wind speed", {
    met <- .na_met()
    met$wind_speed_80m <- 8
    met$wind_speed_80m[3] <- -2
    expect_error(odour_risk(met, .na_odour_site()), "wind_speed_80m")
  })

  it("accepts boundary values and still treats NA as missing", {
    met <- .na_met()
    met$relative_humidity_2m[1] <- 100
    met$cloud_cover[2] <- 0
    met$wind_direction_10m[3] <- 360
    met$soil_moisture_0_to_1cm[4] <- 1
    expect_no_error(odour_risk(met, .na_odour_site()))
    met$relative_humidity_2m[5] <- NA
    res <- .count_missing(odour_risk(met, .na_odour_site()))
    expect_equal(res$n, 1L)
    expect_true(all(is.na(res$value[5, ])))
  })
})
