# litter wetness state: shortwave_radiation derived from direct + diffuse
# radiation when absent (or NA). Fixtures in helper-na-policy.R.

describe("shortwave_radiation derived from direct + diffuse", {
  it("litter_hazard(use_wetness_state = TRUE) derives it when absent", {
    met <- .na_met()
    with_sw <- met
    with_sw$shortwave_radiation <- met$direct_radiation + met$diffuse_radiation
    expect_message(
      derived <- litter_hazard(met, use_wetness_state = TRUE),
      class = "meteoHazard_derived_input"
    )
    expect_identical(derived, litter_hazard(with_sw, use_wetness_state = TRUE))
    expect_no_message(litter_hazard(met, use_wetness_state = TRUE, verbose = FALSE))
  })

  it("litter_wetness() and litter_risk() derive it too", {
    met <- .na_met()
    with_sw <- met
    with_sw$shortwave_radiation <- met$direct_radiation + met$diffuse_radiation
    expect_identical(
      suppressMessages(litter_wetness(met)), litter_wetness(with_sw)
    )
    expect_identical(
      suppressMessages(litter_risk(met, .na_litter_site(), use_wetness_state = TRUE)),
      litter_risk(with_sw, .na_litter_site(), use_wetness_state = TRUE)
    )
  })

  it("fills NA shortwave hours from direct + diffuse where available", {
    met <- .na_met()
    met$shortwave_radiation <- met$direct_radiation + met$diffuse_radiation
    ref <- litter_wetness(met)
    met$shortwave_radiation[3] <- NA
    expect_message(out <- litter_wetness(met), class = "meteoHazard_derived_input")
    expect_identical(out, ref)
  })

  it("still errors when neither shortwave nor direct + diffuse is present", {
    met <- .na_met()
    met$diffuse_radiation <- NULL
    expect_error(litter_hazard(met, use_wetness_state = TRUE),
                 class = "meteoHazard_input_error")
  })
})

