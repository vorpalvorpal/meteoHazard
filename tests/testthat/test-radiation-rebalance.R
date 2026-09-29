# Negative radiation components: Open-Meteo's dawn hours can carry a negative
# diffuse_radiation beside a positive direct_radiation (real values below,
# Blaxland / Katoomba 2026-10-03 07:00, shortwave 23 and 32 W/m^2). The
# components are rebalanced so the global horizontal total direct + diffuse is
# kept, with one meteoHazard_input_adjusted message per call; only a clearly
# negative total (< -5 W/m^2) errors. Fixtures in helper-na-policy.R.

# Real Open-Meteo rows (direct, diffuse) and their rebalanced equivalents.
.dawn_met <- function() {
  met <- .na_met()
  met$direct_radiation[2:4]  <- c(85, 3, 43)
  met$diffuse_radiation[2:4] <- c(-62, -2, -11)
  met
}
.dawn_met_fixed <- function() {
  met <- .na_met()
  met$direct_radiation[2:4]  <- c(23, 1, 32)
  met$diffuse_radiation[2:4] <- c(0, 0, 0)
  met
}

.count_adjusted <- function(expr) {
  n <- 0L
  value <- withCallingHandlers(
    expr,
    meteoHazard_input_adjusted = function(m) {
      n <<- n + 1L
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, n = n)
}

describe("negative radiation rebalance", {
  twl <- function(met) {
    generate_twl(
      datetime = met$datetime, latitude = -33.7, longitude = 150.3,
      temp = met$temperature_2m, wind_speed = met$wind_speed_10m / 3,
      wind_height = 10, RH = met$relative_humidity_2m,
      direct_solar = met$direct_radiation, diffuse_solar = met$diffuse_radiation,
      pressure = met$surface_pressure, verbose = FALSE
    )
  }

  it("generate_twl() rebalances Open-Meteo's negative diffuse, keeping the total", {
    res <- .count_adjusted(twl(.dawn_met()))
    expect_equal(res$n, 1L)
    expect_identical(res$value, twl(.dawn_met_fixed()))
  })

  it("moves a small negative direct into diffuse, keeping the total", {
    met <- .na_met()
    met$direct_radiation[3]  <- -3
    met$diffuse_radiation[3] <- 40
    ref <- .na_met()
    ref$direct_radiation[3]  <- 0
    ref$diffuse_radiation[3] <- 37
    res <- .count_adjusted(twl(met))
    expect_equal(res$n, 1L)
    expect_identical(res$value, twl(ref))
  })

  it("clamps a total just below zero (within 5 W/m^2) to zero", {
    met <- .na_met()
    met$direct_radiation[3]  <- 0
    met$diffuse_radiation[3] <- -3
    ref <- .na_met()
    ref$direct_radiation[3]  <- 0
    ref$diffuse_radiation[3] <- 0
    res <- .count_adjusted(twl(met))
    expect_equal(res$n, 1L)
    expect_identical(res$value, twl(ref))
  })

  it("still errors, naming the row, when direct + diffuse is clearly negative", {
    met <- .na_met()
    met$direct_radiation[3]  <- 10
    met$diffuse_radiation[3] <- -20
    expect_error(twl(met), class = "meteoHazard_input_error")
    expect_error(twl(met), "row 3")
  })

  it("leaves non-negative radiation untouched, with no message", {
    expect_no_message(twl(.na_met()))
  })

  it("odour_risk() and odour_hazard() rebalance direct_radiation the same way", {
    site <- .na_odour_site()
    res <- .count_adjusted(odour_risk(.dawn_met(), site))
    expect_equal(res$n, 1L)
    expect_identical(res$value, odour_risk(.dawn_met_fixed(), site))
    res <- .count_adjusted(odour_hazard(.dawn_met()))
    expect_equal(res$n, 1L)
    expect_identical(res$value, odour_hazard(.dawn_met_fixed()))
  })

  it("odour without a diffuse column clamps a small negative direct only", {
    met <- .na_met()
    met$diffuse_radiation <- NULL
    ref <- met
    met$direct_radiation[3] <- -2
    ref$direct_radiation[3] <- 0
    res <- .count_adjusted(odour_hazard(met))
    expect_equal(res$n, 1L)
    expect_identical(res$value, odour_hazard(ref))
    met$direct_radiation[3] <- -50
    expect_error(odour_hazard(met), "direct_radiation")
  })

  it("litter's derived shortwave is the preserved total (23 W/m^2), no rebalance needed", {
    met <- .dawn_met()
    ref <- met
    ref$shortwave_radiation <- met$direct_radiation + met$diffuse_radiation
    expect_equal(ref$shortwave_radiation[2:4], c(23, 1, 32))
    out <- suppressMessages(litter_wetness(met))
    expect_identical(out, litter_wetness(ref))
    expect_identical(
      suppressMessages(litter_risk(met, .na_litter_site(), use_wetness_state = TRUE)),
      litter_risk(ref, .na_litter_site(), use_wetness_state = TRUE)
    )
  })

  it("litter's derived shortwave clamps a total just below zero and errors on a clearly negative one", {
    met <- .na_met()
    met$direct_radiation[3]  <- 0
    met$diffuse_radiation[3] <- -3
    ref <- met
    ref$shortwave_radiation <- met$direct_radiation + met$diffuse_radiation
    ref$shortwave_radiation[3] <- 0
    res <- .count_adjusted(suppressMessages(litter_wetness(met), classes = "meteoHazard_derived_input"))
    expect_equal(res$n, 1L)
    expect_identical(res$value, litter_wetness(ref))
    met$diffuse_radiation[3] <- -20
    expect_error(suppressMessages(litter_wetness(met)), class = "meteoHazard_input_error")
  })
})
