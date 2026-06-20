# Behaviour spec for the odour ventilation state (C2, issue #15):
# ventilation_state() and its pool_top / cbl_growth / residual_wind / generation
# components. pool_top resolves D3 (terrain-modulated heat-deficit estimate).
# Pending specs (skipped) are the checklist /implement turns green.

describe("ventilation_state()", {
  it("returns per-hour u_eff, h_mix, s, is_calm, is_day, PM, W_rain, pool_top, cbl_growth and residual_wind")
  it("reproduces the current dispersion-state u_eff, h_mix and s exactly")
})

describe("ventilation_state(): pool_top (terrain-modulated heat deficit)", {
  it("computes net longwave cooling from the Brunt formula for given T, RH and cloud")
  it("gives stronger cooling under clear skies than overcast")
  it("increases monotonically with cooling hours until saturation, then plateaus")
  it("is floored by the Venkatram mechanical depth 2400*u*^1.5")
  it("matches the Venkatram floor for a given friction velocity")
  it("is capped by the basin sill (valley_depth / basin_capacity)")
  it("is amplified, not reduced, by a topographic amplification factor >= 1")
  it("falls back to the mechanical floor / h_mix outside the clear-calm-stable regime")
  it("yields no accumulated pool for an input window containing no night")
  it("is taken as the overnight maximum, frozen at release, on the AGL datum")
})

describe("ventilation_state(): cbl_growth", {
  it("is positive across the morning transition and ~0 at night")
  it("returns 0 (not an NA error) for a flat or missing h_mix series")
})

describe("ventilation_state(): residual_wind", {
  it("is the overnight circular-mean direction at each available level")
  it("falls back to the next lower level when an upper level is missing")
  it("returns NA when all winds are missing")
})

describe(".odour_generation()", {
  it("reproduces the current generation modifier G exactly")
})

describe("ventilation_state(): validation", {
  it("errors (classed) on missing required met columns, including the multi-level winds")
})
