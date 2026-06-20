# Behaviour spec for the geo-referenced site model: mh_site() and mh_terrain().
# Plan: GitHub issue #14 (C1). Pending specs (skipped) are the checklist that
# /implement turns green; concrete assertions, tolerances and seeds are filled
# in then. See the C1 plan and the Stage-1 terrain-schema spike findings (#14).

describe("mh_site()", {
  it("bundles features, roles and terrain into one validated object")
  it("requires every roles$feature_id to exist in features$id")
  it("requires a single CRS across all features")
  it("requires a projected metric CRS")
  it("reprojects a geographic (lon/lat) CRS to the supplied metric epsg")
  it("matches sf::st_transform for a known point under reprojection")
  it("enforces a shared AGL datum between source and terrain heights")
  it("allows a feature that carries no roles")
  it("allows a hazard with no source or no receptor at construction")
  it("errors (classed) when no metric epsg is available")
  it("errors (classed) on an unknown hazard or role value")
})

describe("mh_terrain()", {
  it("holds the scalar terrain descriptors with units and a _meta scale block")
  it("records relief as a height above local base in metres")          # spike: not standardised DEV
  it("carries _meta scale fields (relief_radius, channel_threshold, fetch_L)")
  it("rejects a negative relief, valley_depth or basin_capacity")
  it("rejects a topographic amplification factor below 1")
  it("rejects a drainage bearing outside [0, 360)")
  it("requires flow_convergence and shelter_index to be finite")
  it("permits NA descriptors as 'no terrain effect' (flat site)")
})

# The existing package units contract, exercised through the new constructors.
describe("mh_site() / mh_terrain() units handling", {
  it("accepts distances/elevations as bare numerics in the documented unit")
  it("accepts units-tagged distances/elevations and converts them")
  it("errors (classed) on a dimensionally incompatible unit")
})

# === C4 addition (issue #18): dust receptor roles on the unified model ===
describe("mh_site(): dust receptor roles", {
  it("validates a (dust, receptor) role")
  it("does not require a (dust, barrier) role")
})
