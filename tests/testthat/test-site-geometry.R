# Behaviour spec for the shared site geometry helpers (C1, issue #14):
# .bearing_distance(), .relative_elevation(), .role_features(),
# .crosswind_halfwidth(). All computed in a projected metric CRS.
# Pending specs (skipped) are the checklist /implement turns green.

describe(".bearing_distance()", {
  it("returns a bearing in degrees from north (clockwise) and a distance in metres")
  it("gives bearing 90 and the coordinate separation for a receptor due east")
  it("gives bearing 0 due north and 180 due south")
  it("is symmetric in distance between two features")
  it("gives a reverse bearing that differs by 180 degrees")
  it("returns distance 0 and bearing NA for coincident points")
})

describe(".relative_elevation()", {
  it("returns receptor elevation minus the source base on the AGL datum")
  it("propagates NA when an elevation is missing")
})

describe(".role_features()", {
  it("returns exactly the features whose roles match the hazard and role")
  it("returns a zero-row sf for an absent hazard or role")
})

describe(".crosswind_halfwidth()", {
  it("returns half the footprint extent perpendicular to the wind")
  it("gives half the side for a square aligned N-S with wind from the north")
  it("gives half of side*sqrt(2) for that square with wind from the north-east")
  it("is non-negative")
  it("is symmetric under a 180-degree wind reversal")
  it("returns 0 for a degenerate (point) geometry")
})
