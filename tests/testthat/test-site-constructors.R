# Behaviour spec for the site convenience constructors (C4, issue #18):
# site_from_sectors() — compass-sector descriptions -> features geometry.
# Pending specs (skipped) are the checklist /implement turns green.

describe("site_from_sectors()", {
  it("places a barrier from arc [NE, SE] spanning bearings 45-135 from the centroid")
  it("carries permeability and sensitive onto the feature roles")
  it("wraps an arc that crosses north (e.g. [NW, NE]) correctly")
  it("requires an explicit origin (centroid / working face)")
})
