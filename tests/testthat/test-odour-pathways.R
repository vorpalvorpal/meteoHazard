# Behaviour spec for the odour terrain morning-pulse physics (C3b, issue #17):
# .pool_partition(), .cw_venting() [1a], .cw_fumigation() [1b], .morning_release().
# Pending specs (skipped) are the checklist /implement turns green.

describe(".pool_partition()", {
  it("returns within-pool (1a) and above-pool (1b) fractions that sum to 1")
  it("increases the above-pool fraction monotonically with (emit-height - pool_top)")
  it("sends ~all emission to 1a for a low source deep in the pool")
  it("sends ~all emission to 1b for a high source above a shallow pool")
  it("varies continuously as pool_top sweeps an E(z) band edge")
})

describe(".cw_venting() [pathway 1a]", {
  it("confines at night, lowering an up-slope receptor relative to flat terrain")
  it("vents up-slope at the morning transition, raising an aligned receptor")
  it("uses radial drainage for a lone mound and channelled drainage along the drainage bearing")
})

describe(".cw_fumigation() [pathway 1b]", {
  it("directs the morning burst downwind along the residual-layer wind, not the drainage axis")
  it("barely raises a receptor up-drainage but crosswind to the residual wind")
})

describe("both pathways over a tall mound", {
  it("contributes to both an up-slope and a downwind receptor in the same morning")
})

describe(".morning_release()", {
  it("conserves the accumulated mass A over the release window")
  it("gives a shorter, taller pulse with the same integral when CBL growth doubles")
  it("saturates the accumulation A with cooling hours")
  it("is idempotent: the same input window yields the same result")
  it("produces no pulse (and no error) for a window with no preceding night")
  it("draws 1a confinement, 1a venting and 1b fumigation from a single budget A")
})

describe("odour_exposure() terrain-backend edges", {
  it("falls back to surface advection when pool_top or residual_wind is NA")
  it("reproduces the C3a flat result with terrain_backend = 'none'")
})
