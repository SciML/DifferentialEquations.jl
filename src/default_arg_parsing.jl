const LOW_TOLERANCE = 1e-6
const HIGH_TOLERANCE = 1e-1

get_alg_hints(o) = :alg_hints ∈ keys(o) ? alg_hints = o[:alg_hints] : alg_hints = Symbol[:nonstiff]
