get_alg_hints(o) = :alg_hints ∈ keys(o) ? alg_hints = o[:alg_hints] : alg_hints = Symbol[:nonstiff]

const LOW_TOL = 1e-6
const MED_TOL = 1e-2
const EXTREME_TOL = 1e-9

function get_tolerance_level(o)
  :reltol ∈ keys(o) ? reltol = o[:reltol] : reltol = 1e-3
  if reltol > MED_TOL
    return :high_tol
  elseif reltol > LOW_TOL
    return :med_tol
  elseif reltol > EXTREME_TOL
    return :low_tol
  else
    return :extreme_tol
  end
end
