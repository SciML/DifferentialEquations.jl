"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190
"""
function hermite3_interpolant(Θ,Δt,y₀,y₁,k₀,k₁)
  (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*Δt*k₀ + Θ*Δt*k₁)
end

"""
hermite3_interpolate(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function hermite3_interpolate(tvals::AbstractArray,ts,timeseries,ks)
  if typeof(tvals) <: StepRange
  else
    sort!(tvals) # Make sure it's sorted
  end
  i = 2 # Start the search thinking it's between ts[1] and ts[2]
  vals = Vector{eltype(timeseries)}(0)
  for t in tvals
    i = findfirst((x)->x>=t,ts[i:end])+i-1 # It's in the interval ts[i-1] to ts[i]
    if ts[i] == t
      push!(vals,timeseries[i])
    else
      Δt = ts[i] - ts[i-1]
      Θ = (t-ts[i-1])/Δt
      push!(vals,hermite3_interpolant(Θ,Δt,timeseries[i-1],timeseries[i],ks[i-1],ks[i]))
    end
  end
  vals
end

"""
hermite3_interpolate(tval::Number,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function hermite3_interpolate(tval::Number,ts,timeseries,ks)
  i = findfirst((x)->x>=tval,ts) # It's in the interval ts[i-1] to ts[i]
  if ts[i] == tval
    val = timeseries[i]
  else
    Δt = ts[i] - ts[i-1]
    Θ = (tval-ts[i-1])/Δt
    val = hermite3_interpolant(Θ,Δt,timeseries[i-1],timeseries[i],ks[i-1],ks[i])
  end
  val
end

function ode_interpolation(tvals,ts,timeseries,ks,alg) # Default to Hermite
  if alg ∉ DIFFERENTIALEQUATIONSJL_SPECIALDENSEALGS
    vals =  hermite3_interpolate(tvals,ts,timeseries,ks)
  end
  return vals
end
