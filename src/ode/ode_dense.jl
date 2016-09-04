ode_interpolant(Θ,Δt,y₀,y₁,k₀,k₁,alg::Symbol) = ode_interpolant(Θ,Δt,y₀,y₁,k₀,k₁,Val{alg}) # Dispatch interpolant by alg
ode_interpolation(tvals,ts,timeseries,ks,alg::Symbol) = ode_interpolation(tvals,ts,timeseries,ks,Val{alg}) # Dispatch by alg

"""
ode_interpolation(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function ode_interpolation{alg}(tvals,ts,timeseries,ks,T::Type{Val{alg}})
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
      push!(vals,ode_interpolant(Θ,Δt,timeseries[i-1],timeseries[i],ks[i-1],ks[i],alg))
    end
  end
  vals
end

"""
ode_interpolation(tval::Number,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function ode_interpolation{alg}(tval::Number,ts,timeseries,ks,T::Type{Val{alg}})
  i = findfirst((x)->x>=tval,ts) # It's in the interval ts[i-1] to ts[i]
  if ts[i] == tval
    val = timeseries[i]
  else
    Δt = ts[i] - ts[i-1]
    Θ = (tval-ts[i-1])/Δt
    val = ode_interpolant(Θ,Δt,timeseries[i-1],timeseries[i],ks[i-1],ks[i],alg)
  end
  val
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190

Herimte Interpolation, chosen if no other dispatch for ode_interpolant
"""
function ode_interpolant{alg}(Θ,Δt,y₀,y₁,k₀,k₁,T::Type{Val{alg}}) # Default interpolant is Hermite
  (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*Δt*k₀ + Θ*Δt*k₁)
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190
"""
function ode_interpolant(Θ,Δt,y₀,y₁,kprevious,k,T::Type{Val{:DP5}})
  b1,b3,b4,b5,b6,b7 = DP5_dense_bs(eltype(y₀/y₀)) # Divide away the units
  b7Θ = Θ^2 * (Θ-1) + Θ^2 * (Θ-1)^2 *10*(7414447 - 829305Θ)/29380423
  b1Θ = Θ^2 * (3-2Θ)*b1 + Θ*(Θ-1)^2 - Θ^2*(Θ-1)^2 *5*(2558722523 - 31403016Θ)/11282082432
  b3Θ = Θ^2 * (3-2Θ)*b3 + Θ^2 * (Θ-1)^2 * 100   * (882725551 - 15701508Θ)/32700410799
  b4Θ = Θ^2 * (3-2Θ)*b4 - Θ^2 * (Θ-1)^2 * 25    * (443332067 - 31403016Θ)/1880347072
  b5Θ = Θ^2 * (3-2Θ)*b5 + Θ^2 * (Θ-1)^2 * 32805 * (23143187  - 3489224Θ )/199316789632
  b6Θ = Θ^2 * (3-2Θ)*b6 - Θ^2 * (Θ-1)^2 * 55    * (29972135  - 7076736Θ )/822651844
  y₀ + (k[1]*b1Θ + k[2]*b3Θ + k[3]*b4Θ + k[4]*b5Θ + k[5]*b6Θ + k[6]*b7Θ) # No k2
end

function DP5_dense_bs(T)
  b1  = T(5179//57600)
  b3  = T(7571//16695)
  b4  = T(393//640)
  b5  = T(-92097//339200)
  b6  = T(187//2100)
  b7  = T(1//40)
  return b1,b3,b4,b5,b6,b7
end
