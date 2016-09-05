ode_interpolant(Θ,Δt,y₀,y₁,k₀,k₁,alg::Symbol) = ode_interpolant(Θ,Δt,y₀,y₁,k₀,k₁,Val{alg}) # Dispatch interpolant by alg
ode_interpolation(tvals,ts,timeseries,ks,alg::Symbol,f) = ode_interpolation(tvals,ts,timeseries,ks,Val{alg},f) # Dispatch by alg
ode_addsteps!(k,t,u,Δt,alg::Symbol,f) = ode_addsteps!(k,t,u,Δt,f,Val{alg},typeof(k[1]./Δt),eltype(k[1]./k[1]))

"""
ode_interpolation(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function ode_interpolation{alg}(tvals,ts,timeseries,ks,T::Type{Val{alg}},f)
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
    elseif ts[i-1] == t # Can happen if it's the first value!
      push!(vals,timeseries[i-1])
    else
      Δt = ts[i] - ts[i-1]
      Θ = (t-ts[i-1])/Δt
      ode_addsteps!(ks[i],ts[i-1],timeseries[i-1],Δt,alg,f) # update the kcurrent, since kprevious not used in special algs
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
function ode_interpolation{alg}(tval::Number,ts,timeseries,ks,T::Type{Val{alg}},f)
  i = findfirst((x)->x>=tval,ts) # It's in the interval ts[i-1] to ts[i]
  if ts[i] == tval
    val = timeseries[i]
  elseif ts[i-1] == tval # Can happen if it's the first value!
    push!(vals,timeseries[i-1])
  else
    Δt = ts[i] - ts[i-1]
    Θ = (tval-ts[i-1])/Δt
    ode_addsteps!(ks[i],ts[i-1],timeseries[i-1],Δt,alg,f) # update the kcurrent, since kprevious not used in special algs
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
By default, never add steps (no op)
"""
function ode_addsteps!{rateType,uEltypeNoUnits,alg}(k,t,u,Δt,f,T::Type{Val{alg}},T2::Type{rateType},T3::Type{uEltypeNoUnits})
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

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190
"""
function ode_interpolant(Θ,Δt,y₀,y₁,kprevious,k,T::Type{Val{:DP5Vectorized}})
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

"""
Runge–Kutta pairs of order 5(4) satisfying only the first column
simplifying assumption

Ch. Tsitouras
"""
function ode_interpolant(Θ,Δt,y₀,y₁,kprevious,k,T::Type{Val{:Tsit5}})
  b1Θ = -1.0530884977290216Θ * (Θ - 1.3299890189751412)*(Θ^2 - 1.4364028541716351Θ + 0.7139816917074209)
  b2Θ = 0.1017Θ^2 * (Θ^2 - 2.1966568338249754Θ + 1.2949852507374631)
  b3Θ = 2.490627285651252793Θ^2 * (Θ^2 - 2.38535645472061657Θ + 1.57803468208092486)
  b4Θ = -16.54810288924490272*(Θ - 1.21712927295533244)*(Θ - 0.61620406037800089)*Θ^2
  b5Θ = 47.37952196281928122*(Θ - 1.203071208372362603)*(Θ - 0.658047292653547382)*Θ^2
  b6Θ = -34.87065786149660974*(Θ - 1.2)*(Θ - 0.666666666666666667)*Θ^2
  b7Θ = 2.5*(Θ - 1)*(Θ - 0.6)*Θ^2
  y₀ + (k[1]*b1Θ + k[2]*b2Θ + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ)
end

"""
Runge–Kutta pairs of order 5(4) satisfying only the first column
simplifying assumption

Ch. Tsitouras
"""
function ode_interpolant(Θ,Δt,y₀,y₁,kprevious,k,T::Type{Val{:Tsit5Vectorized}})
  b1Θ = -1.0530884977290216Θ * (Θ - 1.3299890189751412)*(Θ^2 - 1.4364028541716351Θ + 0.7139816917074209)
  b2Θ = 0.1017Θ^2 * (Θ^2 - 2.1966568338249754Θ + 1.2949852507374631)
  b3Θ = 2.490627285651252793Θ^2 * (Θ^2 - 2.38535645472061657Θ + 1.57803468208092486)
  b4Θ = -16.54810288924490272*(Θ - 1.21712927295533244)*(Θ - 0.61620406037800089)*Θ^2
  b5Θ = 47.37952196281928122*(Θ - 1.203071208372362603)*(Θ - 0.658047292653547382)*Θ^2
  b6Θ = -34.87065786149660974*(Θ - 1.2)*(Θ - 0.666666666666666667)*Θ^2
  b7Θ = 2.5*(Θ - 1)*(Θ - 0.6)*Θ^2
  y₀ + (k[1]*b1Θ + k[2]*b2Θ + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ)
end

"""
Coefficients taken from RKSuite
"""
function ode_interpolant(Θ,Δt,y₀,y₁,kprevious,k,T::Type{Val{:BS5}})
  r016,r015,r014,r013,r012,r011,r036,r035,r034,r033,r032,r046,r045,r044,r043,r042,r056,r055,r054,r053,r052,r066,r065,r064,r063,r062,r076,r075,r074,r073,r072,r086,r085,r084,r083,r082,r096,r095,r094,r093,r092,r106,r105,r104,r103,r102,r116,r115,r114,r113,r112 = BS5Interp_polyweights(eltype(y₀))
  Θ2 = Θ^2
  Θ3 = Θ2*Θ
  Θ4 = Θ3*Θ
  Θ5 = Θ4*Θ
  Θ6 = Θ5*Θ
  b1Θ = r011*Θ  + r012*Θ2 + r013*Θ3 + r014*Θ4 + r015*Θ5 + r016*Θ6
  b3Θ =           r032*Θ2 + r033*Θ3 + r034*Θ4 + r035*Θ5 + r036*Θ6
  b4Θ =           r042*Θ2 + r043*Θ3 + r044*Θ4 + r045*Θ5 + r046*Θ6
  b5Θ =           r052*Θ2 + r053*Θ3 + r054*Θ4 + r055*Θ5 + r056*Θ6
  b6Θ =           r062*Θ2 + r063*Θ3 + r064*Θ4 + r065*Θ5 + r066*Θ6
  b7Θ =           r072*Θ2 + r073*Θ3 + r074*Θ4 + r075*Θ5 + r076*Θ6
  b8Θ =           r082*Θ2 + r083*Θ3 + r084*Θ4 + r085*Θ5 + r086*Θ6
  b9Θ =           r092*Θ2 + r093*Θ3 + r094*Θ4 + r095*Θ5 + r096*Θ6
  b10Θ=           r102*Θ2 + r103*Θ3 + r104*Θ4 + r105*Θ5 + r106*Θ6
  b11Θ=           r112*Θ2 + r113*Θ3 + r114*Θ4 + r115*Θ5 + r116*Θ6
  y₀ + (k[1]*b1Θ + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ)
end

"""
Coefficients taken from RKSuite
"""
function ode_interpolant(Θ,Δt,y₀,y₁,kprevious,k,T::Type{Val{:BS5Vectorized}})
  r016,r015,r014,r013,r012,r011,r036,r035,r034,r033,r032,r046,r045,r044,r043,r042,r056,r055,r054,r053,r052,r066,r065,r064,r063,r062,r076,r075,r074,r073,r072,r086,r085,r084,r083,r082,r096,r095,r094,r093,r092,r106,r105,r104,r103,r102,r116,r115,r114,r113,r112 = BS5Interp_polyweights(eltype(y₀))
  Θ2 = Θ^2
  Θ3 = Θ2*Θ
  Θ4 = Θ3*Θ
  Θ5 = Θ4*Θ
  Θ6 = Θ5*Θ
  b1Θ = r011*Θ  + r012*Θ2 + r013*Θ3 + r014*Θ4 + r015*Θ5 + r016*Θ6
  b3Θ =           r032*Θ2 + r033*Θ3 + r034*Θ4 + r035*Θ5 + r036*Θ6
  b4Θ =           r042*Θ2 + r043*Θ3 + r044*Θ4 + r045*Θ5 + r046*Θ6
  b5Θ =           r052*Θ2 + r053*Θ3 + r054*Θ4 + r055*Θ5 + r056*Θ6
  b6Θ =           r062*Θ2 + r063*Θ3 + r064*Θ4 + r065*Θ5 + r066*Θ6
  b7Θ =           r072*Θ2 + r073*Θ3 + r074*Θ4 + r075*Θ5 + r076*Θ6
  b8Θ =           r082*Θ2 + r083*Θ3 + r084*Θ4 + r085*Θ5 + r086*Θ6
  b9Θ =           r092*Θ2 + r093*Θ3 + r094*Θ4 + r095*Θ5 + r096*Θ6
  b10Θ=           r102*Θ2 + r103*Θ3 + r104*Θ4 + r105*Θ5 + r106*Θ6
  b11Θ=           r112*Θ2 + r113*Θ3 + r114*Θ4 + r115*Θ5 + r116*Θ6
  y₀ + (k[1]*b1Θ + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ)
end

"""

"""
function ode_interpolant(Θ,Δt,y₀,y₁,kprevious,k,T::Type{Val{:Vern6}})
  r011,r012,r013,r014,r015,r016,r042,r043,r044,r045,r046,r052,r053,r054,r055,r056,r062,r063,r064,r065,r066,r072,r073,r074,r075,r076,r082,r083,r084,r085,r086,r092,r093,r094,r095,r096,r102,r103,r104,r105,r106,r112,r113,r114,r115,r116,r122,r123,r124,r125,r126 = Vern6Interp_polyweights(eltype(y₀))
  Θ2 = Θ^2
  Θ3 = Θ2*Θ
  Θ4 = Θ3*Θ
  Θ5 = Θ4*Θ
  Θ6 = Θ5*Θ
  b1Θ = r011*Θ + r012*Θ2 + r013*Θ3 + r014*Θ4 + r015*Θ5 + r016*Θ6
  b4Θ =          r042*Θ2 + r043*Θ3 + r044*Θ4 + r045*Θ5 + r046*Θ6
  b5Θ =          r052*Θ2 + r053*Θ3 + r054*Θ4 + r055*Θ5 + r056*Θ6
  b6Θ =          r062*Θ2 + r063*Θ3 + r064*Θ4 + r065*Θ5 + r066*Θ6
  b7Θ =          r072*Θ2 + r073*Θ3 + r074*Θ4 + r075*Θ5 + r076*Θ6
  b8Θ =          r082*Θ2 + r083*Θ3 + r084*Θ4 + r085*Θ5 + r086*Θ6
  b9Θ =          r092*Θ2 + r093*Θ3 + r094*Θ4 + r095*Θ5 + r096*Θ6
  b10Θ=          r102*Θ2 + r103*Θ3 + r104*Θ4 + r105*Θ5 + r106*Θ6
  b11Θ=          r112*Θ2 + r113*Θ3 + r114*Θ4 + r115*Θ5 + r116*Θ6
  b12Θ=          r122*Θ2 + r123*Θ3 + r124*Θ4 + r125*Θ5 + r126*Θ6
  y₀ + k[1]*b1Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ + k[12]*b12Θ
end



"""
An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28

Called to add the extra k9, k10, k11 steps for the Order 5 interpolation when needed
"""
function ode_addsteps!{rateType<:Number,uEltypeNoUnits}(k,t,u,Δt,f,T::Type{Val{:BS5}},T2::Type{rateType},T3::Type{uEltypeNoUnits})
  if length(k) < 11 # Have not added the extra stages yet
    c6,c7,c8,a91,a92,a93,a94,a95,a96,a97,a98,a101,a102,a103,a104,a105,a106,a107,a108,a109,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110 = BS5Interp(uEltypeNoUnits)
    push!(k,Δt*f(t+c6*Δt,u+a91*k[1]+a92*k[2]+a93*k[3]+a94*k[4]+a95*k[5]+a96*k[6]+a97*k[7]+a98*k[8]))
    push!(k,Δt*f(t+c7*Δt,u+a101*k[1]+a102*k[2]+a103*k[3]+a104*k[4]+a105*k[5]+a106*k[6]+a107*k[7]+a108*k[8]+a109*k[9]))
    push!(k,Δt*f(t+c8*Δt,u+a111*k[1]+a112*k[2]+a113*k[3]+a114*k[4]+a115*k[5]+a116*k[6]+a117*k[7]+a118*k[8]+a119*k[9]+a1110*k[10]))
  end
  nothing
end

"""
An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28

Called to add the extra k9, k10, k11 steps for the Order 5 interpolation when needed
"""
function ode_addsteps!{rateType<:AbstractArray,uEltypeNoUnits}(k,t,u,Δt,f,T::Type{Val{:BS5}},T2::Type{rateType},T3::Type{uEltypeNoUnits})
  if length(k) < 11 # Have not added the extra stages yet
    c6,c7,c8,a91,a92,a93,a94,a95,a96,a97,a98,a101,a102,a103,a104,a105,a106,a107,a108,a109,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110 = BS5Interp(uEltypeNoUnits)
    rtmp = rateType(size(k[1]))
    f(t+c6*Δt,u+a91*k[1]+a92*k[2]+a93*k[3]+a94*k[4]+a95*k[5]+a96*k[6]+a97*k[7]+a98*k[8],rtmp); push!(k,Δt*rtmp)
    f(t+c7*Δt,u+a101*k[1]+a102*k[2]+a103*k[3]+a104*k[4]+a105*k[5]+a106*k[6]+a107*k[7]+a108*k[8]+a109*k[9],rtmp); push!(k,Δt*rtmp)
    f(t+c8*Δt,u+a111*k[1]+a112*k[2]+a113*k[3]+a114*k[4]+a115*k[5]+a116*k[6]+a117*k[7]+a118*k[8]+a119*k[9]+a1110*k[10],rtmp); push!(k,Δt*rtmp)
  end
  nothing
end

"""
An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28

Called to add the extra k9, k10, k11 steps for the Order 5 interpolation when needed
"""
function ode_addsteps!{rateType<:AbstractArray,uEltypeNoUnits}(k,t,u,Δt,f,T::Type{Val{:BS5Vectorized}},T2::Type{rateType},T3::Type{uEltypeNoUnits})
  if length(k) < 11 # Have not added the extra stages yet
    c6,c7,c8,a91,a92,a93,a94,a95,a96,a97,a98,a101,a102,a103,a104,a105,a106,a107,a108,a109,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110 = BS5Interp(uEltypeNoUnits)
    rtmp = rateType(size(k[1]))
    f(t+c6*Δt,u+a91*k[1]+a92*k[2]+a93*k[3]+a94*k[4]+a95*k[5]+a96*k[6]+a97*k[7]+a98*k[8],rtmp); push!(k,Δt*rtmp)
    f(t+c7*Δt,u+a101*k[1]+a102*k[2]+a103*k[3]+a104*k[4]+a105*k[5]+a106*k[6]+a107*k[7]+a108*k[8]+a109*k[9],rtmp); push!(k,Δt*rtmp)
    f(t+c8*Δt,u+a111*k[1]+a112*k[2]+a113*k[3]+a114*k[4]+a115*k[5]+a116*k[6]+a117*k[7]+a118*k[8]+a119*k[9]+a1110*k[10],rtmp); push!(k,Δt*rtmp)
  end
  nothing
end

"""

"""
function ode_addsteps!{rateType<:Number,uEltypeNoUnits}(k,t,u,Δt,f,T::Type{Val{:Vern6}},T2::Type{rateType},T3::Type{uEltypeNoUnits})
  if length(k) < 12 # Have not added the extra stages yet
    c10,a1001,a1004,a1005,a1006,a1007,a1008,a1009,c11,a1101,a1102,a1103,a1104,a1105,a1106,a1107,a1108,a1109,a1110,c12,a1201,a1202,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = Vern6Interp(uEltypeNoUnits)
    push!(k,Δt*f(t+c10*Δt,u+a1001*k[1]+a1004*k[4]+a1005*k[5]+a1006*k[6]+a1007*k[7]+a1008*k[8]+a1009*k[9]))
    push!(k,Δt*f(t+c11*Δt,u+a1101*k[1]+a1102*k[2]+a1103*k[3]+a1104*k[4]+a1105*k[5]+a1106*k[6]+a1107*k[7]+a1108*k[8]+a1109*k[9]+a1110*k[10]))
    push!(k,Δt*f(t+c12*Δt,u+a1201*k[1]+a1202*k[2]+a1203*k[3]+a1204*k[4]+a1205*k[5]+a1206*k[6]+a1207*k[7]+a1208*k[8]+a1209*k[9]+a1210*k[10]+a1211*k[11]))
  end
  nothing
end

"""

"""
function ode_addsteps!{rateType<:AbstractArray,uEltypeNoUnits}(k,t,u,Δt,f,T::Type{Val{:Vern6}},T2::Type{rateType},T3::Type{uEltypeNoUnits})
  if length(k) < 12 # Have not added the extra stages yet
    c10,a1001,a1004,a1005,a1006,a1007,a1008,a1009,c11,a1101,a1102,a1103,a1104,a1105,a1106,a1107,a1108,a1109,a1110,c12,a1201,a1202,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = Vern6Interp(uEltypeNoUnits)
    rtmp = rateType(size(k[1]))
    f(t+c10*Δt,u+a1001*k[1]+a1004*k[4]+a1005*k[5]+a1006*k[6]+a1007*k[7]+a1008*k[8]+a1009*k[9],rtmp); push!(k,Δt*rtmp)
    f(t+c11*Δt,u+a1101*k[1]+a1102*k[2]+a1103*k[3]+a1104*k[4]+a1105*k[5]+a1106*k[6]+a1107*k[7]+a1108*k[8]+a1109*k[9]+a1110*k[10],rtmp); push!(k,Δt*rtmp)
    f(t+c12*Δt,u+a1201*k[1]+a1202*k[2]+a1203*k[3]+a1204*k[4]+a1205*k[5]+a1206*k[6]+a1207*k[7]+a1208*k[8]+a1209*k[9]+a1210*k[10]+a1211*k[11],rtmp); push!(k,Δt*rtmp)
  end
  nothing
end

"""

"""
function ode_addsteps!{rateType<:AbstractArray,uEltypeNoUnits}(k,t,u,Δt,f,T::Type{Val{:Vern6Vectorized}},T2::Type{rateType},T3::Type{uEltypeNoUnits})
  if length(k) < 12 # Have not added the extra stages yet
    c10,a1001,a1004,a1005,a1006,a1007,a1008,a1009,c11,a1101,a1102,a1103,a1104,a1105,a1106,a1107,a1108,a1109,a1110,c12,a1201,a1202,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = Vern6Interp(uEltypeNoUnits)
    rtmp = rateType(size(k[1]))
    f(t+c10*Δt,u+a1001*k[1]+a1004*k[4]+a1005*k[5]+a1006*k[6]+a1007*k[7]+a1008*k[8]+a1009*k[9],rtmp); push!(k,Δt*rtmp)
    f(t+c11*Δt,u+a1101*k[1]+a1102*k[2]+a1103*k[3]+a1104*k[4]+a1105*k[5]+a1106*k[6]+a1107*k[7]+a1108*k[8]+a1109*k[9]+a1110*k[10],rtmp); push!(k,Δt*rtmp)
    f(t+c12*Δt,u+a1201*k[1]+a1202*k[2]+a1203*k[3]+a1204*k[4]+a1205*k[5]+a1206*k[6]+a1207*k[7]+a1208*k[8]+a1209*k[9]+a1210*k[10]+a1211*k[11],rtmp); push!(k,Δt*rtmp)
  end
  nothing
end
