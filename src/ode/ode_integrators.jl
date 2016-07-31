@def ode_loopheader begin
  iter += 1
  if iter > maxiters
    warn("Max Iters Reached. Aborting")
    # u = map((x)->oftype(x,NaN),u)
    break
  end
end

@def ode_savevalues begin
  if save_timeseries && iter%timeseries_steps==0
    push!(timeseries,u)
    push!(ts,t)
  end
end

@def ode_loopfooter begin
  if adaptive
    standard = abs(1/(γ*EEst))^(1/order)
    if isinf(standard)
        q = qmax
    else
       q = min(qmax,max(standard,eps()))
    end
    if q > 1
      t = t + Δt
      u = utmp
      @ode_savevalues
    end
    Δtpropose = min(Δtmax,q*Δt)
    Δt = max(min(Δtpropose,abs(T-t)),Δtmin) #abs to fix complex sqrt issue at end
  else #Not adaptive
    t = t + Δt
    @ode_savevalues
  end
  (progressbar && atomloaded && iter%progress_steps==0) ? Main.Atom.progress(t/T) : nothing #Use Atom's progressbar if loaded
end

function ode_euler(f::Function,u,t,Δt,T,iter,maxiters,
                    timeseries,ts,timeseries_steps,save_timeseries,adaptive,progressbar)
  while t < T
    @ode_loopheader
    u = u + Δt.*f(u,t)
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_midpoint(f::Function,u::Number,t,Δt,T,iter,
                      maxiters,timeseries,ts,timeseries_steps,save_timeseries,adaptive,progressbar)
  halfΔt = Δt/2
  while t < T
    @ode_loopheader
    u = u + Δt.*f(u+halfΔt.*f(u,t),t+halfΔt)
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_midpoint(f::Function,u::AbstractArray,t,Δt,T,iter,
                      maxiters,timeseries,ts,timeseries_steps,save_timeseries,adaptive,progressbar)
  halfΔt = Δt/2
  utilde = similar(u)
  while t < T
    @ode_loopheader
    utilde[:] = u+halfΔt.*f(u,t)
    u = u + Δt.*f(utilde,t+halfΔt)
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_rk4(f::Function,u::Number,t,Δt,T,iter,maxiters,
                timeseries,ts,timeseries_steps,save_timeseries,adaptive,progressbar)
  halfΔt = Δt/2
  while t < T
    @ode_loopheader
    k₁ = f(u,t)
    ttmp = t+halfΔt
    k₂ = f(u+halfΔt*k₁,ttmp)
    k₃ = f(u+halfΔt*k₂,ttmp)
    k₄ = f(u+Δt*k₃,t+Δt)
    u = u + Δt*(k₁ + 2k₂ + 2k₃ + k₄)/6
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_rk4(f::Function,u::AbstractArray,t,Δt,T,
                iter,maxiters,timeseries,ts,timeseries_steps,save_timeseries,adaptive,progressbar)
  halfΔt = Δt/2
  k₁ = similar(u)
  k₂ = similar(u)
  k₃ = similar(u)
  k₄ = similar(u)
  while t < T
    @ode_loopheader
    k₁[:] = f(u,t)
    ttmp = t+halfΔt
    k₂[:] = f(u+halfΔt*k₁,ttmp)
    k₃[:] = f(u+halfΔt*k₂,ttmp)
    k₄[:] = f(u+Δt*k₃,t+Δt)
    u = u + Δt*(k₁ + 2k₂ + 2k₃ + k₄)/6
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_explicitrk(f::Function,u::Number,t,Δt,T,iter,maxiters,timeseries,ts,
                        timeseries_steps,save_timeseries,A,c,α,αEEst,stages,order,γ,adaptive,
                        abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  ks = Array{typeof(u)}(stages)
  while t < T
    @ode_loopheader
    for i = 1:stages
      utilde = 0
      for j = 1:i-1
        utilde += A[i,j]*ks[j]
      end
      ks[i] = f(u+Δt*utilde,t+c[i]*Δt)
    end
    utilde = α[1]*ks[1]
    for i = 2:stages
      utilde += α[i]*ks[i]
    end
    if adaptive
      utmp = u + Δt*utilde
      uEEst = αEEst[1]*ks[1]
      for i = 2:stages
        uEEst += αEEst[i]*ks[i]
      end
      EEst = norm((utilde-uEEst)./(abstol+u*reltol),internalnorm)
    else
      u = u + Δt*utilde
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_explicitrk(f::Function,u::AbstractArray,t,Δt,T,iter,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        A,c,α,αEEst,stages,order,γ,adaptive,
                        abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  ks = Array{eltype(u)}(size(u)...,stages)
  while t < T
    @ode_loopheader
    for i = 1:stages
      utilde = zeros(u)
      for j = 1:i-1
        utilde += A[i,j]*ks[..,j]
      end
      ks[..,i] = f(u+Δt*utilde,t+c[i]*Δt)
    end
    utilde = α[1]*ks[..,1]
    for i = 2:stages
      utilde += α[i]*ks[..,i]
    end
    if adaptive
      utmp = u + Δt*utilde
      uEEst = αEEst[1]*ks[..,1]
      for i = 2:stages
        uEEst += αEEst[i]*ks[..,i]
      end
      EEst = norm((utilde-uEEst)./(abstol+u*reltol),internalnorm)
    else
      u = u + Δt*utilde
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_feagin10(f::Function,u::AbstractArray,t,Δt,T,iter,order,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b,c = constructFeagin10(eltype(u))
  k = Vector{typeof(u)}(17)
  sumIdx = [collect(1:3);5;7;collect(9:17)]
  if adaptive
    utmp = similar(u)
  end
  while t < T
    @ode_loopheader
    k[1]  = Δt*f(u,t)
    k[2]  = Δt*f(u + a0100*k[1],t + c[1]*Δt)
    k[3]  = Δt*f(u + a0200*k[1] + a0201*k[2],t + c[2]*Δt )
    k[4]  = Δt*f(u + a0300*k[1]              + a0302*k[3],t + c[3]*Δt)
    k[5]  = Δt*f(u + a0400*k[1]              + a0402*k[3] + a0403*k[4],t + c[4]*Δt)
    k[6]  = Δt*f(u + a0500*k[1]                           + a0503*k[4] + a0504*k[5],t + c[5]*Δt)
    k[7]  = Δt*f(u + a0600*k[1]                           + a0603*k[4] + a0604*k[5] + a0605*k[6],t + c[6]*Δt)
    k[8]  = Δt*f(u + a0700*k[1]                                        + a0704*k[5] + a0705*k[6] + a0706*k[7],t + c[7]*Δt)
    k[9]  = Δt*f(u + a0800*k[1]                                                     + a0805*k[6] + a0806*k[7] + a0807*k[8],t + c[8]*Δt)
    k[10] = Δt*f(u + a0900*k[1]                                                     + a0905*k[6] + a0906*k[7] + a0907*k[8] + a0908*k[9],t + c[9]*Δt)
    k[11] = Δt*f(u + a1000*k[1]                                                     + a1005*k[6] + a1006*k[7] + a1007*k[8] + a1008*k[9] + a1009*k[10],t + c[10]*Δt)
    k[12] = Δt*f(u + a1100*k[1]                                                     + a1105*k[6] + a1106*k[7] + a1107*k[8] + a1108*k[9] + a1109*k[10] + a1110*k[11],t + c[11]*Δt)
    k[13] = Δt*f(u + a1200*k[1]                           + a1203*k[4] + a1204*k[5] + a1205*k[6] + a1206*k[7] + a1207*k[8] + a1208*k[9] + a1209*k[10] + a1210*k[11] + a1211*k[12],t + c[12]*Δt)
    k[14] = Δt*f(u + a1300*k[1]              + a1302*k[3] + a1303*k[4]              + a1305*k[6] + a1306*k[7] + a1307*k[8] + a1308*k[9] + a1309*k[10] + a1310*k[11] + a1311*k[12] + a1312*k[13],t + c[13]*Δt)
    k[15] = Δt*f(u + a1400*k[1] + a1401*k[2]                           + a1404*k[5]              + a1406*k[7] +                                                                     a1412*k[13] + a1413*k[14],t + c[14]*Δt)
    k[16] = Δt*f(u + a1500*k[1]              + a1502*k[3]                                                                                                                                                     + a1514*k[15],t + c[15]*Δt)
    k[17] = Δt*f(u + a1600*k[1] + a1601*k[2] + a1602*k[3]              + a1604*k[5] + a1605*k[6] + a1606*k[7] + a1607*k[8] + a1608*k[9] + a1609*k[10] + a1610*k[11] + a1611*k[12] + a1612*k[13] + a1613*k[14] + a1614*k[15] + a1615*k[16],t + c[16]*Δt)
    if adaptive
      utmp = copy(u)
      for i=sumIdx
        utmp += b[i]*k[i]
      end
      EEst = norm(((k[2] - k[16]) / 360)./(abstol+u*reltol),internalnorm)
    else #no chance of rejecting, so in-place
      #=
      for i=1:17
        u[:]+=vec(b[i]*k[i])
      end
      =#
      #=
      for i in eachindex(u)
        u[i] = u[i] + b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i] + b[7]*k[7][i] + b[9]*k[9][i] + b[10]*k[10][i] + b[11]*k[11][i] + b[12]*k[12][i] + b[13]*k[13][i] + b[14]*k[14][i] + b[15]*k[15][i] + b[16]*k[16][i] + b[17]*k[17][i]
      end
      =#
      u = u + b[1]*k[1] + b[2]*k[2] + b[3]*k[3] + b[5]*k[5] + b[7]*k[7] + b[9]*k[9] + b[10]*k[10] + b[11]*k[11] + b[12]*k[12] + b[13]*k[13] + b[14]*k[14] + b[15]*k[15] + b[16]*k[16] + b[17]*k[17]
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_impliciteuler(f::Function,u,t,Δt,T,iter,maxiters,
                            timeseries,ts,timeseries_steps,save_timeseries,adaptive,sizeu,progressbar,autodiff)
  function rhsIE(u,resid,uOld,t,Δt)
    u = reshape(u,sizeu...)
    resid = reshape(resid,sizeu...)
    resid[:] = u - uOld - Δt*f(u,t+Δt)
    vec(u)
    vec(resid)
  end
  while t < T
    @ode_loopheader
    uOld = copy(u)
    u = vec(u)
    nlres = NLsolve.nlsolve((u,resid)->rhsIE(u,resid,uOld,t,Δt),u,autodiff=autodiff)
    u = reshape(nlres.zero,sizeu...)
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_trapezoid(f::Function,u,t,Δt,T,iter,maxiters,
                      timeseries,ts,timeseries_steps,save_timeseries,adaptive,sizeu,progressbar,autodiff)
  function rhsTrap(u,resid,uOld,t,Δt)
    u = reshape(u,sizeu...)
    resid = reshape(resid,sizeu...)
    resid[:] = u - uOld - Δt*(f(u,t+Δt)+f(uOld,t))/2
    u = vec(u)
    resid = vec(resid)
  end
  while t < T
    @ode_loopheader
    uOld = copy(u)
    u = vec(u)
    nlres = NLsolve.nlsolve((u,resid)->rhsTrap(u,resid,uOld,t,Δt),u,autodiff=autodiff)
    u = reshape(nlres.zero,sizeu...)
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_rosenbrock32(f::Function,u::AbstractArray,t,Δt,T,iter,
                          maxiters,timeseries,ts,timeseries_steps,save_timeseries,adaptive,
                          sizeu,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar,γ)
  order = 2
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  k₁ = similar(u)
  k₂ = similar(u)
  k₃ = similar(u)
  function vecf(u,t)
    return(vec(f(reshape(u,sizeu...),t)))
  end
  while t < T
    @ode_loopheader
    # Time derivative
    dT = ForwardDiff.derivative((t)->f(u,t),t)
    J = ForwardDiff.jacobian((u)->vecf(u,t),vec(u))
    W = one(J)-Δt*d*J
    f₀ = f(u,t)
    k₁[:] = reshape(W\vec(f₀ + Δt*d*dT),sizeu...)
    f₁ = f(u+Δt*k₁/2,t+Δt/2)
    k₂[:] = reshape(W\vec(f₁-k₁),sizeu...) + k₁
    if adaptive
      utmp = u + Δt*k₂
      f₂ = f(utmp,t+Δt)
      k₃[:] = reshape(W\vec(f₂ - c₃₂*(k₂-f₁)-2(k₁-f₀)+Δt*d*T),sizeu...)
      EEst = norm((Δt(k₁ - 2k₂ + k₃)/6)./(abstol+u*reltol),internalnorm)
    else
      u = u + Δt*k₂
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_rosenbrock32(f::Function,u::Number,t,Δt,T,iter,
                          maxiters,timeseries,ts,timeseries_steps,save_timeseries,adaptive,
                          sizeu,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar,γ)
  order = 2
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  function vecf(u,t)
    return(vec(f(reshape(u,sizeu...),t)))
  end
  while t < T
    @ode_loopheader
    # Time derivative
    dT = ForwardDiff.derivative((t)->f(u,t),t)
    J = ForwardDiff.jacobian((u)->vecf(u,t),vec(u))
    W = one(J)-Δt*d*J
    f₀ = f(u,t)
    k₁ = reshape(W\vec(f₀ + Δt*d*dT),sizeu...)
    f₁ = f(u+Δt*k₁/2,t+Δt/2)
    k₂ = reshape(W\vec(f₁-k₁),sizeu...) + k₁
    if adaptive
      utmp = u + Δt*k₂
      f₂ = f(utmp,t+Δt)
      k₃ = reshape(W\vec(f₂ - c₃₂*(k₂-f₁)-2(k₁-f₀)+Δt*d*T),sizeu...)
      EEst = norm((Δt(k₁ - 2k₂ + k₃)/6)./(abstol+u*reltol),internalnorm)
    else
      u = u + Δt*k₂
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end
