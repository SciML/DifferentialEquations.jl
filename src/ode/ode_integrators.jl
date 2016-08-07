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

@def ode_implicitsavevalues begin
  if save_timeseries && iter%timeseries_steps==0
    push!(timeseries,reshape(u,sizeu...))
    push!(ts,t)
  end
end

@def ode_numberimplicitsavevalues begin
  if save_timeseries && iter%timeseries_steps==0
    push!(timeseries,u[1])
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
      u = copy_if_possible!(u, utmp)
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

@def ode_implicitloopfooter begin
  if adaptive
    standard = abs(1/(γ*EEst))^(1/order)
    if isinf(standard)
        q = qmax
    else
       q = min(qmax,max(standard,eps()))
    end
    if q > 1
      t = t + Δt
      u = copy_if_possible!(u, utmp)
      @ode_implicitsavevalues
    end
    Δtpropose = min(Δtmax,q*Δt)
    Δt = max(min(Δtpropose,abs(T-t)),Δtmin) #abs to fix complex sqrt issue at end
  else #Not adaptive
    t = t + Δt
    @ode_implicitsavevalues
  end
  (progressbar && atomloaded && iter%progress_steps==0) ? Main.Atom.progress(t/T) : nothing #Use Atom's progressbar if loaded
end

@def ode_numberimplicitloopfooter begin
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
      @ode_numberimplicitsavevalues
    end
    Δtpropose = min(Δtmax,q*Δt)
    Δt = max(min(Δtpropose,abs(T-t)),Δtmin) #abs to fix complex sqrt issue at end
  else #Not adaptive
    t = t + Δt
    @ode_numberimplicitsavevalues
  end
  (progressbar && atomloaded && iter%progress_steps==0) ? Main.Atom.progress(t/T) : nothing #Use Atom's progressbar if loaded
end

function ode_euler(f::Function,u::Number,t,Δt,T,iter,maxiters,
                    timeseries,ts,timeseries_steps,save_timeseries,adaptive,progressbar)
  @inbounds while t < T
    @ode_loopheader
    u = u + Δt.*f(u,t)
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_euler(f::Function,u::AbstractArray,t,Δt,T,iter,maxiters,
                    timeseries,ts,timeseries_steps,save_timeseries,adaptive,progressbar)
  du = similar(u)
  @inbounds while t < T
    @ode_loopheader
    f(du,u,t)
    for i in eachindex(u)
      u[i] = u[i] + Δt*du[i]
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_midpoint(f::Function,u::Number,t,Δt,T,iter,
                      maxiters,timeseries,ts,timeseries_steps,save_timeseries,adaptive,progressbar)
  halfΔt = Δt/2
  @inbounds while t < T
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
  du = similar(u)
  @inbounds while t < T
    @ode_loopheader
    f(du,u,t)
    for i in eachindex(u)
      utilde[i] = u[i]+halfΔt*du[i]
    end
    f(du,utilde,t+halfΔt)
    for i in eachindex(u)
      u[i] = u[i] + Δt*du[i]
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_rk4(f::Function,u::Number,t,Δt,T,iter,maxiters,
                timeseries,ts,timeseries_steps,save_timeseries,adaptive,progressbar)
  halfΔt = Δt/2
  @inbounds while t < T
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
  tmp = similar(u)
  @inbounds while t < T
    @ode_loopheader
    f(k₁,u,t)
    ttmp = t+halfΔt
    for i in eachindex(u)
      tmp[i] = u[i]+halfΔt*k₁[i]
    end
    f(k₂,tmp,ttmp)
    for i in eachindex(u)
      tmp[i] = u[i]+halfΔt*k₂[i]
    end
    f(k₃,tmp,ttmp)
    for i in eachindex(u)
      tmp[i] = u[i]+Δt*k₃[i]
    end
    f(k₄,tmp,t+Δt)
    for i in eachindex(u)
      u[i] = u[i] + Δt*(k₁[i] + 2k₂[i] + 2k₃[i] + k₄[i])/6
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_explicitrk(f::Function,u::Number,t,Δt,T,iter,maxiters,timeseries,ts,
                        timeseries_steps,save_timeseries,A,c,α,αEEst,stages,order,γ,adaptive,
                        abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  ks = Array{typeof(u)}(stages)
  @inbounds while t < T
    @ode_loopheader
    for i = 1:stages
      utilde = zero(u)
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
  ks = Vector{typeof(u)}(0)
  for i = 1:stages
    push!(ks,similar(u))
  end
  utilde = similar(u)
  tmp = similar(u)
  utmp = zeros(u)
  uEEst = similar(u)
  @inbounds while t < T
    @ode_loopheader
    for i = 1:stages
      utilde[:] = zero(eltype(u))
      for j = 1:i-1
        for k in eachindex(u)
          utilde[k] += A[i,j]*ks[j][k]
        end
      end
      for k in eachindex(u)
        tmp[k] = u[k]+Δt*utilde[k]
      end
      f(ks[i],tmp,t+c[i]*Δt)
    end
    utilde[:] = α[1]*ks[1]
    for i = 2:stages
      for k in eachindex(u)
        utilde[k] += α[i]*ks[i][k]
      end
    end
    if adaptive
      for i in eachindex(u)
        utmp[i] = u[i] + Δt*utilde[i]
      end
      uEEst[:] = αEEst[1]*ks[1]
      for i = 2:stages
        for j in eachindex(u)
          uEEst[j] += αEEst[i]*ks[i][j]
        end
      end
      EEst = norm((utilde-uEEst)./(abstol+u*reltol),internalnorm)
    else
      for i in eachindex(u)
        u[i] = u[i] + Δt*utilde[i]
      end
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_feagin10_vectorized(f::Function,u::AbstractArray,t,Δt,T,iter,order,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b,c = constructFeagin10(eltype(u))
  k = Vector{typeof(u)}(0)
  for i = 1:17
    push!(k,similar(u))
  end
  update = similar(u)
  utmp = similar(u)
  @inbounds while t < T
    @ode_loopheader
    f(k[1],u,t); k[1]*=Δt
    f(k[2],u + a0100*k[1],t + c[1]*Δt); k[2]*=Δt
    f(k[3],u + a0200*k[1] + a0201*k[2],t + c[2]*Δt ); k[3]*=Δt
    f(k[4],u + a0300*k[1]              + a0302*k[3],t + c[3]*Δt); k[4]*=Δt
    f(k[5],u + a0400*k[1]              + a0402*k[3] + a0403*k[4],t + c[4]*Δt); k[5]*=Δt
    f(k[6],u + a0500*k[1]                           + a0503*k[4] + a0504*k[5],t + c[5]*Δt); k[6]*=Δt
    f(k[7],u + a0600*k[1]                           + a0603*k[4] + a0604*k[5] + a0605*k[6],t + c[6]*Δt); k[7]*=Δt
    f(k[8],u + a0700*k[1]                                        + a0704*k[5] + a0705*k[6] + a0706*k[7],t + c[7]*Δt); k[8]*=Δt
    f(k[9],u + a0800*k[1]                                                     + a0805*k[6] + a0806*k[7] + a0807*k[8],t + c[8]*Δt); k[9]*=Δt
    f(k[10],u + a0900*k[1]                                                     + a0905*k[6] + a0906*k[7] + a0907*k[8] + a0908*k[9],t + c[9]*Δt); k[10]*=Δt
    f(k[11],u + a1000*k[1]                                                     + a1005*k[6] + a1006*k[7] + a1007*k[8] + a1008*k[9] + a1009*k[10],t + c[10]*Δt); k[11]*=Δt
    f(k[12],u + a1100*k[1]                                                     + a1105*k[6] + a1106*k[7] + a1107*k[8] + a1108*k[9] + a1109*k[10] + a1110*k[11],t + c[11]*Δt); k[12]*=Δt
    f(k[13],u + a1200*k[1]                           + a1203*k[4] + a1204*k[5] + a1205*k[6] + a1206*k[7] + a1207*k[8] + a1208*k[9] + a1209*k[10] + a1210*k[11] + a1211*k[12],t + c[12]*Δt); k[13]*=Δt
    f(k[14],u + a1300*k[1]              + a1302*k[3] + a1303*k[4]              + a1305*k[6] + a1306*k[7] + a1307*k[8] + a1308*k[9] + a1309*k[10] + a1310*k[11] + a1311*k[12] + a1312*k[13],t + c[13]*Δt); k[14]*=Δt
    f(k[15],u + a1400*k[1] + a1401*k[2]                           + a1404*k[5]              + a1406*k[7] +                                                                     a1412*k[13] + a1413*k[14],t + c[14]*Δt); k[15]*=Δt
    f(k[16],u + a1500*k[1]              + a1502*k[3]                                                                                                                                                     + a1514*k[15],t + c[15]*Δt); k[16]*=Δt
    f(k[17],u + a1600*k[1] + a1601*k[2] + a1602*k[3]              + a1604*k[5] + a1605*k[6] + a1606*k[7] + a1607*k[8] + a1608*k[9] + a1609*k[10] + a1610*k[11] + a1611*k[12] + a1612*k[13] + a1613*k[14] + a1614*k[15] + a1615*k[16],t + c[16]*Δt); k[17]*=Δt
    for i in eachindex(u)
      update[i] = (b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i]) + (b[7]*k[7][i] + b[9]*k[9][i] + b[10]*k[10][i] + b[11]*k[11][i]) + (b[12]*k[12][i] + b[13]*k[13][i] + b[14]*k[14][i] + b[15]*k[15][i]) + (b[16]*k[16][i] + b[17]*k[17][i])
    end
    if adaptive
      for i in eachindex(u)
        utmp[i] = u[i] + update[i]
      end
      EEst = norm(((k[2] - k[16]) / 360)./(abstol+u*reltol),internalnorm)
    else #no chance of rejecting, so in-place
      for i in eachindex(u)
        u[i] = u[i] + update[i]
      end
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_feagin10(f::Function,u::AbstractArray,t,Δt,T,iter,order,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b,c = constructFeagin10(eltype(u))
  k = Vector{typeof(u)}(0)
  for i = 1:17
    push!(k,similar(u))
  end
  tmp = similar(u)
  utmp = similar(u)
  @inbounds while t < T
    @ode_loopheader
    f(k[1],u,t); k[1]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0100*k[1][i]
    end
    f(k[2],tmp,t + c[1]*Δt); k[2]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0200*k[1][i] + a0201*k[2][i]
    end
    f(k[3],tmp,t + c[2]*Δt ); k[3]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0300*k[1][i] + a0302*k[3][i]
    end
    f(k[4],tmp,t + c[3]*Δt); k[4]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0400*k[1][i] + a0402*k[3][i] + a0403*k[4][i]
    end
    f(k[5],tmp,t + c[4]*Δt); k[5]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0500*k[1][i] + a0503*k[4][i] + a0504*k[5][i]
    end
    f(k[6],tmp,t + c[5]*Δt); k[6]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0600*k[1][i] + a0603*k[4][i] + a0604*k[5][i] + a0605*k[6][i]
    end
    f(k[7],tmp,t + c[6]*Δt); k[7]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a0700*k[1][i] + a0704*k[5][i] + a0705*k[6][i]) + a0706*k[7][i]
    end
    f(k[8],tmp,t + c[7]*Δt); k[8]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a0800*k[1][i] + a0805*k[6][i] + a0806*k[7][i]) + a0807*k[8][i]
    end
    f(k[9],tmp,t + c[8]*Δt); k[9]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a0900*k[1][i] + a0905*k[6][i] + a0906*k[7][i]) + a0907*k[8][i] + a0908*k[9][i]
    end
    f(k[10],tmp,t + c[9]*Δt); k[10]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1000*k[1][i] + a1005*k[6][i] + a1006*k[7][i]) + a1007*k[8][i] + a1008*k[9][i] + a1009*k[10][i]
    end
    f(k[11],tmp,t + c[10]*Δt); k[11]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1100*k[1][i] + a1105*k[6][i] + a1106*k[7][i]) + (a1107*k[8][i] + a1108*k[9][i] + a1109*k[10][i] + a1110*k[11][i])
    end
    f(k[12],tmp,t + c[11]*Δt); k[12]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1200*k[1][i] + a1203*k[4][i] + a1204*k[5][i]) + (a1205*k[6][i] + a1206*k[7][i] + a1207*k[8][i] + a1208*k[9][i]) + (a1209*k[10][i] + a1210*k[11][i] + a1211*k[12][i])
    end
    f(k[13],tmp,t + c[12]*Δt); k[13]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1300*k[1][i] + a1302*k[3][i] + a1303*k[4][i]) + (a1305*k[6][i] + a1306*k[7][i] + a1307*k[8][i] + a1308*k[9][i]) + (a1309*k[10][i] + a1310*k[11][i] + a1311*k[12][i] + a1312*k[13][i])
    end
    f(k[14],tmp,t + c[13]*Δt); k[14]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1400*k[1][i] + a1401*k[2][i] + a1404*k[5][i]) + (a1406*k[7][i] + a1412*k[13][i] + a1413*k[14][i])
    end
    f(k[15],tmp,t + c[14]*Δt); k[15]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a1500*k[1][i] + a1502*k[3][i] + a1514*k[15][i]
    end
    f(k[16],tmp,t + c[15]*Δt); k[16]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1600*k[1][i] + a1601*k[2][i] + a1602*k[3][i]) + (a1604*k[5][i] + a1605*k[6][i] + a1606*k[7][i] + a1607*k[8][i]) + (a1608*k[9][i] + a1609*k[10][i] + a1610*k[11][i] + a1611*k[12][i]) + (a1612*k[13][i] + a1613*k[14][i] + a1614*k[15][i] + a1615*k[16][i])
    end
    f(k[17],tmp,t + c[16]*Δt); k[17]*=Δt
    for i in eachindex(u)
      tmp[i] = (b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i]) + (b[7]*k[7][i] + b[9]*k[9][i] + b[10]*k[10][i] + b[11]*k[11][i]) + (b[12]*k[12][i] + b[13]*k[13][i] + b[14]*k[14][i] + b[15]*k[15][i]) + (b[16]*k[16][i] + b[17]*k[17][i])
    end
    if adaptive
      for i in eachindex(u)
        utmp[i] = u[i] + tmp[i]
        tmp[i] = ((k[2][i] - k[16][i]) / 360)./(abstol+u[i]*reltol)
      end
      EEst = norm(tmp,internalnorm)
    else #no chance of rejecting, so in-place
      for i in eachindex(u)
        u[i] = u[i] + tmp[i]
      end
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_feagin10(f::Function,u::Number,t,Δt,T,iter,order,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b,c = constructFeagin10(eltype(u))
  k = Vector{typeof(u)}(17)
  sumIdx = [collect(1:3);5;7;collect(9:17)]
  @inbounds while t < T
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
    update = (b[1]*k[1] + b[2]*k[2] + b[3]*k[3] + b[5]*k[5]) + (b[7]*k[7] + b[9]*k[9] + b[10]*k[10] + b[11]*k[11]) + (b[12]*k[12] + b[13]*k[13] + b[14]*k[14] + b[15]*k[15]) + (b[16]*k[16] + b[17]*k[17])
    if adaptive
      utmp = u + update
      EEst = norm(((k[2] - k[16]) / 360)./(abstol+u*reltol),internalnorm)
    else
      u = u + update
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_feagin12_vectorized(f::Function,u::AbstractArray,t,Δt,T,iter,order,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,b,c = constructFeagin12(eltype(u))
  k = Vector{typeof(u)}(0)
  for i = 1:25
    push!(k,similar(u))
  end
  update = similar(u)
  utmp = similar(u)
  adaptiveConst = 49/640
  @inbounds while t < T
    @ode_loopheader
    f(k[1] ,u,t); k[1]*=Δt
    f(k[2] ,u + a0100*k[1],t + c[1]*Δt); k[2]*=Δt
    f(k[3] ,u + a0200*k[1] + a0201*k[2],t + c[2]*Δt ); k[3]*=Δt
    f(k[4] ,u + a0300*k[1]              + a0302*k[3],t + c[3]*Δt); k[4]*=Δt
    f(k[5] ,u + a0400*k[1]              + a0402*k[3] + a0403*k[4],t + c[4]*Δt); k[5]*=Δt
    f(k[6] ,u + a0500*k[1]                           + a0503*k[4] + a0504*k[5],t + c[5]*Δt); k[6]*=Δt
    f(k[7] ,u + a0600*k[1]                           + a0603*k[4] + a0604*k[5] + a0605*k[6],t + c[6]*Δt); k[7]*=Δt
    f(k[8] ,u + a0700*k[1]                                        + a0704*k[5] + a0705*k[6] + a0706*k[7],t + c[7]*Δt); k[8]*=Δt
    f(k[9] ,u + a0800*k[1]                                                     + a0805*k[6] + a0806*k[7] + a0807*k[8],t + c[8]*Δt); k[9]*=Δt
    f(k[10],u + a0900*k[1]                                                     + a0905*k[6] + a0906*k[7] + a0907*k[8] + a0908*k[9],t + c[9]*Δt); k[10]*=Δt
    f(k[11],u + a1000*k[1]                                                     + a1005*k[6] + a1006*k[7] + a1007*k[8] + a1008*k[9] + a1009*k[10],t + c[10]*Δt); k[11]*=Δt
    f(k[12],u + a1100*k[1]                                                     + a1105*k[6] + a1106*k[7] + a1107*k[8] + a1108*k[9] + a1109*k[10] + a1110*k[11],t + c[11]*Δt); k[12]*=Δt
    f(k[13],u + a1200*k[1]                                                                                            + a1208*k[9] + a1209*k[10] + a1210*k[11] + a1211*k[12],t + c[12]*Δt); k[13]*=Δt
    f(k[14],u + a1300*k[1]                                                                                            + a1308*k[9] + a1309*k[10] + a1310*k[11] + a1311*k[12] + a1312*k[13],t + c[13]*Δt); k[14]*=Δt
    f(k[15],u + a1400*k[1]                                                                                            + a1408*k[9] + a1409*k[10] + a1410*k[11] + a1411*k[12] + a1412*k[13] + a1413*k[14],t + c[14]*Δt); k[15]*=Δt
    f(k[16],u + a1500*k[1]                                                                                            + a1508*k[9] + a1509*k[10] + a1510*k[11] + a1511*k[12] + a1512*k[13] + a1513*k[14] + a1514*k[15],t + c[15]*Δt); k[16]*=Δt
    f(k[17],u + a1600*k[1]                                                                                            + a1608*k[9] + a1609*k[10] + a1610*k[11] + a1611*k[12] + a1612*k[13] + a1613*k[14] + a1614*k[15] + a1615*k[16],t + c[16]*Δt); k[17]*=Δt
    f(k[18],u + a1700*k[1]                                                     + a1705*k[6] + a1706*k[7] + a1707*k[8] + a1708*k[9] + a1709*k[10] + a1710*k[11] + a1711*k[12] + a1712*k[13] + a1713*k[14] + a1714*k[15] + a1715*k[16] + a1716*k[17],t + c[17]*Δt); k[18]*=Δt
    f(k[19],u + a1800*k[1]                                                     + a1805*k[6] + a1806*k[7] + a1807*k[8] + a1808*k[9] + a1809*k[10] + a1810*k[11] + a1811*k[12] + a1812*k[13] + a1813*k[14] + a1814*k[15] + a1815*k[16] + a1816*k[17] + a1817*k[18],t + c[18]*Δt); k[19]*=Δt
    f(k[20],u + a1900*k[1]                                        + a1904*k[5] + a1905*k[6] + a1906*k[7]              + a1908*k[9] + a1909*k[10] + a1910*k[11] + a1911*k[12] + a1912*k[13] + a1913*k[14] + a1914*k[15] + a1915*k[16] + a1916*k[17] + a1917*k[18] + a1918*k[19],t + c[19]*Δt); k[20]*=Δt
    f(k[21],u + a2000*k[1]                           + a2003*k[4] + a2004*k[5] + a2005*k[6]              + a2007*k[8]              + a2009*k[10] + a2010*k[11]                                                                                     + a2017*k[18] + a2018*k[19] + a2019*k[20],t + c[20]*Δt); k[21]*=Δt
    f(k[22],u + a2100*k[1]              + a2102*k[3] + a2103*k[4]                           + a2106*k[7] + a2107*k[8]              + a2109*k[10] + a2110*k[11]                                                                                     + a2117*k[18] + a2118*k[19] + a2119*k[20] + a2120*k[21],t + c[21]*Δt); k[22]*=Δt
    f(k[23],u + a2200*k[1] + a2201*k[2]                           + a2204*k[5]              + a2206*k[7]                                                                                                                                                                                     + a2220*k[21] + a2221*k[22],t + c[22]*Δt); k[23]*=Δt
    f(k[24],u + a2300*k[1]              + a2302*k[3]                                                                                                                                                                                                                                                                     + a2322*k[23],t + c[23]*Δt); k[24]*=Δt
    f(k[25],u + a2400*k[1] + a2401*k[2] + a2402*k[3]              + a2404*k[5]              + a2406*k[7] + a2407*k[8] + a2408*k[9] + a2409*k[10] + a2410*k[11] + a2411*k[12] + a2412*k[13] + a2413*k[14] + a2414*k[15] + a2415*k[16] + a2416*k[17] + a2417*k[18] + a2418*k[19] + a2419*k[20] + a2420*k[21] + a2421*k[22] + a2422*k[23] + a2423*k[24],t + c[24]*Δt); k[25]*=Δt

    for i in eachindex(u)
      update[i] = (b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i]) + (b[7]*k[7][i] + b[8]*k[8][i] + b[10]*k[10][i] + b[11]*k[11][i]) + (b[13]*k[13][i] + b[14]*k[14][i] + b[15]*k[15][i] + b[16]*k[16][i]) + (b[17]*k[17][i] + b[18]*k[18][i] + b[19]*k[19][i] + b[20]*k[20][i]) + (b[21]*k[21][i] + b[22]*k[22][i] + b[23]*k[23][i] + b[24]*k[24][i]) + b[25]*k[25][i]
    end
    if adaptive
      for i in eachindex(u)
        utmp[i] = u[i] + update[i]
      end
      EEst = norm(((k[1] - k[23]) * adaptiveConst)./(abstol+u*reltol),internalnorm)
    else #no chance of rejecting so in-place
      for i in eachindex(u)
        u[i] = u[i] + update[i]
      end
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_feagin12(f::Function,u::AbstractArray,t,Δt,T,iter,order,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,b,c = constructFeagin12(eltype(u))
  k = Vector{typeof(u)}(0)
  for i = 1:25
    push!(k,similar(u))
  end
  update = similar(u)
  utmp = similar(u)
  tmp = similar(u)
  adaptiveConst = 49/640
  @inbounds while t < T
    @ode_loopheader
    f(k[1] ,u,t); k[1]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0100*k[1][i]
    end
    f(k[2] ,tmp,t + c[1]*Δt); k[2]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0200*k[1][i] + a0201*k[2][i]
    end
    f(k[3] ,tmp,t + c[2]*Δt ); k[3]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0300*k[1][i] + a0302*k[3][i]
    end
    f(k[4] ,tmp,t + c[3]*Δt); k[4]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0400*k[1][i] + a0402*k[3][i] + a0403*k[4][i]
    end
    f(k[5] ,tmp,t + c[4]*Δt); k[5]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0500*k[1][i] + a0503*k[4][i] + a0504*k[5][i]
    end
    f(k[6] ,tmp,t + c[5]*Δt); k[6]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a0600*k[1][i] + a0603*k[4][i] + a0604*k[5][i]) + a0605*k[6][i]
    end
    f(k[7] ,tmp,t + c[6]*Δt); k[7]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a0700*k[1][i] + a0704*k[5][i] + a0705*k[6][i]) + a0706*k[7][i]
    end
    f(k[8] ,tmp,t + c[7]*Δt); k[8]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a0800*k[1][i] + a0805*k[6][i] + a0806*k[7][i]) + a0807*k[8][i]
    end
    f(k[9] ,tmp,t + c[8]*Δt); k[9]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a0900*k[1][i] + a0905*k[6][i] + a0906*k[7][i]) + (a0907*k[8][i] + a0908*k[9][i])
    end
    f(k[10],tmp,t + c[9]*Δt); k[10]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1000*k[1][i] + a1005*k[6][i] + a1006*k[7][i]) + (a1007*k[8][i] + a1008*k[9][i] + a1009*k[10][i])
    end
    f(k[11],tmp,t + c[10]*Δt); k[11]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1100*k[1][i] + a1105*k[6][i] + a1106*k[7][i]) + (a1107*k[8][i] + a1108*k[9][i] + a1109*k[10][i] + a1110*k[11][i])
    end
    f(k[12],tmp,t + c[11]*Δt); k[12]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1200*k[1][i] + a1208*k[9][i] + a1209*k[10][i]) + (a1210*k[11][i] + a1211*k[12][i])
    end
    f(k[13],tmp,t + c[12]*Δt); k[13]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1300*k[1][i] + a1308*k[9][i] + a1309*k[10][i]) + (a1310*k[11][i] + a1311*k[12][i] + a1312*k[13][i])
    end
    f(k[14],tmp,t + c[13]*Δt); k[14]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1400*k[1][i] + a1408*k[9][i] + a1409*k[10][i]) + (a1410*k[11][i] + a1411*k[12][i] + a1412*k[13][i] + a1413*k[14][i])
    end
    f(k[15],tmp,t + c[14]*Δt); k[15]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1500*k[1][i] + a1508*k[9][i] + a1509*k[10][i]) + (a1510*k[11][i] + a1511*k[12][i] + a1512*k[13][i] + a1513*k[14][i]) + a1514*k[15][i]
    end
    f(k[16],tmp,t + c[15]*Δt); k[16]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1600*k[1][i] + a1608*k[9][i] + a1609*k[10][i]) + (a1610*k[11][i] + a1611*k[12][i] + a1612*k[13][i] + a1613*k[14][i]) + (a1614*k[15][i] + a1615*k[16][i])
    end
    f(k[17],tmp,t + c[16]*Δt); k[17]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1700*k[1][i] + a1705*k[6][i] + a1706*k[7][i]) + (a1707*k[8][i] + a1708*k[9][i] + a1709*k[10][i] + a1710*k[11][i]) + (a1711*k[12][i] + a1712*k[13][i] + a1713*k[14][i] + a1714*k[15][i]) + (a1715*k[16][i] + a1716*k[17][i])
    end
    f(k[18],tmp,t + c[17]*Δt); k[18]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1800*k[1][i] + a1805*k[6][i] + a1806*k[7][i]) + (a1807*k[8][i] + a1808*k[9][i] + a1809*k[10][i] + a1810*k[11][i]) + (a1811*k[12][i] + a1812*k[13][i] + a1813*k[14][i] + a1814*k[15][i]) + (a1815*k[16][i] + a1816*k[17][i] + a1817*k[18][i])
    end
    f(k[19],tmp,t + c[18]*Δt); k[19]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1900*k[1][i] + a1904*k[5][i] + a1905*k[6][i]) + (a1906*k[7][i] + a1908*k[9][i] + a1909*k[10][i] + a1910*k[11][i]) + (a1911*k[12][i] + a1912*k[13][i] + a1913*k[14][i] + a1914*k[15][i]) + (a1915*k[16][i] + a1916*k[17][i] + a1917*k[18][i] + a1918*k[19][i])
    end
    f(k[20],tmp,t + c[19]*Δt); k[20]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2000*k[1][i] + a2003*k[4][i] + a2004*k[5][i]) + (a2005*k[6][i] + a2007*k[8][i] + a2009*k[10][i] + a2010*k[11][i]) + (a2017*k[18][i] + a2018*k[19][i] + a2019*k[20][i])
    end
    f(k[21],tmp,t + c[20]*Δt); k[21]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2100*k[1][i] + a2102*k[3][i] + a2103*k[4][i]) + (a2106*k[7][i] + a2107*k[8][i] + a2109*k[10][i] + a2110*k[11][i]) + (a2117*k[18][i] + a2118*k[19][i] + a2119*k[20][i] + a2120*k[21][i])
    end
    f(k[22],tmp,t + c[21]*Δt); k[22]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2200*k[1][i] + a2201*k[2][i] + a2204*k[5][i]) + (a2206*k[7][i] + a2220*k[21][i] + a2221*k[22][i])
    end
    f(k[23],tmp,t + c[22]*Δt); k[23]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a2300*k[1][i] + a2302*k[3][i] + a2322*k[23][i]
    end
    f(k[24],tmp,t + c[23]*Δt); k[24]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2400*k[1][i] + a2401*k[2][i] + a2402*k[3][i]) + (a2404*k[5][i] + a2406*k[7][i] + a2407*k[8][i] + a2408*k[9][i]) + (a2409*k[10][i] + a2410*k[11][i] + a2411*k[12][i] + a2412*k[13][i]) + (a2413*k[14][i] + a2414*k[15][i] + a2415*k[16][i] + a2416*k[17][i]) + (a2417*k[18][i] + a2418*k[19][i] + a2419*k[20][i] + a2420*k[21][i]) + (a2421*k[22][i] + a2422*k[23][i] + a2423*k[24][i])
    end
    f(k[25],tmp,t + c[24]*Δt); k[25]*=Δt

    for i in eachindex(u)
      update[i] = (b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i]) + (b[7]*k[7][i] + b[8]*k[8][i] + b[10]*k[10][i] + b[11]*k[11][i]) + (b[13]*k[13][i] + b[14]*k[14][i] + b[15]*k[15][i] + b[16]*k[16][i]) + (b[17]*k[17][i] + b[18]*k[18][i] + b[19]*k[19][i] + b[20]*k[20][i]) + (b[21]*k[21][i] + b[22]*k[22][i] + b[23]*k[23][i] + b[24]*k[24][i]) + b[25]*k[25][i]
    end
    if adaptive
      for i in eachindex(u)
        utmp[i] = u[i] + update[i]
        tmp[i] = ((k[1][i] - k[23][i]) * adaptiveConst)/(abstol+u[i]*reltol)
      end
      EEst = norm(tmp,internalnorm)
    else #no chance of rejecting so in-place
      for i in eachindex(u)
        u[i] = u[i] + update[i]
      end
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_feagin12(f::Function,u::Number,t,Δt,T,iter,order,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,b,c = constructFeagin12(eltype(u))
  k = Vector{typeof(u)}(25)
  adaptiveConst = 49/640
  @inbounds while t < T
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
    k[13] = Δt*f(u + a1200*k[1]                                                                                            + a1208*k[9] + a1209*k[10] + a1210*k[11] + a1211*k[12],t + c[12]*Δt)
    k[14] = Δt*f(u + a1300*k[1]                                                                                            + a1308*k[9] + a1309*k[10] + a1310*k[11] + a1311*k[12] + a1312*k[13],t + c[13]*Δt)
    k[15] = Δt*f(u + a1400*k[1]                                                                                            + a1408*k[9] + a1409*k[10] + a1410*k[11] + a1411*k[12] + a1412*k[13] + a1413*k[14],t + c[14]*Δt)
    k[16] = Δt*f(u + a1500*k[1]                                                                                            + a1508*k[9] + a1509*k[10] + a1510*k[11] + a1511*k[12] + a1512*k[13] + a1513*k[14] + a1514*k[15],t + c[15]*Δt)
    k[17] = Δt*f(u + a1600*k[1]                                                                                            + a1608*k[9] + a1609*k[10] + a1610*k[11] + a1611*k[12] + a1612*k[13] + a1613*k[14] + a1614*k[15] + a1615*k[16],t + c[16]*Δt)
    k[18] = Δt*f(u + a1700*k[1]                                                     + a1705*k[6] + a1706*k[7] + a1707*k[8] + a1708*k[9] + a1709*k[10] + a1710*k[11] + a1711*k[12] + a1712*k[13] + a1713*k[14] + a1714*k[15] + a1715*k[16] + a1716*k[17],t + c[17]*Δt)
    k[19] = Δt*f(u + a1800*k[1]                                                     + a1805*k[6] + a1806*k[7] + a1807*k[8] + a1808*k[9] + a1809*k[10] + a1810*k[11] + a1811*k[12] + a1812*k[13] + a1813*k[14] + a1814*k[15] + a1815*k[16] + a1816*k[17] + a1817*k[18],t + c[18]*Δt)
    k[20] = Δt*f(u + a1900*k[1]                                        + a1904*k[5] + a1905*k[6] + a1906*k[7]              + a1908*k[9] + a1909*k[10] + a1910*k[11] + a1911*k[12] + a1912*k[13] + a1913*k[14] + a1914*k[15] + a1915*k[16] + a1916*k[17] + a1917*k[18] + a1918*k[19],t + c[19]*Δt)
    k[21] = Δt*f(u + a2000*k[1]                           + a2003*k[4] + a2004*k[5] + a2005*k[6]              + a2007*k[8]              + a2009*k[10] + a2010*k[11]                                                                                     + a2017*k[18] + a2018*k[19] + a2019*k[20],t + c[20]*Δt)
    k[22] = Δt*f(u + a2100*k[1]              + a2102*k[3] + a2103*k[4]                           + a2106*k[7] + a2107*k[8]              + a2109*k[10] + a2110*k[11]                                                                                     + a2117*k[18] + a2118*k[19] + a2119*k[20] + a2120*k[21],t + c[21]*Δt)
    k[23] = Δt*f(u + a2200*k[1] + a2201*k[2]                           + a2204*k[5]              + a2206*k[7]                                                                                                                                                                                     + a2220*k[21] + a2221*k[22],t + c[22]*Δt)
    k[24] = Δt*f(u + a2300*k[1]              + a2302*k[3]                                                                                                                                                                                                                                                                     + a2322*k[23],t + c[23]*Δt)
    k[25] = Δt*f(u + a2400*k[1] + a2401*k[2] + a2402*k[3]              + a2404*k[5]              + a2406*k[7] + a2407*k[8] + a2408*k[9] + a2409*k[10] + a2410*k[11] + a2411*k[12] + a2412*k[13] + a2413*k[14] + a2414*k[15] + a2415*k[16] + a2416*k[17] + a2417*k[18] + a2418*k[19] + a2419*k[20] + a2420*k[21] + a2421*k[22] + a2422*k[23] + a2423*k[24],t + c[24]*Δt)

    update = (b[1]*k[1] + b[2]*k[2] + b[3]*k[3] + b[5]*k[5]) + (b[7]*k[7] + b[8]*k[8] + b[10]*k[10] + b[11]*k[11]) + (b[13]*k[13] + b[14]*k[14] + b[15]*k[15] + b[16]*k[16]) + (b[17]*k[17] + b[18]*k[18] + b[19]*k[19] + b[20]*k[20]) + (b[21]*k[21] + b[22]*k[22] + b[23]*k[23] + b[24]*k[24]) + (b[25]*k[25])
    if adaptive
      utmp = u + update
      EEst = norm(((k[1] - k[23]) * adaptiveConst)./(abstol+u*reltol),internalnorm)
    else #no chance of rejecting so in-place
      u = u + update
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_feagin14(f::Function,u::AbstractArray,t,Δt,T,iter,order,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,b,c = constructFeagin14(eltype(u))
  k = Vector{typeof(u)}(0)
  for i = 1:35
    push!(k,similar(u))
  end
  update = similar(u)
  utmp = similar(u)
  tmp = similar(u)
  adaptiveConst = 1/1000
  @inbounds while t < T
    @ode_loopheader
    f(k[1] ,u,t); k[1]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0100*k[1][i]
    end
    f(k[2] ,tmp,t + c[1]*Δt); k[2]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0200*k[1][i] + a0201*k[2][i]
    end
    f(k[3] ,tmp,t + c[2]*Δt ); k[3]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0300*k[1][i] + a0302*k[3][i]
    end
    f(k[4] ,tmp,t + c[3]*Δt); k[4]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0400*k[1][i] + a0402*k[3][i] + a0403*k[4][i]
    end
    f(k[5] ,tmp,t + c[4]*Δt); k[5]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a0500*k[1][i] + a0503*k[4][i] + a0504*k[5][i]
    end
    f(k[6] ,tmp,t + c[5]*Δt); k[6]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a0600*k[1][i] + a0603*k[4][i] + a0604*k[5][i]) + a0605*k[6][i]
    end
    f(k[7] ,tmp,t + c[6]*Δt); k[7]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a0700*k[1][i] + a0704*k[5][i] + a0705*k[6][i]) + a0706*k[7][i]
    end
    f(k[8] ,tmp,t + c[7]*Δt); k[8]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a0800*k[1][i] + a0805*k[6][i] + a0806*k[7][i]) + a0807*k[8][i]
    end
    f(k[9] ,tmp,t + c[8]*Δt); k[9]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a0900*k[1][i] + a0905*k[6][i] + a0906*k[7][i]) + a0907*k[8][i] + a0908*k[9][i]
    end
    f(k[10],tmp,t + c[9]*Δt); k[10]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1000*k[1][i] + a1005*k[6][i] + a1006*k[7][i]) + (a1007*k[8][i] + a1008*k[9][i] + a1009*k[10][i])
    end
    f(k[11],tmp,t + c[10]*Δt); k[11]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1100*k[1][i] + a1105*k[6][i] + a1106*k[7][i]) + (a1107*k[8][i] + a1108*k[9][i] + a1109*k[10][i] + a1110*k[11][i])
    end
    f(k[12],tmp,t + c[11]*Δt); k[12]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1200*k[1][i] + a1208*k[9][i] + a1209*k[10][i]) + (a1210*k[11][i] + a1211*k[12][i])
    end
    f(k[13],tmp,t + c[12]*Δt); k[13]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1300*k[1][i] + a1308*k[9][i] + a1309*k[10][i]) + (a1310*k[11][i] + a1311*k[12][i] + a1312*k[13][i])
    end
    f(k[14],tmp,t + c[13]*Δt); k[14]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1400*k[1][i] + a1408*k[9][i] + a1409*k[10][i]) + (a1410*k[11][i] + a1411*k[12][i] + a1412*k[13][i] + a1413*k[14][i])
    end
    f(k[15],tmp,t + c[14]*Δt); k[15]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1500*k[1][i] + a1508*k[9][i] + a1509*k[10][i]) + (a1510*k[11][i] + a1511*k[12][i] + a1512*k[13][i] + a1513*k[14][i]) + a1514*k[15][i]
    end
    f(k[16],tmp,t + c[15]*Δt); k[16]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1600*k[1][i] + a1608*k[9][i] + a1609*k[10][i]) + (a1610*k[11][i] + a1611*k[12][i] + a1612*k[13][i] + a1613*k[14][i]) + a1614*k[15][i] + a1615*k[16][i]
    end
    f(k[17],tmp,t + c[16]*Δt); k[17]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1700*k[1][i] + a1712*k[13][i] + a1713*k[14][i]) + (a1714*k[15][i] + a1715*k[16][i] + a1716*k[17][i])
    end
    f(k[18],tmp,t + c[17]*Δt); k[18]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1800*k[1][i] + a1812*k[13][i] + a1813*k[14][i]) + (a1814*k[15][i] + a1815*k[16][i] + a1816*k[17][i] + a1817*k[18][i])
    end
    f(k[19],tmp,t + c[18]*Δt); k[19]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a1900*k[1][i] + a1912*k[13][i] + a1913*k[14][i]) + (a1914*k[15][i] + a1915*k[16][i] + a1916*k[17][i] + a1917*k[18][i]) + a1918*k[19][i]
    end
    f(k[20],tmp,t + c[19]*Δt); k[20]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2000*k[1][i] + a2012*k[13][i] + a2013*k[14][i]) + (a2014*k[15][i] + a2015*k[16][i] + a2016*k[17][i] + a2017*k[18][i]) + (a2018*k[19][i] + a2019*k[20][i])
    end
    f(k[21],tmp,t + c[20]*Δt); k[21]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2100*k[1][i] + a2112*k[13][i] + a2113*k[14][i]) + (a2114*k[15][i] + a2115*k[16][i] + a2116*k[17][i] + a2117*k[18][i]) + (a2118*k[19][i] + a2119*k[20][i] + a2120*k[21][i])
    end
    f(k[22],tmp,t + c[21]*Δt); k[22]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2200*k[1][i] + a2212*k[13][i] + a2213*k[14][i]) + (a2214*k[15][i] + a2215*k[16][i] + a2216*k[17][i] + a2217*k[18][i]) + (a2218*k[19][i] + a2219*k[20][i] + a2220*k[21][i] + a2221*k[22][i])
    end
    f(k[23],tmp,t + c[22]*Δt); k[23]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2300*k[1][i] + a2308*k[9][i] + a2309*k[10][i]) + (a2310*k[11][i] + a2311*k[12][i] + a2312*k[13][i] + a2313*k[14][i]) + (a2314*k[15][i] + a2315*k[16][i] + a2316*k[17][i] + a2317*k[18][i]) + (a2318*k[19][i] + a2319*k[20][i] + a2320*k[21][i] + a2321*k[22][i]) + (a2322*k[23][i])
    end
    f(k[24],tmp,t + c[23]*Δt); k[24]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2400*k[1][i] + a2408*k[9][i] + a2409*k[10][i]) + (a2410*k[11][i] + a2411*k[12][i] + a2412*k[13][i] + a2413*k[14][i]) + (a2414*k[15][i] + a2415*k[16][i] + a2416*k[17][i] + a2417*k[18][i]) + (a2418*k[19][i] + a2419*k[20][i] + a2420*k[21][i] + a2421*k[22][i]) + (a2422*k[23][i] + a2423*k[24][i])
    end
    f(k[25],tmp,t + c[24]*Δt); k[25]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2500*k[1][i] + a2508*k[9][i] + a2509*k[10][i]) + (a2510*k[11][i] + a2511*k[12][i] + a2512*k[13][i] + a2513*k[14][i]) + (a2514*k[15][i] + a2515*k[16][i] + a2516*k[17][i] + a2517*k[18][i]) + (a2518*k[19][i] + a2519*k[20][i] + a2520*k[21][i] + a2521*k[22][i]) + (a2522*k[23][i] + a2523*k[24][i] + a2524*k[25][i])
    end
    f(k[26],tmp,t + c[25]*Δt); k[26]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2600*k[1][i] + a2605*k[6][i] + a2606*k[7][i]) + (a2607*k[8][i] + a2608*k[9][i] + a2609*k[10][i] + a2610*k[11][i]) + (a2612*k[13][i] + a2613*k[14][i] + a2614*k[15][i] + a2615*k[16][i]) + (a2616*k[17][i] + a2617*k[18][i] + a2618*k[19][i] + a2619*k[20][i]) + (a2620*k[21][i] + a2621*k[22][i] + a2622*k[23][i] + a2623*k[24][i]) + (a2624*k[25][i] + a2625*k[26][i])
    end
    f(k[27],tmp,t + c[26]*Δt); k[27]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2700*k[1][i] + a2705*k[6][i] + a2706*k[7][i]) + (a2707*k[8][i] + a2708*k[9][i] + a2709*k[10][i] + a2711*k[12][i]) + (a2712*k[13][i] + a2713*k[14][i] + a2714*k[15][i] + a2715*k[16][i]) + (a2716*k[17][i] + a2717*k[18][i] + a2718*k[19][i] + a2719*k[20][i]) + (a2720*k[21][i] + a2721*k[22][i] + a2722*k[23][i] + a2723*k[24][i]) + (a2724*k[25][i] + a2725*k[26][i] + a2726*k[27][i])
    end
    f(k[28],tmp,t + c[27]*Δt); k[28]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2800*k[1][i] + a2805*k[6][i] + a2806*k[7][i]) + (a2807*k[8][i] + a2808*k[9][i] + a2810*k[11][i] + a2811*k[12][i]) + (a2813*k[14][i] + a2814*k[15][i] + a2815*k[16][i] + a2823*k[24][i]) + (a2824*k[25][i] + a2825*k[26][i] + a2826*k[27][i] + a2827*k[28][i])
    end
    f(k[29],tmp,t + c[28]*Δt); k[29]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a2900*k[1][i] + a2904*k[5][i] + a2905*k[6][i]) + (a2906*k[7][i] + a2909*k[10][i] + a2910*k[11][i] + a2911*k[12][i]) + (a2913*k[14][i] + a2914*k[15][i] + a2915*k[16][i] + a2923*k[24][i]) + (a2924*k[25][i] + a2925*k[26][i] + a2926*k[27][i] + a2927*k[28][i]) + (a2928*k[29][i])
    end
    f(k[30],tmp,t + c[29]*Δt); k[30]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a3000*k[1][i] + a3003*k[4][i] + a3004*k[5][i]) + (a3005*k[6][i] + a3007*k[8][i] + a3009*k[10][i] + a3010*k[11][i]) + (a3013*k[14][i] + a3014*k[15][i] + a3015*k[16][i] + a3023*k[24][i]) + (a3024*k[25][i] + a3025*k[26][i] + a3027*k[28][i] + a3028*k[29][i]) + (a3029*k[30][i])
    end
    f(k[31],tmp,t + c[30]*Δt); k[31]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a3100*k[1][i] + a3102*k[3][i] + a3103*k[4][i]) + (a3106*k[7][i] + a3107*k[8][i] + a3109*k[10][i] + a3110*k[11][i]) + (a3113*k[14][i] + a3114*k[15][i] + a3115*k[16][i] + a3123*k[24][i]) + (a3124*k[25][i] + a3125*k[26][i] + a3127*k[28][i] + a3128*k[29][i]) + (a3129*k[30][i] + a3130*k[31][i])
    end
    f(k[32],tmp,t + c[31]*Δt); k[32]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a3200*k[1][i] + a3201*k[2][i] + a3204*k[5][i]) + (a3206*k[7][i] + a3230*k[31][i] + a3231*k[32][i])
    end
    f(k[33],tmp,t + c[32]*Δt); k[33]*=Δt
    for i in eachindex(u)
      tmp[i] = u[i] + a3300*k[1][i] + a3302*k[3][i] + a3332*k[33][i]
    end
    f(k[34],tmp,t + c[33]*Δt); k[34]*=Δt
    for i in eachindex(u)
      tmp[i] = (u[i] + a3400*k[1][i] + a3401*k[2][i] + a3402*k[3][i]) + (a3404*k[5][i] + a3406*k[7][i] + a3407*k[8][i] + a3409*k[10][i]) + (a3410*k[11][i] + a3411*k[12][i] + a3412*k[13][i] + a3413*k[14][i]) + (a3414*k[15][i] + a3415*k[16][i] + a3416*k[17][i] + a3417*k[18][i]) + (a3418*k[19][i] + a3419*k[20][i] + a3420*k[21][i] + a3421*k[22][i]) + (a3422*k[23][i] + a3423*k[24][i] + a3424*k[25][i] + a3425*k[26][i]) + (a3426*k[27][i] + a3427*k[28][i] + a3428*k[29][i] + a3429*k[30][i]) + (a3430*k[31][i] + a3431*k[32][i] + a3432*k[33][i] + a3433*k[34][i])
    end
    f(k[35],tmp,t + c[34]*Δt); k[35]*=Δt
    for i in eachindex(u)
      update[i] = (b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i]) + (b[7]*k[7][i] + b[8]*k[8][i] + b[10]*k[10][i] + b[11]*k[11][i]) + (b[12]*k[12][i] + b[14]*k[14][i] + b[15]*k[15][i] + b[16]*k[16][i]) + (b[18]*k[18][i] + b[19]*k[19][i] + b[20]*k[20][i] + b[21]*k[21][i]) + (b[22]*k[22][i] + b[23]*k[23][i] + b[24]*k[24][i] + b[25]*k[25][i]) + (b[26]*k[26][i] + b[27]*k[27][i] + b[28]*k[28][i] + b[29]*k[29][i]) + (b[30]*k[30][i] + b[31]*k[31][i] + b[32]*k[32][i] + b[33]*k[33][i]) + (b[34]*k[34][i] + b[35]*k[35][i])
    end
    if adaptive
      for i in eachindex(u)
        utmp[i] = u[i] + update[i]
        tmp[i] = ((k[1][i] - k[33][i]) * adaptiveConst)./(abstol+u[i]*reltol)
      end
      EEst = norm(tmp,internalnorm)
    else #no chance of rejecting, so in-place
      for i in eachindex(u)
        u[i] = u[i] + update[i]
      end
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_feagin14_vectorized(f::Function,u::AbstractArray,t,Δt,T,iter,order,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,b,c = constructFeagin14(eltype(u))
  k = Vector{typeof(u)}(0)
  for i = 1:35
    push!(k,similar(u))
  end
  update = similar(u)
  utmp = similar(u)
  adaptiveConst = 1/1000
  @inbounds while t < T
    @ode_loopheader
    f(k[1] ,u,t); k[1]*=Δt
    f(k[2] ,u + a0100*k[1],t + c[1]*Δt); k[2]*=Δt
    f(k[3] ,u + a0200*k[1] + a0201*k[2],t + c[2]*Δt ); k[3]*=Δt
    f(k[4] ,u + a0300*k[1]              + a0302*k[3],t + c[3]*Δt); k[4]*=Δt
    f(k[5] ,u + a0400*k[1]              + a0402*k[3] + a0403*k[4],t + c[4]*Δt); k[5]*=Δt
    f(k[6] ,u + a0500*k[1]                           + a0503*k[4] + a0504*k[5],t + c[5]*Δt); k[6]*=Δt
    f(k[7] ,u + a0600*k[1]                           + a0603*k[4] + a0604*k[5] + a0605*k[6],t + c[6]*Δt); k[7]*=Δt
    f(k[8] ,u + a0700*k[1]                                        + a0704*k[5] + a0705*k[6] + a0706*k[7],t + c[7]*Δt); k[8]*=Δt
    f(k[9] ,u + a0800*k[1]                                                     + a0805*k[6] + a0806*k[7] + a0807*k[8],t + c[8]*Δt); k[9]*=Δt
    f(k[10],u + a0900*k[1]                                                     + a0905*k[6] + a0906*k[7] + a0907*k[8] + a0908*k[9],t + c[9]*Δt); k[10]*=Δt
    f(k[11],u + a1000*k[1]                                                     + a1005*k[6] + a1006*k[7] + a1007*k[8] + a1008*k[9] + a1009*k[10],t + c[10]*Δt); k[11]*=Δt
    f(k[12],u + a1100*k[1]                                                     + a1105*k[6] + a1106*k[7] + a1107*k[8] + a1108*k[9] + a1109*k[10] + a1110*k[11],t + c[11]*Δt); k[12]*=Δt
    f(k[13],u + a1200*k[1]                                                                                            + a1208*k[9] + a1209*k[10] + a1210*k[11] + a1211*k[12],t + c[12]*Δt); k[13]*=Δt
    f(k[14],u + a1300*k[1]                                                                                            + a1308*k[9] + a1309*k[10] + a1310*k[11] + a1311*k[12] + a1312*k[13],t + c[13]*Δt); k[14]*=Δt
    f(k[15],u + a1400*k[1]                                                                                            + a1408*k[9] + a1409*k[10] + a1410*k[11] + a1411*k[12] + a1412*k[13] + a1413*k[14],t + c[14]*Δt); k[15]*=Δt
    f(k[16],u + a1500*k[1]                                                                                            + a1508*k[9] + a1509*k[10] + a1510*k[11] + a1511*k[12] + a1512*k[13] + a1513*k[14] + a1514*k[15],t + c[15]*Δt); k[16]*=Δt
    f(k[17],u + a1600*k[1]                                                                                            + a1608*k[9] + a1609*k[10] + a1610*k[11] + a1611*k[12] + a1612*k[13] + a1613*k[14] + a1614*k[15] + a1615*k[16],t + c[16]*Δt); k[17]*=Δt
    f(k[18],u + a1700*k[1]                                                                                                                                                   + a1712*k[13] + a1713*k[14] + a1714*k[15] + a1715*k[16] + a1716*k[17],t + c[17]*Δt); k[18]*=Δt
    f(k[19],u + a1800*k[1]                                                                                                                                                   + a1812*k[13] + a1813*k[14] + a1814*k[15] + a1815*k[16] + a1816*k[17] + a1817*k[18],t + c[18]*Δt); k[19]*=Δt
    f(k[20],u + a1900*k[1]                                                                                                                                                   + a1912*k[13] + a1913*k[14] + a1914*k[15] + a1915*k[16] + a1916*k[17] + a1917*k[18] + a1918*k[19],t + c[19]*Δt); k[20]*=Δt
    f(k[21],u + a2000*k[1]                                                                                                                                                   + a2012*k[13] + a2013*k[14] + a2014*k[15] + a2015*k[16] + a2016*k[17] + a2017*k[18] + a2018*k[19] + a2019*k[20],t + c[20]*Δt); k[21]*=Δt
    f(k[22],u + a2100*k[1]                                                                                                                                                   + a2112*k[13] + a2113*k[14] + a2114*k[15] + a2115*k[16] + a2116*k[17] + a2117*k[18] + a2118*k[19] + a2119*k[20] + a2120*k[21],t + c[21]*Δt); k[22]*=Δt
    f(k[23],u + a2200*k[1]                                                                                                                                                   + a2212*k[13] + a2213*k[14] + a2214*k[15] + a2215*k[16] + a2216*k[17] + a2217*k[18] + a2218*k[19] + a2219*k[20] + a2220*k[21] + a2221*k[22],t + c[22]*Δt); k[23]*=Δt
    f(k[24],u + a2300*k[1]                                                                                            + a2308*k[9] + a2309*k[10] + a2310*k[11] + a2311*k[12] + a2312*k[13] + a2313*k[14] + a2314*k[15] + a2315*k[16] + a2316*k[17] + a2317*k[18] + a2318*k[19] + a2319*k[20] + a2320*k[21] + a2321*k[22] + a2322*k[23],t + c[23]*Δt); k[24]*=Δt
    f(k[25],u + a2400*k[1]                                                                                            + a2408*k[9] + a2409*k[10] + a2410*k[11] + a2411*k[12] + a2412*k[13] + a2413*k[14] + a2414*k[15] + a2415*k[16] + a2416*k[17] + a2417*k[18] + a2418*k[19] + a2419*k[20] + a2420*k[21] + a2421*k[22] + a2422*k[23] + a2423*k[24],t + c[24]*Δt); k[25]*=Δt
    f(k[26],u + a2500*k[1]                                                                                            + a2508*k[9] + a2509*k[10] + a2510*k[11] + a2511*k[12] + a2512*k[13] + a2513*k[14] + a2514*k[15] + a2515*k[16] + a2516*k[17] + a2517*k[18] + a2518*k[19] + a2519*k[20] + a2520*k[21] + a2521*k[22] + a2522*k[23] + a2523*k[24] + a2524*k[25],t + c[25]*Δt); k[26]*=Δt
    f(k[27],u + a2600*k[1]                                                     + a2605*k[6] + a2606*k[7] + a2607*k[8] + a2608*k[9] + a2609*k[10] + a2610*k[11]               + a2612*k[13] + a2613*k[14] + a2614*k[15] + a2615*k[16] + a2616*k[17] + a2617*k[18] + a2618*k[19] + a2619*k[20] + a2620*k[21] + a2621*k[22] + a2622*k[23] + a2623*k[24] + a2624*k[25] + a2625*k[26],t + c[26]*Δt); k[27]*=Δt
    f(k[28],u + a2700*k[1]                                                     + a2705*k[6] + a2706*k[7] + a2707*k[8] + a2708*k[9] + a2709*k[10]               + a2711*k[12] + a2712*k[13] + a2713*k[14] + a2714*k[15] + a2715*k[16] + a2716*k[17] + a2717*k[18] + a2718*k[19] + a2719*k[20] + a2720*k[21] + a2721*k[22] + a2722*k[23] + a2723*k[24] + a2724*k[25] + a2725*k[26] + a2726*k[27],t + c[27]*Δt); k[28]*=Δt
    f(k[29],u + a2800*k[1]                                                     + a2805*k[6] + a2806*k[7] + a2807*k[8] + a2808*k[9]               + a2810*k[11] + a2811*k[12]               + a2813*k[14] + a2814*k[15] + a2815*k[16]                                                                                                   + a2823*k[24] + a2824*k[25] + a2825*k[26] + a2826*k[27] + a2827*k[28],t + c[28]*Δt); k[29]*=Δt
    f(k[30],u + a2900*k[1]                                        + a2904*k[5] + a2905*k[6] + a2906*k[7]                           + a2909*k[10] + a2910*k[11] + a2911*k[12]               + a2913*k[14] + a2914*k[15] + a2915*k[16]                                                                                                   + a2923*k[24] + a2924*k[25] + a2925*k[26] + a2926*k[27] + a2927*k[28] + a2928*k[29],t + c[29]*Δt); k[30]*=Δt
    f(k[31],u + a3000*k[1]                           + a3003*k[4] + a3004*k[5] + a3005*k[6]              + a3007*k[8]              + a3009*k[10] + a3010*k[11]                             + a3013*k[14] + a3014*k[15] + a3015*k[16]                                                                                                   + a3023*k[24] + a3024*k[25] + a3025*k[26]               + a3027*k[28] + a3028*k[29] + a3029*k[30],t + c[30]*Δt); k[31]*=Δt
    f(k[32],u + a3100*k[1]              + a3102*k[3] + a3103*k[4]                           + a3106*k[7] + a3107*k[8]              + a3109*k[10] + a3110*k[11]                             + a3113*k[14] + a3114*k[15] + a3115*k[16]                                                                                                   + a3123*k[24] + a3124*k[25] + a3125*k[26]               + a3127*k[28] + a3128*k[29] + a3129*k[30] + a3130*k[31],t + c[31]*Δt); k[32]*=Δt
    f(k[33],u + a3200*k[1] + a3201*k[2]                           + a3204*k[5]              + a3206*k[7]                                                                                                                                                                                                                                                                                                                                 + a3230*k[31] + a3231*k[32],t + c[32]*Δt); k[33]*=Δt
    f(k[34],u + a3300*k[1]              + a3302*k[3]                                                                                                                                                                                                                                                                                                                                                                                                                 + a3332*k[33],t + c[33]*Δt); k[34]*=Δt
    f(k[35],u + a3400*k[1] + a3401*k[2] + a3402*k[3]              + a3404*k[5]              + a3406*k[7] + a3407*k[8]              + a3409*k[10] + a3410*k[11] + a3411*k[12] + a3412*k[13] + a3413*k[14] + a3414*k[15] + a3415*k[16] + a3416*k[17] + a3417*k[18] + a3418*k[19] + a3419*k[20] + a3420*k[21] + a3421*k[22] + a3422*k[23] + a3423*k[24] + a3424*k[25] + a3425*k[26] + a3426*k[27] + a3427*k[28] + a3428*k[29] + a3429*k[30] + a3430*k[31] + a3431*k[32] + a3432*k[33] + a3433*k[34],t + c[34]*Δt); k[35]*=Δt
    for i in eachindex(u)
      update[i] = (b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i]) + (b[7]*k[7][i] + b[8]*k[8][i] + b[10]*k[10][i] + b[11]*k[11][i]) + (b[12]*k[12][i] + b[14]*k[14][i] + b[15]*k[15][i] + b[16]*k[16][i]) + (b[18]*k[18][i] + b[19]*k[19][i] + b[20]*k[20][i] + b[21]*k[21][i]) + (b[22]*k[22][i] + b[23]*k[23][i] + b[24]*k[24][i] + b[25]*k[25][i]) + (b[26]*k[26][i] + b[27]*k[27][i] + b[28]*k[28][i] + b[29]*k[29][i]) + (b[30]*k[30][i] + b[31]*k[31][i] + b[32]*k[32][i] + b[33]*k[33][i]) + (b[34]*k[34][i] + b[35]*k[35][i])
    end
    if adaptive
      for i in eachindex(u)
        utmp[i] = u[i] + update[i]
      end
      EEst = norm(((k[1] - k[33]) * adaptiveConst)./(abstol+u*reltol),internalnorm)
    else #no chance of rejecting, so in-place
      for i in eachindex(u)
        u[i] = u[i] + update[i]
      end
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_feagin14(f::Function,u::Number,t,Δt,T,iter,order,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,b,c = constructFeagin14(eltype(u))
  k = Vector{typeof(u)}(35)
  adaptiveConst = 1/1000
  @inbounds while t < T
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
    k[13] = Δt*f(u + a1200*k[1]                                                                                            + a1208*k[9] + a1209*k[10] + a1210*k[11] + a1211*k[12],t + c[12]*Δt)
    k[14] = Δt*f(u + a1300*k[1]                                                                                            + a1308*k[9] + a1309*k[10] + a1310*k[11] + a1311*k[12] + a1312*k[13],t + c[13]*Δt)
    k[15] = Δt*f(u + a1400*k[1]                                                                                            + a1408*k[9] + a1409*k[10] + a1410*k[11] + a1411*k[12] + a1412*k[13] + a1413*k[14],t + c[14]*Δt)
    k[16] = Δt*f(u + a1500*k[1]                                                                                            + a1508*k[9] + a1509*k[10] + a1510*k[11] + a1511*k[12] + a1512*k[13] + a1513*k[14] + a1514*k[15],t + c[15]*Δt)
    k[17] = Δt*f(u + a1600*k[1]                                                                                            + a1608*k[9] + a1609*k[10] + a1610*k[11] + a1611*k[12] + a1612*k[13] + a1613*k[14] + a1614*k[15] + a1615*k[16],t + c[16]*Δt)
    k[18] = Δt*f(u + a1700*k[1]                                                                                                                                                   + a1712*k[13] + a1713*k[14] + a1714*k[15] + a1715*k[16] + a1716*k[17],t + c[17]*Δt)
    k[19] = Δt*f(u + a1800*k[1]                                                                                                                                                   + a1812*k[13] + a1813*k[14] + a1814*k[15] + a1815*k[16] + a1816*k[17] + a1817*k[18],t + c[18]*Δt)
    k[20] = Δt*f(u + a1900*k[1]                                                                                                                                                   + a1912*k[13] + a1913*k[14] + a1914*k[15] + a1915*k[16] + a1916*k[17] + a1917*k[18] + a1918*k[19],t + c[19]*Δt)
    k[21] = Δt*f(u + a2000*k[1]                                                                                                                                                   + a2012*k[13] + a2013*k[14] + a2014*k[15] + a2015*k[16] + a2016*k[17] + a2017*k[18] + a2018*k[19] + a2019*k[20],t + c[20]*Δt)
    k[22] = Δt*f(u + a2100*k[1]                                                                                                                                                   + a2112*k[13] + a2113*k[14] + a2114*k[15] + a2115*k[16] + a2116*k[17] + a2117*k[18] + a2118*k[19] + a2119*k[20] + a2120*k[21],t + c[21]*Δt)
    k[23] = Δt*f(u + a2200*k[1]                                                                                                                                                   + a2212*k[13] + a2213*k[14] + a2214*k[15] + a2215*k[16] + a2216*k[17] + a2217*k[18] + a2218*k[19] + a2219*k[20] + a2220*k[21] + a2221*k[22],t + c[22]*Δt)
    k[24] = Δt*f(u + a2300*k[1]                                                                                            + a2308*k[9] + a2309*k[10] + a2310*k[11] + a2311*k[12] + a2312*k[13] + a2313*k[14] + a2314*k[15] + a2315*k[16] + a2316*k[17] + a2317*k[18] + a2318*k[19] + a2319*k[20] + a2320*k[21] + a2321*k[22] + a2322*k[23],t + c[23]*Δt)
    k[25] = Δt*f(u + a2400*k[1]                                                                                            + a2408*k[9] + a2409*k[10] + a2410*k[11] + a2411*k[12] + a2412*k[13] + a2413*k[14] + a2414*k[15] + a2415*k[16] + a2416*k[17] + a2417*k[18] + a2418*k[19] + a2419*k[20] + a2420*k[21] + a2421*k[22] + a2422*k[23] + a2423*k[24],t + c[24]*Δt)
    k[26] = Δt*f(u + a2500*k[1]                                                                                            + a2508*k[9] + a2509*k[10] + a2510*k[11] + a2511*k[12] + a2512*k[13] + a2513*k[14] + a2514*k[15] + a2515*k[16] + a2516*k[17] + a2517*k[18] + a2518*k[19] + a2519*k[20] + a2520*k[21] + a2521*k[22] + a2522*k[23] + a2523*k[24] + a2524*k[25],t + c[25]*Δt)
    k[27] = Δt*f(u + a2600*k[1]                                                     + a2605*k[6] + a2606*k[7] + a2607*k[8] + a2608*k[9] + a2609*k[10] + a2610*k[11]               + a2612*k[13] + a2613*k[14] + a2614*k[15] + a2615*k[16] + a2616*k[17] + a2617*k[18] + a2618*k[19] + a2619*k[20] + a2620*k[21] + a2621*k[22] + a2622*k[23] + a2623*k[24] + a2624*k[25] + a2625*k[26],t + c[26]*Δt)
    k[28] = Δt*f(u + a2700*k[1]                                                     + a2705*k[6] + a2706*k[7] + a2707*k[8] + a2708*k[9] + a2709*k[10]               + a2711*k[12] + a2712*k[13] + a2713*k[14] + a2714*k[15] + a2715*k[16] + a2716*k[17] + a2717*k[18] + a2718*k[19] + a2719*k[20] + a2720*k[21] + a2721*k[22] + a2722*k[23] + a2723*k[24] + a2724*k[25] + a2725*k[26] + a2726*k[27],t + c[27]*Δt)
    k[29] = Δt*f(u + a2800*k[1]                                                     + a2805*k[6] + a2806*k[7] + a2807*k[8] + a2808*k[9]               + a2810*k[11] + a2811*k[12]               + a2813*k[14] + a2814*k[15] + a2815*k[16]                                                                                                   + a2823*k[24] + a2824*k[25] + a2825*k[26] + a2826*k[27] + a2827*k[28],t + c[28]*Δt)
    k[30] = Δt*f(u + a2900*k[1]                                        + a2904*k[5] + a2905*k[6] + a2906*k[7]                           + a2909*k[10] + a2910*k[11] + a2911*k[12]               + a2913*k[14] + a2914*k[15] + a2915*k[16]                                                                                                   + a2923*k[24] + a2924*k[25] + a2925*k[26] + a2926*k[27] + a2927*k[28] + a2928*k[29],t + c[29]*Δt)
    k[31] = Δt*f(u + a3000*k[1]                           + a3003*k[4] + a3004*k[5] + a3005*k[6]              + a3007*k[8]              + a3009*k[10] + a3010*k[11]                             + a3013*k[14] + a3014*k[15] + a3015*k[16]                                                                                                   + a3023*k[24] + a3024*k[25] + a3025*k[26]               + a3027*k[28] + a3028*k[29] + a3029*k[30],t + c[30]*Δt)
    k[32] = Δt*f(u + a3100*k[1]              + a3102*k[3] + a3103*k[4]                           + a3106*k[7] + a3107*k[8]              + a3109*k[10] + a3110*k[11]                             + a3113*k[14] + a3114*k[15] + a3115*k[16]                                                                                                   + a3123*k[24] + a3124*k[25] + a3125*k[26]               + a3127*k[28] + a3128*k[29] + a3129*k[30] + a3130*k[31],t + c[31]*Δt)
    k[33] = Δt*f(u + a3200*k[1] + a3201*k[2]                           + a3204*k[5]              + a3206*k[7]                                                                                                                                                                                                                                                                                                                                 + a3230*k[31] + a3231*k[32],t + c[32]*Δt)
    k[34] = Δt*f(u + a3300*k[1]              + a3302*k[3]                                                                                                                                                                                                                                                                                                                                                                                                                 + a3332*k[33],t + c[33]*Δt)
    k[35] = Δt*f(u + a3400*k[1] + a3401*k[2] + a3402*k[3]              + a3404*k[5]              + a3406*k[7] + a3407*k[8]              + a3409*k[10] + a3410*k[11] + a3411*k[12] + a3412*k[13] + a3413*k[14] + a3414*k[15] + a3415*k[16] + a3416*k[17] + a3417*k[18] + a3418*k[19] + a3419*k[20] + a3420*k[21] + a3421*k[22] + a3422*k[23] + a3423*k[24] + a3424*k[25] + a3425*k[26] + a3426*k[27] + a3427*k[28] + a3428*k[29] + a3429*k[30] + a3430*k[31] + a3431*k[32] + a3432*k[33] + a3433*k[34],t + c[34]*Δt)
    update = (b[1]*k[1] + b[2]*k[2] + b[3]*k[3] + b[5]*k[5]) + (b[7]*k[7] + b[8]*k[8] + b[10]*k[10] + b[11]*k[11]) + (b[12]*k[12] + b[14]*k[14] + b[15]*k[15] + b[16]*k[16]) + (b[18]*k[18] + b[19]*k[19] + b[20]*k[20] + b[21]*k[21]) + (b[22]*k[22] + b[23]*k[23] + b[24]*k[24] + b[25]*k[25]) + (b[26]*k[26] + b[27]*k[27] + b[28]*k[28] + b[29]*k[29]) + (b[30]*k[30] + b[31]*k[31] + b[32]*k[32] + b[33]*k[33]) + (b[34]*k[34] + b[35]*k[35])
    if adaptive
      utmp = u + update
      EEst = norm(((k[1] - k[33]) * adaptiveConst)./(abstol+u*reltol),internalnorm)
    else #no chance of rejecting, so in-place
      u = u + update
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_impliciteuler(f::Function,u::Number,t,Δt,T,iter,maxiters,
                            timeseries,ts,timeseries_steps,save_timeseries,adaptive,sizeu,progressbar,autodiff)
  function rhs_ie(u,resid,u_old,t,Δt)
    resid[1] = u[1] - u_old[1] - Δt*f(u,t+Δt)[1]
  end
  u = [u]
  u_old = similar(u)
  @inbounds while t < T
    @ode_loopheader
    u_old[1] = u[1]
    nlres = NLsolve.nlsolve((u,resid)->rhs_ie(u,resid,u_old,t,Δt),u,autodiff=autodiff)
    u[1] = nlres.zero[1]
    @ode_numberimplicitloopfooter
  end
  u = u[1]
  return u,t,timeseries,ts
end

function ode_impliciteuler(f::Function,u::AbstractArray,t,Δt,T,iter,maxiters,
                            timeseries,ts,timeseries_steps,save_timeseries,adaptive,sizeu,progressbar,autodiff)
  if autodiff
    cache = DiffCache(u)
    rhs_ie = (u,resid,u_old,t,Δt,cache) -> begin
      du = get_du(cache, eltype(u))
      f(du,reshape(u,sizeu),t+Δt)
      for i in eachindex(u)
        resid[i] = u[i] - u_old[i] - Δt*du[i]
      end
    end
  else
    cache = similar(u)
    rhs_ie = (u,resid,u_old,t,Δt,du) -> begin
      f(du,reshape(u,sizeu),t+Δt)
      for i in eachindex(u)
        resid[i] = u[i] - u_old[i] - Δt*du[i]
      end
    end
  end

  u = vec(u); u_old = similar(u)
  @inbounds while t < T
    @ode_loopheader
    copy!(u_old,u)
    nlres = NLsolve.nlsolve((u,resid)->rhs_ie(u,resid,u_old,t,Δt,cache),u,autodiff=autodiff)
    u[:] = nlres.zero
    @ode_implicitloopfooter
  end
  u = reshape(u,sizeu...)
  return u,t,timeseries,ts
end

function ode_trapezoid(f::Function,u::AbstractArray,t,Δt,T,iter,maxiters,
                            timeseries,ts,timeseries_steps,save_timeseries,adaptive,sizeu,progressbar,autodiff)
  if autodiff
    cache1 = DiffCache(u)
    cache2 = DiffCache(u)
    Δto2 = Δt/2
    rhs_trap = (u,resid,u_old,t,Δt,cache1,cache2) -> begin
      du1 = get_du(cache1, eltype(u)); du2 = get_du(cache2, eltype(u_old))
      f(du2,reshape(u_old,sizeu),t)
      f(du1,reshape(u,sizeu),t+Δt)
      for i in eachindex(u)
        resid[i] = u[i] - u_old[i] - Δto2*(du1[i]+du2[i])
      end
    end
  else
    cache1 = similar(u)
    cache2 = similar(u)
    Δto2 = Δt/2
    rhs_trap = (u,resid,u_old,t,Δt,du1,du2) -> begin
      f(du2,reshape(u_old,sizeu),t)
      f(du1,reshape(u,sizeu),t+Δt)
      for i in eachindex(u)
        resid[i] = u[i] - u_old[i] - Δto2*(du1[i]+du2[i])
      end
    end
  end
  u = vec(u); u_old = similar(u)
  @inbounds while t < T
    @ode_loopheader
    copy!(u_old,u)
    nlres = NLsolve.nlsolve((u,resid)->rhs_trap(u,resid,u_old,t,Δt,cache1,cache2),u,autodiff=autodiff)
    u[:] = nlres.zero
    @ode_implicitloopfooter
  end
  u = reshape(u,sizeu...)
  return u,t,timeseries,ts
end

function ode_trapezoid(f::Function,u::Number,t,Δt,T,iter,maxiters,
                            timeseries,ts,timeseries_steps,save_timeseries,adaptive,sizeu,progressbar,autodiff)
  Δto2 = Δt/2
  function rhs_trap(u,resid,u_old,t,Δt)
    resid[1] = u[1] - u_old[1] - Δto2*(f(u,t+Δt)[1] + f(u_old,t)[1])
  end
  u = [u]
  u_old = similar(u)
  @inbounds while t < T
    @ode_loopheader
    u_old[1] = u[1]
    nlres = NLsolve.nlsolve((u,resid)->rhs_trap(u,resid,u_old,t,Δt),u,autodiff=autodiff)
    u[1] = nlres.zero[1]
    @ode_numberimplicitloopfooter
  end
  u = u[1]
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
  function vecf(du,u,t)
    return(vec(f(reshape(du,sizeu...),reshape(u,sizeu...),t)))
  end
  du1 = similar(u)
  du2 = similar(u)
  f₀ = similar(u)
  f₁ = similar(u)
  f₂ = similar(u)
  utmp = similar(u)
  @inbounds while t < T
    @ode_loopheader
    dT = ForwardDiff.derivative((t)->f(du2,u,t),t) # Time derivative
    J = ForwardDiff.jacobian((du1,u)->vecf(du1,u,t),du1,vec(u))
    W = one(J)-Δt*d*J
    f(f₀,u,t)
    k₁[:] = reshape(W\vec(f₀ + Δt*d*dT),sizeu...)
    for i in eachindex(u)
      utmp[i]=u[i]+Δt*k₁[i]/2
    end
    f(f₁,utmp,t+Δt/2)
    k₂[:] = reshape(W\vec(f₁-k₁),sizeu...) + k₁
    if adaptive
      for i in eachindex(u)
        utmp[i] = u[i] + Δt*k₂[i]
      end
      f(f₂,utmp,t+Δt)
      k₃[:] = reshape(W\vec(f₂ - c₃₂*(k₂-f₁)-2(k₁-f₀)+Δt*d*T),sizeu...)
      EEst = norm((Δt(k₁ - 2k₂ + k₃)/6)./(abstol+u*reltol),internalnorm)
    else
      for i in eachindex(u)
        u[i] = u[i] + Δt*k₂[i]
      end
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
  @inbounds while t < T
    @ode_loopheader
    # Time derivative
    dT = ForwardDiff.derivative((t)->f(u,t),t)
    J = ForwardDiff.derivative((u)->f(u,t),u)
    W = one(J)-Δt*d*J
    f₀ = f(u,t)
    k₁ = W\(f₀ + Δt*d*dT)
    f₁ = f(u+Δt*k₁/2,t+Δt/2)
    k₂ = W\(f₁-k₁) + k₁
    if adaptive
      utmp = u + Δt*k₂
      f₂ = f(utmp,t+Δt)
      k₃ = W\(f₂ - c₃₂*(k₂-f₁)-2(k₁-f₀)+Δt*d*T)
      EEst = norm((Δt*(k₁ - 2k₂ + k₃)/6)./(abstol+u*reltol),internalnorm)
    else
      u = u + Δt*k₂
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end
