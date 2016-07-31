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
    push!(timeseries,copy(u))
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
      utmp = u + b[1]*k[1] + b[2]*k[2] + b[3]*k[3] + b[5]*k[5] + b[7]*k[7] + b[9]*k[9] + b[10]*k[10] + b[11]*k[11] + b[12]*k[12] + b[13]*k[13] + b[14]*k[14] + b[15]*k[15] + b[16]*k[16] + b[17]*k[17]
      EEst = norm(((k[2] - k[16]) / 360)./(abstol+u*reltol),internalnorm)
    else #no chance of rejecting, so in-place
      #=
      #why doesn't this work?
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

function ode_feagin12(f::Function,u::AbstractArray,t,Δt,T,iter,order,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,b,c = constructFeagin12(eltype(u))
  k = Vector{typeof(u)}(25)
  sumIdx = [collect(1:3);5;7;collect(9:17)]
  adaptiveConst = 49/640
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

    if adaptive
      utmp = u + b[1]*k[1] + b[2]*k[2] + b[3]*k[3] + b[5]*k[5] + b[7]*k[7] + b[8]*k[8] + b[10]*k[10] + b[11]*k[11] + b[13]*k[13] + b[14]*k[14] + b[15]*k[15] + b[16]*k[16] + b[17]*k[17] + b[18]*k[18] + b[19]*k[19] + b[20]*k[20] + b[21]*k[21] + b[22]*k[22] + b[23]*k[23] + b[24]*k[24] + b[25]*k[25]
      EEst = norm(((k[1] - k[23]) * adaptiveConst)./(abstol+u*reltol),internalnorm)
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
      u = u + b[1]*k[1] + b[2]*k[2] + b[3]*k[3] + b[5]*k[5] + b[7]*k[7] + b[8]*k[8] + b[10]*k[10] + b[11]*k[11] + b[13]*k[13] + b[14]*k[14] + b[15]*k[15] + b[16]*k[16] + b[17]*k[17] + b[18]*k[18] + b[19]*k[19] + b[20]*k[20] + b[21]*k[21] + b[22]*k[22] + b[23]*k[23] + b[24]*k[24] + b[25]*k[25]
    end
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end

function ode_feagin14(f::Function,u::AbstractArray,t,Δt,T,iter,order,
                        maxiters,timeseries,ts,timeseries_steps,save_timeseries,
                        γ,adaptive,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,progressbar)
  a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,b,c = constructFeagin14(eltype(u))
  k = Vector{typeof(u)}(35)
  sumIdx = [collect(1:3);5;7;collect(9:17)]
  adaptiveConst = 1/1000
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
    if adaptive
      utmp = u + b[1]*k[1] + b[2]*k[2] + b[3]*k[3] + b[5]*k[5] + b[7]*k[7] + b[8]*k[8] + b[10]*k[10] + b[11]*k[11] + b[12]*k[12] + b[14]*k[14] + b[15]*k[15] + b[16]*k[16] + b[18]*k[18] + b[19]*k[19] + b[20]*k[20] + b[21]*k[21] + b[22]*k[22] + b[23]*k[23] + b[24]*k[24] + b[25]*k[25] + b[26]*k[26] + b[27]*k[27] + b[28]*k[28] + b[29]*k[29] + b[30]*k[30] + b[31]*k[31] + b[32]*k[32] + b[33]*k[33] + b[34]*k[34] + b[35]*k[35]
      EEst = norm(((k[1] - k[33]) * adaptiveConst)./(abstol+u*reltol),internalnorm)
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
      u = u + b[1]*k[1] + b[2]*k[2] + b[3]*k[3] + b[5]*k[5] + b[7]*k[7] + b[8]*k[8] + b[10]*k[10] + b[11]*k[11] + b[12]*k[12] + b[14]*k[14] + b[15]*k[15] + b[16]*k[16] + b[18]*k[18] + b[19]*k[19] + b[20]*k[20] + b[21]*k[21] + b[22]*k[22] + b[23]*k[23] + b[24]*k[24] + b[25]*k[25] + b[26]*k[26] + b[27]*k[27] + b[28]*k[28] + b[29]*k[29] + b[30]*k[30] + b[31]*k[31] + b[32]*k[32] + b[33]*k[33] + b[34]*k[34] + b[35]*k[35]
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
