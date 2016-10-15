function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:ExplicitRKVectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  local A::Matrix{uEltypeNoUnits}
  local c::Vector{uEltypeNoUnits}
  local α::Vector{uEltypeNoUnits}
  local αEEst::Vector{uEltypeNoUnits}
  local stages::Int
  uidx = eachindex(u)
  @unpack A,c,α,αEEst,stages,fsal = integrator.tableau
  kk = Vector{rateType}(0)
  for i = 1:stages
    push!(kk,rateType(sizeu))
  end
  A = A' # Transpose A to column major looping
  utilde = rateType(sizeu)
  tmp = similar(u)
  utmp = zeros(u)
  uEEst = rateType(sizeu)
  if fsal
    f(t,u,fsalfirst)
  end
  kk[1] = fsalfirst
  kk[end] = fsallast
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      #First
      if !fsal
        f(t,u,kk[1])
      end
      #Middle
      for i = 2:stages-1
        for l in uidx
          utilde[l] = zero(kk[1][1])
        end
        for j = 1:i-1
          utilde += A[j,i]*kk[j]
        end
        tmp = u+Δt*utilde
        f(t+c[i]*Δt,tmp,kk[i])
      end
      # Last
      for l in uidx
        utilde[l] = zero(kk[1][1])
      end
      for j = 1:stages-1
        utilde += A[j,end]*kk[j]
      end
      tmp = u+Δt*utilde
      f(t+c[end]*Δt,tmp,fsallast)
      #Accumulate
      utilde[:] = α[1]*kk[1]
      for i = 2:stages
        utilde += α[i]*kk[i]
      end
      if adaptive
        utmp = u + Δt*utilde
        uEEst[:] = αEEst[1]*kk[1]
        for i = 2:stages
          uEEst += αEEst[i]*kk[i]
        end
        EEst = sqrt( sum((Δt*(utilde-uEEst)./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor)
      else
        u = u + Δt*utilde
      end
      if calck
        k = kk[end]
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:BS3Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  a21,a32,a41,a42,a43,c1,c2,b1,b2,b3,b4  = constructBS3(uEltypeNoUnits)
  k1 = rateType(sizeu)
  k2 = rateType(sizeu)
  k3 = rateType(sizeu)
  k4 = rateType(sizeu)
  local utilde::uType
  f(t,u,fsalfirst) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = fsalfirst
      f(t+c1*Δt,u+Δt*a21*k1,k2)
      f(t+c2*Δt,u+Δt*a32*k2,k3)
      utmp = u+Δt*(a41*k1+a42*k2+a43*k3)
      f(t+Δt,utmp,k4); fsallast = k4
      if adaptive
        utilde = u + Δt*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
        EEst = sqrt( sum(((utilde-utmp)./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if calck
        k = fsallast
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:BS5Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8   = constructBS5(uEltypeNoUnits)
  k1::rateType = rateType(sizeu)
  k2::rateType = rateType(sizeu)
  k3::rateType = rateType(sizeu)
  k4::rateType = rateType(sizeu)
  k5::rateType = rateType(sizeu)
  k6::rateType = rateType(sizeu)
  k7::rateType = rateType(sizeu)
  k8::rateType = rateType(sizeu)
  const kshortsize = 8
  local utilde::uType
  local uhat::uType
  local EEst2::uEltypeNoUnits
  if calck
    k = ksEltype()
    for i in 1:kshortsize
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    if calcprevs
      kprev = deepcopy(k)
      for i in 1:3 # Make it full-sized
        push!(kprev,rateType(sizeu))
      end
    end
  end
  f(t,u,fsalfirst) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = fsalfirst
      f(t+c1*Δt,u+Δt*(a21*k1),k2)
      f(t+c2*Δt,u+Δt*(a31*k1+a32*k2),k3)
      f(t+c3*Δt,u+Δt*(a41*k1+a42*k2+a43*k3),k4)
      f(t+c4*Δt,u+Δt*(a51*k1+a52*k2+a53*k3+a54*k4),k5)
      f(t+c5*Δt,u+Δt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5),k6)
      f(t+Δt,u+Δt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6),k7)
      utmp = u+Δt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
      f(t+Δt,utmp,fsallast); k8 = fsallast
      if adaptive
        uhat   = Δt*(bhat1*k1 + bhat3*k3 + bhat4*k4 + bhat5*k5 + bhat6*k6)
        utilde = u + Δt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7 + btilde8*k8)
        EEst1 = sqrt( sum(((uhat)./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor)
        EEst2 = sqrt( sum(((utilde-utmp)./(abstol+max(abs.(u),abs(utmp))*reltol)).^2) * normfactor)
        EEst = max(EEst1,EEst2)
      else
        u = utmp
      end
      if calck
        k[1]=k1; k[2]=k2; k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Tsit5Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,b1,b2,b3,b4,b5,b6,b7 = constructTsit5(uEltypeNoUnits)
  k1::rateType = rateType(sizeu)
  k2::rateType = rateType(sizeu)
  k3::rateType = rateType(sizeu)
  k4::rateType = rateType(sizeu)
  k5::rateType = rateType(sizeu)
  k6::rateType = rateType(sizeu)
  k7::rateType = rateType(sizeu)
  utilde::uType = similar(u)
  const kshortsize = 7
  if calck
    k = ksEltype()
    for i in 1:kshortsize
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    if calcprevs
      kprev = deepcopy(k)
    end
  end
  f(t,u,fsalfirst) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = fsalfirst
      f(t+c1*Δt,u+Δt*(a21*k1),k2)
      f(t+c2*Δt,u+Δt*(a31*k1+a32*k2),k3)
      f(t+c3*Δt,u+Δt*(a41*k1+a42*k2+a43*k3),k4)
      f(t+c4*Δt,u+Δt*(a51*k1+a52*k2+a53*k3+a54*k4),k5)
      f(t+Δt,u+Δt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5),k6)
      utmp = u+Δt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
      f(t+Δt,utmp,k7); fsallast = k7
      if adaptive
        utilde = u + Δt*(b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7)
        EEst = sqrt( sum(((utilde-utmp)./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if calck
        k[1] = k1
        k[2] = k2
        k[3] = k3
        k[4] = k4
        k[5] = k5
        k[6] = k6
        k[7] = k7
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:DP5Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = constructDP5(uEltypeNoUnits)
  k1 = rateType(sizeu)
  k2 = rateType(sizeu)
  k3 = rateType(sizeu)
  k4 = rateType(sizeu)
  k5 = rateType(sizeu)
  k6 = rateType(sizeu)
  k7 = rateType(sizeu)
  update = rateType(sizeu)
  utilde = similar(u)
  const kshortsize = 4
  if calck
    d1,d3,d4,d5,d6,d7 = DP5_dense_ds(uEltypeNoUnits)
    k = ksEltype()
    for i in 1:kshortsize
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    if calcprevs
      kprev = deepcopy(k)
    end
  end
  f(t,u,fsalfirst); #Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1=fsalfirst
      f(t+c1*Δt,u+Δt*(a21*k1),k2)
      f(t+c2*Δt,u+Δt*(a31*k1+a32*k2),k3)
      f(t+c3*Δt,u+Δt*(a41*k1+a42*k2+a43*k3),k4)
      f(t+c4*Δt,u+Δt*(a51*k1+a52*k2+a53*k3+a54*k4),k5)
      f(t+Δt,u+Δt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5),k6)
      update = a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
      utmp = u+Δt*update
      f(t+Δt,utmp,fsallast); k7=fsallast
      if adaptive
        utilde = u + Δt*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7)
        EEst = sqrt( sum(((utilde-utmp)./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if calck
        k[1] = update
        bspl = k1 - update
        k[2] = bspl
        k[3] = update - k7 - bspl
        k[4] = (d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Vern6Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,b1,b4,b5,b6,b7,b8,b9= constructVern6(uEltypeNoUnits)
  k1 = rateType(sizeu); k2 = rateType(sizeu) ; k3 = rateType(sizeu); k4 = rateType(sizeu)
  k5 = rateType(sizeu); k6 = rateType(sizeu) ; k7 = rateType(sizeu); k8 = rateType(sizeu)
  k9 = rateType(sizeu)
  const kshortsize = 9
  if calck
    k = ksEltype()
    for i in 1:kshortsize
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    if calcprevs
      kprev = deepcopy(k)
      for i in 1:3
        push!(kprev,rateType(sizeu))
      end
    end
  end
  utilde = similar(u)
  f(t,u,fsalfirst) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = fsalfirst
      f(t+c1*Δt,u+Δt*(a21*k1),k2)
      f(t+c2*Δt,u+Δt*(a31*k1+a32*k2),k3)
      f(t+c3*Δt,u+Δt*(a41*k1       +a43*k3),k4)
      f(t+c4*Δt,u+Δt*(a51*k1       +a53*k3+a54*k4),k5)
      f(t+c5*Δt,u+Δt*(a61*k1       +a63*k3+a64*k4+a65*k5),k6)
      f(t+c6*Δt,u+Δt*(a71*k1       +a73*k3+a74*k4+a75*k5+a76*k6),k7)
      f(t+Δt,u+Δt*(a81*k1       +a83*k3+a84*k4+a85*k5+a86*k6+a87*k7),k8)
      utmp=u+Δt*(a91*k1              +a94*k4+a95*k5+a96*k6+a97*k7+a98*k8)
      f(t+Δt,utmp,fsallast); k9 = fsallast
      if adaptive
        utilde = u + Δt*(b1*k1 + b4*k4 + b5*k5 + b6*k6 + b7*k7 + b8*k8 + b9*k9)
        EEst = sqrt( sum(((utilde-utmp)./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if calck
        k[1]=k1; k[2]=k2; k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8; k[9]=k9
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Vern7Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c2,c3,c4,c5,c6,c7,c8,a021,a031,a032,a041,a043,a051,a053,a054,a061,a063,a064,a065,a071,a073,a074,a075,a076,a081,a083,a084,a085,a086,a087,a091,a093,a094,a095,a096,a097,a098,a101,a103,a104,a105,a106,a107,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat10= constructVern7(uEltypeNoUnits)
  k1 = rateType(sizeu); k2 = rateType(sizeu); k3 = rateType(sizeu); k4 = rateType(sizeu);
  k5 = rateType(sizeu); k6 = rateType(sizeu); k7 = rateType(sizeu); k8 = rateType(sizeu);
  k9 = rateType(sizeu); k10 = rateType(sizeu); utilde = similar(u); update = similar(u)
  const kshortsize = 10
  if calck
    k = ksEltype()
    for i in 1:kshortsize
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    if calcprevs
      kprev = deepcopy(k)
      for i in 1:6 # Make it full-sized
        push!(kprev,rateType(sizeu))
      end
    end
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k1)
      f(t+c2*Δt,u+Δt*(a021*k1),k2)
      f(t+c3*Δt,u+Δt*(a031*k1+a032*k2),k3)
      f(t+c4*Δt,u+Δt*(a041*k1       +a043*k3),k4)
      f(t+c5*Δt,u+Δt*(a051*k1       +a053*k3+a054*k4),k5)
      f(t+c6*Δt,u+Δt*(a061*k1       +a063*k3+a064*k4+a065*k5),k6)
      f(t+c7*Δt,u+Δt*(a071*k1       +a073*k3+a074*k4+a075*k5+a076*k6),k7)
      f(t+c8*Δt,u+Δt*(a081*k1       +a083*k3+a084*k4+a085*k5+a086*k6+a087*k7),k8)
      f(t+Δt,u+Δt*(a091*k1          +a093*k3+a094*k4+a095*k5+a096*k6+a097*k7+a098*k8),k9)
      f(t+Δt,u+Δt*(a101*k1          +a103*k3+a104*k4+a105*k5+a106*k6+a107*k7),k10)
      update = Δt*(k1*b1 + k4*b4 + k5*b5 + k6*b6 + k7*b7 + k8*b8 + k9*b9)
      utmp = u + update
      if adaptive
        EEst = sqrt( sum(((update - Δt*(bhat1*k1 + bhat4*k4 + bhat5*k5 + bhat6*k6 + bhat7*k7 + bhat10*k10))./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if calck
        k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Vern8Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13= constructVern8(uEltypeNoUnits)
  k1 = rateType(sizeu); k2 = rateType(sizeu); k3 = rateType(sizeu); k4 = rateType(sizeu);
  k5 = rateType(sizeu); k6 = rateType(sizeu); k7 = rateType(sizeu); k8 = rateType(sizeu);
  k9 = rateType(sizeu); k10 = rateType(sizeu); k11 = rateType(sizeu); k12 = rateType(sizeu); k13 = rateType(sizeu)
  utilde = similar(u); update = similar(u);
  const kshortsize = 13
  if calck
    k = ksEltype()
    for i in 1:kshortsize
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    if calcprevs
      kprev = deepcopy(k)
      for i in 1:8 # Make it full-sized
        push!(kprev,rateType(sizeu))
      end
    end
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k1)
      f(t+c2*Δt ,u+Δt*(a0201*k1),k2)
      f(t+c3*Δt ,u+Δt*(a0301*k1+a0302*k2),k3)
      f(t+c4*Δt ,u+Δt*(a0401*k1       +a0403*k3),k4)
      f(t+c5*Δt ,u+Δt*(a0501*k1       +a0503*k3+a0504*k4),k5)
      f(t+c6*Δt ,u+Δt*(a0601*k1                +a0604*k4+a0605*k5),k6)
      f(t+c7*Δt ,u+Δt*(a0701*k1                +a0704*k4+a0705*k5+a0706*k6),k7)
      f(t+c8*Δt ,u+Δt*(a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7),k8)
      f(t+c9*Δt ,u+Δt*(a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8),k9)
      f(t+c10*Δt,u+Δt*(a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9),k10)
      f(t+c11*Δt,u+Δt*(a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10),k11)
      f(t+    Δt,u+Δt*(a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11),k12)
      f(t+    Δt,u+Δt*(a1301*k1                +a1304*k4+a1305*k5+a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10),k13)
      update = Δt*(k1*b1 + k6*b6 + k7*b7 + k8*b8 + k9*b9 + k10*b10 + k11*b11 + k12*b12)
      utmp = u + update
      if adaptive
        EEst = sqrt( sum(((update - Δt*(bhat1*k1 + bhat6*k6 + bhat7*k7 + bhat8*k8 + bhat9*k9 + bhat10*k10 + bhat13*k13))./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if calck
        k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10;k[11]=k11;k[12]=k12;k[13]=k13
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:TanYam7Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a93,a94,a95,a96,a97,a98,a101,a103,a104,a105,a106,a107,a108,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat8,bhat10 = constructTanYam7(uEltypeNoUnits)
  k1 = rateType(sizeu); k2 = rateType(sizeu) ; k3 = rateType(sizeu); k4 = rateType(sizeu)
  k5 = rateType(sizeu); k6 = rateType(sizeu) ; k7 = rateType(sizeu); k8 = rateType(sizeu)
  k9 = rateType(sizeu); k10= rateType(sizeu) ; k  = rateType(sizeu)
  utilde = similar(u);
  if calck
    pop!(ks) # Get rid of the one it starts with
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k); k1=k
      f(t+c1*Δt,u+Δt*(a21*k1),k2)
      f(t+c2*Δt,u+Δt*(a31*k1+a32*k2),k3)
      f(t+c3*Δt,u+Δt*(a41*k1       +a43*k3),k4)
      f(t+c4*Δt,u+Δt*(a51*k1       +a53*k3+a54*k4),k5)
      f(t+c5*Δt,u+Δt*(a61*k1       +a63*k3+a64*k4+a65*k5),k6)
      f(t+c6*Δt,u+Δt*(a71*k1       +a73*k3+a74*k4+a75*k5+a76*k6),k7)
      f(t+c7*Δt,u+Δt*(a81*k1       +a83*k3+a84*k4+a85*k5+a86*k6+a87*k7),k8)
      f(t+Δt,u+Δt*(a91*k1       +a93*k3+a94*k4+a95*k5+a96*k6+a97*k7+a98*k8),k9)
      f(t+Δt,u+Δt*(a101*k1      +a103*k3+a104*k4+a105*k5+a106*k6+a107*k7+a108*k8),k10)
      utmp = u + Δt*(k1*b1+k4*b4+k5*b5+k6*b6+k7*b7+k8*b8+k9*b9)
      if adaptive
        utilde = u + Δt*(k1*bhat1+k4*bhat4+k5*bhat5+k6*bhat6+k7*bhat7+k8*bhat8+k10*bhat10)
        EEst = sqrt( sum(((utilde-utmp)./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  if calck
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end


function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:DP8Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,b1,b6,b7,b8,b9,b10,b11,b12,bhh1,bhh2,bhh3,er1,er6,er7,er8,er9,er10,er11,er12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = constructDP8(uEltypeNoUnits)
  k1 = rateType(sizeu); k2  = rateType(sizeu); k3  = rateType(sizeu);  k4 = rateType(sizeu)
  k5 = rateType(sizeu); k6  = rateType(sizeu); k7  = rateType(sizeu);  k8 = rateType(sizeu)
  k9 = rateType(sizeu); k10 = rateType(sizeu); k11 = rateType(sizeu); k12 = rateType(sizeu)
  kupdate = rateType(sizeu); local k14::rateType; local k15::rateType; local k16::rateType;
  local udiff::rateType; local bspl::rateType
  utilde = similar(u); err5 = similar(u); err3 = similar(u)
  const kshortsize = 7
  if calck
    c14,c15,c16,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = DP8Interp(uEltypeNoUnits)
    d401,d406,d407,d408,d409,d410,d411,d412,d413,d414,d415,d416,d501,d506,d507,d508,d509,d510,d511,d512,d513,d514,d515,d516,d601,d606,d607,d608,d609,d610,d611,d612,d613,d614,d615,d616,d701,d706,d707,d708,d709,d710,d711,d712,d713,d714,d715,d716 = DP8Interp_polyweights(uEltypeNoUnits)
    if calck
      k = ksEltype()
      for i in 1:kshortsize
        push!(k,rateType(sizeu))
      end
      push!(ks,deepcopy(k)) #Initialize ks
      if calcprevs
        kprev = deepcopy(k)
      end
    end
    k13 = rateType(sizeu)
    k14 = rateType(sizeu)
    k15 = rateType(sizeu)
    k16 = rateType(sizeu)
    udiff = rateType(sizeu)
    bspl = rateType(sizeu)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k1)
      f(t+c2*Δt,u+Δt*(a0201*k1),k2)
      f(t+c3*Δt,u+Δt*(a0301*k1+a0302*k2),k3)
      f(t+c4*Δt,u+Δt*(a0401*k1       +a0403*k3),k4)
      f(t+c5*Δt,u+Δt*(a0501*k1       +a0503*k3+a0504*k4),k5)
      f(t+c6*Δt,u+Δt*(a0601*k1                +a0604*k4+a0605*k5),k6)
      f(t+c7*Δt,u+Δt*(a0701*k1                +a0704*k4+a0705*k5+a0706*k6),k7)
      f(t+c8*Δt,u+Δt*(a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7),k8)
      f(t+c9*Δt,u+Δt*(a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8),k9)
      f(t+c10*Δt,u+Δt*(a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9),k10)
      f(t+c11*Δt,u+Δt*(a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10),k11)
      f(t+Δt,u+Δt*(a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11),k12)
      kupdate = b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12
      update = Δt*kupdate
      utmp = u + update
      if adaptive
        err5 = sqrt(sum((Δt*(k1*er1 + k6*er6 + k7*er7 + k8*er8 + k9*er9 + k10*er10 + k11*er11 + k12*er12)./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor) # Order 5
        err3 = sqrt(sum(((update - Δt*(bhh1*k1 + bhh2*k9 + bhh3*k12))./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor) # Order 3
        err52 = err5*err5
        EEst = err52/sqrt(err52 + 0.01*err3*err3)
      else
        u = utmp
      end
      if calck
        f(t+Δt,utmp,k13)
        f(t+c14*Δt,u+Δt*(a1401*k1         +a1407*k7+a1408*k8+a1409*k9+a1410*k10+a1411*k11+a1412*k12+a1413*k13),k14)
        f(t+c15*Δt,u+Δt*(a1501*k1+a1506*k6+a1507*k7+a1508*k8                   +a1511*k11+a1512*k12+a1513*k13+a1514*k14),k15)
        f(t+c16*Δt,u+Δt*(a1601*k1+a1606*k6+a1607*k7+a1608*k8+a1609*k9                              +a1613*k13+a1614*k14+a1615*k15),k16)
        udiff = kupdate
        k[1] = udiff
        bspl = k1 - udiff
        k[2] = bspl
        k[3] = udiff - k13 - bspl
        k[4] = (d401*k1+d406*k6+d407*k7+d408*k8+d409*k9+d410*k10+d411*k11+d412*k12+d413*k13+d414*k14+d415*k15+d416*k16)
        k[5] = (d501*k1+d506*k6+d507*k7+d508*k8+d509*k9+d510*k10+d511*k11+d512*k12+d513*k13+d514*k14+d515*k15+d516*k16)
        k[6] = (d601*k1+d606*k6+d607*k7+d608*k8+d609*k9+d610*k10+d611*k11+d612*k12+d613*k13+d614*k14+d615*k15+d616*k16)
        k[7] = (d701*k1+d706*k6+d707*k7+d708*k8+d709*k9+d710*k10+d711*k11+d712*k12+d713*k13+d714*k14+d715*k15+d716*k16)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end


function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:TsitPap8Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13 = constructTsitPap8(uEltypeNoUnits)
  k1 = rateType(sizeu); k2 = rateType(sizeu); k3 = rateType(sizeu); k4 = rateType(sizeu)
  k5 = rateType(sizeu); k6 = rateType(sizeu); k7 = rateType(sizeu); k8 = rateType(sizeu)
  k9 = rateType(sizeu); k10 = rateType(sizeu); k11 = rateType(sizeu); k12 = rateType(sizeu)
  k13::rateType = rateType(sizeu); utilde = similar(u);
  k = rateType(sizeu)
  if calcprevs
    kprev = rateType(sizeu)
  end
  if calck
    pop!(ks)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k); k1=k
      f(t+c1*Δt,u+Δt*(a0201*k1),k2)
      f(t+c2*Δt,u+Δt*(a0301*k1+a0302*k2),k3)
      f(t+c3*Δt,u+Δt*(a0401*k1       +a0403*k3),k4)
      f(t+c4*Δt,u+Δt*(a0501*k1       +a0503*k3+a0504*k4),k5)
      f(t+c5*Δt,u+Δt*(a0601*k1                +a0604*k4+a0605*k5),k6)
      f(t+c6*Δt,u+Δt*(a0701*k1                +a0704*k4+a0705*k5+a0706*k6),k7)
      f(t+c7*Δt,u+Δt*(a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7),k8)
      f(t+c8*Δt,u+Δt*(a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8),k9)
      f(t+c9*Δt,u+Δt*(a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9),k10)
      f(t+c10*Δt,u+Δt*(a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10),k11)
      f(t+Δt,u+Δt*(a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11),k12)
      f(t+Δt,u+Δt*(a1301*k1                +a1304*k4+a1305*k5+a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10),k13)
      update = Δt*(b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12)
      utmp = u + update
      if adaptive
        EEst = sqrt(sum(((update - Δt*(k1*bhat1 + k6*bhat6 + k7*bhat7 + k8*bhat8 + k9*bhat9 + k10*bhat10 + k13*bhat13))./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor) # Order 5
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  if calck
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end



function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Vern9Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0806,a0807,a0901,a0906,a0907,a0908,a1001,a1006,a1007,a1008,a1009,a1101,a1106,a1107,a1108,a1109,a1110,a1201,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1401,a1406,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,b1,b8,b9,b10,b11,b12,b13,b14,b15,bhat1,bhat8,bhat9,bhat10,bhat11,bhat12,bhat13,bhat16 = constructVern9(uEltypeNoUnits)
  k1 = rateType(sizeu); k2 = rateType(sizeu);k3 = rateType(sizeu); k4 = rateType(sizeu);
  k5 = rateType(sizeu); k6 = rateType(sizeu);k7 = rateType(sizeu); k8 = rateType(sizeu);
  k9 = rateType(sizeu); k10 = rateType(sizeu); k11 = rateType(sizeu); k12 = rateType(sizeu);
  k13 = rateType(sizeu); k14 = rateType(sizeu); k15 = rateType(sizeu); k16 =rateType(sizeu);
  utilde = similar(u); update = similar(u)
  const kshortsize = 16
  if calck
    k = ksEltype()
    for i in 1:kshortsize
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    if calcprevs
      kprev = deepcopy(k)
      for i in 1:3 # Make it full-sized
        push!(kprev,rateType(sizeu))
      end
    end
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k1)
      f(t+c1*Δt,u+Δt*(a0201*k1),k2)
      f(t+c2*Δt,u+Δt*(a0301*k1+a0302*k2),k3)
      f(t+c3*Δt,u+Δt*(a0401*k1       +a0403*k3),k4)
      f(t+c4*Δt,u+Δt*(a0501*k1       +a0503*k3+a0504*k4),k5)
      f(t+c5*Δt,u+Δt*(a0601*k1                +a0604*k4+a0605*k5),k6)
      f(t+c6*Δt,u+Δt*(a0701*k1                +a0704*k4+a0705*k5+a0706*k6),k7)
      f(t+c7*Δt,u+Δt*(a0801*k1                                  +a0806*k6+a0807*k7),k8)
      f(t+c8*Δt,u+Δt*(a0901*k1                                  +a0906*k6+a0907*k7+a0908*k8),k9)
      f(t+c9*Δt,u+Δt*(a1001*k1                                  +a1006*k6+a1007*k7+a1008*k8+a1009*k9),k10)
      f(t+c10*Δt,u+Δt*(a1101*k1                                  +a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10),k11)
      f(t+c11*Δt,u+Δt*(a1201*k1                                  +a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11),k12)
      f(t+c12*Δt,u+Δt*(a1301*k1                                  +a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10+a1311*k11+a1312*k12),k13)
      f(t+c13*Δt,u+Δt*(a1401*k1                                  +a1406*k6+a1407*k7+a1408*k8+a1409*k9+a1410*k10+a1411*k11+a1412*k12+a1413*k13),k14)
      f(t+Δt,u+Δt*(a1501*k1                                  +a1506*k6+a1507*k7+a1508*k8+a1509*k9+a1510*k10+a1511*k11+a1512*k12+a1513*k13+a1514*k14),k15)
      f(t+Δt,u+Δt*(a1601*k1                                  +a1606*k6+a1607*k7+a1608*k8+a1609*k9+a1610*k10+a1611*k11+a1612*k12+a1613*k13),k16)
      update = Δt*(k1*b1+k8*b8+k9*b9+k10*b10+k11*b11+k12*b12+k13*b13+k14*b14+k15*b15)
      utmp = u + update
      if adaptive
        EEst = sqrt(sum(((update - Δt*(k1*bhat1 + k8*bhat8 + k9*bhat9 + k10*bhat10 + k11*bhat11 + k12*bhat12 + k13*bhat13 + k16*bhat16))./(abstol+max(abs.(u),abs.(utmp))*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if calck
        k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10;k[11]=k11;k[12]=k12;k[13]=k13;k[14]=k14;k[15]=k15;k[16]=k16
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end


function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Feagin10Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16 = constructFeagin10(uEltypeNoUnits)
  k1 = rateType(sizeu); k2 = rateType(sizeu); k3 = rateType(sizeu); k4 = rateType(sizeu); k5 = rateType(sizeu)
  k6 = rateType(sizeu); k7 = rateType(sizeu); k8 = rateType(sizeu); k9 = rateType(sizeu); k10 = rateType(sizeu)
  k11 = rateType(sizeu); k12 = rateType(sizeu); k13 = rateType(sizeu); k14 = rateType(sizeu)
  k15 = rateType(sizeu); k16 = rateType(sizeu); k17 = rateType(sizeu);
  k = rateType(sizeu)
  if calcprevs
    kprev = rateType(sizeu)
  end
  if calck
    pop!(ks)
  end
  update = rateType(sizeu)
  utmp = similar(u)
  uidx = eachindex(u)

  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k); k1=k
      f(t + c1*Δt,u + Δt*(a0100*k1),k2)
      f(t + c2*Δt,u + Δt*(a0200*k1 + a0201*k2),k3)
      f(t + c3*Δt,u + Δt*(a0300*k1              + a0302*k3),k4)
      f(t + c4*Δt,u + Δt*(a0400*k1              + a0402*k3 + a0403*k4),k5)
      f(t + c5*Δt,u + Δt*(a0500*k1                           + a0503*k4 + a0504*k5),k6)
      f(t + c6*Δt,u + Δt*(a0600*k1                           + a0603*k4 + a0604*k5 + a0605*k6),k7)
      f(t + c7*Δt,u + Δt*(a0700*k1                                        + a0704*k5 + a0705*k6 + a0706*k7),k8)
      f(t + c8*Δt,u + Δt*(a0800*k1                                                     + a0805*k6 + a0806*k7 + a0807*k8),k9)
      f(t + c9*Δt,u + Δt*(a0900*k1                                                     + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9),k10)
      f(t + c10*Δt,u + Δt*(a1000*k1                                                     + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10),k11)
      f(t + c11*Δt,u + Δt*(a1100*k1                                                     + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11),k12)
      f(t + c12*Δt,u + Δt*(a1200*k1                           + a1203*k4 + a1204*k5 + a1205*k6 + a1206*k7 + a1207*k8 + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12),k13)
      f(t + c13*Δt,u + Δt*(a1300*k1              + a1302*k3 + a1303*k4              + a1305*k6 + a1306*k7 + a1307*k8 + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13),k14)
      f(t + c14*Δt,u + Δt*(a1400*k1 + a1401*k2                           + a1404*k5              + a1406*k7 +                                                                     a1412*k13 + a1413*k14),k15)
      f(t + c15*Δt,u + Δt*(a1500*k1              + a1502*k3                                                                                                                                                     + a1514*k15),k16)
      f(t + c16*Δt,u + Δt*(a1600*k1 + a1601*k2 + a1602*k3              + a1604*k5 + a1605*k6 + a1606*k7 + a1607*k8 + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16),k17)
      for i in uidx
        update[i] = (b1*k1[i] + b2*k2[i] + b3*k3[i] + b5*k5[i]) + (b7*k7[i] + b9*k9[i] + b10*k10[i] + b11*k11[i]) + (b12*k12[i] + b13*k13[i] + b14*k14[i] + b15*k15[i]) + (b16*k16[i] + b17*k17[i])
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + Δt*update[i]
        end
        EEst = norm((Δt*(k2 - k16) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting, so in-place
        for i in uidx
          u[i] = u[i] + Δt*update[i]
        end
      end
      @ode_loopfooter
    end
  end
  if calck
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end


function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Feagin12Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25 = constructFeagin12(uEltypeNoUnits)
  k1 = rateType(sizeu); k2 = rateType(sizeu); k3 = rateType(sizeu); k4 = rateType(sizeu); k5 = rateType(sizeu)
  k6 = rateType(sizeu); k7 = rateType(sizeu); k8 = rateType(sizeu); k9 = rateType(sizeu); k10 = rateType(sizeu)
  k11 = rateType(sizeu); k12 = rateType(sizeu); k13 = rateType(sizeu); k14 = rateType(sizeu)
  k15 = rateType(sizeu); k16 = rateType(sizeu); k17 = rateType(sizeu); k18 = rateType(sizeu)
  k19 = rateType(sizeu); k20 = rateType(sizeu); k21 = rateType(sizeu); k22 = rateType(sizeu)
  k23 = rateType(sizeu); k24 = rateType(sizeu); k25 = rateType(sizeu)
  update = similar(u);
  utmp = similar(u)
  uidx = eachindex(u)
  k = rateType(sizeu)
  if calcprevs
    kprev = rateType(sizeu)
  end
  if calck
    pop!(ks)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k); k1=k
      f(t + c1*Δt,u + Δt*(a0100*k1),k2)
      f(t + c2*Δt ,u + Δt*(a0200*k1 + a0201*k2),k3)
      f(t + c3*Δt,u + Δt*(a0300*k1              + a0302*k3),k4)
      f(t + c4*Δt,u + Δt*(a0400*k1              + a0402*k3 + a0403*k4),k5)
      f(t + c5*Δt,u + Δt*(a0500*k1                           + a0503*k4 + a0504*k5),k6)
      f(t + c6*Δt,u + Δt*(a0600*k1                           + a0603*k4 + a0604*k5 + a0605*k6),k7)
      f(t + c7*Δt,u + Δt*(a0700*k1                                        + a0704*k5 + a0705*k6 + a0706*k7),k8)
      f(t + c8*Δt,u + Δt*(a0800*k1                                                     + a0805*k6 + a0806*k7 + a0807*k8),k9)
      f(t + c9*Δt,u + Δt*(a0900*k1                                                     + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9),k10)
      f(t + c10*Δt,u + Δt*(a1000*k1                                                     + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10),k11)
      f(t + c11*Δt,u + Δt*(a1100*k1                                                     + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11),k12)
      f(t + c12*Δt,u + Δt*(a1200*k1                                                                                            + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12),k13)
      f(t + c13*Δt,u + Δt*(a1300*k1                                                                                            + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13),k14)
      f(t + c14*Δt,u + Δt*(a1400*k1                                                                                            + a1408*k9 + a1409*k10 + a1410*k11 + a1411*k12 + a1412*k13 + a1413*k14),k15)
      f(t + c15*Δt,u + Δt*(a1500*k1                                                                                            + a1508*k9 + a1509*k10 + a1510*k11 + a1511*k12 + a1512*k13 + a1513*k14 + a1514*k15),k16)
      f(t + c16*Δt,u + Δt*(a1600*k1                                                                                            + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16),k17)
      f(t + c17*Δt,u + Δt*(a1700*k1                                                     + a1705*k6 + a1706*k7 + a1707*k8 + a1708*k9 + a1709*k10 + a1710*k11 + a1711*k12 + a1712*k13 + a1713*k14 + a1714*k15 + a1715*k16 + a1716*k17),k18)
      f(t + c18*Δt,u + Δt*(a1800*k1                                                     + a1805*k6 + a1806*k7 + a1807*k8 + a1808*k9 + a1809*k10 + a1810*k11 + a1811*k12 + a1812*k13 + a1813*k14 + a1814*k15 + a1815*k16 + a1816*k17 + a1817*k18),k19)
      f(t + c19*Δt,u + Δt*(a1900*k1                                        + a1904*k5 + a1905*k6 + a1906*k7              + a1908*k9 + a1909*k10 + a1910*k11 + a1911*k12 + a1912*k13 + a1913*k14 + a1914*k15 + a1915*k16 + a1916*k17 + a1917*k18 + a1918*k19),k20)
      f(t + c20*Δt,u + Δt*(a2000*k1                           + a2003*k4 + a2004*k5 + a2005*k6              + a2007*k8              + a2009*k10 + a2010*k11                                                                                     + a2017*k18 + a2018*k19 + a2019*k20),k21)
      f(t + c21*Δt,u + Δt*(a2100*k1              + a2102*k3 + a2103*k4                           + a2106*k7 + a2107*k8              + a2109*k10 + a2110*k11                                                                                     + a2117*k18 + a2118*k19 + a2119*k20 + a2120*k21),k22)
      f(t + c22*Δt,u + Δt*(a2200*k1 + a2201*k2                           + a2204*k5              + a2206*k7                                                                                                                                                                                     + a2220*k21 + a2221*k22),k23)
      f(t + c23*Δt,u + Δt*(a2300*k1              + a2302*k3                                                                                                                                                                                                                                                                     + a2322*k23),k24)
      f(t + c24*Δt,u + Δt*(a2400*k1 + a2401*k2 + a2402*k3              + a2404*k5              + a2406*k7 + a2407*k8 + a2408*k9 + a2409*k10 + a2410*k11 + a2411*k12 + a2412*k13 + a2413*k14 + a2414*k15 + a2415*k16 + a2416*k17 + a2417*k18 + a2418*k19 + a2419*k20 + a2420*k21 + a2421*k22 + a2422*k23 + a2423*k24),k25)

      for i in uidx
        update[i] = Δt*((b1*k1[i] + b2*k2[i] + b3*k3[i] + b5*k5[i]) + (b7*k7[i] + b8*k8[i] + b10*k10[i] + b11*k11[i]) + (b13*k13[i] + b14*k14[i] + b15*k15[i] + b16*k16[i]) + (b17*k17[i] + b18*k18[i] + b19*k19[i] + b20*k20[i]) + (b21*k21[i] + b22*k22[i] + b23*k23[i] + b24*k24[i]) + b25*k25[i])
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + update[i]
        end
        EEst = norm((Δt*(k2 - k24) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting so in-place
        for i in uidx
          u[i] = u[i] + update[i]
        end
      end
      @ode_loopfooter
    end
  end
  if calck
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end


function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Feagin14Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,b35 = constructFeagin14(uEltypeNoUnits)
  k1 = rateType(sizeu); k2 = rateType(sizeu); k3 = rateType(sizeu); k4 = rateType(sizeu); k5 = rateType(sizeu)
  k6 = rateType(sizeu); k7 = rateType(sizeu); k8 = rateType(sizeu); k9 = rateType(sizeu); k10 = rateType(sizeu)
  k11 = rateType(sizeu); k12 = rateType(sizeu); k13 = rateType(sizeu); k14 = rateType(sizeu)
  k15 = rateType(sizeu); k16 = rateType(sizeu); k17 = rateType(sizeu); k18 = rateType(sizeu)
  k19 = rateType(sizeu); k20 = rateType(sizeu); k21 = rateType(sizeu); k22 = rateType(sizeu)
  k23 = rateType(sizeu); k24 = rateType(sizeu); k25 = rateType(sizeu)
  k26 = rateType(sizeu); k27 = rateType(sizeu); k28 = rateType(sizeu)
  k29 = rateType(sizeu); k30 = rateType(sizeu); k31 = rateType(sizeu); k32 = rateType(sizeu)
  k33 = rateType(sizeu); k34 = rateType(sizeu); k35 = rateType(sizeu)
  update = similar(u)
  utmp = similar(u)
  uidx = eachindex(u)

  k = rateType(sizeu)
  if calcprevs
    kprev = rateType(sizeu)
  end
  if calck
    pop!(ks)
  end
  k = k1
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k1)
      f(t + c1*Δt,u + Δt*(a0100*k1),k2)
      f(t + c2*Δt ,u + Δt*(a0200*k1 + a0201*k2),k3)
      f(t + c3*Δt,u + Δt*(a0300*k1              + a0302*k3),k4)
      f(t + c4*Δt,u + Δt*(a0400*k1              + a0402*k3 + a0403*k4),k5)
      f(t + c5*Δt,u + Δt*(a0500*k1                           + a0503*k4 + a0504*k5),k6)
      f(t + c6*Δt,u + Δt*(a0600*k1                           + a0603*k4 + a0604*k5 + a0605*k6),k7)
      f(t + c7*Δt,u + Δt*(a0700*k1                                        + a0704*k5 + a0705*k6 + a0706*k7),k8)
      f(t + c8*Δt,u + Δt*(a0800*k1                                                     + a0805*k6 + a0806*k7 + a0807*k8),k9)
      f(t + c9*Δt,u + Δt*(a0900*k1                                                     + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9),k10)
      f(t + c10*Δt,u + Δt*(a1000*k1                                                     + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10),k11)
      f(t + c11*Δt,u + Δt*(a1100*k1                                                     + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11),k12)
      f(t + c12*Δt,u + Δt*(a1200*k1                                                                                            + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12),k13)
      f(t + c13*Δt,u + Δt*(a1300*k1                                                                                            + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13),k14)
      f(t + c14*Δt,u + Δt*(a1400*k1                                                                                            + a1408*k9 + a1409*k10 + a1410*k11 + a1411*k12 + a1412*k13 + a1413*k14),k15)
      f(t + c15*Δt,u + Δt*(a1500*k1                                                                                            + a1508*k9 + a1509*k10 + a1510*k11 + a1511*k12 + a1512*k13 + a1513*k14 + a1514*k15),k16)
      f(t + c16*Δt,u + Δt*(a1600*k1                                                                                            + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16),k17)
      f(t + c17*Δt,u + Δt*(a1700*k1                                                                                                                                                   + a1712*k13 + a1713*k14 + a1714*k15 + a1715*k16 + a1716*k17),k18)
      f(t + c18*Δt,u + Δt*(a1800*k1                                                                                                                                                   + a1812*k13 + a1813*k14 + a1814*k15 + a1815*k16 + a1816*k17 + a1817*k18),k19)
      f(t + c19*Δt,u + Δt*(a1900*k1                                                                                                                                                   + a1912*k13 + a1913*k14 + a1914*k15 + a1915*k16 + a1916*k17 + a1917*k18 + a1918*k19),k20)
      f(t + c20*Δt,u + Δt*(a2000*k1                                                                                                                                                   + a2012*k13 + a2013*k14 + a2014*k15 + a2015*k16 + a2016*k17 + a2017*k18 + a2018*k19 + a2019*k20),k21)
      f(t + c21*Δt,u + Δt*(a2100*k1                                                                                                                                                   + a2112*k13 + a2113*k14 + a2114*k15 + a2115*k16 + a2116*k17 + a2117*k18 + a2118*k19 + a2119*k20 + a2120*k21),k22)
      f(t + c22*Δt,u + Δt*(a2200*k1                                                                                                                                                   + a2212*k13 + a2213*k14 + a2214*k15 + a2215*k16 + a2216*k17 + a2217*k18 + a2218*k19 + a2219*k20 + a2220*k21 + a2221*k22),k23)
      f(t + c23*Δt,u + Δt*(a2300*k1                                                                                            + a2308*k9 + a2309*k10 + a2310*k11 + a2311*k12 + a2312*k13 + a2313*k14 + a2314*k15 + a2315*k16 + a2316*k17 + a2317*k18 + a2318*k19 + a2319*k20 + a2320*k21 + a2321*k22 + a2322*k23),k24)
      f(t + c24*Δt,u + Δt*(a2400*k1                                                                                            + a2408*k9 + a2409*k10 + a2410*k11 + a2411*k12 + a2412*k13 + a2413*k14 + a2414*k15 + a2415*k16 + a2416*k17 + a2417*k18 + a2418*k19 + a2419*k20 + a2420*k21 + a2421*k22 + a2422*k23 + a2423*k24),k25)
      f(t + c25*Δt,u + Δt*(a2500*k1                                                                                            + a2508*k9 + a2509*k10 + a2510*k11 + a2511*k12 + a2512*k13 + a2513*k14 + a2514*k15 + a2515*k16 + a2516*k17 + a2517*k18 + a2518*k19 + a2519*k20 + a2520*k21 + a2521*k22 + a2522*k23 + a2523*k24 + a2524*k25),k26)
      f(t + c26*Δt,u + Δt*(a2600*k1                                                     + a2605*k6 + a2606*k7 + a2607*k8 + a2608*k9 + a2609*k10 + a2610*k11               + a2612*k13 + a2613*k14 + a2614*k15 + a2615*k16 + a2616*k17 + a2617*k18 + a2618*k19 + a2619*k20 + a2620*k21 + a2621*k22 + a2622*k23 + a2623*k24 + a2624*k25 + a2625*k26),k27)
      f(t + c27*Δt,u + Δt*(a2700*k1                                                     + a2705*k6 + a2706*k7 + a2707*k8 + a2708*k9 + a2709*k10               + a2711*k12 + a2712*k13 + a2713*k14 + a2714*k15 + a2715*k16 + a2716*k17 + a2717*k18 + a2718*k19 + a2719*k20 + a2720*k21 + a2721*k22 + a2722*k23 + a2723*k24 + a2724*k25 + a2725*k26 + a2726*k27),k28)
      f(t + c28*Δt,u + Δt*(a2800*k1                                                     + a2805*k6 + a2806*k7 + a2807*k8 + a2808*k9               + a2810*k11 + a2811*k12               + a2813*k14 + a2814*k15 + a2815*k16                                                                                                   + a2823*k24 + a2824*k25 + a2825*k26 + a2826*k27 + a2827*k28),k29)
      f(t + c29*Δt,u + Δt*(a2900*k1                                        + a2904*k5 + a2905*k6 + a2906*k7                           + a2909*k10 + a2910*k11 + a2911*k12               + a2913*k14 + a2914*k15 + a2915*k16                                                                                                   + a2923*k24 + a2924*k25 + a2925*k26 + a2926*k27 + a2927*k28 + a2928*k29),k30)
      f(t + c30*Δt,u + Δt*(a3000*k1                           + a3003*k4 + a3004*k5 + a3005*k6              + a3007*k8              + a3009*k10 + a3010*k11                             + a3013*k14 + a3014*k15 + a3015*k16                                                                                                   + a3023*k24 + a3024*k25 + a3025*k26               + a3027*k28 + a3028*k29 + a3029*k30),k31)
      f(t + c31*Δt,u + Δt*(a3100*k1              + a3102*k3 + a3103*k4                           + a3106*k7 + a3107*k8              + a3109*k10 + a3110*k11                             + a3113*k14 + a3114*k15 + a3115*k16                                                                                                   + a3123*k24 + a3124*k25 + a3125*k26               + a3127*k28 + a3128*k29 + a3129*k30 + a3130*k31),k32)
      f(t + c32*Δt,u + Δt*(a3200*k1 + a3201*k2                           + a3204*k5              + a3206*k7                                                                                                                                                                                                                                                                                                                                 + a3230*k31 + a3231*k32),k33)
      f(t + c33*Δt,u + Δt*(a3300*k1              + a3302*k3                                                                                                                                                                                                                                                                                                                                                                                                                 + a3332*k33),k34)
      f(t + c34*Δt,u + Δt*(a3400*k1 + a3401*k2 + a3402*k3              + a3404*k5              + a3406*k7 + a3407*k8              + a3409*k10 + a3410*k11 + a3411*k12 + a3412*k13 + a3413*k14 + a3414*k15 + a3415*k16 + a3416*k17 + a3417*k18 + a3418*k19 + a3419*k20 + a3420*k21 + a3421*k22 + a3422*k23 + a3423*k24 + a3424*k25 + a3425*k26 + a3426*k27 + a3427*k28 + a3428*k29 + a3429*k30 + a3430*k31 + a3431*k32 + a3432*k33 + a3433*k34),k35)
      for i in uidx
        update[i] = Δt*((b1*k1[i] + b2*k2[i] + b3*k3[i] + b5*k5[i]) + (b7*k7[i] + b8*k8[i] + b10*k10[i] + b11*k11[i]) + (b12*k12[i] + b14*k14[i] + b15*k15[i] + b16*k16[i]) + (b18*k18[i] + b19*k19[i] + b20*k20[i] + b21*k21[i]) + (b22*k22[i] + b23*k23[i] + b24*k24[i] + b25*k25[i]) + (b26*k26[i] + b27*k27[i] + b28*k28[i] + b29*k29[i]) + (b30*k30[i] + b31*k31[i] + b32*k32[i] + b33*k33[i]) + (b34*k34[i] + b35*k35[i]))
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + update[i]
        end
        EEst = norm((Δt*(k2 - k34) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting, so in-place
        for i in uidx
          u[i] = u[i] + update[i]
        end
      end
      @ode_loopfooter
    end
  end
  if calck
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end


"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 192
"""
function ode_interpolant(Θ,Δt,y₀,y₁,kprevious,k,T::Type{Val{:DP5Vectorized}})
  Θ1 = 1-Θ
  y₀ + Δt*Θ*(k[1]+Θ1*(k[2]+Θ*(k[3]+Θ1*k[4])))
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
  y₀ + Δt*(k[1]*b1Θ + k[2]*b2Θ + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ)
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
  b1Θ =           r012*Θ2 + r013*Θ3 + r014*Θ4 + r015*Θ5 + r016*Θ6
  b3Θ =           r032*Θ2 + r033*Θ3 + r034*Θ4 + r035*Θ5 + r036*Θ6
  b4Θ =           r042*Θ2 + r043*Θ3 + r044*Θ4 + r045*Θ5 + r046*Θ6
  b5Θ =           r052*Θ2 + r053*Θ3 + r054*Θ4 + r055*Θ5 + r056*Θ6
  b6Θ =           r062*Θ2 + r063*Θ3 + r064*Θ4 + r065*Θ5 + r066*Θ6
  b7Θ =           r072*Θ2 + r073*Θ3 + r074*Θ4 + r075*Θ5 + r076*Θ6
  b8Θ =           r082*Θ2 + r083*Θ3 + r084*Θ4 + r085*Θ5 + r086*Θ6
  b9Θ =           r092*Θ2 + r093*Θ3 + r094*Θ4 + r095*Θ5 + r096*Θ6
  b10Θ=           r102*Θ2 + r103*Θ3 + r104*Θ4 + r105*Θ5 + r106*Θ6
  b11Θ=           r112*Θ2 + r113*Θ3 + r114*Θ4 + r115*Θ5 + r116*Θ6
  y₀ + Δt*(k[1]*b1Θ + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ)
end


"""

"""
function ode_interpolant(Θ,Δt,y₀,y₁,kprevious,k,T::Type{Val{:DP8Vectorized}})
  Θ1 = 1-Θ
  conpar = k[4] + Θ*(k[5] + Θ1*(k[6]+Θ*k[7]))
  y₀ + Δt*Θ*(k[1] + Θ1*(k[2] + Θ*(k[3]+Θ1*conpar)))
end



"""
An Efficient Runge-Kutta (4,5) Pair by P.Bogacki and L.F.Shampine
 Computers and Mathematics with Applications, Vol. 32, No. 6, 1996, pages 15 to 28

Called to add the extra k9, k10, k11 steps for the Order 5 interpolation when needed
"""
function ode_addsteps!{rateType<:AbstractArray,uEltypeNoUnits}(k,t,u,Δt,f,T::Type{Val{:BS5Vectorized}},T2::Type{rateType},T3::Type{uEltypeNoUnits})
  if length(k) < 8
    c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = constructBS5(uEltypeNoUnits)
    rtmp = rateType(size(u))
    f(t,u,rtmp); push!(k,copy(rtmp))
    f(t+c1*Δt,u+Δt*a21*k[1],rtmp); push!(k,copy(rtmp))
    f(t+c2*Δt,u+Δt*(a31*k[1]+a32*k[2]),rtmp); push!(k,copy(rtmp))
    f(t+c3*Δt,u+Δt*(a41*k[1]+a42*k[2]+a43*k[3]),rtmp); push!(k,copy(rtmp))
    f(t+c4*Δt,u+Δt*(a51*k[1]+a52*k[2]+a53*k[3]+a54*k[4]),rtmp); push!(k,copy(rtmp))
    f(t+c5*Δt,u+Δt*(a61*k[1]+a62*k[2]+a63*k[3]+a64*k[4]+a65*k[5]),rtmp); push!(k,copy(rtmp))
    f(t+Δt,u+Δt*(a71*k[1]+a72*k[2]+a73*k[3]+a74*k[4]+a75*k[5]+a76*k[6]),rtmp); push!(k,copy(rtmp))
    f(t+Δt,u+Δt*(a81*k[1]+a83*k[3]+a84*k[4]+a85*k[5]+a86*k[6]+a87*k[7]),rtmp); push!(k,copy(rtmp))
  end
  if length(k) < 11 # Have not added the extra stages yet
    c6,c7,c8,a91,a92,a93,a94,a95,a96,a97,a98,a101,a102,a103,a104,a105,a106,a107,a108,a109,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110 = BS5Interp(uEltypeNoUnits)
    rtmp = rateType(size(u))
    f(t+c6*Δt,u+Δt*(a91*k[1]+a92*k[2]+a93*k[3]+a94*k[4]+a95*k[5]+a96*k[6]+a97*k[7]+a98*k[8]),rtmp); push!(k,copy(rtmp))
    f(t+c7*Δt,u+Δt*(a101*k[1]+a102*k[2]+a103*k[3]+a104*k[4]+a105*k[5]+a106*k[6]+a107*k[7]+a108*k[8]+a109*k[9]),rtmp); push!(k,copy(rtmp))
    f(t+c8*Δt,u+Δt*(a111*k[1]+a112*k[2]+a113*k[3]+a114*k[4]+a115*k[5]+a116*k[6]+a117*k[7]+a118*k[8]+a119*k[9]+a1110*k[10]),rtmp); push!(k,copy(rtmp))
  end
  nothing
end


"""

"""
function ode_addsteps!{rateType<:AbstractArray,uEltypeNoUnits}(k,t,u,Δt,f,T::Type{Val{:Vern6Vectorized}},T2::Type{rateType},T3::Type{uEltypeNoUnits})
  if length(k) < 9
    c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,b1,b4,b5,b6,b7,b8,b9= constructVern6(uEltypeNoUnits)
    rtmp = rateType(size(u))
    f(t,u,rtmp); push!(k,copy(rtmp))
    f(t+c1*Δt,u+Δt*(a21*k[1]),rtmp); push!(k,copy(rtmp))
    f(t+c2*Δt,u+Δt*(a31*k[1]+a32*k[2]),rtmp); push!(k,copy(rtmp))
    f(t+c3*Δt,u+Δt*(a41*k[1]       +a43*k[3]),rtmp); push!(k,copy(rtmp))
    f(t+c4*Δt,u+Δt*(a51*k[1]       +a53*k[3]+a54*k[4]),rtmp); push!(k,copy(rtmp))
    f(t+c5*Δt,u+Δt*(a61*k[1]       +a63*k[3]+a64*k[4]+a65*k[5]),rtmp); push!(k,copy(rtmp))
    f(t+c6*Δt,u+Δt*(a71*k[1]       +a73*k[3]+a74*k[4]+a75*k[5]+a76*k[6]),rtmp); push!(k,copy(rtmp))
    f(t+Δt,u+Δt*(a81*k[1]       +a83*k[3]+a84*k[4]+a85*k[5]+a86*k[6]+a87*k[7]),rtmp); push!(k,copy(rtmp))
    f(t+Δt,u+Δt*(a91*k[1]              +a94*k[4]+a95*k[5]+a96*k[6]+a97*k[7]+a98*k[8]),rtmp); push!(k,copy(rtmp))
  end
  if length(k) < 12 # Have not added the extra stages yet
    c10,a1001,a1004,a1005,a1006,a1007,a1008,a1009,c11,a1101,a1102,a1103,a1104,a1105,a1106,a1107,a1108,a1109,a1110,c12,a1201,a1202,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = Vern6Interp(uEltypeNoUnits)
    rtmp = rateType(size(u))
    f(t+c10*Δt,u+Δt*(a1001*k[1]+a1004*k[4]+a1005*k[5]+a1006*k[6]+a1007*k[7]+a1008*k[8]+a1009*k[9]),rtmp); push!(k,copy(rtmp))
    f(t+c11*Δt,u+Δt*(a1101*k[1]+a1102*k[2]+a1103*k[3]+a1104*k[4]+a1105*k[5]+a1106*k[6]+a1107*k[7]+a1108*k[8]+a1109*k[9]+a1110*k[10]),rtmp); push!(k,copy(rtmp))
    f(t+c12*Δt,u+Δt*(a1201*k[1]+a1202*k[2]+a1203*k[3]+a1204*k[4]+a1205*k[5]+a1206*k[6]+a1207*k[7]+a1208*k[8]+a1209*k[9]+a1210*k[10]+a1211*k[11]),rtmp); push!(k,copy(rtmp))
  end
  nothing
end


"""

"""
function ode_addsteps!{rateType<:AbstractArray,uEltypeNoUnits}(k,t,u,Δt,f,T::Type{Val{:Vern7Vectorized}},T2::Type{rateType},T3::Type{uEltypeNoUnits})
  if length(k) < 10
    c2,c3,c4,c5,c6,c7,c8,a021,a031,a032,a041,a043,a051,a053,a054,a061,a063,a064,a065,a071,a073,a074,a075,a076,a081,a083,a084,a085,a086,a087,a091,a093,a094,a095,a096,a097,a098,a101,a103,a104,a105,a106,a107,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat10= constructVern7(uEltypeNoUnit)
    rtmp = rateType(size(u))
    f(t,u,rtmp); push!(k,copy(rtmp))
    f(t+c2*Δt,u+Δt*(a021*k[1]),rtmp); push!(k,copy(rtmp))
    f(t+c3*Δt,u+Δt*(a031*k[1]+a032*k[2]),rtmp); push!(k,copy(rtmp))
    f(t+c4*Δt,u+Δt*(a041*k[1]       +a043*k[3]),rtmp); push!(k,copy(rtmp))
    f(t+c5*Δt,u+Δt*(a051*k[1]       +a053*k[3]+a054*k[4]),rtmp); push!(k,copy(rtmp))
    f(t+c6*Δt,u+Δt*(a061*k[1]       +a063*k[3]+a064*k[4]+a065*k[5]),rtmp); push!(k,copy(rtmp))
    f(t+c7*Δt,u+Δt*(a071*k[1]       +a073*k[3]+a074*k[4]+a075*k[5]+a076*k[6]),rtmp); push!(k,copy(rtmp))
    f(t+c8*Δt,u+Δt*(a081*k[1]       +a083*k[3]+a084*k[4]+a085*k[5]+a086*k[6]+a087*k[7]),rtmp); push!(k,copy(rtmp))
    f(t+Δt,u+Δt*(a091*k[1]          +a093*k[3]+a094*k[4]+a095*k[5]+a096*k[6]+a097*k[7]+a098*k[8]),rtmp); push!(k,copy(rtmp))
    f(t+Δt,u+Δt*(a101*k[1]          +a103*k[3]+a104*k[4]+a105*k[5]+a106*k[6]+a107*k[7]),rtmp); push!(k,copy(rtmp))
  end
  if length(k) < 16 # Have not added the extra stages yet
    c11,a1101,a1104,a1105,a1106,a1107,a1108,a1109,c12,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1211,c13,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1311,a1312,c14,a1401,a1404,a1405,a1406,a1407,a1408,a1409,a1411,a1412,a1413,c15,a1501,a1504,a1505,a1506,a1507,a1508,a1509,a1511,a1512,a1513,c16,a1601,a1604,a1605,a1606,a1607,a1608,a1609,a1611,a1612,a1613 = Vern7Interp(uEltypeNoUnits)
    rtmp = rateType(size(u))
    f(t+c11*Δt,u+Δt*(a1101*k[1]+a1104*k[4]+a1105*k[5]+a1106*k[6]+a1107*k[7]+a1108*k[8]+a1109*k[9]),rtmp); push!(k,copy(rtmp))
    f(t+c12*Δt,u+Δt*(a1201*k[1]+a1204*k[4]+a1205*k[5]+a1206*k[6]+a1207*k[7]+a1208*k[8]+a1209*k[9]+a1211*k[11]),rtmp); push!(k,copy(rtmp))
    f(t+c13*Δt,u+Δt*(a1301*k[1]+a1304*k[4]+a1305*k[5]+a1306*k[6]+a1307*k[7]+a1308*k[8]+a1309*k[9]+a1311*k[11]+a1312*k[12]),rtmp); push!(k,copy(rtmp))
    f(t+c14*Δt,u+Δt*(a1401*k[1]+a1404*k[4]+a1405*k[5]+a1406*k[6]+a1407*k[7]+a1408*k[8]+a1409*k[9]+a1411*k[11]+a1412*k[12]+a1413*k[13]),rtmp); push!(k,copy(rtmp))
    f(t+c15*Δt,u+Δt*(a1501*k[1]+a1504*k[4]+a1505*k[5]+a1506*k[6]+a1507*k[7]+a1508*k[8]+a1509*k[9]+a1511*k[11]+a1512*k[12]+a1513*k[13]),rtmp); push!(k,copy(rtmp))
    f(t+c16*Δt,u+Δt*(a1601*k[1]+a1604*k[4]+a1605*k[5]+a1606*k[6]+a1607*k[7]+a1608*k[8]+a1609*k[9]+a1611*k[11]+a1612*k[12]+a1613*k[13]),rtmp); push!(k,copy(rtmp))
  end
  nothing
end


"""

"""
function ode_addsteps!{rateType<:AbstractArray,uEltypeNoUnits}(k,t,u,Δt,f,T::Type{Val{:Vern8Vectorized}},T2::Type{rateType},T3::Type{uEltypeNoUnits})
  if length(k) <13
    c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13= constructVern8(uEltypeNoUnits)
    rtmp = rateType(size(u))
    f(t,u,rtmp); push!(k,copy(rtmp))
    f(t+c2*Δt ,u+Δt*(a0201*k[1]),rtmp); push!(k,copy(rtmp))
    f(t+c3*Δt ,u+Δt*(a0301*k[1]+a0302*k[2]),rtmp); push!(k,copy(rtmp))
    f(t+c4*Δt ,u+Δt*(a0401*k[1]       +a0403*k[3]),rtmp); push!(k,copy(rtmp))
    f(t+c5*Δt ,u+Δt*(a0501*k[1]       +a0503*k[3]+a0504*k[4]),rtmp); push!(k,copy(rtmp))
    f(t+c6*Δt ,u+Δt*(a0601*k[1]                +a0604*k[4]+a0605*k[5]),rtmp); push!(k,copy(rtmp))
    f(t+c7*Δt ,u+Δt*(a0701*k[1]                +a0704*k[4]+a0705*k[5]+a0706*k[6]),rtmp); push!(k,copy(rtmp))
    f(t+c8*Δt ,u+Δt*(a0801*k[1]                +a0804*k[4]+a0805*k[5]+a0806*k[6]+a0807*k[7]),rtmp); push!(k,copy(rtmp))
    f(t+c9*Δt ,u+Δt*(a0901*k[1]                +a0904*k[4]+a0905*k[5]+a0906*k[6]+a0907*k[7]+a0908*k[8]),rtmp); push!(k,copy(rtmp))
    f(t+c10*Δt,u+Δt*(a1001*k[1]                +a1004*k[4]+a1005*k[5]+a1006*k[6]+a1007*k[7]+a1008*k[8]+a1009*k[9]),rtmp); push!(k,copy(rtmp))
    f(t+c11*Δt,u+Δt*(a1101*k[1]                +a1104*k[4]+a1105*k[5]+a1106*k[6]+a1107*k[7]+a1108*k[8]+a1109*k[9]+a1110*k[10]),rtmp); push!(k,copy(rtmp))
    f(t+    Δt,u+Δt*(a1201*k[1]                +a1204*k[4]+a1205*k[5]+a1206*k[6]+a1207*k[7]+a1208*k[8]+a1209*k[9]+a1210*k[10]+a1211*k[11]),rtmp); push!(k,copy(rtmp))
    f(t+    Δt,u+Δt*(a1301*k[1]                +a1304*k[4]+a1305*k[5]+a1306*k[6]+a1307*k[7]+a1308*k[8]+a1309*k[9]+a1310*k[10]),rtmp); push!(k,copy(rtmp))
  end
  if length(k) < 26 # Have not added the extra stages yet
    c14,a1401,a1406,a1407,a1408,a1409,a1410,a1411,a1412,c15,a1501,a1506,a1507,a1508,a1509,a1510,a1511,a1512,a1514,c16,a1601,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1614,a1615,c17,a1701,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1714,a1715,a1716,c18,a1801,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1814,a1815,a1816,a1817,c19,a1901,a1906,a1907,a1908,a1909,a1910,a1911,a1912,a1914,a1915,a1916,a1917,c20,a2001,a2006,a2007,a2008,a2009,a2010,a2011,a2012,a2014,a2015,a2016,a2017,c21,a2101,a2106,a2107,a2108,a2109,a2110,a2111,a2112,a2114,a2115,a2116,a2117 = Vern8Interp(uEltypeNoUnits)
    rtmp = rateType(size(u))
    f(t+c14*Δt,u+Δt*(a1401*k[1]+a1406*k[6]+a1407*k[7]+a1408*k[8]+a1409*k[9]+a1410*k[10]+a1411*k[11]+a1412*k[12]),rtmp); push!(k,copy(rtmp))
    f(t+c15*Δt,u+Δt*(a1501*k[1]+a1506*k[6]+a1507*k[7]+a1508*k[8]+a1509*k[9]+a1510*k[10]+a1511*k[11]+a1512*k[12]+a1514*k[14]),rtmp); push!(k,copy(rtmp))
    f(t+c16*Δt,u+Δt*(a1601*k[1]+a1606*k[6]+a1607*k[7]+a1608*k[8]+a1609*k[9]+a1610*k[10]+a1611*k[11]+a1612*k[12]+a1614*k[14]+a1615*k[15]),rtmp); push!(k,copy(rtmp))
    f(t+c17*Δt,u+Δt*(a1701*k[1]+a1706*k[6]+a1707*k[7]+a1708*k[8]+a1709*k[9]+a1710*k[10]+a1711*k[11]+a1712*k[12]+a1714*k[14]+a1715*k[15]+a1716*k[16]),rtmp); push!(k,copy(rtmp))
    f(t+c18*Δt,u+Δt*(a1801*k[1]+a1806*k[6]+a1807*k[7]+a1808*k[8]+a1809*k[9]+a1810*k[10]+a1811*k[11]+a1812*k[12]+a1814*k[14]+a1815*k[15]+a1816*k[16]+a1817*k[17]),rtmp); push!(k,copy(rtmp))
    f(t+c19*Δt,u+Δt*(a1901*k[1]+a1906*k[6]+a1907*k[7]+a1908*k[8]+a1909*k[9]+a1910*k[10]+a1911*k[11]+a1912*k[12]+a1914*k[14]+a1915*k[15]+a1916*k[16]+a1917*k[17]),rtmp); push!(k,copy(rtmp))
    f(t+c20*Δt,u+Δt*(a2001*k[1]+a2006*k[6]+a2007*k[7]+a2008*k[8]+a2009*k[9]+a2010*k[10]+a2011*k[11]+a2012*k[12]+a2014*k[14]+a2015*k[15]+a2016*k[16]+a2017*k[17]),rtmp); push!(k,copy(rtmp))
    f(t+c21*Δt,u+Δt*(a2101*k[1]+a2106*k[6]+a2107*k[7]+a2108*k[8]+a2109*k[9]+a2110*k[10]+a2111*k[11]+a2112*k[12]+a2114*k[14]+a2115*k[15]+a2116*k[16]+a2117*k[17]),rtmp); push!(k,copy(rtmp))
  end
  nothing
end


"""

"""
function ode_addsteps!{rateType<:AbstractArray,uEltypeNoUnits}(k,t,u,Δt,f,T::Type{Val{:Vern9Vectorized}},T2::Type{rateType},T3::Type{uEltypeNoUnits})
  if length(k) < 16
    c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0806,a0807,a0901,a0906,a0907,a0908,a1001,a1006,a1007,a1008,a1009,a1101,a1106,a1107,a1108,a1109,a1110,a1201,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1401,a1406,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,b1,b8,b9,b10,b11,b12,b13,b14,b15,bhat1,bhat8,bhat9,bhat10,bhat11,bhat12,bhat13,bhat16 = constructVern9(uEltypeNoUnits)
    rtmp = rateType(size(u))
    f(t,u,rtmp); push!(k,copy(rtmp))
    f(t+c1*Δt,u+Δt*(a0201*k[1]),rtmp); push!(k,copy(rtmp))
    f(t+c2*Δt,u+Δt*(a0301*k[1]+a0302*k[2]),rtmp); push!(k,copy(rtmp))
    f(t+c3*Δt,u+Δt*(a0401*k[1]       +a0403*k[3]),rtmp); push!(k,copy(rtmp))
    f(t+c4*Δt,u+Δt*(a0501*k[1]       +a0503*k[3]+a0504*k[4]),rtmp); push!(k,copy(rtmp))
    f(t+c5*Δt,u+Δt*(a0601*k[1]                +a0604*k[4]+a0605*k[5]),rtmp); push!(k,copy(rtmp))
    f(t+c6*Δt,u+Δt*(a0701*k[1]                +a0704*k[4]+a0705*k[5]+a0706*k[6]),rtmp); push!(k,copy(rtmp))
    f(t+c7*Δt,u+Δt*(a0801*k[1]                                  +a0806*k[6]+a0807*k[7]),rtmp); push!(k,copy(rtmp))
    f(t+c8*Δt,u+Δt*(a0901*k[1]                                  +a0906*k[6]+a0907*k[7]+a0908*k[8]),rtmp); push!(k,copy(rtmp))
    f(t+c9*Δt,u+Δt*(a1001*k[1]                                  +a1006*k[6]+a1007*k[7]+a1008*k[8]+a1009*k[9]),rtmp); push!(k,copy(rtmp))
    f(t+c10*Δt,u+Δt*(a1101*k[1]                                  +a1106*k[6]+a1107*k[7]+a1108*k[8]+a1109*k[9]+a1110*k[10]),rtmp); push!(k,copy(rtmp))
    f(t+c11*Δt,u+Δt*(a1201*k[1]                                  +a1206*k[6]+a1207*k[7]+a1208*k[8]+a1209*k[9]+a1210*k[10]+a1211*k[11]),rtmp); push!(k,copy(rtmp))
    f(t+c12*Δt,u+Δt*(a1301*k[1]                                  +a1306*k[6]+a1307*k[7]+a1308*k[8]+a1309*k[9]+a1310*k[10]+a1311*k[11]+a1312*k[12]),rtmp); push!(k,copy(rtmp))
    f(t+c13*Δt,u+Δt*(a1401*k[1]                                  +a1406*k[6]+a1407*k[7]+a1408*k[8]+a1409*k[9]+a1410*k[10]+a1411*k[11]+a1412*k[12]+a1413*k[13]),rtmp); push!(k,copy(rtmp))
    f(t+Δt,u+Δt*(a1501*k[1]                                  +a1506*k[6]+a1507*k[7]+a1508*k[8]+a1509*k[9]+a1510*k[10]+a1511*k[11]+a1512*k[12]+a1513*k[13]+a1514*k[14]),rtmp); push!(k,copy(rtmp))
    f(t+Δt,u+Δt*(a1601*k[1]                                  +a1606*k[6]+a1607*k[7]+a1608*k[8]+a1609*k[9]+a1610*k[10]+a1611*k[11]+a1612*k[12]+a1613*k[13]),rtmp); push!(k,copy(rtmp))
  end
  if length(k) < 26 # Have not added the extra stages yet
    rtmp = rateType(size(u))
    c17,a1701,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,c18,a1801,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1817,c19,a1901,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1917,a1918,c20,a2001,a2008,a2009,a2010,a2011,a2012,a2013,a2014,a2015,a2017,a2018,a2019,c21,a2101,a2108,a2109,a2110,a2111,a2112,a2113,a2114,a2115,a2117,a2118,a2119,a2120,c22,a2201,a2208,a2209,a2210,a2211,a2212,a2213,a2214,a2215,a2217,a2218,a2219,a2220,a2221,c23,a2301,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2317,a2318,a2319,a2320,a2321,c24,a2401,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2417,a2418,a2419,a2420,a2421,c25,a2501,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2517,a2518,a2519,a2520,a2521,c26,a2601,a2608,a2609,a2610,a2611,a2612,a2613,a2614,a2615,a2617,a2618,a2619,a2620,a2621 = Vern9Interp(uEltypeNoUnits)
    f(t+c17*Δt,u+Δt*(a1701*k[1]+a1708*k[8]+a1709*k[9]+a1710*k[10]+a1711*k[11]+a1712*k[12]+a1713*k[13]+a1714*k[14]+a1715*k[15]),rtmp); push!(k,copy(rtmp))
    f(t+c18*Δt,u+Δt*(a1801*k[1]+a1808*k[8]+a1809*k[9]+a1810*k[10]+a1811*k[11]+a1812*k[12]+a1813*k[13]+a1814*k[14]+a1815*k[15]+a1817*k[17]),rtmp); push!(k,copy(rtmp))
    f(t+c19*Δt,u+Δt*(a1901*k[1]+a1908*k[8]+a1909*k[9]+a1910*k[10]+a1911*k[11]+a1912*k[12]+a1913*k[13]+a1914*k[14]+a1915*k[15]+a1917*k[17]+a1918*k[18]),rtmp); push!(k,copy(rtmp))
    f(t+c20*Δt,u+Δt*(a2001*k[1]+a2008*k[8]+a2009*k[9]+a2010*k[10]+a2011*k[11]+a2012*k[12]+a2013*k[13]+a2014*k[14]+a2015*k[15]+a2017*k[17]+a2018*k[18]+a2019*k[19]),rtmp); push!(k,copy(rtmp))
    f(t+c21*Δt,u+Δt*(a2101*k[1]+a2108*k[8]+a2109*k[9]+a2110*k[10]+a2111*k[11]+a2112*k[12]+a2113*k[13]+a2114*k[14]+a2115*k[15]+a2117*k[17]+a2118*k[18]+a2119*k[19]+a2120*k[20]),rtmp); push!(k,copy(rtmp))
    f(t+c22*Δt,u+Δt*(a2201*k[1]+a2208*k[8]+a2209*k[9]+a2210*k[10]+a2211*k[11]+a2212*k[12]+a2213*k[13]+a2214*k[14]+a2215*k[15]+a2217*k[17]+a2218*k[18]+a2219*k[19]+a2220*k[20]+a2221*k[21]),rtmp); push!(k,copy(rtmp))
    f(t+c23*Δt,u+Δt*(a2301*k[1]+a2308*k[8]+a2309*k[9]+a2310*k[10]+a2311*k[11]+a2312*k[12]+a2313*k[13]+a2314*k[14]+a2315*k[15]+a2317*k[17]+a2318*k[18]+a2319*k[19]+a2320*k[20]+a2321*k[21]),rtmp); push!(k,copy(rtmp))
    f(t+c24*Δt,u+Δt*(a2401*k[1]+a2408*k[8]+a2409*k[9]+a2410*k[10]+a2411*k[11]+a2412*k[12]+a2413*k[13]+a2414*k[14]+a2415*k[15]+a2417*k[17]+a2418*k[18]+a2419*k[19]+a2420*k[20]+a2421*k[21]),rtmp); push!(k,copy(rtmp))
    f(t+c25*Δt,u+Δt*(a2501*k[1]+a2508*k[8]+a2509*k[9]+a2510*k[10]+a2511*k[11]+a2512*k[12]+a2513*k[13]+a2514*k[14]+a2515*k[15]+a2517*k[17]+a2518*k[18]+a2519*k[19]+a2520*k[20]+a2521*k[21]),rtmp); push!(k,copy(rtmp))
    f(t+c26*Δt,u+Δt*(a2601*k[1]+a2608*k[8]+a2609*k[9]+a2610*k[10]+a2611*k[11]+a2612*k[12]+a2613*k[13]+a2614*k[14]+a2615*k[15]+a2617*k[17]+a2618*k[18]+a2619*k[19]+a2620*k[20]+a2621*k[21]),rtmp); push!(k,copy(rtmp))
 end
  nothing
end
