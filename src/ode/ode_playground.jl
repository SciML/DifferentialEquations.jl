
function ode_solve{uType<:AbstractArray,uEltype<:Float64,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:DP5Threaded,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = constructDP5(uEltypeNoUnits)
  k2::rateType = rateType(sizeu)
  k3::rateType = rateType(sizeu)
  k4::rateType = rateType(sizeu)
  k5::rateType = rateType(sizeu)
  k6::rateType = rateType(sizeu)
  update::rateType = rateType(sizeu)
  bspl::rateType = rateType(sizeu)
  utilde = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  uidx::Base.OneTo{Int64} = eachindex(u)
  const kshortsize = 4
  if calck
    d1,d3,d4,d5,d6,d7 = DP5_dense_ds(uEltypeNoUnits)
    k = ksEltype()
    for i in 1:kshortsize
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    # Setup k pointers
    k[1] = update
    if calcprevs
      kprev = deepcopy(k)
    end
  end
  k1 = fsalfirst; k7 = fsallast # Setup pointers
  f(t,u,fsalfirst);  # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      dp5threaded_loop1(Δt,tmp,u,k1,uidx)
      f(t+c1*Δt,tmp,k2)
      dp5threaded_loop2(Δt,tmp,u,k1,k2,uidx)
      f(t+c2*Δt,tmp,k3)
      dp5threaded_loop3(Δt,tmp,u,k1,k2,k3,uidx)
      f(t+c3*Δt,tmp,k4)
      dp5threaded_loop4(Δt,tmp,u,k1,k2,k3,k4,uidx)
      f(t+c4*Δt,tmp,k5)
      dp5threaded_loop5(Δt,tmp,u,k1,k2,k3,k4,k5,uidx)
      f(t+Δt,tmp,k6)
      dp5threaded_loop6(Δt,utmp,u,k1,k3,k4,k5,k6,update,uidx)
      f(t+Δt,utmp,k7)
      if adaptive
        dp5threaded_adaptiveloop(Δt,utilde,u,k1,k3,k4,k5,k6,k7,atmp,utmp,abstol,reltol,uidx)
        EEst = sqrt( sum(atmp) * normfactor)
      else
        recursivecopy!(u, utmp)
      end
      if calck
        dp5threaded_denseloop(bspl,update,k1,k3,k4,k5,k6,k7,k,uidx)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end


@noinline function dp5threaded_loop1{T<:Float64,N}(Δt,tmp::Array{T,N},u,k1,uidx)
  Threads.@threads for i in uidx
    tmp[i] = u[i]+Δt*(0.2*k1[i])
  end
end

@noinline function dp5threaded_loop2{T<:Float64,N}(Δt,tmp::Array{T,N},u,k1,k2,uidx)
  Threads.@threads for i in uidx
    tmp[i] = u[i]+Δt*(0.075*k1[i]+0.225*k2[i])
  end
end

@noinline function dp5threaded_loop3{T<:Float64,N}(Δt,tmp::Array{T,N},u,k1,k2,k3,uidx)
  Threads.@threads for i in uidx
    tmp[i] = u[i]+Δt*(0.9777777777777777*k1[i]-3.7333333333333334*k2[i]+3.5555555555555554*k3[i])
  end
end

@noinline function dp5threaded_loop4{T<:Float64,N}(Δt,tmp::Array{T,N},u,k1,k2,k3,k4,uidx)
  Threads.@threads for i in uidx
    tmp[i] =u[i]+ Δt*(2.9525986892242035*k1[i]+-11.595793324188385*k2[i]+9.822892851699436*k3[i]-0.2908093278463649*k4[i])
  end
end

@noinline function dp5threaded_loop5{T<:Float64,N}(Δt,tmp::Array{T,N},u,k1,k2,k3,k4,k5,uidx)
  Threads.@threads for i in uidx
    tmp[i] = u[i]+Δt*(2.8462752525252526*k1[i]-10.757575757575758*k2[i]+8.906422717743473*k3[i]+0.2784090909090909*k4[i]-0.2735313036020583*k5[i])
  end
end

@noinline function dp5threaded_loop6{T<:Float64,N}(Δt,utmp::Array{T,N},u,k1,k3,k4,k5,k6,update,uidx)
  Threads.@threads for i in uidx
    update[i] = 0.09114583333333333*k1[i]+0.44923629829290207*k3[i]+0.6510416666666666*k4[i]-0.322376179245283*k5[i]+0.13095238095238096*k6[i]
    utmp[i] = u[i]+Δt*update[i]
  end
end

@noinline function dp5threaded_adaptiveloop{T<:Float64,N}(Δt,utilde::Array{T,N},u,k1,k3,k4,k5,k6,k7,tmp,utmp,abstol,reltol,uidx)
  Threads.@threads for i in uidx
    utilde[i] = u[i] + Δt*(0.08991319444444444*k1[i] + 0.4534890685834082*k3[i] + 0.6140625*k4[i] - 0.2715123820754717*k5[i] + 0.08904761904761904*k6[i] + 0.025*k7[i])
    tmp[i] = ((utilde[i]-utmp[i])/(abstol+max(abs(u[i]),abs(utmp[i]))*reltol))^2
  end
end

@noinline function dp5threaded_denseloop{T<:Float64,N}(bspl::Array{T,N},update,k1,k3,k4,k5,k6,k7,k,uidx)
  Threads.@threads for i in uidx
    bspl[i] = k1[i] - update[i]
    k[2][i] = bspl[i]
    k[3][i] = update[i] - k7[i] - bspl[i]
    k[4][i] = (-1.1270175653862835k1[i]+2.675424484351598k3[i]+-5.685526961588504k4[i]+3.5219323679207912k5[i]+-1.7672812570757455k6[i]+2.382468931778144k7[i])
  end
end
