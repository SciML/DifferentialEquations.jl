immutable SDEIntegrator{T1,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType}
  f::Function
  g::Function
  u::uType
  t::tType
  Δt::tType
  T::tType
  maxiters::Int
  timeseries::Vector{uType}
  Ws::Vector{randType}
  ts::Vector{tType}
  timeseries_steps::Int
  save_timeseries::Bool
  adaptive::Bool
  adaptivealg::Symbol
  δ::uEltypeNoUnits
  γ::uEltypeNoUnits
  abstol::uEltype
  reltol::uEltypeNoUnits
  qmax::uEltypeNoUnits
  Δtmax::tType
  Δtmin::tType
  internalnorm::Int
  numvars::Int
  discard_length::tType
  progressbar::Bool
  atomloaded::Bool
  progress_steps::Int
  rands::ChunkedArray{uEltypeNoUnits,Nm1,N}
  sqΔt::tType
  W::randType
  Z::randType
  tableau::tableauType
end

@def sde_preamble begin
  local u::uType
  local t::tType
  local Δt::tType
  local T::tType
  local ΔW::randType
  local ΔZ::randType
  @unpack f,g,u,t,Δt,T,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,discard_length,progressbar,atomloaded,progress_steps,rands,sqΔt,W,Z,tableau = integrator
  sizeu = size(u)
  iter = 0
  max_stack_size = 0
  max_stack_size2 = 0
  ΔW = sqΔt*next(rands) # Take one first
  ΔZ = sqΔt*next(rands) # Take one first
end

@def sde_sritableaupreamble begin
  local c₀::Vector{uEltypeNoUnits}
  local c₁::Vector{uEltypeNoUnits}
  local A₀::Matrix{uEltypeNoUnits}
  local A₁::Matrix{uEltypeNoUnits}
  local B₀::Matrix{uEltypeNoUnits}
  local B₁::Matrix{uEltypeNoUnits}
  local α::Vector{uEltypeNoUnits}
  local β₁::Vector{uEltypeNoUnits}
  local β₂::Vector{uEltypeNoUnits}
  local β₃::Vector{uEltypeNoUnits}
  local β₄::Vector{uEltypeNoUnits}
end

@def sde_sratableaupreamble begin
  local c₀::Vector{uEltypeNoUnits}
  local c₁::Vector{uEltypeNoUnits}
  local A₀::Matrix{uEltypeNoUnits}
  local B₀::Matrix{uEltypeNoUnits}
  local α::Vector{uEltypeNoUnits}
  local β₁::Vector{uEltypeNoUnits}
  local β₂::Vector{uEltypeNoUnits}
end

@def sde_loopheader begin
  iter += 1
  if iter > maxiters
    warn("Max Iters Reached. Aborting")
    # u = map((x)->oftype(x,NaN),u)
    break
  end
  if any(isnan,u)
    warn("NaNs detected. Aborting")
    break
  end
end

@def sde_savevalues begin
  if save_timeseries && iter%timeseries_steps==0
    push!(timeseries,copy(u))
    push!(ts,t)
    push!(Ws,copy(W))
  end
end

@def sde_loopfooter begin
  t = t + Δt
  W = W + ΔW
  Z = Z + ΔZ
  ΔW = sqΔt*next(rands)
  ΔZ = sqΔt*next(rands)
  @sde_savevalues
  (atomloaded && progressbar && iter%progress_steps==0) ? Main.Atom.progress(t/T) : nothing #Use Atom's progressbar if loaded
end

@def sde_adaptiveprelim begin
max_stack_size = 0
max_stack_size2 = 0
end

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:Number,rateType<:Number}(integrator::SDEIntegrator{:EM,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  @fastmath @inbounds while t<T
    @sde_loopheader

    u = u + Δt.*f(t,u) + g(t,u).*ΔW

    t = t + Δt
    W = W + ΔW
    ΔW = sqΔt*next(rands)
    @sde_savevalues
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:AbstractArray,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:AbstractArray,rateType<:AbstractArray}(integrator::SDEIntegrator{:EM,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  utmp1 = similar(u); utmp2 = similar(u)
  @fastmath @inbounds while t<T
    @sde_loopheader
    f(t,u,utmp1)
    g(t,u,utmp2)
    for i in eachindex(u)
      u[i] = u[i] + Δt*utmp1[i] + utmp2[i]*ΔW[i]
      W[i] = W[i] + ΔW[i]
    end
    t = t + Δt
    ΔW = sqΔt*next(rands)
    @sde_savevalues
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:AbstractArray,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:AbstractArray,rateType<:AbstractArray}(integrator::SDEIntegrator{:SRI,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  @sde_sritableaupreamble
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄ = tableau
  stages = length(α)
  H0 = Vector{typeof(u)}(0)
  H1 = Vector{typeof(u)}(0)
  for i = 1:stages
    push!(H0,similar(u))
    push!(H1,similar(u))
  end
  #TODO Reduce memory
  A0temp::uType = similar(u); A1temp::uType = similar(u)
  B0temp::uType = similar(u); B1temp::uType = similar(u)
  A0temp2::uType = similar(u); A1temp2::uType = similar(u)
  B0temp2::uType = similar(u); B1temp2::uType = similar(u)
  atemp::uType = similar(u); btemp::uType = similar(u)
  E₁::uType = similar(u); E₂::uType = similar(u); E₁temp::uType = similar(u)
  ftemp::uType = similar(u); gtemp::uType = similar(u)
  chi1::randType = similar(ΔW); chi2::randType = similar(ΔW); chi3::randType = similar(ΔW)
  @sde_adaptiveprelim
  @fastmath @inbounds while t<T
    @sde_loopheader

    for i in eachindex(u)
      chi1[i] = .5*(ΔW[i].^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
      chi2[i] = .5*(ΔW[i] + ΔZ[i]/sqrt(3)) #I_(1,0)/h
      chi3[i] = 1/6 * (ΔW[i].^3 - 3*ΔW[i]*Δt)/Δt #I_(1,1,1)/h
    end
    for i=1:stages
      H0[i][:]=zero(uEltype)
      H1[i][:]=zero(uEltype)
    end
    for i = 1:stages
      A0temp[:]=zero(uEltype)
      B0temp[:]=zero(uEltype)
      A1temp[:]=zero(uEltype)
      B1temp[:]=zero(uEltype)
      for j = 1:i-1
        f(t + c₀[j]*Δt,H0[j],ftemp)
        g(t + c₁[j]*Δt,H1[j],gtemp)
        for k in eachindex(u)
          A0temp[k] += A₀[i,j]*ftemp[k]
          B0temp[k] += B₀[i,j]*gtemp[k]
          A1temp[k] += A₁[i,j]*ftemp[k]
          B1temp[k] += B₁[i,j]*gtemp[k]
        end
      end
      H0[i] = u + A0temp*Δt + B0temp.*chi2
      H1[i] = u + A1temp*Δt + B1temp*sqΔt
    end
    atemp[:]=zero(uEltype)
    btemp[:]=zero(uEltype)
    E₂[:]=zero(uEltype)
    E₁temp[:]=zero(uEltype)
    for i = 1:stages
      f(t+c₀[i]*Δt,H0[i],ftemp)
      g(t+c₁[i]*Δt,H1[i],gtemp)
      for j in eachindex(u)
        atemp[j] += α[i]*ftemp[j]
        btemp[j] += (β₁[i]*ΔW[j] + β₂[i]*chi1[j])*gtemp[j]
        E₂[j]    += (β₃[i]*chi2[j] + β₄[i]*chi3[j])*gtemp[j]
      end
      if i<3 #1 or 2
        for j in eachindex(u)
          E₁temp[j] += ftemp[j]
        end
      end
    end

    for i in eachindex(u)
      E₁[i] = Δt*E₁temp[i]
      u[i] = u[i] + Δt*atemp[i] + btemp[i] + E₂[i]
    end

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:AbstractArray,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:AbstractArray,rateType<:AbstractArray}(integrator::SDEIntegrator{:SRIW1Optimized,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  chi1::randType = similar(ΔW)
  chi2::randType = similar(ΔW)
  chi3::randType = similar(ΔW)
  fH01o4::uType = similar(u)
  g₁o2::uType = similar(u)
  H0::uType = similar(u)
  H11::uType = similar(u)
  H12::uType = similar(u)
  H13::uType = similar(u)
  g₂o3::uType = similar(u)
  Fg₂o3::uType = similar(u)
  g₃o3::uType = similar(u)
  Tg₃o3::uType = similar(u)
  mg₁::uType = similar(u)
  E₁::uType = similar(u)
  E₂::uType = similar(u)
  fH01::uType = similar(u); fH02::uType = similar(u)
  g₁::uType = similar(u); g₂::uType = similar(u); g₃::uType = similar(u); g₄::uType = similar(u)
  @sde_adaptiveprelim
  @fastmath @inbounds while t<T
    @sde_loopheader

    for i in eachindex(u)
      chi1[i] = (ΔW[i].^2 - Δt)/2sqΔt #I_(1,1)/sqrt(h)
      chi2[i] = (ΔW[i] + ΔZ[i]/sqrt(3))/2 #I_(1,0)/h
      chi3[i] = (ΔW[i].^3 - 3ΔW[i]*Δt)/6Δt #I_(1,1,1)/h
    end
    f(t,u,fH01);fH01*=Δt
    g(t,u,g₁)
    Δto4 = Δt/4
    for i in eachindex(u)
      fH01o4[i] = fH01[i]/4
      g₁o2[i] = g₁[i]/2
      H0[i] =  u[i] + 3*(fH01o4[i]  + chi2[i]*g₁o2[i])
      H11[i] = u[i] + fH01o4[i]   + sqΔt*g₁o2[i]
      H12[i] = u[i] + fH01[i]     - sqΔt*g₁[i]
    end
    g(t+Δto4,H11,g₂)
    g(t+Δt,H12,g₃)
    for i in eachindex(u)
      H13[i] = u[i] + fH01o4[i] + sqΔt*(-5g₁[i] + 3g₂[i] + g₃[i]/2)
    end

    g(t+Δto4,H13,g₄)
    f(t+3Δto4,H0,fH02); fH02*=Δt
    for i in eachindex(u)
      g₂o3[i] = g₂[i]/3
      Fg₂o3[i] = 4g₂o3[i]
      g₃o3[i] = g₃[i]/3
      Tg₃o3[i] = 2g₃o3[i]
      mg₁[i] = -g₁[i]
      E₁[i] = fH01[i]+fH02[i]
      E₂[i] = chi2[i]*(2g₁[i] - Fg₂o3[i] - Tg₃o3[i]) + chi3[i]*(2mg₁[i] + 5g₂o3[i] - Tg₃o3[i] + g₄[i])
      u[i] = u[i] +  (fH01[i] + 2fH02[i])/3 + ΔW[i]*(mg₁[i] + Fg₂o3[i] + Tg₃o3[i]) + chi1[i]*(mg₁[i] + Fg₂o3[i] - g₃o3[i]) + E₂[i]
    end

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:Number,rateType<:Number}(integrator::SDEIntegrator{:SRIW1Optimized,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  local H0::uType
  @sde_adaptiveprelim
  local fH01::uType; local g₁::uType
  local fH01o4::uType; local g₁o2::uType
  local H11::uType; local H12::uType
  local g₂::uType; local g₃::uType; local g₄::uType
  local H13::uType; local fH02::uType
  local g₂o3::uType; local Fg₂o3::uType
  local g₃o3::uType; local Tg₃o3::uType
  local mg₁::uType; local E₁::uType; local E₂::uType
  @fastmath @inbounds while t<T
    @sde_loopheader

    chi1 = (ΔW.^2 - Δt)/2sqΔt #I_(1,1)/sqrt(h)
    chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
    chi3 = (ΔW.^3 - 3ΔW*Δt)/6Δt #I_(1,1,1)/h
    fH01 = Δt*f(t,u)

    g₁ = g(t,u)
    fH01o4 = fH01/4
    Δto4 = Δt/4
    g₁o2 = g₁/2
    H0 =  u + 3*(fH01o4  + chi2.*g₁o2)
    H11 = u + fH01o4   + sqΔt*g₁o2
    H12 = u + fH01     - sqΔt*g₁
    g₂ = g(t+Δto4,H11)
    g₃ = g(t+Δt,H12)
    H13 = u + fH01o4 + sqΔt*(-5g₁ + 3g₂ + g₃/2)


    g₄ = g(t+Δto4,H13)
    fH02 = Δt*f(t+3Δto4,H0)

    g₂o3 = g₂/3
    Fg₂o3 = 4g₂o3
    g₃o3 = g₃/3
    Tg₃o3 = 2g₃o3
    mg₁ = -g₁
    E₁ = fH01+fH02
    E₂ = chi2.*(2g₁ - Fg₂o3 - Tg₃o3) + chi3.*(2mg₁ + 5g₂o3 - Tg₃o3 + g₄)

    u = u + (fH01 + 2fH02)/3 + ΔW.*(mg₁ + Fg₂o3 + Tg₃o3) + chi1.*(mg₁ + Fg₂o3 - g₃o3) + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:Number,rateType<:Number}(integrator::SDEIntegrator{:SRI,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  @sde_sritableaupreamble
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄ = tableau
  stages::Int = length(α)
  H0 = Array{typeof(u)}(stages)
  H1 = Array{typeof(u)}(stages)
  local A0temp::uType; local A1temp::uType
  local B0temp::uType; local B1temp::uType
  local atemp::uType;  local btemp::uType
  local E₁::uType; local E₂::uType
  local E₁temp::uType; local ftemp::uType
  @sde_adaptiveprelim
  @fastmath @inbounds while t<T
    @sde_loopheader

    chi1 = .5*(ΔW.^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    chi3 = 1/6 * (ΔW.^3 - 3*ΔW*Δt)/Δt #I_(1,1,1)/h

    H0[:]=zero(typeof(u))
    H1[:]=zero(typeof(u))
    for i = 1:stages
      A0temp = zero(u)
      B0temp = zero(u)
      A1temp = zero(u)
      B1temp = zero(u)
      for j = 1:i-1
        A0temp += A₀[i,j]*f(t + c₀[j]*Δt,H0[j])
        B0temp += B₀[i,j]*g(t + c₁[j]*Δt,H1[j])
        A1temp += A₁[i,j]*f(t + c₀[j]*Δt,H0[j])
        B1temp += B₁[i,j]*g(t + c₁[j]*Δt,H1[j])
      end
      H0[i] = u + A0temp*Δt + B0temp.*chi2
      H1[i] = u + A1temp*Δt + B1temp*sqΔt
    end
    atemp = zero(u)
    btemp = zero(u)
    E₂    = zero(u)
    E₁temp= zero(u)
    for i = 1:stages
      ftemp = f(t+c₀[i]*Δt,H0[i])
      atemp += α[i]*ftemp
      btemp += (β₁[i]*ΔW + β₂[i]*chi1).*g(t+c₁[i]*Δt,H1[i])
      E₂    += (β₃[i]*chi2 + β₄[i]*chi3).*g(t+c₁[i]*Δt,H1[i])
      if i<3 #1 or 2
        E₁temp += ftemp
      end
    end
    E₁ = Δt*E₁temp

    u = u + Δt*atemp + btemp + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:Number,rateType<:Number}(integrator::SDEIntegrator{:SRIVectorized,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  @sde_sritableaupreamble
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄ = tableau
  stages::Int = length(α)
  H0 = Array{uEltype}(stages)
  H1 = Array{uEltype}(stages)
  @sde_adaptiveprelim
  @fastmath @inbounds while t<T
    @sde_loopheader

    chi1 = .5*(ΔW.^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    chi3 = 1/6 * (ΔW.^3 - 3*ΔW*Δt)/Δt #I_(1,1,1)/h
    H0[:]=zeros(uType,4)
    H1[:]=zeros(uType,4)
    for i = 1:stages
      H0temp = u + Δt*dot(vec(A₀[i,:]),f(t + c₀*Δt,H0)) + chi2*dot(vec(B₀[i,:]),g(t+c₁*Δt,H1))
      H1[i]  = u + Δt*dot(vec(A₁[i,:]),f(t + c₀*Δt,H0)) + sqΔt*dot(vec(B₁[i,:]),g(t+c₁*Δt,H1))
      H0[i] = H0temp
    end
    fVec = Δt*f(t+c₀*Δt,H0)
    E₁ = fVec[1]+fVec[2]
    E₂ = dot(β₃*chi2 + β₄*chi3,g(t+c₁*Δt,H1))

    u = u + dot(α,fVec) + dot(β₁*ΔW + β₂*chi1,g(t+c₁*Δt,H1)) + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:AbstractArray,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:AbstractArray,rateType<:AbstractArray}(integrator::SDEIntegrator{:RKMil,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  du1::uType = similar(u); du2::uType = similar(u)
  K::uType = similar(u); utilde::uType = similar(u); L::uType = similar(u)
  @fastmath @inbounds while t<T
    @sde_loopheader
    f(t,u,du1)
    g(t,u,L)
    for i in eachindex(u)
      K[i] = u[i] + Δt*du1[i]
      utilde[i] = K[i] + L[i]*sqΔt
    end
    g(t,utilde,du2)
    for i in eachindex(u)
      u[i] = K[i]+L[i]*ΔW[i]+(du2[i]-L[i])./(2sqΔt).*(ΔW[i].^2 - Δt)
      W[i] = W[i] + ΔW[i]
    end
    t = t + Δt
    ΔW = sqΔt*next(rands)
    @sde_savevalues
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:Number,rateType<:Number}(integrator::SDEIntegrator{:RKMil,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  local L::uType; local K::uType; local utilde::uType
  @fastmath @inbounds while t<T
    @sde_loopheader

    K = u + Δt.*f(t,u)
    L = g(t,u)
    utilde = K + L*sqΔt
    u = K+L*ΔW+(g(t,utilde)-g(t,u))/(2sqΔt)*(ΔW^2 - Δt)

    t = t + Δt
    W = W + ΔW
    ΔW = sqΔt*next(rands)
    @sde_savevalues
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end


function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:Number,rateType<:Number}(integrator::SDEIntegrator{:SRA1Optimized,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  H0 = Array{uEltype}(size(u)...,2)
  local k₁::uType; local k₂::uType; local E₁::uType; local E₂::uType
  @sde_adaptiveprelim
  @fastmath @inbounds while t<T
    @sde_loopheader

    chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
    k₁ = Δt*f(t,u)
    k₂ = Δt*f(t+3Δt/4,u+3k₁/4 + 3chi2*g(t+Δt,u)/2)
    E₁ = k₁ + k₂
    E₂ = chi2.*(g(t,u)-g(t+Δt,u)) #Only for additive!

    u = u + k₁/3 + 2k₂/3 + E₂ + ΔW*g(t+Δt,u)

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:AbstractArray,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:AbstractArray,rateType<:AbstractArray}(integrator::SDEIntegrator{:SRA1Optimized,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble

  H0 = Array{uEltype}(size(u)...,2)
  chi2::randType = similar(ΔW)
  tmp1::uType = similar(u)
  E₁::uType = similar(u); gt::uType = similar(u); gpΔt::uType = similar(u)
  E₂::uType = similar(u); k₁::uType = similar(u); k₂::uType = similar(u)
  @sde_adaptiveprelim
  @fastmath @inbounds while t<T
    @sde_loopheader
    g(t,u,gt)
    g(t+Δt,u,gpΔt)
    f(t,u,k₁); k₁*=Δt
    for i in eachindex(u)
      chi2[i] = (ΔW[i] + ΔZ[i]/sqrt(3))/2 #I_(1,0)/h
      tmp1[i] = u[i]+3k₁[i]/4 + 3chi2[i]*gpΔt[i]/2
    end

    f(t+3Δt/4,tmp1,k₂); k₂*=Δt
    for i in eachindex(u)
      E₁[i] = k₁[i] + k₂[i]
      E₂[i] = chi2[i]*(gt[i]-gpΔt[i]) #Only for additive!
      u[i] = u[i] + k₁[i]/3 + 2k₂[i]/3 + E₂[i] + ΔW[i]*gpΔt[i]
    end
    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:AbstractArray,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:AbstractArray,rateType<:AbstractArray}(integrator::SDEIntegrator{:SRA,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  @sde_sratableaupreamble
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂ = tableau
  stages::Int = length(α)
  H0 = Vector{typeof(u)}(0)
  for i = 1:stages
    push!(H0,similar(u))
  end
  A0temp::uType = similar(u); B0temp::uType = similar(u)
  ftmp::uType = similar(u); gtmp::uType = similar(u); chi2::uType = similar(u)
  atemp::uType = similar(u); btemp::uType = similar(u); E₂::uType = similar(u); E₁temp::uType = similar(u)
  E₁::uType = similar(u)
  @sde_adaptiveprelim
  @fastmath @inbounds while t<T
    @sde_loopheader
    for i in eachindex(u)
      chi2[i] = .5*(ΔW[i] + ΔZ[i]/sqrt(3)) #I_(1,0)/h
    end
    for i in 1:stages
      H0[i][:]=zero(uEltype)
    end
    for i = 1:stages
      A0temp[:] = zero(uEltype)
      B0temp[:] = zero(uEltype)
      for j = 1:i-1
        f(t + c₀[j]*Δt,H0[j],ftmp)
        g(t + c₁[j]*Δt,H0[j],gtmp)
        for k in eachindex(u)
          A0temp[k] += A₀[i,j]*ftmp[k]
          B0temp[k] += B₀[i,j]*gtmp[k]
        end
      end
      for j in eachindex(u)
        H0[i][j] = u[j] + A0temp[j]*Δt + B0temp[j]*chi2[j]
      end
    end
    atemp[:] = zero(uEltype)
    btemp[:] = zero(uEltype)
    E₂[:]    = zero(uEltype)
    E₁temp[:]= zero(uEltype)

    for i = 1:stages
      f(t+c₀[i]*Δt,H0[i],ftmp)
      g(t+c₁[i]*Δt,H0[i],gtmp)
      for j in eachindex(u)
        atemp[j] += α[i]*ftmp[j]
        btemp[j] += (β₁[i]*ΔW[j])*gtmp[j]
        E₂[j]    += (β₂[i]*chi2[j])*gtmp[j]
        E₁temp[j] += ftmp[j]
      end
    end
    for i in eachindex(u)
      E₁[i] = Δt*E₁temp[i]
      u[i] = u[i] + Δt*atemp[i] + btemp[i] + E₂[i]
    end

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:Number,rateType<:Number}(integrator::SDEIntegrator{:SRA,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  @sde_sratableaupreamble
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂ = tableau
  stages::Int = length(α)
  H0 = Array{uEltype}(stages)
  local atemp::uType; local btemp::uType
  local E₂::uType; local E₁::uType; local E₁temp::uType
  local ftemp::uType; local A0temp::uType; local B0temp::uType
  @sde_adaptiveprelim
  @fastmath @inbounds while t<T
    @sde_loopheader

    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    H0[:]=zeros(stages)
    for i = 1:stages
      A0temp = zero(u)
      B0temp = zero(u)
      for j = 1:i-1
        A0temp += A₀[i,j]*f(t + c₀[j]*Δt,H0[j])
        B0temp += B₀[i,j]*g(t + c₁[j]*Δt,H0[j]) #H0[..,i] argument ignored
      end
      H0[i] = u + A0temp*Δt + B0temp.*chi2
    end

    atemp = zero(u)
    btemp = zero(u)
    E₂    = zero(u)
    E₁temp= zero(u)

    for i = 1:stages
      ftemp = f(t+c₀[i]*Δt,H0[i])
      atemp += α[i]*ftemp
      btemp += (β₁[i]*ΔW ).*g(t+c₁[i]*Δt,H0[i]) #H0[..,i] argument ignored
      E₂    += (β₂[i]*chi2).*g(t+c₁[i]*Δt,H0[i]) #H0[..,i] argument ignored

      E₁temp += ftemp
    end
    E₁ = Δt*E₁temp
    u = u + Δt*atemp + btemp + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau,uEltypeNoUnits<:Number,randType<:Number,rateType<:Number}(integrator::SDEIntegrator{:SRAVectorized,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType})
  @sde_preamble
  @sde_sratableaupreamble
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂ = tableau
  stages::Int = length(α)
  @sde_adaptiveprelim
  H0 = Array{uEltype}(stages)

  @fastmath @inbounds while t<T
    @sde_loopheader

    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    H0[:]=zeros(stages)
    for i = 1:stages
      H0[i] = u + Δt*dot(vec(A₀[i,:]),f(t + c₀*Δt,H0)) + chi2*dot(vec(B₀[i,:]),g(t+c₁*Δt,H0))
    end
    fVec = f(t+c₀*Δt,H0)
    E₁ = Δt*(fVec[1]+fVec[2])
    E₂ = dot(β₂*chi2,g(t+c₁*Δt,H0))

    u = u + Δt*dot(α,f(t+c₀*Δt,H0)) + dot(β₁*ΔW,g(t+c₁*Δt,H0)) + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end
