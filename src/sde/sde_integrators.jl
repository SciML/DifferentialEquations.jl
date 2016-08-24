immutable SDEIntegrator{T1,uType,uEltype,Nm1,N,tType,tableauType}
  f::Function
  σ::Function
  u::uType
  t::tType
  Δt::tType
  T::tType
  maxiters::Int
  timeseries::GrowableArray{uEltype,uType,N}
  Ws::GrowableArray{uEltype,uType,N}
  ts::Vector{tType}
  timeseries_steps::Int
  save_timeseries::Bool
  adaptive::Bool
  adaptivealg::Symbol
  δ::uEltype
  γ::uEltype
  abstol::uEltype
  reltol::uEltype
  qmax::uEltype
  Δtmax::tType
  Δtmin::tType
  internalnorm::Int
  numvars::Int
  discard_length::tType
  progressbar::Bool
  atomloaded::Bool
  progress_steps::Int
  rands::ChunkedArray{uEltype,Nm1,N}
  sqΔt::tType
  W::uType
  Z::uType
  tableau::tableauType
end

@def sde_preamble begin
  local u::uType
  local t::tType
  local Δt::tType
  local T::tType
  local ΔW::uType
  local ΔZ::uType
  @unpack integrator: f,σ,u,t,Δt,T,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,discard_length,progressbar,atomloaded,progress_steps,rands,sqΔt,W,Z,tableau
  sizeu = size(u)
  iter = 0
  max_stack_size = 0
  max_stack_size2 = 0
  ΔW = sqΔt*next(rands) # Take one first
  ΔZ = sqΔt*next(rands) # Take one first
end

@def sde_sritableaupreamble begin
  local c₀::Vector{uEltype}
  local c₁::Vector{uEltype}
  local A₀::Matrix{uEltype}
  local A₁::Matrix{uEltype}
  local B₀::Matrix{uEltype}
  local B₁::Matrix{uEltype}
  local α::Vector{uEltype}
  local β₁::Vector{uEltype}
  local β₂::Vector{uEltype}
  local β₃::Vector{uEltype}
  local β₄::Vector{uEltype}
end

@def sde_sratableaupreamble begin
  local c₀::Vector{uEltype}
  local c₁::Vector{uEltype}
  local A₀::Matrix{uEltype}
  local B₀::Matrix{uEltype}
  local α::Vector{uEltype}
  local β₁::Vector{uEltype}
  local β₂::Vector{uEltype}
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
    push!(timeseries,u)
    push!(ts,t)
    push!(Ws,W)
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

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:EM,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  @inbounds while t<T
    @sde_loopheader

    u = u + Δt.*f(t,u) + σ(t,u).*ΔW

    t = t + Δt
    W = W + ΔW
    ΔW = sqΔt*next(rands)
    @sde_savevalues
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:AbstractArray,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:EM,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  utmp1 = similar(u); utmp2 = similar(u)
  @inbounds while t<T
    @sde_loopheader
    f(t,u,utmp1)
    σ(t,u,utmp2)
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

function sde_solve{uType<:AbstractArray,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:SRI,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  @sde_sritableaupreamble
  @unpack tableau: c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄
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
  ftemp::uType = similar(u); σtemp::uType = similar(u)
  chi1::uType = similar(u); chi2::uType = similar(u); chi3::uType = similar(u)
  @sde_adaptiveprelim
  @inbounds while t<T
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
        σ(t + c₁[j]*Δt,H1[j],σtemp)
        for k in eachindex(u)
          A0temp[k] += A₀[i,j]*ftemp[k]
          B0temp[k] += B₀[i,j]*σtemp[k]
          A1temp[k] += A₁[i,j]*ftemp[k]
          B1temp[k] += B₁[i,j]*σtemp[k]
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
      σ(t+c₁[i]*Δt,H1[i],σtemp)
      for j in eachindex(u)
        atemp[j] += α[i]*ftemp[j]
        btemp[j] += (β₁[i]*ΔW[j] + β₂[i]*chi1[j])*σtemp[j]
        E₂[j]    += (β₃[i]*chi2[j] + β₄[i]*chi3[j])*σtemp[j]
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

function sde_solve{uType<:AbstractArray,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:SRIW1Optimized,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  chi1::uType = similar(u)
  chi2::uType = similar(u)
  chi3::uType = similar(u)
  fH01o4::uType = similar(u)
  σ₁o2::uType = similar(u)
  H0::uType = similar(u)
  H11::uType = similar(u)
  H12::uType = similar(u)
  H13::uType = similar(u)
  σ₂o3::uType = similar(u)
  Fσ₂o3::uType = similar(u)
  σ₃o3::uType = similar(u)
  Tσ₃o3::uType = similar(u)
  mσ₁::uType = similar(u)
  E₁::uType = similar(u)
  E₂::uType = similar(u)
  fH01::uType = similar(u); fH02::uType = similar(u)
  σ₁::uType = similar(u); σ₂::uType = similar(u); σ₃::uType = similar(u); σ₄::uType = similar(u)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    for i in eachindex(u)
      chi1[i] = (ΔW[i].^2 - Δt)/2sqΔt #I_(1,1)/sqrt(h)
      chi2[i] = (ΔW[i] + ΔZ[i]/sqrt(3))/2 #I_(1,0)/h
      chi3[i] = (ΔW[i].^3 - 3ΔW[i]*Δt)/6Δt #I_(1,1,1)/h
    end
    f(t,u,fH01);fH01*=Δt
    σ(t,u,σ₁)
    Δto4 = Δt/4
    for i in eachindex(u)
      fH01o4[i] = fH01[i]/4
      σ₁o2[i] = σ₁[i]/2
      H0[i] =  u[i] + 3*(fH01o4[i]  + chi2[i]*σ₁o2[i])
      H11[i] = u[i] + fH01o4[i]   + sqΔt*σ₁o2[i]
      H12[i] = u[i] + fH01[i]     - sqΔt*σ₁[i]
    end
    σ(t+Δto4,H11,σ₂)
    σ(t+Δt,H12,σ₃)
    for i in eachindex(u)
      H13[i] = u[i] + fH01o4[i] + sqΔt*(-5σ₁[i] + 3σ₂[i] + σ₃[i]/2)
    end

    σ(t+Δto4,H13,σ₄)
    f(t+3Δto4,H0,fH02); fH02*=Δt
    for i in eachindex(u)
      σ₂o3[i] = σ₂[i]/3
      Fσ₂o3[i] = 4σ₂o3[i]
      σ₃o3[i] = σ₃[i]/3
      Tσ₃o3[i] = 2σ₃o3[i]
      mσ₁[i] = -σ₁[i]
      E₁[i] = fH01[i]+fH02[i]
      E₂[i] = chi2[i]*(2σ₁[i] - Fσ₂o3[i] - Tσ₃o3[i]) + chi3[i]*(2mσ₁[i] + 5σ₂o3[i] - Tσ₃o3[i] + σ₄[i])
      u[i] = u[i] +  (fH01[i] + 2fH02[i])/3 + ΔW[i]*(mσ₁[i] + Fσ₂o3[i] + Tσ₃o3[i]) + chi1[i]*(mσ₁[i] + Fσ₂o3[i] - σ₃o3[i]) + E₂[i]
    end

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:SRIW1Optimized,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  local H0::uType
  @sde_adaptiveprelim
  local fH01::uType; local σ₁::uType
  local fH01o4::uType; local σ₁o2::uType
  local H11::uType; local H12::uType
  local σ₂::uType; local σ₃::uType; local σ₄::uType
  local H13::uType; local fH02::uType
  local σ₂o3::uType; local Fσ₂o3::uType
  local σ₃o3::uType; local Tσ₃o3::uType
  local mσ₁::uType; local E₁::uType; local E₂::uType
  @inbounds while t<T
    @sde_loopheader

    chi1 = (ΔW.^2 - Δt)/2sqΔt #I_(1,1)/sqrt(h)
    chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
    chi3 = (ΔW.^3 - 3ΔW*Δt)/6Δt #I_(1,1,1)/h
    fH01 = Δt*f(t,u)

    σ₁ = σ(t,u)
    fH01o4 = fH01/4
    Δto4 = Δt/4
    σ₁o2 = σ₁/2
    H0 =  u + 3*(fH01o4  + chi2.*σ₁o2)
    H11 = u + fH01o4   + sqΔt*σ₁o2
    H12 = u + fH01     - sqΔt*σ₁
    σ₂ = σ(t+Δto4,H11)
    σ₃ = σ(t+Δt,H12)
    H13 = u + fH01o4 + sqΔt*(-5σ₁ + 3σ₂ + σ₃/2)


    σ₄ = σ(t+Δto4,H13)
    fH02 = Δt*f(t+3Δto4,H0)

    σ₂o3 = σ₂/3
    Fσ₂o3 = 4σ₂o3
    σ₃o3 = σ₃/3
    Tσ₃o3 = 2σ₃o3
    mσ₁ = -σ₁
    E₁ = fH01+fH02
    E₂ = chi2.*(2σ₁ - Fσ₂o3 - Tσ₃o3) + chi3.*(2mσ₁ + 5σ₂o3 - Tσ₃o3 + σ₄)

    u = u + (fH01 + 2fH02)/3 + ΔW.*(mσ₁ + Fσ₂o3 + Tσ₃o3) + chi1.*(mσ₁ + Fσ₂o3 - σ₃o3) + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:SRI,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  @sde_sritableaupreamble
  @unpack tableau: c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄
  stages::Int = length(α)
  H0 = Array{typeof(u)}(stages)
  H1 = Array{typeof(u)}(stages)
  local A0temp::uType; local A1temp::uType
  local B0temp::uType; local B1temp::uType
  local atemp::uType;  local btemp::uType
  local E₁::uType; local E₂::uType
  local E₁temp::uType; local ftemp::uType
  @sde_adaptiveprelim
  @inbounds while t<T
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
        A0temp += A₀[i,j]*f(t + c₀[j]*Δt,H0[..,j])
        B0temp += B₀[i,j]*σ(t + c₁[j]*Δt,H1[..,j])
        A1temp += A₁[i,j]*f(t + c₀[j]*Δt,H0[..,j])
        B1temp += B₁[i,j]*σ(t + c₁[j]*Δt,H1[..,j])
      end
      H0[..,i] = u + A0temp*Δt + B0temp.*chi2
      H1[..,i] = u + A1temp*Δt + B1temp*sqΔt
    end
    atemp = zero(u)
    btemp = zero(u)
    E₂    = zero(u)
    E₁temp= zero(u)
    for i = 1:stages
      ftemp = f(t+c₀[i]*Δt,H0[..,i])
      atemp += α[i]*ftemp
      btemp += (β₁[i]*ΔW + β₂[i]*chi1).*σ(t+c₁[i]*Δt,H1[..,i])
      E₂    += (β₃[i]*chi2 + β₄[i]*chi3).*σ(t+c₁[i]*Δt,H1[..,i])
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

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:SRIVectorized,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  @sde_sritableaupreamble
  @unpack tableau: c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄
  stages::Int = length(α)
  H0 = Array{uEltype}(stages)
  H1 = Array{uEltype}(stages)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    chi1 = .5*(ΔW.^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    chi3 = 1/6 * (ΔW.^3 - 3*ΔW*Δt)/Δt #I_(1,1,1)/h
    H0[:]=zeros(uType,4)
    H1[:]=zeros(uType,4)
    for i = 1:stages
      H0temp = u + Δt*dot(vec(A₀[i,:]),f(t + c₀*Δt,H0)) + chi2*dot(vec(B₀[i,:]),σ(t+c₁*Δt,H1))
      H1[i]  = u + Δt*dot(vec(A₁[i,:]),f(t + c₀*Δt,H0)) + sqΔt*dot(vec(B₁[i,:]),σ(t+c₁*Δt,H1))
      H0[i] = H0temp
    end
    fVec = Δt*f(t+c₀*Δt,H0)
    E₁ = fVec[1]+fVec[2]
    E₂ = dot(β₃*chi2 + β₄*chi3,σ(t+c₁*Δt,H1))

    u = u + dot(α,fVec) + dot(β₁*ΔW + β₂*chi1,σ(t+c₁*Δt,H1)) + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:AbstractArray,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:RKMil,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  du1::uType = similar(u); du2::uType = similar(u)
  K::uType = similar(u); utilde::uType = similar(u); L::uType = similar(u)
  @inbounds while t<T
    @sde_loopheader
    f(t,u,du1)
    L = σ(t,u,du2)
    for i in eachindex(u)
      K[i] = u[i] + Δt*du1[i]
      utilde[i] = K[i] + L[i]*sqΔt
    end
    σ(t,utilde,du1)
    for i in eachindex(u)
      u[i] = K[i]+L[i]*ΔW[i]+(du1[i]-du2[i])./(2sqΔt).*(ΔW[i].^2 - Δt)
      W[i] = W[i] + ΔW[i]
    end
    t = t + Δt
    ΔW = sqΔt*next(rands)
    @sde_savevalues
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:RKMil,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  local L::uType; local K::uType; local utilde::uType
  @inbounds while t<T
    @sde_loopheader

    K = u + Δt.*f(t,u)
    L = σ(t,u)
    utilde = K + L.*sqΔt
    u = K+L.*ΔW+(σ(t,utilde)-σ(t,u))./(2sqΔt).*(ΔW.^2 - Δt)

    t = t + Δt
    W = W + ΔW
    ΔW = sqΔt*next(rands)
    @sde_savevalues
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end


function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:SRA1Optimized,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  H0 = Array{uEltype}(size(u)...,2)
  local k₁::uType; local k₂::uType; local E₁::uType; local E₂::uType
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
    k₁ = Δt*f(t,u)
    k₂ = Δt*f(t+3Δt/4,u+3k₁/4 + 3chi2*σ(t+Δt,u)/2)
    E₁ = k₁ + k₂
    E₂ = chi2.*(σ(t,u)-σ(t+Δt,u)) #Only for additive!

    u = u + k₁/3 + 2k₂/3 + E₂ + ΔW*σ(t+Δt,u)

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:AbstractArray,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:SRA1Optimized,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble

  H0 = Array{uEltype}(size(u)...,2)
  chi2::uType = similar(u)
  tmp1::uType = similar(u)
  E₁::uType = similar(u); σt::uType = similar(u); σpΔt::uType = similar(u)
  E₂::uType = similar(u); k₁::uType = similar(u); k₂::uType = similar(u)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader
    σ(t,u,σt)
    σ(t+Δt,u,σpΔt)
    f(t,u,k₁); k₁*=Δt
    for i in eachindex(u)
      chi2[i] = (ΔW[i] + ΔZ[i]/sqrt(3))/2 #I_(1,0)/h
      tmp1[i] = u[i]+3k₁[i]/4 + 3chi2[i]*σpΔt[i]/2
    end

    f(t+3Δt/4,tmp1,k₂); k₂*=Δt
    for i in eachindex(u)
      E₁[i] = k₁[i] + k₂[i]
      E₂[i] = chi2[i]*(σt[i]-σpΔt[i]) #Only for additive!
      u[i] = u[i] + k₁[i]/3 + 2k₂[i]/3 + E₂[i] + ΔW[i]*σpΔt[i]
    end
    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:AbstractArray,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:SRA,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  @sde_sratableaupreamble
  @unpack tableau: c₀,c₁,A₀,B₀,α,β₁,β₂
  stages::Int = length(α)
  H0 = Vector{typeof(u)}(0)
  for i = 1:stages
    push!(H0,similar(u))
  end
  A0temp::uType = similar(u); B0temp::uType = similar(u)
  ftmp::uType = similar(u); σtmp::uType = similar(u); chi2::uType = similar(u)
  atemp::uType = similar(u); btemp::uType = similar(u); E₂::uType = similar(u); E₁temp::uType = similar(u)
  E₁::uType = similar(u)
  @sde_adaptiveprelim
  @inbounds while t<T
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
        σ(t + c₁[j]*Δt,H0[j],σtmp)
        for k in eachindex(u)
          A0temp[k] += A₀[i,j]*ftmp[k]
          B0temp[k] += B₀[i,j]*σtmp[k]
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
      f(t+c₀[i]*Δt,H0[..,i],ftmp)
      σ(t+c₁[i]*Δt,H0[..,i],σtmp)
      for j in eachindex(u)
        atemp[j] += α[i]*ftmp[j]
        btemp[j] += (β₁[i]*ΔW[j])*σtmp[j]
        E₂[j]    += (β₂[i]*chi2[j])*σtmp[j]
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

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:SRA,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  @sde_sratableaupreamble
  @unpack tableau: c₀,c₁,A₀,B₀,α,β₁,β₂
  stages::Int = length(α)
  H0 = Array{uEltype}(stages)
  local atemp::uType; local btemp::uType
  local E₂::uType; local E₁::uType; local E₁temp::uType
  local ftemp::uType; local A0temp::uType; local B0temp::uType
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    H0[:]=zeros(stages)
    for i = 1:stages
      A0temp = zero(u)
      B0temp = zero(u)
      for j = 1:i-1
        A0temp += A₀[i,j]*f(t + c₀[j]*Δt,H0[..,j])
        B0temp += B₀[i,j]*σ(t + c₁[j]*Δt,H0[..,j]) #H0[..,i] argument ignored
      end
      H0[..,i] = u + A0temp*Δt + B0temp.*chi2
    end

    atemp = zero(u)
    btemp = zero(u)
    E₂    = zero(u)
    E₁temp= zero(u)

    for i = 1:stages
      ftemp = f(t+c₀[i]*Δt,H0[..,i])
      atemp += α[i]*ftemp
      btemp += (β₁[i]*ΔW ).*σ(t+c₁[i]*Δt,H0[..,i]) #H0[..,i] argument ignored
      E₂    += (β₂[i]*chi2).*σ(t+c₁[i]*Δt,H0[..,i]) #H0[..,i] argument ignored

      E₁temp += ftemp
    end
    E₁ = Δt*E₁temp
    u = u + Δt*atemp + btemp + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:Number,uEltype<:Number,Nm1,N,tType<:Number,tableauType<:Tableau}(integrator::SDEIntegrator{:SRAVectorized,uType,uEltype,Nm1,N,tType,tableauType})
  @sde_preamble
  @sde_sratableaupreamble
  @unpack tableau: c₀,c₁,A₀,B₀,α,β₁,β₂
  stages::Int = length(α)
  @sde_adaptiveprelim
  H0 = Array{uEltype}(stages)

  @inbounds while t<T
    @sde_loopheader

    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    H0[:]=zeros(stages)
    for i = 1:stages
      H0[i] = u + Δt*dot(vec(A₀[i,:]),f(t + c₀*Δt,H0)) + chi2*dot(vec(B₀[i,:]),σ(t+c₁*Δt,H0))
    end
    fVec = f(t+c₀*Δt,H0)
    E₁ = Δt*(fVec[1]+fVec[2])
    E₂ = dot(β₂*chi2,σ(t+c₁*Δt,H0))

    u = u + Δt*dot(α,f(t+c₀*Δt,H0)) + dot(β₁*ΔW,σ(t+c₁*Δt,H0)) + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end
