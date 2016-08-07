immutable SDEIntegrator{T1,uType,tType,tableauType}
  f::Function
  σ::Function
  u::uType
  t::tType
  Δt::tType
  T::tType
  maxiters::Int
  timeseries::AbstractArray
  Ws::GrowableArray
  ts::AbstractArray
  timeseries_steps::Int
  save_timeseries::Bool
  adaptive::Bool
  adaptivealg::Symbol
  δ::uType
  γ::tType
  abstol::uType
  reltol::uType
  qmax::tType
  Δtmax::tType
  Δtmin::tType
  internalnorm::Int
  numvars::Int
  discard_length::tType
  progressbar::Bool
  atomloaded::Bool
  progress_steps::Int
  rands::ChunkedArray
  sqΔt::tType
  W::uType
  Z::uType
  tableau::tableauType
end

@def sde_preamble begin
  @unpack integrator: f,σ,u,t,Δt,T,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,atomloaded,progress_steps,rands,sqΔt,W,Z,tableau
  sizeu = size(u)
  iter = 0
  max_stack_size = 0
  max_stack_size2 = 0
  ΔW = sqΔt*next(rands) # Take one first
  ΔZ = sqΔt*next(rands) # Take one first
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

function sde_solve(integrator::SDEIntegrator{:EM,:Number})
  @sde_preamble
  @inbounds while t<T
    @sde_loopheader

    u = u + Δt.*f(u,t) + σ(u,t).*ΔW

    t = t + Δt
    W = W + ΔW
    ΔW = sqΔt*next(rands)
    @sde_savevalues
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve(integrator::SDEIntegrator{:EM,:AbstractArray})
  @sde_preamble
  utmp1 = similar(u); utmp2 = similar(u)
  @inbounds while t<T
    @sde_loopheader
    f(utmp1,u,t)
    σ(utmp2,u,t)
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

function sde_solve(integrator::SDEIntegrator{:SRI,:AbstractArray})
  @sde_preamble
  @unpack tableau: c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄
  stages = length(α)
  H0 = Vector{typeof(u)}(0)
  H1 = Vector{typeof(u)}(0)
  for i = 1:stages
    push!(H0,similar(u))
    push!(H1,similar(u))
  end
  #TODO Reduce memory
  A0temp = similar(u); A1temp = similar(u)
  B0temp = similar(u); B1temp = similar(u)
  A0temp2 = similar(u); A1temp2 = similar(u)
  B0temp2 = similar(u); B1temp2 = similar(u)
  atemp = similar(u); btemp = similar(u)
  E₁ = similar(u); E₂ = similar(u); E₁temp = similar(u)
  ftemp = similar(u); σtemp = similar(u)
  chi1 = similar(u); chi2 = similar(u); chi3 = similar(u)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    for i in eachindex(u)
      chi1[i] = .5*(ΔW[i].^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
      chi2[i] = .5*(ΔW[i] + ΔZ[i]/sqrt(3)) #I_(1,0)/h
      chi3[i] = 1/6 * (ΔW[i].^3 - 3*ΔW[i]*Δt)/Δt #I_(1,1,1)/h
    end
    for i=1:stages
      H0[i][:]=zero(eltype(u))
      H1[i][:]=zero(eltype(u))
    end
    for i = 1:stages
      A0temp[:]=zero(eltype(u))
      B0temp[:]=zero(eltype(u))
      A1temp[:]=zero(eltype(u))
      B1temp[:]=zero(eltype(u))
      for j = 1:i-1
        f(ftemp,H0[j],t + c₀[j]*Δt)
        σ(σtemp,H1[j],t + c₁[j]*Δt)
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
    atemp[:]=zero(eltype(u))
    btemp[:]=zero(eltype(u))
    E₂[:]=zero(eltype(u))
    E₁temp[:]=zero(eltype(u))
    for i = 1:stages
      f(ftemp,H0[i],t+c₀[i]*Δt)
      σ(σtemp,H1[i],t+c₁[i]*Δt)
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

function sde_solve(integrator::SDEIntegrator{:SRIW1Optimized,:AbstractArray})
  @sde_preamble
  H0 = Array{eltype(u)}(size(u)...,4)
  H1 = Array{eltype(u)}(size(u)...,4)
  chi1 = similar(u)
  chi2 = similar(u)
  chi3 = similar(u)
  fH01o4 = similar(u)
  σ₁o2 = similar(u)
  H0 = similar(u)
  H11 = similar(u)
  H12 = similar(u)
  H13 = similar(u)
  σ₂o3 = similar(u)
  Fσ₂o3 = similar(u)
  σ₃o3 = similar(u)
  Tσ₃o3 = similar(u)
  mσ₁ = similar(u)
  E₁ = similar(u)
  E₂ = similar(u)
  fH01 = similar(u); fH02 = similar(u)
  σ₁ = similar(u); σ₂ = similar(u); σ₃ = similar(u); σ₄ = similar(u)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    for i in eachindex(u)
      chi1[i] = (ΔW[i].^2 - Δt)/2sqΔt #I_(1,1)/sqrt(h)
      chi2[i] = (ΔW[i] + ΔZ[i]/sqrt(3))/2 #I_(1,0)/h
      chi3[i] = (ΔW[i].^3 - 3ΔW[i]*Δt)/6Δt #I_(1,1,1)/h
    end
    f(fH01,u,t);fH01*=Δt
    σ(σ₁,u,t)
    Δto4 = Δt/4
    for i in eachindex(u)
      fH01o4[i] = fH01[i]/4
      σ₁o2[i] = σ₁[i]/2
      H0[i] =  u[i] + 3*(fH01o4[i]  + chi2[i]*σ₁o2[i])
      H11[i] = u[i] + fH01o4[i]   + sqΔt*σ₁o2[i]
      H12[i] = u[i] + fH01[i]     - sqΔt*σ₁[i]
    end
    σ(σ₂,H11,t+Δto4)
    σ(σ₃,H12,t+Δt)
    for i in eachindex(u)
      H13[i] = u[i] + fH01o4[i] + sqΔt*(-5σ₁[i] + 3σ₂[i] + σ₃[i]/2)
    end

    σ(σ₄,H13,t+Δto4)
    f(fH02,H0,t+3Δto4); fH02*=Δt
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

function sde_solve(integrator::SDEIntegrator{:SRIW1Optimized,:Number})
  @sde_preamble
  H0 = Array{eltype(u)}(size(u)...,4)
  H1 = Array{eltype(u)}(size(u)...,4)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    chi1 = (ΔW.^2 - Δt)/2sqΔt #I_(1,1)/sqrt(h)
    chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
    chi3 = (ΔW.^3 - 3ΔW*Δt)/6Δt #I_(1,1,1)/h
    fH01 = Δt*f(u,t)

    σ₁ = σ(u,t)
    fH01o4 = fH01/4
    Δto4 = Δt/4
    σ₁o2 = σ₁/2
    H0 =  u + 3*(fH01o4  + chi2.*σ₁o2)
    H11 = u + fH01o4   + sqΔt*σ₁o2
    H12 = u + fH01     - sqΔt*σ₁
    σ₂ = σ(H11,t+Δto4)
    σ₃ = σ(H12,t+Δt)
    H13 = u + fH01o4 + sqΔt*(-5σ₁ + 3σ₂ + σ₃/2)


    σ₄ = σ(H13,t+Δto4)
    fH02 = Δt*f(H0,t+3Δto4)

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

function sde_solve(integrator::SDEIntegrator{:SRI,:Number})
  @sde_preamble
  @unpack tableau: c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄
  stages = length(α)
  H0 = Array{typeof(u)}(stages)
  H1 = Array{typeof(u)}(stages)
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
        A0temp += A₀[i,j]*f(H0[..,j],t + c₀[j]*Δt)
        B0temp += B₀[i,j]*σ(H1[..,j],t + c₁[j]*Δt)
        A1temp += A₁[i,j]*f(H0[..,j],t + c₀[j]*Δt)
        B1temp += B₁[i,j]*σ(H1[..,j],t + c₁[j]*Δt)
      end
      H0[..,i] = u + A0temp*Δt + B0temp.*chi2
      H1[..,i] = u + A1temp*Δt + B1temp*sqΔt
    end
    atemp = zero(u)
    btemp = zero(u)
    E₂    = zero(u)
    E₁temp= zero(u)
    for i = 1:stages
      ftemp = f(H0[..,i],t+c₀[i]*Δt)
      atemp += α[i]*ftemp
      btemp += (β₁[i]*ΔW + β₂[i]*chi1).*σ(H1[..,i],t+c₁[i]*Δt)
      E₂    += (β₃[i]*chi2 + β₄[i]*chi3).*σ(H1[..,i],t+c₁[i]*Δt)
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

function sde_solve(integrator::SDEIntegrator{:SRIVectorized,:Number})
  @sde_preamble
  uType = typeof(u)
  @unpack tableau: c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄
  stages = length(α)
  H0 = Array{eltype(u)}(stages)
  H1 = Array{eltype(u)}(stages)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    chi1 = .5*(ΔW.^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    chi3 = 1/6 * (ΔW.^3 - 3*ΔW*Δt)/Δt #I_(1,1,1)/h
    H0[:]=zeros(uType,4)
    H1[:]=zeros(uType,4)
    for i = 1:stages
      H0temp = u + Δt*dot(vec(A₀[i,:]),f(H0,t + c₀*Δt)) + chi2*dot(vec(B₀[i,:]),σ(H1,t+c₁*Δt))
      H1[i]  = u + Δt*dot(vec(A₁[i,:]),f(H0,t + c₀*Δt)) + sqΔt*dot(vec(B₁[i,:]),σ(H1,t+c₁*Δt))
      H0[i] = H0temp
    end
    fVec = Δt*f(H0,t+c₀*Δt)
    E₁ = fVec[1]+fVec[2]
    E₂ = dot(β₃*chi2 + β₄*chi3,σ(H1,t+c₁*Δt))

    u = u + dot(α,fVec) + dot(β₁*ΔW + β₂*chi1,σ(H1,t+c₁*Δt)) + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve(integrator::SDEIntegrator{:RKMil,:AbstractArray})
  @sde_preamble
  du1 = similar(u); du2 = similar(u)
  K = similar(u); utilde = similar(u); L = similar(u)
  @inbounds while t<T
    @sde_loopheader
    f(du1,u,t)
    L = σ(du2,u,t)
    for i in eachindex(u)
      K[i] = u[i] + Δt*du1[i]
      utilde[i] = K[i] + L[i]*sqΔt
    end
    σ(du1,utilde,t)
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

function sde_solve(integrator::SDEIntegrator{:RKMil,:Number})
  @sde_preamble
  @inbounds while t<T
    @sde_loopheader

    K = u + Δt.*f(u,t)
    L = σ(u,t)
    utilde = K + L.*sqΔt
    u = K+L.*ΔW+(σ(utilde,t)-σ(u,t))./(2sqΔt).*(ΔW.^2 - Δt)

    t = t + Δt
    W = W + ΔW
    ΔW = sqΔt*next(rands)
    @sde_savevalues
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end


function sde_solve(integrator::SDEIntegrator{:SRA1Optimized,:Number})
  @sde_preamble
  uType = typeof(u)
  H0 = Array{eltype(u)}(size(u)...,2)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
    k₁ = Δt*f(u,t)
    k₂ = Δt*f(u+3k₁/4 + 3chi2*σ(u,t+Δt)/2,t+3Δt/4)
    E₁ = k₁ + k₂
    E₂ = chi2.*(σ(u,t)-σ(u,t+Δt)) #Only for additive!

    u = u + k₁/3 + 2k₂/3 + E₂ + ΔW*σ(u,t+Δt)

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve(integrator::SDEIntegrator{:SRA1Optimized,:AbstractArray})
  @sde_preamble
  uType = typeof(u)
  H0 = Array{eltype(u)}(size(u)...,2)
  chi2 = similar(u)
  tmp1 = similar(u)
  E₁ = similar(u); σt = similar(u); σpΔt = similar(u)
  E₂ = similar(u); k₁ = similar(u); k₂ = similar(u)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader
    σ(σt,u,t)
    σ(σpΔt,u,t+Δt)
    f(k₁,u,t); k₁*=Δt
    for i in eachindex(u)
      chi2[i] = (ΔW[i] + ΔZ[i]/sqrt(3))/2 #I_(1,0)/h
      tmp1[i] = u[i]+3k₁[i]/4 + 3chi2[i]*σpΔt[i]/2
    end

    f(k₂,tmp1,t+3Δt/4); k₂*=Δt
    for i in eachindex(u)
      E₁[i] = k₁[i] + k₂[i]
      E₂[i] = chi2[i]*(σt[i]-σpΔt[i]) #Only for additive!
      u[i] = u[i] + k₁[i]/3 + 2k₂[i]/3 + E₂[i] + ΔW[i]*σpΔt[i]
    end
    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve(integrator::SDEIntegrator{:SRA,:AbstractArray})
  @sde_preamble
  uType = typeof(u)
  @unpack tableau: c₀,c₁,A₀,B₀,α,β₁,β₂
  stages = length(α)
  H0 = Vector{typeof(u)}(0)
  for i = 1:stages
    push!(H0,similar(u))
  end
  A0temp = similar(u); B0temp = similar(u)
  ftmp = similar(u); σtmp = similar(u); chi2 = similar(u)
  atemp = similar(u); btemp = similar(u); E₂ = similar(u); E₁temp = similar(u)
  E₁ = similar(u)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader
    for i in eachindex(u)
      chi2[i] = .5*(ΔW[i] + ΔZ[i]/sqrt(3)) #I_(1,0)/h
    end
    for i in 1:stages
      H0[i][:]=zero(eltype(u))
    end
    for i = 1:stages
      A0temp[:] = zero(eltype(u))
      B0temp[:] = zero(eltype(u))
      for j = 1:i-1
        f(ftmp,H0[j],t + c₀[j]*Δt)
        σ(σtmp,H0[j],t + c₁[j]*Δt)
        for k in eachindex(u)
          A0temp[k] += A₀[i,j]*ftmp[k]
          B0temp[k] += B₀[i,j]*σtmp[k]
        end
      end
      for j in eachindex(u)
        H0[i][j] = u[j] + A0temp[j]*Δt + B0temp[j]*chi2[j]
      end
    end
    atemp[:] = zero(eltype(u))
    btemp[:] = zero(eltype(u))
    E₂[:]    = zero(eltype(u))
    E₁temp[:]= zero(eltype(u))

    for i = 1:stages
      f(ftmp,H0[..,i],t+c₀[i]*Δt)
      σ(σtmp,H0[..,i],t+c₁[i]*Δt)
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

function sde_solve(integrator::SDEIntegrator{:SRA,:Number})
  @sde_preamble
  uType = typeof(u)
  @unpack tableau: c₀,c₁,A₀,B₀,α,β₁,β₂
  stages = length(α)
  H0 = Array{eltype(u)}(stages)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    H0[:]=zeros(stages)
    for i = 1:stages
      A0temp = zero(u)
      B0temp = zero(u)
      for j = 1:i-1
        A0temp += A₀[i,j]*f(H0[..,j],t + c₀[j]*Δt)
        B0temp += B₀[i,j]*σ(H0[..,j],t + c₁[j]*Δt) #H0[..,i] argument ignored
      end
      H0[..,i] = u + A0temp*Δt + B0temp.*chi2
    end

    atemp = zero(u)
    btemp = zero(u)
    E₂    = zero(u)
    E₁temp= zero(u)

    for i = 1:stages
      ftemp = f(H0[..,i],t+c₀[i]*Δt)
      atemp += α[i]*ftemp
      btemp += (β₁[i]*ΔW ).*σ(H0[..,i],t+c₁[i]*Δt) #H0[..,i] argument ignored
      E₂    += (β₂[i]*chi2).*σ(H0[..,i],t+c₁[i]*Δt) #H0[..,i] argument ignored
      E₁temp += ftemp
    end
    E₁ = Δt*E₁temp
    u = u + Δt*atemp + btemp + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve(integrator::SDEIntegrator{:SRAVectorized,:Number})
  @sde_preamble
  uType = typeof(u)
  @unpack tableau: c₀,c₁,A₀,B₀,α,β₁,β₂
  stages = length(α)
  @sde_adaptiveprelim
  H0 = Array{eltype(u)}(stages)
  @inbounds while t<T
    @sde_loopheader

    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    H0[:]=zeros(stages)
    for i = 1:stages
      H0[i] = u + Δt*dot(vec(A₀[i,:]),f(H0,t + c₀*Δt)) + chi2*dot(vec(B₀[i,:]),σ(H0,t+c₁*Δt))
    end
    fVec = f(H0,t+c₀*Δt)
    E₁ = Δt*(fVec[1]+fVec[2])
    E₂ = dot(β₂*chi2,σ(H0,t+c₁*Δt))

    u = u + Δt*dot(α,f(H0,t+c₀*Δt)) + dot(β₁*ΔW,σ(H0,t+c₁*Δt)) + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end
