immutable ODEIntegrator{Alg,uType<:Union{AbstractArray,Number},uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Union{AbstractArray,Number},ksEltype<:Union{AbstractArray,Number}} <: DEIntegrator
  f::Function
  u::uType
  t::tType
  Δt::tType
  Ts::Vector{tType}
  maxiters::Int
  timeseries_steps::Int
  save_timeseries::Bool
  adaptive::Bool
  abstol::uEltype
  reltol::uEltypeNoUnits
  γ::uEltypeNoUnits
  qmax::uEltypeNoUnits
  qmin::uEltypeNoUnits
  Δtmax::tType
  Δtmin::tType
  internalnorm::Int
  progressbar::Bool
  tableau::ExplicitRKTableau
  autodiff::Bool
  adaptiveorder::Int
  order::Int
  atomloaded::Bool
  progress_steps::Int
  β::uEltypeNoUnits
  timechoicealg::Symbol
  qoldinit::uEltypeNoUnits
  normfactor::uEltypeNoUnits
  fsal::Bool
  dense::Bool
  saveat::Vector{tType}
end

@def ode_preamble begin
  local u::uType
  local t::tType
  local Δt::tType
  local Ts::Vector{tType}
  local adaptiveorder::Int
  @unpack integrator: f,u,t,Δt,Ts,maxiters,timeseries_steps,γ,qmax,qmin,save_timeseries,adaptive,progressbar,autodiff,adaptiveorder,order,atomloaded,progress_steps,β,timechoicealg,qoldinit,normfactor,fsal, dense, saveat

  issimple_dense = (ksEltype==rateType) # Means ks[i] = f(t[i],timeseries[i]), for Hermite

  # Need to initiate ks in the method

  Tfinal = Ts[end]
  local iter::Int = 0
  sizeu = size(u)
  local utmp::uType
  local k::ksEltype

  # Setup FSAL
  if uType <: Number
    utmp = zero(uType)
    fsallast = zero(rateType)
    if fsal
      fsalfirst = f(t,u)
    else
      fsalfirst = zero(rateType)
    end
  else
    utmp = zeros(u)
    fsallast = rateType(sizeu)
    fsalfirst = rateType(sizeu)
    if fsal
      f(t,u,fsalfirst)
    end
  end

  local standard::uEltype = 0
  local q::uEltypeNoUnits = 0
  local Δtpropose::tType = 0
  local q11::uEltypeNoUnits = 0
  #local k1::uType; local k7::uType
  local qold::uEltypeNoUnits = 0

  expo1 = 1/order - 0.75β
  qminc = inv(qmin)
  qmaxc = inv(qmax)
  local Eest::uEltypeNoUnits = zero(uEltypeNoUnits)
  if adaptive
    @unpack integrator: abstol,reltol,qmax,Δtmax,Δtmin,internalnorm
  end

  timeseries = Vector{uType}(0)

  push!(timeseries,copy(u))
  ts = Vector{tType}(0)
  push!(ts,t)
  ks = Vector{ksEltype}(0)
  if dense && issimple_dense #If issimple_dense, then ks[1]=f(ts[1],timeseries[1])
    if ksEltype <: AbstractArray
      k = rateType(sizeu)
    end
    if fsal
      push!(ks,copy(fsalfirst)) # already computed
    elseif ksEltype <: Number
      push!(ks,f(t,u))
    elseif ksEltype <: AbstractArray
      f(t,u,k)
      push!(ks,deepcopy(k))
    end
  end ## if not simple_dense, you have to initialize k and push the ks[1]!

  (progressbar && atomloaded && iter%progress_steps==0) ? Main.Atom.progress(0) : nothing #Use Atom's progressbar if loaded
end

@def ode_loopheader begin
  iter += 1
  if iter > maxiters
    warn("Max Iters Reached. Aborting")
    # u = map((x)->oftype(x,NaN),u)
    return u,t,timeseries,ts,ks
  end
  Δt = min(Δt,abs(T-t))
end

@def ode_savevalues begin
  if save_timeseries && iter%timeseries_steps==0
    push!(timeseries,copy(u))
    push!(ts,t)
    if dense
      push!(ks,deepcopy(k))
    end
  end
end

@def ode_implicitsavevalues begin
  if save_timeseries && iter%timeseries_steps==0
    push!(timeseries,copy(u))
    push!(ts,t)
    if dense
      push!(ks,deepcopy(k))
    end
  end
end

@def ode_numberimplicitsavevalues begin
  if save_timeseries && iter%timeseries_steps==0
    push!(timeseries,uhold[1])
    push!(ts,t)
    if dense
      push!(ks,deepcopy(k))
    end
  end
end

@def ode_loopfooter begin
  if adaptive
    if timechoicealg == :Lund #Lund stabilization of q
      q11 = EEst^expo1
      q = q11/(qold^β)
      q = max(qmaxc,min(qminc,q/γ))
      Δtnew = Δt/q
      if EEst < 1.0 # Accept
        t = t + Δt
        qold = max(EEst,qoldinit)
        copy!(u, utmp)
        @ode_savevalues
        Δtpropose = min(Δtmax,Δtnew)
        Δt = max(Δtpropose,Δtmin) #abs to fix complex sqrt issue at end
        if fsal
          copy!(fsalfirst,fsallast)
        end
      else # Reject
        Δt = Δt/min(qminc,q11/γ)
      end
    elseif timechoicealg == :Simple
      standard = γ*abs(1/(EEst))^(1/(adaptiveorder))
      if isinf(standard)
          q = qmax
      else
         q = min(qmax,max(standard,eps()))
      end
      if q > 1 # Accept
        t = t + Δt
        copy!(u, utmp)
        @ode_savevalues
        if fsal
          copy!(fsalfirst,fsallast)
        end
      end
      Δtpropose = min(Δtmax,q*Δt)
      Δt = max(min(Δtpropose,abs(T-t)),Δtmin) #abs to fix complex sqrt issue at end
    end
  else #Not adaptive
    t += Δt
    @ode_savevalues
    if fsal
      copy!(fsalfirst,fsallast)
    end
  end
  (progressbar && atomloaded && iter%progress_steps==0) ? Main.Atom.progress(t/Tfinal) : nothing #Use Atom's progressbar if loaded
end

@def ode_numberloopfooter begin
  if adaptive
    if timechoicealg == :Lund #Lund stabilization of q
      q11 = EEst^expo1
      q = q11/(qold^β)
      q = max(qmaxc,min(qminc,q/γ))
      Δtnew = Δt/q
      if EEst < 1.0 # Accept
        t = t + Δt
        qold = max(EEst,qoldinit)
        u = utmp
        @ode_savevalues
        Δtpropose = min(Δtmax,Δtnew)
        Δt = max(Δtpropose,Δtmin) #abs to fix complex sqrt issue at end
        if fsal
          fsalfirst = fsallast
        end
      else # Reject
        Δt = Δt/min(qminc,q11/γ)
      end
    elseif timechoicealg == :Simple
      standard = γ*abs(1/(EEst))^(1/(adaptiveorder))
      if isinf(standard)
          q = qmax
      else
         q = min(qmax,max(standard,eps()))
      end
      if q > 1 # Accept
        t = t + Δt
        u = utmp
        @ode_savevalues
        if fsal
          fsalfirst = fsallast
        end
      end
      Δtpropose = min(Δtmax,q*Δt)
      Δt = max(min(Δtpropose,abs(T-t)),Δtmin) #abs to fix complex sqrt issue at end
    end
  else #Not adaptive
    t += Δt
    @ode_savevalues
    if fsal
      fsalfirst = fsallast
    end
  end
  (progressbar && atomloaded && iter%progress_steps==0) ? Main.Atom.progress(t/Tfinal) : nothing #Use Atom's progressbar if loaded
end

@def ode_implicitloopfooter begin
  if adaptive
    if timechoicealg == :Lund #Lund stabilization of q
      q11 = EEst^expo1
      q = q11/(qold^β)
      q = max(qmaxc,min(qminc,q/γ))
      Δtnew = Δt/q
      if EEst < 1.0 # Accept
        t = t + Δt
        qold = max(EEst,qoldinit)
        copy!(uhold, utmp)
        @ode_implicitsavevalues
        Δtpropose = min(Δtmax,Δtnew)
        Δt = max(Δtpropose,Δtmin) #abs to fix complex sqrt issue at end
      else # Reject
        Δt = Δt/min(qminc,q11/γ)
      end
    elseif timechoicealg == :Simple
      standard = γ*abs(1/(EEst))^(1/(adaptiveorder))
      if isinf(standard)
          q = qmax
      else
         q = min(qmax,max(standard,eps()))
      end
      if q > 1 # Accept
        t = t + Δt
        copy!(uhold, utmp)
        @ode_implicitsavevalues
      end
      Δtpropose = min(Δtmax,q*Δt)
      Δt = max(min(Δtpropose,abs(T-t)),Δtmin) #abs to fix complex sqrt issue at end
    end  else #Not adaptive
    t = t + Δt
    @ode_implicitsavevalues
  end
  (progressbar && atomloaded && iter%progress_steps==0) ? Main.Atom.progress(t/Tfinal) : nothing #Use Atom's progressbar if loaded
end

@def ode_numberimplicitloopfooter begin
  if adaptive
    if timechoicealg == :Lund #Lund stabilization of q
      q11 = EEst^expo1
      q = q11/(qold^β)
      q = max(qmaxc,min(qminc,q/γ))
      Δtnew = Δt/q
      if EEst < 1.0 # Accept
        qold = max(EEst,qoldinit)
        t = t + Δt
        uhold = utmp
        @ode_numberimplicitsavevalues
        Δtpropose = min(Δtmax,Δtnew)
        Δt = max(Δtpropose,Δtmin) #abs to fix complex sqrt issue at end
      else # Reject
        Δt = Δt/min(qminc,q11/γ)
      end
    elseif timechoicealg == :Simple
      standard = γ*abs(1/(EEst))^(1/(adaptiveorder))
      if isinf(standard)
          q = qmax
      else
         q = min(qmax,max(standard,eps()))
      end
      if q > 1
        t = t + Δt
        uhold = utmp
        @ode_numberimplicitsavevalues
      end
      Δtpropose = min(Δtmax,q*Δt)
      Δt = max(min(Δtpropose,abs(T-t)),Δtmin) #abs to fix complex sqrt issue at end
    end
  else #Not adaptive
    t = t + Δt
    @ode_numberimplicitsavevalues
  end
  (progressbar && atomloaded && iter%progress_steps==0) ? Main.Atom.progress(t/Tfinal) : nothing #Use Atom's progressbar if loaded
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:Euler,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  k = f(t,u) # For the interpolation, needs k at the updated point
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      u = u + Δt*k
      k = f(t,u) # For the interpolation, needs k at the updated point
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Euler,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  uidx = eachindex(u)
  if !dense
    k = rateType(sizeu) # Not initialized if not dense
  end
  f(t,u,k) # For the interpolation, needs k at the updated point
  @inbounds for T in Ts
      while t < T
      @ode_loopheader
      for i in uidx
        u[i] = u[i] + Δt*k[i]
      end
      f(t,u,k) # For the interpolation, needs k at the updated point
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:Midpoint,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  halfΔt::tType = Δt/2
  local du::rateType
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k = f(t+halfΔt,u+halfΔt*f(t,u))
      u = u + Δt*k
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Midpoint,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  halfΔt::tType = Δt/2
  utilde::uType = similar(u)
  uidx = eachindex(u)
  if !dense
    k = rateType(sizeu) # Not initialized if not dense
  end
  du = rateType(sizeu)
  @inbounds for T in Ts
      while t < T
      @ode_loopheader
      f(t,u,du)
      for i in uidx
        utilde[i] = u[i]+halfΔt*du[i]
      end
      f(t+halfΔt,utilde,k)
      for i in uidx
        u[i] = u[i] + Δt*k[i]
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:RK4,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  halfΔt::tType = Δt/2
  local k₁::rateType
  local k₂::rateType
  local k₃::rateType
  local k₄::rateType
  local ttmp::tType
  @inbounds for T in Ts
      while t < T
      @ode_loopheader
      k₁ = f(t,u)
      ttmp = t+halfΔt
      k₂ = f(ttmp,u+halfΔt*k₁)
      k₃ = f(ttmp,u+halfΔt*k₂)
      k₄ = f(t+Δt,u+Δt*k₃)
      u = u + Δt*(k₁ + 2(k₂ + k₃) + k₄)/6
      if dense
        k=k₄
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:RK4,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  halfΔt::tType = Δt/2
  k₁ = rateType(sizeu)
  k₂ = rateType(sizeu)
  k₃ = rateType(sizeu)
  k₄ = rateType(sizeu)
  tmp = similar(u)
  uidx = eachindex(u)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k₁)
      ttmp = t+halfΔt
      for i in uidx
        tmp[i] = u[i]+halfΔt*k₁[i]
      end
      f(ttmp,tmp,k₂)
      for i in uidx
        tmp[i] = u[i]+halfΔt*k₂[i]
      end
      f(ttmp,tmp,k₃)
      for i in uidx
        tmp[i] = u[i]+Δt*k₃[i]
      end
      f(t+Δt,tmp,k₄)
      for i in uidx
        u[i] = u[i] + Δt*(k₁[i] + 2k₂[i] + 2k₃[i] + k₄[i])/6
      end
      if dense
        k=k₄
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:ExplicitRK,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  local A::Matrix{uEltypeNoUnits}
  local c::Vector{uEltypeNoUnits}
  local α::Vector{uEltypeNoUnits}
  local αEEst::Vector{uEltypeNoUnits}
  local stages::Int
  local rtmp::ratetype
  @unpack integrator.tableau: A,c,α,αEEst,stages
  A = A' # Transpose A to column major looping
  kk = Array{typeof(u)}(stages) # Not ks since that's for dense
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      # Calc First
      if fsal
        kk[1] = Δt*fsalfirst
      else
        kk[1] = Δt*f(t,u)
      end
      # Calc Middle
      for i = 2:stages-1
        utilde = zero(u)
        for j = 1:i-1
          utilde += A[j,i]*kk[j]
        end
        kk[i] = Δt*f(t+c[i]*Δt,u+utilde);
      end
      #Calc Last
      utilde = zero(u)
      for j = 1:stages-1
        utilde += A[j,end]*kk[j]
      end
      fsallast = f(t+c[end]*Δt,u+utilde); kk[end]=Δt*fsallast # Uses fsallast as temp even if not fsal
      # Accumulate Result
      utilde = α[1]*kk[1]
      for i = 2:stages
        utilde += α[i]*kk[i]
      end
      if adaptive
        utmp = u + utilde
        uEEst = αEEst[1]*kk[1]
        for i = 2:stages
          uEEst += αEEst[i]*kk[i]
        end
        EEst = abs( (utilde-uEEst)/(abstol+max(u,utmp)*reltol))
      else
        u = u + utilde
      end
      if dense
        k = kk[end]
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:ExplicitRK,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  local A::Matrix{uEltypeNoUnits}
  local c::Vector{uEltypeNoUnits}
  local α::Vector{uEltypeNoUnits}
  local αEEst::Vector{uEltypeNoUnits}
  local stages::Int
  uidx = eachindex(u)
  @unpack integrator.tableau: A,c,α,αEEst,stages
  A = A' # Transpose A to column major looping
  kk = Vector{typeof(u)}(0)
  for i = 1:stages
    push!(kk,similar(u))
  end
  utilde = similar(u)
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  rtmp = rateType(sizeu)
  utmp = zeros(u)
  uEEst = similar(u)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      # First
      if fsal
        for l in uidx
          kk[1][l] = Δt*fsalfirst[l]
        end
      else
        f(t,u,rtmp)
        for l in uidx
          kk[1][l]=Δt*rtmp[l]
        end
      end
      # Middle
      for i = 2:stages-1
        utilde[:] = zero(eltype(u))
        for j = 1:i-1
          for l in uidx
            utilde[l] += A[j,i]*kk[j][l]
          end
        end
        for l in uidx
          tmp[l] = u[l]+utilde[l]
        end
        f(t+c[i]*Δt,tmp,rtmp)
        for l in uidx
          kk[i][l]=Δt*rtmp[l]
        end
      end
      #Last
      utilde[:] = zero(eltype(u))
      for j = 1:stages-1
        for l in uidx
          utilde[l] += A[j,end]*kk[j][l]
        end
      end
      for l in uidx
        tmp[l] = u[l]+utilde[l]
      end
      f(t+c[end]*Δt,tmp,fsallast); #fsallast is tmp even if not fsal
      for l in uidx
        kk[end][l]= Δt*fsallast[l]
      end
      #Accumulate
      utilde[:] = α[1]*kk[1]
      for i = 2:stages
        for l in uidx
          utilde[l] += α[i]*kk[i][l]
        end
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + utilde[i]
        end
        uEEst[:] = αEEst[1]*kk[1]
        for i = 2:stages
          for j in uidx
            uEEst[j] += αEEst[i]*kk[i][j]
          end
        end
        for i in uidx
          atmp[i] = ((utilde[i]-uEEst[i])/(abstol+max(u[i],utmp[i])*reltol))^2
        end
        EEst = sqrt( sum(atmp) * normfactor)
      else
        for i in uidx
          u[i] = u[i] + utilde[i]
        end
      end
      if dense
        k = kk[end]
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:ExplicitRKVectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  local A::Matrix{uEltypeNoUnits}
  local c::Vector{uEltypeNoUnits}
  local α::Vector{uEltypeNoUnits}
  local αEEst::Vector{uEltypeNoUnits}
  local stages::Int
  rtmp = rateType(sizeu)
  uidx = eachindex(u)
  @unpack integrator.tableau: A,c,α,αEEst,stages
  ks = Vector{typeof(u)}(0)
  for i = 1:stages
    push!(ks,similar(u))
  end
  A = A' # Transpose A to column major looping
  utilde = similar(u)
  tmp = similar(u)
  utmp = zeros(u)
  uEEst = similar(u)
  if fsal
    f(t,u,fsalfirst)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      #First
      if fsal
        ks[1] = Δt*fsalfirst
      else
        f(t,u,rtmp); ks[1]=Δt*rtmp
      end
      #Middle
      for i = 2:stages-1
        utilde[:] = zero(eltype(u))
        for j = 1:i-1
          utilde += A[j,i]*ks[j]
        end
        tmp = u+utilde
        f(t+c[i]*Δt,tmp,rtmp); ks[i]=Δt*rtmp
      end
      # Last
      utilde[:] = zero(eltype(u))
      for j = 1:stages-1
        utilde += A[j,end]*ks[j]
      end
      tmp = u+utilde
      f(t+c[end]*Δt,tmp,fsallast); ks[end]=fsallast*Δt
      #Accumulate
      utilde[:] = α[1]*ks[1]
      for i = 2:stages
        utilde += α[i]*ks[i]
      end
      if adaptive
        utmp = u + utilde
        uEEst[:] = αEEst[1]*ks[1]
        for i = 2:stages
          uEEst += αEEst[i]*ks[i]
        end
        EEst = sqrt( sum(((utilde-uEEst)./(abstol+max(u,utmp)*reltol)).^2) * normfactor)
      else
        u = u + utilde
      end
      if dense
        k = kk[end]
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:BS3,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  a21,a32,a41,a42,a43,c1,c2,b1,b2,b3,b4  = constructBS3(uEltypeNoUnits)
  local k1::uType
  local k2::uType
  local k3::uType
  local k4::uType
  local utilde::uType
  fsalfirst = f(t,u) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*fsalfirst
      k2 = Δt*f(t+c1*Δt,u+a21*k1)
      k3 = Δt*f(t+c2*Δt,u+a32*k2)
      utmp = u+a41*k1+a42*k2+a43*k3
      fsallast = f(t+Δt,utmp); k4 = Δt*fsallast
      if adaptive
        utilde = u + b1*k1 + b2*k2 + b3*k3 + b4*k4
        EEst = sqrt( sum(((utilde-utmp)/(abstol+max(u,utmp)*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if dense
        k = fsallast
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:BS3Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  a21,a32,a41,a42,a43,c1,c2,b1,b2,b3,b4  = constructBS3(uEltypeNoUnits)
  k1 = similar(u)
  k2 = similar(u)
  k3 = similar(u)
  k4 = similar(u)
  local utilde::uType
  local rtmp::rateType = rateType(sizeu)
  f(t,u,fsalfirst) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*fsalfirst
      f(t+c1*Δt,u+a21*k1,rtmp); k2=Δt*rtmp
      f(t+c2*Δt,u+a32*k2,rtmp); k3=Δt*rtmp
      utmp = u+a41*k1+a42*k2+a43*k3
      f(t+Δt,utmp,fsallast); k4 = Δt*fsallast
      if adaptive
        utilde = u + b1*k1 + b2*k2 + b3*k3 + b4*k4
        EEst = sqrt( sum(((utilde-utmp)./(abstol+max(u,utmp)*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if dense
        k = fsallast
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:BS3,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  a21,a32,a41,a42,a43,c1,c2,b1,b2,b3,b4  = constructBS3(uEltypeNoUnits)
  k1 = similar(u)
  k2 = similar(u)
  k3 = similar(u)
  k4 = similar(u)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  local rtmp = rateType(sizeu)
  uidx = eachindex(u)
  tmp = similar(u)
  k = fsallast
  f(t,u,fsalfirst) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      for i in uidx
        k1[i] = Δt*fsalfirst[i]
        tmp[i] = u[i]+a21*k1[i]
      end
      f(t+c1*Δt,tmp,rtmp)
      for i in uidx
        k2[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a32*k2[i]
      end
      f(t+c2*Δt,tmp,rtmp)
      for i in uidx
        k3[i] = Δt*rtmp[i]
        utmp[i] = u[i]+a41*k1[i]+a42*k2[i]+a43*k3[i]
      end
      f(t+Δt,utmp,fsallast);
      for i in uidx
        k4[i] = Δt*fsallast[i]
      end
      if adaptive
        for i in uidx
          utilde[i] = u[i] + b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i]
          atmp[i] = ((utilde[i]-utmp[i])/(abstol+max(u[i],utmp[i])*reltol))^2
        end
        EEst = sqrt( sum(atmp) * normfactor)
      else
        copy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:BS5,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat2,bhat3,bhat4,bhat5,bhat6,bhat7,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = constructBS5(uEltypeNoUnits)
  local k1::uType
  local k2::uType
  local k3::uType
  local k4::uType
  local k5::uType
  local k6::uType
  local k7::uType
  local k8::uType
  local utilde::uType
  local EEst2::uEltypeNoUnits
  if dense
    k = ksEltype()
    for i in 1:8
      push!(k,zero(rateType))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  fsalfirst = f(t,u) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*fsalfirst
      k2 = Δt*f(t+c1*Δt,u+a21*k1)
      k3 = Δt*f(t+c2*Δt,u+a31*k1+a32*k2)
      k4 = Δt*f(t+c3*Δt,u+a41*k1+a42*k2+a43*k3)
      k5 = Δt*f(t+c4*Δt,u+a51*k1+a52*k2+a53*k3+a54*k4)
      k6 = Δt*f(t+c5*Δt,u+a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
      k7 = Δt*f(t+Δt,u+a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
      utmp = u+a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7
      fsallast = f(t+Δt,utmp); k8 = Δt*fsallast
      if adaptive
        uhat   = u + bhat1*k1 + bhat2*k2 + bhat3*k3 + bhat4*k4 + bhat5*k5 + bhat6*k6 + bhat7*k7
        utilde = u + btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7 + btilde8*k8
        EEst1 = sqrt( sum(((uhat-utmp)./(abstol+max(u,utmp)*reltol)).^2) * normfactor)
        EEst2 = sqrt( sum(((utilde-utmp)./(abstol+max(u,utmp)*reltol)).^2) * normfactor)
        EEst = max(EEst1,EEst2)
      else
        u = utmp
      end
      if dense
        k[1]=k1; k[2]=k2; k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:BS5Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat2,bhat3,bhat4,bhat5,bhat6,bhat7,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8   = constructBS5(uEltypeNoUnits)
  k1::uType = similar(u)
  k2::uType = similar(u)
  k3::uType = similar(u)
  k4::uType = similar(u)
  k5::uType = similar(u)
  k6::uType = similar(u)
  k7::uType = similar(u)
  k8::uType = similar(u)
  rtmp = rateType(sizeu)
  local utilde::uType
  local uhat::uType
  local EEst2::uEltypeNoUnits
  if dense
    k = ksEltype()
    for i in 1:8
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  f(t,u,fsalfirst) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*fsalfirst
      f(t+c1*Δt,u+a21*k1,rtmp); k2=Δt*rtmp
      f(t+c2*Δt,u+a31*k1+a32*k2,rtmp); k3=Δt*rtmp
      f(t+c3*Δt,u+a41*k1+a42*k2+a43*k3,rtmp); k4=Δt*rtmp
      f(t+c4*Δt,u+a51*k1+a52*k2+a53*k3+a54*k4,rtmp); k5=Δt*rtmp
      f(t+c5*Δt,u+a61*k1+a62*k2+a63*k3+a64*k4+a65*k5,rtmp); k6=Δt*rtmp
      f(t+Δt,u+a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6,rtmp); k7=Δt*rtmp
      utmp = u+a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7
      f(t+Δt,utmp,fsallast); k8 = Δt*fsallast
      if adaptive
        uhat   = u + bhat1*k1 + bhat2*k2 + bhat3*k3 + bhat4*k4 + bhat5*k5 + bhat6*k6 + bhat7*k7
        utilde = u + btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7 + btilde8*k8
        EEst1 = sqrt( sum(((uhat-utmp)./(abstol+max(u,utmp)*reltol)).^2) * normfactor)
        EEst2 = sqrt( sum(((utilde-utmp)./(abstol+max(u,utmp)*reltol)).^2) * normfactor)
        EEst = max(EEst1,EEst2)
      else
        u = utmp
      end
      if dense
        k[1]=k1; k[2]=k2; k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:BS5,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat2,bhat3,bhat4,bhat5,bhat6,bhat7,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = constructBS5(uEltypeNoUnits)
  k1::uType = similar(u)
  k2::uType = similar(u)
  k3::uType = similar(u)
  k4::uType = similar(u)
  k5::uType = similar(u)
  k6::uType = similar(u)
  k7::uType = similar(u)
  k8::uType = similar(u); rtmp = rateType(sizeu)
  utilde = similar(u)
  uhat   = similar(u)
  local EEst2::uEltypeNoUnits
  uidx = eachindex(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); atmptilde = similar(u,uEltypeNoUnits)
  if dense
    k = ksEltype()
    for i in 1:8
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  if dense
    k[1]=k1; k[2]=k2; k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8
  end
  f(t,u,fsalfirst) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      for i in uidx
        k1[i] = Δt*fsalfirst[i]
        tmp[i] = u[i]+a21*k1[i]
      end
      f(t+c1*Δt,tmp,rtmp)
      for i in uidx
        k2[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a31*k1[i]+a32*k2[i]
      end
      f(t+c2*Δt,tmp,rtmp)
      for i in uidx
        k3[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a41*k1[i]+a42*k2[i]+a43*k3[i]
      end
      f(t+c3*Δt,tmp,rtmp)
      for i in uidx
        k4[i] = Δt*rtmp[i]
        tmp[i] = (u[i]+a51*k1[i]+a52*k2[i])+(a53*k3[i]+a54*k4[i])
      end
      f(t+c4*Δt,tmp,rtmp)
      for i in uidx
        k5[i] = Δt*rtmp[i]
        tmp[i] = (u[i]+a61*k1[i]+a62*k2[i])+(a63*k3[i]+a64*k4[i]+a65*k5[i])
      end
      f(t+c5*Δt,tmp,rtmp)
      for i in uidx
        k6[i] = Δt*rtmp[i]
        tmp[i] = (u[i]+a71*k1[i]+a72*k2[i]+a73*k3[i])+(a74*k4[i]+a75*k5[i]+a76*k6[i])
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k7[i] = Δt*rtmp[i]
        utmp[i] = (u[i]+a81*k1[i]+a83*k3[i])+(a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i])
      end
      f(t+Δt,utmp,fsallast)
      for i in uidx
        k8[i] = Δt*fsallast[i]
      end
      if adaptive
        for i in uidx
          uhat[i]   = u[i] + bhat1*k1[i] + bhat2*k2[i] + bhat3*k3[i] + bhat4*k4[i] + bhat5*k5[i] + bhat6*k6[i] + bhat7*k7[i]
          utilde[i] = u[i] + btilde1*k1[i] + btilde2*k2[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i] + btilde8*k8[i]
          atmp[i] = ((uhat[i]-utmp[i])./(abstol+max(u[i],utmp[i])*reltol))^2
          atmptilde[i] = ((utilde[i]-utmp[i])./(abstol+max(u[i],utmp[i])*reltol))^2
        end
        EEst1 = sqrt( sum(atmp) * normfactor)
        EEst2 = sqrt( sum(atmptilde) * normfactor)
        EEst = max(EEst1,EEst2)
      else
        copy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:Tsit5,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,b1,b2,b3,b4,b5,b6,b7 = constructTsit5(uEltypeNoUnits)
  local k1::uType
  local k2::uType
  local k3::uType
  local k4::uType
  local k5::uType
  local k6::uType
  local k7::uType
  local utilde::uType
  if dense
    k = ksEltype()
    for i in 1:7
      push!(k,zero(rateType))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  fsalfirst = f(t,u) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*fsalfirst
      k2 = Δt*f(t+c1*Δt,u+a21*k1)
      k3 = Δt*f(t+c2*Δt,u+a31*k1+a32*k2)
      k4 = Δt*f(t+c3*Δt,u+a41*k1+a42*k2+a43*k3)
      k5 = Δt*f(t+c4*Δt,u+a51*k1+a52*k2+a53*k3+a54*k4)
      k6 = Δt*f(t+Δt,u+a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
      utmp = u+a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6
      fsallast = f(t+Δt,utmp); k7 = Δt*fsallast
      if adaptive
        utilde = u + b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7
        EEst = abs(((utilde-utmp)/(abstol+max(u,utmp)*reltol)) * normfactor)
      else
        u = utmp
      end
      if dense
        k[1] = k1
        k[2] = k2
        k[3] = k3
        k[4] = k4
        k[5] = k5
        k[6] = k6
        k[7] = k7
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Tsit5Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,b1,b2,b3,b4,b5,b6,b7 = constructTsit5(uEltypeNoUnits)
  k1::uType = similar(u)
  k2::uType = similar(u)
  k3::uType = similar(u)
  k4::uType = similar(u)
  k5::uType = similar(u)
  k6::uType = similar(u)
  k7::uType = similar(u)
  utilde::uType = similar(u)
  rtmp = rateType(sizeu)
  if dense
    k = ksEltype()
    for i in 1:7
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  f(t,u,fsalfirst) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*fsalfirst
      f(t+c1*Δt,u+a21*k1,rtmp); k2=Δt*rtmp
      f(t+c2*Δt,u+a31*k1+a32*k2,rtmp); k3=Δt*rtmp
      f(t+c3*Δt,u+a41*k1+a42*k2+a43*k3,rtmp); k4=Δt*rtmp
      f(t+c4*Δt,u+(a51*k1+a52*k2+a53*k3+a54*k4),rtmp); k5=Δt*rtmp
      f(t+Δt,(u+a61*k1+a62*k2+a63*k3)+(a64*k4+a65*k5),rtmp); k6=Δt*rtmp
      utmp = (u+a71*k1+a72*k2+a73*k3)+(a74*k4+a75*k5+a76*k6)
      f(t+Δt,utmp,fsallast); k7 = Δt*fsallast
      if adaptive
        utilde = (u + b1*k1 + b2*k2 + b3*k3) + (b4*k4 + b5*k5 + b6*k6 + b7*k7)
        EEst = sqrt( sum(((utilde-utmp)./(abstol+max(u,utmp)*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if dense
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

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Tsit5,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,b1,b2,b3,b4,b5,b6,b7 = constructTsit5(uEltypeNoUnits)
  k1::uType = similar(u)
  k2::uType = similar(u)
  k3::uType = similar(u)
  k4::uType = similar(u)
  k5::uType = similar(u)
  k6::uType = similar(u)
  k7::uType = similar(u)
  utilde::uType = similar(u)
  uidx = eachindex(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  rtmp = rateType(sizeu)
  if dense
    k = ksEltype()
    for i in 1:7
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    # Setup k pointers
    k[1] = k1
    k[2] = k2
    k[3] = k3
    k[4] = k4
    k[5] = k5
    k[6] = k6
    k[7] = k7
  end
  f(t,u,fsalfirst) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      for i in uidx
        k1[i] = Δt*fsalfirst[i]
        tmp[i] = u[i]+a21*k1[i]
      end
      f(t+c1*Δt,tmp,rtmp)
      for i in uidx
        k2[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a31*k1[i]+a32*k2[i]
      end
      f(t+c2*Δt,tmp,rtmp)
      for i in uidx
        k3[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a41*k1[i]+a42*k2[i]+a43*k3[i]
      end
      f(t+c3*Δt,tmp,rtmp)
      for i in uidx
        k4[i] = Δt*rtmp[i]
        tmp[i] = u[i]+(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
      end
      f(t+c4*Δt,tmp,rtmp)
      for i in uidx
        k5[i] = Δt*rtmp[i]
        tmp[i] = (u[i]+a61*k1[i]+a62*k2[i]+a63*k3[i])+(a64*k4[i]+a65*k5[i])
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k6[i] = Δt*rtmp[i]
        utmp[i] = (u[i]+a71*k1[i]+a72*k2[i]+a73*k3[i])+(a74*k4[i]+a75*k5[i]+a76*k6[i])
      end
      f(t+Δt,utmp,fsallast)
      for i in uidx
        k7[i] = Δt*fsallast[i]
      end
      if adaptive
        for i in uidx
          utilde[i] = (u[i] + b1*k1[i] + b2*k2[i] + b3*k3[i]) + (b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i])
          atmp[i] = ((utilde[i]-utmp[i])./(abstol+max(u[i],utmp[i])*reltol)).^2
        end
        EEst = sqrt( sum(atmp) * normfactor)
      else
        copy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:DP5,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = constructDP5(uEltypeNoUnits)
  local k1::uType
  local k2::uType
  local k3::uType
  local k4::uType
  local k5::uType
  local k6::uType
  local k7::uType
  if dense
    k = ksEltype()
    for i in 1:6
      push!(k,zero(rateType))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  local utilde::uType
  fsalfirst = f(t,u) # Pre-start fsal
  @inbounds for T in Ts
    while t<T
      @ode_loopheader
      k1 = Δt*fsalfirst
      k2 = Δt*f(t+c1*Δt,u+a21*k1)
      k3 = Δt*f(t+c2*Δt,u+a31*k1+a32*k2)
      k4 = Δt*f(t+c3*Δt,u+a41*k1+a42*k2+a43*k3)
      k5 = Δt*f(t+c4*Δt,u+a51*k1+a52*k2+a53*k3+a54*k4)
      k6 = Δt*f(t+Δt,u+a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
      utmp = u+a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
      fsallast = f(t+Δt,utmp); k7 = Δt*fsallast
      if adaptive
        utilde = u + b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7
        EEst = abs( ((utilde-utmp)/(abstol+max(u,utmp)*reltol)) * normfactor)
      else
        u = utmp
      end
      if dense
        k[1] = k1
        # No 2
        k[2] = k3
        k[3] = k4
        k[4] = k5
        k[5] = k6
        k[6] = k7
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:DP5Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = constructDP5(uEltypeNoUnits)
  k1 = similar(u)
  k2 = similar(u)
  k3 = similar(u)
  k4 = similar(u)
  k5 = similar(u)
  k6 = similar(u)
  k7 = similar(u)
  utilde = similar(u)
  rtmp = rateType(sizeu)
  if dense
    k = ksEltype()
    for i in 1:6
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  f(t,u,fsalfirst); #Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1=Δt*fsalfirst
      f(t+c1*Δt,u+a21*k1,rtmp); k2=Δt*rtmp
      f(t+c2*Δt,u+a31*k1+a32*k2,rtmp); k3=Δt*rtmp
      f(t+c3*Δt,u+a41*k1+a42*k2+a43*k3,rtmp); k4=Δt*rtmp
      f(t+c4*Δt,u+(a51*k1+a52*k2+a53*k3+a54*k4),rtmp); k5=Δt*rtmp
      f(t+Δt,(u+a61*k1+a62*k2+a63*k3)+(a64*k4+a65*k5),rtmp); k6=Δt*rtmp
      utmp = (u+a71*k1+a73*k3)+(a74*k4+a75*k5+a76*k6)
      f(t+Δt,utmp,fsallast); k7=Δt*fsallast
      if adaptive
        utilde = (u + b1*k1 + b3*k3) + (b4*k4 + b5*k5 + b6*k6 + b7*k7)
        EEst = sqrt( sum(((utilde-utmp)./(abstol+max(u,utmp)*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if dense
        k[1] = k1
        # No 2
        k[2] = k3
        k[3] = k4
        k[4] = k5
        k[5] = k6
        k[6] = k7
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:DP5,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = constructDP5(uEltypeNoUnits)
  k1 = similar(u)
  k2 = similar(u)
  k3 = similar(u)
  k4 = similar(u)
  k5 = similar(u)
  k6 = similar(u)
  k7 = similar(u)
  utilde = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  uidx = eachindex(u)
  rtmp = rateType(sizeu)
  if dense
    k = ksEltype()
    for i in 1:6
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    # Setup k pointers
    k[1] = k1
    # No 2
    k[2] = k3
    k[3] = k4
    k[4] = k5
    k[5] = k6
    k[6] = k7
  end
  f(t,u,fsalfirst);  # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      for i in uidx
        k1[i] = Δt*fsalfirst[i]
        tmp[i] = u[i]+a21*k1[i]
      end
      f(t+c1*Δt,tmp,rtmp)
      for i in uidx
        k2[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a31*k1[i]+a32*k2[i]
      end
      f(t+c2*Δt,tmp,rtmp)
      for i in uidx
        k3[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a41*k1[i]+a42*k2[i]+a43*k3[i]
      end
      f(t+c3*Δt,tmp,rtmp)
      for i in uidx
        k4[i] = Δt*rtmp[i]
        tmp[i] =(u[i]+a51*k1[i])+(a52*k2[i]+a53*k3[i]+a54*k4[i])
      end
      f(t+c4*Δt,tmp,rtmp)
      for i in uidx
        k5[i] = Δt*rtmp[i]
        tmp[i] = (u[i]+a61*k1[i]+a62*k2[i]+a63*k3[i])+(a64*k4[i]+a65*k5[i])
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k6[i] = Δt*rtmp[i]
        utmp[i] = (u[i]+a71*k1[i]+a73*k3[i])+(a74*k4[i]+a75*k5[i]+a76*k6[i])
      end
      f(t+Δt,utmp,fsallast);
      for i in uidx
        k7[i]=Δt*fsallast[i]
      end
      if adaptive
        for i in uidx
          utilde[i] = (u[i] + b1*k1[i] + b3*k3[i]) + (b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i])
          atmp[i] = ((utilde[i]-utmp[i])/(abstol+max(u[i],utmp[i])*reltol))^2
        end
        EEst = sqrt( sum(atmp) * normfactor)
      else
        copy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:DP5Threaded,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = constructDP5(uEltypeNoUnits)
  k1::uType = similar(u)
  k2::uType = similar(u)
  k3::uType = similar(u)
  k4::uType = similar(u)
  k5::uType = similar(u)
  k6::uType = similar(u)
  k7::uType = similar(u)
  utilde = similar(u); rtmp = rateType(sizeu)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  uidx::Base.OneTo{Int64} = eachindex(u)
  if dense
    k = ksEltype()
    for i in 1:6
      push!(k,zero(rateType))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    # Setup k pointers
    k[1] = k1
    # No 2
    k[2] = k3
    k[3] = k4
    k[4] = k5
    k[5] = k6
    k[6] = k7
  end
  f(t,u,fsalfirst);  # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      dp5threaded_loop1(k1,fsalfirst,u,a21,tmp,Δt,uidx)
      f(t+c1*Δt,tmp,rtmp)
      dp5threaded_loop2(rtmp,Δt,tmp,u,a31,k1,a32,k2,uidx)
      f(t+c2*Δt,tmp,rtmp)
      dp5threaded_loop3(rtmp,Δt,tmp,u,a41,k1,a42,k2,a43,k3,uidx)
      f(t+c3*Δt,tmp,rtmp)
      dp5threaded_loop4(rtmp,Δt,tmp,u,a51,k1,a52,k2,a53,k3,a54,k4,uidx)
      f(t+c4*Δt,tmp,rtmp)
      dp5threaded_loop5(rtmp,Δt,tmp,u,a61,k1,a62,k2,a63,k3,a64,k4,a65,k5,uidx)
      f(t+Δt,tmp,rtmp)
      dp5threaded_loop6(rtmp,Δt,utmp,u,a71,k1,a73,k3,a74,k4,a75,k5,a76,k6,uidx)
      f(t+Δt,utmp,fsallast)
      dp5threaded_loop7(k7,Δt,fsallast,uidx)
      if adaptive
        dp5threaded_adaptiveloop(utilde,u,b1,k1,b3,k3,b4,k4,b5,k5,b6,k6,b7,k7,atmp,utmp,abstol,reltol,uidx)
        EEst = sqrt( sum(atmp) * normfactor)
      else
        copy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Float64,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:DP5Threaded,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = constructDP5(uEltypeNoUnits)
  k1::uType = similar(u)
  k2::uType = similar(u)
  k3::uType = similar(u)
  k4::uType = similar(u)
  k5::uType = similar(u)
  k6::uType = similar(u)
  k7::uType = similar(u)
  utilde = similar(u); rtmp = rateType(sizeu)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  uidx::Base.OneTo{Int64} = eachindex(u)
  if dense
    k = ksEltype()
    for i in 1:6
      push!(k,zero(rateType))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    # Setup k pointers
    k[1] = k1
    # No 2
    k[2] = k3
    k[3] = k4
    k[4] = k5
    k[5] = k6
    k[6] = k7
  end
  f(t,u,fsalfirst);  # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      dp5threaded_loop1(k1,fsalfirst,u,tmp,Δt,uidx)
      f(t+c1*Δt,tmp,rtmp);
      dp5threaded_loop2(rtmp,Δt,tmp,u,k1,k2,uidx)
      f(t+c2*Δt,tmp,rtmp);
      dp5threaded_loop3(rtmp,Δt,tmp,u,k1,k2,k3,uidx)
      f(t+c3*Δt,tmp,rtmp);
      dp5threaded_loop4(rtmp,Δt,tmp,u,k1,k2,k3,k4,uidx)
      f(t+c4*Δt,tmp,rtmp);
      dp5threaded_loop5(rtmp,Δt,tmp,u,k1,k2,k3,k4,k5,uidx)
      f(t+Δt,tmp,rtmp);
      dp5threaded_loop6(rtmp,Δt,utmp,u,k1,k3,k4,k5,k6,uidx)
      f(t+Δt,utmp,fsallast);
      dp5threaded_loop7(k7,Δt,fsallast,uidx)
      if adaptive
        dp5threaded_adaptiveloop(utilde,u,k1,k3,k4,k5,k6,k7,atmp,utmp,abstol,reltol,uidx)
        EEst = sqrt( sum(atmp) * normfactor)
      else
        copy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

@noinline function dp5threaded_loop1(k1,fsalfirst,u,a21,tmp,Δt,uidx)
  Threads.@threads for i in uidx
    k1[i] = Δt*fsalfirst[i]
    tmp[i] = u[i]+a21*k1[i]
  end
end

@noinline function dp5threaded_loop2(rtmp,Δt,tmp,u,a31,k1,a32,k2,uidx)
  Threads.@threads for i in uidx
    k2[i] = Δt*rtmp[i]
    tmp[i] = u[i]+a31*k1[i]+a32*k2[i]
  end
end

@noinline function dp5threaded_loop3(rtmp,Δt,tmp,u,a41,k1,a42,k2,a43,k3,uidx)
  Threads.@threads for i in uidx
    k3[i] = Δt*rtmp[i]
    tmp[i] = u[i]+a41*k1[i]+a42*k2[i]+a43*k3[i]
  end
end

@noinline function dp5threaded_loop4(rtmp,Δt,tmp,u,a51,k1,a52,k2,a53,k3,a54,k4,uidx)
  Threads.@threads for i in uidx
    k4[i] = Δt*rtmp[i]
    tmp[i] =(u[i]+a51*k1[i])+(a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
end

@noinline function dp5threaded_loop5(rtmp,Δt,tmp,u,a61,k1,a62,k2,a63,k3,a64,k4,a65,k5,uidx)
  Threads.@threads for i in uidx
    k5[i] = Δt*rtmp[i]
    tmp[i] = (u[i]+a61*k1[i]+a62*k2[i]+a63*k3[i])+(a64*k4[i]+a65*k5[i])
  end
end

@noinline function dp5threaded_loop6(rtmp,Δt,utmp,u,a71,k1,a73,k3,a74,k4,a75,k5,a76,k6,uidx)
  Threads.@threads for i in uidx
    k6[i] = Δt*rtmp[i]
    utmp[i] = (u[i]+a71*k1[i]+a73*k3[i])+(a74*k4[i]+a75*k5[i]+a76*k6[i])
  end
end

@noinline function dp5threaded_loop7(k7,Δt,fsallast,uidx)
  Threads.@threads for i in uidx
    k7[i]=Δt*fsallast[i]
  end
end

@noinline function dp5threaded_adaptiveloop(utilde,u,b1,k1,b3,k3,b4,k4,b5,k5,b6,k6,b7,k7,tmp,utmp,abstol,reltol,uidx)
  Threads.@threads for i in uidx
    utilde[i] = (u[i] + b1*k1[i] + b3*k3[i]) + (b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i])
    tmp[i] = ((utilde[i]-utmp[i])/(abstol+max(u[i],utmp[i])*reltol))^2
  end
end

@noinline function dp5threaded_loop1{T<:Float64,N}(k1::Array{T,N},fsalfirst,u,tmp,Δt,uidx)
  Threads.@threads for i in uidx
    k1[i] = Δt*fsalfirst[i]
    tmp[i] = u[i]+0.2*k1[i]
  end
end

@noinline function dp5threaded_loop2{T<:Float64,N}(rtmp,Δt,tmp::Array{T,N},u,k1,k2,uidx)
  Threads.@threads for i in uidx
    k2[i] = Δt*rtmp[i]
    tmp[i] = u[i]+0.075*k1[i]+0.225*k2[i]
  end
end

@noinline function dp5threaded_loop3{T<:Float64,N}(rtmp,Δt,tmp::Array{T,N},u,k1,k2,k3,uidx)
  Threads.@threads for i in uidx
    k3[i] = Δt*rtmp[i]
    tmp[i] = u[i]+0.9777777777777777*k1[i]-3.7333333333333334*k2[i]+3.5555555555555554*k3[i]
  end
end

@noinline function dp5threaded_loop4{T<:Float64,N}(rtmp,Δt,tmp::Array{T,N},u,k1,k2,k3,k4,uidx)
  Threads.@threads for i in uidx
    k4[i] = Δt*rtmp[i]
    tmp[i] =(u[i]+ 2.9525986892242035*k1[i])+(-11.595793324188385*k2[i]+9.822892851699436*k3[i]-0.2908093278463649*k4[i])
  end
end

@noinline function dp5threaded_loop5{T<:Float64,N}(rtmp,Δt,tmp::Array{T,N},u,k1,k2,k3,k4,k5,uidx)
  Threads.@threads for i in uidx
    k5[i]  = Δt*rtmp[i]
    tmp[i] = (u[i]+2.8462752525252526*k1[i]-10.757575757575758*k2[i]+8.906422717743473*k3[i])+(0.2784090909090909*k4[i]-0.2735313036020583*k5[i])
  end
end

@noinline function dp5threaded_loop6{T<:Float64,N}(rtmp,Δt,utmp::Array{T,N},u,k1,k3,k4,k5,k6,uidx)
  Threads.@threads for i in uidx
    k6[i] = Δt*rtmp[i]
    utmp[i] = (u[i]+0.09114583333333333*k1[i]+0.44923629829290207*k3[i])+(0.6510416666666666*k4[i]-0.322376179245283*k5[i]+0.13095238095238096*k6[i])
  end
end

@noinline function dp5threaded_loop7{T<:Float64,N}(k7::Array{T,N},Δt,fsallast,uidx)
  Threads.@threads for i in uidx
    k7[i]=Δt*fsallast[i]
  end
end

@noinline function dp5threaded_adaptiveloop{T<:Float64,N}(utilde::Array{T,N},u,k1,k3,k4,k5,k6,k7,tmp,utmp,abstol,reltol,uidx)
  Threads.@threads for i in uidx
    utilde[i] = (u[i] + 0.08991319444444444*k1[i] + 0.4534890685834082*k3[i]) + (0.6140625*k4[i] - 0.2715123820754717*k5[i] + 0.08904761904761904*k6[i] + 0.025*k7[i])
    tmp[i] = ((utilde[i]-utmp[i])/(abstol+max(u[i],utmp[i])*reltol))^2
  end
end


function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:Vern6,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,b1,b4,b5,b6,b7,b8,b9= constructVern6(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType;
  local k5::uType; local k6::uType; local k7::uType; local k8::uType;
  local k9::uType;
  local utilde::uType; fsalfirst = f(t,u) # Pre-start fsal
  if dense
    k = ksEltype()
    for i in 1:9
      push!(k,zero(rateType))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*fsalfirst
      k2 = Δt*f(t+c1*Δt,u+a21*k1)
      k3 = Δt*f(t+c2*Δt,u+a31*k1+a32*k2)
      k4 = Δt*f(t+c3*Δt,u+a41*k1       +a43*k3)
      k5 = Δt*f(t+c4*Δt,u+a51*k1       +a53*k3+a54*k4)
      k6 = Δt*f(t+c5*Δt,u+a61*k1       +a63*k3+a64*k4+a65*k5)
      k7 = Δt*f(t+c6*Δt,u+a71*k1       +a73*k3+a74*k4+a75*k5+a76*k6)
      k8 = Δt*f(t+Δt,u+a81*k1       +a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
      utmp =    u+a91*k1              +a94*k4+a95*k5+a96*k6+a97*k7+a98*k8
      fsallast = f(t+Δt,utmp); k9 = Δt*fsallast
      if adaptive
        utilde = u + b1*k1 + b4*k4 + b5*k5 + b6*k6 + b7*k7 + b8*k8 + b9*k9
        EEst = abs( ((utilde-utmp)/(abstol+max(u,utmp)*reltol)) * normfactor)
      else
        u = utmp
      end
      if dense
        k[1]=k1; k[2]=k2; k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Vern6Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,b1,b4,b5,b6,b7,b8,b9= constructVern6(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u) ; k3 = similar(u); k4 = similar(u)
  k5 = similar(u); k6 = similar(u) ; k7 = similar(u); k8 = similar(u)
  k9 = similar(u)
  if dense
    k = ksEltype()
    for i in 1:9
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  utilde = similar(u); rtmp = rateType(sizeu)
  f(t,u,fsalfirst) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*fsalfirst
      f(t+c1*Δt,u+a21*k1,rtmp); k2=Δt*rtmp
      f(t+c2*Δt,u+a31*k1+a32*k2,rtmp); k3=Δt*rtmp
      f(t+c3*Δt,u+a41*k1       +a43*k3,rtmp); k4=Δt*rtmp
      f(t+c4*Δt,u+a51*k1       +a53*k3+a54*k4,rtmp); k5=Δt*rtmp
      f(t+c5*Δt,u+a61*k1       +a63*k3+a64*k4+a65*k5,rtmp); k6=Δt*rtmp
      f(t+c6*Δt,u+a71*k1       +a73*k3+a74*k4+a75*k5+a76*k6,rtmp); k7=Δt*rtmp
      f(t+Δt,u+a81*k1       +a83*k3+a84*k4+a85*k5+a86*k6+a87*k7,rtmp); k8=Δt*rtmp
      utmp=u+a91*k1              +a94*k4+a95*k5+a96*k6+a97*k7+a98*k8
      f(t+Δt,utmp,fsallast); k9 = Δt*fsallast
      if adaptive
        utilde = u + b1*k1 + b4*k4 + b5*k5 + b6*k6 + b7*k7 + b8*k8 + b9*k9
        EEst = sqrt( sum(((utilde-utmp)./(abstol+max(u,utmp)*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if dense
        k[1]=k1; k[2]=k2; k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8; k[9]=k9
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Vern6,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,b1,b4,b5,b6,b7,b8,b9= constructVern6(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u); rtmp = rateType(sizeu)
  k5 = similar(u); k6 = similar(u); k7 = similar(u); k8 = similar(u);
  k9 = similar(u);
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits); uidx = eachindex(u)
  f(t,u,fsalfirst) # Pre-start fsal
  if dense
    k = ksEltype()
    for i in 1:9
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    k[1]=k1; k[2]=k2; k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9 # Set the pointers
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      for i in uidx
        k1[i] = Δt*fsalfirst[i]
        tmp[i] = u[i]+a21*k1[i]
      end
      f(t+c1*Δt,tmp,rtmp)
      for i in uidx
        k2[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a31*k1[i]+a32*k2[i]
      end
      f(t+c2*Δt,tmp,rtmp)
      for i in uidx
        k3[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a41*k1[i]+a43*k3[i]
      end
      f(t+c3*Δt,tmp,rtmp)
      for i in uidx
        k4[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a51*k1[i]+a53*k3[i]+a54*k4[i]
      end
      f(t+c4*Δt,tmp,rtmp)
      for i in uidx
        k5[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a61*k1[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]
      end
      f(t+c5*Δt,tmp,rtmp)
      for i in uidx
        k6[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]
      end
      f(t+c6*Δt,tmp,rtmp)
      for i in uidx
        k7[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i]
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k8[i]=Δt*rtmp[i]
        utmp[i]=u[i]+a91*k1[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k8[i]
      end
      f(t+Δt,utmp,fsallast)
      for i in uidx
        k9[i] = Δt*fsallast[i]
      end
      if adaptive
        for i in uidx
          utilde[i] = u[i] + b1*k1[i] + b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i] + b8*k8[i] + b9*k9[i]
          atmp[i] = ((utilde[i]-utmp[i])/(abstol+max(u[i],utmp[i])*reltol))^2
        end
        EEst = sqrt( sum(atmp) * normfactor)
      else
        copy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:Vern7,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c2,c3,c4,c5,c6,c7,c8,a021,a031,a032,a041,a043,a051,a053,a054,a061,a063,a064,a065,a071,a073,a074,a075,a076,a081,a083,a084,a085,a086,a087,a091,a093,a094,a095,a096,a097,a098,a101,a103,a104,a105,a106,a107,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat10= constructVern7(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType;
  local k5::uType; local k6::uType; local k7::uType; local k8::uType;
  local k9::uType; local k10::uType; local utilde::uType; local update::uType
  if dense
    k = ksEltype()
    for i in 1:10
      push!(k,zero(rateType))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*f(t,u)
      k2 = Δt*f(t+c2*Δt,u+a021*k1)
      k3 = Δt*f(t+c3*Δt,u+a031*k1+a032*k2)
      k4 = Δt*f(t+c4*Δt,u+a041*k1       +a043*k3)
      k5 = Δt*f(t+c5*Δt,u+a051*k1       +a053*k3+a054*k4)
      k6 = Δt*f(t+c6*Δt,u+a061*k1       +a063*k3+a064*k4+a065*k5)
      k7 = Δt*f(t+c7*Δt,u+a071*k1       +a073*k3+a074*k4+a075*k5+a076*k6)
      k8 = Δt*f(t+c8*Δt,u+a081*k1       +a083*k3+a084*k4+a085*k5+a086*k6+a087*k7)
      k9 = Δt*f(t+Δt,u+a091*k1          +a093*k3+a094*k4+a095*k5+a096*k6+a097*k7+a098*k8)
      k10= Δt*f(t+Δt,u+a101*k1          +a103*k3+a104*k4+a105*k5+a106*k6+a107*k7)
      update = k1*b1 + k4*b4 + k5*b5 + k6*b6 + k7*b7 + k8*b8 + k9*b9
      utmp = u + update
      if adaptive
        EEst = abs( ((update - bhat1*k1 - bhat4*k4 - bhat5*k5 - bhat6*k6 - bhat7*k7 - bhat10*k10)/(abstol+max(u,utmp)*reltol)) * normfactor)
      else
        u = utmp
      end
      if dense
        k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Vern7Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c2,c3,c4,c5,c6,c7,c8,a021,a031,a032,a041,a043,a051,a053,a054,a061,a063,a064,a065,a071,a073,a074,a075,a076,a081,a083,a084,a085,a086,a087,a091,a093,a094,a095,a096,a097,a098,a101,a103,a104,a105,a106,a107,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat10= constructVern7(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u);
  k5 = similar(u); k6 = similar(u); k7 = similar(u); k8 = similar(u);
  rtmp = rateType(sizeu)
  k9 = similar(u); k10 = similar(u); utilde = similar(u); update = similar(u)
  if dense
    k = ksEltype()
    for i in 1:10
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp); k1 = Δt*rtmp
      f(t+c2*Δt,u+a021*k1,rtmp); k2 = Δt*rtmp
      f(t+c3*Δt,u+a031*k1+a032*k2,rtmp); k3 = Δt*rtmp
      f(t+c4*Δt,u+a041*k1       +a043*k3,rtmp); k4 = Δt*rtmp
      f(t+c5*Δt,u+a051*k1       +a053*k3+a054*k4,rtmp); k5 = Δt*rtmp
      f(t+c6*Δt,u+a061*k1       +a063*k3+a064*k4+a065*k5,rtmp); k6 = Δt*rtmp
      f(t+c7*Δt,u+a071*k1       +a073*k3+a074*k4+a075*k5+a076*k6,rtmp); k7 = Δt*rtmp
      f(t+c8*Δt,u+a081*k1       +a083*k3+a084*k4+a085*k5+a086*k6+a087*k7,rtmp); k8 = Δt*rtmp
      f(t+Δt,u+a091*k1          +a093*k3+a094*k4+a095*k5+a096*k6+a097*k7+a098*k8,rtmp); k9 = Δt*rtmp
      f(t+Δt,u+a101*k1          +a103*k3+a104*k4+a105*k5+a106*k6+a107*k7,rtmp); k10 = Δt*rtmp
      update = k1*b1 + k4*b4 + k5*b5 + k6*b6 + k7*b7 + k8*b8 + k9*b9
      utmp = u + update
      if adaptive
        EEst = sqrt( sum(((update - bhat1*k1 - bhat4*k4 - bhat5*k5 - bhat6*k6 - bhat7*k7 - bhat10*k10)/(abstol+max(u,utmp)*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if dense
        k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Vern7,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c2,c3,c4,c5,c6,c7,c8,a021,a031,a032,a041,a043,a051,a053,a054,a061,a063,a064,a065,a071,a073,a074,a075,a076,a081,a083,a084,a085,a086,a087,a091,a093,a094,a095,a096,a097,a098,a101,a103,a104,a105,a106,a107,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat10= constructVern7(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u);
  k5 = similar(u); k6 = similar(u); k7 = similar(u); k8 = similar(u);
  k9 = similar(u); k10 = similar(u); utilde = similar(u); update = similar(u)
  uidx = eachindex(u); tmp = similar(u); rtmp = rateType(sizeu); atmp = similar(u,uEltypeNoUnits)
  if dense
    k = ksEltype()
    for i in 1:10
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10 # Setup pointers
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp)
      for i in uidx
        k1[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a021*k1[i]
      end
      f(t+c2*Δt,tmp,rtmp)
      for i in uidx
        k2[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a031*k1[i]+a032*k2[i]
      end
      f(t+c3*Δt,tmp,rtmp)
      for i in uidx
        k3[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a041*k1[i]+a043*k3[i]
      end
      f(t+c4*Δt,tmp,rtmp)
      for i in uidx
        k4[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a051*k1[i]+a053*k3[i]+a054*k4[i]
      end
      f(t+c5*Δt,tmp,rtmp)
      for i in uidx
        k5[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a061*k1[i]+a063*k3[i]+a064*k4[i]+a065*k5[i]
      end
      f(t+c6*Δt,tmp,rtmp)
      for i in uidx
        k6[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a071*k1[i]+a073*k3[i]+a074*k4[i]+a075*k5[i]+a076*k6[i]
      end
      f(t+c7*Δt,tmp,rtmp)
      for i in uidx
        k7[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a081*k1[i]+a083*k3[i]+a084*k4[i]+a085*k5[i]+a086*k6[i]+a087*k7[i]
      end
      f(t+c8*Δt,tmp,rtmp)
      for i in uidx
        k8[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a091*k1[i]+a093*k3[i]+a094*k4[i]+a095*k5[i]+a096*k6[i]+a097*k7[i]+a098*k8[i]
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k9[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a101*k1[i]+a103*k3[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+a107*k7[i]
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k10[i] = Δt*rtmp[i]
        update[i] = k1[i]*b1 + k4[i]*b4 + k5[i]*b5 + k6[i]*b6 + k7[i]*b7 + k8[i]*b8 + k9[i]*b9
        utmp[i] = u[i] + update[i]
      end
      if adaptive
        for i in uidx
          atmp[i] = ((update[i] - bhat1*k1[i] - bhat4*k4[i] - bhat5*k5[i] - bhat6*k6[i] - bhat7*k7[i] - bhat10*k10[i])/(abstol+max(u[i],utmp[i])*reltol))^2
        end
        EEst = sqrt(sum(atmp)  * normfactor)
      else
        copy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:Vern8,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13= constructVern8(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType;
  local k5::uType; local k6::uType; local k7::uType; local k8::uType;
  local k9::uType; local k10::uType; local utilde::uType; local update::uType
  if dense
    k = ksEltype()
    for i in 1:13
      push!(k,zero(rateType))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*f(t,u)
      k2 = Δt*f(t+c2*Δt ,u+a0201*k1)
      k3 = Δt*f(t+c3*Δt ,u+a0301*k1+a0302*k2)
      k4 = Δt*f(t+c4*Δt ,u+a0401*k1       +a0403*k3)
      k5 = Δt*f(t+c5*Δt ,u+a0501*k1       +a0503*k3+a0504*k4)
      k6 = Δt*f(t+c6*Δt ,u+a0601*k1                +a0604*k4+a0605*k5)
      k7 = Δt*f(t+c7*Δt ,u+a0701*k1                +a0704*k4+a0705*k5+a0706*k6)
      k8 = Δt*f(t+c8*Δt ,u+a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7)
      k9 = Δt*f(t+c9*Δt ,u+a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8)
      k10= Δt*f(t+c10*Δt,u+a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9)
      k11= Δt*f(t+c11*Δt,u+a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10)
      k12= Δt*f(t+    Δt,u+a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11)
      k13= Δt*f(t+    Δt,u+a1301*k1                +a1304*k4+a1305*k5+a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10)
      update = k1*b1 + k6*b6 + k7*b7 + k8*b8 + k9*b9 + k10*b10 + k11*b11 + k12*b12
      utmp = u + update
      if adaptive
        EEst = abs( ((update - bhat1*k1 - bhat6*k6 - bhat7*k7 - bhat8*k8 - bhat9*k9 - bhat10*k10 - bhat13*k13)/(abstol+max(u,utmp)*reltol)) * normfactor)
      else
        u = utmp
      end
      if dense
        k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10;k[11]=k11;k[12]=k12;k[13]=k13
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Vern8Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13= constructVern8(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u);
  k5 = similar(u); k6 = similar(u); k7 = similar(u); k8 = similar(u);
  k9 = similar(u); k10 = similar(u); k11 = similar(u); k12 = similar(u); k13 = similar(u)
  utilde = similar(u); update = similar(u);
  rtmp = rateType(sizeu)
  if dense
    k = ksEltype()
    for i in 1:13
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp); k1 = Δt*rtmp
      f(t+c2*Δt ,u+a0201*k1,rtmp); k2 = Δt*rtmp
      f(t+c3*Δt ,u+a0301*k1+a0302*k2,rtmp); k3 = Δt*rtmp
      f(t+c4*Δt ,u+a0401*k1       +a0403*k3,rtmp); k4 = Δt*rtmp
      f(t+c5*Δt ,u+a0501*k1       +a0503*k3+a0504*k4,rtmp); k5 = Δt*rtmp
      f(t+c6*Δt ,u+a0601*k1                +a0604*k4+a0605*k5,rtmp); k6 = Δt*rtmp
      f(t+c7*Δt ,u+a0701*k1                +a0704*k4+a0705*k5+a0706*k6,rtmp); k7 = Δt*rtmp
      f(t+c8*Δt ,u+a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7,rtmp); k8 = Δt*rtmp
      f(t+c9*Δt ,u+a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8,rtmp); k9 = Δt*rtmp
      f(t+c10*Δt,u+a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9,rtmp); k10 = Δt*rtmp
      f(t+c11*Δt,u+a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10,rtmp); k11 = Δt*rtmp
      f(t+    Δt,u+a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11,rtmp); k12 = Δt*rtmp
      f(t+    Δt,u+a1301*k1                +a1304*k4+a1305*k5+a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10,rtmp); k13 = Δt*rtmp
      update = k1*b1 + k6*b6 + k7*b7 + k8*b8 + k9*b9 + k10*b10 + k11*b11 + k12*b12
      utmp = u + update
      if adaptive
        EEst = sqrt( sum(((update - bhat1*k1 - bhat6*k6 - bhat7*k7 - bhat8*k8 - bhat9*k9 - bhat10*k10 - bhat13*k13)/(abstol+max(u,utmp)*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if dense
        k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10;k[11]=k11;k[12]=k12;k[13]=k13
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Vern8,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13= constructVern8(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u);
  k5 = similar(u); k6 = similar(u); k7 = similar(u); k8 = similar(u); tmp = similar(u)
  k9 = similar(u); k10 = similar(u); k11 = similar(u); k12 = similar(u); k13 = similar(u)
  utilde = similar(u); update = similar(u);
  rtmp = rateType(sizeu); uidx = eachindex(u); atmp = similar(u,uEltypeNoUnits)
  if dense
    k = ksEltype()
    for i in 1:13
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10;k[11]=k11;k[12]=k12;k[13]=k13 # Setup pointers
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp)
      for i in uidx
        k1[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a0201*k1[i]
      end
      f(t+c2*Δt,tmp,rtmp)
      for i in uidx
        k2[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a0301*k1[i]+a0302*k2[i]
      end
      f(t+c3*Δt,tmp,rtmp)
      for i in uidx
        k3[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a0401*k1[i]+a0403*k3[i]
      end
      f(t+c4*Δt,tmp,rtmp)
      for i in uidx
        k4[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a0501*k1[i]+a0503*k3[i]+a0504*k4[i]
      end
      f(t+c5*Δt,tmp,rtmp)
      for i in uidx
        k5[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a0601*k1[i]+a0604*k4[i]+a0605*k5[i]
      end
      f(t+c6*Δt,tmp,rtmp)
      for i in uidx
        k6[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a0701*k1[i]+a0704*k4[i]+a0705*k5[i]+a0706*k6[i]
      end
      f(t+c7*Δt,tmp,rtmp)
      for i in uidx
        k7[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a0801*k1[i]+a0804*k4[i]+a0805*k5[i]+a0806*k6[i]+a0807*k7[i]
      end
      f(t+c8*Δt,tmp,rtmp)
      for i in uidx
        k8[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a0901*k1[i]+a0904*k4[i]+a0905*k5[i]+a0906*k6[i]+a0907*k7[i]+a0908*k8[i]
      end
      f(t+c9*Δt,tmp,rtmp)
      for i in uidx
        k9[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a1001*k1[i]+a1004*k4[i]+a1005*k5[i]+a1006*k6[i]+a1007*k7[i]+a1008*k8[i]+a1009*k9[i]
      end
      f(t+c10*Δt,tmp,rtmp)
      for i in uidx
        k10[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a1101*k1[i]+a1104*k4[i]+a1105*k5[i]+a1106*k6[i]+a1107*k7[i]+a1108*k8[i]+a1109*k9[i]+a1110*k10[i]
      end
      f(t+c11*Δt,tmp,rtmp)
      for i in uidx
        k11[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a1201*k1[i]+a1204*k4[i]+a1205*k5[i]+a1206*k6[i]+a1207*k7[i]+a1208*k8[i]+a1209*k9[i]+a1210*k10[i]+a1211*k11[i]
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k12[i] = Δt*rtmp[i]
        tmp[i] = u[i]+a1301*k1[i]+a1304*k4[i]+a1305*k5[i]+a1306*k6[i]+a1307*k7[i]+a1308*k8[i]+a1309*k9[i]+a1310*k10[i]
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k13[i] = Δt*rtmp[i]
        update[i] = k1[i]*b1 + k6[i]*b6 + k7[i]*b7 + k8[i]*b8 + k9[i]*b9 + k10[i]*b10 + k11[i]*b11 + k12[i]*b12
        utmp[i] = u[i] + update[i]
      end
      if adaptive
        for i in uidx
          atmp[i] = ((update[i] - bhat1*k1[i] - bhat6*k6[i] - bhat7*k7[i] - bhat8*k8[i] - bhat9*k9[i] - bhat10*k10[i] - bhat13*k13[i])/(abstol+max(u[i],utmp[i])*reltol)).^2
        end
        EEst = sqrt( sum(atmp) * normfactor)
      else
        copy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:TanYam7,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a93,a94,a95,a96,a97,a98,a101,a103,a104,a105,a106,a107,a108,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat8,bhat10 = constructTanYam7(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType;
  local k5::uType; local k6::uType; local k7::uType; local k8::uType;
  local k9::uType; local k10::uType;
  local utilde::uType;
  if dense
    pop!(ks) # Get rid of the one it starts with
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k = f(t,u)
      k1 = Δt*k
      k2 = Δt*f(t+c1*Δt,u+a21*k1)
      k3 = Δt*f(t+c2*Δt,u+a31*k1+a32*k2)
      k4 = Δt*f(t+c3*Δt,u+a41*k1       +a43*k3)
      k5 = Δt*f(t+c4*Δt,u+a51*k1       +a53*k3+a54*k4)
      k6 = Δt*f(t+c5*Δt,u+a61*k1       +a63*k3+a64*k4+a65*k5)
      k7 = Δt*f(t+c6*Δt,u+a71*k1       +a73*k3+a74*k4+a75*k5+a76*k6)
      k8 = Δt*f(t+c7*Δt,u+a81*k1       +a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
      k9 = Δt*f(t+Δt,u+a91*k1       +a93*k3+a94*k4+a95*k5+a96*k6+a97*k7+a98*k8)
      k10= Δt*f(t+Δt,u+a101*k1      +a103*k3+a104*k4+a105*k5+a106*k6+a107*k7+a108*k8)
      utmp = u + k1*b1+k4*b4+k5*b5+k6*b6+k7*b7+k8*b8+k9*b9
      if adaptive
        utilde = u + k1*bhat1+k4*bhat4+k5*bhat5+k6*bhat6+k7*bhat7+k8*bhat8+k10*bhat10
        EEst = abs( ((utilde-utmp)/(abstol+max(u,utmp)*reltol)) * normfactor)
      else
        u = utmp
      end
      @ode_numberloopfooter
    end
  end
  if dense
    k = f(t,u)
    push!(ks,k)
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:TanYam7Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a93,a94,a95,a96,a97,a98,a101,a103,a104,a105,a106,a107,a108,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat8,bhat10 = constructTanYam7(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u) ; k3 = similar(u); k4 = similar(u)
  k5 = similar(u); k6 = similar(u) ; k7 = similar(u); k8 = similar(u)
  k9 = similar(u); k10= similar(u) ; k  = rateType(sizeu)
  utilde = similar(u); rtmp = rateType(sizeu);
  if dense
    pop!(ks) # Get rid of the one it starts with
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k); k1=Δt*k
      f(t+c1*Δt,u+a21*k1,rtmp); k2=Δt*rtmp
      f(t+c2*Δt,u+a31*k1+a32*k2,rtmp); k3=Δt*rtmp
      f(t+c3*Δt,u+a41*k1       +a43*k3,rtmp); k4=Δt*rtmp
      f(t+c4*Δt,u+a51*k1       +a53*k3+a54*k4,rtmp); k5=Δt*rtmp
      f(t+c5*Δt,u+a61*k1       +a63*k3+a64*k4+a65*k5,rtmp); k6=Δt*rtmp
      f(t+c6*Δt,u+a71*k1       +a73*k3+a74*k4+a75*k5+a76*k6,rtmp); k7=Δt*rtmp
      f(t+c7*Δt,u+a81*k1       +a83*k3+a84*k4+a85*k5+a86*k6+a87*k7,rtmp); k8=Δt*rtmp
      f(t+Δt,u+a91*k1       +a93*k3+a94*k4+a95*k5+a96*k6+a97*k7+a98*k8,rtmp); k9=Δt*rtmp
      f(t+Δt,u+a101*k1      +a103*k3+a104*k4+a105*k5+a106*k6+a107*k7+a108*k8,rtmp); k10=Δt*rtmp
      utmp = u + k1*b1+k4*b4+k5*b5+k6*b6+k7*b7+k8*b8+k9*b9
      if adaptive
        utilde = u + k1*bhat1+k4*bhat4+k5*bhat5+k6*bhat6+k7*bhat7+k8*bhat8+k10*bhat10
        EEst = sqrt( sum(((utilde-utmp)./(abstol+max(u,utmp)*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  if dense
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:TanYam7,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a93,a94,a95,a96,a97,a98,a101,a103,a104,a105,a106,a107,a108,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat8,bhat10 = constructTanYam7(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u) ; k3 = similar(u); k4 = similar(u)
  k5 = similar(u); k6 = similar(u) ; k7 = similar(u); k8 = similar(u)
  k9 = similar(u); k10= similar(u) ; k = rateType(sizeu)
  utilde = similar(u); uidx = eachindex(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  rtmp = rateType(sizeu)
  if dense
    pop!(ks) # Get rid of the one it starts with
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k)
      for i in uidx
        k1[i]=Δt*k[i]
        tmp[i] = u[i]+a21*k1[i]
      end
      f(t+c1*Δt,tmp,rtmp)
      for i in uidx
        k2[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a31*k1[i]+a32*k2[i]
      end
      f(t+c2*Δt,tmp,rtmp)
      for i in uidx
        k3[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a41*k1[i]+a43*k3[i]
      end
      f(t+c3*Δt,tmp,rtmp)
      for i in uidx
        k4[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a51*k1[i]+a53*k3[i]+a54*k4[i]
      end
      f(t+c4*Δt,tmp,rtmp)
      for i in uidx
        k5[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a61*k1[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]
      end
      f(t+c5*Δt,tmp,rtmp)
      for i in uidx
        k6[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]
      end
      f(t+c6*Δt,tmp,rtmp)
      for i in uidx
        k7[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i]
      end
      f(t+c7*Δt,tmp,rtmp)
      for i in uidx
        k8[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a91*k1[i]+a93*k3[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k8[i]
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k9[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a101*k1[i]+a103*k3[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+a107*k7[i]+a108*k8[i]
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k10[i]=Δt*rtmp[i]
        utmp[i] = u[i] + k1[i]*b1+k4[i]*b4+k5[i]*b5+k6[i]*b6+k7[i]*b7+k8[i]*b8+k9[i]*b9
      end
      if adaptive
        for i in uidx
          utilde[i] = u[i] + k1[i]*bhat1+k4[i]*bhat4+k5[i]*bhat5+k6[i]*bhat6+k7[i]*bhat7+k8[i]*bhat8+k10[i]*bhat10
          atmp[i] = ((utilde[i]-utmp[i])/(abstol+max(u[i],utmp[i])*reltol))^2
        end
        EEst = sqrt( sum(atmp) * normfactor)
      else
        copy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  if dense
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:DP8,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,b1,b6,b7,b8,b9,b10,b11,b12,bhh1,bhh2,bhh3,er1,er6,er7,er8,er9,er10,er11,er12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = constructDP8(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType;
  local k5::uType; local k6::uType; local k7::uType; local k8::uType;
  local k9::uType; local k10::uType; local k11::uType; local k12::uType;
  local k13::uType; local utilde::uType; local udiff::uType; local bspl::uType
  if dense
    c14,c15,c16,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = DP8Interp(uEltypeNoUnits)
    d401,d406,d407,d408,d409,d410,d411,d412,d413,d414,d415,d416,d501,d506,d507,d508,d509,d510,d511,d512,d513,d514,d515,d516,d601,d606,d607,d608,d609,d610,d611,d612,d613,d614,d615,d616,d701,d706,d707,d708,d709,d710,d711,d712,d713,d714,d715,d716 = DP8Interp_polyweights(uEltypeNoUnits)
    if dense
      k = ksEltype()
      for i in 1:7
        push!(k,zero(rateType))
      end
      push!(ks,deepcopy(k)) #Initialize ks
    end
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*f(t,u)
      k2 = Δt*f(t+c2*Δt,u+a0201*k1)
      k3 = Δt*f(t+c3*Δt,u+a0301*k1+a0302*k2)
      k4 = Δt*f(t+c4*Δt,u+a0401*k1       +a0403*k3)
      k5 = Δt*f(t+c5*Δt,u+a0501*k1       +a0503*k3+a0504*k4)
      k6 = Δt*f(t+c6*Δt,u+a0601*k1                +a0604*k4+a0605*k5)
      k7 = Δt*f(t+c7*Δt,u+a0701*k1                +a0704*k4+a0705*k5+a0706*k6)
      k8 = Δt*f(t+c8*Δt,u+a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7)
      k9 = Δt*f(t+c9*Δt,u+a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8)
      k10 =Δt*f(t+c10*Δt,u+a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9)
      k11= Δt*f(t+c11*Δt,u+a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10)
      k12= Δt*f(t+Δt,u+a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11)
      update = b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12
      utmp = u + update
      if adaptive
        err5 = abs((k1*er1 + k6*er6 + k7*er7 + k8*er8 + k9*er9 + k10*er10 + k11*er11 + k12*er12)/(abstol+max(u,utmp)*reltol) * normfactor) # Order 5
        err3 = abs((update - bhh1*k1 - bhh2*k9 - bhh3*k12)/(abstol+max(u,utmp)*reltol) * normfactor) # Order 3
        err52 = err5*err5
        EEst = err52/sqrt(err52 + 0.01*err3*err3)
      else
        u = utmp
      end
      if dense
        k13 = update
        k14 = Δt*f(t+c14*Δt,u+a1401*k1         +a1407*k7+a1408*k8+a1409*k9+a1410*k10+a1411*k11+a1412*k12+a1413*k13)
        k15 = Δt*f(t+c15*Δt,u+a1501*k1+a1506*k6+a1507*k7+a1508*k8                   +a1511*k11+a1512*k12+a1513*k13+a1514*k14)
        k16 = Δt*f(t+c16*Δt,u+a1601*k1+a1606*k6+a1607*k7+a1608*k8+a1609*k9                              +a1613*k13+a1614*k14+a1615*k15)
        udiff = utmp - u
        k[1] = udiff
        bspl = k1 - udiff
        k[2] = bspl
        k[3] = udiff - update - bspl
        k[4] = d401*k1+d406*k6+d407*k7+d408*k8+d409*k9+d410*k10+d411*k11+d412*k12+d413*k13+d414*k14+d415*k15+d416*k16
        k[5] = d501*k1+d506*k6+d507*k7+d508*k8+d509*k9+d510*k10+d511*k11+d512*k12+d513*k13+d514*k14+d515*k15+d516*k16
        k[6] = d601*k1+d606*k6+d607*k7+d608*k8+d609*k9+d610*k10+d611*k11+d612*k12+d613*k13+d614*k14+d615*k15+d616*k16
        k[7] = d701*k1+d706*k6+d707*k7+d708*k8+d709*k9+d710*k10+d711*k11+d712*k12+d713*k13+d714*k14+d715*k15+d716*k16
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:DP8Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,b1,b6,b7,b8,b9,b10,b11,b12,bhh1,bhh2,bhh3,er1,er6,er7,er8,er9,er10,er11,er12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = constructDP8(uEltypeNoUnits)
  k1 = similar(u); k2  = similar(u); k3  = similar(u);  k4 = similar(u)
  k5 = similar(u); k6  = similar(u); k7  = similar(u);  k8 = similar(u)
  k9 = similar(u); k10 = similar(u); k11 = similar(u); k12 = similar(u)
  local k13::uType; local k14::uType; local k15::uType; local k16::uType;
  local udiff::uType; local bspl::uType
  utilde = similar(u); err5 = similar(u); err3 = similar(u)
  rtmp = rateType(sizeu)
  if dense
    c14,c15,c16,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = DP8Interp(uEltypeNoUnits)
    d401,d406,d407,d408,d409,d410,d411,d412,d413,d414,d415,d416,d501,d506,d507,d508,d509,d510,d511,d512,d513,d514,d515,d516,d601,d606,d607,d608,d609,d610,d611,d612,d613,d614,d615,d616,d701,d706,d707,d708,d709,d710,d711,d712,d713,d714,d715,d716 = DP8Interp_polyweights(uEltypeNoUnits)
    if dense
      k = ksEltype()
      for i in 1:7
        push!(k,rateType(sizeu))
      end
      push!(ks,deepcopy(k)) #Initialize ks
    end
    k14 = similar(u)
    k15 = similar(u)
    k16 = similar(u)
    udiff = similar(u)
    bspl = similar(u)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t, u,rtmp); k1=Δt*rtmp
      f(t+c2*Δt, u+a0201*k1,rtmp); k2=Δt*rtmp
      f(t+c3*Δt, u+a0301*k1+a0302*k2,rtmp); k3=Δt*rtmp
      f(t+c4*Δt, u+a0401*k1       +a0403*k3,rtmp); k4=Δt*rtmp
      f(t+c5*Δt, u+a0501*k1       +a0503*k3+a0504*k4,rtmp); k5=Δt*rtmp
      f(t+c6*Δt, u+a0601*k1                +a0604*k4+a0605*k5,rtmp); k6=Δt*rtmp
      f(t+c7*Δt, u+a0701*k1                +a0704*k4+a0705*k5+a0706*k6,rtmp); k7=Δt*rtmp
      f(t+c8*Δt, u+a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7,rtmp); k8=Δt*rtmp
      f(t+c9*Δt, u+a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8,rtmp); k9=Δt*rtmp
      f(t+c10*Δt,u+a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9,rtmp); k10=Δt*rtmp
      f(t+c11*Δt,u+a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10,rtmp); k11=Δt*rtmp
      f(t+Δt,u+a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11,rtmp); k12=Δt*rtmp
      update = b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12
      utmp = u + update
      if adaptive
        err5 = sqrt(sum(((k1*er1 + k6*er6 + k7*er7 + k8*er8 + k9*er9 + k10*er10 + k11*er11 + k12*er12)./(abstol+max(u,utmp)*reltol)).^2) * normfactor) # Order 5
        err3 = sqrt(sum(((update - bhh1*k1 - bhh2*k9 - bhh3*k12)./(abstol+max(u,utmp)*reltol)).^2) * normfactor) # Order 3
        err52 = err5*err5
        EEst = err52/sqrt(err52 + 0.01*err3*err3)
      else
        u = utmp
      end
      if dense
        k13 = update
        f(t+c14*Δt,u+a1401*k1         +a1407*k7+a1408*k8+a1409*k9+a1410*k10+a1411*k11+a1412*k12+a1413*k13,rtmp); k14 = Δt*rtmp
        f(t+c15*Δt,u+a1501*k1+a1506*k6+a1507*k7+a1508*k8                   +a1511*k11+a1512*k12+a1513*k13+a1514*k14,rtmp); k15 = Δt*rtmp
        f(t+c16*Δt,u+a1601*k1+a1606*k6+a1607*k7+a1608*k8+a1609*k9                              +a1613*k13+a1614*k14+a1615*k15,rtmp); k16 = Δt*rtmp
        udiff = utmp - u
        k[1] = udiff
        bspl = k1 - udiff
        k[2] = bspl
        k[3] = udiff - update - bspl
        k[4] = d401*k1+d406*k6+d407*k7+d408*k8+d409*k9+d410*k10+d411*k11+d412*k12+d413*k13+d414*k14+d415*k15+d416*k16
        k[5] = d501*k1+d506*k6+d507*k7+d508*k8+d509*k9+d510*k10+d511*k11+d512*k12+d513*k13+d514*k14+d515*k15+d516*k16
        k[6] = d601*k1+d606*k6+d607*k7+d608*k8+d609*k9+d610*k10+d611*k11+d612*k12+d613*k13+d614*k14+d615*k15+d616*k16
        k[7] = d701*k1+d706*k6+d707*k7+d708*k8+d709*k9+d710*k10+d711*k11+d712*k12+d713*k13+d714*k14+d715*k15+d716*k16
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:DP8,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,b1,b6,b7,b8,b9,b10,b11,b12,bhh1,bhh2,bhh3,er1,er6,er7,er8,er9,er10,er11,er12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211 = constructDP8(uEltypeNoUnits)
  k1 = similar(u); k2  = similar(u); k3  = similar(u);  k4 = similar(u)
  k5 = similar(u); k6  = similar(u); k7  = similar(u);  k8 = similar(u)
  k9 = similar(u); k10 = similar(u); k11 = similar(u); k12 = similar(u)
  k13 = similar(u); utilde = similar(u); err5 = similar(u); err3 = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); uidx = eachindex(u); atmp2 = similar(u,uEltypeNoUnits); update = similar(u)
  rtmp = rateType(sizeu)
  local k13::uType; local k14::uType; local k15::uType; local k16::uType;
  local udiff::uType; local bspl::uType
  if dense
    c14,c15,c16,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = DP8Interp(uEltypeNoUnits)
    d401,d406,d407,d408,d409,d410,d411,d412,d413,d414,d415,d416,d501,d506,d507,d508,d509,d510,d511,d512,d513,d514,d515,d516,d601,d606,d607,d608,d609,d610,d611,d612,d613,d614,d615,d616,d701,d706,d707,d708,d709,d710,d711,d712,d713,d714,d715,d716 = DP8Interp_polyweights(uEltypeNoUnits)
    if dense
      k = ksEltype()
      for i in 1:7
        push!(k,rateType(sizeu))
      end
      push!(ks,deepcopy(k)) #Initialize ks
    end
    k13 = update
    k14 = similar(u)
    k15 = similar(u)
    k16 = similar(u)
    udiff = similar(u)
    bspl = similar(u)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t, u,rtmp)
      for i in uidx
        k1[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0201*k1[i]
      end
      f(t+c2*Δt,tmp,rtmp)
      for i in uidx
        k2[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0301*k1[i]+a0302*k2[i]
      end
      f(t+c3*Δt,tmp,rtmp)
      for i in uidx
        k3[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0401*k1[i]+a0403*k3[i]
      end
      f(t+c4*Δt,tmp,rtmp)
      for i in uidx
        k4[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0501*k1[i]+a0503*k3[i]+a0504*k4[i]
      end
      f(t+c5*Δt,tmp,rtmp)
      for i in uidx
        k5[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0601*k1[i]+a0604*k4[i]+a0605*k5[i]
      end
      f(t+c6*Δt,tmp,rtmp)
      for i in uidx
        k6[i]=Δt*rtmp[i]
        tmp[i]=u[i]+a0701*k1[i]+a0704*k4[i]+a0705*k5[i]+a0706*k6[i]
      end
      f(t+c7*Δt,tmp,rtmp)
      for i in uidx
        k7[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0801*k1[i]+a0804*k4[i]+a0805*k5[i]+a0806*k6[i]+a0807*k7[i]
      end
      f(t+c8*Δt,tmp,rtmp)
      for i in uidx
        k8[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0901*k1[i]+a0904*k4[i]+a0905*k5[i]+a0906*k6[i]+a0907*k7[i]+a0908*k8[i]
      end
      f(t+c9*Δt,tmp,rtmp)
      for i in uidx
        k9[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1001*k1[i]+a1004*k4[i]+a1005*k5[i]+a1006*k6[i]+a1007*k7[i]+a1008*k8[i]+a1009*k9[i]
      end
      f(t+c10*Δt,tmp,rtmp)
      for i in uidx
        k10[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1101*k1[i]+a1104*k4[i]+a1105*k5[i]+a1106*k6[i]+a1107*k7[i]+a1108*k8[i]+a1109*k9[i]+a1110*k10[i]
      end
      f(t+c11*Δt,tmp,rtmp)
      for i in uidx
        k11[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1201*k1[i]+a1204*k4[i]+a1205*k5[i]+a1206*k6[i]+a1207*k7[i]+a1208*k8[i]+a1209*k9[i]+a1210*k10[i]+a1211*k11[i]
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k12[i]=Δt*rtmp[i]
        update[i] = b1*k1[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+b11*k11[i]+b12*k12[i]
        utmp[i] = u[i] + update[i]
      end
      if adaptive
        for i in uidx
          atmp[i] = ((k1[i]*er1 + k6[i]*er6 + k7[i]*er7 + k8[i]*er8 + k9[i]*er9 + k10[i]*er10 + k11[i]*er11 + k12[i]*er12)/(abstol+max(u[i],utmp[i])*reltol))^2
          atmp2[i]= ((update[i] - bhh1*k1[i] - bhh2*k9[i] - bhh3*k12[i])/(abstol+max(u[i],utmp[i])*reltol))^2
        end
        err5 = sqrt( sum(atmp)  * normfactor) # Order 5
        err3 = sqrt( sum(atmp2) * normfactor) # Order 3
        err52 = err5*err5
        EEst = err52/sqrt(err52 + 0.01*err3*err3)
      else
        copy!(u, utmp)
      end
      if dense
        for i in uidx
          tmp[i] = u[i]+a1401*k1[i]+a1407*k7[i]+a1408*k8[i]+a1409*k9[i]+a1410*k10[i]+a1411*k11[i]+a1412*k12[i]+a1413*k13[i]
        end
        f(t+c14*Δt,tmp,rtmp)
        for i in uidx
          k14[i] = Δt*rtmp[i]
          tmp[i] = u[i]+a1501*k1[i]+a1506*k6[i]+a1507*k7[i]+a1508*k8[i]+a1511*k11[i]+a1512*k12[i]+a1513*k13[i]+a1514*k14[i]
        end
        f(t+c15*Δt,tmp,rtmp)
        for i in uidx
          k15[i] = Δt*rtmp[i]
          tmp[i] = u[i]+a1601*k1[i]+a1606*k6[i]+a1607*k7[i]+a1608*k8[i]+a1609*k9[i]+a1613*k13[i]+a1614*k14[i]+a1615*k15[i]
        end
        f(t+c16*Δt,tmp,rtmp)
        for i in uidx
          k16[i] = Δt*rtmp[i]
          udiff[i] = utmp[i] - u[i]
          k[1][i] = udiff[i]
          bspl[i] = k1[i] - udiff[i]
          k[2][i] = bspl[i]
          k[3][i] = udiff[i] - update[i] - bspl[i]
          k[4][i] = d401*k1[i]+d406*k6[i]+d407*k7[i]+d408*k8[i]+d409*k9[i]+d410*k10[i]+d411*k11[i]+d412*k12[i]+d413*k13[i]+d414*k14[i]+d415*k15[i]+d416*k16[i]
          k[5][i] = d501*k1[i]+d506*k6[i]+d507*k7[i]+d508*k8[i]+d509*k9[i]+d510*k10[i]+d511*k11[i]+d512*k12[i]+d513*k13[i]+d514*k14[i]+d515*k15[i]+d516*k16[i]
          k[6][i] = d601*k1[i]+d606*k6[i]+d607*k7[i]+d608*k8[i]+d609*k9[i]+d610*k10[i]+d611*k11[i]+d612*k12[i]+d613*k13[i]+d614*k14[i]+d615*k15[i]+d616*k16[i]
          k[7][i] = d701*k1[i]+d706*k6[i]+d707*k7[i]+d708*k8[i]+d709*k9[i]+d710*k10[i]+d711*k11[i]+d712*k12[i]+d713*k13[i]+d714*k14[i]+d715*k15[i]+d716*k16[i]
        end
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:TsitPap8,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13 = constructTsitPap8(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType;
  local k5::uType; local k6::uType; local k7::uType; local k8::uType;
  local k9::uType; local k10::uType; local k11::uType; local k12::uType;
  local k13::uType; local utilde::uType;
  if dense
    pop!(ks) # Take out the initial
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k = f(t,u)
      k1 = Δt*k
      k2 = Δt*f(t+c1*Δt,u+a0201*k1)
      k3 = Δt*f(t+c2*Δt,u+a0301*k1+a0302*k2)
      k4 = Δt*f(t+c3*Δt,u+a0401*k1       +a0403*k3)
      k5 = Δt*f(t+c4*Δt,u+a0501*k1       +a0503*k3+a0504*k4)
      k6 = Δt*f(t+c5*Δt,u+a0601*k1                +a0604*k4+a0605*k5)
      k7 = Δt*f(t+c6*Δt,u+a0701*k1                +a0704*k4+a0705*k5+a0706*k6)
      k8 = Δt*f(t+c7*Δt,u+a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7)
      k9 = Δt*f(t+c8*Δt,u+a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8)
      k10 =Δt*f(t+c9*Δt,u+a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9)
      k11= Δt*f(t+c10*Δt,u+a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10)
      k12= Δt*f(t+Δt,u+a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11)
      k13= Δt*f(t+Δt,u+a1301*k1                +a1304*k4+a1305*k5+a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10)
      update = b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12
      utmp = u + update
      if adaptive
        EEst = abs((update - k1*bhat1 - k6*bhat6 - k7*bhat7 - k8*bhat8 - k9*bhat9 - k10*bhat10 - k13*bhat13)/(abstol+max(u,utmp)*reltol) * normfactor) # Order 5
      else
        u = utmp
      end
      @ode_numberloopfooter
    end
  end
  if dense
    k = f(t,u)
    push!(ks,k)
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:TsitPap8Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13 = constructTsitPap8(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u)
  k5 = similar(u); k6 = similar(u); k7 = similar(u); k8 = similar(u)
  k9 = similar(u); k10 = similar(u); k11 = similar(u); k12 = similar(u)
  k13::uType = similar(u); utilde = similar(u); k = rateType(sizeu)
  if dense
    pop!(ks)
  end
  rtmp = rateType(sizeu)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k); k1=Δt*k
      f(t+c1*Δt,u+a0201*k1,rtmp); k2=Δt*rtmp
      f(t+c2*Δt,u+a0301*k1+a0302*k2,rtmp); k3=Δt*rtmp
      f(t+c3*Δt,u+a0401*k1       +a0403*k3,rtmp); k4=Δt*rtmp
      f(t+c4*Δt,u+a0501*k1       +a0503*k3+a0504*k4,rtmp); k5=Δt*rtmp
      f(t+c5*Δt,u+a0601*k1                +a0604*k4+a0605*k5,rtmp); k6=Δt*rtmp
      f(t+c6*Δt,u+a0701*k1                +a0704*k4+a0705*k5+a0706*k6,rtmp); k7=Δt*rtmp
      f(t+c7*Δt,u+a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7,rtmp); k8=Δt*rtmp
      f(t+c8*Δt,u+a0901*k1                +a0904*k4+a0905*k5+a0906*k6+a0907*k7+a0908*k8,rtmp); k9=Δt*rtmp
      f(t+c9*Δt,u+a1001*k1                +a1004*k4+a1005*k5+a1006*k6+a1007*k7+a1008*k8+a1009*k9,rtmp); k10=Δt*rtmp
      f(t+c10*Δt,u+a1101*k1                +a1104*k4+a1105*k5+a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10,rtmp); k11=Δt*rtmp
      f(t+Δt,u+a1201*k1                +a1204*k4+a1205*k5+a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11,rtmp); k12=Δt*rtmp
      f(t+Δt,u+a1301*k1                +a1304*k4+a1305*k5+a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10,rtmp); k13=Δt*rtmp
      update = b1*k1+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k11+b12*k12
      utmp = u + update
      if adaptive
        EEst = sqrt(sum(((update - k1*bhat1 - k6*bhat6 - k7*bhat7 - k8*bhat8 - k9*bhat9 - k10*bhat10 - k13*bhat13)./(abstol+max(u,utmp)*reltol)).^2) * normfactor) # Order 5
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  if dense
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:TsitPap8,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13 = constructTsitPap8(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u)
  k5 = similar(u); k6 = similar(u); k7 = similar(u); k8 = similar(u)
  k9 = similar(u); k10 = similar(u); k11 = similar(u); k12 = similar(u)
  k13 = similar(u); rtmp = rateType(sizeu); update = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); uidx = eachindex(u)
  k13::uType; utilde = similar(u); k = rateType(sizeu)
  if dense
    pop!(ks)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k)
      for i in uidx
        k1[i]=Δt*k[i]
        tmp[i] = u[i]+a0201*k1[i]
      end
      f(t+c1*Δt,tmp,rtmp)
      for i in uidx
        k2[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0301*k1[i]+a0302*k2[i]
      end
      f(t+c2*Δt,tmp,rtmp)
      for i in uidx
        k3[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0401*k1[i]+a0403*k3[i]
      end
      f(t+c3*Δt,tmp,rtmp)
      for i in uidx
        k4[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0501*k1[i]+a0503*k3[i]+a0504*k4[i]
      end
      f(t+c4*Δt,tmp,rtmp)
      for i in uidx
        k5[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0601*k1[i]+a0604*k4[i]+a0605*k5[i]
      end
      f(t+c5*Δt,tmp,rtmp)
      for i in uidx
        k6[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0701*k1[i]+a0704*k4[i]+a0705*k5[i]+a0706*k6[i]
      end
      f(t+c6*Δt,tmp,rtmp)
      for i in uidx
        k7[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0801*k1[i]+a0804*k4[i]+a0805*k5[i]+a0806*k6[i]+a0807*k7[i]
      end
      f(t+c7*Δt,tmp,rtmp)
      for i in uidx
        k8[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0901*k1[i]+a0904*k4[i]+a0905*k5[i]+a0906*k6[i]+a0907*k7[i]+a0908*k8[i]
      end
      f(t+c8*Δt,tmp,rtmp)
      for i in uidx
        k9[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1001*k1[i]+a1004*k4[i]+a1005*k5[i]+a1006*k6[i]+a1007*k7[i]+a1008*k8[i]+a1009*k9[i]
      end
      f(t+c9*Δt,tmp,rtmp)
      for i in uidx
        k10[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1101*k1[i]+a1104*k4[i]+a1105*k5[i]+a1106*k6[i]+a1107*k7[i]+a1108*k8[i]+a1109*k9[i]+a1110*k10[i]
      end
      f(t+c10*Δt,tmp,rtmp)
      for i in uidx
        k11[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1201*k1[i]+a1204*k4[i]+a1205*k5[i]+a1206*k6[i]+a1207*k7[i]+a1208*k8[i]+a1209*k9[i]+a1210*k10[i]+a1211*k11[i]
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k12[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1301*k1[i]+a1304*k4[i]+a1305*k5[i]+a1306*k6[i]+a1307*k7[i]+a1308*k8[i]+a1309*k9[i]+a1310*k10[i]
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k13[i]=Δt*rtmp[i]
        update[i] = b1*k1[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+b11*k11[i]+b12*k12[i]
        utmp[i] = u[i] + update[i]
      end
      if adaptive
        for i in uidx
          atmp[i] = ((update[i] - k1[i]*bhat1 - k6[i]*bhat6 - k7[i]*bhat7 - k8[i]*bhat8 - k9[i]*bhat9 - k10[i]*bhat10 - k13[i]*bhat13)/(abstol+max(u[i],utmp[i])*reltol))^2
        end
        EEst = sqrt(sum(atmp) * normfactor) # Order 5
      else
        copy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  if dense
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:Vern9,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0806,a0807,a0901,a0906,a0907,a0908,a1001,a1006,a1007,a1008,a1009,a1101,a1106,a1107,a1108,a1109,a1110,a1201,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1401,a1406,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,b1,b8,b9,b10,b11,b12,b13,b14,b15,bhat1,bhat8,bhat9,bhat10,bhat11,bhat12,bhat13,bhat16 = constructVern9(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType;
  local k5::uType; local k6::uType; local k7::uType; local k8::uType;
  local k9::uType; local k10::uType; local k11::uType; local k12::uType;
  local k13::uType; local k14::uType; local k15::uType; local k16::uType;
  local utilde::uType; local update::uType
  if dense
    k = ksEltype()
    for i in 1:16
      push!(k,zero(rateType))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*f(t,u)
      k2 = Δt*f(t+c1*Δt,u+a0201*k1)
      k3 = Δt*f(t+c2*Δt,u+a0301*k1+a0302*k2)
      k4 = Δt*f(t+c3*Δt,u+a0401*k1       +a0403*k3)
      k5 = Δt*f(t+c4*Δt,u+a0501*k1       +a0503*k3+a0504*k4)
      k6 = Δt*f(t+c5*Δt,u+a0601*k1                +a0604*k4+a0605*k5)
      k7 = Δt*f(t+c6*Δt,u+a0701*k1                +a0704*k4+a0705*k5+a0706*k6)
      k8 = Δt*f(t+c7*Δt,u+a0801*k1                                  +a0806*k6+a0807*k7)
      k9 = Δt*f(t+c8*Δt,u+a0901*k1                                  +a0906*k6+a0907*k7+a0908*k8)
      k10 =Δt*f(t+c9*Δt,u+a1001*k1                                  +a1006*k6+a1007*k7+a1008*k8+a1009*k9)
      k11= Δt*f(t+c10*Δt,u+a1101*k1                                  +a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10)
      k12= Δt*f(t+c11*Δt,u+a1201*k1                                  +a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11)
      k13= Δt*f(t+c12*Δt,u+a1301*k1                                  +a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10+a1311*k11+a1312*k12)
      k14= Δt*f(t+c13*Δt,u+a1401*k1                                  +a1406*k6+a1407*k7+a1408*k8+a1409*k9+a1410*k10+a1411*k11+a1412*k12+a1413*k13)
      k15= Δt*f(t+Δt,u+a1501*k1                                  +a1506*k6+a1507*k7+a1508*k8+a1509*k9+a1510*k10+a1511*k11+a1512*k12+a1513*k13+a1514*k14)
      k16= Δt*f(t+Δt,u+a1601*k1                                  +a1606*k6+a1607*k7+a1608*k8+a1609*k9+a1610*k10+a1611*k11+a1612*k12+a1613*k13)
      update = k1*b1+k8*b8+k9*b9+k10*b10+k11*b11+k12*b12+k13*b13+k14*b14+k15*b15
      utmp = u + update
      if adaptive
        EEst = abs((update - k1*bhat1 - k8*bhat8 - k9*bhat9 - k10*bhat10 - k11*bhat11 - k12*bhat12 - k13*bhat13 - k16*bhat16)/(abstol+max(u,utmp)*reltol) * normfactor)
      else
        u = utmp
      end
      if dense
        k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10;k[11]=k11;k[12]=k12;k[13]=k13;k[14]=k14;k[15]=k15;k[16]=k16
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Vern9Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0806,a0807,a0901,a0906,a0907,a0908,a1001,a1006,a1007,a1008,a1009,a1101,a1106,a1107,a1108,a1109,a1110,a1201,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1401,a1406,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,b1,b8,b9,b10,b11,b12,b13,b14,b15,bhat1,bhat8,bhat9,bhat10,bhat11,bhat12,bhat13,bhat16 = constructVern9(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u);k3 = similar(u); k4 = similar(u);
  k5 = similar(u); k6 = similar(u);k7 = similar(u); k8 = similar(u);
  k9 = similar(u); k10 = similar(u); k11 = similar(u); k12 = similar(u);
  k13 = similar(u); k14 = similar(u); k15 = similar(u); k16 =similar(u);
  utilde = similar(u); update = similar(u)
  rtmp = rateType(sizeu)
  if dense
    k = ksEltype()
    for i in 1:16
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp); k1=Δt*rtmp
      f(t+c1*Δt,u+a0201*k1,rtmp); k2=Δt*rtmp
      f(t+c2*Δt,u+a0301*k1+a0302*k2,rtmp); k3=Δt*rtmp
      f(t+c3*Δt,u+a0401*k1       +a0403*k3,rtmp); k4=Δt*rtmp
      f(t+c4*Δt,u+a0501*k1       +a0503*k3+a0504*k4,rtmp); k5=Δt*rtmp
      f(t+c5*Δt,u+a0601*k1                +a0604*k4+a0605*k5,rtmp); k6=Δt*rtmp
      f(t+c6*Δt,u+a0701*k1                +a0704*k4+a0705*k5+a0706*k6,rtmp); k7=Δt*rtmp
      f(t+c7*Δt,u+a0801*k1                                  +a0806*k6+a0807*k7,rtmp); k8=Δt*rtmp
      f(t+c8*Δt,u+a0901*k1                                  +a0906*k6+a0907*k7+a0908*k8,rtmp); k9=Δt*rtmp
      f(t+c9*Δt,u+a1001*k1                                  +a1006*k6+a1007*k7+a1008*k8+a1009*k9,rtmp); k10=Δt*rtmp
      f(t+c10*Δt,u+a1101*k1                                  +a1106*k6+a1107*k7+a1108*k8+a1109*k9+a1110*k10,rtmp); k11=Δt*rtmp
      f(t+c11*Δt,u+a1201*k1                                  +a1206*k6+a1207*k7+a1208*k8+a1209*k9+a1210*k10+a1211*k11,rtmp); k12=Δt*rtmp
      f(t+c12*Δt,u+a1301*k1                                  +a1306*k6+a1307*k7+a1308*k8+a1309*k9+a1310*k10+a1311*k11+a1312*k12,rtmp); k13=Δt*rtmp
      f(t+c13*Δt,u+a1401*k1                                  +a1406*k6+a1407*k7+a1408*k8+a1409*k9+a1410*k10+a1411*k11+a1412*k12+a1413*k13,rtmp); k14=Δt*rtmp
      f(t+Δt,u+a1501*k1                                  +a1506*k6+a1507*k7+a1508*k8+a1509*k9+a1510*k10+a1511*k11+a1512*k12+a1513*k13+a1514*k14,rtmp); k15=Δt*rtmp
      f(t+Δt,u+a1601*k1                                  +a1606*k6+a1607*k7+a1608*k8+a1609*k9+a1610*k10+a1611*k11+a1612*k12+a1613*k13,rtmp); k16=Δt*rtmp
      update = k1*b1+k8*b8+k9*b9+k10*b10+k11*b11+k12*b12+k13*b13+k14*b14+k15*b15
      utmp = u + update
      if adaptive
        EEst = sqrt(sum(((update - k1*bhat1 - k8*bhat8 - k9*bhat9 - k10*bhat10 - k11*bhat11 - k12*bhat12 - k13*bhat13 - k16*bhat16)./(abstol+max(u,utmp)*reltol)).^2) * normfactor)
      else
        u = utmp
      end
      if dense
        k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10;k[11]=k11;k[12]=k12;k[13]=k13;k[14]=k14;k[15]=k15;k[16]=k16
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Vern9,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0806,a0807,a0901,a0906,a0907,a0908,a1001,a1006,a1007,a1008,a1009,a1101,a1106,a1107,a1108,a1109,a1110,a1201,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1401,a1406,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,b1,b8,b9,b10,b11,b12,b13,b14,b15,bhat1,bhat8,bhat9,bhat10,bhat11,bhat12,bhat13,bhat16 = constructVern9(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u);k3 = similar(u); k4 = similar(u);
  k5 = similar(u); k6 = similar(u);k7 = similar(u); k8 = similar(u);
  k9 = similar(u); k10 = similar(u); k11 = similar(u); k12 = similar(u); update = similar(u)
  k13 = similar(u); k14 = similar(u); k15 = similar(u); k16 =similar(u); rtmp = rateType(sizeu)
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits); uidx = eachindex(u)
  if dense
    k = ksEltype()
    for i in 1:16
      push!(k,rateType(sizeu))
    end
    push!(ks,deepcopy(k)) #Initialize ks
    k[1]=k1;k[2]=k2;k[3]=k3;k[4]=k4;k[5]=k5;k[6]=k6;k[7]=k7;k[8]=k8;k[9]=k9;k[10]=k10;k[11]=k11;k[12]=k12;k[13]=k13;k[14]=k14;k[15]=k15;k[16]=k16 # Setup pointers
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp)
      for i in uidx
        k1[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0201*k1[i]
      end
      f(t+c1*Δt,tmp,rtmp)
      for i in uidx
        k2[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0301*k1[i]+a0302*k2[i]
      end
      f(t+c2*Δt,tmp,rtmp)
      for i in uidx
        k3[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0401*k1[i]+a0403*k3[i]
      end
      f(t+c3*Δt,tmp,rtmp)
      for i in uidx
        k4[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0501*k1[i]+a0503*k3[i]+a0504*k4[i]
      end
      f(t+c4*Δt,tmp,rtmp)
      for i in uidx
        k5[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0601*k1[i]+a0604*k4[i]+a0605*k5[i]
      end
      f(t+c5*Δt,tmp,rtmp)
      for i in uidx
        k6[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0701*k1[i]+a0704*k4[i]+a0705*k5[i]+a0706*k6[i]
      end
      f(t+c6*Δt,tmp,rtmp)
      for i in uidx
        k7[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0801*k1[i]+a0806*k6[i]+a0807*k7[i]
      end
      f(t+c7*Δt,tmp,rtmp)
      for i in uidx
        k8[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a0901*k1[i]+a0906*k6[i]+a0907*k7[i]+a0908*k8[i]
      end
      f(t+c8*Δt,tmp,rtmp)
      for i in uidx
        k9[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1001*k1[i]+a1006*k6[i]+a1007*k7[i]+a1008*k8[i]+a1009*k9[i]
      end
      f(t+c9*Δt,tmp,rtmp)
      for i in uidx
        k10[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1101*k1[i]+a1106*k6[i]+a1107*k7[i]+a1108*k8[i]+a1109*k9[i]+a1110*k10[i]
      end
      f(t+c10*Δt,tmp,rtmp)
      for i in uidx
        k11[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1201*k1[i]+a1206*k6[i]+a1207*k7[i]+a1208*k8[i]+a1209*k9[i]+a1210*k10[i]+a1211*k11[i]
      end
      f(t+c11*Δt,tmp,rtmp)
      for i in uidx
        k12[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1301*k1[i]+a1306*k6[i]+a1307*k7[i]+a1308*k8[i]+a1309*k9[i]+a1310*k10[i]+a1311*k11[i]+a1312*k12[i]
      end
      f(t+c12*Δt,tmp,rtmp)
      for i in uidx
        k13[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1401*k1[i]+a1406*k6[i]+a1407*k7[i]+a1408*k8[i]+a1409*k9[i]+a1410*k10[i]+a1411*k11[i]+a1412*k12[i]+a1413*k13[i]
      end
      f(t+c13*Δt,tmp,rtmp)
      for i in uidx
        k14[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1501*k1[i]+a1506*k6[i]+a1507*k7[i]+a1508*k8[i]+a1509*k9[i]+a1510*k10[i]+a1511*k11[i]+a1512*k12[i]+a1513*k13[i]+a1514*k14[i]
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k15[i]=Δt*rtmp[i]
        tmp[i] = u[i]+a1601*k1[i]+a1606*k6[i]+a1607*k7[i]+a1608*k8[i]+a1609*k9[i]+a1610*k10[i]+a1611*k11[i]+a1612*k12[i]+a1613*k13[i]
      end
      f(t+Δt,tmp,rtmp)
      for i in uidx
        k16[i]=Δt*rtmp[i]
        update[i] = k1[i]*b1+k8[i]*b8+k9[i]*b9+k10[i]*b10+k11[i]*b11+k12[i]*b12+k13[i]*b13+k14[i]*b14+k15[i]*b15
        utmp[i] = u[i] + update[i]
      end
      if adaptive
        for i in uidx
          atmp[i] = ((update[i] - k1[i]*bhat1 - k8[i]*bhat8 - k9[i]*bhat9 - k10[i]*bhat10 - k11[i]*bhat11 - k12[i]*bhat12 - k13[i]*bhat13 - k16[i]*bhat16)/(abstol+max(u[i],utmp[i])*reltol))^2
        end
        EEst = sqrt(sum(atmp) * normfactor) # Order 5
      else
        copy!(u, utmp)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Feagin10Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16 = constructFeagin10(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u); k5 = similar(u)
  k6 = similar(u); k7 = similar(u); k8 = similar(u); k9 = similar(u); k10 = similar(u)
  k11 = similar(u); k12 = similar(u); k13 = similar(u); k14 = similar(u)
  k15 = similar(u); k16 = similar(u); k17 = similar(u); k = rateType(sizeu)
  if dense
    pop!(ks)
  end
  update = similar(u)
  utmp = similar(u)
  uidx = eachindex(u)
  rtmp = rateType(sizeu)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k); k1=Δt*k
      f(t + c1*Δt,u + a0100*k1,rtmp); k2=Δt*rtmp
      f(t + c2*Δt ,u + a0200*k1 + a0201*k2,rtmp); k3=Δt*rtmp
      f(t + c3*Δt,u + a0300*k1              + a0302*k3,rtmp); k4=Δt*rtmp
      f(t + c4*Δt,u + a0400*k1              + a0402*k3 + a0403*k4,rtmp); k5=Δt*rtmp
      f(t + c5*Δt,u + a0500*k1                           + a0503*k4 + a0504*k5,rtmp); k6=Δt*rtmp
      f(t + c6*Δt,u + a0600*k1                           + a0603*k4 + a0604*k5 + a0605*k6,rtmp); k7=Δt*rtmp
      f(t + c7*Δt,u + a0700*k1                                        + a0704*k5 + a0705*k6 + a0706*k7,rtmp); k8=Δt*rtmp
      f(t + c8*Δt,u + a0800*k1                                                     + a0805*k6 + a0806*k7 + a0807*k8,rtmp); k9=Δt*rtmp
      f(t + c9*Δt,u + a0900*k1                                                     + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9,rtmp); k10=Δt*rtmp
      f(t + c10*Δt,u + a1000*k1                                                     + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10,rtmp); k11=Δt*rtmp
      f(t + c11*Δt,u + a1100*k1                                                     + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11,rtmp); k12=Δt*rtmp
      f(t + c12*Δt,u + a1200*k1                           + a1203*k4 + a1204*k5 + a1205*k6 + a1206*k7 + a1207*k8 + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12,rtmp); k13=Δt*rtmp
      f(t + c13*Δt,u + a1300*k1              + a1302*k3 + a1303*k4              + a1305*k6 + a1306*k7 + a1307*k8 + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13,rtmp); k14=Δt*rtmp
      f(t + c14*Δt,u + a1400*k1 + a1401*k2                           + a1404*k5              + a1406*k7 +                                                                     a1412*k13 + a1413*k14,rtmp); k15=Δt*rtmp
      f(t + c15*Δt,u + a1500*k1              + a1502*k3                                                                                                                                                     + a1514*k15,rtmp); k16=Δt*rtmp
      f(t + c16*Δt,u + a1600*k1 + a1601*k2 + a1602*k3              + a1604*k5 + a1605*k6 + a1606*k7 + a1607*k8 + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16,rtmp); k17=Δt*rtmp
      for i in uidx
        update[i] = (b1*k1[i] + b2*k2[i] + b3*k3[i] + b5*k5[i]) + (b7*k7[i] + b9*k9[i] + b10*k10[i] + b11*k11[i]) + (b12*k12[i] + b13*k13[i] + b14*k14[i] + b15*k15[i]) + (b16*k16[i] + b17*k17[i])
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + update[i]
        end
        EEst = norm(((k2 - k16) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting, so in-place
        for i in uidx
          u[i] = u[i] + update[i]
        end
      end
      @ode_loopfooter
    end
  end
  if dense
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Feagin10,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16 = constructFeagin10(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u); k5 = similar(u)
  k6 = similar(u); k7 = similar(u); k8 = similar(u); k9 = similar(u); k10 = similar(u)
  k11 = similar(u); k12 = similar(u); k13 = similar(u); k14 = similar(u)
  k15 = similar(u); k16 = similar(u); k17 = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  utmp = similar(u); rtmp = rateType(sizeu)
  uidx = eachindex(u)
  k = rateType(sizeu)
  if dense
    pop!(ks)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k)
      for i in uidx
        k1[i]=Δt*k[i]
        tmp[i] = u[i] + a0100*k1[i]
      end
      f(t + c1*Δt,tmp,rtmp)
      for i in uidx
        k2[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0200*k1[i] + a0201*k2[i]
      end
      f(t + c2*Δt ,tmp,rtmp)
      for i in uidx
        k3[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0300*k1[i] + a0302*k3[i]
      end
      f(t + c3*Δt,tmp,rtmp)
      for i in uidx
        k4[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0400*k1[i] + a0402*k3[i] + a0403*k4[i]
      end
      f(t + c4*Δt,tmp,rtmp)
      for i in uidx
        k5[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0500*k1[i] + a0503*k4[i] + a0504*k5[i]
      end
      f(t + c5*Δt,tmp,rtmp)
      for i in uidx
        k6[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0600*k1[i] + a0603*k4[i] + a0604*k5[i] + a0605*k6[i]
      end
      f(t + c6*Δt,tmp,rtmp)
      for i in uidx
        k7[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0700*k1[i] + a0704*k5[i] + a0705*k6[i]) + a0706*k7[i]
      end
      f(t + c7*Δt,tmp,rtmp)
      for i in uidx
        k8[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0800*k1[i] + a0805*k6[i] + a0806*k7[i]) + a0807*k8[i]
      end
      f(t + c8*Δt,tmp,rtmp)
      for i in uidx
        k9[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0900*k1[i] + a0905*k6[i] + a0906*k7[i]) + a0907*k8[i] + a0908*k9[i]
      end
      f(t + c9*Δt,tmp,rtmp)
      for i in uidx
        k10[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1000*k1[i] + a1005*k6[i] + a1006*k7[i]) + a1007*k8[i] + a1008*k9[i] + a1009*k10[i]
      end
      f(t + c10*Δt,tmp,rtmp)
      for i in uidx
        k11[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1100*k1[i] + a1105*k6[i] + a1106*k7[i]) + (a1107*k8[i] + a1108*k9[i] + a1109*k10[i] + a1110*k11[i])
      end
      f(t + c11*Δt,tmp,rtmp)
      for i in uidx
        k12[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1200*k1[i] + a1203*k4[i] + a1204*k5[i]) + (a1205*k6[i] + a1206*k7[i] + a1207*k8[i] + a1208*k9[i]) + (a1209*k10[i] + a1210*k11[i] + a1211*k12[i])
      end
      f(t + c12*Δt,tmp,rtmp)
      for i in uidx
        k13[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1300*k1[i] + a1302*k3[i] + a1303*k4[i]) + (a1305*k6[i] + a1306*k7[i] + a1307*k8[i] + a1308*k9[i]) + (a1309*k10[i] + a1310*k11[i] + a1311*k12[i] + a1312*k13[i])
      end
      f(t + c13*Δt,tmp,rtmp)
      for i in uidx
        k14[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1400*k1[i] + a1401*k2[i] + a1404*k5[i]) + (a1406*k7[i] + a1412*k13[i] + a1413*k14[i])
      end
      f(t + c14*Δt,tmp,rtmp)
      for i in uidx
        k15[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a1500*k1[i] + a1502*k3[i] + a1514*k15[i]
      end
      f(t + c15*Δt,tmp,rtmp)
      for i in uidx
        k16[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1600*k1[i] + a1601*k2[i] + a1602*k3[i]) + (a1604*k5[i] + a1605*k6[i] + a1606*k7[i] + a1607*k8[i]) + (a1608*k9[i] + a1609*k10[i] + a1610*k11[i] + a1611*k12[i]) + (a1612*k13[i] + a1613*k14[i] + a1614*k15[i] + a1615*k16[i])
      end
      f(t + c16*Δt,tmp,rtmp)
      for i in uidx
        k17[i]=Δt*rtmp[i]
        tmp[i] = (b1*k1[i] + b2*k2[i] + b3*k3[i] + b5*k5[i]) + (b7*k7[i] + b9*k9[i] + b10*k10[i] + b11*k11[i]) + (b12*k12[i] + b13*k13[i] + b14*k14[i] + b15*k15[i]) + (b16*k16[i] + b17*k17[i])
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + tmp[i]
          atmp[i] = ((k2[i] - k16[i]) * adaptiveConst)./(abstol+u[i]*reltol)
        end
        EEst = norm(atmp,internalnorm)
      else #no chance of rejecting, so in-place
        for i in uidx
          u[i] = u[i] + tmp[i]
        end
      end
      @ode_loopfooter
    end
  end
  if dense
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:Feagin10,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16 = constructFeagin10(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType; local k5::uType
  local k6::uType; local k7::uType; local k8::uType; local k9::uType; local k10::uType
  local k11::uType; local k12::uType; local k13::uType; local k14::uType
  local k15::uType; local k16::uType; local k17::uType
  if dense
    pop!(ks)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k   = f(t,u)
      k1  = Δt*k
      k2  = Δt*f(t + c1*Δt,u + a0100*k1)
      k3  = Δt*f(t + c2*Δt ,u + a0200*k1 + a0201*k2)
      k4  = Δt*f(t + c3*Δt,u + a0300*k1              + a0302*k3)
      k5  = Δt*f(t + c4*Δt,u + a0400*k1              + a0402*k3 + a0403*k4)
      k6  = Δt*f(t + c5*Δt,u + a0500*k1                           + a0503*k4 + a0504*k5)
      k7  = Δt*f(t + c6*Δt,u + a0600*k1                           + a0603*k4 + a0604*k5 + a0605*k6)
      k8  = Δt*f(t + c7*Δt,u + a0700*k1                                        + a0704*k5 + a0705*k6 + a0706*k7)
      k9  = Δt*f(t + c8*Δt,u + a0800*k1                                                     + a0805*k6 + a0806*k7 + a0807*k8)
      k10 = Δt*f(t + c9*Δt,u + a0900*k1                                                     + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9)
      k11 = Δt*f(t + c10*Δt,u + a1000*k1                                                     + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10)
      k12 = Δt*f(t + c11*Δt,u + a1100*k1                                                     + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11)
      k13 = Δt*f(t + c12*Δt,u + a1200*k1                           + a1203*k4 + a1204*k5 + a1205*k6 + a1206*k7 + a1207*k8 + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12)
      k14 = Δt*f(t + c13*Δt,u + a1300*k1              + a1302*k3 + a1303*k4              + a1305*k6 + a1306*k7 + a1307*k8 + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13)
      k15 = Δt*f(t + c14*Δt,u + a1400*k1 + a1401*k2                           + a1404*k5              + a1406*k7 +                                                                     a1412*k13 + a1413*k14)
      k16 = Δt*f(t + c15*Δt,u + a1500*k1              + a1502*k3                                                                                                                                                     + a1514*k15)
      k17 = Δt*f(t + c16*Δt,u + a1600*k1 + a1601*k2 + a1602*k3              + a1604*k5 + a1605*k6 + a1606*k7 + a1607*k8 + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16)
      update = (b1*k1 + b2*k2 + b3*k3 + b5*k5) + (b7*k7 + b9*k9 + b10*k10 + b11*k11) + (b12*k12 + b13*k13 + b14*k14 + b15*k15) + (b16*k16 + b17*k17)
      if adaptive
        utmp = u + update
        EEst = norm(((k2 - k16) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else
        u = u + update
      end
      @ode_numberloopfooter
    end
  end
  if dense
    k = f(t,u)
    push!(ks,k)
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Feagin12Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25 = constructFeagin12(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u); k5 = similar(u)
  k6 = similar(u); k7 = similar(u); k8 = similar(u); k9 = similar(u); k10 = similar(u)
  k11 = similar(u); k12 = similar(u); k13 = similar(u); k14 = similar(u)
  k15 = similar(u); k16 = similar(u); k17 = similar(u); k18 = similar(u)
  k19 = similar(u); k20 = similar(u); k21 = similar(u); k22 = similar(u)
  k23 = similar(u); k24 = similar(u); k25 = similar(u)
  update = similar(u); rtmp = rateType(sizeu)
  utmp = similar(u)
  uidx = eachindex(u)
  k = rateType(sizeu)
  if dense
    pop!(ks)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k); k1=Δt*k
      f(t + c1*Δt,u + a0100*k1,rtmp); k2=Δt*rtmp
      f(t + c2*Δt ,u + a0200*k1 + a0201*k2,rtmp); k3=Δt*rtmp
      f(t + c3*Δt,u + a0300*k1              + a0302*k3,rtmp); k4=Δt*rtmp
      f(t + c4*Δt,u + a0400*k1              + a0402*k3 + a0403*k4,rtmp); k5=Δt*rtmp
      f(t + c5*Δt,u + a0500*k1                           + a0503*k4 + a0504*k5,rtmp); k6=Δt*rtmp
      f(t + c6*Δt,u + a0600*k1                           + a0603*k4 + a0604*k5 + a0605*k6,rtmp); k7=Δt*rtmp
      f(t + c7*Δt,u + a0700*k1                                        + a0704*k5 + a0705*k6 + a0706*k7,rtmp); k8=Δt*rtmp
      f(t + c8*Δt,u + a0800*k1                                                     + a0805*k6 + a0806*k7 + a0807*k8,rtmp); k9=Δt*rtmp
      f(t + c9*Δt,u + a0900*k1                                                     + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9,rtmp); k10=Δt*rtmp
      f(t + c10*Δt,u + a1000*k1                                                     + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10,rtmp); k11=Δt*rtmp
      f(t + c11*Δt,u + a1100*k1                                                     + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11,rtmp); k12=Δt*rtmp
      f(t + c12*Δt,u + a1200*k1                                                                                            + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12,rtmp); k13=Δt*rtmp
      f(t + c13*Δt,u + a1300*k1                                                                                            + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13,rtmp); k14=Δt*rtmp
      f(t + c14*Δt,u + a1400*k1                                                                                            + a1408*k9 + a1409*k10 + a1410*k11 + a1411*k12 + a1412*k13 + a1413*k14,rtmp); k15=Δt*rtmp
      f(t + c15*Δt,u + a1500*k1                                                                                            + a1508*k9 + a1509*k10 + a1510*k11 + a1511*k12 + a1512*k13 + a1513*k14 + a1514*k15,rtmp); k16=Δt*rtmp
      f(t + c16*Δt,u + a1600*k1                                                                                            + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16,rtmp); k17=Δt*rtmp
      f(t + c17*Δt,u + a1700*k1                                                     + a1705*k6 + a1706*k7 + a1707*k8 + a1708*k9 + a1709*k10 + a1710*k11 + a1711*k12 + a1712*k13 + a1713*k14 + a1714*k15 + a1715*k16 + a1716*k17,rtmp); k18=Δt*rtmp
      f(t + c18*Δt,u + a1800*k1                                                     + a1805*k6 + a1806*k7 + a1807*k8 + a1808*k9 + a1809*k10 + a1810*k11 + a1811*k12 + a1812*k13 + a1813*k14 + a1814*k15 + a1815*k16 + a1816*k17 + a1817*k18,rtmp); k19=Δt*rtmp
      f(t + c19*Δt,u + a1900*k1                                        + a1904*k5 + a1905*k6 + a1906*k7              + a1908*k9 + a1909*k10 + a1910*k11 + a1911*k12 + a1912*k13 + a1913*k14 + a1914*k15 + a1915*k16 + a1916*k17 + a1917*k18 + a1918*k19,rtmp); k20=Δt*rtmp
      f(t + c20*Δt,u + a2000*k1                           + a2003*k4 + a2004*k5 + a2005*k6              + a2007*k8              + a2009*k10 + a2010*k11                                                                                     + a2017*k18 + a2018*k19 + a2019*k20,rtmp); k21=Δt*rtmp
      f(t + c21*Δt,u + a2100*k1              + a2102*k3 + a2103*k4                           + a2106*k7 + a2107*k8              + a2109*k10 + a2110*k11                                                                                     + a2117*k18 + a2118*k19 + a2119*k20 + a2120*k21,rtmp); k22=Δt*rtmp
      f(t + c22*Δt,u + a2200*k1 + a2201*k2                           + a2204*k5              + a2206*k7                                                                                                                                                                                     + a2220*k21 + a2221*k22,rtmp); k23=Δt*rtmp
      f(t + c23*Δt,u + a2300*k1              + a2302*k3                                                                                                                                                                                                                                                                     + a2322*k23,rtmp); k24=Δt*rtmp
      f(t + c24*Δt,u + a2400*k1 + a2401*k2 + a2402*k3              + a2404*k5              + a2406*k7 + a2407*k8 + a2408*k9 + a2409*k10 + a2410*k11 + a2411*k12 + a2412*k13 + a2413*k14 + a2414*k15 + a2415*k16 + a2416*k17 + a2417*k18 + a2418*k19 + a2419*k20 + a2420*k21 + a2421*k22 + a2422*k23 + a2423*k24,rtmp); k25=Δt*rtmp

      for i in uidx
        update[i] = (b1*k1[i] + b2*k2[i] + b3*k3[i] + b5*k5[i]) + (b7*k7[i] + b8*k8[i] + b10*k10[i] + b11*k11[i]) + (b13*k13[i] + b14*k14[i] + b15*k15[i] + b16*k16[i]) + (b17*k17[i] + b18*k18[i] + b19*k19[i] + b20*k20[i]) + (b21*k21[i] + b22*k22[i] + b23*k23[i] + b24*k24[i]) + b25*k25[i]
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + update[i]
        end
        EEst = norm(((k2 - k24) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting so in-place
        for i in uidx
          u[i] = u[i] + update[i]
        end
      end
      @ode_loopfooter
    end
  end
  if dense
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Feagin12,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25 = constructFeagin12(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u); k5 = similar(u)
  k6 = similar(u); k7 = similar(u); k8 = similar(u); k9 = similar(u); k10 = similar(u)
  k11 = similar(u); k12 = similar(u); k13 = similar(u); k14 = similar(u)
  k15 = similar(u); k16 = similar(u); k17 = similar(u); k18 = similar(u)
  k19 = similar(u); k20 = similar(u); k21 = similar(u); k22 = similar(u)
  k23 = similar(u); k24 = similar(u); k25 = similar(u)
  update = similar(u)
  utmp = similar(u); rtmp = rateType(sizeu)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  uidx = eachindex(u)
  k = rateType(sizeu)
  if dense
    pop!(ks)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k)
      for i in uidx
        k1[i]=Δt*k[i]
        tmp[i] = u[i] + a0100*k1[i]
      end
      f(t + c1*Δt,tmp,rtmp)
      for i in uidx
        k2[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0200*k1[i] + a0201*k2[i]
      end
      f(t + c2*Δt ,tmp,rtmp)
      for i in uidx
        k3[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0300*k1[i] + a0302*k3[i]
      end
      f(t + c3*Δt,tmp,rtmp)
      for i in uidx
        k4[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0400*k1[i] + a0402*k3[i] + a0403*k4[i]
      end
      f(t + c4*Δt,tmp,rtmp)
      for i in uidx
        k5[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0500*k1[i] + a0503*k4[i] + a0504*k5[i]
      end
      f(t + c5*Δt,tmp,rtmp)
      for i in uidx
        k6[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0600*k1[i] + a0603*k4[i] + a0604*k5[i]) + a0605*k6[i]
      end
      f(t + c6*Δt,tmp,rtmp)
      for i in uidx
        k7[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0700*k1[i] + a0704*k5[i] + a0705*k6[i]) + a0706*k7[i]
      end
      f(t + c7*Δt,tmp,rtmp)
      for i in uidx
        k8[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0800*k1[i] + a0805*k6[i] + a0806*k7[i]) + a0807*k8[i]
      end
      f(t + c8*Δt,tmp,rtmp)
      for i in uidx
        k9[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0900*k1[i] + a0905*k6[i] + a0906*k7[i]) + (a0907*k8[i] + a0908*k9[i])
      end
      f(t + c9*Δt,tmp,rtmp)
      for i in uidx
        k10[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1000*k1[i] + a1005*k6[i] + a1006*k7[i]) + (a1007*k8[i] + a1008*k9[i] + a1009*k10[i])
      end
      f(t + c10*Δt,tmp,rtmp)
      for i in uidx
        k11[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1100*k1[i] + a1105*k6[i] + a1106*k7[i]) + (a1107*k8[i] + a1108*k9[i] + a1109*k10[i] + a1110*k11[i])
      end
      f(t + c11*Δt,tmp,rtmp)
      for i in uidx
        k12[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1200*k1[i] + a1208*k9[i] + a1209*k10[i]) + (a1210*k11[i] + a1211*k12[i])
      end
      f(t + c12*Δt,tmp,rtmp)
      for i in uidx
        k13[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1300*k1[i] + a1308*k9[i] + a1309*k10[i]) + (a1310*k11[i] + a1311*k12[i] + a1312*k13[i])
      end
      f(t + c13*Δt,tmp,rtmp)
      for i in uidx
        k14[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1400*k1[i] + a1408*k9[i] + a1409*k10[i]) + (a1410*k11[i] + a1411*k12[i] + a1412*k13[i] + a1413*k14[i])
      end
      f(t + c14*Δt,tmp,rtmp)
      for i in uidx
        k15[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1500*k1[i] + a1508*k9[i] + a1509*k10[i]) + (a1510*k11[i] + a1511*k12[i] + a1512*k13[i] + a1513*k14[i]) + a1514*k15[i]
      end
      f(t + c15*Δt,tmp,rtmp)
      for i in uidx
        k16[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1600*k1[i] + a1608*k9[i] + a1609*k10[i]) + (a1610*k11[i] + a1611*k12[i] + a1612*k13[i] + a1613*k14[i]) + (a1614*k15[i] + a1615*k16[i])
      end
      f(t + c16*Δt,tmp,rtmp)
      for i in uidx
        k17[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1700*k1[i] + a1705*k6[i] + a1706*k7[i]) + (a1707*k8[i] + a1708*k9[i] + a1709*k10[i] + a1710*k11[i]) + (a1711*k12[i] + a1712*k13[i] + a1713*k14[i] + a1714*k15[i]) + (a1715*k16[i] + a1716*k17[i])
      end
      f(t + c17*Δt,tmp,rtmp)
      for i in uidx
        k18[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1800*k1[i] + a1805*k6[i] + a1806*k7[i]) + (a1807*k8[i] + a1808*k9[i] + a1809*k10[i] + a1810*k11[i]) + (a1811*k12[i] + a1812*k13[i] + a1813*k14[i] + a1814*k15[i]) + (a1815*k16[i] + a1816*k17[i] + a1817*k18[i])
      end
      f(t + c18*Δt,tmp,rtmp)
      for i in uidx
        k19[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1900*k1[i] + a1904*k5[i] + a1905*k6[i]) + (a1906*k7[i] + a1908*k9[i] + a1909*k10[i] + a1910*k11[i]) + (a1911*k12[i] + a1912*k13[i] + a1913*k14[i] + a1914*k15[i]) + (a1915*k16[i] + a1916*k17[i] + a1917*k18[i] + a1918*k19[i])
      end
      f(t + c19*Δt,tmp,rtmp)
      for i in uidx
        k20[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2000*k1[i] + a2003*k4[i] + a2004*k5[i]) + (a2005*k6[i] + a2007*k8[i] + a2009*k10[i] + a2010*k11[i]) + (a2017*k18[i] + a2018*k19[i] + a2019*k20[i])
      end
      f(t + c20*Δt,tmp,rtmp)
      for i in uidx
        k21[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2100*k1[i] + a2102*k3[i] + a2103*k4[i]) + (a2106*k7[i] + a2107*k8[i] + a2109*k10[i] + a2110*k11[i]) + (a2117*k18[i] + a2118*k19[i] + a2119*k20[i] + a2120*k21[i])
      end
      f(t + c21*Δt,tmp,rtmp)
      for i in uidx
        k22[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2200*k1[i] + a2201*k2[i] + a2204*k5[i]) + (a2206*k7[i] + a2220*k21[i] + a2221*k22[i])
      end
      f(t + c22*Δt,tmp,rtmp)
      for i in uidx
        k23[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a2300*k1[i] + a2302*k3[i] + a2322*k23[i]
      end
      f(t + c23*Δt,tmp,rtmp)
      for i in uidx
        k24[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2400*k1[i] + a2401*k2[i] + a2402*k3[i]) + (a2404*k5[i] + a2406*k7[i] + a2407*k8[i] + a2408*k9[i]) + (a2409*k10[i] + a2410*k11[i] + a2411*k12[i] + a2412*k13[i]) + (a2413*k14[i] + a2414*k15[i] + a2415*k16[i] + a2416*k17[i]) + (a2417*k18[i] + a2418*k19[i] + a2419*k20[i] + a2420*k21[i]) + (a2421*k22[i] + a2422*k23[i] + a2423*k24[i])
      end
      f(t + c24*Δt,tmp,rtmp)
      for i in uidx
        k25[i] = Δt*rtmp[i]
        update[i] = (b1*k1[i] + b2*k2[i] + b3*k3[i] + b5*k5[i]) + (b7*k7[i] + b8*k8[i] + b10*k10[i] + b11*k11[i]) + (b13*k13[i] + b14*k14[i] + b15*k15[i] + b16*k16[i]) + (b17*k17[i] + b18*k18[i] + b19*k19[i] + b20*k20[i]) + (b21*k21[i] + b22*k22[i] + b23*k23[i] + b24*k24[i]) + b25*k25[i]
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + update[i]
          atmp[i] = ((k2[i] - k24[i]) * adaptiveConst)/(abstol+u[i]*reltol)
        end
        EEst = norm(atmp,internalnorm)
      else #no chance of rejecting so in-place
        for i in uidx
          u[i] = u[i] + update[i]
        end
      end
      @ode_loopfooter
    end
  end
  if dense
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:Feagin12,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25 = constructFeagin12(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType; local k5::uType
  local k6::uType; local k7::uType; local k8::uType; local k9::uType; local k10::uType
  local k11::uType; local k12::uType; local k13::uType; local k14::uType
  local k15::uType; local k16::uType; local k17::uType; local k18::uType
  local k19::uType; local k20::uType; local k21::uType; local k22::uType
  local k23::uType; local k24::uType; local k25::uType
  if dense
    pop!(ks)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k   = f(t,u)
      k1  = Δt*k
      k2  = Δt*f(t + c1*Δt,u + a0100*k1)
      k3  = Δt*f(t + c2*Δt ,u + a0200*k1 + a0201*k2)
      k4  = Δt*f(t + c3*Δt,u + a0300*k1              + a0302*k3)
      k5  = Δt*f(t + c4*Δt,u + a0400*k1              + a0402*k3 + a0403*k4)
      k6  = Δt*f(t + c5*Δt,u + a0500*k1                           + a0503*k4 + a0504*k5)
      k7  = Δt*f(t + c6*Δt,u + a0600*k1                           + a0603*k4 + a0604*k5 + a0605*k6)
      k8  = Δt*f(t + c7*Δt,u + a0700*k1                                        + a0704*k5 + a0705*k6 + a0706*k7)
      k9  = Δt*f(t + c8*Δt,u + a0800*k1                                                     + a0805*k6 + a0806*k7 + a0807*k8)
      k10 = Δt*f(t + c9*Δt,u + a0900*k1                                                     + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9)
      k11 = Δt*f(t + c10*Δt,u + a1000*k1                                                     + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10)
      k12 = Δt*f(t + c11*Δt,u + a1100*k1                                                     + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11)
      k13 = Δt*f(t + c12*Δt,u + a1200*k1                                                                                            + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12)
      k14 = Δt*f(t + c13*Δt,u + a1300*k1                                                                                            + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13)
      k15 = Δt*f(t + c14*Δt,u + a1400*k1                                                                                            + a1408*k9 + a1409*k10 + a1410*k11 + a1411*k12 + a1412*k13 + a1413*k14)
      k16 = Δt*f(t + c15*Δt,u + a1500*k1                                                                                            + a1508*k9 + a1509*k10 + a1510*k11 + a1511*k12 + a1512*k13 + a1513*k14 + a1514*k15)
      k17 = Δt*f(t + c16*Δt,u + a1600*k1                                                                                            + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16)
      k18 = Δt*f(t + c17*Δt,u + a1700*k1                                                     + a1705*k6 + a1706*k7 + a1707*k8 + a1708*k9 + a1709*k10 + a1710*k11 + a1711*k12 + a1712*k13 + a1713*k14 + a1714*k15 + a1715*k16 + a1716*k17)
      k19 = Δt*f(t + c18*Δt,u + a1800*k1                                                     + a1805*k6 + a1806*k7 + a1807*k8 + a1808*k9 + a1809*k10 + a1810*k11 + a1811*k12 + a1812*k13 + a1813*k14 + a1814*k15 + a1815*k16 + a1816*k17 + a1817*k18)
      k20 = Δt*f(t + c19*Δt,u + a1900*k1                                        + a1904*k5 + a1905*k6 + a1906*k7              + a1908*k9 + a1909*k10 + a1910*k11 + a1911*k12 + a1912*k13 + a1913*k14 + a1914*k15 + a1915*k16 + a1916*k17 + a1917*k18 + a1918*k19)
      k21 = Δt*f(t + c20*Δt,u + a2000*k1                           + a2003*k4 + a2004*k5 + a2005*k6              + a2007*k8              + a2009*k10 + a2010*k11                                                                                     + a2017*k18 + a2018*k19 + a2019*k20)
      k22 = Δt*f(t + c21*Δt,u + a2100*k1              + a2102*k3 + a2103*k4                           + a2106*k7 + a2107*k8              + a2109*k10 + a2110*k11                                                                                     + a2117*k18 + a2118*k19 + a2119*k20 + a2120*k21)
      k23 = Δt*f(t + c22*Δt,u + a2200*k1 + a2201*k2                           + a2204*k5              + a2206*k7                                                                                                                                                                                     + a2220*k21 + a2221*k22)
      k24 = Δt*f(t + c23*Δt,u + a2300*k1              + a2302*k3                                                                                                                                                                                                                                                                     + a2322*k23)
      k25 = Δt*f(t + c24*Δt,u + a2400*k1 + a2401*k2 + a2402*k3              + a2404*k5              + a2406*k7 + a2407*k8 + a2408*k9 + a2409*k10 + a2410*k11 + a2411*k12 + a2412*k13 + a2413*k14 + a2414*k15 + a2415*k16 + a2416*k17 + a2417*k18 + a2418*k19 + a2419*k20 + a2420*k21 + a2421*k22 + a2422*k23 + a2423*k24)

      update = (b1*k1 + b2*k2 + b3*k3 + b5*k5) + (b7*k7 + b8*k8 + b10*k10 + b11*k11) + (b13*k13 + b14*k14 + b15*k15 + b16*k16) + (b17*k17 + b18*k18 + b19*k19 + b20*k20) + (b21*k21 + b22*k22 + b23*k23 + b24*k24) + (b25*k25)
      if adaptive
        utmp = u + update
        EEst = norm(((k2 - k24) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting so in-place
        u = u + update
      end
      @ode_numberloopfooter
    end
  end
  if dense
    k = f(t,u)
    push!(ks,k)
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Feagin14,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,b35 = constructFeagin14(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u); k5 = similar(u)
  k6 = similar(u); k7 = similar(u); k8 = similar(u); k9 = similar(u); k10 = similar(u)
  k11 = similar(u); k12 = similar(u); k13 = similar(u); k14 = similar(u)
  k15 = similar(u); k16 = similar(u); k17 = similar(u); k18 = similar(u)
  k19 = similar(u); k20 = similar(u); k21 = similar(u); k22 = similar(u)
  k23 = similar(u); k24 = similar(u); k25 = similar(u)
  k26 = similar(u); k27 = similar(u); k28 = similar(u)
  k29 = similar(u); k30 = similar(u); k31 = similar(u); k32 = similar(u)
  k33 = similar(u); k34 = similar(u); k35 = similar(u)
  update = similar(u)
  utmp = similar(u);
  rtmp = rateType(sizeu)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  uidx = eachindex(u)
  k = rateType(sizeu)
  if dense
    pop!(ks)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k)
      for i in uidx
        k1[i]=Δt*k[i]
        tmp[i] = u[i] + a0100*k1[i]
      end
      f(t + c1*Δt,tmp,rtmp)
      for i in uidx
        k2[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0200*k1[i] + a0201*k2[i]
      end
      f(t + c2*Δt ,tmp,rtmp)
      for i in uidx
        k3[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0300*k1[i] + a0302*k3[i]
      end
      f(t + c3*Δt,tmp,rtmp)
      for i in uidx
        k4[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0400*k1[i] + a0402*k3[i] + a0403*k4[i]
      end
      f(t + c4*Δt,tmp,rtmp)
      for i in uidx
        k5[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0500*k1[i] + a0503*k4[i] + a0504*k5[i]
      end
      f(t + c5*Δt,tmp,rtmp)
      for i in uidx
        k6[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0600*k1[i] + a0603*k4[i] + a0604*k5[i]) + a0605*k6[i]
      end
      f(t + c6*Δt,tmp,rtmp)
      for i in uidx
        k7[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0700*k1[i] + a0704*k5[i] + a0705*k6[i]) + a0706*k7[i]
      end
      f(t + c7*Δt,tmp,rtmp)
      for i in uidx
        k8[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0800*k1[i] + a0805*k6[i] + a0806*k7[i]) + a0807*k8[i]
      end
      f(t + c8*Δt,tmp,rtmp)
      for i in uidx
        k9[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0900*k1[i] + a0905*k6[i] + a0906*k7[i]) + a0907*k8[i] + a0908*k9[i]
      end
      f(t + c9*Δt,tmp,rtmp)
      for i in uidx
        k10[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1000*k1[i] + a1005*k6[i] + a1006*k7[i]) + (a1007*k8[i] + a1008*k9[i] + a1009*k10[i])
      end
      f(t + c10*Δt,tmp,rtmp)
      for i in uidx
        k11[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1100*k1[i] + a1105*k6[i] + a1106*k7[i]) + (a1107*k8[i] + a1108*k9[i] + a1109*k10[i] + a1110*k11[i])
      end
      f(t + c11*Δt,tmp,rtmp)
      for i in uidx
        k12[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1200*k1[i] + a1208*k9[i] + a1209*k10[i]) + (a1210*k11[i] + a1211*k12[i])
      end
      f(t + c12*Δt,tmp,rtmp)
      for i in uidx
        k13[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1300*k1[i] + a1308*k9[i] + a1309*k10[i]) + (a1310*k11[i] + a1311*k12[i] + a1312*k13[i])
      end
      f(t + c13*Δt,tmp,rtmp)
      for i in uidx
        k14[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1400*k1[i] + a1408*k9[i] + a1409*k10[i]) + (a1410*k11[i] + a1411*k12[i] + a1412*k13[i] + a1413*k14[i])
      end
      f(t + c14*Δt,tmp,rtmp)
      for i in uidx
        k15[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1500*k1[i] + a1508*k9[i] + a1509*k10[i]) + (a1510*k11[i] + a1511*k12[i] + a1512*k13[i] + a1513*k14[i]) + a1514*k15[i]
      end
      f(t + c15*Δt,tmp,rtmp)
      for i in uidx
        k16[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1600*k1[i] + a1608*k9[i] + a1609*k10[i]) + (a1610*k11[i] + a1611*k12[i] + a1612*k13[i] + a1613*k14[i]) + a1614*k15[i] + a1615*k16[i]
      end
      f(t + c16*Δt,tmp,rtmp)
      for i in uidx
        k17[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1700*k1[i] + a1712*k13[i] + a1713*k14[i]) + (a1714*k15[i] + a1715*k16[i] + a1716*k17[i])
      end
      f(t + c17*Δt,tmp,rtmp)
      for i in uidx
        k18[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1800*k1[i] + a1812*k13[i] + a1813*k14[i]) + (a1814*k15[i] + a1815*k16[i] + a1816*k17[i] + a1817*k18[i])
      end
      f(t + c18*Δt,tmp,rtmp)
      for i in uidx
        k19[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1900*k1[i] + a1912*k13[i] + a1913*k14[i]) + (a1914*k15[i] + a1915*k16[i] + a1916*k17[i] + a1917*k18[i]) + a1918*k19[i]
      end
      f(t + c19*Δt,tmp,rtmp)
      for i in uidx
        k20[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2000*k1[i] + a2012*k13[i] + a2013*k14[i]) + (a2014*k15[i] + a2015*k16[i] + a2016*k17[i] + a2017*k18[i]) + (a2018*k19[i] + a2019*k20[i])
      end
      f(t + c20*Δt,tmp,rtmp)
      for i in uidx
        k21[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2100*k1[i] + a2112*k13[i] + a2113*k14[i]) + (a2114*k15[i] + a2115*k16[i] + a2116*k17[i] + a2117*k18[i]) + (a2118*k19[i] + a2119*k20[i] + a2120*k21[i])
      end
      f(t + c21*Δt,tmp,rtmp)
      for i in uidx
        k22[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2200*k1[i] + a2212*k13[i] + a2213*k14[i]) + (a2214*k15[i] + a2215*k16[i] + a2216*k17[i] + a2217*k18[i]) + (a2218*k19[i] + a2219*k20[i] + a2220*k21[i] + a2221*k22[i])
      end
      f(t + c22*Δt,tmp,rtmp)
      for i in uidx
        k23[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2300*k1[i] + a2308*k9[i] + a2309*k10[i]) + (a2310*k11[i] + a2311*k12[i] + a2312*k13[i] + a2313*k14[i]) + (a2314*k15[i] + a2315*k16[i] + a2316*k17[i] + a2317*k18[i]) + (a2318*k19[i] + a2319*k20[i] + a2320*k21[i] + a2321*k22[i]) + (a2322*k23[i])
      end
      f(t + c23*Δt,tmp,rtmp)
      for i in uidx
        k24[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2400*k1[i] + a2408*k9[i] + a2409*k10[i]) + (a2410*k11[i] + a2411*k12[i] + a2412*k13[i] + a2413*k14[i]) + (a2414*k15[i] + a2415*k16[i] + a2416*k17[i] + a2417*k18[i]) + (a2418*k19[i] + a2419*k20[i] + a2420*k21[i] + a2421*k22[i]) + (a2422*k23[i] + a2423*k24[i])
      end
      f(t + c24*Δt,tmp,rtmp)
      for i in uidx
        k25[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2500*k1[i] + a2508*k9[i] + a2509*k10[i]) + (a2510*k11[i] + a2511*k12[i] + a2512*k13[i] + a2513*k14[i]) + (a2514*k15[i] + a2515*k16[i] + a2516*k17[i] + a2517*k18[i]) + (a2518*k19[i] + a2519*k20[i] + a2520*k21[i] + a2521*k22[i]) + (a2522*k23[i] + a2523*k24[i] + a2524*k25[i])
      end
      f(t + c25*Δt,tmp,rtmp)
      for i in uidx
        k26[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2600*k1[i] + a2605*k6[i] + a2606*k7[i]) + (a2607*k8[i] + a2608*k9[i] + a2609*k10[i] + a2610*k11[i]) + (a2612*k13[i] + a2613*k14[i] + a2614*k15[i] + a2615*k16[i]) + (a2616*k17[i] + a2617*k18[i] + a2618*k19[i] + a2619*k20[i]) + (a2620*k21[i] + a2621*k22[i] + a2622*k23[i] + a2623*k24[i]) + (a2624*k25[i] + a2625*k26[i])
      end
      f(t + c26*Δt,tmp,rtmp)
      for i in uidx
        k27[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2700*k1[i] + a2705*k6[i] + a2706*k7[i]) + (a2707*k8[i] + a2708*k9[i] + a2709*k10[i] + a2711*k12[i]) + (a2712*k13[i] + a2713*k14[i] + a2714*k15[i] + a2715*k16[i]) + (a2716*k17[i] + a2717*k18[i] + a2718*k19[i] + a2719*k20[i]) + (a2720*k21[i] + a2721*k22[i] + a2722*k23[i] + a2723*k24[i]) + (a2724*k25[i] + a2725*k26[i] + a2726*k27[i])
      end
      f(t + c27*Δt,tmp,rtmp)
      for i in uidx
        k28[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2800*k1[i] + a2805*k6[i] + a2806*k7[i]) + (a2807*k8[i] + a2808*k9[i] + a2810*k11[i] + a2811*k12[i]) + (a2813*k14[i] + a2814*k15[i] + a2815*k16[i] + a2823*k24[i]) + (a2824*k25[i] + a2825*k26[i] + a2826*k27[i] + a2827*k28[i])
      end
      f(t + c28*Δt,tmp,rtmp)
      for i in uidx
        k29[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2900*k1[i] + a2904*k5[i] + a2905*k6[i]) + (a2906*k7[i] + a2909*k10[i] + a2910*k11[i] + a2911*k12[i]) + (a2913*k14[i] + a2914*k15[i] + a2915*k16[i] + a2923*k24[i]) + (a2924*k25[i] + a2925*k26[i] + a2926*k27[i] + a2927*k28[i]) + (a2928*k29[i])
      end
      f(t + c29*Δt,tmp,rtmp)
      for i in uidx
        k30[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a3000*k1[i] + a3003*k4[i] + a3004*k5[i]) + (a3005*k6[i] + a3007*k8[i] + a3009*k10[i] + a3010*k11[i]) + (a3013*k14[i] + a3014*k15[i] + a3015*k16[i] + a3023*k24[i]) + (a3024*k25[i] + a3025*k26[i] + a3027*k28[i] + a3028*k29[i]) + (a3029*k30[i])
      end
      f(t + c30*Δt,tmp,rtmp)
      for i in uidx
        k31[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a3100*k1[i] + a3102*k3[i] + a3103*k4[i]) + (a3106*k7[i] + a3107*k8[i] + a3109*k10[i] + a3110*k11[i]) + (a3113*k14[i] + a3114*k15[i] + a3115*k16[i] + a3123*k24[i]) + (a3124*k25[i] + a3125*k26[i] + a3127*k28[i] + a3128*k29[i]) + (a3129*k30[i] + a3130*k31[i])
      end
      f(t + c31*Δt,tmp,rtmp)
      for i in uidx
        k32[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a3200*k1[i] + a3201*k2[i] + a3204*k5[i]) + (a3206*k7[i] + a3230*k31[i] + a3231*k32[i])
      end
      f(t + c32*Δt,tmp,rtmp)
      for i in uidx
        k33[i]=Δt*rtmp[i]
        tmp[i] = u[i] + a3300*k1[i] + a3302*k3[i] + a3332*k33[i]
      end
      f(t + c33*Δt,tmp,rtmp)
      for i in uidx
        k34[i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a3400*k1[i] + a3401*k2[i] + a3402*k3[i]) + (a3404*k5[i] + a3406*k7[i] + a3407*k8[i] + a3409*k10[i]) + (a3410*k11[i] + a3411*k12[i] + a3412*k13[i] + a3413*k14[i]) + (a3414*k15[i] + a3415*k16[i] + a3416*k17[i] + a3417*k18[i]) + (a3418*k19[i] + a3419*k20[i] + a3420*k21[i] + a3421*k22[i]) + (a3422*k23[i] + a3423*k24[i] + a3424*k25[i] + a3425*k26[i]) + (a3426*k27[i] + a3427*k28[i] + a3428*k29[i] + a3429*k30[i]) + (a3430*k31[i] + a3431*k32[i] + a3432*k33[i] + a3433*k34[i])
      end
      f(t + c34*Δt,tmp,rtmp)
      for i in uidx
        k35[i]=Δt*rtmp[i]
        update[i] = (b1*k1[i] + b2*k2[i] + b3*k3[i] + b5*k5[i]) + (b7*k7[i] + b8*k8[i] + b10*k10[i] + b11*k11[i]) + (b12*k12[i] + b14*k14[i] + b15*k15[i] + b16*k16[i]) + (b18*k18[i] + b19*k19[i] + b20*k20[i] + b21*k21[i]) + (b22*k22[i] + b23*k23[i] + b24*k24[i] + b25*k25[i]) + (b26*k26[i] + b27*k27[i] + b28*k28[i] + b29*k29[i]) + (b30*k30[i] + b31*k31[i] + b32*k32[i] + b33*k33[i]) + (b34*k34[i] + b35*k35[i])
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + update[i]
          atmp[i] = ((k2[i] - k34[i]) * adaptiveConst)./(abstol+u[i]*reltol)
        end
        EEst = norm(atmp,internalnorm)
      else #no chance of rejecting, so in-place
        for i in uidx
          u[i] = u[i] + update[i]
        end
      end
      @ode_loopfooter
    end
  end
  if dense
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Feagin14Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,b35 = constructFeagin14(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u); k5 = similar(u)
  k6 = similar(u); k7 = similar(u); k8 = similar(u); k9 = similar(u); k10 = similar(u)
  k11 = similar(u); k12 = similar(u); k13 = similar(u); k14 = similar(u)
  k15 = similar(u); k16 = similar(u); k17 = similar(u); k18 = similar(u)
  k19 = similar(u); k20 = similar(u); k21 = similar(u); k22 = similar(u)
  k23 = similar(u); k24 = similar(u); k25 = similar(u)
  k26 = similar(u); k27 = similar(u); k28 = similar(u)
  k29 = similar(u); k30 = similar(u); k31 = similar(u); k32 = similar(u)
  k33 = similar(u); k34 = similar(u); k35 = similar(u)
  update = similar(u)
  utmp = similar(u)
  uidx = eachindex(u)
  rtmp = rateType(sizeu)
  k = rateType(sizeu)
  if dense
    pop!(ks)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,k); k1=Δt*k
      f(t + c1*Δt,u + a0100*k1,rtmp); k2=Δt*rtmp
      f(t + c2*Δt ,u + a0200*k1 + a0201*k2,rtmp); k3=Δt*rtmp
      f(t + c3*Δt,u + a0300*k1              + a0302*k3,rtmp); k4=Δt*rtmp
      f(t + c4*Δt,u + a0400*k1              + a0402*k3 + a0403*k4,rtmp); k5=Δt*rtmp
      f(t + c5*Δt,u + a0500*k1                           + a0503*k4 + a0504*k5,rtmp); k6=Δt*rtmp
      f(t + c6*Δt,u + a0600*k1                           + a0603*k4 + a0604*k5 + a0605*k6,rtmp); k7=Δt*rtmp
      f(t + c7*Δt,u + a0700*k1                                        + a0704*k5 + a0705*k6 + a0706*k7,rtmp); k8=Δt*rtmp
      f(t + c8*Δt,u + a0800*k1                                                     + a0805*k6 + a0806*k7 + a0807*k8,rtmp); k9=Δt*rtmp
      f(t + c9*Δt,u + a0900*k1                                                     + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9,rtmp); k10=Δt*rtmp
      f(t + c10*Δt,u + a1000*k1                                                     + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10,rtmp); k11=Δt*rtmp
      f(t + c11*Δt,u + a1100*k1                                                     + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11,rtmp); k12=Δt*rtmp
      f(t + c12*Δt,u + a1200*k1                                                                                            + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12,rtmp); k13=Δt*rtmp
      f(t + c13*Δt,u + a1300*k1                                                                                            + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13,rtmp); k14=Δt*rtmp
      f(t + c14*Δt,u + a1400*k1                                                                                            + a1408*k9 + a1409*k10 + a1410*k11 + a1411*k12 + a1412*k13 + a1413*k14,rtmp); k15=Δt*rtmp
      f(t + c15*Δt,u + a1500*k1                                                                                            + a1508*k9 + a1509*k10 + a1510*k11 + a1511*k12 + a1512*k13 + a1513*k14 + a1514*k15,rtmp); k16=Δt*rtmp
      f(t + c16*Δt,u + a1600*k1                                                                                            + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16,rtmp); k17=Δt*rtmp
      f(t + c17*Δt,u + a1700*k1                                                                                                                                                   + a1712*k13 + a1713*k14 + a1714*k15 + a1715*k16 + a1716*k17,rtmp); k18=Δt*rtmp
      f(t + c18*Δt,u + a1800*k1                                                                                                                                                   + a1812*k13 + a1813*k14 + a1814*k15 + a1815*k16 + a1816*k17 + a1817*k18,rtmp); k19=Δt*rtmp
      f(t + c19*Δt,u + a1900*k1                                                                                                                                                   + a1912*k13 + a1913*k14 + a1914*k15 + a1915*k16 + a1916*k17 + a1917*k18 + a1918*k19,rtmp); k20=Δt*rtmp
      f(t + c20*Δt,u + a2000*k1                                                                                                                                                   + a2012*k13 + a2013*k14 + a2014*k15 + a2015*k16 + a2016*k17 + a2017*k18 + a2018*k19 + a2019*k20,rtmp); k21=Δt*rtmp
      f(t + c21*Δt,u + a2100*k1                                                                                                                                                   + a2112*k13 + a2113*k14 + a2114*k15 + a2115*k16 + a2116*k17 + a2117*k18 + a2118*k19 + a2119*k20 + a2120*k21,rtmp); k22=Δt*rtmp
      f(t + c22*Δt,u + a2200*k1                                                                                                                                                   + a2212*k13 + a2213*k14 + a2214*k15 + a2215*k16 + a2216*k17 + a2217*k18 + a2218*k19 + a2219*k20 + a2220*k21 + a2221*k22,rtmp); k23=Δt*rtmp
      f(t + c23*Δt,u + a2300*k1                                                                                            + a2308*k9 + a2309*k10 + a2310*k11 + a2311*k12 + a2312*k13 + a2313*k14 + a2314*k15 + a2315*k16 + a2316*k17 + a2317*k18 + a2318*k19 + a2319*k20 + a2320*k21 + a2321*k22 + a2322*k23,rtmp); k24=Δt*rtmp
      f(t + c24*Δt,u + a2400*k1                                                                                            + a2408*k9 + a2409*k10 + a2410*k11 + a2411*k12 + a2412*k13 + a2413*k14 + a2414*k15 + a2415*k16 + a2416*k17 + a2417*k18 + a2418*k19 + a2419*k20 + a2420*k21 + a2421*k22 + a2422*k23 + a2423*k24,rtmp); k25=Δt*rtmp
      f(t + c25*Δt,u + a2500*k1                                                                                            + a2508*k9 + a2509*k10 + a2510*k11 + a2511*k12 + a2512*k13 + a2513*k14 + a2514*k15 + a2515*k16 + a2516*k17 + a2517*k18 + a2518*k19 + a2519*k20 + a2520*k21 + a2521*k22 + a2522*k23 + a2523*k24 + a2524*k25,rtmp); k26=Δt*rtmp
      f(t + c26*Δt,u + a2600*k1                                                     + a2605*k6 + a2606*k7 + a2607*k8 + a2608*k9 + a2609*k10 + a2610*k11               + a2612*k13 + a2613*k14 + a2614*k15 + a2615*k16 + a2616*k17 + a2617*k18 + a2618*k19 + a2619*k20 + a2620*k21 + a2621*k22 + a2622*k23 + a2623*k24 + a2624*k25 + a2625*k26,rtmp); k27=Δt*rtmp
      f(t + c27*Δt,u + a2700*k1                                                     + a2705*k6 + a2706*k7 + a2707*k8 + a2708*k9 + a2709*k10               + a2711*k12 + a2712*k13 + a2713*k14 + a2714*k15 + a2715*k16 + a2716*k17 + a2717*k18 + a2718*k19 + a2719*k20 + a2720*k21 + a2721*k22 + a2722*k23 + a2723*k24 + a2724*k25 + a2725*k26 + a2726*k27,rtmp); k28=Δt*rtmp
      f(t + c28*Δt,u + a2800*k1                                                     + a2805*k6 + a2806*k7 + a2807*k8 + a2808*k9               + a2810*k11 + a2811*k12               + a2813*k14 + a2814*k15 + a2815*k16                                                                                                   + a2823*k24 + a2824*k25 + a2825*k26 + a2826*k27 + a2827*k28,rtmp); k29=Δt*rtmp
      f(t + c29*Δt,u + a2900*k1                                        + a2904*k5 + a2905*k6 + a2906*k7                           + a2909*k10 + a2910*k11 + a2911*k12               + a2913*k14 + a2914*k15 + a2915*k16                                                                                                   + a2923*k24 + a2924*k25 + a2925*k26 + a2926*k27 + a2927*k28 + a2928*k29,rtmp); k30=Δt*rtmp
      f(t + c30*Δt,u + a3000*k1                           + a3003*k4 + a3004*k5 + a3005*k6              + a3007*k8              + a3009*k10 + a3010*k11                             + a3013*k14 + a3014*k15 + a3015*k16                                                                                                   + a3023*k24 + a3024*k25 + a3025*k26               + a3027*k28 + a3028*k29 + a3029*k30,rtmp); k31=Δt*rtmp
      f(t + c31*Δt,u + a3100*k1              + a3102*k3 + a3103*k4                           + a3106*k7 + a3107*k8              + a3109*k10 + a3110*k11                             + a3113*k14 + a3114*k15 + a3115*k16                                                                                                   + a3123*k24 + a3124*k25 + a3125*k26               + a3127*k28 + a3128*k29 + a3129*k30 + a3130*k31,rtmp); k32=Δt*rtmp
      f(t + c32*Δt,u + a3200*k1 + a3201*k2                           + a3204*k5              + a3206*k7                                                                                                                                                                                                                                                                                                                                 + a3230*k31 + a3231*k32,rtmp); k33=Δt*rtmp
      f(t + c33*Δt,u + a3300*k1              + a3302*k3                                                                                                                                                                                                                                                                                                                                                                                                                 + a3332*k33,rtmp); k34=Δt*rtmp
      f(t + c34*Δt,u + a3400*k1 + a3401*k2 + a3402*k3              + a3404*k5              + a3406*k7 + a3407*k8              + a3409*k10 + a3410*k11 + a3411*k12 + a3412*k13 + a3413*k14 + a3414*k15 + a3415*k16 + a3416*k17 + a3417*k18 + a3418*k19 + a3419*k20 + a3420*k21 + a3421*k22 + a3422*k23 + a3423*k24 + a3424*k25 + a3425*k26 + a3426*k27 + a3427*k28 + a3428*k29 + a3429*k30 + a3430*k31 + a3431*k32 + a3432*k33 + a3433*k34,rtmp); k35=Δt*rtmp
      for i in uidx
        update[i] = (b1*k1[i] + b2*k2[i] + b3*k3[i] + b5*k5[i]) + (b7*k7[i] + b8*k8[i] + b10*k10[i] + b11*k11[i]) + (b12*k12[i] + b14*k14[i] + b15*k15[i] + b16*k16[i]) + (b18*k18[i] + b19*k19[i] + b20*k20[i] + b21*k21[i]) + (b22*k22[i] + b23*k23[i] + b24*k24[i] + b25*k25[i]) + (b26*k26[i] + b27*k27[i] + b28*k28[i] + b29*k29[i]) + (b30*k30[i] + b31*k31[i] + b32*k32[i] + b33*k33[i]) + (b34*k34[i] + b35*k35[i])
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + update[i]
        end
        EEst = norm(((k2 - k34) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting, so in-place
        for i in uidx
          u[i] = u[i] + update[i]
        end
      end
      @ode_loopfooter
    end
  end
  if dense
    f(t,u,k)
    push!(ks,deepcopy(k))
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:Feagin14,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,b35 = constructFeagin14(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType; local k5::uType
  local k6::uType; local k7::uType; local k8::uType; local k9::uType; local k10::uType
  local k11::uType; local k12::uType; local k13::uType; local k14::uType
  local k15::uType; local k16::uType; local k17::uType; local k18::uType
  local k19::uType; local k20::uType; local k21::uType; local k22::uType
  local k23::uType; local k24::uType; local k25::uType
  local k26::uType; local k27::uType; local k28::uType
  local k29::uType; local k30::uType; local k31::uType; local k32::uType
  local k33::uType; local k34::uType; local k35::uType
  if dense
    pop!(ks)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k   = f(t,u)
      k1  = Δt*k
      k2  = Δt*f(t + c1*Δt,u + a0100*k1)
      k3  = Δt*f(t + c2*Δt ,u + a0200*k1 + a0201*k2)
      k4  = Δt*f(t + c3*Δt,u + a0300*k1              + a0302*k3)
      k5  = Δt*f(t + c4*Δt,u + a0400*k1              + a0402*k3 + a0403*k4)
      k6  = Δt*f(t + c5*Δt,u + a0500*k1                           + a0503*k4 + a0504*k5)
      k7  = Δt*f(t + c6*Δt,u + a0600*k1                           + a0603*k4 + a0604*k5 + a0605*k6)
      k8  = Δt*f(t + c7*Δt,u + a0700*k1                                        + a0704*k5 + a0705*k6 + a0706*k7)
      k9  = Δt*f(t + c8*Δt,u + a0800*k1                                                     + a0805*k6 + a0806*k7 + a0807*k8)
      k10 = Δt*f(t + c9*Δt,u + a0900*k1                                                     + a0905*k6 + a0906*k7 + a0907*k8 + a0908*k9)
      k11 = Δt*f(t + c10*Δt,u + a1000*k1                                                     + a1005*k6 + a1006*k7 + a1007*k8 + a1008*k9 + a1009*k10)
      k12 = Δt*f(t + c11*Δt,u + a1100*k1                                                     + a1105*k6 + a1106*k7 + a1107*k8 + a1108*k9 + a1109*k10 + a1110*k11)
      k13 = Δt*f(t + c12*Δt,u + a1200*k1                                                                                            + a1208*k9 + a1209*k10 + a1210*k11 + a1211*k12)
      k14 = Δt*f(t + c13*Δt,u + a1300*k1                                                                                            + a1308*k9 + a1309*k10 + a1310*k11 + a1311*k12 + a1312*k13)
      k15 = Δt*f(t + c14*Δt,u + a1400*k1                                                                                            + a1408*k9 + a1409*k10 + a1410*k11 + a1411*k12 + a1412*k13 + a1413*k14)
      k16 = Δt*f(t + c15*Δt,u + a1500*k1                                                                                            + a1508*k9 + a1509*k10 + a1510*k11 + a1511*k12 + a1512*k13 + a1513*k14 + a1514*k15)
      k17 = Δt*f(t + c16*Δt,u + a1600*k1                                                                                            + a1608*k9 + a1609*k10 + a1610*k11 + a1611*k12 + a1612*k13 + a1613*k14 + a1614*k15 + a1615*k16)
      k18 = Δt*f(t + c17*Δt,u + a1700*k1                                                                                                                                                   + a1712*k13 + a1713*k14 + a1714*k15 + a1715*k16 + a1716*k17)
      k19 = Δt*f(t + c18*Δt,u + a1800*k1                                                                                                                                                   + a1812*k13 + a1813*k14 + a1814*k15 + a1815*k16 + a1816*k17 + a1817*k18)
      k20 = Δt*f(t + c19*Δt,u + a1900*k1                                                                                                                                                   + a1912*k13 + a1913*k14 + a1914*k15 + a1915*k16 + a1916*k17 + a1917*k18 + a1918*k19)
      k21 = Δt*f(t + c20*Δt,u + a2000*k1                                                                                                                                                   + a2012*k13 + a2013*k14 + a2014*k15 + a2015*k16 + a2016*k17 + a2017*k18 + a2018*k19 + a2019*k20)
      k22 = Δt*f(t + c21*Δt,u + a2100*k1                                                                                                                                                   + a2112*k13 + a2113*k14 + a2114*k15 + a2115*k16 + a2116*k17 + a2117*k18 + a2118*k19 + a2119*k20 + a2120*k21)
      k23 = Δt*f(t + c22*Δt,u + a2200*k1                                                                                                                                                   + a2212*k13 + a2213*k14 + a2214*k15 + a2215*k16 + a2216*k17 + a2217*k18 + a2218*k19 + a2219*k20 + a2220*k21 + a2221*k22)
      k24 = Δt*f(t + c23*Δt,u + a2300*k1                                                                                            + a2308*k9 + a2309*k10 + a2310*k11 + a2311*k12 + a2312*k13 + a2313*k14 + a2314*k15 + a2315*k16 + a2316*k17 + a2317*k18 + a2318*k19 + a2319*k20 + a2320*k21 + a2321*k22 + a2322*k23)
      k25 = Δt*f(t + c24*Δt,u + a2400*k1                                                                                            + a2408*k9 + a2409*k10 + a2410*k11 + a2411*k12 + a2412*k13 + a2413*k14 + a2414*k15 + a2415*k16 + a2416*k17 + a2417*k18 + a2418*k19 + a2419*k20 + a2420*k21 + a2421*k22 + a2422*k23 + a2423*k24)
      k26 = Δt*f(t + c25*Δt,u + a2500*k1                                                                                            + a2508*k9 + a2509*k10 + a2510*k11 + a2511*k12 + a2512*k13 + a2513*k14 + a2514*k15 + a2515*k16 + a2516*k17 + a2517*k18 + a2518*k19 + a2519*k20 + a2520*k21 + a2521*k22 + a2522*k23 + a2523*k24 + a2524*k25)
      k27 = Δt*f(t + c26*Δt,u + a2600*k1                                                     + a2605*k6 + a2606*k7 + a2607*k8 + a2608*k9 + a2609*k10 + a2610*k11               + a2612*k13 + a2613*k14 + a2614*k15 + a2615*k16 + a2616*k17 + a2617*k18 + a2618*k19 + a2619*k20 + a2620*k21 + a2621*k22 + a2622*k23 + a2623*k24 + a2624*k25 + a2625*k26)
      k28 = Δt*f(t + c27*Δt,u + a2700*k1                                                     + a2705*k6 + a2706*k7 + a2707*k8 + a2708*k9 + a2709*k10               + a2711*k12 + a2712*k13 + a2713*k14 + a2714*k15 + a2715*k16 + a2716*k17 + a2717*k18 + a2718*k19 + a2719*k20 + a2720*k21 + a2721*k22 + a2722*k23 + a2723*k24 + a2724*k25 + a2725*k26 + a2726*k27)
      k29 = Δt*f(t + c28*Δt,u + a2800*k1                                                     + a2805*k6 + a2806*k7 + a2807*k8 + a2808*k9               + a2810*k11 + a2811*k12               + a2813*k14 + a2814*k15 + a2815*k16                                                                                                   + a2823*k24 + a2824*k25 + a2825*k26 + a2826*k27 + a2827*k28)
      k30 = Δt*f(t + c29*Δt,u + a2900*k1                                        + a2904*k5 + a2905*k6 + a2906*k7                           + a2909*k10 + a2910*k11 + a2911*k12               + a2913*k14 + a2914*k15 + a2915*k16                                                                                                   + a2923*k24 + a2924*k25 + a2925*k26 + a2926*k27 + a2927*k28 + a2928*k29)
      k31 = Δt*f(t + c30*Δt,u + a3000*k1                           + a3003*k4 + a3004*k5 + a3005*k6              + a3007*k8              + a3009*k10 + a3010*k11                             + a3013*k14 + a3014*k15 + a3015*k16                                                                                                   + a3023*k24 + a3024*k25 + a3025*k26               + a3027*k28 + a3028*k29 + a3029*k30)
      k32 = Δt*f(t + c31*Δt,u + a3100*k1              + a3102*k3 + a3103*k4                           + a3106*k7 + a3107*k8              + a3109*k10 + a3110*k11                             + a3113*k14 + a3114*k15 + a3115*k16                                                                                                   + a3123*k24 + a3124*k25 + a3125*k26               + a3127*k28 + a3128*k29 + a3129*k30 + a3130*k31)
      k33 = Δt*f(t + c32*Δt,u + a3200*k1 + a3201*k2                           + a3204*k5              + a3206*k7                                                                                                                                                                                                                                                                                                                                 + a3230*k31 + a3231*k32)
      k34 = Δt*f(t + c33*Δt,u + a3300*k1              + a3302*k3                                                                                                                                                                                                                                                                                                                                                                                                                 + a3332*k33)
      k35 = Δt*f(t + c34*Δt,u + a3400*k1 + a3401*k2 + a3402*k3              + a3404*k5              + a3406*k7 + a3407*k8              + a3409*k10 + a3410*k11 + a3411*k12 + a3412*k13 + a3413*k14 + a3414*k15 + a3415*k16 + a3416*k17 + a3417*k18 + a3418*k19 + a3419*k20 + a3420*k21 + a3421*k22 + a3422*k23 + a3423*k24 + a3424*k25 + a3425*k26 + a3426*k27 + a3427*k28 + a3428*k29 + a3429*k30 + a3430*k31 + a3431*k32 + a3432*k33 + a3433*k34)
      update = (b1*k1 + b2*k2 + b3*k3 + b5*k5) + (b7*k7 + b8*k8 + b10*k10 + b11*k11) + (b12*k12 + b14*k14 + b15*k15 + b16*k16) + (b18*k18 + b19*k19 + b20*k20 + b21*k21) + (b22*k22 + b23*k23 + b24*k24 + b25*k25) + (b26*k26 + b27*k27 + b28*k28 + b29*k29) + (b30*k30 + b31*k31 + b32*k32 + b33*k33) + (b34*k34 + b35*k35)
      if adaptive
        utmp = u + update
        EEst = norm(((k2 - k34) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting, so in-place
        u = u + update
      end
      @ode_numberloopfooter
    end
  end
  if dense
    k = f(t,u)
    push!(ks,k)
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:ImplicitEuler,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  local nlres::NLsolve.SolverResults{uEltype}
  function rhs_ie(u,resid,u_old,t,Δt)
    resid[1] = u[1] - u_old[1] - Δt*f(t+Δt,u)[1]
  end
  uhold::Vector{uType} = Vector{uType}(1)
  u_old::Vector{uType} = Vector{uType}(1)
  uhold[1] = u; u_old[1] = u
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      u_old[1] = uhold[1]
      nlres = NLsolve.nlsolve((uhold,resid)->rhs_ie(uhold,resid,u_old,t,Δt),uhold,autodiff=autodiff)
      uhold[1] = nlres.zero[1]
      if dense
        k = f(t+Δt,uhold[1])
      end
      @ode_numberimplicitloopfooter
    end
  end
  u = uhold[1]
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:ImplicitEuler,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  local nlres::NLsolve.SolverResults{uEltype}
  uidx = eachindex(u)
  if autodiff
    cache = DiffCache(u)
    rhs_ie = (u,resid,u_old,t,Δt,cache) -> begin
      du = get_du(cache, eltype(u))
      f(t+Δt,reshape(u,sizeu),du)
      for i in uidx
        resid[i] = u[i] - u_old[i] - Δt*du[i]
      end
    end
  else
    cache = similar(u)
    rhs_ie = (u,resid,u_old,t,Δt,du) -> begin
      f(t+Δt,reshape(u,sizeu),du)
      for i in uidx
        resid[i] = u[i] - u_old[i] - Δt*du[i]
      end
    end
  end

  uhold = vec(u); u_old = similar(u)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      copy!(u_old,uhold)
      nlres = NLsolve.nlsolve((uhold,resid)->rhs_ie(uhold,resid,u_old,t,Δt,cache),uhold,autodiff=autodiff)
      uhold[:] = nlres.zero
      if dense
        f(t+Δt,u,k)
      end
      @ode_implicitloopfooter
    end
  end
  u = reshape(uhold,sizeu...)
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Trapezoid,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  local nlres::NLsolve.SolverResults{uEltype}
  uidx = eachindex(u)
  if autodiff
    cache1 = DiffCache(u)
    cache2 = DiffCache(u)
    Δto2 = Δt/2
    rhs_trap = (u,resid,u_old,t,Δt,cache1,cache2) -> begin
      du1 = get_du(cache1, eltype(u)); du2 = get_du(cache2, eltype(u_old))
      f(t,reshape(u_old,sizeu),du2)
      f(t+Δt,reshape(u,sizeu),du1)
      for i in uidx
        resid[i] = u[i] - u_old[i] - Δto2*(du1[i]+du2[i])
      end
    end
  else
    cache1 = similar(u)
    cache2 = similar(u)
    Δto2 = Δt/2
    rhs_trap = (u,resid,u_old,t,Δt,du1,du2) -> begin
      f(t,reshape(u_old,sizeu),du2)
      f(t+Δt,reshape(u,sizeu),du1)
      for i in uidx
        resid[i] = u[i] - u_old[i] - Δto2*(du1[i]+du2[i])
      end
    end
  end
  uhold = vec(u); u_old = similar(u)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      copy!(u_old,uhold)
      nlres = NLsolve.nlsolve((uhold,resid)->rhs_trap(uhold,resid,u_old,t,Δt,cache1,cache2),uhold,autodiff=autodiff)
      uhold[:] = nlres.zero
      if dense
        f(t+Δt,u,k)
      end
      @ode_implicitloopfooter
    end
  end
  u = reshape(uhold,sizeu...)
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:Trapezoid,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  Δto2::tType = Δt/2
  function rhs_trap(u,resid,u_old,t,Δt)
    resid[1] = u[1] - u_old[1] - Δto2*(f(t,u_old)[1] + f(t+Δt,u)[1])
  end
  local nlres::NLsolve.SolverResults{uEltype}
  uhold::Vector{uType} = Vector{uType}(1)
  u_old::Vector{uType} = Vector{uType}(1)
  uhold[1] = u; u_old[1] = u
  @inbounds for T in Ts
      while t < T
      @ode_loopheader
      u_old[1] = uhold[1]
      nlres = NLsolve.nlsolve((uhold,resid)->rhs_trap(uhold,resid,u_old,t,Δt),uhold,autodiff=autodiff)
      uhold[1] = nlres.zero[1]
      if dense
        k = f(t+Δt,uhold[1])
      end
      @ode_numberimplicitloopfooter
    end
  end
  u = uhold[1]
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray,ksEltype}(integrator::ODEIntegrator{:Rosenbrock32,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  local k₁::uType = similar(u)
  local k₂ = similar(u)
  local k₃::uType = similar(u)
  local tmp::uType
  function vecf(t,u,du)
    f(t,reshape(u,sizeu...),reshape(du,sizeu...))
    u = vec(u)
    du = vec(du)
  end
  function vecfreturn(t,u,du)
    f(t,reshape(u,sizeu...),reshape(du,sizeu...))
    return vec(du)
  end
  du1 = similar(u)
  du2 = similar(u)
  f₀ = similar(u)
  f₁ = similar(u)
  f₂ = similar(u); vectmp3 = similar(vec(u))
  utmp = similar(u); vectmp2 = similar(vec(u))
  dT = similar(u); vectmp = similar(vec(u))
  J = Matrix{uEltype}(length(u),length(u))
  W = similar(J); tmp2 = similar(u)
  uidx = eachindex(u)
  jidx = eachindex(J)
  f(t,u,fsalfirst)
  if dense
    k = fsalfirst
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      ForwardDiff.derivative!(dT,(t)->vecfreturn(t,u,du2),t) # Time derivative
      ForwardDiff.jacobian!(J,(du1,u)->vecf(t,u,du1),vec(du1),vec(u))

      W[:] = one(J)-Δt*d*J # Can an allocation be cut here?
      @into! vectmp = W\vec(fsalfirst + Δt*d*dT)
      k₁ = reshape(vectmp,sizeu...)
      for i in uidx
        utmp[i]=u[i]+Δt*k₁[i]/2
      end
      f(t+Δt/2,utmp,f₁)
      @into! vectmp2 = W\vec(f₁-k₁)
      tmp = reshape(vectmp2,sizeu...)
      for i in uidx
        k₂[i] = tmp[i] + k₁[i]
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + Δt*k₂[i]
        end
        f(t+Δt,utmp,fsallast)
        @into! vectmp3 = W\vec(f₂ - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+Δt*d*T)
        k₃ = reshape(vectmp3,sizeu...)
        for i in uidx
          tmp2[i] = (Δt*(k₁[i] - 2k₂[i] + k₃[i])/6)./(abstol+u[i]*reltol)
        end
        EEst = norm(tmp2,internalnorm)
      else
        for i in uidx
          u[i] = u[i] + Δt*k₂[i]
        end
        f(t,u,fsalfirst)
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts,ks
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number,ksEltype}(integrator::ODEIntegrator{:Rosenbrock32,uType,uEltype,N,tType,uEltypeNoUnits,rateType,ksEltype})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  local dT::uType
  local J::uType
  #f₀ = fsalfirst
  local k₁::uType
  local f₁::uType
  #f₂ = fsallast
  local k₂::uType
  local k₃::uType
  fsalfirst = f(t,u)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      # Time derivative
      dT = ForwardDiff.derivative((t)->f(t,u),t)
      J = ForwardDiff.derivative((u)->f(t,u),u)
      W = one(J)-Δt*d*J
      #f₀ = f(t,u)
      if dense
        k = fsalfirst
      end
      k₁ = W\(fsalfirst + Δt*d*dT)
      f₁ = f(t+Δt/2,u+Δt*k₁/2)
      k₂ = W\(f₁-k₁) + k₁
      if adaptive
        utmp = u + Δt*k₂
        fsallast = f(t+Δt,utmp)
        k₃ = W\(fsallast - c₃₂*(k₂-f₁)-2(k₁-fsalfirst)+Δt*d*T)
        EEst = norm((Δt*(k₁ - 2k₂ + k₃)/6)./(abstol+u*reltol),internalnorm)
      else
        u = u + Δt*k₂
        fsalfirst = f(t,u)
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts,ks
end
