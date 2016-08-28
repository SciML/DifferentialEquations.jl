immutable ODEIntegrator{Alg,uType<:Union{AbstractArray,Number},uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Union{AbstractArray,Number}} <: DEIntegrator
  f::Function
  u::uType
  t::tType
  Δt::tType
  Ts::Vector{tType}
  maxiters::Int
  timeseries::Vector{uType}
  ts::Vector{tType}
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
end

@def ode_preamble begin
  local u::uType
  local t::tType
  local Δt::tType
  local Ts::Vector{tType}
  local adaptiveorder::Int
  @unpack integrator: f,u,t,Δt,Ts,maxiters,timeseries,ts,timeseries_steps,γ,qmax,qmin,save_timeseries,adaptive,progressbar,autodiff,adaptiveorder,order,atomloaded,progress_steps,β,timechoicealg,qoldinit,normfactor,fsal
  Tfinal = Ts[end]
  local iter::Int = 0
  sizeu = size(u)
  local utmp::uType
  if uType <: Number
    utmp = zero(uType)
    fsalfirst = zero(rateType)
    fsallast = zero(rateType)
  else
    utmp = zeros(u)
    fsalfirst = rateType(sizeu)
    fsallast = rateType(sizeu)
  end
  local standard::uEltype = 0
  local q::uEltypeNoUnits = 0
  local Δtpropose::tType = 0
  local q11::uEltypeNoUnits = 0
  #local k1::uType; local k7::uType
  local qold::uEltypeNoUnits = 0

  if typeof(β) <: SIUnits.SIQuantity
    β = β.val
  end
  expo1 = 1/order - 0.75β
  qminc = inv(qmin)
  qmaxc = inv(qmax)
  local Eest::uEltypeNoUnits = zero(uEltypeNoUnits)
  if adaptive
    @unpack integrator: abstol,reltol,qmax,Δtmax,Δtmin,internalnorm
  end
  (progressbar && atomloaded && iter%progress_steps==0) ? Main.Atom.progress(0) : nothing #Use Atom's progressbar if loaded
end

@def ode_loopheader begin
  iter += 1
  if iter > maxiters
    warn("Max Iters Reached. Aborting")
    # u = map((x)->oftype(x,NaN),u)
    return u,t,timeseries,ts
  end
  Δt = min(Δt,abs(T-t))
end

@def ode_savevalues begin
  if save_timeseries && iter%timeseries_steps==0
    push!(timeseries,copy(u))
    push!(ts,t)
  end
end

@def ode_implicitsavevalues begin
  if save_timeseries && iter%timeseries_steps==0
    push!(timeseries,copy(reshape(u,sizeu...)))
    push!(ts,t)
  end
end

@def ode_numberimplicitsavevalues begin
  if save_timeseries && iter%timeseries_steps==0
    push!(timeseries,uhold[1])
    push!(ts,t)
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
      println(EEst)
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

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:Euler,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      u = u + Δt.*f(t,u)
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Euler,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  du::rateType = rateType(sizeu)
  uidx = eachindex(u)
  @inbounds for T in Ts
      while t < T
      @ode_loopheader
      f(t,u,du)
      for i in uidx
        u[i] = u[i] + Δt*du[i]
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:Midpoint,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  halfΔt::tType = Δt/2
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      u = u + Δt.*f(t+halfΔt,u+halfΔt.*f(t,u))
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Midpoint,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  halfΔt::tType = Δt/2
  utilde::uType = similar(u)
  du::rateType = rateType(sizeu)
  uidx = eachindex(u)
  @inbounds for T in Ts
      while t < T
      @ode_loopheader
      f(t,u,du)
      for i in uidx
        utilde[i] = u[i]+halfΔt*du[i]
      end
      f(t+halfΔt,utilde,du)
      for i in uidx
        u[i] = u[i] + Δt*du[i]
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:RK4,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
      u = u + Δt*(k₁ + 2k₂ + 2k₃ + k₄)/6
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:RK4,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  halfΔt::tType = Δt/2
  k₁ = rateType(sizeu)
  k₂ = rateType(sizeu)
  k₃ = rateType(sizeu)
  k₄ = rateType(sizeu)
  println(k₁)
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
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:ExplicitRK,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  local A::Matrix{uEltypeNoUnits}
  local c::Vector{uEltypeNoUnits}
  local α::Vector{uEltypeNoUnits}
  local αEEst::Vector{uEltypeNoUnits}
  local stages::Int
  local rtmp::ratetype
  @unpack integrator.tableau: A,c,α,αEEst,stages
  A = A' # Transpose A to column major looping
  ks = Array{typeof(u)}(stages)
  if fsal
    fsalfirst = f(t,u)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      # Calc First
      if fsal
        ks[1] = Δt*fsalfirst
      else
        ks[1] = Δt*f(t,u)
      end
      # Calc Middle
      for i = 2:stages-1
        utilde = zero(u)
        for j = 1:i-1
          utilde += A[j,i]*ks[j]
        end
        ks[i] = Δt*f(t+c[i]*Δt,u+utilde);
      end
      #Calc Last
      utilde = zero(u)
      for j = 1:stages-1
        utilde += A[j,end]*ks[j]
      end
      fsallast = f(t+c[end]*Δt,u+utilde); ks[end]=Δt*fsallast # Uses fsallast as temp even if not fsal
      # Accumulate Result
      utilde = α[1]*ks[1]
      for i = 2:stages
        utilde += α[i]*ks[i]
      end
      if adaptive
        utmp = u + utilde
        uEEst = αEEst[1]*ks[1]
        for i = 2:stages
          uEEst += αEEst[i]*ks[i]
        end
        EEst = abs( (utilde-uEEst)/(abstol+max(u,utmp)*reltol))
      else
        u = u + utilde
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:ExplicitRK,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  local A::Matrix{uEltypeNoUnits}
  local c::Vector{uEltypeNoUnits}
  local α::Vector{uEltypeNoUnits}
  local αEEst::Vector{uEltypeNoUnits}
  local stages::Int
  uidx = eachindex(u)
  @unpack integrator.tableau: A,c,α,αEEst,stages
  A = A' # Transpose A to column major looping
  ks = Vector{typeof(u)}(0)
  for i = 1:stages
    push!(ks,similar(u))
  end
  utilde = similar(u)
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  rtmp = rateType(sizeu)
  utmp = zeros(u)
  uEEst = similar(u)
  if fsal
    f(t,u,fsalfirst)
  end
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      # First
      if fsal
        for k in uidx
          ks[1][k] = Δt*fsalfirst[k]
        end
      else
        f(t,u,rtmp)
        for k in uidx
          ks[1][k]=Δt*rtmp[k]
        end
      end
      # Middle
      for i = 2:stages-1
        utilde[:] = zero(eltype(u))
        for j = 1:i-1
          for k in uidx
            utilde[k] += A[j,i]*ks[j][k]
          end
        end
        for k in uidx
          tmp[k] = u[k]+utilde[k]
        end
        f(t+c[i]*Δt,tmp,rtmp)
        for k in uidx
          ks[i][k]=Δt*rtmp[k]
        end
      end
      #Last
      utilde[:] = zero(eltype(u))
      for j = 1:stages-1
        for k in uidx
          utilde[k] += A[j,end]*ks[j][k]
        end
      end
      for k in uidx
        tmp[k] = u[k]+utilde[k]
      end
      f(t+c[end]*Δt,tmp,fsallast); #fsallast is tmp even if not fsal
      for k in uidx
        ks[end][k]= Δt*fsallast[k]
      end
      #Accumulate
      utilde[:] = α[1]*ks[1]
      for i = 2:stages
        for k in uidx
          utilde[k] += α[i]*ks[i][k]
        end
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + utilde[i]
        end
        uEEst[:] = αEEst[1]*ks[1]
        for i = 2:stages
          for j in uidx
            uEEst[j] += αEEst[i]*ks[i][j]
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
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:ExplicitRKVectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:BS3,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:BS3Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:BS3,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  a21,a32,a41,a42,a43,c1,c2,b1,b2,b3,b4  = constructBS3(uEltypeNoUnits)
  k1 = similar(u)
  k2 = similar(u)
  k3 = similar(u)
  k4 = similar(u)
  utilde = similar(u)
  local rtmp = rateType(sizeu)
  uidx = eachindex(u)
  tmp = similar(u,uEltypeNoUnits)
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
        k3[i] = Δt*rtmp
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
        u = utmp
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:BS5,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat2,bhat3,bhat4,bhat5,bhat6,bhat7,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8  = constructBS5(uEltypeNoUnits)
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
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:BS5Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat2,bhat3,bhat4,bhat5,bhat6,bhat7,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8  = constructBS5(uEltypeNoUnits)
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
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:BS5,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat2,bhat3,bhat4,bhat5,bhat6,bhat7,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8  = constructBS5(uEltypeNoUnits)
  k1::uType = similar(u)
  k2::uType = similar(u)
  k3::uType = similar(u)
  k4::uType = similar(u)
  k5::uType = similar(u)
  k6::uType = similar(u)
  k7::uType = similar(u)
  k8::uType = similar(u)
  utilde = similar(u)
  uhat   = similar(u)
  local EEst2::uEltypeNoUnits
  uidx = eachindex(u)
  tmp = similar(u); atmp = similar(u); atmp = similar(u,uEltypeNoUnits); tmptilde = similar(u)
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
          tmp[i] = ((uhat[i]-utmp[i])./(abstol+max(u[i],utmp[i])*reltol))^2
          tmptilde[i] = ((utilde[i]-utmp[i])./(abstol+max(u[i],utmp[i])*reltol))^2
        end
        EEst1 = sqrt( sum(tmp) * normfactor)
        EEst2 = sqrt( sum(tmptilde) * normfactor)
        EEst = max(EEst1,EEst2)
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:Tsit5,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Tsit5Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Tsit5,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
        u = utmp
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:DP5,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = constructDP5(uEltypeNoUnits)
  local k1::uType
  local k2::uType
  local k3::uType
  local k4::uType
  local k5::uType
  local k6::uType
  local k7::uType
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
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:DP5Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:DP5,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
        u = utmp
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:DP5Threaded,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  a21::uEltype,a31::uEltype,a32::uEltype,a41::uEltype,a42::uEltype,a43::uEltype,a51::uEltype,a52::uEltype,a53::uEltype,a54::uEltype,a61::uEltype,a62::uEltype,a63::uEltype,a64::uEltype,a65::uEltype,a71::uEltype,a73::uEltype,a74::uEltype,a75::uEltype,a76::uEltype,b1::uEltype,b3::uEltype,b4::uEltype,b5::uEltype,b6::uEltype,b7::uEltype,c1::uEltype,c2::uEltype,c3::uEltype,c4::uEltype,c5::uEltype,c6::uEltype = constructDP5(uEltypeNoUnits)
  k1::uType = similar(u)
  k2::uType = similar(u)
  k3::uType = similar(u)
  k4::uType = similar(u)
  k5::uType = similar(u)
  k6::uType = similar(u)
  k7::uType = similar(u)
  utilde = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  uidx::Base.OneTo{Int64} = eachindex(u)
  f(t,u,fsalfirst);  # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      dp5threaded_loop1(rtmp,Δt,k1,fsalfirst,u,a21,tmp,Δt,uidx)
      f(t+c1*Δt,tmp,rtmp);
      dp5threaded_loop2(rtmp,Δt,tmp,u,a31,k1,a32,k2,uidx)
      f(t+c2*Δt,tmp,rtmp);
      dp5threaded_loop3(rtmp,Δt,tmp,u,a41,k1,a42,k2,a43,k3,uidx)
      f(t+c3*Δt,tmp,rtmp);
      dp5threaded_loop4(rtmp,Δt,tmp,u,a51,k1,a52,k2,a53,k3,a54,k4,uidx)
      f(t+c4*Δt,tmp,rtmp);
      dp5threaded_loop5(rtmp,Δt,tmp,u,a61,k1,a62,k2,a63,k3,a64,k4,a65,k5,uidx)
      f(t+Δt,tmp,rtmp);
      dp5threaded_loop6(rtmp,Δt,utmp,u,a71,k1,a73,k3,a74,k4,a75,k5,a76,k6,uidx)
      f(t+Δt,utmp,fsallast);
      dp5threaded_loop7(k7,Δt,fsallast,uidx)
      if adaptive
        dp5threaded_adaptiveloop(utilde,u,b1,k1,b3,k3,b4,k4,b5,k5,b6,k6,b7,k7,atmp,utmp,abstol,reltol,uidx)
        EEst = sqrt( sum(atmp) * normfactor)
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Float64,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:DP5Threaded,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  a21::uEltype,a31::uEltype,a32::uEltype,a41::uEltype,a42::uEltype,a43::uEltype,a51::uEltype,a52::uEltype,a53::uEltype,a54::uEltype,a61::uEltype,a62::uEltype,a63::uEltype,a64::uEltype,a65::uEltype,a71::uEltype,a73::uEltype,a74::uEltype,a75::uEltype,a76::uEltype,b1::uEltype,b3::uEltype,b4::uEltype,b5::uEltype,b6::uEltype,b7::uEltype,c1::uEltype,c2::uEltype,c3::uEltype,c4::uEltype,c5::uEltype,c6::uEltype = constructDP5(uEltypeNoUnits)
  k1::uType = similar(u)
  k2::uType = similar(u)
  k3::uType = similar(u)
  k4::uType = similar(u)
  k5::uType = similar(u)
  k6::uType = similar(u)
  k7::uType = similar(u)
  utilde = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  uidx::Base.OneTo{Int64} = eachindex(u)
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
        u = utmp
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
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


function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:Vern6,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,b1,b4,b5,b6,b7,b8,b9= constructVern6(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType;
  local k5::uType; local k6::uType; local k7::uType; local k8::uType;
  local utilde::uType; fsalfirst = f(t,u) # Pre-start fsal
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
      fsallast = Δt*f(t+Δt,utmp); k9 = Δt*fsallast
      if adaptive
        utilde = u + b1*k1 + b4*k4 + b5*k5 + b6*k6 + b7*k7 + b8*k8 + b9*k9
        EEst = abs( ((utilde-utmp)/(abstol+max(u,utmp)*reltol)) * normfactor)
      else
        u = utmp
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Vern6Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,b1,b4,b5,b6,b7,b8,b9= constructVern6(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u) ; k3 = similar(u); k4 = similar(u)
  k5 = similar(u); k6 = similar(u) ; k7 = similar(u); k8 = similar(u)
  utilde = similar(u); rtmp = rateType(sizeu)
  fsalfirst = f(t,u,du) # Pre-start fsal
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
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Vern6,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,b1,b4,b5,b6,b7,b8,b9= constructVern6(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u)
  k5 = similar(u); k6 = similar(u); k7 = similar(u); k8 = similar(u)
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits); uidx = eachindex(u)
  fsalfirst = f(t,u,du) # Pre-start fsal
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      for i in uidx
        k1[i] = Δt*fsalfirst[i]
        tmp[i] = u[i]+a21*k1[i]
      end
      f(t+c1*Δt,tmp,k2)
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
        u[i]+a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]
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
      f(t+Δt,utmp,fsallast); k9 = Δt*fsallast
      if adaptive
        for i in uidx
          utilde[i] = u[i] + b1*k1[i] + b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i] + b8*k8[i] + b9*k9[i]
          atmp[i] = ((utilde[i]-utmp[i])/(abstol+max(u[i],utmp[i])*reltol))^2
        end
        EEst = sqrt( sum(atmp) * normfactor)
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:TanYam7,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a93,a94,a95,a96,a97,a98,a101,a103,a104,a105,a106,a107,a108,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat8,bhat10 = constructVern6(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType;
  local k5::uType; local k6::uType; local k7::uType; local k8::uType;
  local utilde::uType; @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k1 = Δt*f(t,u)
      k2 = Δt*f(t+c1*Δt,u+a21*k1)
      k3 = Δt*f(t+c2*Δt,u+a31*k1+a32*k2)
      k4 = Δt*f(t+c3*Δt,u+a41*k1       +a43*k3)
      k5 = Δt*f(t+c4*Δt,u+a51*k1       +a53*k3+a54*k4)
      k6 = Δt*f(t+c5*Δt,u+a61*k1       +a63*k3+a64*k4+a65*k5)
      k6 = Δt*f(t+c6*Δt,u+a71*k1       +a73*k3+a74*k4+a75*k5+a76*k6)
      k7 = Δt*f(t+c7*Δt,u+a81*k1       +a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
      k8 = Δt*f(t+Δt,u+a91*k1       +a93*k3+a94*k4+a95*k5+a96*k6+a97*k7+a98*k8)
      k9 = Δt*f(t+Δt,u+a101*k1      +a103*k3+a104*k4+a105*k5+a106*k6+a107*k7+a108*k8)
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
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:TanYam7Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a93,a94,a95,a96,a97,a98,a101,a103,a104,a105,a106,a107,a108,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat8,bhat10 = constructVern6(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u) ; k3 = similar(u); k4 = similar(u)
  k5 = similar(u); k6 = similar(u) ; k7 = similar(u); k8 = similar(u)
  utilde = similar(u); rtmp = rateType(sizeu)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp); k1=Δt*rtmp
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
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:TanYam7,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,a21,a31,a32,a41,a43,a51,a53,a54,a61,a63,a64,a65,a71,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,a91,a93,a94,a95,a96,a97,a98,a101,a103,a104,a105,a106,a107,a108,b1,b4,b5,b6,b7,b8,b9,bhat1,bhat4,bhat5,bhat6,bhat7,bhat8,bhat10 = constructVern6(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u) ; k3 = similar(u); k4 = similar(u)
  k5 = similar(u); k6 = similar(u) ; k7 = similar(u); k8 = similar(u)
  utilde = similar(u); uidx = eachindex(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  rtmp = rateType(sizeu)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp)
      for i in uidx
        k1[i]=Δt*rtmp[i]
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
        u = utmp
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:DP8,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,c14,c15,c16,b1,b6,b7,b8,b9,b10,b11,b12,bhh1,bhh2,bhh3,er1,er6,er7,er8,er9,er10,er11,er12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = constructDP8(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType;
  local k5::uType; local k6::uType; local k7::uType; local k8::uType;
  local k9::uType; local k10::uType; local k11::uType; local k12::uType;
  local k13::uType; local utilde::uType;
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
        err3 = abs((update - bhh1*k1 - bhh2*k9 - bhh3*k3)/(abstol+max(u,utmp)*reltol) * normfactor) # Order 3
        err52 = err5*err5
        EEst = err52/sqrt(err52 + 0.01*err3*err3)
      else
        u = utmp
      end
      @ode_numberloopfooter
    end
  end
  # Dense output: a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:DP8Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,c14,c15,c16,b1,b6,b7,b8,b9,b10,b11,b12,bhh1,bhh2,bhh3,er1,er6,er7,er8,er9,er10,er11,er12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = constructDP8(uEltypeNoUnits)
  k1 = similar(u); k2  = similar(u); k3  = similar(u);  k4 = similar(u)
  k5 = similar(u); k6  = similar(u); k7  = similar(u);  k8 = similar(u)
  k9 = similar(u); k10 = similar(u); k11 = similar(u); k12 = similar(u)
  k13 = similar(u); utilde = similar(u); err5 = similar(u); err3 = similar(u)
  rtmp = rateType(sizeu)
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
        err3 = sqrt(sum(((update - bhh1*k1 - bhh2*k9 - bhh3*k3)./(abstol+max(u,utmp)*reltol)).^2) * normfactor) # Order 3
        err52 = err5*err5
        EEst = err52/sqrt(err52 + 0.01*err3*err3)
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  # Dense output: a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:DP8,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c7,c8,c9,c10,c11,c6,c5,c4,c3,c2,c14,c15,c16,b1,b6,b7,b8,b9,b10,b11,b12,bhh1,bhh2,bhh3,er1,er6,er7,er8,er9,er10,er11,er12,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615 = constructDP8(uEltypeNoUnits)
  k1 = similar(u); k2  = similar(u); k3  = similar(u);  k4 = similar(u)
  k5 = similar(u); k6  = similar(u); k7  = similar(u);  k8 = similar(u)
  k9 = similar(u); k10 = similar(u); k11 = similar(u); k12 = similar(u)
  k13 = similar(u); utilde = similar(u); err5 = similar(u); err3 = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); uidx = eachindex(u); atmp2 = similar(u,uEltypeNoUnits); update = similar(u)
  rtmp = simlar(u,rateType)
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
          atmp2[i]= ((update[i] - bhh1*k1[i] - bhh2*k9[i] - bhh3*k3[i])/(abstol+max(u[i],utmp[i])*reltol))^2
        end
        err5 = sqrt( sum(atmp)  * normfactor) # Order 5
        err3 = sqrt( sum(atmp2) * normfactor) # Order 3
        err52 = err5*err5
        EEst = err52/sqrt(err52 + 0.01*err3*err3)
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  # Dense output: a1401,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1613,a1614,a1615
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:TsitPap8,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13 = constructTsitPap8(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType;
  local k5::uType; local k6::uType; local k7::uType; local k8::uType;
  local k9::uType; local k10::uType; local k11::uType; local k12::uType;
  local k13::uType; local utilde::uType;
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
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:TsitPap8Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13 = constructTsitPap8(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u)
  k5 = similar(u); k6 = similar(u); k7 = similar(u); k8 = similar(u)
  k9 = similar(u); k10 = similar(u); k11 = similar(u); k12 = similar(u)
  k13::uType; utilde = similar(u);
  rtmp = rateType(sizeu)
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
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:TsitPap8,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0804,a0805,a0806,a0807,a0901,a0904,a0905,a0906,a0907,a0908,a1001,a1004,a1005,a1006,a1007,a1008,a1009,a1101,a1104,a1105,a1106,a1107,a1108,a1109,a1110,a1201,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1304,a1305,a1306,a1307,a1308,a1309,a1310,b1,b6,b7,b8,b9,b10,b11,b12,bhat1,bhat6,bhat7,bhat8,bhat9,bhat10,bhat13 = constructTsitPap8(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u); k3 = similar(u); k4 = similar(u)
  k5 = similar(u); k6 = similar(u); k7 = similar(u); k8 = similar(u)
  k9 = similar(u); k10 = similar(u); k11 = similar(u); k12 = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); uidx = eachindex(u)
  k13::uType; utilde = similar(u);
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
        tmp[i] = u+a0801*k1[i]+a0804*k4[i]+a0805*k5[i]+a0806*k6[i]+a0807*k7[i]
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
        tmp[i] = u+a1201*k1[i]+a1204*k4[i]+a1205*k5[i]+a1206*k6[i]+a1207*k7[i]+a1208*k8[i]+a1209*k9[i]+a1210*k10[i]+a1211*k11[i]
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
        u = utmp
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:Vern9,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0806,a0807,a0901,a0906,a0907,a0908,a1001,a1006,a1007,a1008,a1009,a1101,a1106,a1107,a1108,a1109,a1110,a1201,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1401,a1406,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,b1,b8,b9,b10,b11,b12,b13,b14,b15,bhat1,bhat8,bhat9,bhat10,bhat11,bhat12,bhat13,bhat16 = constructVern9(uEltypeNoUnits)
  local k1::uType; local k2::uType; local k3::uType; local k4::uType;
  local k5::uType; local k6::uType; local k7::uType; local k8::uType;
  local k9::uType; local k10::uType; local k11::uType; local k12::uType;
  local k13::uType; local k14::uType; local k15::uType; local k16::uType;
  local utilde::uType;
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
      k8 = Δt*f(t+c7*Δt,u+a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7)
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
        EEst = abs((update - k1*bhat1 - k8*bhat8 - k9*bhat9 - k10*bhat10 - k11*bhat11 - k12*bhat12 - k13*bhat13 - k16*bhat16)/(abstol+max(u,utmp)*reltol) * normfactor) # Order 5
      else
        u = utmp
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Vern9Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0806,a0807,a0901,a0906,a0907,a0908,a1001,a1006,a1007,a1008,a1009,a1101,a1106,a1107,a1108,a1109,a1110,a1201,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1401,a1406,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,b1,b8,b9,b10,b11,b12,b13,b14,b15,bhat1,bhat8,bhat9,bhat10,bhat11,bhat12,bhat13,bhat16 = constructVern9(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u);k3 = similar(u); k4 = similar(u);
  k5 = similar(u); k6 = similar(u);k7 = similar(u); k8 = similar(u);
  k9 = similar(u); k10 = simlar(u); k11 = similar(u); k12 = similar(u);
  k13 = similar(u); k14 = similar(u); k15 = similar(u); k16 =similar(u);
  utilde = similar(u);
  rtmp = rateType(sizeu)
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
      f(t+c7*Δt,u+a0801*k1                +a0804*k4+a0805*k5+a0806*k6+a0807*k7,rtmp); k8=Δt*rtmp
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
        EEst = sqrt(sum(((update - k1*bhat1 - k8*bhat8 - k9*bhat9 - k10*bhat10 - k11*bhat11 - k12*bhat12 - k13*bhat13 - k16*bhat16)./(abstol+max(u,utmp)*reltol)).^2) * normfactor) # Order 5
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Vern9,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,a0201,a0301,a0302,a0401,a0403,a0501,a0503,a0504,a0601,a0604,a0605,a0701,a0704,a0705,a0706,a0801,a0806,a0807,a0901,a0906,a0907,a0908,a1001,a1006,a1007,a1008,a1009,a1101,a1106,a1107,a1108,a1109,a1110,a1201,a1206,a1207,a1208,a1209,a1210,a1211,a1301,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1401,a1406,a1407,a1408,a1409,a1410,a1411,a1412,a1413,a1501,a1506,a1507,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1601,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,b1,b8,b9,b10,b11,b12,b13,b14,b15,bhat1,bhat8,bhat9,bhat10,bhat11,bhat12,bhat13,bhat16 = constructVern9(uEltypeNoUnits)
  k1 = similar(u); k2 = similar(u);k3 = similar(u); k4 = similar(u);
  k5 = similar(u); k6 = similar(u);k7 = similar(u); k8 = similar(u);
  k9 = similar(u); k10 = simlar(u); k11 = similar(u); k12 = similar(u);
  k13 = similar(u); k14 = similar(u); k15 = similar(u); k16 =similar(u);
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits); uidx = eachindex(u)
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
        tmp[i] = u[i]+a0801*k1[i]+a0804*k4[i]+a0805*k5[i]+a0806*k6[i]+a0807*k7[i]
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
          atmp[i] = ((update[i] - k1[i]*bhat1 - k8[i]*bhat8 - k9[i]*bhat9 - k10[i]*bhat10 - k11[i]*bhat11 - k12[i]*bhat12 - k13[i]*bhat13 - k16[i]*bhat16)/(abstol+max(u,utmp)*reltol))^2
        end
        EEst = sqrt(sum(atmp) * normfactor) # Order 5
      else
        u = utmp
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Feagin10Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b,c = constructFeagin10(uEltypeNoUnits)
  k = Vector{typeof(u)}(0)
  for i = 1:17
    push!(k,similar(u))
  end
  update = similar(u)
  utmp = similar(u)
  uidx = eachindex(u)
  rtmp = rateType(sizeu)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp); k[1]=Δt*rtmp
      f(t + c[1]*Δt,u + a0100*k[1],rtmp); k[2]=Δt*rtmp
      f(t + c[2]*Δt ,u + a0200*k[1] + a0201*k[2],rtmp); k[3]=Δt*rtmp
      f(t + c[3]*Δt,u + a0300*k[1]              + a0302*k[3],rtmp); k[4]=Δt*rtmp
      f(t + c[4]*Δt,u + a0400*k[1]              + a0402*k[3] + a0403*k[4],rtmp); k[5]=Δt*rtmp
      f(t + c[5]*Δt,u + a0500*k[1]                           + a0503*k[4] + a0504*k[5],rtmp); k[6]=Δt*rtmp
      f(t + c[6]*Δt,u + a0600*k[1]                           + a0603*k[4] + a0604*k[5] + a0605*k[6],rtmp); k[7]=Δt*rtmp
      f(t + c[7]*Δt,u + a0700*k[1]                                        + a0704*k[5] + a0705*k[6] + a0706*k[7],rtmp); k[8]=Δt*rtmp
      f(t + c[8]*Δt,u + a0800*k[1]                                                     + a0805*k[6] + a0806*k[7] + a0807*k[8],rtmp); k[9]=Δt*rtmp
      f(t + c[9]*Δt,u + a0900*k[1]                                                     + a0905*k[6] + a0906*k[7] + a0907*k[8] + a0908*k[9],rtmp); k[10]=Δt*rtmp
      f(t + c[10]*Δt,u + a1000*k[1]                                                     + a1005*k[6] + a1006*k[7] + a1007*k[8] + a1008*k[9] + a1009*k[10],rtmp); k[11]=Δt*rtmp
      f(t + c[11]*Δt,u + a1100*k[1]                                                     + a1105*k[6] + a1106*k[7] + a1107*k[8] + a1108*k[9] + a1109*k[10] + a1110*k[11],rtmp); k[12]=Δt*rtmp
      f(t + c[12]*Δt,u + a1200*k[1]                           + a1203*k[4] + a1204*k[5] + a1205*k[6] + a1206*k[7] + a1207*k[8] + a1208*k[9] + a1209*k[10] + a1210*k[11] + a1211*k[12],rtmp); k[13]=Δt*rtmp
      f(t + c[13]*Δt,u + a1300*k[1]              + a1302*k[3] + a1303*k[4]              + a1305*k[6] + a1306*k[7] + a1307*k[8] + a1308*k[9] + a1309*k[10] + a1310*k[11] + a1311*k[12] + a1312*k[13],rtmp); k[14]=Δt*rtmp
      f(t + c[14]*Δt,u + a1400*k[1] + a1401*k[2]                           + a1404*k[5]              + a1406*k[7] +                                                                     a1412*k[13] + a1413*k[14],rtmp); k[15]=Δt*rtmp
      f(t + c[15]*Δt,u + a1500*k[1]              + a1502*k[3]                                                                                                                                                     + a1514*k[15],rtmp); k[16]=Δt*rtmp
      f(t + c[16]*Δt,u + a1600*k[1] + a1601*k[2] + a1602*k[3]              + a1604*k[5] + a1605*k[6] + a1606*k[7] + a1607*k[8] + a1608*k[9] + a1609*k[10] + a1610*k[11] + a1611*k[12] + a1612*k[13] + a1613*k[14] + a1614*k[15] + a1615*k[16],rtmp); k[17]=Δt*rtmp
      for i in uidx
        update[i] = (b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i]) + (b[7]*k[7][i] + b[9]*k[9][i] + b[10]*k[10][i] + b[11]*k[11][i]) + (b[12]*k[12][i] + b[13]*k[13][i] + b[14]*k[14][i] + b[15]*k[15][i]) + (b[16]*k[16][i] + b[17]*k[17][i])
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + update[i]
        end
        EEst = norm(((k[2] - k[16]) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting, so in-place
        for i in uidx
          u[i] = u[i] + update[i]
        end
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Feagin10,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b,c = constructFeagin10(uEltypeNoUnits)
  k = Vector{typeof(u)}(0)
  for i = 1:17
    push!(k,similar(u))
  end
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  utmp = similar(u)
  uidx = eachindex(u)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp)
      for i in uidx
        k[1][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0100*k[1][i]
      end
      f(t + c[1]*Δt,tmp,rtmp)
      for i in uidx
        k[2][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0200*k[1][i] + a0201*k[2][i]
      end
      f(t + c[2]*Δt ,tmp,rtmp)
      for i in uidx
        k[3][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0300*k[1][i] + a0302*k[3][i]
      end
      f(t + c[3]*Δt,tmp,rtmp)
      for i in uidx
        k[4][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0400*k[1][i] + a0402*k[3][i] + a0403*k[4][i]
      end
      f(t + c[4]*Δt,tmp,rtmp)
      for i in uidx
        k[5][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0500*k[1][i] + a0503*k[4][i] + a0504*k[5][i]
      end
      f(t + c[5]*Δt,tmp,rtmp)
      for i in uidx
        k[6][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0600*k[1][i] + a0603*k[4][i] + a0604*k[5][i] + a0605*k[6][i]
      end
      f(t + c[6]*Δt,tmp,rtmp)
      for i in uidx
        k[7][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0700*k[1][i] + a0704*k[5][i] + a0705*k[6][i]) + a0706*k[7][i]
      end
      f(t + c[7]*Δt,tmp,rtmp)
      for i in uidx
        k[8][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0800*k[1][i] + a0805*k[6][i] + a0806*k[7][i]) + a0807*k[8][i]
      end
      f(t + c[8]*Δt,tmp,rtmp)
      for i in uidx
        k[9][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0900*k[1][i] + a0905*k[6][i] + a0906*k[7][i]) + a0907*k[8][i] + a0908*k[9][i]
      end
      f(t + c[9]*Δt,tmp,rtmp)
      for i in uidx
        k[10][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1000*k[1][i] + a1005*k[6][i] + a1006*k[7][i]) + a1007*k[8][i] + a1008*k[9][i] + a1009*k[10][i]
      end
      f(t + c[10]*Δt,tmp,rtmp)
      for i in uidx
        k[11][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1100*k[1][i] + a1105*k[6][i] + a1106*k[7][i]) + (a1107*k[8][i] + a1108*k[9][i] + a1109*k[10][i] + a1110*k[11][i])
      end
      f(t + c[11]*Δt,tmp,rtmp)
      for i in uidx
        k[12][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1200*k[1][i] + a1203*k[4][i] + a1204*k[5][i]) + (a1205*k[6][i] + a1206*k[7][i] + a1207*k[8][i] + a1208*k[9][i]) + (a1209*k[10][i] + a1210*k[11][i] + a1211*k[12][i])
      end
      f(t + c[12]*Δt,tmp,rtmp)
      for i in uidx
        k[13][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1300*k[1][i] + a1302*k[3][i] + a1303*k[4][i]) + (a1305*k[6][i] + a1306*k[7][i] + a1307*k[8][i] + a1308*k[9][i]) + (a1309*k[10][i] + a1310*k[11][i] + a1311*k[12][i] + a1312*k[13][i])
      end
      f(t + c[13]*Δt,tmp,rtmp)
      for i in uidx
        k[14][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1400*k[1][i] + a1401*k[2][i] + a1404*k[5][i]) + (a1406*k[7][i] + a1412*k[13][i] + a1413*k[14][i])
      end
      f(t + c[14]*Δt,tmp,rtmp)
      for i in uidx
        k[15][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a1500*k[1][i] + a1502*k[3][i] + a1514*k[15][i]
      end
      f(t + c[15]*Δt,tmp,rtmp)
      for i in uidx
        k[16][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1600*k[1][i] + a1601*k[2][i] + a1602*k[3][i]) + (a1604*k[5][i] + a1605*k[6][i] + a1606*k[7][i] + a1607*k[8][i]) + (a1608*k[9][i] + a1609*k[10][i] + a1610*k[11][i] + a1611*k[12][i]) + (a1612*k[13][i] + a1613*k[14][i] + a1614*k[15][i] + a1615*k[16][i])
      end
      f(t + c[16]*Δt,tmp,rtmp)
      for i in uidx
        k[17][i]=Δt*rtmp[i]
        tmp[i] = (b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i]) + (b[7]*k[7][i] + b[9]*k[9][i] + b[10]*k[10][i] + b[11]*k[11][i]) + (b[12]*k[12][i] + b[13]*k[13][i] + b[14]*k[14][i] + b[15]*k[15][i]) + (b[16]*k[16][i] + b[17]*k[17][i])
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + tmp[i]
          tmp[i] = ((k[2][i] - k[16][i]) * adaptiveConst)./(abstol+u[i]*reltol)
        end
        EEst = norm(atmp,internalnorm)
      else #no chance of rejecting, so in-place
        for i in uidx
          u[i] = u[i] + atmp[i]
        end
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:Feagin10,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1203,a1204,a1205,a1206,a1207,a1208,a1209,a1210,a1211,a1300,a1302,a1303,a1305,a1306,a1307,a1308,a1309,a1310,a1311,a1312,a1400,a1401,a1404,a1406,a1412,a1413,a1500,a1502,a1514,a1600,a1601,a1602,a1604,a1605,a1606,a1607,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,b,c = constructFeagin10(uEltypeNoUnits)
  k = Vector{typeof(u)}(17)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k[1]  = Δt*f(t,u)
      k[2]  = Δt*f(t + c[1]*Δt,u + a0100*k[1])
      k[3]  = Δt*f(t + c[2]*Δt ,u + a0200*k[1] + a0201*k[2])
      k[4]  = Δt*f(t + c[3]*Δt,u + a0300*k[1]              + a0302*k[3])
      k[5]  = Δt*f(t + c[4]*Δt,u + a0400*k[1]              + a0402*k[3] + a0403*k[4])
      k[6]  = Δt*f(t + c[5]*Δt,u + a0500*k[1]                           + a0503*k[4] + a0504*k[5])
      k[7]  = Δt*f(t + c[6]*Δt,u + a0600*k[1]                           + a0603*k[4] + a0604*k[5] + a0605*k[6])
      k[8]  = Δt*f(t + c[7]*Δt,u + a0700*k[1]                                        + a0704*k[5] + a0705*k[6] + a0706*k[7])
      k[9]  = Δt*f(t + c[8]*Δt,u + a0800*k[1]                                                     + a0805*k[6] + a0806*k[7] + a0807*k[8])
      k[10] = Δt*f(t + c[9]*Δt,u + a0900*k[1]                                                     + a0905*k[6] + a0906*k[7] + a0907*k[8] + a0908*k[9])
      k[11] = Δt*f(t + c[10]*Δt,u + a1000*k[1]                                                     + a1005*k[6] + a1006*k[7] + a1007*k[8] + a1008*k[9] + a1009*k[10])
      k[12] = Δt*f(t + c[11]*Δt,u + a1100*k[1]                                                     + a1105*k[6] + a1106*k[7] + a1107*k[8] + a1108*k[9] + a1109*k[10] + a1110*k[11])
      k[13] = Δt*f(t + c[12]*Δt,u + a1200*k[1]                           + a1203*k[4] + a1204*k[5] + a1205*k[6] + a1206*k[7] + a1207*k[8] + a1208*k[9] + a1209*k[10] + a1210*k[11] + a1211*k[12])
      k[14] = Δt*f(t + c[13]*Δt,u + a1300*k[1]              + a1302*k[3] + a1303*k[4]              + a1305*k[6] + a1306*k[7] + a1307*k[8] + a1308*k[9] + a1309*k[10] + a1310*k[11] + a1311*k[12] + a1312*k[13])
      k[15] = Δt*f(t + c[14]*Δt,u + a1400*k[1] + a1401*k[2]                           + a1404*k[5]              + a1406*k[7] +                                                                     a1412*k[13] + a1413*k[14])
      k[16] = Δt*f(t + c[15]*Δt,u + a1500*k[1]              + a1502*k[3]                                                                                                                                                     + a1514*k[15])
      k[17] = Δt*f(t + c[16]*Δt,u + a1600*k[1] + a1601*k[2] + a1602*k[3]              + a1604*k[5] + a1605*k[6] + a1606*k[7] + a1607*k[8] + a1608*k[9] + a1609*k[10] + a1610*k[11] + a1611*k[12] + a1612*k[13] + a1613*k[14] + a1614*k[15] + a1615*k[16])
      update = (b[1]*k[1] + b[2]*k[2] + b[3]*k[3] + b[5]*k[5]) + (b[7]*k[7] + b[9]*k[9] + b[10]*k[10] + b[11]*k[11]) + (b[12]*k[12] + b[13]*k[13] + b[14]*k[14] + b[15]*k[15]) + (b[16]*k[16] + b[17]*k[17])
      if adaptive
        utmp = u + update
        EEst = norm(((k[2] - k[16]) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else
        u = u + update
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Feagin12Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,b,c = constructFeagin12(uEltypeNoUnits)
  k = Vector{uType}(0)
  for i = 1:25
    push!(k,similar(u))
  end
  update = similar(u)
  utmp = similar(u)
  uidx = eachindex(u)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp); k[1]=Δt*rtmp
      f(t + c[1]*Δt,u + a0100*k[1],rtmp); k[2]=Δt*rtmp
      f(t + c[2]*Δt ,u + a0200*k[1] + a0201*k[2],rtmp); k[3]=Δt*rtmp
      f(t + c[3]*Δt,u + a0300*k[1]              + a0302*k[3],rtmp); k[4]=Δt*rtmp
      f(t + c[4]*Δt,u + a0400*k[1]              + a0402*k[3] + a0403*k[4],rtmp); k[5]=Δt*rtmp
      f(t + c[5]*Δt,u + a0500*k[1]                           + a0503*k[4] + a0504*k[5],rtmp); k[6]=Δt*rtmp
      f(t + c[6]*Δt,u + a0600*k[1]                           + a0603*k[4] + a0604*k[5] + a0605*k[6],rtmp); k[7]=Δt*rtmp
      f(t + c[7]*Δt,u + a0700*k[1]                                        + a0704*k[5] + a0705*k[6] + a0706*k[7],rtmp); k[8]=Δt*rtmp
      f(t + c[8]*Δt,u + a0800*k[1]                                                     + a0805*k[6] + a0806*k[7] + a0807*k[8],rtmp); k[9]=Δt*rtmp
      f(t + c[9]*Δt,u + a0900*k[1]                                                     + a0905*k[6] + a0906*k[7] + a0907*k[8] + a0908*k[9],rtmp); k[10]=Δt*rtmp
      f(t + c[10]*Δt,u + a1000*k[1]                                                     + a1005*k[6] + a1006*k[7] + a1007*k[8] + a1008*k[9] + a1009*k[10],rtmp); k[11]=Δt*rtmp
      f(t + c[11]*Δt,u + a1100*k[1]                                                     + a1105*k[6] + a1106*k[7] + a1107*k[8] + a1108*k[9] + a1109*k[10] + a1110*k[11],rtmp); k[12]=Δt*rtmp
      f(t + c[12]*Δt,u + a1200*k[1]                                                                                            + a1208*k[9] + a1209*k[10] + a1210*k[11] + a1211*k[12],rtmp); k[13]=Δt*rtmp
      f(t + c[13]*Δt,u + a1300*k[1]                                                                                            + a1308*k[9] + a1309*k[10] + a1310*k[11] + a1311*k[12] + a1312*k[13],rtmp); k[14]=Δt*rtmp
      f(t + c[14]*Δt,u + a1400*k[1]                                                                                            + a1408*k[9] + a1409*k[10] + a1410*k[11] + a1411*k[12] + a1412*k[13] + a1413*k[14],rtmp); k[15]=Δt*rtmp
      f(t + c[15]*Δt,u + a1500*k[1]                                                                                            + a1508*k[9] + a1509*k[10] + a1510*k[11] + a1511*k[12] + a1512*k[13] + a1513*k[14] + a1514*k[15],rtmp); k[16]=Δt*rtmp
      f(t + c[16]*Δt,u + a1600*k[1]                                                                                            + a1608*k[9] + a1609*k[10] + a1610*k[11] + a1611*k[12] + a1612*k[13] + a1613*k[14] + a1614*k[15] + a1615*k[16],rtmp); k[17]=Δt*rtmp
      f(t + c[17]*Δt,u + a1700*k[1]                                                     + a1705*k[6] + a1706*k[7] + a1707*k[8] + a1708*k[9] + a1709*k[10] + a1710*k[11] + a1711*k[12] + a1712*k[13] + a1713*k[14] + a1714*k[15] + a1715*k[16] + a1716*k[17],rtmp); k[18]=Δt*rtmp
      f(t + c[18]*Δt,u + a1800*k[1]                                                     + a1805*k[6] + a1806*k[7] + a1807*k[8] + a1808*k[9] + a1809*k[10] + a1810*k[11] + a1811*k[12] + a1812*k[13] + a1813*k[14] + a1814*k[15] + a1815*k[16] + a1816*k[17] + a1817*k[18],rtmp); k[19]=Δt*rtmp
      f(t + c[19]*Δt,u + a1900*k[1]                                        + a1904*k[5] + a1905*k[6] + a1906*k[7]              + a1908*k[9] + a1909*k[10] + a1910*k[11] + a1911*k[12] + a1912*k[13] + a1913*k[14] + a1914*k[15] + a1915*k[16] + a1916*k[17] + a1917*k[18] + a1918*k[19],rtmp); k[20]=Δt*rtmp
      f(t + c[20]*Δt,u + a2000*k[1]                           + a2003*k[4] + a2004*k[5] + a2005*k[6]              + a2007*k[8]              + a2009*k[10] + a2010*k[11]                                                                                     + a2017*k[18] + a2018*k[19] + a2019*k[20],rtmp); k[21]=Δt*rtmp
      f(t + c[21]*Δt,u + a2100*k[1]              + a2102*k[3] + a2103*k[4]                           + a2106*k[7] + a2107*k[8]              + a2109*k[10] + a2110*k[11]                                                                                     + a2117*k[18] + a2118*k[19] + a2119*k[20] + a2120*k[21],rtmp); k[22]=Δt*rtmp
      f(t + c[22]*Δt,u + a2200*k[1] + a2201*k[2]                           + a2204*k[5]              + a2206*k[7]                                                                                                                                                                                     + a2220*k[21] + a2221*k[22],rtmp); k[23]=Δt*rtmp
      f(t + c[23]*Δt,u + a2300*k[1]              + a2302*k[3]                                                                                                                                                                                                                                                                     + a2322*k[23],rtmp); k[24]=Δt*rtmp
      f(t + c[24]*Δt,u + a2400*k[1] + a2401*k[2] + a2402*k[3]              + a2404*k[5]              + a2406*k[7] + a2407*k[8] + a2408*k[9] + a2409*k[10] + a2410*k[11] + a2411*k[12] + a2412*k[13] + a2413*k[14] + a2414*k[15] + a2415*k[16] + a2416*k[17] + a2417*k[18] + a2418*k[19] + a2419*k[20] + a2420*k[21] + a2421*k[22] + a2422*k[23] + a2423*k[24],rtmp); k[25]=Δt*rtmp

      for i in uidx
        update[i] = (b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i]) + (b[7]*k[7][i] + b[8]*k[8][i] + b[10]*k[10][i] + b[11]*k[11][i]) + (b[13]*k[13][i] + b[14]*k[14][i] + b[15]*k[15][i] + b[16]*k[16][i]) + (b[17]*k[17][i] + b[18]*k[18][i] + b[19]*k[19][i] + b[20]*k[20][i]) + (b[21]*k[21][i] + b[22]*k[22][i] + b[23]*k[23][i] + b[24]*k[24][i]) + b[25]*k[25][i]
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + update[i]
        end
        EEst = norm(((k[2] - k[24]) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting so in-place
        for i in uidx
          u[i] = u[i] + update[i]
        end
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Feagin12,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,b,c = constructFeagin12(uEltypeNoUnits)
  k = Vector{typeof(u)}(0)
  for i = 1:25
    push!(k,similar(u))
  end
  update = similar(u)
  utmp = similar(u); rtmp = rateType(sizeu)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  uidx = eachindex(u)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp)
      for i in uidx
        k[1][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0100*k[1][i]
      end
      f(t + c[1]*Δt,tmp,rtmp)
      for i in uidx
        k[2][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0200*k[1][i] + a0201*k[2][i]
      end
      f(t + c[2]*Δt ,tmp,rtmp)
      for i in uidx
        k[3][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0300*k[1][i] + a0302*k[3][i]
      end
      f(t + c[3]*Δt,tmp,rtmp)
      for i in uidx
        k[4][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0400*k[1][i] + a0402*k[3][i] + a0403*k[4][i]
      end
      f(t + c[4]*Δt,tmp,rtmp)
      for i in uidx
        k[5][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0500*k[1][i] + a0503*k[4][i] + a0504*k[5][i]
      end
      f(t + c[5]*Δt,tmp,rtmp)
      for i in uidx
        k[6][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0600*k[1][i] + a0603*k[4][i] + a0604*k[5][i]) + a0605*k[6][i]
      end
      f(t + c[6]*Δt,tmp,rtmp)
      for i in uidx
        k[7][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0700*k[1][i] + a0704*k[5][i] + a0705*k[6][i]) + a0706*k[7][i]
      end
      f(t + c[7]*Δt,tmp,rtmp)
      for i in uidx
        k[8][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0800*k[1][i] + a0805*k[6][i] + a0806*k[7][i]) + a0807*k[8][i]
      end
      f(t + c[8]*Δt,tmp,rtmp)
      for i in uidx
        k[9][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0900*k[1][i] + a0905*k[6][i] + a0906*k[7][i]) + (a0907*k[8][i] + a0908*k[9][i])
      end
      f(t + c[9]*Δt,tmp,rtmp)
      for i in uidx
        k[10][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1000*k[1][i] + a1005*k[6][i] + a1006*k[7][i]) + (a1007*k[8][i] + a1008*k[9][i] + a1009*k[10][i])
      end
      f(t + c[10]*Δt,tmp,rtmp)
      for i in uidx
        k[11][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1100*k[1][i] + a1105*k[6][i] + a1106*k[7][i]) + (a1107*k[8][i] + a1108*k[9][i] + a1109*k[10][i] + a1110*k[11][i])
      end
      f(t + c[11]*Δt,tmp,rtmp)
      for i in uidx
        k[12][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1200*k[1][i] + a1208*k[9][i] + a1209*k[10][i]) + (a1210*k[11][i] + a1211*k[12][i])
      end
      f(t + c[12]*Δt,tmp,rtmp)
      for i in uidx
        k[13][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1300*k[1][i] + a1308*k[9][i] + a1309*k[10][i]) + (a1310*k[11][i] + a1311*k[12][i] + a1312*k[13][i])
      end
      f(t + c[13]*Δt,tmp,rtmp)
      for i in uidx
        k[14][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1400*k[1][i] + a1408*k[9][i] + a1409*k[10][i]) + (a1410*k[11][i] + a1411*k[12][i] + a1412*k[13][i] + a1413*k[14][i])
      end
      f(t + c[14]*Δt,tmp,rtmp)
      for i in uidx
        k[15][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1500*k[1][i] + a1508*k[9][i] + a1509*k[10][i]) + (a1510*k[11][i] + a1511*k[12][i] + a1512*k[13][i] + a1513*k[14][i]) + a1514*k[15][i]
      end
      f(t + c[15]*Δt,tmp,rtmp)
      for i in uidx
        k[16][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1600*k[1][i] + a1608*k[9][i] + a1609*k[10][i]) + (a1610*k[11][i] + a1611*k[12][i] + a1612*k[13][i] + a1613*k[14][i]) + (a1614*k[15][i] + a1615*k[16][i])
      end
      f(t + c[16]*Δt,tmp,rtmp)
      for i in uidx
        k[17][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1700*k[1][i] + a1705*k[6][i] + a1706*k[7][i]) + (a1707*k[8][i] + a1708*k[9][i] + a1709*k[10][i] + a1710*k[11][i]) + (a1711*k[12][i] + a1712*k[13][i] + a1713*k[14][i] + a1714*k[15][i]) + (a1715*k[16][i] + a1716*k[17][i])
      end
      f(t + c[17]*Δt,tmp,rtmp)
      for i in uidx
        k[18][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1800*k[1][i] + a1805*k[6][i] + a1806*k[7][i]) + (a1807*k[8][i] + a1808*k[9][i] + a1809*k[10][i] + a1810*k[11][i]) + (a1811*k[12][i] + a1812*k[13][i] + a1813*k[14][i] + a1814*k[15][i]) + (a1815*k[16][i] + a1816*k[17][i] + a1817*k[18][i])
      end
      f(t + c[18]*Δt,tmp,rtmp)
      for i in uidx
        k[19][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1900*k[1][i] + a1904*k[5][i] + a1905*k[6][i]) + (a1906*k[7][i] + a1908*k[9][i] + a1909*k[10][i] + a1910*k[11][i]) + (a1911*k[12][i] + a1912*k[13][i] + a1913*k[14][i] + a1914*k[15][i]) + (a1915*k[16][i] + a1916*k[17][i] + a1917*k[18][i] + a1918*k[19][i])
      end
      f(t + c[19]*Δt,tmp,rtmp)
      for i in uidx
        k[20][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2000*k[1][i] + a2003*k[4][i] + a2004*k[5][i]) + (a2005*k[6][i] + a2007*k[8][i] + a2009*k[10][i] + a2010*k[11][i]) + (a2017*k[18][i] + a2018*k[19][i] + a2019*k[20][i])
      end
      f(t + c[20]*Δt,tmp,rtmp)
      for i in uidx
        k[21][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2100*k[1][i] + a2102*k[3][i] + a2103*k[4][i]) + (a2106*k[7][i] + a2107*k[8][i] + a2109*k[10][i] + a2110*k[11][i]) + (a2117*k[18][i] + a2118*k[19][i] + a2119*k[20][i] + a2120*k[21][i])
      end
      f(t + c[21]*Δt,tmp,rtmp)
      for i in uidx
        k[22][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2200*k[1][i] + a2201*k[2][i] + a2204*k[5][i]) + (a2206*k[7][i] + a2220*k[21][i] + a2221*k[22][i])
      end
      f(t + c[22]*Δt,tmp,rtmp)
      for i in uidx
        k[23][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a2300*k[1][i] + a2302*k[3][i] + a2322*k[23][i]
      end
      f(t + c[23]*Δt,tmp,rtmp)
      for i in uidx
        k[24][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2400*k[1][i] + a2401*k[2][i] + a2402*k[3][i]) + (a2404*k[5][i] + a2406*k[7][i] + a2407*k[8][i] + a2408*k[9][i]) + (a2409*k[10][i] + a2410*k[11][i] + a2411*k[12][i] + a2412*k[13][i]) + (a2413*k[14][i] + a2414*k[15][i] + a2415*k[16][i] + a2416*k[17][i]) + (a2417*k[18][i] + a2418*k[19][i] + a2419*k[20][i] + a2420*k[21][i]) + (a2421*k[22][i] + a2422*k[23][i] + a2423*k[24][i])
      end
      f(t + c[24]*Δt,tmp,rtmp)
      for i in uidx
        k[25][i] = Δt*rtmp[i]
        update[i] = (b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i]) + (b[7]*k[7][i] + b[8]*k[8][i] + b[10]*k[10][i] + b[11]*k[11][i]) + (b[13]*k[13][i] + b[14]*k[14][i] + b[15]*k[15][i] + b[16]*k[16][i]) + (b[17]*k[17][i] + b[18]*k[18][i] + b[19]*k[19][i] + b[20]*k[20][i]) + (b[21]*k[21][i] + b[22]*k[22][i] + b[23]*k[23][i] + b[24]*k[24][i]) + b[25]*k[25][i]
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + update[i]
          atmp[i] = ((k[2][i] - k[24][i]) * adaptiveConst)/(abstol+u[i]*reltol)
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
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:Feagin12,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1705,a1706,a1707,a1708,a1709,a1710,a1711,a1712,a1713,a1714,a1715,a1716,a1800,a1805,a1806,a1807,a1808,a1809,a1810,a1811,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1904,a1905,a1906,a1908,a1909,a1910,a1911,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2003,a2004,a2005,a2007,a2009,a2010,a2017,a2018,a2019,a2100,a2102,a2103,a2106,a2107,a2109,a2110,a2117,a2118,a2119,a2120,a2200,a2201,a2204,a2206,a2220,a2221,a2300,a2302,a2322,a2400,a2401,a2402,a2404,a2406,a2407,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,b,c = constructFeagin12(uEltypeNoUnits)
  k = Vector{typeof(u)}(25)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k[1]  = Δt*f(t,u)
      k[2]  = Δt*f(t + c[1]*Δt,u + a0100*k[1])
      k[3]  = Δt*f(t + c[2]*Δt ,u + a0200*k[1] + a0201*k[2])
      k[4]  = Δt*f(t + c[3]*Δt,u + a0300*k[1]              + a0302*k[3])
      k[5]  = Δt*f(t + c[4]*Δt,u + a0400*k[1]              + a0402*k[3] + a0403*k[4])
      k[6]  = Δt*f(t + c[5]*Δt,u + a0500*k[1]                           + a0503*k[4] + a0504*k[5])
      k[7]  = Δt*f(t + c[6]*Δt,u + a0600*k[1]                           + a0603*k[4] + a0604*k[5] + a0605*k[6])
      k[8]  = Δt*f(t + c[7]*Δt,u + a0700*k[1]                                        + a0704*k[5] + a0705*k[6] + a0706*k[7])
      k[9]  = Δt*f(t + c[8]*Δt,u + a0800*k[1]                                                     + a0805*k[6] + a0806*k[7] + a0807*k[8])
      k[10] = Δt*f(t + c[9]*Δt,u + a0900*k[1]                                                     + a0905*k[6] + a0906*k[7] + a0907*k[8] + a0908*k[9])
      k[11] = Δt*f(t + c[10]*Δt,u + a1000*k[1]                                                     + a1005*k[6] + a1006*k[7] + a1007*k[8] + a1008*k[9] + a1009*k[10])
      k[12] = Δt*f(t + c[11]*Δt,u + a1100*k[1]                                                     + a1105*k[6] + a1106*k[7] + a1107*k[8] + a1108*k[9] + a1109*k[10] + a1110*k[11])
      k[13] = Δt*f(t + c[12]*Δt,u + a1200*k[1]                                                                                            + a1208*k[9] + a1209*k[10] + a1210*k[11] + a1211*k[12])
      k[14] = Δt*f(t + c[13]*Δt,u + a1300*k[1]                                                                                            + a1308*k[9] + a1309*k[10] + a1310*k[11] + a1311*k[12] + a1312*k[13])
      k[15] = Δt*f(t + c[14]*Δt,u + a1400*k[1]                                                                                            + a1408*k[9] + a1409*k[10] + a1410*k[11] + a1411*k[12] + a1412*k[13] + a1413*k[14])
      k[16] = Δt*f(t + c[15]*Δt,u + a1500*k[1]                                                                                            + a1508*k[9] + a1509*k[10] + a1510*k[11] + a1511*k[12] + a1512*k[13] + a1513*k[14] + a1514*k[15])
      k[17] = Δt*f(t + c[16]*Δt,u + a1600*k[1]                                                                                            + a1608*k[9] + a1609*k[10] + a1610*k[11] + a1611*k[12] + a1612*k[13] + a1613*k[14] + a1614*k[15] + a1615*k[16])
      k[18] = Δt*f(t + c[17]*Δt,u + a1700*k[1]                                                     + a1705*k[6] + a1706*k[7] + a1707*k[8] + a1708*k[9] + a1709*k[10] + a1710*k[11] + a1711*k[12] + a1712*k[13] + a1713*k[14] + a1714*k[15] + a1715*k[16] + a1716*k[17])
      k[19] = Δt*f(t + c[18]*Δt,u + a1800*k[1]                                                     + a1805*k[6] + a1806*k[7] + a1807*k[8] + a1808*k[9] + a1809*k[10] + a1810*k[11] + a1811*k[12] + a1812*k[13] + a1813*k[14] + a1814*k[15] + a1815*k[16] + a1816*k[17] + a1817*k[18])
      k[20] = Δt*f(t + c[19]*Δt,u + a1900*k[1]                                        + a1904*k[5] + a1905*k[6] + a1906*k[7]              + a1908*k[9] + a1909*k[10] + a1910*k[11] + a1911*k[12] + a1912*k[13] + a1913*k[14] + a1914*k[15] + a1915*k[16] + a1916*k[17] + a1917*k[18] + a1918*k[19])
      k[21] = Δt*f(t + c[20]*Δt,u + a2000*k[1]                           + a2003*k[4] + a2004*k[5] + a2005*k[6]              + a2007*k[8]              + a2009*k[10] + a2010*k[11]                                                                                     + a2017*k[18] + a2018*k[19] + a2019*k[20])
      k[22] = Δt*f(t + c[21]*Δt,u + a2100*k[1]              + a2102*k[3] + a2103*k[4]                           + a2106*k[7] + a2107*k[8]              + a2109*k[10] + a2110*k[11]                                                                                     + a2117*k[18] + a2118*k[19] + a2119*k[20] + a2120*k[21])
      k[23] = Δt*f(t + c[22]*Δt,u + a2200*k[1] + a2201*k[2]                           + a2204*k[5]              + a2206*k[7]                                                                                                                                                                                     + a2220*k[21] + a2221*k[22])
      k[24] = Δt*f(t + c[23]*Δt,u + a2300*k[1]              + a2302*k[3]                                                                                                                                                                                                                                                                     + a2322*k[23])
      k[25] = Δt*f(t + c[24]*Δt,u + a2400*k[1] + a2401*k[2] + a2402*k[3]              + a2404*k[5]              + a2406*k[7] + a2407*k[8] + a2408*k[9] + a2409*k[10] + a2410*k[11] + a2411*k[12] + a2412*k[13] + a2413*k[14] + a2414*k[15] + a2415*k[16] + a2416*k[17] + a2417*k[18] + a2418*k[19] + a2419*k[20] + a2420*k[21] + a2421*k[22] + a2422*k[23] + a2423*k[24])

      update = (b[1]*k[1] + b[2]*k[2] + b[3]*k[3] + b[5]*k[5]) + (b[7]*k[7] + b[8]*k[8] + b[10]*k[10] + b[11]*k[11]) + (b[13]*k[13] + b[14]*k[14] + b[15]*k[15] + b[16]*k[16]) + (b[17]*k[17] + b[18]*k[18] + b[19]*k[19] + b[20]*k[20]) + (b[21]*k[21] + b[22]*k[22] + b[23]*k[23] + b[24]*k[24]) + (b[25]*k[25])
      if adaptive
        utmp = u + update
        EEst = norm(((k[2] - k[24]) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting so in-place
        u = u + update
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Feagin14,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,b,c = constructFeagin14(uEltypeNoUnits)
  k = Vector{typeof(u)}(0)
  for i = 1:35
    push!(k,similar(u))
  end
  update = similar(u)
  utmp = similar(u);
  rtmp = rateType(sizeu)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  uidx = eachindex(u)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp)
      for i in uidx
        k[1][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0100*k[1][i]
      end
      f(t + c[1]*Δt,tmp,rtmp)
      for i in uidx
        k[2][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0200*k[1][i] + a0201*k[2][i]
      end
      f(t + c[2]*Δt ,tmp,rtmp)
      for i in uidx
        k[3][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0300*k[1][i] + a0302*k[3][i]
      end
      f(t + c[3]*Δt,tmp,rtmp)
      for i in uidx
        k[4][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0400*k[1][i] + a0402*k[3][i] + a0403*k[4][i]
      end
      f(t + c[4]*Δt,tmp,rtmp)
      for i in uidx
        k[5][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a0500*k[1][i] + a0503*k[4][i] + a0504*k[5][i]
      end
      f(t + c[5]*Δt,tmp,rtmp)
      for i in uidx
        k[6][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0600*k[1][i] + a0603*k[4][i] + a0604*k[5][i]) + a0605*k[6][i]
      end
      f(t + c[6]*Δt,tmp,rtmp)
      for i in uidx
        k[7][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0700*k[1][i] + a0704*k[5][i] + a0705*k[6][i]) + a0706*k[7][i]
      end
      f(t + c[7]*Δt,tmp,rtmp)
      for i in uidx
        k[8][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0800*k[1][i] + a0805*k[6][i] + a0806*k[7][i]) + a0807*k[8][i]
      end
      f(t + c[8]*Δt,tmp,rtmp)
      for i in uidx
        k[9][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a0900*k[1][i] + a0905*k[6][i] + a0906*k[7][i]) + a0907*k[8][i] + a0908*k[9][i]
      end
      f(t + c[9]*Δt,tmp,rtmp)
      for i in uidx
        k[10][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1000*k[1][i] + a1005*k[6][i] + a1006*k[7][i]) + (a1007*k[8][i] + a1008*k[9][i] + a1009*k[10][i])
      end
      f(t + c[10]*Δt,tmp,rtmp)
      for i in uidx
        k[11][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1100*k[1][i] + a1105*k[6][i] + a1106*k[7][i]) + (a1107*k[8][i] + a1108*k[9][i] + a1109*k[10][i] + a1110*k[11][i])
      end
      f(t + c[11]*Δt,tmp,rtmp)
      for i in uidx
        k[12][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1200*k[1][i] + a1208*k[9][i] + a1209*k[10][i]) + (a1210*k[11][i] + a1211*k[12][i])
      end
      f(t + c[12]*Δt,tmp,rtmp)
      for i in uidx
        k[13][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1300*k[1][i] + a1308*k[9][i] + a1309*k[10][i]) + (a1310*k[11][i] + a1311*k[12][i] + a1312*k[13][i])
      end
      f(t + c[13]*Δt,tmp,rtmp)
      for i in uidx
        k[14][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1400*k[1][i] + a1408*k[9][i] + a1409*k[10][i]) + (a1410*k[11][i] + a1411*k[12][i] + a1412*k[13][i] + a1413*k[14][i])
      end
      f(t + c[14]*Δt,tmp,rtmp)
      for i in uidx
        k[15][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1500*k[1][i] + a1508*k[9][i] + a1509*k[10][i]) + (a1510*k[11][i] + a1511*k[12][i] + a1512*k[13][i] + a1513*k[14][i]) + a1514*k[15][i]
      end
      f(t + c[15]*Δt,tmp,rtmp)
      for i in uidx
        k[16][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1600*k[1][i] + a1608*k[9][i] + a1609*k[10][i]) + (a1610*k[11][i] + a1611*k[12][i] + a1612*k[13][i] + a1613*k[14][i]) + a1614*k[15][i] + a1615*k[16][i]
      end
      f(t + c[16]*Δt,tmp,rtmp)
      for i in uidx
        k[17][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1700*k[1][i] + a1712*k[13][i] + a1713*k[14][i]) + (a1714*k[15][i] + a1715*k[16][i] + a1716*k[17][i])
      end
      f(t + c[17]*Δt,tmp,rtmp)
      for i in uidx
        k[18][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1800*k[1][i] + a1812*k[13][i] + a1813*k[14][i]) + (a1814*k[15][i] + a1815*k[16][i] + a1816*k[17][i] + a1817*k[18][i])
      end
      f(t + c[18]*Δt,tmp,rtmp)
      for i in uidx
        k[19][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a1900*k[1][i] + a1912*k[13][i] + a1913*k[14][i]) + (a1914*k[15][i] + a1915*k[16][i] + a1916*k[17][i] + a1917*k[18][i]) + a1918*k[19][i]
      end
      f(t + c[19]*Δt,tmp,rtmp)
      for i in uidx
        k[20][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2000*k[1][i] + a2012*k[13][i] + a2013*k[14][i]) + (a2014*k[15][i] + a2015*k[16][i] + a2016*k[17][i] + a2017*k[18][i]) + (a2018*k[19][i] + a2019*k[20][i])
      end
      f(t + c[20]*Δt,tmp,rtmp)
      for i in uidx
        k[21][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2100*k[1][i] + a2112*k[13][i] + a2113*k[14][i]) + (a2114*k[15][i] + a2115*k[16][i] + a2116*k[17][i] + a2117*k[18][i]) + (a2118*k[19][i] + a2119*k[20][i] + a2120*k[21][i])
      end
      f(t + c[21]*Δt,tmp,rtmp)
      for i in uidx
        k[22][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2200*k[1][i] + a2212*k[13][i] + a2213*k[14][i]) + (a2214*k[15][i] + a2215*k[16][i] + a2216*k[17][i] + a2217*k[18][i]) + (a2218*k[19][i] + a2219*k[20][i] + a2220*k[21][i] + a2221*k[22][i])
      end
      f(t + c[22]*Δt,tmp,rtmp)
      for i in uidx
        k[23][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2300*k[1][i] + a2308*k[9][i] + a2309*k[10][i]) + (a2310*k[11][i] + a2311*k[12][i] + a2312*k[13][i] + a2313*k[14][i]) + (a2314*k[15][i] + a2315*k[16][i] + a2316*k[17][i] + a2317*k[18][i]) + (a2318*k[19][i] + a2319*k[20][i] + a2320*k[21][i] + a2321*k[22][i]) + (a2322*k[23][i])
      end
      f(t + c[23]*Δt,tmp,rtmp)
      for i in uidx
        k[24][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2400*k[1][i] + a2408*k[9][i] + a2409*k[10][i]) + (a2410*k[11][i] + a2411*k[12][i] + a2412*k[13][i] + a2413*k[14][i]) + (a2414*k[15][i] + a2415*k[16][i] + a2416*k[17][i] + a2417*k[18][i]) + (a2418*k[19][i] + a2419*k[20][i] + a2420*k[21][i] + a2421*k[22][i]) + (a2422*k[23][i] + a2423*k[24][i])
      end
      f(t + c[24]*Δt,tmp,rtmp)
      for i in uidx
        k[25][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2500*k[1][i] + a2508*k[9][i] + a2509*k[10][i]) + (a2510*k[11][i] + a2511*k[12][i] + a2512*k[13][i] + a2513*k[14][i]) + (a2514*k[15][i] + a2515*k[16][i] + a2516*k[17][i] + a2517*k[18][i]) + (a2518*k[19][i] + a2519*k[20][i] + a2520*k[21][i] + a2521*k[22][i]) + (a2522*k[23][i] + a2523*k[24][i] + a2524*k[25][i])
      end
      f(t + c[25]*Δt,tmp,rtmp)
      for i in uidx
        k[26][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2600*k[1][i] + a2605*k[6][i] + a2606*k[7][i]) + (a2607*k[8][i] + a2608*k[9][i] + a2609*k[10][i] + a2610*k[11][i]) + (a2612*k[13][i] + a2613*k[14][i] + a2614*k[15][i] + a2615*k[16][i]) + (a2616*k[17][i] + a2617*k[18][i] + a2618*k[19][i] + a2619*k[20][i]) + (a2620*k[21][i] + a2621*k[22][i] + a2622*k[23][i] + a2623*k[24][i]) + (a2624*k[25][i] + a2625*k[26][i])
      end
      f(t + c[26]*Δt,tmp,rtmp)
      for i in uidx
        k[27][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2700*k[1][i] + a2705*k[6][i] + a2706*k[7][i]) + (a2707*k[8][i] + a2708*k[9][i] + a2709*k[10][i] + a2711*k[12][i]) + (a2712*k[13][i] + a2713*k[14][i] + a2714*k[15][i] + a2715*k[16][i]) + (a2716*k[17][i] + a2717*k[18][i] + a2718*k[19][i] + a2719*k[20][i]) + (a2720*k[21][i] + a2721*k[22][i] + a2722*k[23][i] + a2723*k[24][i]) + (a2724*k[25][i] + a2725*k[26][i] + a2726*k[27][i])
      end
      f(t + c[27]*Δt,tmp,rtmp)
      for i in uidx
        k[28][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2800*k[1][i] + a2805*k[6][i] + a2806*k[7][i]) + (a2807*k[8][i] + a2808*k[9][i] + a2810*k[11][i] + a2811*k[12][i]) + (a2813*k[14][i] + a2814*k[15][i] + a2815*k[16][i] + a2823*k[24][i]) + (a2824*k[25][i] + a2825*k[26][i] + a2826*k[27][i] + a2827*k[28][i])
      end
      f(t + c[28]*Δt,tmp,rtmp)
      for i in uidx
        k[29][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a2900*k[1][i] + a2904*k[5][i] + a2905*k[6][i]) + (a2906*k[7][i] + a2909*k[10][i] + a2910*k[11][i] + a2911*k[12][i]) + (a2913*k[14][i] + a2914*k[15][i] + a2915*k[16][i] + a2923*k[24][i]) + (a2924*k[25][i] + a2925*k[26][i] + a2926*k[27][i] + a2927*k[28][i]) + (a2928*k[29][i])
      end
      f(t + c[29]*Δt,tmp,rtmp)
      for i in uidx
        k[30][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a3000*k[1][i] + a3003*k[4][i] + a3004*k[5][i]) + (a3005*k[6][i] + a3007*k[8][i] + a3009*k[10][i] + a3010*k[11][i]) + (a3013*k[14][i] + a3014*k[15][i] + a3015*k[16][i] + a3023*k[24][i]) + (a3024*k[25][i] + a3025*k[26][i] + a3027*k[28][i] + a3028*k[29][i]) + (a3029*k[30][i])
      end
      f(t + c[30]*Δt,tmp,rtmp)
      for i in uidx
        k[31][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a3100*k[1][i] + a3102*k[3][i] + a3103*k[4][i]) + (a3106*k[7][i] + a3107*k[8][i] + a3109*k[10][i] + a3110*k[11][i]) + (a3113*k[14][i] + a3114*k[15][i] + a3115*k[16][i] + a3123*k[24][i]) + (a3124*k[25][i] + a3125*k[26][i] + a3127*k[28][i] + a3128*k[29][i]) + (a3129*k[30][i] + a3130*k[31][i])
      end
      f(t + c[31]*Δt,tmp,rtmp)
      for i in uidx
        k[32][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a3200*k[1][i] + a3201*k[2][i] + a3204*k[5][i]) + (a3206*k[7][i] + a3230*k[31][i] + a3231*k[32][i])
      end
      f(t + c[32]*Δt,tmp,rtmp)
      for i in uidx
        k[33][i]=Δt*rtmp[i]
        tmp[i] = u[i] + a3300*k[1][i] + a3302*k[3][i] + a3332*k[33][i]
      end
      f(t + c[33]*Δt,tmp,rtmp)
      for i in uidx
        k[34][i]=Δt*rtmp[i]
        tmp[i] = (u[i] + a3400*k[1][i] + a3401*k[2][i] + a3402*k[3][i]) + (a3404*k[5][i] + a3406*k[7][i] + a3407*k[8][i] + a3409*k[10][i]) + (a3410*k[11][i] + a3411*k[12][i] + a3412*k[13][i] + a3413*k[14][i]) + (a3414*k[15][i] + a3415*k[16][i] + a3416*k[17][i] + a3417*k[18][i]) + (a3418*k[19][i] + a3419*k[20][i] + a3420*k[21][i] + a3421*k[22][i]) + (a3422*k[23][i] + a3423*k[24][i] + a3424*k[25][i] + a3425*k[26][i]) + (a3426*k[27][i] + a3427*k[28][i] + a3428*k[29][i] + a3429*k[30][i]) + (a3430*k[31][i] + a3431*k[32][i] + a3432*k[33][i] + a3433*k[34][i])
      end
      f(t + c[34]*Δt,tmp,rtmp)
      for i in uidx
        k[35][i]=Δt*rtmp[i]
        update[i] = (b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i]) + (b[7]*k[7][i] + b[8]*k[8][i] + b[10]*k[10][i] + b[11]*k[11][i]) + (b[12]*k[12][i] + b[14]*k[14][i] + b[15]*k[15][i] + b[16]*k[16][i]) + (b[18]*k[18][i] + b[19]*k[19][i] + b[20]*k[20][i] + b[21]*k[21][i]) + (b[22]*k[22][i] + b[23]*k[23][i] + b[24]*k[24][i] + b[25]*k[25][i]) + (b[26]*k[26][i] + b[27]*k[27][i] + b[28]*k[28][i] + b[29]*k[29][i]) + (b[30]*k[30][i] + b[31]*k[31][i] + b[32]*k[32][i] + b[33]*k[33][i]) + (b[34]*k[34][i] + b[35]*k[35][i])
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + update[i]
          atmp[i] = ((k[2][i] - k[34][i]) * adaptiveConst)./(abstol+u[i]*reltol)
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
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Feagin14Vectorized,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,b,c = constructFeagin14(uEltypeNoUnits)
  k = Vector{typeof(u)}(0)
  for i = 1:35
    push!(k,similar(u))
  end
  update = similar(u)
  utmp = similar(u)
  uidx = eachindex(u)
  rtmp = rateType(sizeu)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      f(t,u,rtmp); k[1]=Δt*rtmp
      f(t + c[1]*Δt,u + a0100*k[1],rtmp); k[2]=Δt*rtmp
      f(t + c[2]*Δt ,u + a0200*k[1] + a0201*k[2],rtmp); k[3]=Δt*rtmp
      f(t + c[3]*Δt,u + a0300*k[1]              + a0302*k[3],rtmp); k[4]=Δt*rtmp
      f(t + c[4]*Δt,u + a0400*k[1]              + a0402*k[3] + a0403*k[4],rtmp); k[5]=Δt*rtmp
      f(t + c[5]*Δt,u + a0500*k[1]                           + a0503*k[4] + a0504*k[5],rtmp); k[6]=Δt*rtmp
      f(t + c[6]*Δt,u + a0600*k[1]                           + a0603*k[4] + a0604*k[5] + a0605*k[6],rtmp); k[7]=Δt*rtmp
      f(t + c[7]*Δt,u + a0700*k[1]                                        + a0704*k[5] + a0705*k[6] + a0706*k[7],rtmp); k[8]=Δt*rtmp
      f(t + c[8]*Δt,u + a0800*k[1]                                                     + a0805*k[6] + a0806*k[7] + a0807*k[8],rtmp); k[9]=Δt*rtmp
      f(t + c[9]*Δt,u + a0900*k[1]                                                     + a0905*k[6] + a0906*k[7] + a0907*k[8] + a0908*k[9],rtmp); k[10]=Δt*rtmp
      f(t + c[10]*Δt,u + a1000*k[1]                                                     + a1005*k[6] + a1006*k[7] + a1007*k[8] + a1008*k[9] + a1009*k[10],rtmp); k[11]=Δt*rtmp
      f(t + c[11]*Δt,u + a1100*k[1]                                                     + a1105*k[6] + a1106*k[7] + a1107*k[8] + a1108*k[9] + a1109*k[10] + a1110*k[11],rtmp); k[12]=Δt*rtmp
      f(t + c[12]*Δt,u + a1200*k[1]                                                                                            + a1208*k[9] + a1209*k[10] + a1210*k[11] + a1211*k[12],rtmp); k[13]=Δt*rtmp
      f(t + c[13]*Δt,u + a1300*k[1]                                                                                            + a1308*k[9] + a1309*k[10] + a1310*k[11] + a1311*k[12] + a1312*k[13],rtmp); k[14]=Δt*rtmp
      f(t + c[14]*Δt,u + a1400*k[1]                                                                                            + a1408*k[9] + a1409*k[10] + a1410*k[11] + a1411*k[12] + a1412*k[13] + a1413*k[14],rtmp); k[15]=Δt*rtmp
      f(t + c[15]*Δt,u + a1500*k[1]                                                                                            + a1508*k[9] + a1509*k[10] + a1510*k[11] + a1511*k[12] + a1512*k[13] + a1513*k[14] + a1514*k[15],rtmp); k[16]=Δt*rtmp
      f(t + c[16]*Δt,u + a1600*k[1]                                                                                            + a1608*k[9] + a1609*k[10] + a1610*k[11] + a1611*k[12] + a1612*k[13] + a1613*k[14] + a1614*k[15] + a1615*k[16],rtmp); k[17]=Δt*rtmp
      f(t + c[17]*Δt,u + a1700*k[1]                                                                                                                                                   + a1712*k[13] + a1713*k[14] + a1714*k[15] + a1715*k[16] + a1716*k[17],rtmp); k[18]=Δt*rtmp
      f(t + c[18]*Δt,u + a1800*k[1]                                                                                                                                                   + a1812*k[13] + a1813*k[14] + a1814*k[15] + a1815*k[16] + a1816*k[17] + a1817*k[18],rtmp); k[19]=Δt*rtmp
      f(t + c[19]*Δt,u + a1900*k[1]                                                                                                                                                   + a1912*k[13] + a1913*k[14] + a1914*k[15] + a1915*k[16] + a1916*k[17] + a1917*k[18] + a1918*k[19],rtmp); k[20]=Δt*rtmp
      f(t + c[20]*Δt,u + a2000*k[1]                                                                                                                                                   + a2012*k[13] + a2013*k[14] + a2014*k[15] + a2015*k[16] + a2016*k[17] + a2017*k[18] + a2018*k[19] + a2019*k[20],rtmp); k[21]=Δt*rtmp
      f(t + c[21]*Δt,u + a2100*k[1]                                                                                                                                                   + a2112*k[13] + a2113*k[14] + a2114*k[15] + a2115*k[16] + a2116*k[17] + a2117*k[18] + a2118*k[19] + a2119*k[20] + a2120*k[21],rtmp); k[22]=Δt*rtmp
      f(t + c[22]*Δt,u + a2200*k[1]                                                                                                                                                   + a2212*k[13] + a2213*k[14] + a2214*k[15] + a2215*k[16] + a2216*k[17] + a2217*k[18] + a2218*k[19] + a2219*k[20] + a2220*k[21] + a2221*k[22],rtmp); k[23]=Δt*rtmp
      f(t + c[23]*Δt,u + a2300*k[1]                                                                                            + a2308*k[9] + a2309*k[10] + a2310*k[11] + a2311*k[12] + a2312*k[13] + a2313*k[14] + a2314*k[15] + a2315*k[16] + a2316*k[17] + a2317*k[18] + a2318*k[19] + a2319*k[20] + a2320*k[21] + a2321*k[22] + a2322*k[23],rtmp); k[24]=Δt*rtmp
      f(t + c[24]*Δt,u + a2400*k[1]                                                                                            + a2408*k[9] + a2409*k[10] + a2410*k[11] + a2411*k[12] + a2412*k[13] + a2413*k[14] + a2414*k[15] + a2415*k[16] + a2416*k[17] + a2417*k[18] + a2418*k[19] + a2419*k[20] + a2420*k[21] + a2421*k[22] + a2422*k[23] + a2423*k[24],rtmp); k[25]=Δt*rtmp
      f(t + c[25]*Δt,u + a2500*k[1]                                                                                            + a2508*k[9] + a2509*k[10] + a2510*k[11] + a2511*k[12] + a2512*k[13] + a2513*k[14] + a2514*k[15] + a2515*k[16] + a2516*k[17] + a2517*k[18] + a2518*k[19] + a2519*k[20] + a2520*k[21] + a2521*k[22] + a2522*k[23] + a2523*k[24] + a2524*k[25],rtmp); k[26]=Δt*rtmp
      f(t + c[26]*Δt,u + a2600*k[1]                                                     + a2605*k[6] + a2606*k[7] + a2607*k[8] + a2608*k[9] + a2609*k[10] + a2610*k[11]               + a2612*k[13] + a2613*k[14] + a2614*k[15] + a2615*k[16] + a2616*k[17] + a2617*k[18] + a2618*k[19] + a2619*k[20] + a2620*k[21] + a2621*k[22] + a2622*k[23] + a2623*k[24] + a2624*k[25] + a2625*k[26],rtmp); k[27]=Δt*rtmp
      f(t + c[27]*Δt,u + a2700*k[1]                                                     + a2705*k[6] + a2706*k[7] + a2707*k[8] + a2708*k[9] + a2709*k[10]               + a2711*k[12] + a2712*k[13] + a2713*k[14] + a2714*k[15] + a2715*k[16] + a2716*k[17] + a2717*k[18] + a2718*k[19] + a2719*k[20] + a2720*k[21] + a2721*k[22] + a2722*k[23] + a2723*k[24] + a2724*k[25] + a2725*k[26] + a2726*k[27],rtmp); k[28]=Δt*rtmp
      f(t + c[28]*Δt,u + a2800*k[1]                                                     + a2805*k[6] + a2806*k[7] + a2807*k[8] + a2808*k[9]               + a2810*k[11] + a2811*k[12]               + a2813*k[14] + a2814*k[15] + a2815*k[16]                                                                                                   + a2823*k[24] + a2824*k[25] + a2825*k[26] + a2826*k[27] + a2827*k[28],rtmp); k[29]=Δt*rtmp
      f(t + c[29]*Δt,u + a2900*k[1]                                        + a2904*k[5] + a2905*k[6] + a2906*k[7]                           + a2909*k[10] + a2910*k[11] + a2911*k[12]               + a2913*k[14] + a2914*k[15] + a2915*k[16]                                                                                                   + a2923*k[24] + a2924*k[25] + a2925*k[26] + a2926*k[27] + a2927*k[28] + a2928*k[29],rtmp); k[30]=Δt*rtmp
      f(t + c[30]*Δt,u + a3000*k[1]                           + a3003*k[4] + a3004*k[5] + a3005*k[6]              + a3007*k[8]              + a3009*k[10] + a3010*k[11]                             + a3013*k[14] + a3014*k[15] + a3015*k[16]                                                                                                   + a3023*k[24] + a3024*k[25] + a3025*k[26]               + a3027*k[28] + a3028*k[29] + a3029*k[30],rtmp); k[31]=Δt*rtmp
      f(t + c[31]*Δt,u + a3100*k[1]              + a3102*k[3] + a3103*k[4]                           + a3106*k[7] + a3107*k[8]              + a3109*k[10] + a3110*k[11]                             + a3113*k[14] + a3114*k[15] + a3115*k[16]                                                                                                   + a3123*k[24] + a3124*k[25] + a3125*k[26]               + a3127*k[28] + a3128*k[29] + a3129*k[30] + a3130*k[31],rtmp); k[32]=Δt*rtmp
      f(t + c[32]*Δt,u + a3200*k[1] + a3201*k[2]                           + a3204*k[5]              + a3206*k[7]                                                                                                                                                                                                                                                                                                                                 + a3230*k[31] + a3231*k[32],rtmp); k[33]=Δt*rtmp
      f(t + c[33]*Δt,u + a3300*k[1]              + a3302*k[3]                                                                                                                                                                                                                                                                                                                                                                                                                 + a3332*k[33],rtmp); k[34]=Δt*rtmp
      f(t + c[34]*Δt,u + a3400*k[1] + a3401*k[2] + a3402*k[3]              + a3404*k[5]              + a3406*k[7] + a3407*k[8]              + a3409*k[10] + a3410*k[11] + a3411*k[12] + a3412*k[13] + a3413*k[14] + a3414*k[15] + a3415*k[16] + a3416*k[17] + a3417*k[18] + a3418*k[19] + a3419*k[20] + a3420*k[21] + a3421*k[22] + a3422*k[23] + a3423*k[24] + a3424*k[25] + a3425*k[26] + a3426*k[27] + a3427*k[28] + a3428*k[29] + a3429*k[30] + a3430*k[31] + a3431*k[32] + a3432*k[33] + a3433*k[34],rtmp); k[35]=Δt*rtmp
      for i in uidx
        update[i] = (b[1]*k[1][i] + b[2]*k[2][i] + b[3]*k[3][i] + b[5]*k[5][i]) + (b[7]*k[7][i] + b[8]*k[8][i] + b[10]*k[10][i] + b[11]*k[11][i]) + (b[12]*k[12][i] + b[14]*k[14][i] + b[15]*k[15][i] + b[16]*k[16][i]) + (b[18]*k[18][i] + b[19]*k[19][i] + b[20]*k[20][i] + b[21]*k[21][i]) + (b[22]*k[22][i] + b[23]*k[23][i] + b[24]*k[24][i] + b[25]*k[25][i]) + (b[26]*k[26][i] + b[27]*k[27][i] + b[28]*k[28][i] + b[29]*k[29][i]) + (b[30]*k[30][i] + b[31]*k[31][i] + b[32]*k[32][i] + b[33]*k[33][i]) + (b[34]*k[34][i] + b[35]*k[35][i])
      end
      if adaptive
        for i in uidx
          utmp[i] = u[i] + update[i]
        end
        EEst = norm(((k[2] - k[34]) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting, so in-place
        for i in uidx
          u[i] = u[i] + update[i]
        end
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:Feagin14,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  adaptiveConst,a0100,a0200,a0201,a0300,a0302,a0400,a0402,a0403,a0500,a0503,a0504,a0600,a0603,a0604,a0605,a0700,a0704,a0705,a0706,a0800,a0805,a0806,a0807,a0900,a0905,a0906,a0907,a0908,a1000,a1005,a1006,a1007,a1008,a1009,a1100,a1105,a1106,a1107,a1108,a1109,a1110,a1200,a1208,a1209,a1210,a1211,a1300,a1308,a1309,a1310,a1311,a1312,a1400,a1408,a1409,a1410,a1411,a1412,a1413,a1500,a1508,a1509,a1510,a1511,a1512,a1513,a1514,a1600,a1608,a1609,a1610,a1611,a1612,a1613,a1614,a1615,a1700,a1712,a1713,a1714,a1715,a1716,a1800,a1812,a1813,a1814,a1815,a1816,a1817,a1900,a1912,a1913,a1914,a1915,a1916,a1917,a1918,a2000,a2012,a2013,a2014,a2015,a2016,a2017,a2018,a2019,a2100,a2112,a2113,a2114,a2115,a2116,a2117,a2118,a2119,a2120,a2200,a2212,a2213,a2214,a2215,a2216,a2217,a2218,a2219,a2220,a2221,a2300,a2308,a2309,a2310,a2311,a2312,a2313,a2314,a2315,a2316,a2317,a2318,a2319,a2320,a2321,a2322,a2400,a2408,a2409,a2410,a2411,a2412,a2413,a2414,a2415,a2416,a2417,a2418,a2419,a2420,a2421,a2422,a2423,a2500,a2508,a2509,a2510,a2511,a2512,a2513,a2514,a2515,a2516,a2517,a2518,a2519,a2520,a2521,a2522,a2523,a2524,a2600,a2605,a2606,a2607,a2608,a2609,a2610,a2612,a2613,a2614,a2615,a2616,a2617,a2618,a2619,a2620,a2621,a2622,a2623,a2624,a2625,a2700,a2705,a2706,a2707,a2708,a2709,a2711,a2712,a2713,a2714,a2715,a2716,a2717,a2718,a2719,a2720,a2721,a2722,a2723,a2724,a2725,a2726,a2800,a2805,a2806,a2807,a2808,a2810,a2811,a2813,a2814,a2815,a2823,a2824,a2825,a2826,a2827,a2900,a2904,a2905,a2906,a2909,a2910,a2911,a2913,a2914,a2915,a2923,a2924,a2925,a2926,a2927,a2928,a3000,a3003,a3004,a3005,a3007,a3009,a3010,a3013,a3014,a3015,a3023,a3024,a3025,a3027,a3028,a3029,a3100,a3102,a3103,a3106,a3107,a3109,a3110,a3113,a3114,a3115,a3123,a3124,a3125,a3127,a3128,a3129,a3130,a3200,a3201,a3204,a3206,a3230,a3231,a3300,a3302,a3332,a3400,a3401,a3402,a3404,a3406,a3407,a3409,a3410,a3411,a3412,a3413,a3414,a3415,a3416,a3417,a3418,a3419,a3420,a3421,a3422,a3423,a3424,a3425,a3426,a3427,a3428,a3429,a3430,a3431,a3432,a3433,b,c = constructFeagin14(uEltypeNoUnits)
  k = Vector{typeof(u)}(35)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      k[1]  = Δt*f(t,u)
      k[2]  = Δt*f(t + c[1]*Δt,u + a0100*k[1])
      k[3]  = Δt*f(t + c[2]*Δt ,u + a0200*k[1] + a0201*k[2])
      k[4]  = Δt*f(t + c[3]*Δt,u + a0300*k[1]              + a0302*k[3])
      k[5]  = Δt*f(t + c[4]*Δt,u + a0400*k[1]              + a0402*k[3] + a0403*k[4])
      k[6]  = Δt*f(t + c[5]*Δt,u + a0500*k[1]                           + a0503*k[4] + a0504*k[5])
      k[7]  = Δt*f(t + c[6]*Δt,u + a0600*k[1]                           + a0603*k[4] + a0604*k[5] + a0605*k[6])
      k[8]  = Δt*f(t + c[7]*Δt,u + a0700*k[1]                                        + a0704*k[5] + a0705*k[6] + a0706*k[7])
      k[9]  = Δt*f(t + c[8]*Δt,u + a0800*k[1]                                                     + a0805*k[6] + a0806*k[7] + a0807*k[8])
      k[10] = Δt*f(t + c[9]*Δt,u + a0900*k[1]                                                     + a0905*k[6] + a0906*k[7] + a0907*k[8] + a0908*k[9])
      k[11] = Δt*f(t + c[10]*Δt,u + a1000*k[1]                                                     + a1005*k[6] + a1006*k[7] + a1007*k[8] + a1008*k[9] + a1009*k[10])
      k[12] = Δt*f(t + c[11]*Δt,u + a1100*k[1]                                                     + a1105*k[6] + a1106*k[7] + a1107*k[8] + a1108*k[9] + a1109*k[10] + a1110*k[11])
      k[13] = Δt*f(t + c[12]*Δt,u + a1200*k[1]                                                                                            + a1208*k[9] + a1209*k[10] + a1210*k[11] + a1211*k[12])
      k[14] = Δt*f(t + c[13]*Δt,u + a1300*k[1]                                                                                            + a1308*k[9] + a1309*k[10] + a1310*k[11] + a1311*k[12] + a1312*k[13])
      k[15] = Δt*f(t + c[14]*Δt,u + a1400*k[1]                                                                                            + a1408*k[9] + a1409*k[10] + a1410*k[11] + a1411*k[12] + a1412*k[13] + a1413*k[14])
      k[16] = Δt*f(t + c[15]*Δt,u + a1500*k[1]                                                                                            + a1508*k[9] + a1509*k[10] + a1510*k[11] + a1511*k[12] + a1512*k[13] + a1513*k[14] + a1514*k[15])
      k[17] = Δt*f(t + c[16]*Δt,u + a1600*k[1]                                                                                            + a1608*k[9] + a1609*k[10] + a1610*k[11] + a1611*k[12] + a1612*k[13] + a1613*k[14] + a1614*k[15] + a1615*k[16])
      k[18] = Δt*f(t + c[17]*Δt,u + a1700*k[1]                                                                                                                                                   + a1712*k[13] + a1713*k[14] + a1714*k[15] + a1715*k[16] + a1716*k[17])
      k[19] = Δt*f(t + c[18]*Δt,u + a1800*k[1]                                                                                                                                                   + a1812*k[13] + a1813*k[14] + a1814*k[15] + a1815*k[16] + a1816*k[17] + a1817*k[18])
      k[20] = Δt*f(t + c[19]*Δt,u + a1900*k[1]                                                                                                                                                   + a1912*k[13] + a1913*k[14] + a1914*k[15] + a1915*k[16] + a1916*k[17] + a1917*k[18] + a1918*k[19])
      k[21] = Δt*f(t + c[20]*Δt,u + a2000*k[1]                                                                                                                                                   + a2012*k[13] + a2013*k[14] + a2014*k[15] + a2015*k[16] + a2016*k[17] + a2017*k[18] + a2018*k[19] + a2019*k[20])
      k[22] = Δt*f(t + c[21]*Δt,u + a2100*k[1]                                                                                                                                                   + a2112*k[13] + a2113*k[14] + a2114*k[15] + a2115*k[16] + a2116*k[17] + a2117*k[18] + a2118*k[19] + a2119*k[20] + a2120*k[21])
      k[23] = Δt*f(t + c[22]*Δt,u + a2200*k[1]                                                                                                                                                   + a2212*k[13] + a2213*k[14] + a2214*k[15] + a2215*k[16] + a2216*k[17] + a2217*k[18] + a2218*k[19] + a2219*k[20] + a2220*k[21] + a2221*k[22])
      k[24] = Δt*f(t + c[23]*Δt,u + a2300*k[1]                                                                                            + a2308*k[9] + a2309*k[10] + a2310*k[11] + a2311*k[12] + a2312*k[13] + a2313*k[14] + a2314*k[15] + a2315*k[16] + a2316*k[17] + a2317*k[18] + a2318*k[19] + a2319*k[20] + a2320*k[21] + a2321*k[22] + a2322*k[23])
      k[25] = Δt*f(t + c[24]*Δt,u + a2400*k[1]                                                                                            + a2408*k[9] + a2409*k[10] + a2410*k[11] + a2411*k[12] + a2412*k[13] + a2413*k[14] + a2414*k[15] + a2415*k[16] + a2416*k[17] + a2417*k[18] + a2418*k[19] + a2419*k[20] + a2420*k[21] + a2421*k[22] + a2422*k[23] + a2423*k[24])
      k[26] = Δt*f(t + c[25]*Δt,u + a2500*k[1]                                                                                            + a2508*k[9] + a2509*k[10] + a2510*k[11] + a2511*k[12] + a2512*k[13] + a2513*k[14] + a2514*k[15] + a2515*k[16] + a2516*k[17] + a2517*k[18] + a2518*k[19] + a2519*k[20] + a2520*k[21] + a2521*k[22] + a2522*k[23] + a2523*k[24] + a2524*k[25])
      k[27] = Δt*f(t + c[26]*Δt,u + a2600*k[1]                                                     + a2605*k[6] + a2606*k[7] + a2607*k[8] + a2608*k[9] + a2609*k[10] + a2610*k[11]               + a2612*k[13] + a2613*k[14] + a2614*k[15] + a2615*k[16] + a2616*k[17] + a2617*k[18] + a2618*k[19] + a2619*k[20] + a2620*k[21] + a2621*k[22] + a2622*k[23] + a2623*k[24] + a2624*k[25] + a2625*k[26])
      k[28] = Δt*f(t + c[27]*Δt,u + a2700*k[1]                                                     + a2705*k[6] + a2706*k[7] + a2707*k[8] + a2708*k[9] + a2709*k[10]               + a2711*k[12] + a2712*k[13] + a2713*k[14] + a2714*k[15] + a2715*k[16] + a2716*k[17] + a2717*k[18] + a2718*k[19] + a2719*k[20] + a2720*k[21] + a2721*k[22] + a2722*k[23] + a2723*k[24] + a2724*k[25] + a2725*k[26] + a2726*k[27])
      k[29] = Δt*f(t + c[28]*Δt,u + a2800*k[1]                                                     + a2805*k[6] + a2806*k[7] + a2807*k[8] + a2808*k[9]               + a2810*k[11] + a2811*k[12]               + a2813*k[14] + a2814*k[15] + a2815*k[16]                                                                                                   + a2823*k[24] + a2824*k[25] + a2825*k[26] + a2826*k[27] + a2827*k[28])
      k[30] = Δt*f(t + c[29]*Δt,u + a2900*k[1]                                        + a2904*k[5] + a2905*k[6] + a2906*k[7]                           + a2909*k[10] + a2910*k[11] + a2911*k[12]               + a2913*k[14] + a2914*k[15] + a2915*k[16]                                                                                                   + a2923*k[24] + a2924*k[25] + a2925*k[26] + a2926*k[27] + a2927*k[28] + a2928*k[29])
      k[31] = Δt*f(t + c[30]*Δt,u + a3000*k[1]                           + a3003*k[4] + a3004*k[5] + a3005*k[6]              + a3007*k[8]              + a3009*k[10] + a3010*k[11]                             + a3013*k[14] + a3014*k[15] + a3015*k[16]                                                                                                   + a3023*k[24] + a3024*k[25] + a3025*k[26]               + a3027*k[28] + a3028*k[29] + a3029*k[30])
      k[32] = Δt*f(t + c[31]*Δt,u + a3100*k[1]              + a3102*k[3] + a3103*k[4]                           + a3106*k[7] + a3107*k[8]              + a3109*k[10] + a3110*k[11]                             + a3113*k[14] + a3114*k[15] + a3115*k[16]                                                                                                   + a3123*k[24] + a3124*k[25] + a3125*k[26]               + a3127*k[28] + a3128*k[29] + a3129*k[30] + a3130*k[31])
      k[33] = Δt*f(t + c[32]*Δt,u + a3200*k[1] + a3201*k[2]                           + a3204*k[5]              + a3206*k[7]                                                                                                                                                                                                                                                                                                                                 + a3230*k[31] + a3231*k[32])
      k[34] = Δt*f(t + c[33]*Δt,u + a3300*k[1]              + a3302*k[3]                                                                                                                                                                                                                                                                                                                                                                                                                 + a3332*k[33])
      k[35] = Δt*f(t + c[34]*Δt,u + a3400*k[1] + a3401*k[2] + a3402*k[3]              + a3404*k[5]              + a3406*k[7] + a3407*k[8]              + a3409*k[10] + a3410*k[11] + a3411*k[12] + a3412*k[13] + a3413*k[14] + a3414*k[15] + a3415*k[16] + a3416*k[17] + a3417*k[18] + a3418*k[19] + a3419*k[20] + a3420*k[21] + a3421*k[22] + a3422*k[23] + a3423*k[24] + a3424*k[25] + a3425*k[26] + a3426*k[27] + a3427*k[28] + a3428*k[29] + a3429*k[30] + a3430*k[31] + a3431*k[32] + a3432*k[33] + a3433*k[34])
      update = (b[1]*k[1] + b[2]*k[2] + b[3]*k[3] + b[5]*k[5]) + (b[7]*k[7] + b[8]*k[8] + b[10]*k[10] + b[11]*k[11]) + (b[12]*k[12] + b[14]*k[14] + b[15]*k[15] + b[16]*k[16]) + (b[18]*k[18] + b[19]*k[19] + b[20]*k[20] + b[21]*k[21]) + (b[22]*k[22] + b[23]*k[23] + b[24]*k[24] + b[25]*k[25]) + (b[26]*k[26] + b[27]*k[27] + b[28]*k[28] + b[29]*k[29]) + (b[30]*k[30] + b[31]*k[31] + b[32]*k[32] + b[33]*k[33]) + (b[34]*k[34] + b[35]*k[35])
      if adaptive
        utmp = u + update
        EEst = norm(((k[2] - k[34]) * adaptiveConst)./(abstol+u*reltol),internalnorm)
      else #no chance of rejecting, so in-place
        u = u + update
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:ImplicitEuler,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
      @ode_numberimplicitloopfooter
    end
  end
  u = uhold[1]
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:ImplicitEuler,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
      @ode_implicitloopfooter
    end
  end
  u = reshape(uhold,sizeu...)
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Trapezoid,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
      @ode_implicitloopfooter
    end
  end
  u = reshape(uhold,sizeu...)
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:Trapezoid,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
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
      @ode_numberimplicitloopfooter
    end
  end
  u = uhold[1]
  return u,t,timeseries,ts
end

function ode_solve{uType<:AbstractArray,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:AbstractArray}(integrator::ODEIntegrator{:Rosenbrock32,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  k₁::uType = similar(u)
  k₂ = similar(u)
  k₃::uType = similar(u)
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
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      ForwardDiff.derivative!(dT,(t)->vecfreturn(t,u,du2),t) # Time derivative
      ForwardDiff.jacobian!(J,(du1,u)->vecf(t,u,du1),vec(du1),vec(u))

      W[:] = one(J)-Δt*d*J # Can an allocation be cut here?
      f(t,u,f₀)
      @into! vectmp = W\vec(f₀ + Δt*d*dT)
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
        f(t+Δt,utmp,f₂)
        @into! vectmp3 = W\vec(f₂ - c₃₂*(k₂-f₁)-2(k₁-f₀)+Δt*d*T)
        k₃ = reshape(vectmp3,sizeu...)
        for i in uidx
          tmp2[i] = (Δt*(k₁[i] - 2k₂[i] + k₃[i])/6)./(abstol+u[i]*reltol)
        end
        EEst = norm(tmp2,internalnorm)
      else
        for i in uidx
          u[i] = u[i] + Δt*k₂[i]
        end
      end
      @ode_loopfooter
    end
  end
  return u,t,timeseries,ts
end

function ode_solve{uType<:Number,uEltype<:Number,N,tType<:Number,uEltypeNoUnits<:Number,rateType<:Number}(integrator::ODEIntegrator{:Rosenbrock32,uType,uEltype,N,tType,uEltypeNoUnits,rateType})
  @ode_preamble
  c₃₂ = 6 + sqrt(2)
  d = 1/(2+sqrt(2))
  local dT::uType
  local J::uType
  local f₀::uType
  local k₁::uType
  local f₁::uType
  local f₂::uType
  local k₂::uType
  local k₃::uType
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      # Time derivative
      dT = ForwardDiff.derivative((t)->f(t,u),t)
      J = ForwardDiff.derivative((u)->f(t,u),u)
      W = one(J)-Δt*d*J
      f₀ = f(t,u)
      k₁ = W\(f₀ + Δt*d*dT)
      f₁ = f(t+Δt/2,u+Δt*k₁/2)
      k₂ = W\(f₁-k₁) + k₁
      if adaptive
        utmp = u + Δt*k₂
        f₂ = f(t+Δt,utmp)
        k₃ = W\(f₂ - c₃₂*(k₂-f₁)-2(k₁-f₀)+Δt*d*T)
        EEst = norm((Δt*(k₁ - 2k₂ + k₃)/6)./(abstol+u*reltol),internalnorm)
      else
        u = u + Δt*k₂
      end
      @ode_numberloopfooter
    end
  end
  return u,t,timeseries,ts
end

max{ND,N,m,kg,s,A,K,mol,cd,rad,sr}(u::Array{SIUnits.SIQuantity{N,m,kg,s,A,K,mol,cd,rad,sr},ND},utmp::Array{SIUnits.SIQuantity{N,m,kg,s,A,K,mol,cd,rad,sr},ND}) = map((x,y)->SIUnits.SIQuantity{N,m,kg,s,A,K,mol,cd,rad,sr}(max(x.val,y.val)),u,utmp)
