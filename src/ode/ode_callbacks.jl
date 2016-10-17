macro ode_callback(ex)
  esc(quote
    function (alg,f,t,u,k,tprev,uprev,kprev,ts,timeseries,ks,Δtprev,Δt,saveat,cursaveat,saveiter,iter,save_timeseries,timeseries_steps,uEltype,ksEltype,dense,kshortsize,issimple_dense,fsal,fsalfirst,cache,calck,T,Ts)
      reeval_fsal = false
      event_occured = false
      $(ex)
      cursaveat,saveiter,Δt,t,T,reeval_fsal
    end
  end)
end

macro ode_event(event_f,apply_event!,rootfind_event_loc=true,interp_points=5,terminate_on_event=false,Δt_safety=1)
  esc(quote
    # Event Handling
    if $interp_points!=0
      ode_addsteps!(k,tprev,uprev,Δtprev,alg,f)
      Θs = linspace(0,1,$(interp_points))
    end
    interp_index = 0
    # Check if the event occured
    if $event_f(t,u)<0
      event_occured = true
      interp_index = $interp_points
    elseif $interp_points!=0 # Use the interpolants for safety checking
      for i in 2:length(Θs)-1
        if $event_f(t+Δt*Θs[i],ode_interpolant(Θs[i],Δtprev,uprev,u,kprev,k,alg))<0
          event_occured = true
          interp_index = i
          break
        end
      end
    end

    if event_occured
      if interp_index == $interp_points # If no safety interpolations, start in the middle as well
        initial_Θ = [.5]
      else
        initial_Θ = [Θs[interp_index]] # Start at the closest
      end
      if $rootfind_event_loc
        find_zero = (Θ,val) -> begin
          val[1] = event_f(t+Θ[1]*Δt,ode_interpolant(Θ[1],Δtprev,uprev,u,kprev,k,alg))
        end
        res = nlsolve(find_zero,initial_Θ)
        val = ode_interpolant(res.zero[1],Δtprev,uprev,u,kprev,k,alg)
        copy!(u,val)
        Δtprev *= res.zero[1]
      elseif interp_index != $interp_points
          Δtprev *= Θs[interp_index]
          copy!(u,ode_interpolant(Θs[interp_index],Δtprev,uprev,u,kprev,k,alg))
      end
      # If no solve and no interpolants, just use endpoint

      t = tprev + Δtprev
      if calck
        if alg ∈ DIFFERENTIALEQUATIONSJL_SPECIALDENSEALGS
          resize!(k,kshortsize) # Reset k for next step
          k = typeof(k)() # Make a local blank k for saving
          ode_addsteps!(k,tprev,uprev,Δtprev,alg,f)
        elseif typeof(u) <: Number
          k = f(t,u)
        else
          f(t,u,k)
        end
      end
    end

    @ode_savevalues

    if event_occured
      if $terminate_on_event
        @ode_terminate
      else
        $apply_event!(u,cache)
        if calck
          if alg ∉ DIFFERENTIALEQUATIONSJL_SPECIALDENSEALGS
            if typeof(u) <: Number
              k = f(t,u)
            else
              f(t,u,k)
            end
          end
        end
        @ode_savevalues
        if fsal
          reeval_fsal = true
        end
        Δt *= $Δt_safety # Safety Δt change
      end
    end
  end)
end

macro ode_change_cachesize(cache,resize_ex)
  resize_ex = cache_replace_length(resize_ex)
  esc(quote
    for i in 1:length($cache)
      resize!($cache[i],$resize_ex)
    end
  end)
end

macro ode_change_deleteat(cache,deleteat_ex)
  deleteat_ex = cache_replace_length(deleteat_ex)
  esc(quote
    for i in 1:length($cache)
      deleteat!($cache[i],$deleteat_ex)
    end
  end)
end

function cache_replace_length(ex::Expr)
  for (i,arg) in enumerate(ex.args)
    if isa(arg,Expr)
      cache_replace_length(ex)
    elseif isa(arg,Symbol)
      if arg == :length
        ex.args[i] = :(length(cache[i]))
      end
    end
  end
  ex
end

function cache_replace_length(ex::Symbol)
  if ex == :length
    ex = :(length(cache[i]))
  end
  ex
end

function cache_replace_length(ex::Any)
  ex
end

@def ode_terminate begin
  T = t
  while length(Ts)>1
    pop!(Ts)
  end
end
