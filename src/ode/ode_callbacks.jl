macro ode_callback(ex)
  esc(quote
    function (alg,f,t,u,k,tprev,uprev,kprev,ts,timeseries,ks,Δtprev,Δt,saveat,cursaveat,iter,save_timeseries,timeseries_steps,uEltype,ksEltype,dense,kshortsize,issimple_dense,fsal,fsalfirst,cache)
      reeval_fsal = false
      event_occured = false
      $(ex)
      cursaveat,Δt,t,reeval_fsal
    end
  end)
end

macro ode_event(event_f,apply_event!,interp_points=0,Δt_safety=1)
  esc(quote
    # Event Handling
    ode_addsteps!(k,tprev,uprev,Δtprev,alg,f)
    Θs = linspace(0,1,$(interp_points))
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
      find_zero = (Θ,val) -> begin
        val[1] = event_f(t+Θ[1]*Δt,ode_interpolant(Θ[1],Δtprev,uprev,u,kprev,k,alg))
      end
      res = nlsolve(find_zero,initial_Θ)
      val = ode_interpolant(res.zero[1],Δtprev,uprev,u,kprev,k,alg)
      for i in eachindex(u)
        u[i] = val[i]
      end
      Δtprev *= res.zero[1]
      t = tprev + Δtprev

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

    @ode_savevalues

    if event_occured
      $apply_event!(u,cache)
      if alg ∉ DIFFERENTIALEQUATIONSJL_SPECIALDENSEALGS
        if typeof(u) <: Number
          k = f(t,u)
        else
          f(t,u,k)
        end
      end
      @ode_savevalues
      if fsal
        reeval_fsal = true
      end
      Δt *= $Δt_safety # Safety Δt change
    end
  end)
end
