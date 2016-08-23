"""
animate(sol::FEMSolution)

Plots an animation of the solution. Requires `save_timeseries=true` was enabled in the solver.
"""
function animate(sol::FEMSolution;filename="tmp.gif",fps=15,kw...)
  atomloaded = isdefined(Main,:Atom)
  anim = Plots.Animation()
  for j=1:length(sol.timeseries[1])
    plot(sol,tslocation=j;kw...)
    Plots.frame(anim)
    atomloaded ? Main.Atom.progress(j/length(sol.timeseries[1])) : nothing #Use Atom's progressbar if loaded
  end
  gif(anim,filename,fps=fps)
end

@recipe function f(sol::FEMSolution;plot_analytic=false,tslocation=0)
  if tslocation==0 #Plot solution at end
    out = Any[]
    for i = 1:size(sol.u,2)
      push!(out,sol.u[:,i])
    end
    if plot_analytic
      for i = 1:size(sol.u,2)
        push!(out,sol.u_analytic[:,i])
      end
    end
  else #use timeseries
    out = Any[]
    for i = 1:sol.prob.numvars
      push!(out,sol.timeseries[i][tslocation])
    end
  end
  seriestype --> :surface
  layout --> length(out)
  sol.fem_mesh.node[:,1], sol.fem_mesh.node[:,2], out
end

@recipe function f(sol::SDESolution;plot_analytic=false)
  if ndims(sol.timeseries) > 2
    totaldims = 1
    for i=2:ndims(sol.timeseries)
      totaldims *= size(sol.timeseries,i)
    end
    #println(totaldims)
    vals = reshape(sol.timeseries,size(sol.timeseries,1),totaldims)
  else
    vals = sol.timeseries
  end
  if plot_analytic
    if ndims(sol.timeseries_analytic) > 2
      totaldims = 1
      for i=2:ndims(sol.timeseries)
        totaldims *= size(sol.timeseries_analytic,i)
      end
      vals = [vals reshape(sol.timeseries_analytic,size(sol.timeseries_analytic,1),totaldims)]
    else
      vals = [vals sol.timeseries_analytic]
    end
  end
  #u = Any[sol.timeseries];
  #plot_analytic && push!(u, sol.timeseries_analytic);
  seriestype --> :path
  #layout --> length(u)
  map(Float64,sol.t), map(Float64,vals) #Remove when Tom commits
end

@recipe function f(sol::ODESolution;plot_analytic=false)
  if ndims(sol.timeseries) > 2
    totaldims = 1
    for i=2:ndims(sol.timeseries)
      totaldims *= size(sol.timeseries,i)
    end
    #println(totaldims)
    vals = reshape(sol.timeseries,size(sol.timeseries,1),totaldims)
  else
    vals = sol.timeseries
  end
  if plot_analytic
    if ndims(sol.timeseries_analytic) > 2
      totaldims = 1
      for i=2:ndims(sol.timeseries)
        totaldims *= size(sol.timeseries_analytic,i)
      end
      vals = [vals reshape(sol.timeseries_analytic,size(sol.timeseries_analytic,1),totaldims)]
    else
      vals = [vals sol.timeseries_analytic]
    end
  end
  #u = Any[sol.timeseries];
  #plot_analytic && push!(u, sol.timeseries_analytic);
  seriestype --> :path
  #layout --> length(u)
  sol.t, vals
end

@recipe function f(sim::ConvergenceSimulation)
  if ndims(collect(values(sim.errors))[1])>1 #Monte Carlo
    vals = [mean(x,1)' for x in values(sim.errors)]
  else #Deterministic
    vals = [x for x in values(sim.errors)]
  end
  seriestype --> :path
  label  --> [string(key) for key in keys(sim.errors)]'
  xguide  --> "Convergence Axis"
  yguide  --> "Error"
  xscale --> :log10
  yscale --> :log10
  sim.convergence_axis, vals
end

@recipe function f(shoot::Shootout)
  seriestype --> :bar
  legend := false
  xguide --> "Time (s)"
  yguide --> "Efficiency"
  shoot.names,shoot.effs
end

@recipe function f(wp::WorkPrecision)
  seriestype --> :path
  label -->  wp.names'
  xguide --> "Time (s)"
  yguide --> "Error"
  xscale --> :log10
  yscale --> :log10
  wp.times, wp.errors
end

@recipe function f(wp_set::WorkPrecisionSet)
  seriestype --> :path
  label -->  wp_set.names'
  xguide --> "Time (s)"
  yguide --> "Error"
  xscale --> :log10
  yscale --> :log10
  errors = Vector{Any}(0)
  times = Vector{Any}(0)
  for i in 1:length(wp_set)
    push!(errors,wp_set[i].errors)
    push!(times,wp_set[i].times)
  end
  times,errors
end

@recipe function f(mesh::Mesh)
  seriestype --> :surface #:plot_trisurf
  #triangles  --> mesh.elem-1
  mesh.node[:,1], mesh.node[:,2], ones(mesh.node[:,1])
end

# mesh = meshExample_lakemesh()
# PyPlot.plot_trisurf(mesh.node[:,1],mesh.node[:,2],ones(mesh.node[:,2]),triangles=mesh.elem-1)
