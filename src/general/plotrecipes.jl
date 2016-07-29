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

@recipe function f(sol::FEMSolution;plottrue=false,tslocation=0)
  if tslocation==0 #Plot solution at end
    out = Any[]
    for i = 1:size(sol.u,2)
      push!(out,sol.u[:,i])
    end
    if plottrue
      for i = 1:size(sol.u,2)
        push!(out,sol.uTrue[:,i])
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

@recipe function f(sol::SDESolution;plottrue=false)
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
  if plottrue
    if ndims(sol.analytics) > 2
      totaldims = 1
      for i=2:ndims(sol.timeseries)
        totaldims *= size(sol.analytics,i)
      end
      vals = [vals reshape(sol.analytics,size(sol.analytics,1),totaldims)]
    else
      vals = [vals sol.analytics]
    end
  end
  #u = Any[sol.timeseries];
  #plottrue && push!(u, sol.analytics);
  seriestype --> :path
  #layout --> length(u)
  map(Float64,sol.ts), map(Float64,vals) #Remove when Tom commits
end

@recipe function f(sol::ODESolution;plottrue=false)
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
  if plottrue
    if ndims(sol.analytics) > 2
      totaldims = 1
      for i=2:ndims(sol.timeseries)
        totaldims *= size(sol.analytics,i)
      end
      vals = [vals reshape(sol.analytics,size(sol.analytics,1),totaldims)]
    else
      vals = [vals sol.analytics]
    end
  end
  #u = Any[sol.timeseries];
  #plottrue && push!(u, sol.analytics);
  seriestype --> :path
  #layout --> length(u)
  sol.ts, vals
end

@recipe function f(sim::ConvergenceSimulation)
  if ndims(collect(values(sim.errors))[1])>1 #Monte Carlo
    vals = [mean(x,1)' for x in values(sim.errors)]
  else #Deterministic
    vals = [x for x in values(sim.errors)]
  end
  seriestype --> :path
  label  --> [key for key in keys(sim.errors)]
  xguide  --> "Convergence Axis"
  yguide  --> "Error"
  xscale --> :log10
  yscale --> :log10
  sim.convergence_axis, vals
end


@recipe function f(mesh::Mesh)
  seriestype --> :surface #:plot_trisurf
  #triangles  --> mesh.elem-1
  mesh.node[:,1], mesh.node[:,2], ones(mesh.node[:,1])
end

# mesh = meshExample_lakemesh()
# PyPlot.plot_trisurf(mesh.node[:,1],mesh.node[:,2],ones(mesh.node[:,2]),triangles=mesh.elem-1)
