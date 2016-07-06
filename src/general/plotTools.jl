"""
animate(sol::FEMSolution)

Plots an animation of the solution. Requires `fullSave=true` was enabled in the solver.
"""
function animate(sol::FEMSolution;filename="tmp.gif",fps=15,kw...)
  atomLoaded = isdefined(Main,:Atom)
  anim = Plots.Animation()
  for j=1:length(sol.timeSeries[1])
    plot(sol,tsLocation=j;kw...)
    Plots.frame(anim)
    atomLoaded ? Main.Atom.progress(j/length(sol.timeSeries[1])) : nothing #Use Atom's progressbar if loaded
  end
  gif(anim,filename,fps=fps)
end

@recipe function f(sol::FEMSolution;plottrue=false,tsLocation=0)
  tsLocation = pop!(d,:tsLocation)
  if tsLocation==0 #Plot solution at end
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
    for i = 1:sol.prob.numVars
      push!(out,sol.timeSeries[i][tsLocation])
    end
  end
  seriestype --> :surface
  layout --> length(out)
  sol.femMesh.node[:,1], sol.femMesh.node[:,2], out
end

@recipe function f(sol::SDESolution;plottrue=false)
  if ndims(sol.uFull) > 2
    totaldims = 1
    for i=2:ndims(sol.uFull)
      totaldims *= size(sol.uFull,i)
    end
    #println(totaldims)
    vals = reshape(sol.uFull,size(sol.uFull,1),totaldims)
  else
    vals = sol.uFull
  end
  if plottrue
    if ndims(sol.solFull) > 2
      totaldims = 1
      for i=2:ndims(sol.uFull)
        totaldims *= size(sol.solFull,i)
      end
      vals = [vals reshape(sol.solFull,size(sol.solFull,1),totaldims)]
    else
      vals = [vals sol.solFull]
    end
  end
  #u = Any[sol.uFull];
  #plottrue && push!(u, sol.solFull);
  seriestype --> :path
  #layout --> length(u)
  map(Float64,sol.tFull), map(Float64,vals) #Remove when Tom commits
end

@recipe function f(sol::ODESolution;plottrue=false)
  if ndims(sol.uFull) > 2
    totaldims = 1
    for i=2:ndims(sol.uFull)
      totaldims *= size(sol.uFull,i)
    end
    #println(totaldims)
    vals = reshape(sol.uFull,size(sol.uFull,1),totaldims)
  else
    vals = sol.uFull
  end
  if plottrue
    if ndims(sol.solFull) > 2
      totaldims = 1
      for i=2:ndims(sol.uFull)
        totaldims *= size(sol.solFull,i)
      end
      vals = [vals reshape(sol.solFull,size(sol.solFull,1),totaldims)]
    else
      vals = [vals sol.solFull]
    end
  end
  #u = Any[sol.uFull];
  #plottrue && push!(u, sol.solFull);
  seriestype --> :path
  #layout --> length(u)
  sol.tFull, vals
end

@recipe function f(sim::ConvergenceSimulation)
  if ndims(collect(values(sim.errors))[1])>1 #Monte Carlo
    vals = [mean(x,1) for x in values(sim.errors)]'
  else #Deterministic
    vals = [x for x in values(sim.errors)]
  end
  seriestype --> :path
  label  --> [key for key in keys(sim.errors)]'
  xguide  --> "Convergence Axis"
  yguide  --> "Error"
  xscale --> :log10
  yscale --> :log10
  sim.convergenceAxis, vals
end


@recipe function f(mesh::Mesh)
  seriestype --> :surface #:plot_trisurf
  #triangles  --> mesh.elem-1
  mesh.node[:,1], mesh.node[:,2], ones(mesh.node[:,1])
end

# mesh = meshExample_lakemesh()
# PyPlot.plot_trisurf(mesh.node[:,1],mesh.node[:,2],ones(mesh.node[:,2]),triangles=mesh.elem-1)
