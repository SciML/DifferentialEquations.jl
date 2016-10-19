"""
`animate(sol::FEMSolution)`

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

@recipe function f(sol::AbstractODESolution;sensitivity=false,plot_analytic=false,denseplot=true,plotdensity=100)
  plotseries = Vector{Any}(0)
  if sensitivity
    plott = sol.t
    names = Vector{Symbol}(0)
    for (k,v) in sol.sensitivity.normseries
      push!(plotseries,v)
      push!(names,k)
    end
    label  --> reshape(names,1,length(names))
  else
    if typeof(sol) <: SDESolution; denseplot=false; end

    if denseplot && sol.dense # Generate the points from the plot from dense function
      plott = collect(Ranges.linspace(sol.t[1],sol.t[end],plotdensity))
      plot_timeseries = sol(plott)
      if plot_analytic
        plot_analytic_timeseries = Vector{typeof(sol.u)}(length(plott))
        for i in eachindex(plott)
          tmp[i] = sol.prob.analytic(plott[i],sol.prob.u₀)
        end
      end
    else # Plot for not dense output use the timeseries itself
      plot_timeseries = sol.timeseries
      if plot_analytic
        plot_analytic_timeseries = sol.timeseries_analytic
      end
      plott = sol.t
    end

    # Make component-wise plots
    if typeof(sol.u) <:AbstractArray
      for i in eachindex(sol.u)
        tmp = Vector{eltype(sol.u)}(length(plot_timeseries))
        for j in 1:length(plot_timeseries)
          tmp[j] = plot_timeseries[j][i]
        end
        push!(plotseries,tmp)
      end
    else
      push!(plotseries,plot_timeseries)
    end
    if plot_analytic
      if typeof(sol.u) <: AbstractArray
        for i in eachindex(sol.u)
          tmp = Vector{eltype(sol.u)}(length(plot_timeseries))
          for j in 1:length(plot_timeseries)
            tmp[j] = plot_analytic_timeseries[j][i]
          end
          push!(plotseries,tmp)
        end
      else
        push!(plotseries,plot_analytic_timeseries)
      end
    end
  end
  seriestype --> :path
  lw --> 3
  xtickfont --> font(11)
  ytickfont --> font(11)
  legendfont --> font(11)
  guidefont  --> font(11)
  #layout --> length(u)
  plott, plotseries
end

@recipe function f(sim::ConvergenceSimulation)
  if ndims(collect(values(sim.errors))[1])>1 #Monte Carlo
    vals = [mean(x,1)' for x in values(sim.errors)]
  else #Deterministic
    vals = [x for x in values(sim.errors)]
  end
  seriestype --> :path
  label  --> permuteddims([string(key) for key in keys(sim.errors)],[2,1])
  xguide  --> "Convergence Axis"
  yguide  --> "Error"
  xscale --> :log10
  yscale --> :log10
  sim.convergence_axis, vals
end

@recipe function f(shoot::Shootout)
  seriestype --> :bar
  legend := false
  xguide --> "Algorithms"
  yguide --> "Efficiency"
  shoot.names,shoot.effs
end

@recipe function f(wp::WorkPrecision)
  seriestype --> :path
  label -->  wp.name
  lw --> 3
  yguide --> "Time (s)"
  xguide --> "Error"
  xscale --> :log10
  yscale --> :log10
  wp.errors,wp.times
end

@recipe function f(wp_set::WorkPrecisionSet)
  seriestype --> :path
  label -->  reshape(wp_set.names,1,length(wp_set))
  lw --> 3
  yguide --> "Time (s)"
  xguide --> "Error"
  xscale --> :log10
  yscale --> :log10
  errors = Vector{Any}(0)
  times = Vector{Any}(0)
  for i in 1:length(wp_set)
    push!(errors,wp_set[i].errors)
    push!(times,wp_set[i].times)
  end
  errors,times
end

@recipe function f(mesh::Mesh)
  seriestype --> :surface #:plot_trisurf
  #triangles  --> mesh.elem-1
  mesh.node[:,1], mesh.node[:,2], ones(mesh.node[:,1])
end

# mesh = meshExample_lakemesh()
# PyPlot.plot_trisurf(mesh.node[:,1],mesh.node[:,2],ones(mesh.node[:,2]),triangles=mesh.elem-1)

@recipe function f(tab::ODERKTableau;Δx=1/100,Δy=1/100,xlim=[-6,1],ylim=[-5,5])
  x = xlim[1]:Δx:xlim[2]
  y = ylim[1]:Δy:ylim[2]
  f = (u,v)-> abs(stability_region(u+v*im,tab))<1
  seriestype --> :contour
  fill --> true
  cbar --> false
  x,y,f
end
