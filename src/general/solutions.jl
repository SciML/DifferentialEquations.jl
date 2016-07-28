"""
FEMSolution

Holds the data for the solution to a finite element problem.

### Fields

* `fem_mesh::FEMmesh`: The finite element mesh the problem was solved on.
* `u::Array{Float64}`: The solution (at the final timepoint)
* `trueKnown::Bool`: Boolean flag for if the true solution is given.
* `uTrue::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `appxTrue::Bool`: Boolean flag for if uTrue was an approximation.
* `timeseries`::AbstractArrayOrVoid`: u over time. Only saved if `save_timeseries=true`
is specified in the solver.
* `ts::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
is specified in the solver.
* `prob::DEProblem`: Holds the problem object used to define the problem.
* `save_timeseries::Bool`: True if solver saved the extra timepoints.

"""
type FEMSolution <: DESolution
  fem_mesh::FEMmesh
  u#::Array{Number}
  trueKnown::Bool
  uTrue::AbstractArrayOrVoid
  errors#::Dict{String,Float64}
  appxTrue::Bool
  timeSeries#::GrowableArray
  ts::AbstractArrayOrVoid
  prob::DEProblem
  save_timeseries::Bool
  function FEMSolution(fem_mesh::FEMmesh,u,uTrue,sol,Du,timeSeries,ts,prob;save_timeseries=true)
    errors = Dict("L2"=>getL2error(fem_mesh,sol,u),"H1"=>getH1error(fem_mesh,Du,u),
                  :l∞=> maximum(abs(u-uTrue)), :l2=> norm(u-uTrue,2))
    return(new(fem_mesh,u,true,uTrue,errors,false,timeSeries,ts,prob,true))
  end
  FEMSolution(fem_mesh,u,uTrue,sol,Du,prob) = FEMSolution(fem_mesh::FEMmesh,u,uTrue,sol,Du,nothing,nothing,prob,save_timeseries=false)
  function FEMSolution(fem_mesh::FEMmesh,u::AbstractArray,prob)
    return(FEMSolution(fem_mesh,u,nothing,nothing,prob,save_timeseries=false))
  end
  function FEMSolution(fem_mesh::FEMmesh,u::AbstractArray,timeSeries,ts,prob;save_timeseries=true)
    return(new(fem_mesh,u,false,nothing,Dict{String,Float64},false,timeSeries,ts,prob,save_timeseries))
  end
end

"""
SDESolution

Holds the data for the solution to a SDE problem.

### Fields

* `u::Array{Float64}`: The solution (at the final timepoint)
* `trueKnown::Bool`: Boolean flag for if the true solution is given.
* `uTrue::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `timeseries`::AbstractArrayOrVoid`: u over time. Only saved if `save_timeseries=true`
is specified in the solver.
* `ts::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
is specified in the solver.
* `Ws`: All of the W's in the solution. Only saved if `save_timeseries=true` is specified
in the solver.
* `sols`: If `save_timeseries=true`, saves the solution at each save point.
* `prob::DEProblem`: Holds the problem object used to define the problem.
* `save_timeseries::Bool`: True if solver saved the extra timepoints.
* `appxTrue::Bool`: Boolean flag for if uTrue was an approximation.

"""
type SDESolution <: DESolution
  u#::AbstractArrayOrNumber
  trueKnown::Bool
  uTrue#::AbstractArrayOrNumber
  errors#::Dict{}
  timeseries::AbstractArrayOrVoid
  ts::AbstractArrayOrVoid
  Δts::AbstractArrayOrVoid
  Ws::AbstractArrayOrVoid
  sols::AbstractArrayOrVoid
  appxTrue::Bool
  save_timeseries::Bool
  maxStackSize::Int
  W
  function SDESolution(u;timeseries=nothing,sols=nothing,ts=nothing,Δts=nothing,Ws=nothing,maxStackSize=nothing,W=nothing)
    save_timeseries = timeseries == nothing
    trueKnown = false
    return(new(u,trueKnown,nothing,Dict(),timeseries,ts,Δts,Ws,sols,false,save_timeseries,maxStackSize,W))
  end
  function SDESolution(u,uTrue;timeseries=nothing,sols=nothing,ts=nothing,Δts=nothing,Ws=nothing,maxStackSize=nothing,W=nothing)
    save_timeseries = timeseries != nothing
    trueKnown = true
    errors = Dict(:final=>mean(abs(u-uTrue)))
    if save_timeseries
      errors = Dict(:final=>mean(abs(u-uTrue)),:l∞=>maximum(abs(timeseries-sols)),:l2=>sqrt(mean((timeseries-sols).^2)))
    end
    return(new(u,trueKnown,uTrue,errors,timeseries,ts,Δts,Ws,sols,false,save_timeseries,maxStackSize,W))
  end
  #Required to convert pmap results
  SDESolution(a::Any) = new(a.u,a.trueKnown,a.uTrue,a.errors,a.timeseries,a.ts,a.Δts,a.Ws,a.sols,a.appxTrue,a.save_timeseries,a.maxStackSize,a.W)
end

"""
ODESolution

Holds the data for the solution to an ODE problem.

### Fields

* `u::Array{Float64}`: The solution (at the final timepoint)
* `trueKnown::Bool`: Boolean flag for if the true solution is given.
* `uTrue::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `timeseries`::AbstractArrayOrVoid`: u over time. Only saved if `save_timeseries=true`
is specified in the solver.
* `ts::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
is specified in the solver.
* `sols`: If `save_timeseries=true`, saves the solution at each timestep.
* `prob::DEProblem`: Holds the problem object used to define the problem.
* `save_timeseries::Bool`: True if solver saved the extra timepoints.
* `appxTrue::Bool`: Boolean flag for if uTrue was an approximation.

"""
type ODESolution <: DESolution
  u#::AbstractArrayOrNumber
  trueKnown::Bool
  uTrue#::AbstractArrayOrNumber
  errors#::Dict{}
  timeseries::AbstractArrayOrVoid
  ts::AbstractArrayOrVoid
  sols::AbstractArrayOrVoid
  appxTrue::Bool
  save_timeseries::Bool
  function ODESolution(u;timeseries=nothing,sols=nothing,ts=nothing)
    save_timeseries = timeseries == nothing
    trueKnown = false
    return(new(u,trueKnown,nothing,Dict(),timeseries,ts,sols,false,save_timeseries))
  end
  function ODESolution(u,uTrue;timeseries=nothing,sols=nothing,ts=nothing)
    save_timeseries = timeseries != nothing
    trueKnown = true
    errors = Dict(:final=>mean(abs(u-uTrue)))
    if save_timeseries
      errors = Dict(:final=>mean(abs(u-uTrue)),:l∞=>maximum(abs(timeseries-sols)),:l2=>sqrt(mean((timeseries-sols).^2)))
    end
    return(new(u,trueKnown,uTrue,errors,timeseries,ts,sols,false,save_timeseries))
  end
end

"""
StokesSolution

Holds the data for the solution to a Stokes problem.

### Fields

* u
* v
* p
* uTrue
* vTrue
* pTrue
* mesh
* trueKnown
* errors
* converrors

"""
type StokesSolution <: DESolution
  u
  v
  p
  uTrue
  vTrue
  pTrue
  mesh::FDMMesh
  trueKnown::Bool
  errors
  converrors
  StokesSolution(u,v,p,uTrue,vTrue,pTrue,mesh,trueKnown;errors=nothing,converrors=nothing) = new(u,v,p,uTrue,vTrue,pTrue,mesh,trueKnown,errors,converrors)
end

"""
appxTrue!(res,res2)

Adds the solution from res2 to the FEMSolution object res.
Useful to add a quasi-true solution when none is known by
computing once at a very small time/space step and taking
that solution as the "true" solution
"""
function appxTrue!(res::FEMSolution,res2::FEMSolution)
  res.uTrue = res2.u
  res.errors = Dict(:l∞=>maximum(abs(res.u-res.uTrue)),:l2=>norm(res.u-res.uTrue,2))
  res.appxTrue = true
end

"""
S = FEMSolutionTS(timeseries::GrowableArray,numvars::Int)
S[i][j] => Variable i at time j.
"""
function FEMSolutionTS(timeseries::GrowableArray,numvars::Int)
  G = GrowableArray(timeseries[1][:,1])
  for j = 2:length(timeseries)
    push!(G,timeseries[j][:,1])
  end
  ts = GrowableArray(G)
  if numvars > 1
    for i=2:numvars
      G = GrowableArray(timeseries[1][:,i])
      for j = 2:length(timeseries)
        push!(G,timeseries[j][:,i])
      end
      push!(ts,G)
    end
  end
  return(ts)
end
