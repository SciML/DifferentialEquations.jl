"""
FEMSolution

Holds the data for the solution to a finite element problem.

### Fields

* `femMesh::FEMmesh`: The finite element mesh the problem was solved on.
* `u::Array{Float64}`: The solution (at the final timepoint)
* `trueKnown::Bool`: Boolean flag for if the true solution is given.
* `uTrue::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `l2Err::NumberOrVoid`: The L2 error between u and uTrue.
* `h1Err::NumberOrVoid`: The H1 error between u and uTrue.
* `maxErr::NumberOrVoid`: The nodal maximum error between u and uTrue.
* `nodeErr2::NumberOrVoid`: The nodal l2 error between y abd uTrue.
* `appxTrue::Bool`: Boolean flag for if uTrue was an approximation.
* `uFull`::AbstractArrayOrVoid`: u over time. Only saved if `fullSave=true`
is specified in the solver.
* `tFull::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `fullSave=true`
is specified in the solver.
* `fullSave::Bool`: True if solver saved the extra timepoints.

"""
type FEMSolution <: DESolution
  femMesh::FEMmesh
  u::Array{Float64}
  trueKnown::Bool
  uTrue::AbstractArrayOrVoid
  errors#::Dict{String,Float64}
  appxTrue::Bool
  uFull::AbstractArrayOrVoid
  tFull::AbstractArrayOrVoid
  fullSave::Bool
  function FEMSolution(femMesh::FEMmesh,u,uTrue,sol,Du,uFull,tFull;fullSave=true)
    errors = Dict("L2"=>getL2error(femMesh,sol,u),"H1"=>getH1error(femMesh,Du,u),
                  "l∞"=> maximum(abs(u-uTrue)), "l2"=> norm(u-uTrue,2))
    return(new(femMesh,u,true,uTrue,errors,false,uFull,tFull,true))
  end
  FEMSolution(femMesh,u,uTrue,sol,Du) = FEMSolution(femMesh::FEMmesh,u,uTrue,sol,Du,nothing,nothing,fullSave=false)
  function FEMSolution(femMesh::FEMmesh,u::AbstractArray)
    return(FEMSolution(femMesh,u,nothing,nothing,fullSave=true))
  end
  function FEMSolution(femMesh::FEMmesh,u::AbstractArray,uFull,tFull;fullSave=true)
    return(new(femMesh,u,false,nothing,Dict{String,Float64},false,uFull,tFull,fullSave))
  end
end

type SDESolution <: DESolution
  u::Float64
  trueKnown::Bool
  uTrue::Float64
  errors#::Dict{}
  uFull::AbstractArrayOrVoid
  tFull::AbstractArrayOrVoid
  WFull::AbstractArrayOrVoid
  solFull::AbstractArrayOrVoid
  appxTrue::Bool
  fullSave::Bool
  function SDESolution(u;uFull=nothing,solFull=nothing,tFull=nothing,WFull=nothing)
    fullSave = uFull == nothing
    trueKnown = false
    return(new(u,trueKnown,nothing,Dict(),uFull,solFull,tFull,WFull,false,fullSave))
  end
  function SDESolution(u,uTrue;uFull=nothing,solFull=nothing,tFull=nothing,WFull=nothing)
    fullSave = uFull != nothing
    trueKnown = true
    errors = Dict("final"=>abs(u-uTrue))
    if fullSave
      errors = Dict("final"=>abs(u-uTrue),"l∞"=>maximum(abs(uFull-solFull)),"l2"=>sqrt(mean((uFull-solFull).^2)))
    end
    return(new(u,trueKnown,uTrue,errors,uFull,solFull,tFull,WFull,false,fullSave))
  end
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
  res.errors = Dict("l∞"=>maximum(abs(res.u-res.uTrue)),"l2"=>norm(res.u-res.uTrue,2))
  res.appxTrue = true
end
