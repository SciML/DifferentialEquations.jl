"""
FEMSolution

A type which holds the data for the solution to a finite element problem.
"""
type FEMSolution
  femMesh::FEMmesh
  u::Array{Float64}
  trueKnown::Bool
  uTrue::AbstractArray
  l2Err::Float64
  h1Err::Float64
  maxErr::Float64
  nodeErr2::Float64
  appxTrue::Bool
  uFull::AbstractArray
  tFull::AbstractArray
  fullSave::Bool
  function FEMSolution(femMesh::FEMmesh,u,uTrue,sol,Du,uFull,tFull;fullSave=true)
    l2Err = getL2error(femMesh,sol,u)
    h1Err = getH1error(femMesh,Du,u)
    maxErr = maximum(abs(u-uTrue))
    nodeErr2 = norm(u-uTrue,2)
    return(new(femMesh,u,true,uTrue,l2Err,h1Err,maxErr,nodeErr2,uFull,tFull,true))
  end
  FEMSolution(femMesh,u,uTrue,sol,Du) = FEMSolution(femMesh::FEMmesh,u,uTrue,sol,Du,nothing,nothing,fullSave=false)
  function FEMSolution(femMesh::FEMmesh,u::AbstractArray)
    return(FEMSolution(femMesh,u,nothing))
  end
  function FEMSolution(femMesh::FEMmesh,u::AbstractArray,uFull,tFull)
    return(new(femMesh,u,false,nothing,nothing,nothing,nothing,nothing,uFull,tFull))
  end
end

"""
ConvergenceSimulation

A type which holds the data from a convergence simulation.
"""
type ConvergenceSimulation
  solutions
  h1Errors
  l2Errors
  maxErrors
  node2Errors
  N
  Δts
  Δxs
  μs
  νs
  ConvEst_h1
  ConvEst_l2
  ConvEst_max
  ConvEst_node2
  function ConvergenceSimulation(solutions::Array{FEMSolution})
    N = length(solutions)
    Δts = Vector{Float64}(N)
    Δxs = Vector{Float64}(N)
    μs = Vector{Float64}(N)
    νs = Vector{Float64}(N)
    for i = 1:N
      Δts[i] = solutions[i].femMesh.Δt
      Δxs[i] = solutions[i].femMesh.Δx
      μs[i] = solutions[i].femMesh.μ
      νs[i] = solutions[i].femMesh.ν
    end
    if solutions[1].trueKnown
      h1errors = Vector{Float64}(N)
      l2errors = Vector{Float64}(N)
      maxErrs  = Vector{Float64}(N)
      nodeErr2s = Vector{Float64}(N)
      for i = 1:N
        h1errors[i] = solutions[i].h1Err
        l2errors[i] = solutions[i].l2Err
        maxErrs[i]  = solutions[i].maxErr
        nodeErr2s[i] = solutions[i].nodeErr2
      end
      h1Est    = conv_ests(h1errors)
      l2Est    = conv_ests(l2errors)
      maxEst   = conv_ests(maxErrs)
      node2Est = conv_ests(nodeErr2s)
      return(new(solutions,h1errors,l2errors,maxErrs,nodeErr2s,N,Δts,Δxs,μs,
                 νs,h1Est,l2Est,maxEst,node2Est))
    else # No known solution
      if solutions[1].appxTrue # But appx true solution known
        maxErrs  = Vector{Float64}(N)
        nodeErr2s = Vector{Float64}(N)
        for i = 1:N
          maxErrs[i]  = solutions[i].maxErr
          nodeErr2s[i] = solutions[i].nodeErr2
        end
        maxEst   = conv_ests(maxErrs)
        node2Est = conv_ests(nodeErr2s)
        return(new(solutions,nothing,nothing,maxErrs,nodeErr2s,N,Δts,Δxs,μs,νs,nothing,nothing,maxEst,node2Est))
      else #Nothing true known
        return(new(solutions,nothing,nothing,nothing,nothing,N,Δts,Δxs,μs,νs))
      end
    end

  end
end

"""
length(simres::ConvergenceSimulation)

Returns the number of simultations in the Convergence Simulation
"""
Base.length(simres::ConvergenceSimulation) = simres.N

"""
conv_ests(error::Vector{Number})

Computes the pairwise convergence estimate for a convergence test done by
halving/doubling stepsizes via

log2(error[i+1]/error[i])

Returns the mean of the convergence estimates
"""
function conv_ests(error::Vector{Float64})
  ## Calculate pairwise rate of convergence estimates
  S = Vector{Float64}(length(error)-1)
  for i=1:length(error)-1
    S[i] = log2(error[i+1]/error[i])
  end
  return(abs(mean(S)))
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
  res.maxErr = maximum(abs(res.u-res.uTrue))
  res.nodeErr2 = norm(res.u-res.uTrue,2)
  res.appxTrue = true
end
