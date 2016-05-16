
<a id='The-Solution-Type-1'></a>

# The Solution Type


<a id='Related-Functions-1'></a>

## Related Functions

<a id='DifferentialEquations.FEMSolution' href='#DifferentialEquations.FEMSolution'>#</a>
**`DifferentialEquations.FEMSolution`** &mdash; *Type*.



FEMSolution

Holds the data for the solution to a finite element problem.

**Fields**

  * `femMesh::FEMmesh`: The finite element mesh the problem was solved on.
  * `u::Array{Float64}`: The solution (at the final timepoint)
  * `trueKnown::Bool`: Boolean flag for if the true solution is given.
  * `uTrue::AbstractArrayOrVoid`: The true solution at the final timepoint.
  * `l2Err::NumberOrVoid`: The L2 error between u and uTrue.
  * `h1Err::NumberOrVoid`: The H1 error between u and uTrue.
  * `maxErr::NumberOrVoid`: The nodal maximum error between u and uTrue.
  * `nodeErr2::NumberOrVoid`: The nodal l2 error between y abd uTrue.
  * `appxTrue::Bool`: Boolean flag for if uTrue was an approximation.
  * `uFull`::AbstractArrayOrVoid`: u over time. Only saved if `fullSave=true` is specified in the solver.
  * `tFull::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `fullSave=true` is specified in the solver.
  * `fullSave::Bool`: True if solver saved the extra timepoints.

<a id='DifferentialEquations.appxTrue!' href='#DifferentialEquations.appxTrue!'>#</a>
**`DifferentialEquations.appxTrue!`** &mdash; *Function*.



appxTrue!(res,res2)

Adds the solution from res2 to the FEMSolution object res. Useful to add a quasi-true solution when none is known by computing once at a very small time/space step and taking that solution as the "true" solution

<a id='DifferentialEquations.PdeSolution' href='#DifferentialEquations.PdeSolution'>#</a>
**`DifferentialEquations.PdeSolution`** &mdash; *Type*.



PdeSolution: Wrapper for the objects obtained from a PdeSolver

