
<a id='Information-on-Solvers-1'></a>

# Information on Solvers


<a id='Finite-Element-Method-Solvers-1'></a>

## Finite Element Method Solvers

<a id='DifferentialEquations.fem_solvepoisson' href='#DifferentialEquations.fem_solvepoisson'>#</a>
**`DifferentialEquations.fem_solvepoisson`** &mdash; *Function*.



fem_solvepoisson

<a id='DifferentialEquations.fem_solveheat' href='#DifferentialEquations.fem_solveheat'>#</a>
**`DifferentialEquations.fem_solveheat`** &mdash; *Function*.



fem_solveheat

`fem_solveheat(femMesh::FEMmesh,u0::AbstractArray,gD::Function,f::Function,isLinear::Bool)`

`fem_solveheat(femMesh::FEMmesh,u0::Function,gD::Function,f::Function,isLinear::Bool)`

`fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem)`

Takes in a definition for the heat equation `u_t = Δu + f` on a finite element mesh with initial condtion u0 and returns the solution.

