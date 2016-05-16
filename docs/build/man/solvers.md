
<a id='Information-on-Solvers-1'></a>

# Information on Solvers


<a id='Finite-Element-Method-Solvers-1'></a>

## Finite Element Method Solvers

<a id='DifferentialEquations.fem_solvepoisson' href='#DifferentialEquations.fem_solvepoisson'>#</a>
**`DifferentialEquations.fem_solvepoisson`** &mdash; *Function*.



**Finite Element Poisson Equation Solver**

fem_solvepoisson(femMesh::FEMmesh,pdeProb::PoissonProblem)

Takes in a definition for the heat equation `-Δu = f` on `femMesh` with functions as defined in `pdeProb`. If `σ` is specified in `pdeProb`, then this solves the stochastic Poisson equation `-Δu = f + σdW`.

**Keyword Arguments**

  * `solver` = Linear solver algorithm. This is the algorithm which is chosen for solving the implicit equation `Ax=b`. The default is `LU`. The choices are:
  * `Direct` = Solves `Ax=b` using ``
  * `CG` = Conjugate-Gradient. Best when the space is very large and `I ± ΔtM⁻¹A` is positive definite.
  * `GMRES` = GMRES. Best when the space is very large and `I ± ΔtM⁻¹A` is not positive definite.
  * `saveSteps` = If `fullSave=true`, then this is the number of steps between the saves.
  * `autodiff` = Whether or not autodifferentiation (as provided by AutoDiff.jl) is used for the nonlinear solving. By default autodiff is false.

<a id='DifferentialEquations.fem_solveheat' href='#DifferentialEquations.fem_solveheat'>#</a>
**`DifferentialEquations.fem_solveheat`** &mdash; *Function*.



**Finite Element Heat Equation Solver**

`fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem)`

Takes in a definition for the heat equation `u_t = Δu + f` on `femMesh` with functions as defined in `pdeProb`. If `σ` is specified in `pdeProb`, then this solves the stochastic heat equation `u_t = Δu + f + σdW_t`.

**Keyword Arguments**

  * `alg` = Solution algorithm. Default is Euler. The choices are:
  * Linear     * Euler (Explicit)     * Implicit Euler (Implicit)     * Crank-Nicholson (Implicit)
  * Nonlinear     * Euler (Explicit)     * Implicit Euler (Nonlinear Solve)     * Crank-Nicholson (Nonlinear Solve)     * Semi-Implicit Euler (Implicit)     * Semi-Implicit Crank-Nicholson (Implicit)

Explicit algorithms only require solving matrix multiplications `Au`. Implicit algorithms require solving the linear equation `Ax=b` where `x` is the unknown. Nonlinear Solve algorithms require solving the nonlinear equation f(x)=0 using methods like Newton's method and is provided by NLSolve.jl. Explicit algorithms have the least stability and should be used either small Δt and non-stiff equations. The implicit algorithms have better stability, but for nonlinear equations require costly nonlinear solves in order to be solved exactly. The semi-implicit algorithms discretize with part of the equation implicit and another part explicit in order to allow for the algorithm to not require a nonlinear solve, but at the cost of some stability (though still vastly better at stability than explicit algorithms).

  * `solver` = Linear solver algorithm. This is the algorithm which is chosen for solving the implicit equation `Ax=b`. The default is `LU`. The choices are:
  * `Direct` = Solves using `` (no factorization). Not recommended.
  * `Cholesky` = Cholsky decomposition. Only stable of `I ± ΔtM⁻¹A` is positive definite.     This means that this works best when Δt is small. When applicable, this is the fastest.
  * `LU` = LU-Decomposition. A good mix between fast and stable.
  * `QR` = QR-Decomposition. Less numerical roundoff error than `LU`, but slightly slower.
  * `SVD` = SVD-Decomposition. By far the slowest, but the most robust to roundoff error.
  * `CG` = Conjugate-Gradient. Best when the space is very large and `I ± ΔtM⁻¹A` is positive definite.
  * `GMRES` = GMRES. Best when the space is very large and `I ± ΔtM⁻¹A` is not positive definite.
  * `fullSave` = Makes the algorithm save the output at every `saveSteps` timesteps. By default fullSave is false.
  * `saveSteps` = If `fullSave=true`, then this is the number of steps between the saves.
  * `autodiff` = Whether or not autodifferentiation (as provided by AutoDiff.jl) is used for the nonlinear solving. By default autodiff is false.

