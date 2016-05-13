__precompile__()

@doc """
###DifferentialEquations

This is a package for solving numerically solving differential equations in Julia
by Chris Rackauckas. The purpose of this package is to supply efficient Julia
implementations of solvers for various differential equations. Equations within
the realm of this package include stochastic ordinary differential equations
(SODEs or SDEs), stochastic partial differential equations (SPDEs), partial
differential equations (with both finite difference and finite element methods),
and differential delay equations. For ordinary differential equation solvers,
see [ODE.jl](https://github.com/JuliaLang/ODE.jl)

This package is for efficient and parallel implementations of research-level
algorithms, many of which are quite recent. These algorithms aim to be optimized
for HPC applications, including the use of GPUs, Xeon Phis, and multi-node
parallelism. With the easy to use plot/convergence testing algorithms,
this package also provides a good sandbox for developing novel numerical schemes.
""" ->
module DifferentialEquations

import PyPlot
using LaTeXStrings, Plots, Atom, IterativeSolvers, NLsolve
import Base: length

"PdeProblem: Defines PDE problems via its internal functions"
abstract PdeProblem
"PdeSolution: Wrapper for the objects obtained from a PdeSolver"
abstract PdeSolution

AbstractArrayOrVoid = Union{AbstractArray,Void}
NumberOrVoid = Union{Number,Void}
FunctionOrVoid = Union{Function,Void}

include("fem/meshTools.jl")
include("fem/assemblyTools.jl")
include("fem/boundaryTools.jl")
include("fem/errorTools.jl")
include("general/solutionTools.jl")
include("general/plotTools.jl")
include("general/problemTools.jl")
include("general/stochasticTools.jl")
include("general/miscTools.jl")
include("examples/exampleProblems.jl")
include("fem/femSolvers.jl")

#Types
export PdeProblem, PdeSolution, HeatProblem, PoissonProblem, FEMSolution, ConvergenceSimulation, FEMmesh

#Example Problems
export  heatProblemExample_moving, heatProblemExample_diffuse, heatProblemExample_pure,
        poissonProblemExample_wave, poissonProblemExample_noisyWave, heatProblemExample_birthdeath,
        poissonProblemExample_birthdeath, heatProblemExample_stochasticbirthdeath

#Plot Functions
export  solplot_animation, solplot_appxvstrue, solplot_appx, showmesh, convplot,
        convplot_fullΔt, convplot_fullΔx, convplot_fullΔx, convplot_l2vsΔt, convplot_node2vsΔt,
        convplot_maxvsΔt, convplot_h1vsΔx, convplot_l2vsΔx, convplot_node2vsΔx, convplot_maxvsΔx,
        solplot

#General Functions
export conv_ests, appxTrue!, accumarray

#FEM Functions
export  assemblematrix, findboundary, setboundary, findbdtype, getL2error, quadpts, getH1error,
        gradu, gradbasis, fem_solvepoisson, fem_solveheat, quadfbasis, fem_squaremesh, CFLμ, CFLν,
        meshgrid, notime_squaremesh, parabolic_squaremesh, quadpts1

#Misc Tools
export quadfbasis2, CG2
end # module
