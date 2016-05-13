__precompile__()

@doc """
###DifferentialEquations

Developed by Chris Rackauckas for the solution of various differential equations, including stochastic differential equations (SDEs),
partial differential equations (PDEs), and stochastic partial differential equations (SPDEs). Included are tools for finite difference
and finite element methods, including functions for mesh generation and plotting.
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
        convplot_maxvsΔt, convplot_h1vsΔx, convplot_l2vsΔx, convplot_node2vsΔx, convplot_maxvsΔx

#General Functions
export conv_ests, appxTrue!, accumarray

#FEM Functions
export  assemblematrix, findboundary, setboundary, findbdtype, getL2error, quadpts, getH1error,
        gradu, gradbasis, fem_solvepoisson, fem_solveheat, quadfbasis, fem_squaremesh, CFLμ, CFLν,
        meshgrid, notime_squaremesh, parabolic_squaremesh, quadpts1

end # module
