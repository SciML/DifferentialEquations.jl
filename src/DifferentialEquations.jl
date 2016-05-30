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

using LaTeXStrings, Plots, IterativeSolvers, NLsolve, Parameters, Compat, EllipsesNotation
import Plots: plot
import Base: length, size
import JLD: load
using GrowableArrays, CheckImports, ChunkedArrays

"PdeProblem: Defines differential equation problems via its internal functions"
abstract DEProblem
"PdeSolution: Wrapper for the objects obtained from a solver"
abstract DESolution
"Mesh: An abstract type which holds a (node,elem) pair and other information for a mesh"
abstract Mesh
typealias String AbstractString
AbstractArrayOrVoid = Union{AbstractArray,Void}
NumberOrVoid = Union{Number,Void}
FunctionOrVoid = Union{Function,Void}

#Constants

const TEST_FLOPS_CUTOFF = 1e10

include("fem/meshTools.jl")
include("fem/assemblyTools.jl")
include("fem/boundaryTools.jl")
include("fem/errorTools.jl")
include("general/problemTools.jl")
include("general/solutionTools.jl")
include("general/stochasticTools.jl")
include("general/miscTools.jl")
include("general/convergenceTools.jl")
include("examples/exampleProblems.jl")
include("examples/exampleMeshes.jl")
include("fem/femSolvers.jl")
include("fdm/stokesSolvers.jl")
include("sde/sdeSolvers.jl")
include("general/plotTools.jl")
include("general/coefficientTypes.jl")
include("general/parallelHelpers.jl")

#Types
export DEProblem, DESolution, HeatProblem, PoissonProblem, FEMSolution, Mesh,
       ConvergenceSimulation, FEMmesh, SimpleMesh, SDEProblem, StokesProblem,
       SDESolution, FDMMesh

#SDE Example Problems
export linearSDEExample, cubicSDEExample, waveSDEExample, additiveSDEExample,
       multiDimAdditiveSDEExample, twoDimlinearSDEExample

#FEM Example Problems
export  heatProblemExample_moving, heatProblemExample_diffuse, heatProblemExample_pure,
        poissonProblemExample_wave, poissonProblemExample_noisyWave, heatProblemExample_birthdeath,
        poissonProblemExample_birthdeath, heatProblemExample_stochasticbirthdeath,
        homogeneousStokesExample, dirichletzeroStokesExample

#Example Meshes
export  meshExample_bunny, meshExample_flowpastcylindermesh, meshExample_lakemesh,
        meshExample_Lshapemesh, meshExample_Lshapeunstructure, meshExample_oilpump,
        meshExample_wavymesh, meshExample_wavyperturbmesh

#Plot Functions
export  plot, showmesh, convplot, solplot_animation,
        convplot_fullΔt, convplot_fullΔx, convplot_l2vsΔt, convplot_node2vsΔt,
        convplot_maxvsΔt, convplot_h1vsΔx, convplot_l2vsΔx, convplot_node2vsΔx, convplot_maxvsΔx,
        solplot

#General Functions
export conv_ests, appxTrue!, accumarray, solve, testConvergence, monteCarloSim

#FEM Functions
export  assemblematrix, findboundary, setboundary, findbdtype, getL2error, quadpts, getH1error,
        gradu, gradbasis, quadfbasis, fem_squaremesh, CFLμ, CFLν,
        meshgrid, notime_squaremesh, parabolic_squaremesh, quadpts1

#Misc Tools
export quadfbasis2, CG2, numparameters, checkSRIOrder, checkSRAOrder,
       constructSRIW1, constructSRA1, def

end # module
