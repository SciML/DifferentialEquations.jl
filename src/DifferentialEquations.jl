__precompile__()

module DifferentialEquations

  using IterativeSolvers, Parameters, Plots, GenericSVD, ForwardDiff,
        ChunkedArrays, InplaceOps, SIUnits, Sundials
  import Base: length, size, getindex, endof, show, print, max, linspace
  import Plots: plot
  import ForwardDiff.Dual

  "`DEProblem`: Defines differential equation problems via its internal functions"
  abstract DEProblem
  abstract AbstractODEProblem <: DEProblem
  abstract AbstractSDEProblem <: DEProblem
  abstract AbstractDAEProblem <: DEProblem
  abstract AbstractDDEProblem <: DEProblem
  "`PdeSolution`: Wrapper for the objects obtained from a solver"
  abstract DESolution
  abstract AbstractODESolution <: DESolution
  "`Mesh`: An abstract type which holds a (node,elem) pair and other information for a mesh"
  abstract Mesh
  "`Tableau`: Holds the information for a Runge-Kutta Tableau"
  abstract Tableau
  "`DEIntegrator`: A DifferentialEquations Integrator type, used to initiate a solver."
  abstract DEIntegrator
  "`DEParameters`: Holds the parameters used in a DifferntialEquations model"
  abstract DEParameters
  "`ODERKTableau`: A Runge-Kutta Tableau for an ODE integrator"
  abstract ODERKTableau <: Tableau
  typealias KW Dict{Symbol,Any}
  AbstractArrayOrVoid = Union{AbstractArray,Void}
  NumberOrVoid = Union{Number,Void}
  FunctionOrVoid = Union{Function,Void}

  #Constants

  const TEST_FLOPS_CUTOFF = 1e10
  const initialized_backends = Set{Symbol}()

  include("general/units.jl")
  include("general/backends.jl")
  include("general/replacement_macros.jl")
  include("general/misc_utils.jl")
  include("fem/fem_heatintegrators.jl")
  include("fem/meshes.jl")
  include("fem/fem_assembly.jl")
  include("fem/fem_boundary.jl")
  include("fem/fem_error.jl")
  include("sde/sde_noise_process.jl")
  include("general/problems.jl")
  include("general/solutions.jl")
  include("general/stochastic_utils.jl")
  include("general/convergence.jl")
  include("premades/premade_problems.jl")
  include("premades/premade_meshes.jl")
  include("fem/fem_solve.jl")
  include("fdm/stokes_solve.jl")
  include("sde/sde_solve.jl")
  include("sde/sde_integrators.jl")
  include("ode/ode_tableaus.jl")
  include("ode/ode_integrators.jl")
  include("ode/ode_constants.jl")
  include("ode/ode_callbacks.jl")
  include("ode/ode_solve.jl")
  include("ode/ode_dense.jl")
  include("dae/dae_solve.jl")
  include("sde/sde_tableaus.jl")
  include("general/benchmark.jl")
  include("general/plotrecipes.jl")

  #Types
  export DEProblem, DESolution, DEParameters, HeatProblem, PoissonProblem, FEMSolution, Mesh,
         ConvergenceSimulation, FEMmesh, SimpleMesh, SDEProblem, StokesProblem, DAEProblem,
         DAESolution, SDESolution, ODESolution, ODEProblem, FDMMesh, ExplicitRKTableau,
         MonteCarloSimulation,ImplicitRKTableau, Shootout, ShootoutSet,AbstractODEProblem,
         AbstractSDEProblem, TestSolution

  #SDE Example Problems
  export prob_sde_wave, prob_sde_linear, prob_sde_cubic, prob_sde_2Dlinear, prob_sde_lorenz,
         prob_sde_2Dlinear, prob_sde_additive, prob_sde_additivesystem

  #ODE Example Problems
  export prob_ode_linear, prob_ode_bigfloatlinear, prob_ode_2Dlinear,
         prob_ode_large2Dlinear, prob_ode_bigfloat2Dlinear, prob_ode_rigidbody,
         prob_ode_2Dlinear_notinplace, prob_ode_vanderpol, prob_ode_vanderpol_stiff,
         prob_ode_lorenz, prob_ode_rober, prob_ode_threebody

  #FEM Example Problems
  export  prob_femheat_moving, prob_femheat_pure, prob_femheat_diffuse,
          prob_poisson_wave, prob_poisson_noisywave, prob_femheat_birthdeath,
          prob_poisson_birthdeath, prob_femheat_stochasticbirthdeath,
          prob_stokes_homogenous, prob_stokes_dirichletzero, prob_poisson_birthdeathsystem,
          prob_poisson_birthdeathinteractingsystem, prob_femheat_birthdeathinteractingsystem,
          prob_femheat_birthdeathsystem, prob_femheat_diffusionconstants,
          heatProblemExample_grayscott,heatProblemExample_gierermeinhardt

  #DAE Example Problems
  export prob_dae_resrob

  #Example Meshes
  export  meshExample_bunny, meshExample_flowpastcylindermesh, meshExample_lakemesh,
          meshExample_Lshapemesh, meshExample_Lshapeunstructure, meshExample_oilpump,
          meshExample_wavymesh, meshExample_wavyperturbmesh

  #Benchmark Functions
  export ode_shootout, ode_shootoutset, ode_workprecision, ode_workprecision_set

  #Plot Functions
  export  plot, animate, stability_region

  #General Functions
  export appxtrue!, accumarray, solve, test_convergence

  #Stochastic Utils
  export monteCarloSim, construct_correlated_noisefunc, WHITE_NOISE, NoiseProcess

  #FEM Functions
  export  assemblematrix, findboundary, setboundary, findbdtype, getL2error, quadpts, getH1error,
          gradu, gradbasis, quadfbasis, fem_squaremesh, CFLμ, CFLν,
          meshgrid, notime_squaremesh, parabolic_squaremesh, quadpts1

  #Tableus
  export constructEuler, constructKutta3, constructRK4, constructRK438Rule,
         constructImplicitEuler, constructMidpointRule, constructTrapezoidalRule,
         constructLobattoIIIA4, constructLobattoIIIB2, constructLobattoIIIB4,
         constructLobattoIIIC2, constructLobattoIIIC4, constructLobattoIIICStar2,
         constructLobattoIIICStar4, constructLobattoIIID2, constructLobattoIIID4,
         constructRadauIA3, constructRadauIA5, constructRadauIIA3, constructRadauIIA5,
         constructRalston, constructHeun, constructRKF5, constructBogakiShampine3,
         constructCashKarp, constructDormandPrince, constructRKF8, constructDormandPrince8,
         constructMSRI1,constructFeagin10, constructFeagin12, constructFeagin14,
         constructDormandPrince8_64bit, constructDP8, constructRKF5, constructRungeFirst5,
         constructCassity5, constructLawson5, constructLutherKonen5, constructLutherKonen52,
         constructLutherKonen53, constructPapakostasPapaGeorgiou5, constructPapakostasPapaGeorgiou52,
         constructTsitouras5, constructBogakiShampine5, constructSharpSmart5,
         constructButcher6, constructButcher7, constructDverk, constructClassicVerner6,
         constructClassicVerner7, constructClassicVerner8, constructClassicVerner92,
         constructVernerRobust7, constructEnrightVerner7, constructTanakaYamashitaStable7,
         constructTanakaYamashitaEfficient7, constructSharpSmart7, constructSharpVerner7,
         constructVernerEfficient7, constructCooperVerner8, constructCooperVerner82,
         constructTsitourasPapakostas8, constructdverk78, constructEnrightVerner8,
         constructCurtis8, constructVernerRobust9, constructVernerEfficient9,
         constructSharp9, constructTsitouras9, constructTsitouras92,constructFeagin14Tableau,
         constructFeagin12Tableau, constructOno12, constructCurtis10, constructOno10, constructFeagin10Tableau,
         constructCurtis10, constructBaker10, constructHairer10, constructButcher63,
         constructButcher6, constructButcher62, constructVerner6, constructDormandPrince6,
         constructSharpVerner6, constructVerner9162, constructVerner916, constructVernerRobust6,
         constructVernerEfficient6, constructPapakostas6, constructLawson6,
         constructTsitourasPapakostas6, constructDormandLockyerMcCorriganPrince6,
         constructTanakaKasugaYamashitaYazaki6D, constructTanakaKasugaYamashitaYazaki6C,
         constructTanakaKasugaYamashitaYazaki6B, constructTanakaKasugaYamashitaYazaki6A,
         constructMikkawyEisa, constructChummund6, constructChummund62,
         constructHuta62, constructHuta6, constructRKF4, constructVern6, constructTsit5,
         constructBS5, constructTanYam7, constructTsitPap8, constructVern9,
         constructVerner8, constructVerner7

  #Misc Tools
  export numparameters, checkSRIOrder, checkSRAOrder, @ode_savevalues,
         constructSRIW1, constructSRA1, def, @ode_define, @fem_define

  #Callback Necessary
  export ode_addsteps!, ode_interpolant, DIFFERENTIALEQUATIONSJL_SPECIALDENSEALGS,
         @ode_callback, @ode_event, @ode_change_cachesize, @ode_change_deleteat,
         @ode_terminate, copyat_or_push!

  #=
   include("precompile.jl")
   __precompile__()
   =#
end # module
