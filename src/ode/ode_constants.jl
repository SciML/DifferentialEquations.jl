"""
ODE_DEFAULT_TABLEAU

Sets the default tableau for the ODE solver. Currently Dormand-Prince 4/5.
"""
const ODE_DEFAULT_TABLEAU = constructDormandPrince()

const DIFFERENTIALEQUATIONSJL_ALGORITHMS = Set([:Euler,:Midpoint,:RK4,:ExplicitRK,:ImplicitEuler,:Trapezoid,:Rosenbrock32,:Feagin10])
const ODEINTERFACE_ALGORITHMS = Set([:dopri5,:dop853,:odex,:radau5,:radau,:seulex])
const ODEJL_ALGORITHMS = Set([:ode23,:ode45,:ode78,:ode23s,:ode1,:ode2_midpoint,:ode2_heun,:ode4,:ode45_fe])

const DIFFERENTIALEQUATIONSJL_DEFAULT_OPTIONS = Dict(:Δt => 0,
                                 :save_timeseries => true,
                                 :timeseries_steps => 1,
                                 :tableau => DifferentialEquations.ODE_DEFAULT_TABLEAU,
                                 :adaptive => true,
                                 :γ=>2.0,
                                 :abstol=>1//10^3,
                                 :reltol=>1//10^6,
                                 :qmax=>10,
                                 :maxiters => round(Int,1e9),
                                 :Δtmax=>nothing,
                                 :Δtmin=>nothing,
                                 :autodiff=>true,
                                 :internalnorm => 2,
                                 :progressbar=>false,
                                 :progress_steps=>1000)

const ODEJL_OPTION_LIST = Set([:tout,:tstop,:reltol,:abstol,:minstep,:maxstep,:initstep,:norm,:maxiters,:isoutofdomain])
const ODEJL_ALIASES = Dict{Symbol,Symbol}(:minstep=>:Δtmin,:maxstep=>:Δtmax,:initstep=>:Δt,:tstop=>:T)
const ODEJL_ALIASES_REVERSED = Dict{Symbol,Symbol}([(v,k) for (k,v) in ODEJL_ALIASES])
const ODEINTERFACE_OPTION_LIST = Set([:RTOL,:ATOL,:OUTPUTFCN,:OUTPUTMODE,
                                :MAXSTEPS,:STEST,:EPS,:RHO,:SSMINSEL,
                                :SSMAXSEL,:SSBETA,:MAXSS,:INITIALSS,
                                :MAXEXCOLUMN,:STEPSIZESEQUENCE,:MAXSTABCHECKS,
                                :MAXSTABCHECKLINE,:DENSEOUTPUTWOEE,:INTERPOLDEGRE,
                                :SSREDUCTION,:SSSELECTPAR1,:SSSELECTPAR2,
                                :ORDERDECFRAC,:ORDERINCFRAC,:OPT_RHO,:OPT_RHO2,
                                :RHSAUTONOMOUS,:M1,:M2,:LAMBDADENSE,:TRANSJTOH,
                                :STEPSIZESEQUENCE,:JACRECOMPFACTOR,:MASSMATRIX,
                                :JACOBIMATRIX,:JACOBIBANDSSTRUCT,:WORKFORRHS,
                                :WORKFORJAC,:WORKFORDEC,:WORKFORSOL,:MAXNEWTONITER,
                                :NEWTONSTARTZERO,:NEWTONSTOPCRIT,:DIMFIND1VAR,
                                :MAXSTAGES,:MINSTAGES,:INITSTAGES,:STEPSIZESTRATEGY,
                                :FREEZESSLEFT,:FREEZESSRIGHT,:ORDERDECFACTOR,
                                :ORDERINCFACTOR,:ORDERDECCSTEPFAC1,:ORDERDECSTEPFAC2
                                ])

const ODEINTERFACE_ALIASES = Dict{Symbol,Symbol}(:RTOL=>:reltol,
                                                 :ATOL=>:abstol,
                                                 :MAXSTEPS=> :maxiters,
                                                 :MAXSS=>:Δtmax,
                                                 :INITIALSS=>:Δt,
                                                 #:SSMINSEL=>:qmin,
                                                 :SSMAXSEL=>:qmax)
const ODEINTERFACE_ALIASES_REVERSED = Dict{Symbol,Symbol}([(v,k) for (k,v) in ODEINTERFACE_ALIASES])

const DIFFERENTIALEQUATIONSJL_ORDERS = Dict{Symbol,Int}(:Euler=>1,
                                                        :Midpoint=>2,
                                                        :RK4=>4,
                                                        :ExplicitRK=>4, #Gets overwritten
                                                        :ImplicitEuler=>1,
                                                        :Trapezoid=>2,
                                                        :Rosenbrock32=>3,
                                                        :Feagin10=>10)
const DIFFERENTIALEQUATIONSJL_ADAPTIVEALGS = Set([:ExplicitRK,:Rosenbrock32,:Feagin10])
const DIFFERENTIALEQUATIONSJL_IMPLICITALGS = Set([:ImplicitEuler,:Trapezoid,:Rosenbrock32])
const ODEINTERFACE_STRINGS = Dict{Symbol,String}(
  :LOGIO            => "logio",
  :LOGLEVEL         => "loglevel",
  :RHS_CALLMODE     => "RightHandSideCallMode",

  :RTOL             => "RelTol",
  :ATOL             => "AbsTol",
  :MAXSTEPS         => "MaxNumberOfSteps",
  :EPS              => "eps",

  :OUTPUTFCN        => "OutputFcn",
  :OUTPUTMODE       => "OutputFcnMode",

  :STEST            => "StiffTestAfterStep",
  :RHO              => "rho",
  :SSMINSEL         => "StepSizeMinSelection",
  :SSMAXSEL         => "StepSizeMaxSelection",
  :SSBETA           => "StepSizeBeta",
  :MAXSS            => "MaxStep",
  :INITIALSS        => "InitialStep",


  :MAXEXCOLUMN      => "MaxExtrapolationColumn",
  :MAXSTABCHECKS    => "MaxNumberOfStabilityChecks",
  :MAXSTABCHECKLINE => "MaxLineForStabilityCheck",
  :INTERPOLDEGREE   => "DegreeOfInterpolation",
  :ORDERDECFRAC     => "OrderDecreaseFraction",
  :ORDERINCFRAC     => "OrderIncreaseFraction",
  :STEPSIZESEQUENCE => "StepSizeSequence",
  :SSREDUCTION      => "StepSizeReduction",
  :SSSELECTPAR1     => "StepSizeSelectionParam1",
  :SSSELECTPAR2     => "StepSizeSelectionParam2",
  :RHO2             => "rho2",
  :DENSEOUTPUTWOEE  => "DeactivateErrorEstInDenseOutput",

  :TRANSJTOH        => "TransfromJACtoHess",
  :MAXNEWTONITER    => "MaxNewtonIterations",
  :NEWTONSTARTZERO  => "StartNewtonWithZeros",
  :DIMOFIND1VAR     => "DimensionOfIndex1Vars",
  :DIMOFIND2VAR     => "DimensionOfIndex2Vars",
  :DIMOFIND3VAR     => "DimensionOfIndex3Vars",
  :STEPSIZESTRATEGY => "StepSizeStrategy",
  :M1               => "M1",
  :M2               => "M2",
  :JACRECOMPFACTOR  => "RecomputeJACFactor",
  :NEWTONSTOPCRIT   => "NewtonStopCriterion",
  :FREEZESSLEFT     => "FreezeStepSizeLeftBound",
  :FREEZESSRIGHT    => "FreezeStepSizeRightBound",
  :MASSMATRIX       => "MassMatrix",
  :JACOBIMATRIX     => "JacobiMatrix",
  :JACOBIBANDSTRUCT => "JacobiBandStructure",

  :MAXSTAGES        => "MaximalNumberOfStages",
  :MINSTAGES        => "MinimalNumberOfStages",
  :INITSTAGES       => "InitialNumberOfStages",
  :ORDERINCFACTOR   => "OrderIncreaseFactor",
  :ORDERDECFACTOR   => "OrderDecreaseFactor",
  :ORDERDECSTEPFAC1 => "OrderDecreaseStepFactor1",
  :ORDERDECSTEPFAC2 => "OrderDecreaseStepFactor2",

  :RHSAUTONOMOUS    => "AutonomousRHS",
  :LAMBDADENSE      => "LambdaForDenseOutput",
  :WORKFORRHS       => "WorkForRightHandSide",
  :WORKFORJAC       => "WorkForJacobimatrix",
  :WORKFORDEC       => "WorkForLuDecomposition",
  :WORKFORSOL       => "WorkForSubstitution",

  :BVPCLASS         => "BoundaryValueProblemClass",
  :SOLMETHOD        => "SolutionMethod",
  :IVPOPT           => "OptionsForIVPsolver")
