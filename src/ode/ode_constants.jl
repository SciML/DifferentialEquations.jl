"""
ODE_DEFAULT_TABLEAU

Sets the default tableau for the ODE solver. Currently Dormand-Prince 4/5.
"""
const ODE_DEFAULT_TABLEAU = constructDormandPrince()

const DIFFERENTIALEQUATIONSJL_ALGORITHMS = Set([:Euler,:Midpoint,:RK4,:ExplicitRK,:ExplicitRKVectorized,:BS3,:BS3Vectorized,:DP5,:DP5Vectorized,:DP8,:DP8Vectorized,:ImplicitEuler,:Trapezoid,:Rosenbrock32,:Feagin10,:Feagin12,:Feagin14,:Feagin10Vectorized,:Feagin12Vectorized,:Feagin14Vectorized])

const DIFFERENTIALEQUATIONSJL_FASLALGS = Set([:DP5,:DP5Vectorized,:DP8,:DP8Vectorized,:BS3,:BS3Vectorized])
const ODEINTERFACE_ALGORITHMS = Set([:dopri5,:dop853,:odex,:radau5,:radau,:seulex])
const ODEJL_ALGORITHMS = Set([:ode23,:ode45,:ode78,:ode23s,:ode1,:ode2_midpoint,:ode2_heun,:ode4,:ode45_fe])

const DIFFERENTIALEQUATIONSJL_DEFAULT_OPTIONS = Dict(:Δt => 0,
                                 :save_timeseries => true,
                                 :timeseries_steps => 1,
                                 :tableau => DifferentialEquations.ODE_DEFAULT_TABLEAU,
                                 :adaptive => true,
                                 :γ=>.9,
                                 :abstol=>1//10^6,
                                 :reltol=>1//10^3,
                                 :qmax=>10.0,
                                 :qmin=>0.2,
                                 :qoldinit=>1//10^4,
                                 :fullnormalize=>false,
                                 :β=>nothing,
                                 :timechoicealg=>:Lund,
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
                                                 :SSBETA=>:β,
                                                 :SSMAXSEL=>:qmax)
const ODEINTERFACE_ALIASES_REVERSED = Dict{Symbol,Symbol}([(v,k) for (k,v) in ODEINTERFACE_ALIASES])

const DIFFERENTIALEQUATIONSJL_ORDERS = Dict{Symbol,Int}(:Euler=>1,
                                                        :Midpoint=>2,
                                                        :RK4=>4,
                                                        :ExplicitRK=>4, #Gets overwritten
                                                        :ExplicitRKVectorized=>4,#Gets overwritten
                                                        :BS3=>3,
                                                        :BS3Vectorized=>3,
                                                        :DP5=>5,
                                                        :DP5Vectorized=>5,
                                                        :DP8=>8,
                                                        :DP8Vectorized=>8,
                                                        :ImplicitEuler=>1,
                                                        :Trapezoid=>2,
                                                        :Rosenbrock32=>3,
                                                        :Feagin10=>10,
                                                        :Feagin12=>12,
                                                        :Feagin14=>14,
                                                        :Feagin10Vectorized=>10,
                                                        :Feagin12Vectorized=>12,
                                                        :Feagin14Vectorized=>14)


const DIFFERENTIALEQUATIONSJL_ADAPTIVEORDERS = Dict{Symbol,Int}(:ExplicitRK=>4, #Gets overwritten
                                                                :ExplicitRKVectorized=>4,#Gets overwritten
                                                                :BS3=>2,
                                                                :BS3Vectorized=>3,
                                                                :DP5=>4,
                                                                :DP5Vectorized=>4,
                                                                :DP8=>6,
                                                                :DP8Vectorized=>6,
                                                                :Rosenbrock32=>2,
                                                                :Feagin10=>8,
                                                                :Feagin12=>10,
                                                                :Feagin14=>12,
                                                                :Feagin10Vectorized=>8,
                                                                :Feagin12Vectorized=>10,
                                                                :Feagin14Vectorized=>12)
const DIFFERENTIALEQUATIONSJL_ADAPTIVEALGS = Set([:ExplicitRK,:ExplicitRKVectorized,:BS3,:BS3Vectorized,:DP5,:DP5Vectorized,:DP8,:DP8Vectorized,:Rosenbrock32,:Feagin10,:Feagin12,:Feagin14,:Feagin10Vectorized,:Feagin12Vectorized,:Feagin14Vectorized])
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
