# Ordinary Differential Equation Solvers

```@docs
solve(::AbstractODEProblem,::AbstractArray)
```

## Recommended Methods

Currently, over 100 algorithm choices are available. This guide is to help
you choose the right one.

### Non-Stiff Problems

For non-stiff problems, the native DifferentialEquations.jl algorithms are vastly
more efficient than the other choices (ODEInterface and ODE.jl). For most non-stiff
problems, we recommend `:DP5` (this is the default algorithm). When more robust
error control is required, `:BS5` is a good choice.
For fast solving at lower tolerances, we recommend `:BS3`. For tolerances
which are at about the truncation error of Float64 (1e-16), we recommend
`:DP8` as a robust choice and `:Vern6`, `:Vern7`, or `:Vern8` as an efficient choice.

For high accuracy non-stiff solving, we recommend the `:Feagin12` or `:Feagin14`
methods. These are more robust than Adams-Bashforth methods to discontinuities
and achieve very high precision, and are much more efficient than the extrapolation
methods. Note that the Feagin methods are the only high-order optimized methods
which do not include a high-order interpolant (they do include a 3rd order
Hermite interpolation if needed). If a high-order method is needed with a high
order interpolant, then you should choose `:Vern9` which is Order 9 with an
Order 9 interpolant.

### Stiff Problems

For mildly stiff problems which are not oscillatory, `:Rosenbrock32` is a good choice.
If there exist large oscillations, then this form (local extrapolated) is not L-stable,
and thus it is recommended that you use `:Rosenbrock23` which is L-stable (but Order 2
instead of 3). As a native DifferentialEquations.jl solver, Julia-defined numbers will work.
This method uses ForwardDiff to automatically guess the Jacobian. For faster solving
when the Jacobian is known, use `radau`. For highly stiff problems where Julia-defined
numbers need to be used (SIUnits, Arbs), `:Trapezoid` is the current best choice.
However, for the most efficient highly stiff solvers, use `:radau` or `:cvode_BDF` provided by wrappers
to the ODEInterface and Sundials packages respectively ([see the conditional dependencies documentation](http://juliadiffeq.github.io/DifferentialEquations.jl/latest/man/conditional_dependencies.html))

## Full List of Methods

Choose one of these methods with the `alg` keyword in `solve`.

* DifferentialEquations.jl

Unless otherwise specified, the DifferentialEquations algorithms all come with a
3rd order Hermite polynomial interpolation. The algorithms denoted as having a "free"
interpolation means that no extra steps are required for the interpolation. For
the non-free higher order interpolating functions, the extra steps are computed
lazily (i.e. not during the solve).

  - `:Euler`- The canonical forward Euler method.
  - `:Midpoint` - The second order midpoint method.
  - `:RK4` - The canonical Runge-Kutta Order 4 method.
  - `:BS3` - Bogacki-Shampine 3/2 method.
  - `:DP5` - Dormand-Prince's 5/4 Runge-Kutta method. (free 4th order interpolant)
  - `Tsit5` - Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant)
  - `BS5` - Bogacki-Shampine 5/4 Runge-Kutta method. (5th order interpolant)
  - `Vern6` - Verner's "Most Efficient" 6/5 Runge-Kutta method. (6th order interpolant)
  - `Vern7` - Verner's "Most Efficient" 7/6 Runge-Kutta method. (7th order interpolant)
  - `TanYam7` - Tanaka-Yamashita 7 Runge-Kutta method.
  - `DP8` - Hairer's 8/5/3 adaption of the Dormand-Prince 8
    method Runge-Kutta method. (7th order interpolant)
  - `TsitPap8` - Tsitouras-Papakostas 8/7 Runge-Kutta method.
  - `Vern8` - Verner's "Most Efficient" 8/7 Runge-Kutta method. (8th order interpolant)
  - `Vern9` - Verner's "Most Efficient" 9/8 Runge-Kutta method. (9th order interpolant)
  - `:Feagin10` - Feagin's 10th-order Runge-Kutta method.
  - `:Feagin12` - Feagin's 12th-order Runge-Kutta method.
  - `:Feagin14` - Feagin's 14th-order Runge-Kutta method.
  - `:ExplicitRK` - A general Runge-Kutta solver which takes in a tableau. Can be adaptive. Tableaus
    are specified via the keyword argument `tab=tableau`. The default tableau is
    for Dormand-Prince 4/5. Other supplied tableaus can be found in the Supplied Tableaus section.
  - `:ImplicitEuler` - A 1st order implicit solver. Unconditionally stable.
  - `:Trapezoid` - A second order unconditionally stable implicit solver. Good for highly stiff.
  - `:Rosenbrock23` - An Order 2/3 L-Stable fast solver which is good for mildy stiff equations with oscillations.
  - `:Rosenbrock32` - An Order 3/2 A-Stable fast solver which is good for mildy stiff equations without oscillations.

* ODEInterface.jl

  - `:dopri5` - Hairer's classic implementation of the Dormand-Prince 4/5 method.
  - `:dop853` - Explicit Runge-Kutta 8(5,3) by Dormand-Prince
  - `:odex` - GBS extrapolation-algorithm based on the midpoint rule
  - `:seulex` - extrapolation-algorithm based on the linear implicit Euler method
  - `:radau` - implicit Runge-Kutta (Radau IIA) of variable order between 5 and 13
  - `:radau5` - implicit Runge-Kutta method (Radau IIA) of order 5

* ODE.jl

The ODE.jl algorithms all come with a 3rd order Hermite polynomial interpolation.

  - `:ode23` - Bogakai-Shampine's 2/3 method
  - `:ode45` - Dormand-Prince's 4/5 method
  - `:ode78` - Runge-Kutta-Fehlberg 7/8 method
  - `:ode23s` - Rosenbrock's 2/3 method
  - `:ode1` - Forward Euler
  - `:ode2_midpoint` - Midpoint Method
  - `:ode2_heun` - Heun's Method
  - `:ode4` - RK4
  - `:ode45_fe` - Runge-Kutta-Fehlberg 4/5 method

* Sundials.jl

  - `:cvode_BDF` - CVode Backward Differentiation Formula (BDF) solver.
  - `:cvode_Adams` - CVode Adams-Moulton solver

## List of Supplied Tableaus

A large variety of tableaus have been supplied by default. For the most useful
and common algorithms, a hand-optimized version is supplied and is recommended
for general uses (i.e. use `:DP5` instead of `:ExplicitRK` with `tableau=constructDormandPrince()`).
However, these serve as a good method for comparing between tableaus and understanding
the pros/cons of the methods. Implemented are every published tableau (that I know exist).
Note that user-defined tableaus also are accepted.
To see how to define a tableau, checkout the [premade tableau source code](https://github.com/JuliaDiffEq/DifferentialEquations.jl/blob/master/src/ode/ode_tableaus.jl).
Tableau docstrings should have appropriate citations (if not, file an issue).

A plot recipes is provided which will plot the stability region for a given tableau.

### Explicit Runge-Kutta Methods

* `constructEuler` - Euler's 1st order method.
* `constructHuen()` Huen's order 2 method.
* `constructRalston()` - Ralston's order 2 method.
* `constructKutta3` - Kutta's classic 3rd order method
* `constructRK4` - The classic 4th order "Runge-Kutta" method
* `constructRK438Rule` - The classic 4th order "3/8th's Rule" method
* `constructBogakiShampine3()` - Bogakai-Shampine's 2/3 method.
* `constructRKF4()` - Runge-Kutta-Fehlberg 3/4.
* `constructRKF5()` - Runge-Kutta-Fehlberg 4/5.
* `constructRungeFirst5()` - Runge's first 5th order method.
* `constructCassity5()` - Cassity's 5th order method.
* `constructLawson5()` - Lawson's 5th order method.
* `constructLutherKonen5` - Luther-Konen's first 5th order method.
* `constructLutherKonen52()` - Luther-Konen's second 5th order method.
* `constructLutherKonen53()` - Luther-Konen's third 5th order method.
* `constructPapakostasPapaGeorgiou5()` - Papakostas and PapaGeorgiou more stable order 5 method.
* `constructPapakostasPapaGeorgiou52()` - Papakostas and PapaGeorgiou more efficient order 5 method.
* `constructTsitouras5()` - Tsitouras's order 5 method.
* `constructBogakiShampine5()` - Bogaki and Shampine's Order 5 method.
* `constructSharpSmart5()` - Sharp and Smart's Order 5 method.
* `constructCashKarp()` - Cash-Karp method 4/5.
* `constructDormandPrince()` - Dormand-Prince 4/5.
* `constructButcher6()` - Butcher's first order 6 method.
* `constructButcher62()` - Butcher's second order 6 method.
* `constructButcher63()` - Butcher's third order 6 method.
* `constructDormandPrince6()` - Dormand-Prince's 5/6 method.
* `constructSharpVerner6()` Sharp-Verner's 5/6 method.
* `constructVerner916()` - Verner's more efficient order 6 method (1991).
* `constructVerner9162()` - Verner's second more efficient order 6 method (1991).
* `constructVernerRobust6()` - Verner's "most robust" order 6 method.
* `constructVernerEfficient6()` - Verner's "most efficient" order 6 method.
* `constructPapakostas6()` - Papakostas's order 6 method.
* `constructLawson6()` - Lawson's order 6 method.
* `constructTsitourasPapakostas6()` - Tsitouras and Papakostas's order 6 method.
* `constructDormandLockyerMcCorriganPrince6()` - the Dormand-Lockyer-McCorrigan-Prince order 6 method.
* `constructTanakaKasugaYamashitaYazaki6A()` - Tanaka-Kasuga-Yamashita-Yazaki order 6 method A.
* `constructTanakaKasugaYamashitaYazaki6B()` - Tanaka-Kasuga-Yamashita-Yazaki order 6 method B.
* `constructTanakaKasugaYamashitaYazaki6C()` - Tanaka-Kasuga-Yamashita-Yazaki order 6 method C.
* `constructTanakaKasugaYamashitaYazaki6D()` - Tanaka-Kasuga-Yamashita-Yazaki order 6 method D.
* `constructMikkawyEisa()` - Mikkawy and Eisa's order 6 method.
* `constructChummund6()` - Chummund's first order 6 method.
* `constructChummund62()` - Chummund's second order 6 method.
* `constructHuta6()` - Huta's first order 6 method.
* `constructHuta62()` - Huta's second order 6 method.
* `constructVerner6()` - An old order 6 method attributed to Verner.
* `constructDverk()` - The classic DVERK algorithm attributed to Verner.
* `constructClassicVerner6()` - A classic Verner order 6 algorithm (1978).
* `constructButcher7()` - Butcher's order 7 algorithm.
* `constructClassicVerner7()`- A classic Verner order 7 algorithm (1978).
* `constructVernerRobust7()` - Verner's "most robust" order 7 algorithm.
* `constructTanakaYamashitaStable7()` - Tanaka-Yamashita more stable order 7 algorithm.
* `constructTanakaYamashitaEfficient7()` - Tanaka-Yamashita more efficient order 7 algorithm.
* `constructSharpSmart7()` - Sharp-Smart's order 7 algorithm.
* `constructSharpVerner7()` - Sharp-Verner's order 7 algorithm.
* `constructVerner7()` - Verner's "most efficient" order 7 algorithm.
* `constructVernerEfficient7()` - Verner's "most efficient" order 7 algorithm.
* `constructClassicVerner8()` - A classic Verner order 8 algorithm (1978).
* `constructCooperVerner8()` - Cooper-Verner's first order 8 algorithm.
* `constructCooperVerner82()` - Cooper-Verner's second order 8 algorithm.
* `constructTsitourasPapakostas8()` - Tsitouras-Papakostas order 8 algorithm.
* `constructdverk78()` - The classic order 8 DVERK algorithm.
* `constructEnrightVerner8()` - Enright-Verner order 8 algorithm.
* `constructCurtis8()` - Curtis' order 8 algorithm.
* `constructVerner8()` - Verner's "most efficient" order 8 algorithm.
* `constructRKF8()` - Runge-Kutta-Fehlberg Order 7/8 method.
* `constructDormandPrice8()` - Dormand-Prince Order 7/8 method.
* `constructDormandPrince8_64bit()` - Dormand-Prince Order 7/8 method.
  Coefficients are rational approximations good for 64 bits.
* `constructVernerRobust9()` - Verner's "most robust" order 9 method.
* `constructVernerEfficient9()` - Verner's "most efficient" order 9 method.
* `constructSharp9()` - Sharp's order 9 method.
* `constructTsitouras9()` - Tsitouras's first order 9 method.
* `constructTsitouras92()` - Tsitouras's second order 9 method.
* `constructCurtis10()` - Curtis' order 10 method.
* `constructOno10()` - Ono's order 10 method.
* `constructFeagin10Tableau()` - Feagin's order 10 method.
* `constructCurtis10()` - Curtis' order 10 method.
* `constructBaker10()` - Baker's order 10 method.
* `constructHairer10()` Hairer's order 10 method.
* `constructFeagin12Tableau()` - Feagin's order 12 method.
* `constructOno12()` - Ono's order 12 method.
* `constructFeagin14Tableau()` Feagin's order 14 method.

### Implicit Runge-Kutta Methods

* `constructImplicitEuler` - The 1st order Implicit Euler method.
* `constructMidpointRule` - The 2nd order Midpoint method.
* `constructTrapezoidalRule` - The 2nd order Trapezoidal rule (2nd order LobattoIIIA)
* `constructLobattoIIIA4` - The 4th order LobattoIIIA
* `constructLobattoIIIB2` - The 2nd order LobattoIIIB
* `constructLobattoIIIB4` - The 4th order LobattoIIIB
* `constructLobattoIIIC2` - The 2nd order LobattoIIIC
* `constructLobattoIIIC4` - The 4th order LobattoIIIC
* `constructLobattoIIICStar2` - The 2nd order LobattoIIIC*
* `constructLobattoIIICStar4` - The 4th order LobattoIIIC*
* `constructLobattoIIID2` - The 2nd order LobattoIIID
* `constructLobattoIIID4` - The 4th order LobattoIIID
* `constructRadauIA3` - The 3rd order RadauIA
* `constructRadauIA5` - The 5th order RadauIA
* `constructRadauIIA3` - The 3rd order RadauIIA
* `constructRadauIIA5` - The 5th order RadauIIA

## Analysis of Methods

For an in-depth walkthrough of the various method pros/cons, see [notes on algorithms](http://juliadiffeq.github.io/DifferentialEquations.jl/latest/internals/notes_on_algorithms.html)
