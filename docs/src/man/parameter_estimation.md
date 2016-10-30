# Parameter Estimation

Parameter estimation for ODE models is provided by the DiffEq suite. The current
functionality includes `build_optim_objective` and `lm_fit`. Note these require
that the problem is defined using a [ParameterizedFunction](https://github.com/JuliaDiffEq/ParameterizedFunctions.jl).

### build_optim_objective

`build_optim_objective` builds an objective function to be used with Optim.jl.

```julia
build_optim_objective(prob::DEProblem,tspan,t,data;loss_func = L2DistLoss,kwargs...)
```

The first argument is the DEProblem to solve. Second is the `tspan`. Next is `t`,
the set of timepoints which the data is found at. The last argument which is required
is the data, which are the values where are known, in order to be optimized against.
Optionally, one can choose a loss function from LossFunctions.jl or use the default
of an L2 loss. The keyword arguments are passed to the ODE solver.

### lm_fit

`lm_fit` is a function for fitting the parameters of an ODE using the Levenberg-Marquardt
algorithm. This algorithm is really bad and thus not recommended since, for example,
the Optim.jl algorithms on an L2 loss are more performant and robust. However,
this is provided for completeness as most other differential equation libraries
use an LM-based algorithm, so this allows one to test the increased effectiveness
of not using LM.

```julia
lm_fit(prob::DEProblem,tspan,t,data,p0;kwargs...)
```

The arguments are similar to before, but with `p0` being the initial conditions
for the parameters and the `kwargs` as the args passed to the LsqFit `curve_fit`
function (which is used for the LM solver). This returns the fitted parameters.

#### Examples

We choose to optimize the parameters on the Lotka-Volterra equation. We do so
by defining the function as a [ParmaeterizedFunction](https://github.com/JuliaDiffEq/ParameterizedFunctions.jl):


```julia
f = @ode_def_nohes LotkaVolterraTest begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=1.0 c=3.0 d=1.0

u0 = [1.0;1.0]
tspan = [0;10.0]
prob = ODEProblem(f,u0)
```

Notice that since we only used `=>` for `a`, it's the only free parameter.
We create data using the numerical result with `a=1.5`:

```julia
t = collect(linspace(0,10,200))
randomized = [(sol(t[i]) + .01randn(2)) for i in 1:length(t)]
data = vecvec_to_mat(randomized)
```

Here we used `vecvec_to_mat` from [RecursiveArrayTools.jl](https://github.com/ChrisRackauckas/RecursiveArrayTools.jl)
to turn the result of an ODE to a matrix.

To build the objective function for Optim.jl, we simply call the `build_optim_objective`
funtion:

```julia
cost_function = build_optim_objective(prob,tspan,t,data,alg=:Vern6)
```

Now this cost function can be used with Optim.jl in order to get the parameters.
For example, we can use Brent's algorithm to search for the best solution on
the interval `[0,10]` by:

```julia
result = optimize(cost_function, 0.0, 10.0)
```

This returns `result.minimum[1]==1.5` as the best parameter to match the data.
We can also use the multivariate optimization functions. For example, we can use
the `BFGS` algorithm to optimize the parameter starting at `a=1.42` using

```julia
result = optimize(cost_function, [1.42], BFGS())
```

Note that some of the algorithms may be sensitive to the initial condtion. For more
details on using Optim.jl, see the [documentation for Optim.jl](http://www.juliaopt.org/Optim.jl/latest/).
