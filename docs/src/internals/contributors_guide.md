# Contributor's Guide

So you're looking to help out DifferentialEquations.jl? We'd be happy to have
your help. It is recommended you first discuss with some of the developers
[on the Gitter channel](https://gitter.im/ChrisRackauckas/DifferentialEquations.jl)
to make sure that you're up-to-date with current developments.

## Developing New Solver Algorithms

The easiest way to get started would be to add new solver algorithms. This is a
pretty simple task as there are tools which put you right into the "hot loop".
For example, take a look at the ODE solver code. The mode `solve(::ODEProblem,::AbstractArray)`
is glue code to a bunch of solver algorithms. The algorithms which are coded
in DifferentialEquations.jl can be found in ode_integrators.jl. For example,
take a look at the Midpoint method's implementation:

```julia
function ode_midpoint(f::Function,u::AbstractArray,t,Δt,T,iter,
                      maxiters,timeseries,ts,timeseries_steps,save_timeseries,adaptive,progressbar)
  halfΔt = Δt/2
  utilde = similar(u)
  while t < T
    @ode_loopheader
    utilde[:] = u+halfΔt.*f(u,t)
    u = u + Δt.*f(utilde,t+halfΔt)
    @ode_loopfooter
  end
  return u,t,timeseries,ts
end
```

The parts in the signature are the items you have available. Most are self-explanatory
(they are from the ODE problem). The extra are for parts of the header and footer
for exiting at max iterations, and plugging into the Juno progressbar. These are
done in the `@ode_loopheader` and `@ode_loopfooter` macros, which are defined
using the `@def` macro (they essentially copy-paste the code from the line which
says `@def ode_loopheader begin ... end`). Note that the loopfooter code takes
care of the code for doing the adaptive timestepping. All that is required for
the adaptivity is that the algorithm computes an error estimate `EEst` each time,
save the value `utmp` to be what will replace `u` if the step is not rejected,
and add the algorithm's symbol is added to the dictionary `ODE_DIFFERENTIALEQUATIONSJL_ADAPTIVEALGS`
in `ode_constants.jl`. If implicit solving is needed (via NLsolve),
add the algorithm's symbol to `DIFFERENTIALEQUATIONSJL_IMPLICITALGS` and the
conditional dependency will be supplied. Note that you may need more function
arguments. Use another method as a template.

When the solver is completed, add a call to the solver in the glue code
`solve(::ODEProblem,::AbstractArray)` (you will see all the others),
add the symbol for the algorithm to `DIFFERENTIALEQUATIONSJL_ALGORITHMS`, and
the order to `DIFFERENTIALEQUATIONSJL_ORDERS`. It's that quick! Lastly, add
your method to the convergence tests in the appropriate /test file.  Feel free
to implement any interesting or educational algorithm: they don't have to be
the fastest and it is always is useful to have such algorithms (like Simpson's method)
available for demonstration purposes.

Adding algorithms to the other problems is very similar.

## Adding Conditional Dependencies

If your algorithm requires a conditional dependency (a package, but not one
that everyone who uses DifferentialEquations.jl would need), you can add them
as follows. Before the loop, add the line `initialize_backend(:PkgName)` where
`:PkgName` is the same name as the package you wish to use. Then, in `general/backends.jl`
add a dispatch to `init_package`. A common one would be:

```julia
init_package(b::backend{:PkgName}) = @eval begin
      import PkgName
      export PkgName
    end
```

Now inside your method you can use any function from the package via PkgName.function.
The first time it is used it import the package (or tell the user to install it).

## Developing A New Problem

To develop a new problem, you need to make a new `DEProblem` and a new `DESolution`.
The `DEProblem` type should hold all of the mathematical information about the
problem, and the `DESolution` should hold all of the information for the solution.
Then all that is required is to define a `solve(::DEProblem,*Extra Mesh Things*;kwargs)`
which takes in the problem and returns a solution. To add plotting functionality,
add a plot recipe for the solution type to `/general/plotrecipes`. For testing
that the algorithm works, add a dispatch for `test_convergence` which makes
a `ConvergenceSimulation` type. This type already has a plot recipe, so
plotting functionality will already be embedded. This requires that your
problem can take in a true solution, and has a field `errors` which is a
dictionary of symbols for the different error estimates (L2,L infinity, etc.)

## Other Help

There's always more to be. Improved plot recipes and new series recipes are
always nice to add more default plots. It is always helpful to have benchmarks
between different algorithms to see "which is best". Adding examples IJulia
notebooks to `/examples/` is a good way to share knowledge about DifferentialEquations.jl.
Also, please feel free to comb through the solvers and look for ways to make them
more efficient. Lastly, the documentation could always use improvements. If you
have any questions on how to help, just ask them in the Gitter!
