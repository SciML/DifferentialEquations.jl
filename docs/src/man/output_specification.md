# Output Specification

DifferentialEquations.jl allows for specifying many forms of output. The default
is "as verbose as possible", with items saved to give a continuous interpolating
function as the output for ease of use. However, all of this functionality can
be toned down or turned off in order to improve performance and reduce the memory
usage. This page is to describe the different techniques which can be employed
to change the output specification. It will be described from the top down: the
most powerful is continuous (dense) output, which can instead be used for step-wise
interpolation via `saveat`, to using no interpolations and only save the timeseries
at `timeseries_steps`, to finally turning `save_timeseries=false` to only save the
value at the end.

## Availability

Note that the `dense` and `saveat` functions (the functionality which involves interpolations)
is currently only available for the ODE solvers. The other functionality is available
with all solvers.

## Continuous (Dense) Output

Continuous output is designated by the keyword argument `dense`. This is only available
for the ODE solvers.  By default this is turned on with `dense=true`. At every timepoint
it saves the appropriate derivative approximations `sol.k` in order to produce an
interpolating function which is accessed directly by calling the solution object.
For example, if `sol` is the solution object, the value at time `t` can be found via

```julia
sol(t)
```

Many of the special Runge-Kutta methods include a high-order interpolant which
matches or is one less than the order of the solver. By default the other methods
use an Order 3 Hermite interpolation. Since the `k` array must be stored, this
has the highest memory requirement. Note that for methods which use extra steps
for the high order interpolant that the extra steps are lazy evaluated and thus
only computing when an interpolated value in the appropriate interval is requested

## Choosing Intermediate Locations for the Solution

If `dense` solving is too high of a memory cost, one can specify values to be
interpolated during the solving via the array `saveat`. For example, if we are
solving on the interval `tspan=[0,1]`, we can set `saveat=[0.5]` and the solver
will ensure that an approximate value will be given at `t=0.5`. If this value is
reached by the solver, it will be ignored. If the solver skips over this point,
then an interpolated value is computed and saved for this point. This only requires
saving the derivatives at two timesteps, and thus has a drastically reduced memory
footprint than full dense output. Note that this, being associated with dense output,
is only available for the ODE solvers.

Another way to specify an output location is to add that value to `tspan`. For example,
we can force the solver to solve at `0.5` via `tspan=[0,0.5,1]`. However, notice that
this will require that the solver actually hits `t=0.5`. In some cases this can slow
down the solver by lowering the `Δt` leading to extra steps. In some cases, this may
be advantagous. For example, if you know that there is a large discontinuity at
`t=0.5`, using `tspan=[0,0.5,1]` will force the solver to first solve on `[0,0.5]`,
and then continue to solve on `[0.5,1]`. This will give a much better approximation
by perfectly specifying the moment of discontinuity, and can help the solver through
tough regions.

## Timeseries Specifications

To further reduce the memory usage, we can control the density that the timeseries
is saved at. By default `timeseries_steps=1`, meaning that every step is saved.
Note that `timeseries_steps=1` is required for dense output to work correctly.
If we change this value to `timeseries_steps=n`, then every `n`th step will be
saved. Note that it will always have the first and final steps. We can turn off
the saving of intermediate steps completely via the keyword `save_timeseries=false`.
This can be used to minimize the memory usage.
