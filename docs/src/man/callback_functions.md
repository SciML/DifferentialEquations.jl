# Event Handling and Callback Functions

## Introduction to Callback Function

DifferentialEquations.jl allows for using callback functions to inject user code
into the solver algorithms. This is done by defining a callback function and
passing that function to the solver. After each accepted iteration this function
is called. The standard callback is defined as:

```julia
default_callback = @ode_callback begin
  @ode_savevalues
end
```

This runs the saving routine at every timestep (inside of this saving routine it
  checks the iterations vs `timeseries_steps` etc., so it's quite complicated).
  However, you can add any code you want to this callback. For example, we can
  make it print the value at each timestep by doing:

```julia
my_callback = @ode_callback begin
  println(u)
  @ode_savevalues
end
```

and pass this to the solver:

```julia
sol = solve(prob,tspan,callback=my_callback)
```

Later in the manual the full API for callbacks is given (the callbacks are very
  general and thus the full API is complex enough to handle just about anything),
  but for most users it's recommended that you use the simplified event handling
  DSL described below.

## Event Handling

Since event handling is a very common issue, a special domain-specific language
(DSL) was created to make event handling callbacks simple to define.

### Example 1: Bouncing Ball

First let's look at the bouncing ball. `@ode_def` from
[ParameterizedFunctions.jl](https://github.com/JuliaDiffEq/ParameterizedFunctions.jl)
was to define the problem, where the first variable `y` is the height which changes
by `v` the velocity, where the velocity is always changing at `-g` where is the
gravitational constant. This is the equation:

```julia
f = @ode_def BallBounce begin
  dy =  v
  dv = -g
end g=9.81
```

All we have to do in order specify the event is to have a function which
should always be positive with an event occurring at 0. For now at least
that's how it's specified, if a generalization is needed we can talk about
this (but it needs to be "root-findable"). For here it's clear that we just
want to check if the ball's height ever hits zero:

```julia
function event_f(t,u) # Event when event_f(t,u,k) == 0
  u[1]
end
```

Now we have to say what to do when the event occurs. In this case we just
flip the velocity (the second variable)

```julia
function apply_event!(u,cache)
  u[2] = -u[2]
end
```

That's all you need to specify the callback function with the macro:

```julia
callback = @ode_callback begin
  @ode_event event_f apply_event!
end
```

One thing to note is that by default this will only check at each timestep if
the event condition is satisfied (i.e. if `event_f(t,u)<0`). If your problem
is oscillatory, sometime too large of a timestep will miss the event. In that
case, you will want to specify a number of points in the interval to interpolate
at and check the condition as well. This is done with one more parameter to `@ode_event`.

Lastly, you can also tell the solver to decrease Δt after the event occurs.
This can be helpful if the discontinuity changes the problem immensely.
Using the full power of the macro, we can define an event as

```julia
const Δt_safety = 1 # Multiplier to Δt after an event
const interp_points = 10
callback = @ode_callback begin
  @ode_event event_f apply_event! interp_points Δt_safety
end
```

Then you can solve and plot:

```julia
u0 = [50.0,0.0]
prob = ODEProblem(f,u0)
tspan = [0;15]
sol = solve(prob,tspan,callback=callback)
plot(sol)
```

[BallBounce](../assets/ballbounce.png)

As you can see from the resulting image, DifferentialEquations.jl is smart enough
to use the interpolation to hone in on the time of the event and apply the event
back at the correct time. Thus one does not have to worry about the adaptive timestepping
"overshooting" the event as this is handled for you.

### Example 2: Growing Cell Population

Another interesting issue are models of changing sizes. The ability to handle
such events is a unique feature of DifferentialEquations.jl! The problem we would
like to tackle here is a cell population. We start with 1 cell with a protein `X`
which increases linearly with time with rate parameter `α`. Since we are going
to be changing the size of the population, we write the model in the general form:

```julia
const α = 0.3
f = function (t,u,du)
  for i in 1:length(u)
    du[i] = α*u[i]
  end
end
```

Our model is that, whenever the protein `X` gets to a concentration of 1, it
triggers a cell division. So we check to see if any concentrations hit 1:

```julia
function event_f(t,u) # Event when event_f(t,u,k) == 0
  1-maximum(u)
end
```

Again, recall that this function finds events as switching from positive to negative,
so `1-maximum(u)` is positive until a cell has a concentration of `X` which is
1, which then triggers the event. At the event, we have that the call splits
into two cells, giving a random amount of protein to each one. We can do this
by resizing the cache (adding 1 to the length of all of the caches) and setting
the values of these two cells at the time of the event:

```julia
function apply_event!(u,cache)
  @ode_change_cachesize cache length+1
  maxidx = findmax(u)[2]
  Θ = rand()
  u[maxidx] = Θ
  u[end] = 1-Θ
end
```

`@ode_change_cachesize cache length+1` is used to change the length of all of the
internal caches (which includes `u`) to be their current length + 1, growing the
ODE system. Then the following code sets the new protein concentrations. Now we
can solve:

```julia
const Δt_safety = 1
const interp_points = 10
callback = @ode_callback begin
  @ode_event event_f apply_event! interp_points Δt_safety
end
u0 = [0.2]
prob = ODEProblem(f,u0)
tspan = [0;10]
sol = solve(prob,tspan,callback=callback)
```

The plot recipes do not have a way of handling the changing size, but we can
plot from the solution object directly. For example, let's make a plot of how
many cells there are at each time. Since these are discrete values, we calculate
and plot them directly:

```julia
plot(sol.t,map((x)->length(x),sol[:]),lw=3,
     ylabel="Number of Cells",xlabel="Time")
```

[NumberOfCells](../assets/numcells.png)

Now let's check-in on a cell. We can still use the interpolation to get a nice
plot of the concentration of cell 1 over time. This is done with the command:

```julia
ts = linspace(0,10,100)
plot(ts,map((x)->x[1],sol.(ts)),lw=3,
     ylabel="Amount of X in Cell 1",xlabel="Time")
```

[Cell1](../assets/cell1.png)

Notice that every time it hits 1 the cell divides, giving cell 1 a random amount
of `X` which then grows until the next division.

## Advanced: Callback Function API

The callback functions have access to a lot of the functionality of the solver.
The macro defines a function which is written as follows:

```julia
macro ode_callback(ex)
  esc(quote
    function (alg,f,t,u,k,tprev,uprev,kprev,ts,timeseries,ks,Δtprev,Δt,saveat,cursaveat,iter,save_timeseries,timeseries_steps,uEltype,ksEltype,dense,kshortsize,issimple_dense,fsal,fsalfirst,cache)
      reeval_fsal = false
      event_occured = false
      $(ex)
      cursaveat,Δt,t,reeval_fsal
    end
  end)
end
```

All of the parts of the algorithm are defined in the internal solver documentation.

### Example: Bouncing Ball Without Macros

Here is an example of the defining the ball bouncing callback without the usage
of macros. The entire code in its fully glory is generic enough to handle any
of the implemented DifferentialEquations.jl algorithms, which special differences
depending on the type of interpolant, implementation of FSAL, etc. For these
reasons it's usually recommended to use the event handling macro, though this kind
of code will allow you handle pretty much anything!

```julia
manual_callback = function (alg,f,t,u,k,tprev,uprev,kprev,ts,timeseries,ks,Δtprev,Δt,saveat,cursaveat,iter,save_timeseries,timeseries_steps,uEltype,ksEltype,dense,kshortsize,issimple_dense,fsal,fsalfirst,cache)
  reeval_fsal = false
  event_occured = false
  Δt_safety = 1
  interp_points = 10

  # Event Handling
  ode_addsteps!(k,tprev,uprev,Δtprev,alg,f)
  Θs = linspace(0,1,interp_points)
  interp_index = 0
  # Check if the event occured
  if event_f(t,u)<0
    event_occured = true
    interp_index = interp_points
  elseif interp_points!=0 # Use the interpolants for safety checking
    for i in 2:length(Θs)-1
      if event_f(t+Δt*Θs[i],ode_interpolant(Θs[i],Δtprev,uprev,u,kprev,k,alg))<0
        event_occured = true
        interp_index = i
        break
      end
    end
  end

  if event_occured
    if interp_index == interp_points # If no safety interpolations, start in the middle as well
      initial_Θ = [.5]
    else
      initial_Θ = [Θs[interp_index]] # Start at the closest
    end
    find_zero = (Θ,val) -> begin
      val[1] = event_f(t+Θ[1]*Δt,ode_interpolant(Θ[1],Δtprev,uprev,u,kprev,k,alg))
    end
    res = nlsolve(find_zero,initial_Θ)
    val = ode_interpolant(res.zero[1],Δtprev,uprev,u,kprev,k,alg)
    for i in eachindex(u)
      u[i] = val[i]
    end
    Δtprev *= res.zero[1]
    t = tprev + Δtprev

    if alg ∈ DIFFERENTIALEQUATIONSJL_SPECIALDENSEALGS
      resize!(k,kshortsize) # Reset k for next step
      k = typeof(k)() # Make a local blank k for saving
      ode_addsteps!(k,tprev,uprev,Δtprev,alg,f)
    elseif typeof(u) <: Number
      k = f(t,u)
    else
      f(t,u,k)
    end
  end

  @ode_savevalues

  if event_occured
    apply_event!(u)
    if alg ∉ DIFFERENTIALEQUATIONSJL_SPECIALDENSEALGS
      if typeof(u) <: Number
        k = f(t,u)
      else
        f(t,u,k)
      end
    end
    @ode_savevalues
    if fsal
      reeval_fsal = true
    end
    Δt *= Δt_safety # Safety Δt change
  end

  cursaveat,Δt,t,reeval_fsal
end
```
