# Benchmark Suite

DiffernetialEquations.jl provides a benchmarking suite to be able to test the
difference in error, speed, and efficiency between algorithms. DifferentialEquations.jl includes current benchmarking notebooks to help users
understand the performance of the methods. These benchmarking notebooks use the included benchmarking suite. There are two parts to the benchmarking suite: shootouts and work-precision. The `Shootout` tests methods head-to-head for timing and error on the same problem. A `WorkPrecision` draws a work-precision diagram
for the algorithms in question on the chosen problem.

### Using the Benchmarking Notebooks

To use the benchmarking notebooks, IJulia is required. The commands are as follows:

```julia
using IJulia
notebook(dir = Pkg.dir("DifferentialEquations")*"/benchmarks")
```

### Shootout

 A
shootout is where you compare between algorithms. For example, so see how
different Runge-Kutta algorithms fair against each other, one can define a setup
which is a dictionary of Symbols to Any, where the symbol is the keyword argument.
Then you call `ode_shootout` on that setup. The code is as follows:

```julia
tspan = [0,10]
setups = [Dict(:alg=>:DP5)
          Dict(:abstol=>1e-3,:reltol=>1e-6,:alg=>:ode45) # Fix ODE to be normal
          Dict(:alg=>:dopri5)]
prob = DifferentialEquations.prob_ode_large2Dlinear
names = ["DifferentialEquations";"ODE";"ODEInterface"]
shoot = ode_shootout(prob,tspan,setups;Δt=1/2^(10),names=names)
```

Note that keyword arguments applied to ode_shootout are applie dot every run, so
in this example every run has the same starting timestep.  Here we explicitly chose names.
If you don't, then the algorithm name is the default.
This returns a Shootout type where which holds the times it took for each algorithm
and the errors. Using these, it calculates the efficiency defnied as
1/(error*time), i.e. if the error is low or the run was quick then
it's efficient. `print(shoot)` will show all of this information,
and `plot(shoot)` will show the efficiencies of the algorithms
in comparison to each other.

For every benchmark function there is a special keyword `numruns` which controls
the number of runs used in the time estimate. To be more precise, these functions
by default run the algorithm 20 times on the problem and take the average time.
This amount can be increased and decreased as needed.

A ShootoutSet is a where you define a vector of probs and tspans and run a shootout
on each of these values.

### WorkPrecision

A WorkPrecision calculates the necessary componnets of a work-precision plot. This
shows how time scales with the user chosen tolerances on a given problem. To make
a WorkPrecision, you give it a vector of absolute and relative tolerances:

```julia
abstols = 1./10.^(3:10)
reltols = 1./10.^(3:10)
wp = ode_workprecision(prob,tspan,abstols,reltols;alg=:DP5,name="Dormand-Prince 4/5")
```

If we want to plot many WorkPrecisions together in order to compare between
algorithms, you can make a WorkPrecisionSet. To do so, you pass the setups
into the function as well:

```julia
wp_set = ode_workprecision_set(prob,tspan,abstols,reltols,setups;Δt=1/2^4,numruns=2)
setups = [Dict(:alg=>:RK4);Dict(:alg=>:Euler);Dict(:alg=>:BS3);
          Dict(:alg=>:Midpoint);Dict(:alg=>:BS5);Dict(:alg=>:DP5)]
wp_set = ode_workprecision_set(prob,tspan,abstols,reltols,setups;Δt=1/2^4,numruns=2)
```

Both of these types have a plot recipe to produce a work-precision diagram,
and a print which will show some relevant information.
