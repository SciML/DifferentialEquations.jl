# Juno Progress Bar Integration

DifferentialEquations.jl integrates with the Juno progress bar in order to make
long calculations more manageable. By default this feature is off for ODE and
SDE solvers, but can be turned on via the keyword argument `progressbar=true`.
The progress bar updates every `progress_steps` timesteps, which has a default
value of 1000. Note that making this value really low could cause a performance
hit, though from some basic testing it seems that with updates of at least
 1000 steps on number (the fastest problems) there's no discernable performance degradation,
 giving an high upper bound.
