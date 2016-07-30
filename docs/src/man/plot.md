# Plot Functions

## Standard Plots

Plotting functionality is provided by a recipe to Plots.jl. To
use plot solutions, simply call the `plot(type)` and the plotter will generate
appropriate plots. If `save_timeseries` was used, the plotters can
generate animations of the solutions to evolution equations.
Plots can be customized using [all of the keyword arguments
provided by Plots.jl](https://juliaplots.github.io/supported/). Please see [Plots.jl's documentation](https://juliaplots.github.io/) for more information.

A few extra arguments are provided addition to the Plots.jl keyword arguments. They are as follows:

* plot_analytic: Specifies whether the true solution (if known) should be plotted alongside
  the numerically approximated solution. Default is false.
* filename: Specifies the filename to save an animation. Only applies to the `animate` function.
* fps: Determines the fps in an animation. Only applies to the `animate` function.

## Extra Plot Functions

```@docs
DifferentialEquations.animate
```
