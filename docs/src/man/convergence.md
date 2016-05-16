# Convergence Simulations

The convergence simulation type is useful for deriving order of convergence estimates
from a group of simulations. This object will automatically assemble error vectors
into a more useful manner and provide plotting functionality. Convergence estimates
are also given by pair-wise estimates.

## The ConvergenceSimulation Type

```
{docs}
DifferentialEquations.ConvergenceSimulation
```

## Related Functions

```
{docs}
Base.length(::ConvergenceSimulation)
DifferentialEquations.conv_ests
```

## Plot Functions

```
{docs}
DifferentialEquations.convplot_fullΔt
DifferentialEquations.convplot_fullΔx
DifferentialEquations.convplot_node2vsΔt
DifferentialEquations.convplot_maxvsΔx
DifferentialEquations.convplot_l2vsΔt
DifferentialEquations.convplot_h1vsΔt
DifferentialEquations.convplot_l2vsΔx
DifferentialEquations.convplot_h1vsΔx
DifferentialEquations.convplot_node2vsΔx
DifferentialEquations.convplot_maxvsΔt
```
