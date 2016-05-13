
<a id='Convergence-Simulations-1'></a>

# Convergence Simulations


The convergence simulation


<a id='The-ConvergenceSimulation-Type-1'></a>

## The ConvergenceSimulation Type

<a id='DifferentialEquations.ConvergenceSimulation' href='#DifferentialEquations.ConvergenceSimulation'>#</a>
**`DifferentialEquations.ConvergenceSimulation`** &mdash; *Type*.



ConvergenceSimulation

A type which holds the data from a convergence simulation.


<a id='Related-Functions-1'></a>

## Related Functions

<a id='Base.length-Tuple{DifferentialEquations.ConvergenceSimulation}' href='#Base.length-Tuple{DifferentialEquations.ConvergenceSimulation}'>#</a>
**`Base.length`** &mdash; *Method*.



length(simres::ConvergenceSimulation)

Returns the number of simultations in the Convergence Simulation

<a id='DifferentialEquations.conv_ests' href='#DifferentialEquations.conv_ests'>#</a>
**`DifferentialEquations.conv_ests`** &mdash; *Function*.



conv_ests(error::Vector{Number})

Computes the pairwise convergence estimate for a convergence test done by halving/doubling stepsizes via

log2(error[i+1]/error[i])

Returns the mean of the convergence estimates


<a id='Plot-Functions-1'></a>

## Plot Functions

<a id='DifferentialEquations.convplot_fullΔt' href='#DifferentialEquations.convplot_fullΔt'>#</a>
**`DifferentialEquations.convplot_fullΔt`** &mdash; *Function*.



convplot_fullΔt(simres::ConvergenceSimulation;titleStr="All Convergences",savefile="")

<a id='DifferentialEquations.convplot_fullΔx' href='#DifferentialEquations.convplot_fullΔx'>#</a>
**`DifferentialEquations.convplot_fullΔx`** &mdash; *Function*.



convplot_fullΔx(simres::ConvergenceSimulation;titleStr="All Convergences",savefile="")

<a id='DifferentialEquations.convplot_node2vsΔt' href='#DifferentialEquations.convplot_node2vsΔt'>#</a>
**`DifferentialEquations.convplot_node2vsΔt`** &mdash; *Function*.



convplot_node2vsΔt(simres::ConvergenceSimulation)

<a id='DifferentialEquations.convplot_maxvsΔx' href='#DifferentialEquations.convplot_maxvsΔx'>#</a>
**`DifferentialEquations.convplot_maxvsΔx`** &mdash; *Function*.



convplot_maxvsΔx(simres::ConvergenceSimulation)

<a id='DifferentialEquations.convplot_l2vsΔt' href='#DifferentialEquations.convplot_l2vsΔt'>#</a>
**`DifferentialEquations.convplot_l2vsΔt`** &mdash; *Function*.



convplot_l2vsΔt(simres::ConvergenceSimulation)

<a id='DifferentialEquations.convplot_h1vsΔt' href='#DifferentialEquations.convplot_h1vsΔt'>#</a>
**`DifferentialEquations.convplot_h1vsΔt`** &mdash; *Function*.



convplot_h1vsΔt(simres::ConvergenceSimulation)

<a id='DifferentialEquations.convplot_l2vsΔx' href='#DifferentialEquations.convplot_l2vsΔx'>#</a>
**`DifferentialEquations.convplot_l2vsΔx`** &mdash; *Function*.



convplot_l2vsΔx(simres::ConvergenceSimulation)

<a id='DifferentialEquations.convplot_h1vsΔx' href='#DifferentialEquations.convplot_h1vsΔx'>#</a>
**`DifferentialEquations.convplot_h1vsΔx`** &mdash; *Function*.



convplot_h1vsΔx(simres::ConvergenceSimulation)

<a id='DifferentialEquations.convplot_node2vsΔx' href='#DifferentialEquations.convplot_node2vsΔx'>#</a>
**`DifferentialEquations.convplot_node2vsΔx`** &mdash; *Function*.



convplot_node2vsΔx(simres::ConvergenceSimulation)

<a id='DifferentialEquations.convplot_maxvsΔt' href='#DifferentialEquations.convplot_maxvsΔt'>#</a>
**`DifferentialEquations.convplot_maxvsΔt`** &mdash; *Function*.



convplot_maxvsΔt(simres::ConvergenceSimulation)

