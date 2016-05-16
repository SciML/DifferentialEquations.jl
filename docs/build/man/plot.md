
<a id='Plot-Functions-1'></a>

# Plot Functions


<a id='Related-Functions-1'></a>

## Related Functions

<a id='DifferentialEquations.solplot_appxvstrue' href='#DifferentialEquations.solplot_appxvstrue'>#</a>
**`DifferentialEquations.solplot_appxvstrue`** &mdash; *Function*.



solplot_appxvstrue(res::FEMSolution)

Plots the approximate solution and the true solution.

**Keyword Arguments**

  * `savefile`: Designates a file to save the plot in. Save type is defined by the chosen extension. Default is "" which implies no saving.
  * `title`: The title at the top of the plot. Default is "PDE Solution".
  * `appxTitle`: The title above the approximate solution. Default is "Approximated Solution".
  * `trueTitle`: The title above the true solution. Default is "True Solution".
  * `cmap`: Specifies the color map in the plot. Default is PyPlot.get_cmap("winter").

<a id='DifferentialEquations.solplot_animation' href='#DifferentialEquations.solplot_animation'>#</a>
**`DifferentialEquations.solplot_animation`** &mdash; *Function*.



solplot_animation(res::FEMSolution)

Plots an animation of the solution. Requires `fullSave=true` was enabled in the solver.

**Keyword Arguments**

  * `zlim`: The limits on the z-axis in the simulation. Default nothing.
  * `cbar`: Boolean flag which turns on/off the color bar. Default true.

<a id='DifferentialEquations.convplot' href='#DifferentialEquations.convplot'>#</a>
**`DifferentialEquations.convplot`** &mdash; *Function*.



convplot(measure,err)

Makes a convergence plot of err vs measure in loglog scale.

**Keyword Arguments**

  * `ErrStr`: The y-axis label. Default is "Error".
  * `measureStr`: The x-axis label. Default is "Measure".
  * `titleStr`: The title. Default is ":ErrStr vs :measureStr Convergence Plot".

<a id='DifferentialEquations.solplot' href='#DifferentialEquations.solplot'>#</a>
**`DifferentialEquations.solplot`** &mdash; *Function*.



solplot(res::FEMSolution)

Plots the approximate solution and, if available, the true solution.

**Keyword Arguments**

  * `savefile`: Designates a file to save the plot in. Save type is defined by the chosen extension. Default is "" which implies no saving.
  * `title`: The title at the top of the plot. Default is "PDE Solution".
  * `appxTitle`: The title above the approximate solution. Default is "Approximated Solution".
  * `trueTitle`: The title above the true solution. Default is "True Solution".
  * `cmap`: Specifies the color map in the plot. Default is PyPlot.get_cmap("winter").

<a id='DifferentialEquations.solplot_appx' href='#DifferentialEquations.solplot_appx'>#</a>
**`DifferentialEquations.solplot_appx`** &mdash; *Function*.



solplot_appx(res::FEMSolution)

Plots the approximate solution.

**Keyword Arguments**

  * `savefile`: Designates a file to save the plot in. Save type is defined by the chosen extension. Default is "" which implies no saving.
  * `title`: The title at the top of the plot. Default is "PDE Solution".
  * `cmap`: Specifies the color map in the plot. Default is PyPlot.get_cmap("winter").

<a id='DifferentialEquations.showmesh' href='#DifferentialEquations.showmesh'>#</a>
**`DifferentialEquations.showmesh`** &mdash; *Function*.



showmesh(femMesh::FEMmesh)

Shows the mesh which is defined by the (node,elem) structure.

**Keyword Arguments**

  * `cmap`: Specifies the color map in the plot. Default is PyPlot.get_cmap("winter").

