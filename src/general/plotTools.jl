#plot(z, zlim = (should_override ? zlim_override : default(:zlim)))

"""
solplot_animation(res::FEMSolution)

Plots an animation of the solution. Requires `fullSave=true` was enabled in the solver.

### Keyword Arguments

* `zlim`: The limits on the z-axis in the simulation. Default nothing.
* `cbar`: Boolean flag which turns on/off the color bar. Default true.
"""
function solplot_animation(res::FEMSolution;zlim=nothing,cbar=true)
  Plots.pyplot(reuse=true,size=(750,750))
  if zlim==nothing
    @gif for i=1:size(res.uFull,2)
        surface(res.femMesh.node[:,1],res.femMesh.node[:,2],res.uFull[:,i],cbar=cbar)
    end
  else
    @gif for i=1:size(res.uFull,2)
        surface(res.femMesh.node[:,1],res.femMesh.node[:,2],res.uFull[:,i],zlim=zlim,cbar=cbar)
    end
  end
end

"""
solplot_appxvstrue(res::FEMSolution)

Plots the approximate solution and the true solution.

### Keyword Arguments

* `savefile`: Designates a file to save the plot in. Save type is defined by the chosen
extension. Default is "" which implies no saving.
* `title`: The title at the top of the plot. Default is "PDE Solution".
* `appxTitle`: The title above the approximate solution. Default is "Approximated Solution".
* `trueTitle`: The title above the true solution. Default is "True Solution".
* `cmap`: Specifies the color map in the plot. Default is PyPlot.get_cmap("winter").
"""
function solplot_appxvstrue(res::FEMSolution;savefile="",title="PDE Solution",
      appxTitle="Approximated Solution",trueTitle="True Solution",
      cmap=PyPlot.get_cmap("winter"))
  fig = PyPlot.figure("pyplot_appx_vs_true",figsize=(10,10))
  PyPlot.subplot(211,projection="3d")
  PyPlot.plot_trisurf(res.femMesh.node[:,1],res.femMesh.node[:,2],res.u,cmap=cmap)
  PyPlot.title(appxTitle)
  PyPlot.subplot(212,projection="3d")
  PyPlot.plot_trisurf(res.femMesh.node[:,1],res.femMesh.node[:,2],res.uTrue,cmap=cmap)
  PyPlot.title(trueTitle)
  PyPlot.suptitle(title)
  if savefile!=""
    PyPlot.savefig(savefile)
  end
  #=
  Plots.subplot([node[:,1] node[:,2] u node[:,1] node[:,2] uTrue],n=2,title=title,t=[:surface :surface])
  surface(node[:,1],node[:,2],u,cmap=PyPlot.get_cmap("winter"),title="Approximated Solution")
  surface(node[:,1],node[:,2],uTrue,cmap=PyPlot.get_cmap("winter"),title="True Solution")
  if savefile!=""
    Plots.savefig(savefile)
  end
  =#
end

"""
solplot(res::FEMSolution)

Plots the approximate solution and, if available, the true solution.

### Keyword Arguments

* `savefile`: Designates a file to save the plot in. Save type is defined by the chosen
extension. Default is "" which implies no saving.
* `title`: The title at the top of the plot. Default is "PDE Solution".
* `appxTitle`: The title above the approximate solution. Default is "Approximated Solution".
* `trueTitle`: The title above the true solution. Default is "True Solution".
* `cmap`: Specifies the color map in the plot. Default is PyPlot.get_cmap("winter").
"""
function solplot(res::FEMSolution;savefile="",title="PDE Solution",
        appxTitle="Approximated Solution",trueTitle="True Solution",
        cmap=PyPlot.get_cmap("winter"))
  if res.trueKnown
    solplot_appxvstrue(res,savefile=savefile,title=title,appxTitle=appxTitle,
    trueTitle=trueTitle,cmap=cmap)
  else
    solplot_appx(res,savefile=savefile,title=title,cmap=cmap)
  end
end

"""
solplot_appx(res::FEMSolution)

Plots the approximate solution.

### Keyword Arguments

* `savefile`: Designates a file to save the plot in. Save type is defined by the chosen
extension. Default is "" which implies no saving.
* `title`: The title at the top of the plot. Default is "PDE Solution".
* `cmap`: Specifies the color map in the plot. Default is PyPlot.get_cmap("winter").
"""
function solplot_appx(res::FEMSolution;savefile="",title="Approximated Solution",
  cmap=PyPlot.get_cmap("winter"))
  Plots.surface(res.femMesh.node[:,1],res.femMesh.node[:,2],res.u,cmap=PyPlot.get_cmap("winter"),title=title)
  if savefile!=""
    Plots.savefig(savefile)
  end
end

"""
showmesh(femMesh::FEMmesh)

Shows the mesh which is defined by the (node,elem) structure.

### Keyword Arguments

* `cmap`: Specifies the color map in the plot. Default is PyPlot.get_cmap("winter").
"""
function showmesh(femMesh::Mesh;cmap=PyPlot.get_cmap("ocean"))
  @unpack femMesh: node, elem
  dim = size(node,2)
  nv = size(elem,2)
  if (dim==2) && (nv==3) # planar triangulation

    h = PyPlot.plot_trisurf(node[:,1],node[:,2],zeros(size(node,1)),cmap=cmap)
  end

  #Not supported yet
  if (dim==2) && (nv==4) # planar quadrilateration
    C = 0.5*zeros(length(node[:,1]),length(node[:,2]))
    D = zeros(length(node[:,1]))
    p = pcolormesh(node[:,1],node[:,2],C,edgecolor=D,cmap=cmap)
  end
  if (dim==3)
    if size(elem,2) == 3 # surface meshes
      plot_trisurf(node[:,1],node[:,2],node[:,3],cmap=cmap)
    elseif size(elem,2) == 4
      mxcall(:showmesh3,0,node,elem); # PyPlot and GNUplot do not have tetrahedron plotting
      return
    end
  end
end

"""
convplot(measure,err)

Makes a convergence plot of err vs measure in loglog scale.

### Keyword Arguments

* `ErrStr`: The y-axis label. Default is "Error".
* `measureStr`: The x-axis label. Default is "Measure".
* `titleStr`: The title. Default is "\$ErrStr vs \$measureStr Convergence Plot".
"""
function convplot(measure,err;ErrStr="Error",measureStr="Measure",titleStr="$ErrStr vs $measureStr Convergence Plot")
  PyPlot.loglog(measure,err)
  PyPlot.title(titleStr)
  PyPlot.xlabel(measureStr)
  PyPlot.ylabel(ErrStr)
  #=
  Plots.plot(measure,err)
  Plots.title!(titleStr)
  Plots.xaxis!(measureStr,:log10)
  Plots.yaxis!(ErrStr,:log10)
  =#
end

"""
convplot_fullΔt(simres::ConvergenceSimulation)

Plots a grid which shows the H1, L2, nodal maximum, and nodal l2 error convergence
over changes of Δt.

### Keyword Arguments

* `savefile`: Designates a file to save the plot in. Save type is defined by the chosen
extension. Default is "" which implies no saving.
"""
function convplot_fullΔt(simres::ConvergenceSimulation;titleStr="All Convergences",savefile="")
  fig = PyPlot.figure("pyplot_appx_vs_true",figsize=(10,10))
  PyPlot.subplot(221)
  convplot_h1vsΔt(simres)
  PyPlot.subplot(222)
  convplot_l2vsΔt(simres)
  PyPlot.subplot(223)
  convplot_node2vsΔt(simres)
  PyPlot.subplot(224)
  convplot_maxvsΔt(simres)
  PyPlot.suptitle(titleStr)
  if savefile!=""
    PyPlot.savefig(savefile)
  end
  return(fig)
  #=
  Plots.subplot(convplot_h1vsΔx(simres),
  convplot_l2vsΔx(simres),
  convplot_node2vsΔx(simres),
  convplot_maxvsΔx(simres),
  layout=[2,2])
  if savefile!=""
    Plots.savefig(savefile)
  end
  =#
end

"""
convplot_fullΔx(simres::ConvergenceSimulation)

Plots a grid which shows the H1, L2, nodal maximum, and nodal l2 error convergence
over changes of Δx.

### Keyword Arguments

* `savefile`: Designates a file to save the plot in. Save type is defined by the chosen
extension. Default is "" which implies no saving.
"""
function convplot_fullΔx(simres::ConvergenceSimulation;titleStr="All Convergences",savefile="")
  fig = PyPlot.figure("pyplot_appx_vs_true",figsize=(10,10))
  PyPlot.subplot(221)
  convplot_h1vsΔx(simres)
  PyPlot.subplot(222)
  convplot_l2vsΔx(simres)
  PyPlot.subplot(223)
  convplot_node2vsΔx(simres)
  PyPlot.subplot(224)
  convplot_maxvsΔx(simres)
  PyPlot.suptitle(titleStr)
  if savefile!=""
    PyPlot.savefig(savefile)
  end
  return(fig)
end

"""
convplot_h1vsΔt(simres::ConvergenceSimulation)

Shows the H1 error convergence over changes of Δt.
"""
convplot_h1vsΔt(simres::ConvergenceSimulation) = convplot(simres.Δts,simres.h1Errors;ErrStr="H1 Error",measureStr=L"$\Delta t$")

"""
convplot_l2vsΔt(simres::ConvergenceSimulation)

Shows the L2 error convergence over changes of Δt.
"""
convplot_l2vsΔt(simres::ConvergenceSimulation) = convplot(simres.Δts,simres.l2Errors;ErrStr="L2 Error",measureStr=L"$\Delta t$")

"""
convplot_node2vsΔt(simres::ConvergenceSimulation)

Shows the nodal l2 error convergence over changes of Δt.
"""
convplot_node2vsΔt(simres::ConvergenceSimulation) = convplot(simres.Δts,simres.node2Errors;ErrStr="Nodal L2 Error",measureStr=L"$\Delta t$")

"""
convplot_maxvsΔt(simres::ConvergenceSimulation)

Shows the nodal maximum error convergence over changes of Δt.
"""
convplot_maxvsΔt(simres::ConvergenceSimulation) = convplot(simres.Δts,simres.maxErrors;ErrStr="Nodal Max Error",measureStr=L"$\Delta t$")

"""
convplot_h1vsΔx(simres::ConvergenceSimulation)

Shows the H1 error convergence over changes of Δx.
"""
convplot_h1vsΔx(simres::ConvergenceSimulation) = convplot(simres.Δxs,simres.h1Errors;ErrStr="H1 Error",measureStr=L"$\Delta x$")

"""
convplot_l2vsΔx(simres::ConvergenceSimulation)

Shows the nodal l2 error convergence over changes of Δx.
"""
convplot_l2vsΔx(simres::ConvergenceSimulation) = convplot(simres.Δxs,simres.l2Errors;ErrStr="L2 Error",measureStr=L"$\Delta x$")

"""
convplot_node2vsΔx(simres::ConvergenceSimulation)

Shows the nodal l2 error convergence over changes of Δx.
"""
convplot_node2vsΔx(simres::ConvergenceSimulation) = convplot(simres.Δxs,simres.node2Errors;ErrStr="Nodal L2 Error",measureStr=L"$\Delta x$")

"""
convplot_maxvsΔx(simres::ConvergenceSimulation)

Shows the nodal maximum error convergence over changes of Δx.
"""
convplot_maxvsΔx(simres::ConvergenceSimulation) = convplot(simres.Δxs,simres.maxErrors;ErrStr="Nodal Max Error",measureStr=L"$\Delta x$")
