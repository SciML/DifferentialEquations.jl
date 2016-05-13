#plot(z, zlim = (should_override ? zlim_override : default(:zlim)))

function solplot_animation(node,uFull;zlim=(0,1),vmax=1,cbar=true)
  Plots.pyplot(reuse=true,size=(750,750))
  @gif for i=1:size(uFull,2)
      surface(node[:,1],node[:,2],uFull[:,i],zlim=zlim,vmax=vmax,cbar=cbar)
  end
end

solplot_animation(res::FEMSolution;zlim=(0,1),vmax=1,cbar=true) = solplot_animation(res.femMesh.node,res.uFull,zlim=zlim,vmax=vmax,cbar=cbar)


function solplot_appxvstrue(node,u,uTrue;savefile="",title="PDE Solution",appxTitle="Approximated Solution",trueTitle="True Solution",cmap=PyPlot.get_cmap("winter"))
  fig = PyPlot.figure("pyplot_appx_vs_true",figsize=(10,10))
  PyPlot.subplot(211,projection="3d")
  PyPlot.plot_trisurf(node[:,1],node[:,2],u,cmap=cmap)
  PyPlot.title(appxTitle)
  PyPlot.subplot(212,projection="3d")
  PyPlot.plot_trisurf(node[:,1],node[:,2],uTrue,cmap=cmap)
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

solplot_appxvstrue(res::FEMSolution;savefile="",title="PDE Solution",
      appxTitle="Approximated Solution",trueTitle="True Solution",
      cmap=PyPlot.get_cmap("winter")) = solplot_appxvstrue(res.femMesh.node,
      res.u,res.uTrue,savefile=savefile,title=title,appxTitle=appxTitle,
      trueTitle=trueTitle,cmap=cmap)

function solplot(res::FEMSolution;savefile="",title="PDE Solution",
        appxTitle="Approximated Solution",trueTitle="True Solution",
        cmap=PyPlot.get_cmap("winter"))
  if res.trueKnown
    solplot_appxvstrue(res,savefile=savefile,title=title,appxTitle=appxTitle,
    trueTitle=trueTitle,cmap=cmap)
  else
    solplot_appx(res,savefile=savefile,title=title,appxTitle=appxTitle,
    cmap=cmap)
  end
end

function solplot_appx(node,u;savefile="",title="Approximated Solution",appxTitle="Approximated Solution",
  cmap=PyPlot.get_cmap("winter"))
  Plots.surface(node[:,1],node[:,2],u,cmap=PyPlot.get_cmap("winter"),title=title)
  if savefile!=""
    Plots.savefig(savefile)
  end
end

solplot_appx(res::FEMSolution;savefile="",title="Approximated Solution",appxTitle="Approximated Solution",
  cmap=PyPlot.get_cmap("winter")) = solplot_appx(res.femMesh.node,res.u,savefile=savefile,
  title=title,appxTitle=appxTitle,cmap=cmap)

function showmesh(node,elem,cmap=PyPlot.get_cmap("ocean"))
  dim = size(node,2)
  nv = size(elem,2)
  if (dim==2) && (nv==3) # planar triangulation
    h = surface(node[:,1],node[:,2],zeros(size(node,1)),cmap=cmap)
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

showmesh(femMesh) = showmesh(femMesh.node,femMesh.elem)

"""
convplot(measure,err;ErrStr="Error",measureStr="Measure",titleStr="$ErrStr vs $measureStr Convergence Plot")
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
convplot_fullΔt(simres::ConvergenceSimulation;titleStr="All Convergences",savefile="")
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
convplot_fullΔx(simres::ConvergenceSimulation;titleStr="All Convergences",savefile="")
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
"""
convplot_h1vsΔt(simres::ConvergenceSimulation) = convplot(simres.Δts,simres.h1Errors;ErrStr="H1 Error",measureStr=L"$\Delta t$")

"""
convplot_l2vsΔt(simres::ConvergenceSimulation)
"""
convplot_l2vsΔt(simres::ConvergenceSimulation) = convplot(simres.Δts,simres.l2Errors;ErrStr="L2 Error",measureStr=L"$\Delta t$")

"""
convplot_node2vsΔt(simres::ConvergenceSimulation)
"""
convplot_node2vsΔt(simres::ConvergenceSimulation) = convplot(simres.Δts,simres.node2Errors;ErrStr="Nodal L2 Error",measureStr=L"$\Delta t$")

"""
convplot_maxvsΔt(simres::ConvergenceSimulation)
"""
convplot_maxvsΔt(simres::ConvergenceSimulation) = convplot(simres.Δts,simres.maxErrors;ErrStr="Nodal Max Error",measureStr=L"$\Delta t$")

"""
convplot_h1vsΔx(simres::ConvergenceSimulation)
"""
convplot_h1vsΔx(simres::ConvergenceSimulation) = convplot(simres.Δxs,simres.h1Errors;ErrStr="H1 Error",measureStr=L"$\Delta x$")

"""
convplot_l2vsΔx(simres::ConvergenceSimulation)
"""
convplot_l2vsΔx(simres::ConvergenceSimulation) = convplot(simres.Δxs,simres.l2Errors;ErrStr="L2 Error",measureStr=L"$\Delta x$")

"""
convplot_node2vsΔx(simres::ConvergenceSimulation)
"""
convplot_node2vsΔx(simres::ConvergenceSimulation) = convplot(simres.Δxs,simres.node2Errors;ErrStr="Nodal L2 Error",measureStr=L"$\Delta x$")

"""
convplot_maxvsΔx(simres::ConvergenceSimulation)
"""
convplot_maxvsΔx(simres::ConvergenceSimulation) = convplot(simres.Δxs,simres.maxErrors;ErrStr="Nodal Max Error",measureStr=L"$\Delta x$")
