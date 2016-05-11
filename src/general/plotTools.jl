function solplot_animation(node,animateFrames;zlim=(0,1),vmax=1,cbar=true)
  Plots.pyplot(reuse=true)
  @gif for i=1:size(animateFrames,2)
      surface(node[:,1],node[:,2],animateFrames[:,i],zlim=zlim,vmax=vmax,cbar=cbar)
  end
end

solplot_animation(res::FEMSolution;zlim=(0,1),vmax=1,cbar=true) = solplot_animation(res.femMesh.node,res.animateFrames,zlim=zlim,vmax=vmax,cbar=cbar)

function solplot_appxvstrue(node,u,uTrue;savefile="")
  fig = PyPlot.figure("pyplot_appx_vs_true",figsize=(10,10))
  PyPlot.subplot(211,projection="3d")
  PyPlot.plot_trisurf(node[:,1],node[:,2],u,cmap=PyPlot.get_cmap("winter"))
  PyPlot.title("Approximated Solution")
  PyPlot.subplot(212,projection="3d")
  PyPlot.plot_trisurf(node[:,1],node[:,2],uTrue,cmap=PyPlot.get_cmap("winter"))
  PyPlot.title("True Solution")
  PyPlot.suptitle("Approximated vs True Solution")
  if savefile!=""
    PyPlot.savefig(savefile)
  end
end

solplot_appxvstrue(res::FEMSolution;savefile="") = solplot_appxvstrue(res.femMesh.node,res.u,res.uTrue,savefile=savefile)

function solplot(res::FEMSolution;savefile="")
  if res.trueKnown
    solplot_appxvstrue(res,savefile=savefile)
  else
    solplot_appx(res,savefile=savefile)
  end
end

function solplot_appx(node,u;savefile="")
  fig = PyPlot.figure("pyplot_appx",figsize=(10,10))
  PyPlot.plot_trisurf(node[:,1],node[:,2],u,cmap=PyPlot.get_cmap("winter"))
  PyPlot.title("Approximated Solution")
  if savefile!=""
    PyPlot.savefig(savefile)
  end
end

solplot_appx(res::FEMSolution;savefile="") = solplot_appx(res.femMesh.node,res.u,savefile=savefile)

function showmesh(node,elem)
  dim = size(node,2)
  nv = size(elem,2)
  if (dim==2) && (nv==3) # planar triangulation
    h = plot_trisurf(node[:,1],node[:,2],zeros(size(node,1)),cmap=get_cmap("ocean"))
  end

  if (dim==2) && (nv==4) # planar quadrilateration
    C = 0.5*zeros(length(node[:,1]),length(node[:,2]))
    D = zeros(length(node[:,1]))
    p = pcolormesh(node[:,1],node[:,2],C,edgecolor=D,cmap="ocean")
  end
  if (dim==3)
    if size(elem,2) == 3 # surface meshes
      plot_trisurf(node[:,1],node[:,2],node[:,3],cmap=get_cmap("ocean"))
    elseif size(elem,2) == 4
      mxcall(:showmesh3,0,node,elem); # PyPlot and GNUplot do not have tetrahedron plotting
      return
    end
  end
end

showmesh(femMesh) = showmesh(femMesh.node,femMesh.elem)

function convplot(measure,err;ErrStr="Error",measureStr="Measure",titleStr="$ErrStr vs $measureStr Convergence Plot")
  PyPlot.loglog(measure,err)
  PyPlot.title(titleStr)
  PyPlot.xlabel(measureStr)
  PyPlot.ylabel(ErrStr)
end

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
end

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

convplot_h1vsΔt(simres::ConvergenceSimulation) = convplot(simres.Δts,simres.h1Errors;ErrStr="H1 Error",measureStr=L"$\Delta t$")

convplot_l2vsΔt(simres::ConvergenceSimulation) = convplot(simres.Δts,simres.l2Errors;ErrStr="L2 Error",measureStr=L"$\Delta t$")

convplot_node2vsΔt(simres::ConvergenceSimulation) = convplot(simres.Δts,simres.node2Errors;ErrStr="Nodal L2 Error",measureStr=L"$\Delta t$")

convplot_maxvsΔt(simres::ConvergenceSimulation) = convplot(simres.Δts,simres.maxErrors;ErrStr="Nodal Max Error",measureStr=L"$\Delta t$")

convplot_h1vsΔx(simres::ConvergenceSimulation) = convplot(simres.Δxs,simres.h1Errors;ErrStr="H1 Error",measureStr=L"$\Delta x$")

convplot_l2vsΔx(simres::ConvergenceSimulation) = convplot(simres.Δxs,simres.l2Errors;ErrStr="L2 Error",measureStr=L"$\Delta x$")

convplot_node2vsΔx(simres::ConvergenceSimulation) = convplot(simres.Δxs,simres.node2Errors;ErrStr="Nodal L2 Error",measureStr=L"$\Delta x$")

convplot_maxvsΔx(simres::ConvergenceSimulation) = convplot(simres.Δxs,simres.maxErrors;ErrStr="Nodal Max Error",measureStr=L"$\Delta x$")
