"""
getL2error(node,elem,uexact,uh,quadOrder=[])

getL2error(femMesh::FEMmesh,sol,u)

Estimates the L2 error between uexact and uh on the mesh (node,elem). It
reads the mesh to estimate the element type and uses this to choose a
quadrature order unless specified.
"""
function getL2error(node,elem,uexact,uh,quadOrder=[])
## GETL2ERROR L2 norm of the approximation error.

## Number of vertices, elements, edges, and degrees of freedom
Nu = length(uh);    N = size(node,1);   NT = size(elem,1)
# Euler formula N-NE+NT = c
NE = N + NT - 1;    NP2 = N + NE;   NP3 = N + 2*NE + NT

if Nu > N+NT-5
    elem2dof = dofP2(elem)
    NP2 = max(elem2dof(:))
    NE = NP2 - N
    NP3 = N+2*NE+NT
end


## Default quadrature orders for different elements
if isempty(quadOrder)
    if Nu==NT # piecewise constant function P0
            quadOrder = 2
    elseif Nu==N      # piecewise linear function P1 element
            quadOrder = 3
    elseif Nu==NE     # piecewise linear function CR element
            quadOrder = 3
    elseif Nu==N+NT   # piecewise linear function + constant function
            quadOrder = 3
    elseif Nu==NP2    # piecewise quadratic function
            quadOrder = 4
    elseif Nu==NE+NT  # weak Galerkin element
            quadOrder = 3
    elseif Nu==NP3    # P3 element
            quadOrder = 5
    end
end

## compute L2 error element-wise using quadrature rule with order quadOrder
err = zeros(NT)
lambda,weight = quadpts(quadOrder)
# basis function at quadrature points
if Nu==N # P1 piecewise linear function
        phi = lambda # linear bases
elseif Nu==N+NT # P1+P0
        phi = lambda # linear bases
elseif Nu==NE  # CR nonconforming P1 element
        phi = 1-2*lambda
        elem2edge = elem2dof[:,4:6] - N
elseif Nu==NP2 # P2 piecewise quadratic elements
        phi[:,1] =   lambda[:,1].*(2*lambda[:,1]-1)
        phi[:,2] =   lambda[:,2].*(2*lambda[:,2]-1)
        phi[:,3] =   lambda[:,3].*(2*lambda[:,3]-1)
        phi[:,4] = 4*lambda[:,2].*lambda[:,3]
        phi[:,5] = 4*lambda[:,1].*lambda[:,3]
        phi[:,6] = 4*lambda[:,2].*lambda[:,1]
elseif Nu==NE+NT  # weak Galerkin element
#             uhp = uh(1:NT) # only count the interior part
        phi = 1-2*lambda
        elem2edge = elem2dof[:,4:6] - N + NT
elseif Nu==2*NE+NT+N # P3 piecewise cubic elements
        phi[:,1]  = 0.5*(3*lambda[:,1]-1).*(3*lambda[:,1]-2).*lambda[:,1]
        phi[:,2]  = 0.5*(3*lambda[:,2]-1).*(3*lambda[:,2]-2).*lambda[:,2]
        phi[:,3]  = 0.5*(3*lambda[:,3]-1).*(3*lambda[:,3]-2).*lambda[:,3]
        phi[:,4]  = 9/2*lambda[:,3].*lambda[:,2].*(3*lambda[:,2]-1)
        phi[:,5]  = 9/2*lambda[:,3].*lambda[:,2].*(3*lambda[:,3]-1)
        phi[:,6]  = 9/2*lambda[:,1].*lambda[:,3].*(3*lambda[:,3]-1)
        phi[:,7]  = 9/2*lambda[:,1].*lambda[:,3].*(3*lambda[:,1]-1)
        phi[:,8]  = 9/2*lambda[:,1].*lambda[:,2].*(3*lambda[:,1]-1)
        phi[:,9]  = 9/2*lambda[:,1].*lambda[:,2].*(3*lambda[:,2]-1)
        phi[:,10] = 27* lambda[:,1].*lambda[:,2].*lambda[:,3]
        elem2dof = dofP3(elem) #Not Implemented
end
nQuad = size(lambda,1)
for p = 1:nQuad
    # evaluate uh at quadrature point
    if Nu==NT # P0 piecewise constant function
            uhp = uh
    elseif Nu==N    # P1 piecewise linear function
            uhp = uh[elem[:,1]]*phi[p,1] +
                  uh[elem[:,2]]*phi[p,2] +
                  uh[elem[:,3]]*phi[p,3]
    elseif Nu==N+NT # P1+P0
            uhp = uh[elem[:,1]]*phi[p,1] +
                  uh[elem[:,2]]*phi[p,2] +
                  uh[elem[:,3]]*phi[p,3]
            uhp = uhp + uh[N+1:end]
    elseif Nu==NE  # CR nonconforming P1 element
            uhp = uh[elem2edge[:,1]]*phi[p,1] +
                  uh[elem2edge[:,2]]*phi[p,2] +
                  uh[elem2edge[:,3]]*phi[p,3]
    elseif Nu==NP2 # P2 piecewise quadratic function
            uhp = uh[elem2dof[:,1]].*phi[p,1] +
                  uh[elem2dof[:,2]].*phi[p,2] +
                  uh[elem2dof[:,3]].*phi[p,3] +
                  uh[elem2dof[:,4]].*phi[p,4] +
                  uh[elem2dof[:,5]].*phi[p,5] +
                  uh[elem2dof[:,6]].*phi[p,6]
    elseif Nu==NP3
            uhp = uh[elem2dof[:,1]].*phi[p,1] +
                  uh[elem2dof[:,2]].*phi[p,2] +
                  uh[elem2dof[:,3]].*phi[p,3] +
                  uh[elem2dof[:,4]].*phi[p,4] +
                  uh[elem2dof[:,5]].*phi[p,5] +
                  uh[elem2dof[:,6]].*phi[p,6] +
                  uh[elem2dof[:,7]].*phi[p,7] +
                  uh[elem2dof[:,8]].*phi[p,8] +
                  uh[elem2dof[:,9]].*phi[p,9] +
                  uh[elem2dof[:,10]].*phi[p,10]
    elseif Nu==NE+NT  # weak Galerkin element
            #             uhp = uh(1:NT) # only count the interior part
            uhp = uh[elem2edge[:,1]]*phi[p,1] +
                  uh[elem2edge[:,2]]*phi[p,2] +
                  uh[elem2edge[:,3]]*phi[p,3]
    end
    # quadrature points in the x-y coordinate
    pxy = lambda[p,1]*node[elem[:,1],:] +
        lambda[p,2]*node[elem[:,2],:] +
        lambda[p,3]*node[elem[:,3],:]
    err = err + weight[p]*(uexact(pxy) - uhp).^2
end
# Modification
# area of triangles
ve2 = node[elem[:,1],:]-node[elem[:,3],:]
ve3 = node[elem[:,2],:]-node[elem[:,1],:]
area = 0.5*abs(-ve3[:,1].*ve2[:,2]+ve3[:,2].*ve2[:,1])
err = area.*err
err[isnan(err)] = 0 # singular values, i.e. uexact(p) = infty, are excluded
err = sqrt(sum(err))
return(err)
end

"""
quadpts(order)

Returns the quadrature points and weights for and order ### quadrature in 2D.
"""
function quadpts(order)
## QUADPTS quadrature points in 2-D.
  if order>9
      order = 9
  elseif order==1     # Order 1, nQuad 1
          lambda = [1/3 1/3 1/3]
          weight = 1
  elseif order==2     # Order 2, nQuad 3
          lambda = [2/3 1/6 1/6
                    1/6 2/3 1/6
                    1/6 1/6 2/3]
          weight = [1/3 1/3 1/3]
  elseif order==3     # Order 3, nQuad 4
          lambda = [1/3 1/3 1/3
                    0.6 0.2 0.2
                    0.2 0.6 0.2
                    0.2 0.2 0.6]
          weight = [-27/48 25/48 25/48 25/48]
  elseif order==4     # Order 4, nQuad 6
          lambda = [0.108103018168070 0.445948490915965 0.445948490915965
                    0.445948490915965 0.108103018168070 0.445948490915965
                    0.445948490915965 0.445948490915965 0.108103018168070
                    0.816847572980459 0.091576213509771 0.091576213509771
                    0.091576213509771 0.816847572980459 0.091576213509771
                    0.091576213509771 0.091576213509771 0.816847572980459]
          weight = [0.223381589678011 0.223381589678011 0.223381589678011
                    0.109951743655322 0.109951743655322 0.109951743655322]
  elseif order==5     # Order 5, nQuad 7
          alpha1 = 0.059715871789770 ;     beta1 = 0.470142064105115
          alpha2 = 0.797426985353087 ;     beta2 = 0.101286507323456
          lambda = [   1/3    1/3    1/3
                    alpha1  beta1  beta1
                     beta1 alpha1  beta1
                     beta1  beta1 alpha1
                    alpha2  beta2  beta2
                     beta2 alpha2  beta2
                     beta2  beta2 alpha2]
          weight = [0.225 0.132394152788506 0.132394152788506 0.132394152788506
              0.125939180544827 0.125939180544827 0.125939180544827]
  elseif order==6
          A =[0.249286745170910  0.249286745170910  0.116786275726379
              0.249286745170910  0.501426509658179  0.116786275726379
              0.501426509658179  0.249286745170910  0.116786275726379
              0.063089014491502  0.063089014491502  0.050844906370207
              0.063089014491502  0.873821971016996  0.050844906370207
              0.873821971016996  0.063089014491502  0.050844906370207
              0.310352451033784  0.636502499121399  0.082851075618374
              0.636502499121399  0.053145049844817  0.082851075618374
              0.053145049844817  0.310352451033784  0.082851075618374
              0.636502499121399  0.310352451033784  0.082851075618374
              0.310352451033784  0.053145049844817  0.082851075618374
              0.053145049844817  0.636502499121399  0.082851075618374]
          lambda = [A[:,[1;2]] 1-sum(A[:,[1;2]],2)]
          weight = A[:,3]
  elseif order==7
          A =[0.333333333333333  0.333333333333333 -0.149570044467682
              0.260345966079040  0.260345966079040  0.175615257433208
              0.260345966079040  0.479308067841920  0.175615257433208
              0.479308067841920  0.260345966079040  0.175615257433208
              0.065130102902216  0.065130102902216  0.053347235608838
              0.065130102902216  0.869739794195568  0.053347235608838
              0.869739794195568  0.065130102902216  0.053347235608838
              0.312865496004874  0.638444188569810  0.077113760890257
              0.638444188569810  0.048690315425316  0.077113760890257
              0.048690315425316  0.312865496004874  0.077113760890257
              0.638444188569810  0.312865496004874  0.077113760890257
              0.312865496004874  0.048690315425316  0.077113760890257
              0.048690315425316  0.638444188569810  0.077113760890257]
              lambda = [A[:,[1;2]] 1-sum(A[:,[1;2]],2)]
              weight = A[:,3]
  elseif order==8
          A =[0.333333333333333  0.333333333333333  0.144315607677787
              0.081414823414554  0.459292588292723  0.095091634267285
              0.459292588292723  0.081414823414554  0.095091634267285
              0.459292588292723  0.459292588292723  0.095091634267285
              0.658861384496480  0.170569307751760  0.103217370534718
              0.170569307751760  0.658861384496480  0.103217370534718
              0.170569307751760  0.170569307751760  0.103217370534718
              0.898905543365938  0.050547228317031  0.032458497623198
              0.050547228317031  0.898905543365938  0.032458497623198
              0.050547228317031  0.050547228317031  0.032458497623198
              0.008394777409958  0.263112829634638  0.027230314174435
              0.008394777409958  0.728492392955404  0.027230314174435
              0.263112829634638  0.008394777409958  0.027230314174435
              0.728492392955404  0.008394777409958  0.027230314174435
              0.263112829634638  0.728492392955404  0.027230314174435
              0.728492392955404  0.263112829634638  0.027230314174435]
              lambda = [A[:,[1;2]] 1-sum(A[:,[1;2]],2)]
              weight = A[:,3]
  elseif order==9
          A =[0.333333333333333  0.333333333333333  0.097135796282799
              0.020634961602525  0.489682519198738  0.031334700227139
              0.489682519198738  0.020634961602525  0.031334700227139
              0.489682519198738  0.489682519198738  0.031334700227139
              0.125820817014127  0.437089591492937  0.07782754100474
              0.437089591492937  0.125820817014127  0.07782754100474
              0.437089591492937  0.437089591492937  0.07782754100474
              0.623592928761935  0.188203535619033  0.079647738927210
              0.188203535619033  0.623592928761935  0.079647738927210
              0.188203535619033  0.188203535619033  0.079647738927210
              0.910540973211095  0.044729513394453  0.025577675658698
              0.044729513394453  0.910540973211095  0.025577675658698
              0.044729513394453  0.044729513394453  0.025577675658698
              0.036838412054736  0.221962989160766  0.043283539377289
              0.036838412054736  0.741198598784498  0.043283539377289
              0.221962989160766  0.036838412054736  0.043283539377289
              0.741198598784498  0.036838412054736  0.043283539377289
              0.221962989160766  0.741198598784498  0.043283539377289
              0.741198598784498  0.221962989160766  0.043283539377289]
              lambda = [A[:,[1;2]] 1-sum(A[:,[1;2]],2)]
              weight = A[:,3]
  end
  return(lambda,weight)
end

"""
function getH1error(node,elem,Du,uh,K=[],quadOrder=[])

getH1error(femMesh::FEMmesh,Du,u)

Estimates the H1 error between uexact and uh on the mesh (node,elem). It
reads the mesh to estimate the element type and uses this to choose a
quadrature order unless specified. If K is specified then it is the
diffusion coefficient matrix.
"""
function getH1error(node,elem,Du,uh,K=[],quadOrder=[])

  Nu = size(uh,1);    N = size(node,1);   NT = size(elem,1);
  # Euler formula N-NE+NT = c # rough estimateus using Euler formula
  NE = N + NT;    NP2 = N + NE;   NP3 = N + 2*NE + NT
  if Nu > N+NT-5   # Euler formula N-NE+NT = c
      elem2dof = dofP2(elem)
      NP2 = max(elem2dof(:))
      NE = NP2 - N
      NP3 = N+2*NE+NT
  end

  ## Default quadrature orders for different elements
  if isempty(quadOrder)
      if Nu==NT     # piecewise constant vector (uh is Duh)
              quadOrder = 3
      elseif Nu==N      # piecewise linear function P1 element
              quadOrder = 3
      elseif Nu==NE     # piecewise linear function CR element
              quadOrder = 3
      elseif Nu==NP2    # piecewise quadratic function
              quadOrder = 4
      elseif Nu==NE + NT # WG element
              quadOrder = 3
      elseif Nu==NP3    # P3 element
              quadOrder = 5
      end
  end

  ## compute gradient of finite element function uh
  #Only P1 Implemented
  if (size(uh,2) == 2) && (Nu == NT)      # uh is a piecewise constant vector
      Duh = uh
      area = abs(simplexvolume(node,elem))
  elseif size(uh,2) == 1   # scalar function uh
      if Nu==N      # piecewise linear function P1 element
              Duh,area = gradu(node,elem,uh)
      elseif Nu==NE     # piecewise linear function CR element
              elem2edge = elem2dof(:,4:6) - N
              Duh,area = graduCR(node,elem,elem2edge,uh)
      elseif Nu==NE + NT # weak Galerkin element
              elem2edge = elem2dof(:,4:6) - N
              Duh,area = graduWG(node,elem,elem2edge,uh)
      elseif Nu==NP2    # piecewise quadratic function
              Dlambda,area = gradbasis(node,elem)
      elseif Nu==NP3
              Dlambda,area = gradbasis(node,elem)
              elem2dof = dofP3(elem)
      end
  end

  ## compute H1 error element-wise using quadrature rule with order quadOrder
  lambda,weight = quadpts(quadOrder)
  nQuad = size(lambda,1)
  err = zeros(NT)
  for p = 1:nQuad
      pxy = lambda[p,1]*node[elem[:,1],:] +
            lambda[p,2]*node[elem[:,2],:] +
            lambda[p,3]*node[elem[:,3],:]
      if Nu == NP2 # piecewise quadratic function
          Dphip1 = (4*lambda[p,1]-1).*Dlambda[:,:,1]
          Dphip2 = (4*lambda[p,2]-1).*Dlambda[:,:,2]
          Dphip3 = (4*lambda[p,3]-1).*Dlambda[:,:,3]
          Dphip4 = 4*(lambda[p,2]*Dlambda[:,:,3]+lambda[p,3]*Dlambda[:,:,2])
          Dphip5 = 4*(lambda[p,3]*Dlambda[:,:,1]+lambda[p,1]*Dlambda[:,:,3])
          Dphip6 = 4*(lambda[p,1]*Dlambda[:,:,2]+lambda[p,2]*Dlambda[:,:,1])
          Duh = repmat(uh(elem2dof[:,1]),1,2).*Dphip1 +
                repmat(uh(elem2dof[:,2]),1,2).*Dphip2 +
                repmat(uh(elem2dof[:,3]),1,2).*Dphip3 +
                repmat(uh(elem2dof[:,4]),1,2).*Dphip4 +
                repmat(uh(elem2dof[:,5]),1,2).*Dphip5 +
                repmat(uh(elem2dof[:,6]),1,2).*Dphip6
      end
      if Nu == NP3 # piecewise cubic function
          Dphip1 = (27/2*lambda(p,1)*lambda(p,1)-9*lambda(p,1)+1).*Dlambda(:,:,1)
          Dphip2 = (27/2*lambda(p,2)*lambda(p,2)-9*lambda(p,2)+1).*Dlambda(:,:,2)
          Dphip3 = (27/2*lambda(p,3)*lambda(p,3)-9*lambda(p,3)+1).*Dlambda(:,:,3)
          Dphip4 = 9/2*((3*lambda(p,2)*lambda(p,2)-lambda(p,2)).*Dlambda(:,:,3)+
                  lambda(p,3)*(6*lambda(p,2)-1).*Dlambda(:,:,2))
          Dphip5 = 9/2*((3*lambda(p,3)*lambda(p,3)-lambda(p,3)).*Dlambda(:,:,2)+
                   lambda(p,2)*(6*lambda(p,3)-1).*Dlambda(:,:,3))
          Dphip6 = 9/2*((3*lambda(p,3)*lambda(p,3)-lambda(p,3)).*Dlambda(:,:,1)+
                   lambda(p,1)*(6*lambda(p,3)-1).*Dlambda(:,:,3))
          Dphip7 = 9/2*((3*lambda(p,1)*lambda(p,1)-lambda(p,1)).*Dlambda(:,:,3)+
                   lambda(p,3)*(6*lambda(p,1)-1).*Dlambda(:,:,1))
          Dphip8 = 9/2*((3*lambda(p,1)*lambda(p,1)-lambda(p,1)).*Dlambda(:,:,2)+
                   lambda(p,2)*(6*lambda(p,1)-1).*Dlambda(:,:,1))
          Dphip9 = 9/2*((3*lambda(p,2)*lambda(p,2)-lambda(p,2)).*Dlambda(:,:,1)+
                   lambda(p,1)*(6*lambda(p,2)-1).*Dlambda(:,:,2))
          Dphip10= 27*(lambda(p,1)*lambda(p,2)*Dlambda(:,:,3)+lambda(p,1)*lambda(p,3)*Dlambda(:,:,2)+
                   lambda(p,3)*lambda(p,2)*Dlambda(:,:,1))
          Duh = repmat(uh(elem2dof(:,1)),1,2).*Dphip1 +
                repmat(uh(elem2dof(:,2)),1,2).*Dphip2 +
                repmat(uh(elem2dof(:,3)),1,2).*Dphip3 +
                repmat(uh(elem2dof(:,4)),1,2).*Dphip4 +
                repmat(uh(elem2dof(:,5)),1,2).*Dphip5 +
                repmat(uh(elem2dof(:,6)),1,2).*Dphip6 +
                repmat(uh(elem2dof(:,7)),1,2).*Dphip7 +
                repmat(uh(elem2dof(:,8)),1,2).*Dphip8 +
                repmat(uh(elem2dof(:,9)),1,2).*Dphip9 +
                repmat(uh(elem2dof(:,10)),1,2).*Dphip10
      end
      if ~isempty(K) && ~isnumeric(K) # K is a function
          err = err + weight[p]*K(pxy).*sum((Du(pxy)-Duh).^2,2)
      else
          err = err + weight[p]*sum((Du(pxy)-Duh).^2,2)
      end
  end
  if ~isempty(K) && isnumeric(K) && size(K,1) == NT
      err = K.*err    # K is piecewise constant
  end
  err = area.*err
  err[isnan(err)] = 0 # singular values are excluded
  err = sqrt(sum(err))
  return(err)
end

"""
gradu(node,elem,u,Dlambda=[])

Estimates the gradient of u on the mesh (node,elem)
"""
function gradu(node,elem,u,Dlambda=[])
  ## GRADU gradient of a finite element function.
  if isempty(Dlambda)
      Dlambda,area = gradbasis(node,elem)
  end
  dudx =  u[elem[:,1]].*Dlambda[:,1,1] + u[elem[:,2]].*Dlambda[:,1,2] +
        u[elem[:,3]].*Dlambda[:,1,3]
  dudy =  u[elem[:,1]].*Dlambda[:,2,1] + u[elem[:,2]].*Dlambda[:,2,2] +
        u[elem[:,3]].*Dlambda[:,2,3]
  Du = [dudx dudy];
  return(Du,area,Dlambda)
end

"""
gradbasis(node,elem)

Returns the gradient of the barycentric basis elements.
"""
function gradbasis(node,elem)
  ## GRADBASIS gradient of barycentric basis.
  NT = size(elem,1)
  Dlambda = Array{Float64}(NT,2,3)
  # $\nabla \phi_i = rotation(l_i)/(2|\tau|)$
  ve1 = node[elem[:,3],:]-node[elem[:,2],:]
  ve2 = node[elem[:,1],:]-node[elem[:,3],:]
  ve3 = node[elem[:,2],:]-node[elem[:,1],:]
  area = 0.5*(-ve3[:,1].*ve2[:,2] + ve3[:,2].*ve2[:,1])
  Dlambda[1:NT,:,3] = [-ve3[:,2]./(2*area) ve3[:,1]./(2*area)]
  Dlambda[1:NT,:,1] = [-ve1[:,2]./(2*area) ve1[:,1]./(2*area)]
  Dlambda[1:NT,:,2] = [-ve2[:,2]./(2*area) ve2[:,1]./(2*area)]

  # When the triangle is not positive orientated, we reverse the sign of the
  # area. The sign of Dlambda is always right since signed area is used in
  # the computation.
  idx = (area.<0)
  area[idx,:] = -area[idx,:]
  elemSign = ones(NT)
  elemSign[idx] = -1
  return(Dlambda,area,elemSign)
end

function getL2error(femMesh::FEMmesh,sol,u)
  if femMesh.evolutionEq
    return(getL2error(femMesh.node,femMesh.elem,x->sol(x,femMesh.T),u))
  else
    return(getL2error(femMesh.node,femMesh.elem,sol,u))
  end
end

function getH1error(femMesh::FEMmesh,Du,u)
  if femMesh.evolutionEq
    return(getH1error(femMesh.node,femMesh.elem,x->Du(x,femMesh.T),u))
  else
    return(getH1error(femMesh.node,femMesh.elem,Du,u))
  end
end
