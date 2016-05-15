"""
getL2error(node,elem,uexact,uh,quad𝒪=[])

getL2error(femMesh::FEMmesh,sol,u)

Estimates the L2 error between uexact and uh on the mesh (node,elem). It
reads the mesh to estimate the element type and uses this to choose a
quadrature 𝒪 unless specified.
"""
function getL2error(node,elem,uexact,uh,quad𝒪=[])

Nu = length(uh);    N = size(node,1);   NT = size(elem,1)
# Euler formula: N-NE+NT = c
NE = N + NT - 1;    NP2 = N + NE;   NP3 = N + 2*NE + NT

if Nu > N+NT-5
    elem2dof = dofP2(elem)
    NP2 = max(vec(elem2dof))
    NE = NP2 - N
    NP3 = N+2*NE+NT
end


## Default quadrature 𝒪's for different elements
if isempty(quad𝒪)
    if Nu==NT # piecewise constant function P0
            quad𝒪 = 2
    elseif Nu==N      # piecewise linear function P1 element
            quad𝒪 = 3
    elseif Nu==NE     # piecewise linear function CR element
            quad𝒪 = 3
    elseif Nu==N+NT   # piecewise linear function + constant function
            quad𝒪 = 3
    elseif Nu==NP2    # piecewise quadratic function
            quad𝒪 = 4
    elseif Nu==NE+NT  # weak Galerkin element
            quad𝒪 = 3
    elseif Nu==NP3    # P3 element
            quad𝒪 = 5
    end
end

## compute L2 error element-wise using quadrature rule with 𝒪 quad𝒪
err = zeros(NT)
λ,ω = quadpts(quad𝒪)
if Nu==N # P1 piecewise linear function
        ϕ = λ # linear bases
elseif Nu==N+NT # P1+P0
        ϕ = λ # linear bases
elseif Nu==NE  # CR nonconforming P1 element
        ϕ = 1-2*λ
        elem2edge = elem2dof[:,4:6] - N
elseif Nu==NP2 # P2 piecewise quadratic elements
        ϕ[:,1] =   λ[:,1].*(2*λ[:,1]-1)
        ϕ[:,2] =   λ[:,2].*(2*λ[:,2]-1)
        ϕ[:,3] =   λ[:,3].*(2*λ[:,3]-1)
        ϕ[:,4] = 4*λ[:,2].*λ[:,3]
        ϕ[:,5] = 4*λ[:,1].*λ[:,3]
        ϕ[:,6] = 4*λ[:,2].*λ[:,1]
elseif Nu==NE+NT  # weak Galerkin element
#             uhp = uh(1:NT) # only count the interior part
        ϕ = 1-2*λ
        elem2edge = elem2dof[:,4:6] - N + NT
elseif Nu==2*NE+NT+N # P3 piecewise cubic elements
        ϕ[:,1]  = 0.5*(3*λ[:,1]-1).*(3*λ[:,1]-2).*λ[:,1]
        ϕ[:,2]  = 0.5*(3*λ[:,2]-1).*(3*λ[:,2]-2).*λ[:,2]
        ϕ[:,3]  = 0.5*(3*λ[:,3]-1).*(3*λ[:,3]-2).*λ[:,3]
        ϕ[:,4]  = 9/2*λ[:,3].*λ[:,2].*(3*λ[:,2]-1)
        ϕ[:,5]  = 9/2*λ[:,3].*λ[:,2].*(3*λ[:,3]-1)
        ϕ[:,6]  = 9/2*λ[:,1].*λ[:,3].*(3*λ[:,3]-1)
        ϕ[:,7]  = 9/2*λ[:,1].*λ[:,3].*(3*λ[:,1]-1)
        ϕ[:,8]  = 9/2*λ[:,1].*λ[:,2].*(3*λ[:,1]-1)
        ϕ[:,9]  = 9/2*λ[:,1].*λ[:,2].*(3*λ[:,2]-1)
        ϕ[:,10] = 27* λ[:,1].*λ[:,2].*λ[:,3]
        elem2dof = dofP3(elem) #Not Implemented
end
nQuad = size(λ,1)
for p = 1:nQuad
    # evaluate uh at quadrature point
    if Nu==NT # P0 piecewise constant function
            uhp = uh
    elseif Nu==N    # P1 piecewise linear function
            uhp = uh[elem[:,1]]*ϕ[p,1] +
                  uh[elem[:,2]]*ϕ[p,2] +
                  uh[elem[:,3]]*ϕ[p,3]
    elseif Nu==N+NT # P1+P0
            uhp = uh[elem[:,1]]*ϕ[p,1] +
                  uh[elem[:,2]]*ϕ[p,2] +
                  uh[elem[:,3]]*ϕ[p,3]
            uhp = uhp + uh[N+1:end]
    elseif Nu==NE  # CR nonconforming P1 element
            uhp = uh[elem2edge[:,1]]*ϕ[p,1] +
                  uh[elem2edge[:,2]]*ϕ[p,2] +
                  uh[elem2edge[:,3]]*ϕ[p,3]
    elseif Nu==NP2 # P2 piecewise quadratic function
            uhp = uh[elem2dof[:,1]].*ϕ[p,1] +
                  uh[elem2dof[:,2]].*ϕ[p,2] +
                  uh[elem2dof[:,3]].*ϕ[p,3] +
                  uh[elem2dof[:,4]].*ϕ[p,4] +
                  uh[elem2dof[:,5]].*ϕ[p,5] +
                  uh[elem2dof[:,6]].*ϕ[p,6]
    elseif Nu==NP3
            uhp = uh[elem2dof[:,1]].*ϕ[p,1] +
                  uh[elem2dof[:,2]].*ϕ[p,2] +
                  uh[elem2dof[:,3]].*ϕ[p,3] +
                  uh[elem2dof[:,4]].*ϕ[p,4] +
                  uh[elem2dof[:,5]].*ϕ[p,5] +
                  uh[elem2dof[:,6]].*ϕ[p,6] +
                  uh[elem2dof[:,7]].*ϕ[p,7] +
                  uh[elem2dof[:,8]].*ϕ[p,8] +
                  uh[elem2dof[:,9]].*ϕ[p,9] +
                  uh[elem2dof[:,10]].*ϕ[p,10]
    elseif Nu==NE+NT  # weak Galerkin element
            #             uhp = uh(1:NT) # only count the interior part
            uhp = uh[elem2edge[:,1]]*ϕ[p,1] +
                  uh[elem2edge[:,2]]*ϕ[p,2] +
                  uh[elem2edge[:,3]]*ϕ[p,3]
    end
    # quadrature points in the x-y coordinate
    pxy = λ[p,1]*node[elem[:,1],:] +
        λ[p,2]*node[elem[:,2],:] +
        λ[p,3]*node[elem[:,3],:]
    err = err + ω[p]*(uexact(pxy) - uhp).^2
end
# Modify by area
ve2 = node[elem[:,1],:]-node[elem[:,3],:]
ve3 = node[elem[:,2],:]-node[elem[:,1],:]
area = 0.5*abs(-ve3[:,1].*ve2[:,2]+ve3[:,2].*ve2[:,1])
err = area.*err
err[isnan(err)] = 0 # singular values are excluded
err = sqrt(sum(err))
return(err)
end

"""
quadpts(𝒪)

Returns the quadrature points and ω's for and 𝒪 ### quadrature in 2D.

Reference:
David Dunavant. High degree efficient symmetrical Gaussian
quadrature rules for the triangle. International journal for numerical
methods in engineering. 21(6):1129--1148, 1985.
"""
function quadpts(𝒪)
  if 𝒪>9
      𝒪 = 9
  elseif 𝒪==1
          λ = [1/3 1/3 1/3]
          ω = 1
  elseif 𝒪==2
          λ = [2/3 1/6 1/6
                    1/6 2/3 1/6
                    1/6 1/6 2/3]
          ω = [1/3 1/3 1/3]
  elseif 𝒪==3
          λ = [1/3 1/3 1/3
                    0.6 0.2 0.2
                    0.2 0.6 0.2
                    0.2 0.2 0.6]
          ω = [-27/48 25/48 25/48 25/48]
  elseif 𝒪==4
          λ = [0.108103018168070 0.445948490915965 0.445948490915965
                    0.445948490915965 0.108103018168070 0.445948490915965
                    0.445948490915965 0.445948490915965 0.108103018168070
                    0.816847572980459 0.091576213509771 0.091576213509771
                    0.091576213509771 0.816847572980459 0.091576213509771
                    0.091576213509771 0.091576213509771 0.816847572980459]
          ω = [0.223381589678011 0.223381589678011 0.223381589678011
                    0.109951743655322 0.109951743655322 0.109951743655322]
  elseif 𝒪==5
          α1 = 0.059715871789770 ;     β1 = 0.470142064105115
          α2 = 0.797426985353087 ;     β2 = 0.101286507323456
          λ = [   1/3    1/3    1/3
                    α1  β1  β1
                     β1 α1  β1
                     β1  β1 α1
                    α2  β2  β2
                     β2 α2  β2
                     β2  β2 α2]
          ω = [0.225 0.132394152788506 0.132394152788506 0.132394152788506
              0.125939180544827 0.125939180544827 0.125939180544827]
  elseif 𝒪==6
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
          λ = [A[:,[1;2]] 1-sum(A[:,[1;2]],2)]
          ω = A[:,3]
  elseif 𝒪==7
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
              λ = [A[:,[1;2]] 1-sum(A[:,[1;2]],2)]
              ω = A[:,3]
  elseif 𝒪==8
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
              λ = [A[:,[1;2]] 1-sum(A[:,[1;2]],2)]
              ω = A[:,3]
  elseif 𝒪==9
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
              λ = [A[:,[1;2]] 1-sum(A[:,[1;2]],2)]
              ω = A[:,3]
  end
  return(λ,ω)
end

"""
quadpts1(𝒪)

References: 
Pavel Holoborodko: http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/
"""
function quadpts1(𝒪)
  numPts = ceil((𝒪+1)/2)
  if numPts > 10
     numPts = 10
  end
if numPts==1
          A = [0      2.0000000000000000000000000]
elseif numPts==2
          A = [0.5773502691896257645091488 	1.0000000000000000000000000
              -0.5773502691896257645091488 	1.0000000000000000000000000]
elseif numPts==3
          A = [0 	0.8888888888888888888888889
              0.7745966692414833770358531 	0.5555555555555555555555556
              -0.7745966692414833770358531 	0.5555555555555555555555556]
elseif numPts==4
          A = [0.3399810435848562648026658 	0.6521451548625461426269361
              0.8611363115940525752239465 	0.3478548451374538573730639
              -0.3399810435848562648026658 	0.6521451548625461426269361
              -0.8611363115940525752239465 	0.3478548451374538573730639]
elseif numPts==5
          A = [0 	                            0.5688888888888888888888889
              0.5384693101056830910363144 	0.4786286704993664680412915
              0.9061798459386639927976269 	0.2369268850561890875142640
              -0.5384693101056830910363144 	0.4786286704993664680412915
              -0.9061798459386639927976269 	0.2369268850561890875142640]
elseif numPts==6
          A = [0.2386191860831969086305017 	0.4679139345726910473898703
              0.6612093864662645136613996 	0.3607615730481386075698335
              0.9324695142031520278123016 	0.1713244923791703450402961
              -0.2386191860831969086305017 	0.4679139345726910473898703
              -0.6612093864662645136613996 	0.3607615730481386075698335
              -0.9324695142031520278123016 	0.1713244923791703450402961]
elseif numPts==7
          A = [0 	                            0.4179591836734693877551020
              0.4058451513773971669066064 	0.3818300505051189449503698
              0.7415311855993944398638648 	0.2797053914892766679014678
              0.9491079123427585245261897 	0.1294849661688696932706114
              -0.4058451513773971669066064 	0.3818300505051189449503698
              -0.7415311855993944398638648 	0.2797053914892766679014678
              -0.9491079123427585245261897 	0.1294849661688696932706114]
elseif numPts==8
          A = [0.1834346424956498049394761 	0.3626837833783619829651504
              0.5255324099163289858177390 	0.3137066458778872873379622
              0.7966664774136267395915539 	0.2223810344533744705443560
              0.9602898564975362316835609 	0.1012285362903762591525314
              -0.1834346424956498049394761 	0.3626837833783619829651504
              -0.5255324099163289858177390 	0.3137066458778872873379622
              -0.7966664774136267395915539 	0.2223810344533744705443560
              -0.9602898564975362316835609 	0.1012285362903762591525314]
elseif numPts==9
          A = [0 	                            0.3302393550012597631645251
              0.3242534234038089290385380 	0.3123470770400028400686304
              0.6133714327005903973087020 	0.2606106964029354623187429
              0.8360311073266357942994298 	0.1806481606948574040584720
              0.9681602395076260898355762 	0.0812743883615744119718922
              -0.3242534234038089290385380 	0.3123470770400028400686304
              -0.6133714327005903973087020 	0.2606106964029354623187429
              -0.8360311073266357942994298 	0.1806481606948574040584720
              -0.9681602395076260898355762 	0.0812743883615744119718922]
elseif numPts==10
          A = [0.1488743389816312108848260 	0.2955242247147528701738930
              0.4333953941292471907992659 	0.2692667193099963550912269
              0.6794095682990244062343274 	0.2190863625159820439955349
              0.8650633666889845107320967 	0.1494513491505805931457763
              0.9739065285171717200779640 	0.0666713443086881375935688
              -0.1488743389816312108848260 	0.2955242247147528701738930
              -0.4333953941292471907992659 	0.2692667193099963550912269
              -0.6794095682990244062343274 	0.2190863625159820439955349
              -0.8650633666889845107320967 	0.1494513491505805931457763
              -0.9739065285171717200779640 	0.0666713443086881375935688]
  end
  λ1 = (A[:,1]+1)/2
  λ2 = 1 - λ1
  λ = [λ1 λ2]
  ω = A[:,2]/2
  return(λ,ω)
end
"""
function getH1error(node,elem,Du,uh,K=[],quad𝒪=[])

getH1error(femMesh::FEMmesh,Du,u)

Estimates the H1 error between uexact and uh on the mesh (node,elem). It
reads the mesh to estimate the element type and uses this to choose a
quadrature 𝒪 unless specified. If K is specified then it is the
diffusion coefficient matrix.
"""
function getH1error(node,elem,Du,uh,K=[],quad𝒪=[])

  Nu = size(uh,1);    N = size(node,1);   NT = size(elem,1);
  # Euler formula N-NE+NT = c # rough estimateus using Euler formula
  NE = N + NT;    NP2 = N + NE;   NP3 = N + 2*NE + NT
  if Nu > N+NT-5   # Euler formula N-NE+NT = c
      elem2dof = dofP2(elem)
      NP2 = max(elem2dof(:))
      NE = NP2 - N
      NP3 = N+2*NE+NT
  end

  ## Default quadrature 𝒪s for different elements
  if isempty(quad𝒪)
      if Nu==NT     # piecewise constant vector (uh is Duh)
              quad𝒪 = 3
      elseif Nu==N      # piecewise linear function ℙ1 element
              quad𝒪 = 3
      elseif Nu==NE     # piecewise linear function CR element
              quad𝒪 = 3
      elseif Nu==NP2    # piecewise quadratic function
              quad𝒪 = 4
      elseif Nu==NE + NT # WG element
              quad𝒪 = 3
      elseif Nu==NP3    # ℙ3 element
              quad𝒪 = 5
      end
  end

  ## compute gradient of finite element function uh
  #Only ℙ1 Implemented
  if (size(uh,2) == 2) && (Nu == NT)      # uh is a piecewise constant vector
      Duh = uh
      area = abs(simplexvolume(node,elem))
  elseif size(uh,2) == 1   # scalar function uh
      if Nu==N      # piecewise linear function ℙ1 element
              Duh,area = gradu(node,elem,uh)
      elseif Nu==NE     # piecewise linear function CR element
              elem2edge = elem2dof(:,4:6) - N
              Duh,area = graduCR(node,elem,elem2edge,uh)
      elseif Nu==NE + NT # weak Galerkin element
              elem2edge = elem2dof(:,4:6) - N
              Duh,area = graduWG(node,elem,elem2edge,uh)
      elseif Nu==NP2    # piecewise quadratic function
              Dλ,area = gradbasis(node,elem)
      elseif Nu==NP3
              Dλ,area = gradbasis(node,elem)
              elem2dof = dofP3(elem)
      end
  end

  ## compute H1 error element-wise using quadrature rule with 𝒪 quad𝒪
  λ,ω = quadpts(quad𝒪)
  nQuad = size(λ,1)
  err = zeros(NT)
  for p = 1:nQuad
      pxy = λ[p,1]*node[elem[:,1],:] +
            λ[p,2]*node[elem[:,2],:] +
            λ[p,3]*node[elem[:,3],:]
      if Nu == NP2 # piecewise quadratic function
          Dϕp1 = (4*λ[p,1]-1).*Dλ[:,:,1]
          Dϕp2 = (4*λ[p,2]-1).*Dλ[:,:,2]
          Dϕp3 = (4*λ[p,3]-1).*Dλ[:,:,3]
          Dϕp4 = 4*(λ[p,2]*Dλ[:,:,3]+λ[p,3]*Dλ[:,:,2])
          Dϕp5 = 4*(λ[p,3]*Dλ[:,:,1]+λ[p,1]*Dλ[:,:,3])
          Dϕp6 = 4*(λ[p,1]*Dλ[:,:,2]+λ[p,2]*Dλ[:,:,1])
          Duh = repmat(uh(elem2dof[:,1]),1,2).*Dϕp1 +
                repmat(uh(elem2dof[:,2]),1,2).*Dϕp2 +
                repmat(uh(elem2dof[:,3]),1,2).*Dϕp3 +
                repmat(uh(elem2dof[:,4]),1,2).*Dϕp4 +
                repmat(uh(elem2dof[:,5]),1,2).*Dϕp5 +
                repmat(uh(elem2dof[:,6]),1,2).*Dϕp6
      end
      if Nu == NP3 # piecewise cubic function
          Dϕp1 = (27/2*λ(p,1)*λ(p,1)-9*λ(p,1)+1).*Dλ(:,:,1)
          Dϕp2 = (27/2*λ(p,2)*λ(p,2)-9*λ(p,2)+1).*Dλ(:,:,2)
          Dϕp3 = (27/2*λ(p,3)*λ(p,3)-9*λ(p,3)+1).*Dλ(:,:,3)
          Dϕp4 = 9/2*((3*λ(p,2)*λ(p,2)-λ(p,2)).*Dλ(:,:,3)+
                  λ(p,3)*(6*λ(p,2)-1).*Dλ(:,:,2))
          Dϕp5 = 9/2*((3*λ(p,3)*λ(p,3)-λ(p,3)).*Dλ(:,:,2)+
                   λ(p,2)*(6*λ(p,3)-1).*Dλ(:,:,3))
          Dϕp6 = 9/2*((3*λ(p,3)*λ(p,3)-λ(p,3)).*Dλ(:,:,1)+
                   λ(p,1)*(6*λ(p,3)-1).*Dλ(:,:,3))
          Dϕp7 = 9/2*((3*λ(p,1)*λ(p,1)-λ(p,1)).*Dλ(:,:,3)+
                   λ(p,3)*(6*λ(p,1)-1).*Dλ(:,:,1))
          Dϕp8 = 9/2*((3*λ(p,1)*λ(p,1)-λ(p,1)).*Dλ(:,:,2)+
                   λ(p,2)*(6*λ(p,1)-1).*Dλ(:,:,1))
          Dϕp9 = 9/2*((3*λ(p,2)*λ(p,2)-λ(p,2)).*Dλ(:,:,1)+
                   λ(p,1)*(6*λ(p,2)-1).*Dλ(:,:,2))
          Dϕp10= 27*(λ(p,1)*λ(p,2)*Dλ(:,:,3)+λ(p,1)*λ(p,3)*Dλ(:,:,2)+
                   λ(p,3)*λ(p,2)*Dλ(:,:,1))
          Duh = repmat(uh(elem2dof(:,1)),1,2).*Dϕp1 +
                repmat(uh(elem2dof(:,2)),1,2).*Dϕp2 +
                repmat(uh(elem2dof(:,3)),1,2).*Dϕp3 +
                repmat(uh(elem2dof(:,4)),1,2).*Dϕp4 +
                repmat(uh(elem2dof(:,5)),1,2).*Dϕp5 +
                repmat(uh(elem2dof(:,6)),1,2).*Dϕp6 +
                repmat(uh(elem2dof(:,7)),1,2).*Dϕp7 +
                repmat(uh(elem2dof(:,8)),1,2).*Dϕp8 +
                repmat(uh(elem2dof(:,9)),1,2).*Dϕp9 +
                repmat(uh(elem2dof(:,10)),1,2).*Dϕp10
      end
      if isa(K,Function)
          err = err + ω[p]*K(pxy).*sum((Du(pxy)-Duh).^2,2)
      else
          err = err + ω[p]*sum((Du(pxy)-Duh).^2,2)
      end
  end
  if !isempty(K) && isa(K,Vector{Number}) && size(K,1) == NT
      err = K.*err    # K is piecewise constant
  end
  err = area.*err
  err[isnan(err)] = 0 # singular values are excluded
  err = sqrt(sum(err))
  return(err)
end

"""
gradu(node,elem,u,Dλ=[])

Estimates the gradient of u on the mesh (node,elem)
"""
function gradu(node,elem,u,Dλ=[])
  if isempty(Dλ)
      Dλ,area = gradbasis(node,elem)
  end
  dudx =  u[elem[:,1]].*Dλ[:,1,1] + u[elem[:,2]].*Dλ[:,1,2] +
        u[elem[:,3]].*Dλ[:,1,3]
  dudy =  u[elem[:,1]].*Dλ[:,2,1] + u[elem[:,2]].*Dλ[:,2,2] +
        u[elem[:,3]].*Dλ[:,2,3]
  Du = [dudx dudy];
  return(Du,area,Dλ)
end

"""
gradbasis(node,elem)

Returns the gradient of the barycentric basis elements.
"""
function gradbasis(node,elem)
  NT = size(elem,1)
  Dλ = Array{Float64}(NT,2,3)

  ve1 = node[elem[:,3],:]-node[elem[:,2],:]
  ve2 = node[elem[:,1],:]-node[elem[:,3],:]
  ve3 = node[elem[:,2],:]-node[elem[:,1],:]
  area = 0.5*(-ve3[:,1].*ve2[:,2] + ve3[:,2].*ve2[:,1])
  Dλ[1:NT,:,3] = [-ve3[:,2]./(2*area) ve3[:,1]./(2*area)]
  Dλ[1:NT,:,1] = [-ve1[:,2]./(2*area) ve1[:,1]./(2*area)]
  Dλ[1:NT,:,2] = [-ve2[:,2]./(2*area) ve2[:,1]./(2*area)]

  # When not positive orientated, reverse the sign.
  idx = (area.<0)
  area[idx,:] = -area[idx,:]
  elemSign = ones(NT)
  elemSign[idx] = -1
  return(Dλ,area,elemSign)
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
