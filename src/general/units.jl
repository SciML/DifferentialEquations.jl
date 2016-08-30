max{ND,N,m,kg,s,A,K,mol,cd,rad,sr}(u::Array{SIUnits.SIQuantity{N,m,kg,s,A,K,mol,cd,rad,sr},ND},utmp::Array{SIUnits.SIQuantity{N,m,kg,s,A,K,mol,cd,rad,sr},ND}) = map((x,y)->SIUnits.SIQuantity{N,m,kg,s,A,K,mol,cd,rad,sr}(max(x.val,y.val)),u,utmp)

function fem_squaremesh{N,m,kg,s,A,K,mol,cd,rad,sr}(square,h::SIUnits.SIQuantity{N,m,kg,s,A,K,mol,cd,rad,sr})
  node,elem = fem_squaremesh(square,h.val)
  node = map((x)->SIUnits.SIQuantity{N,m,kg,s,A,K,mol,cd,rad,sr}(x),node)
  return node,elem
end
