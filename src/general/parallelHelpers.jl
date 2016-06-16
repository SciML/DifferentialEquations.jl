function sendtosimple(p::Int, nm, val)
    ref = @spawnat(p, eval(Main, Expr(:(=), nm, val)))
end

macro sendto(p, nm, val)
    return :( sendtosimple($p, $nm, $val) )
end

function sendto(p::Int; args...)
    for (nm, val) in args
        @spawnat(p, eval(Main, Expr(:(=), nm, val)))
    end
end

getfrom(p::Int, nm::Symbol; mod=Main) = fetch(@spawnat(p, getfield(mod, nm)))


function passobj(src::Int, target::Vector{Int}, nm::Symbol;
                 from_mod=Main, to_mod=Main)
    r = RemoteRef(src)
    @spawnat(src, put!(r, getfield(from_mod, nm)))
    for to in target
        @spawnat(to, eval(to_mod, Expr(:(=), nm, fetch(r))))
    end
    nothing
end


function passobj(src::Int, target::Int, nm::Symbol; from_mod=Main, to_mod=Main)
    passobj(src, [target], nm; from_mod=from_mod, to_mod=to_mod)
end


function passobj(src::Int, target, nms::Vector{Symbol};
                 from_mod=Main, to_mod=Main)
    for nm in nms
        passobj(src, target, nm; from_mod=from_mod, to_mod=to_mod)
    end
end

function sendto(ps::Vector{Int}; args...)
    for p in ps
        sendto(p; args...)
    end
end

macro broadcast(nm, val)
    quote
    @sync for p in workers()
        @async sendtosimple(p, $nm, $val)
    end
    end
end

function shapeResult(res)
  #Input is a vector of tuples
  #Output the "columns" as vectors
  out = cell(length(res[1]))
  for i = 1:length(res[1])
    out[i] = Vector{typeof(res[1][i])}(length(res))
  end
  for i = 1:length(res), j = 1:length(out)
    out[j][i] = res[i][j]
  end
  return tuple(out...)
end

function monteCarlo(func::Function,args...)
  res = pmap(func,args...)
  if length(res[1])==1
    return res
  else
    return shapeResult(res)
  end
end

#=
# creates an integer x and Matrix y on processes 1 and 2
sendto([1, 2], x=100, y=rand(2, 3))

# create a variable here, then send it everywhere else
z = randn(10, 10); sendto(workers(), z=z)

# get an object from named x from Main module on process 2. Name it x
x = getfrom(2, :x)

# pass variable named x from process 2 to all other processes
passobj(2, filter(x->x!=2, procs()), :x)

# pass variables t, u, v from process 3 to process 1
passobj(3, 1, [:t, :u, :v])

# Pass a variable from the `Foo` module on process 1 to Main on workers
passobj(1, workers(), [:foo]; from_mod=Foo)
=#
