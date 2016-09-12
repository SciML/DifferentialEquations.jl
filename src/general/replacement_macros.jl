macro def(name, definition)
    return quote
        macro $name()
            esc($(Expr(:quote, definition)))
        end
    end
end

macro ode_define(ex,params...)
  ## Build Symbol dictionary
  dict = Dict{Symbol,Int}()
  for i in 2:2:length(ex.args) #Every odd line is line number
    arg = ex.args[i].args[1] #Get the first thing, should be dsomething
    nodarg = Symbol(string(arg)[2:end]) #Take off the d
    if !haskey(dict,nodarg)
      s = string(arg)
      dict[Symbol(string(arg)[2:end])] = i/2 # and label it the next int if not seen before
    end
  end
  syms = keys(dict)
  pdict = Dict{Symbol,Any}()
  ## Build parameter dictionary
  for i in 1:length(params)
    pdict[params[i].args[1]] = params[i].args[2] # works for k=3, or k=>3
  end
  # Run find replace
  ode_findreplace(ex,dict,syms,pdict)
  # Return the lambda
  return :((t,u,du) -> $(ex))
end

function ode_findreplace(ex,dict,syms,pdict)
  for (i,arg) in enumerate(ex.args)
    if isa(arg,Expr)
      ode_findreplace(arg,dict,syms,pdict)
    elseif isa(arg,Symbol)
      s = string(arg)
      if haskey(dict,arg)
        ex.args[i] = :(u[$(dict[arg])])
      elseif haskey(pdict,arg)
        ex.args[i] = :($(pdict[arg]))
      elseif length(string(arg))>1 && haskey(dict,Symbol(s[nextind(s, 1):end])) && Symbol(s[1])==:d
        tmp = Symbol(s[nextind(s, 1):end]) # Remove the first letter, the d
        ex.args[i] = :(du[$(dict[tmp])])
      end
    end
  end
end

macro fem_define(sig,ex,params...)
  ## Build Symbol dictionary
  dict = Dict{Symbol,Int}()
  for (i,arg) in enumerate(ex.args)
    if i%2 == 0
      dict[Symbol(string(arg.args[1])[2:end])] = i/2 # Change du->u, Fix i counting
    end
  end
  syms = keys(dict)

  pdict = Dict{Symbol,Any}()
  ## Build parameter dictionary
  for i in 1:length(params)
    pdict[params[i].args[1]] = params[i].args[2] # works for k=3, or k=>3
  end
  # Run find replace
  fem_findreplace(ex,dict,syms,pdict)

  funcs = Vector{Expr}(0) # Get all of the functions
  for (i,arg) in enumerate(ex.args)
    if i%2 == 0
      push!(funcs,arg.args[2])
    end
  end
  ex = Expr(:hcat,funcs...)
  # Return the lambda
  :($(sig) -> $(ex))
end

function fem_findreplace(ex,dict,syms,pdict)
  for (i,arg) in enumerate(ex.args)
    if isa(arg,Expr)
      fem_findreplace(arg,dict,syms,pdict)
    elseif isa(arg,Symbol)
      if haskey(dict,arg)
        ex.args[i] = :(u[:,$(dict[arg])])
      elseif haskey(pdict,arg)
        ex.args[i] = :($(pdict[arg]))
      elseif haskey(FEM_SYMBOL_DICT,arg)
        ex.args[i] = FEM_SYMBOL_DICT[arg]
      end
    end
  end
end

FEM_SYMBOL_DICT = Dict{Symbol,Expr}(:x=>:(x[:,1]),:y=>:(x[:,2]),:z=>:(x[:,3]))
