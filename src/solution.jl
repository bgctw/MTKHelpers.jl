"""
    getlast(sol, vars...)
    getlast(sol, vars_vec) 

Apply getindex for each var to sol and return a NamedArray.
    
```jldocstest; output=false
using ModelingToolkit, OrdinaryDiffEq
using MTKHelpers
using NamedArrays
using ModelingToolkit: t_nounits as t, D_nounits as D
function samplesystem(;name) 
    sts = @variables x(t) RHS(t)  # RHS is observed
    ps = @parameters τ       # parameters
    ODESystem([ RHS  ~ (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)
end                     
@named m = samplesystem()
@named sys = embed_system(m)
prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0), [m.τ => 3.0])
sol = solve(prob, Tsit5());
res = getlast(sol, m.x, m.RHS)
res == NamedArray([sol[m.x,end], sol[m.RHS,end]], ([m.x, m.RHS],))   
# output
true
```    
"""
function getlast(sol::SciMLBase.AbstractODESolution, vars...; kwargs...)
    getlast(sol, collect(vars))
end
function getlast(sol::SciMLBase.AbstractODESolution, vars::AbstractVector)
    iend = lastindex(sol, 2)
    names_vars_vec = _names_vars(vars)
    #names_vars_vec = names_vars #convert(Array, names_vars)::Vector{eltype(names_vars)}
    a = getindex.(Ref(sol), names_vars_vec, iend)::Array{eltype(sol)}
    #names = ntuple(i -> symbol_op(vars[i]), length(vars))
    #A = @LArray a names
    NamedArrays.NamedArray(a, (names_vars_vec,))
end
_names_vars(vars_vec) = vars_vec
_names_vars(vars_vec::NamedArrays.NamedArray) = first(names(vars_vec))
