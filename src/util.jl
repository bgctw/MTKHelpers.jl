
"""
    embed_system(m;name)

Embeds system `m` as the single component of a larger system.
This helps to match the naming of the states, parameters, observables to 
the namespace of the system. 

```jldocstest; output=false
using ModelingToolkit, DifferentialEquations
using MTKHelpers
function samplesystem(;name) 
    @variables t 
    D = Differential(t) 
    sts = @variables x(t) RHS(t)  # RHS is observed
    ps = @parameters τ       # parameters
    ODESystem([ RHS  ~ (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)
end                     
@named m = samplesystem()
@named sys = embed_system(m)

# note that keys are `m.x`,`m.τ` or `m.RHS`.
# Hence only m needs to be defined rather than all the states, parameters, 
# and observables of the system.
prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0), [m.τ => 3.0])
sol = solve(prob);
dx1 = sol[m.RHS][1:3] 
dx2 = sol[getproperty(m,:RHS)][1:3]  # access by symbols

# using Plots
# plot(sol, vars=[m.x,m.RHS])    
# plot(sol, vars=getproperty.(Ref(m),[:x, :RHS])) # access by symbols
length(dx2) == 3
# output
true
```
"""
function embed_system(m;name)
    @named _sys_embed = ODESystem(Equation[], ModelingToolkit.get_iv(m))               
    sys = compose(_sys_embed, [m]; name)
    structural_simplify(sys)
end

"""
    symbol(t)

Extract the inner symbol from a Term or Num object.
"""
function symbol(t::Term); symbol(t.f); end
symbol(num::Num) = symbol(num.val)
symbol(s) = Symbol(s)

"""
    strip_namespace(s)
    
Omit the part before the first dot.
"""    
function strip_namespace(s::String); match(r"[^.₊]+$",s).match; end
function strip_namespace(s::Symbol); Symbol(strip_namespace(string(s))); end

"""
    statesyms(sys::ODYSystem)
    parsyms(sys::ODYSystem)

Extract the basic symbols without namespace of system states and system parameters.
"""
function statesyms(sys::ODESystem); symbol.(states(sys)); end
parsyms(sys::ODESystem) = symbol.(parameters(sys))

