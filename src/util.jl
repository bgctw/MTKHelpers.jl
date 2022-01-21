
"""
embed_system(m;name)

Embeds system m as the single component of a larger system.

Then the naming of the states, parameters, observables matches
the namespace of the system.

```julia
#using DifferentialEquations, Plots
function samplesystem(;name) 
@variables t 
D = Differential(t) 
sts = @variables x(t) RHS(t)  # RHS is observed
ps = @parameters τ       # parameters
ODESystem([ RHS  ~ (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)
end                     
@named m = samplesystem()
@named sys = embed_system(m)
# note that state key is m.x, hence only m needs to be defined
prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0), [m.τ => 3.0])
sol = solve(prob);
plot(sol, vars=[m.x,m.RHS])    
plot(sol, vars=getproperty.(Ref(m),[:x, :RHS])) # access by symbols
```
"""
function embed_system(m;name)
@named _sys_embed = ODESystem(Equation[], ModelingToolkit.get_iv(m))               
sys = compose(_sys_embed, [m]; name)
structural_simplify(sys)
end

"Extract the basic symbol"
symbol(t::Term) = symbol(t.f)
symbol(num::Num) = symbol(num.val)
symbol(s) = Symbol(s)

"Omit the part before the first dot"
strip_namespace(s::Symbol) = Symbol(strip_namespace(string(s)))
strip_namespace(s::String) = match(r"[^.₊]+$",s).match

"Extract the basic symbols without namespace of system states"

statesyms(sys::ODESystem) = symbol.(states(sys))

"Extract the basic symbols without namespace of system parameters"

parsyms(sys::ODESystem) = symbol.(parameters(sys))

