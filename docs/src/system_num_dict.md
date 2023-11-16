```@meta
CurrentModule = MTKHelpers
```

# Translating symbols and Nums

ModelingToolkit constructs and AbstractODEProblem from an ODESystem by supplying Dictionaries
that map Nums to values. 
However, it is more convenient to store initial states and parameters
as ComponentVectors with symbolic keys, instead of Dictionaries with Num keys.

The following functions help to translate symbols to Nums of the system.

```@example doc
# setting up a simple example composite system and problem
using ModelingToolkit, OrdinaryDiffEq, ComponentArrays
using MTKHelpers
function samplesystem(;name,τ = 3.0, p1=1.1, p2=1.2) 
    @variables t 
    D = Differential(t) 
    sts = @variables x(t) RHS(t)             # RHS is observed
    ps = @parameters τ=τ p1=p1 p2=p2       # parameters
    ODESystem([ RHS  ~ p1/p2 * (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)
end                     
@named m = samplesystem()
@named sys = embed_system(m)

# setup initial conditions and parameters as ComponentVectors
u0 = ComponentVector(m₊x=0.0)
p = ComponentVector(m₊p1=1.11, m₊τ = 3.1,)

# convert to Dict(Num -> value) in order to create the AbstractODEProblem
p_numdict = system_num_dict(p, m)
prob = ODEProblem(sys, system_num_dict(u0,m), (0,2), p_numdict);

pset = ODEProblemParSetterConcrete(sys, ComponentVector()) 
p_prob = label_par(pset, prob.p)
p_prob.m₊p2 = 1.2 # from default
p_prob.m₊τ == 3.1 # from p
```

## `system_num_dict`

```@docs
system_num_dict
get_system_symbol_dict
```

