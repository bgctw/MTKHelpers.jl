```@meta
CurrentModule = MTKHelpers
```

# Translating symbols and Nums

ModelingToolkit constructs and AbstractODEProblem from an ODESystem by supplying 
Dictionaries that map Nums to values. 
However, it is more convenient to store initial states and parameters
as ComponentVectors with symbolic keys, instead of Dictionaries with Num keys.

First, lets create an example system that can demonstrate symbolic arrays in states and 
parameters.

```@example doc
# setting up a simple example composite system and problem
using ModelingToolkit, OrdinaryDiffEq, ComponentArrays
using MTKHelpers
function samplesystem_vec(; name, τ = 3.0, i=0.1, p = [1.1, 1.2, 1.3])
    n_comp = 2
    @parameters t
    D = Differential(t)
    @variables x(..)[1:n_comp] 
    ps = @parameters τ=τ i=i p[1:3]=p 
    sts = [x(t)[i] for i in 1:n_comp]
    eq = [
        D(x(t)[1]) ~ i - p[1] * x(t)[1] + (p[2] - x(t)[1]^2) / τ, 
        D(x(t)[2]) ~ i - p[3] * x(t)[2], 
     ]
     ODESystem(eq, t, sts, vcat(ps...); name)
end
@named m = samplesystem_vec()
@named sys = embed_system(m)
nothing # hide
```

Next, we want to create the problem and need to specify the initial states
and parameters as Dictionaries Num -> value. However, the Nums are only 
defined within above function `samplesystem_vec()`.

MTKHelpers provides method [`system_num_dict`](@ref).

```@example doc
# setup initial conditions and parameters as ComponentVectors
u0 = ComponentVector(m₊x=[0.1, 0.2])
p = ComponentVector(m₊p=[2.1, 2.2, 2.3], m₊τ = 3.1) # keep i to default

# convert to Dict(Num -> value) in order to create the AbstractODEProblem
p_numdict = system_num_dict(p, sys)
u0_numdict = system_num_dict(u0, sys)
prob = ODEProblem(sys, u0_numdict, (0,2), p_numdict);

# check the parameters of the created problem
pset = ODEProblemParSetterConcrete(sys, ComponentVector()) 
p_prob = label_par(pset, prob.p)
p_prob.m₊i == 0.1     # from default
p_prob.m₊τ == p.m₊τ   # from p
p_prob.m₊p == p.m₊p   # from p
```

MTKHelpers provides a `remake` variant (API below), that allows to set `u0` and `p`
of an `AbstractODEProblem` by providing a ComponentVector of parameters or initial states. 
Note that this is slower than using the position-matching
of `ODEProblemParSetter`.

```@example doc
paropt = ComponentVector(state=(m₊x=[10.1,10.2],), par=(m₊τ = 10.1,))
prob2 = remake(prob, paropt)
p2_prob = label_par(pset, prob2.p); u2_prob = label_state(pset, prob2.u0)
p2_prob.m₊τ == paropt.par.m₊τ     # from paropt
u2_prob.m₊x == paropt.state.m₊x   # from paropt
p2_prob.m₊p == p_prob.m₊p         # from original prob
nothing # hide
```

## API

```@docs
system_num_dict
remake(prob::AbstractODEProblem, paropt::ComponentVector)
get_system
base_num
get_base_num_dict
```

