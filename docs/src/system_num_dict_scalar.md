```@meta
CurrentModule = MTKHelpers
```

# Translating symbols and Nums

ModelingToolkit constructs an AbstractODEProblem from an ODESystem by supplying 
Dictionaries that map Nums to values. 
However, it is more convenient to store initial states and parameters
as ComponentVectors with symbolic keys, instead of Dictionaries with Num keys.

First, lets create an example system that can demonstrate symbolic arrays in states and 
parameters.

```@example doc
# setting up a simple example composite system and problem
using ModelingToolkit, OrdinaryDiffEq, ComponentArrays
using MTKHelpers
using ModelingToolkit: t_nounits as t, D_nounits as D
function samplesystem_scalar(; name, τ = 3.0, i=0.1, p1 = 1.1, p2 = 1.2, p3 = 1.3)
    @variables x1(t) x2(t)
    ps = @parameters τ=τ i=i p1=p1 p2=p2 p3=p3
    sts = [x1, x2]
    eq = [
        D(x1) ~ i - p1 * x1 + (p2 - x1^2) / τ, 
        D(x2) ~ i - p3 * x2, 
     ]
     ODESystem(eq, t, sts, vcat(ps...); name)
end
@named m = samplesystem_scalar()
@named sys = embed_system(m)
nothing # hide
```

Next, we want to create the problem and need to specify the initial states
and parameters as Dictionaries Num -> value. However, the Nums are only 
defined within above function `samplesystem_vec()`.

MTKHelpers provides method [`system_num_dict`](@ref).

```@example doc
# setup initial conditions and parameters as ComponentVectors
u0 = ComponentVector(m₊x1=0.1, m₊x2=0.2)
p = ComponentVector(m₊p1=2.1, m₊p2=2.2, m₊p3=2.3, m₊τ = 3.1) # keep i to default

# convert to Dict(Num -> value) in order to create the AbstractODEProblem
p_numdict = system_num_dict(p, sys)
u0_numdict = system_num_dict(u0, sys)
prob = ODEProblem(sys, merge(u0_numdict, p_numdict), (0,2));

# check the parameters of the created problem
pset = ODEProblemParSetter(sys, ComponentVector()) 
p_prob = get_par_labeled(pset, prob)
p_prob.m₊i == 0.1     # from default
p_prob.m₊τ == p.m₊τ   # from p
p_prob.m₊p1 == p.m₊p1   # from p
nothing # hide
```

MTKHelpers provides a [`remake`](@ref) variant with a [`ODEProblemParSetter`](@ref), 
that allows to update a subset of `u0` and `p`
of an `AbstractODEProblem` by providing a ComponentVector of parameters or initial states. 

```@example doc
paropt = ComponentVector(state=(m₊x1=10.1, m₊x2=10.2), par=(m₊τ = 10.1,))
pset = ODEProblemParSetter(get_system(prob), paropt)
prob2 = remake(prob, paropt, pset)
p2_prob = get_par_labeled(pset, prob2); u2_prob = get_state_labeled(pset, prob2)
p2_prob.m₊τ == paropt.par.m₊τ     # from paropt
p2_prob.m₊p1 == p_prob.m₊p1         # from original prob
u2_prob.m₊x1 == paropt.state.m₊x1
nothing # hide
```

Note that the order of entries symbolic-array components in the state of the problem
may differ from the consecutive order in the ComponentArray.
Functions `remake(prob, paropt, pset)` and `get_state(pset, prob)` take care 
to set and extract the components in correct order.

When updating an ODEProblem, MTKHelpers expects either a ComponentVector
with sub-vectors `state` and `par` or, alternatively, a ComponentVector with keys of 
states before keys of parameters.
Function [`vcat_statesfirst`](@ref) helps to ensure this order when concatenating
ComponentVectors.

## API

```@docs
system_num_dict
get_system
base_num
get_base_num_dict
vcat_statesfirst
```

