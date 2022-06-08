```@meta
CurrentModule = MTKHelpers
```

# Setting Initial state and Parameters of a problem

Oftern one wants to change a subset of the initial
states,`u0`, and a subset of parameters,`p`, of an ODEProblem during an optimization.

Given `u0` and `p` can be expressed as ComponentVectors, 
and `popt` can be expressed as a ComponentVector of optimized parameters, 
which may include initial states,
the following class helps updating the corresponding positions in 
an ODEProblem.

Initial states and parameter components must have different names.

```@example doc
# setting up a simple example composite system and problem
using ModelingToolkit, DifferentialEquations, ComponentArrays
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
prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0))

# setup position matching
popt = ComponentVector(m₊x=0.1, m₊p1=2.1)
pset = ProblemParSetter(sys, popt) # not fully type-inferred

# extract optimized 
get_paropt(pset, prob)          # plain vector
get_paropt_labeled(pset, prob)  # ComponentVector
name_paropt(pset, prob)         # NamedVector 

# update states and parameters
prob2 = update_statepar(pset, popt, prob)
prob2.p # p is still a plain vector
label_par(pset, prob2.p).m₊p1 == popt.m₊p1 # attach labels and access properties
label_state(pset, prob2.u0).m₊x == popt.m₊x # attach labels and access properties
get_paropt_labeled(pset, prob2) == popt
```

Note that ProblemParSetter, `pset`, is only fully type-inferred when constructed with 
three ComponentArrays.Axis objects. This propagates to all ComponentVectors 
constructed by it, e.g. with `label_state`.
Hence, its recommended to pass `pset` across a function barrier for code
where performance matters.

## ProblemParSetter

The following type stores the necessary information, that can be queried.
```@docs
ProblemParSetter
count_state(::AbstractProblemParSetter)
axis_state(::AbstractProblemParSetter)
symbols_state(::AbstractProblemParSetter)
```

## Updating a problem 
```@docs
update_statepar(::AbstractProblemParSetter,popt, prob::ODEProblem)
```

## Extracting optimized parameters
```@docs
get_paropt(::AbstractProblemParSetter,  u0, p)
```

## Labeling 
```@docs
label_state(::AbstractProblemParSetter, u0)
name_state(::AbstractProblemParSetter, state::AbstractVector)
```






