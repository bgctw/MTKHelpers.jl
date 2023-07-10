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

## ProblemUpdater
A [`ProblemParSetter`](@ref) can be combined with a [`KeysProblemParGetter`](@ref)
or other specific implementations of `AbstractProblemParGetter` to 
update an ODEProblem based on information already present in the ODEProblem.

The following example updates parameters `k_R` and `k_P` in the ODEProblem
to the value of `k_L`. This can be useful to ensure that these parameters
are also changed when optimizing parameter `k_L`.

An implementations of `AbstractProblemParGetter` can use any computation of
the source keys to provide its destination keys. It should implement the keys method,
so that when constructing the ProblemUpdater, consistent keys are used,
as in the example below.

```@example doc
f = (u,p,t) -> p[1]*u
u0 = (L=1/2,)
p = (k_L = 1.0, k_R = 2.0, k_P = 3.0)
tspan = (0.0,1.0)
prob = ODEProblem(f,collect(u0),tspan,collect(p))

source = (:k_L,:k_L)
dest = (:k_R,:k_P)
pu = ProblemUpdater(KeysProblemParGetter(source, dest), keys(u0), keys(p))
prob2 = pu(prob)
p2 = label_par(par_setter(pu), prob2.p)
p2.k_R == p.k_L
p2.k_P == p.k_L
```

```@docs
ProblemUpdater
```

## ProblemParGetter

```@docs
KeysProblemParGetter
```

In order to provide computations for parameters to set, declare a
concrete subtype of `KeysProblemParGetter`, and implement a custom 
method `(pg::MyProblemParGetter)(pu::ProblemUpdater, prob)` that returns
a vector of parameter values.

## ProblemParSetter

The following type stores the necessary information, that can be queried.
```@docs
ProblemParSetter
count_state
axis_state
symbols_state
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

## Setting entire state and parameter vectors
```@docs
get_u_map
```
