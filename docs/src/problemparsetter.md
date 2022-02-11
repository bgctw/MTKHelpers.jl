```@meta
CurrentModule = MTKHelpers
```

# Setting Initial state and Parameters of a problem

Oftern one wants to change a subset of the initial
states,`u0`, and a subset of parameters,`p`, of an ODEProblem during an optimization.

Given `u0` and `p` being a sequence of scalars, 
and `popt` is a sequence of optimized parameters, which may include initial states,
the following class helps updating the corresponding positions in 
an ODEProblem.

```@example doc
# setting up a simple example composite system and problem
using ModelingToolkit, DifferentialEquations, LabelledArrays
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
popt = SLVector(m₊x=0.1, m₊p1=2.1)
ps = ProblemParSetter(sys, keys(popt))

# extract optimized 
get_paropt_labeled(ps, prob)

# update states and parameters
prob2 = update_statepar(ps, popt, prob)
get_paropt_labeled(ps, prob2) == popt
```

## ProblemParSetter

The following type stores the necessary information, that can be queried.
```@docs
ProblemParSetter
statesyms(::ProblemParSetter)
count_states
```

## Updating a problem 
```@docs
update_statepar
```

Further a variant or merge is implemented to produce an updated
version of an `SLArray` or `LArray` are produced given a `NamedTuple` of updates.

```@example 
using LabelledArrays, MTKHelpers
popt = SLVector(L = 10.1, k_L = 1.1, k_R = 1/20.1)
popt2 = merge(popt, (k_L = 1.2,))
popt2.k_L == 1.2
```

## Extracting optimized parameters
```@docs
get_paropt
```

## Labeling 
```@docs
label_state
```






