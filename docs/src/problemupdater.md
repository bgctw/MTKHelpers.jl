```@meta
CurrentModule = MTKHelpers
```
# ProblemUpdater

The functioonality of ProblemUpdater can be better modeled since MTK11
suing bindings.




A [`ODEProblemParSetterConcrete`](@ref) can be combined with a [`KeysProblemParGetter`](@ref)
or another specific implementation of [`AbstractProblemParGetter`](@ref) to 
update an AbstractODEProblem based on information already present in the AbstractODEProblem.

The following example updates parameters `k_R` and `k_P` in the AbstractODEProblem
to the value of `k_L`. This can be useful to ensure that these parameters
are also changed when optimizing parameter `k_L`.

An implementations of `AbstractProblemParGetter` can use any computation of
the source keys to provide its destination keys. It should implement the keys method,
so that when constructing the ProblemUpdater, consistent keys are used,
as in the example below.

First, create an example system and example problem.
```@example doc
using ModelingToolkit, OrdinaryDiffEq, ComponentArrays
using MTKHelpers
function get_sys1()
    sts = @variables L(t)
    ps = @parameters k_L, k_R, k_P
    eq = [D(L) ~ 0]
    sys1 = System(eq, t, sts, vcat(ps...); name = :sys1)
end
sys1 = mtkcompile(get_sys1())
u0 = ComponentVector(L = 10.0)
p = ComponentVector(k_L = 1.0, k_R = 1 / 20, k_P = 2.0)
prob = ODEProblem(sys1,
    get_system_symbol_dict(sys1, vcat(u0, p)), (0.0, 1.0))
nothing # hide
```

Next, setup a ProblemUpdater, `pu`, and apply it to the problem via `prob2 = pu(prob)`.
```@example doc
mapping = (:k_L => :k_R, :k_L => :k_P)
pu = get_ode_problemupdater(KeysProblemParGetter(mapping, keys(u0)), get_system(prob))
prob2 = pu(prob)
p2 = get_par_labeled(par_setter(pu), prob2)
p2.k_P == p.k_L
p2.k_R == p.k_L
nothing # hide
```

```@docs
ProblemUpdater
get_ode_problemupdater
NullProblemUpdater
```

# ProblemParGetter

In order to provide computations for parameters to set, declare a
concrete subtype of `AbstractProblemParGetter`, and implement a custom 
method `(pg::MyProblemParGetter)(pu::ProblemUpdater, prob)` that returns
a vector of parameter values. 

```@docs
AbstractProblemParGetter
```

One simple subtype of `AbstractProblemParGetter` is `KeysProblemParGetter`, 
which just extracts variables from the 
original problem to update other parameters.
It can be used to ensure that some parameter of a problem will always equal 
another parameter of the problem. 

```@docs
KeysProblemParGetter
```



# ProblemParSetter
```@docs
AbstractProblemParSetter
```

## Translating between parameters and Problem
```@docs
remake(::SciMLBase.AbstractSciMLProblem, popt, pset::AbstractProblemParSetter)
get_paropt(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem)
```

# Helper functions to access optimized parameters
```@docs
axis_paropt(::AbstractProblemParSetter)
classes_paropt(::AbstractProblemParSetter)
count_paropt(pset::AbstractProblemParSetter) 
keys_paropt(ps::AbstractProblemParSetter) 
symbols_paropt(pset::AbstractProblemParSetter)
```

# Labeling parameter vectors
```@docs
label_paropt(pset::AbstractProblemParSetter, popt)
name_paropt(pset::AbstractProblemParSetter, paropt::AbstractVector)
```

