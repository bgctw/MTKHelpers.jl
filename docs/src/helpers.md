```@meta
CurrentModule = MTKHelpers
```
# Translating between parameters and Problem
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

# ProblemUpdater

The functioonality of ProblemUpdater can be better modeled since MTK version 10
using [`bindings`](https://docs.sciml.ai/ModelingToolkit/dev/tutorials/initialization/#bindings_and_ics).

## Bindings
First, create an example system and example problem.
```@example doc
using ModelingToolkit, OrdinaryDiffEq, ComponentArrays
using ModelingToolkit: ModelingToolkit as MTK
using MTKHelpers
using ModelingToolkit: t_nounits as t, D_nounits as D
function get_sys1()
    sts = @variables L(t)
    ps = @parameters k_L, k_R, k_P
    eq = [D(L) ~ 0]
    sys1 = System(eq, t, sts, vcat(ps...); name = :sys1)
end
sys1 = mtkcompile(get_sys1())
ps1 = [sys1.L => 10.0, sys1.k_L => 1.0, sys1.k_R => 1.0/20.0, sys1.k_P => 2.0]
prob = ODEProblem(sys1, ps1, (0.0, 1.0))
nothing # hide
```

The bindings cannot be changed after the system is created. 
MTKHelper provides function [`override_system`](@ref) that allows
to create a new system with several properties changes. Note
that this is has been developed and tested only for basic ODESystems.

```@example doc
sys2 = mtkcompile(override_system(sys1, 
    bindings = [sys1.k_R => sys1.k_L / 20.0, sys1.k_P => sys1.k_L * 2.0]))
ps2 = [sys2.L => 10.0, sys2.k_L => 1.0]
prob2 = ODEProblem(sys2, ps2, (0.0, 1.0))
prob3 = remake(prob, p = [sys2.k_L => 1.1])
prob.ps[sys2.k_P], prob2.ps[sys2.k_P]
```

If [initial conditions are bound](https://docs.sciml.ai/ModelingToolkit/dev/tutorials/initialization/#Initialization-By-Example:-The-Cartesian-Pendulum), 
then initial conditions could be specified for
the bound property instead of the original one.

## Alternative: Problemupdater

Bindings were introduced in MTK10. Before, a ProblemUpdater could be used
to take care of updating dependent parameters.

A [`ODEProblemParSetterConcrete`](@ref) can be combined with a [`KeysProblemParGetter`](@ref)
or another specific implementation of [`AbstractProblemParGetter`](@ref) to 
update an AbstractODEProblem based on information already present in the AbstractODEProblem.

The following example updates parameters `k_R` and `k_P` in the AbstractODEProblem
to the value of `k_L`. This can be useful to ensure that these parameters
are also changed when optimizing parameter `k_L`.

An implementations of `AbstractProblemParGetter` can use any computation of
the source keys to provide its destination keys. It should implement the keys method,
so that when constructing the `ProblemUpdater`, consistent keys are used,
as in the example below.


Next, setup a `ProblemUpdater`, `pu`, and apply it to the problem via `prob2 = pu(prob)`.
```@example doc
mapping = (:k_L => :k_R, :k_L => :k_P)
pu = get_ode_problemupdater(KeysProblemParGetter(mapping, prob), get_system(prob))
prob2 = pu(prob)
p2 = get_par_labeled(par_setter(pu), prob2)
p2.k_P == p2.k_R == Dict(ps1)[sys1.k_L]
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
Abstract supertype of [`ODEProblemParSetter`](@ref)

```@docs
AbstractProblemParSetter
```

