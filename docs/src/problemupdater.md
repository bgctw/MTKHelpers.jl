```@meta
CurrentModule = MTKHelpers
```
# ProblemUpdater
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

