```@meta
CurrentModule = MTKHelpers
```

# Setting Parameters of a problem



The following type stores the necessary information, that can be queried.
```@docs
AbstractProblemParSetter
remake(::SciMLBase.AbstractSciMLProblem, popt, pset::AbstractProblemParSetter)
get_paropt(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem)
```

## Labeling 
```@docs
label_state(::AbstractProblemParSetter, u0)
name_state(::AbstractProblemParSetter, state::AbstractVector)
```

