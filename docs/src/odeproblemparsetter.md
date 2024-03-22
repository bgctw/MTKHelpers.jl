```@meta
CurrentModule = MTKHelpers
```

# AbstractODEProblemParSetter

The following type implements [Translating between parameters and Problem](@ref)
for AbstractODEProblems.

```@docs
AbstractODEProblemParSetter
```

See the specific implementation [`ODEProblemParSetter`](@ref). 


## Helper functions to access state and parameters
```@docs
axis_state
keys_state
count_state
symbols_state
```

## Labeling 
```@docs
label_state(::AbstractODEProblemParSetter, u0)
name_state(::AbstractODEProblemParSetter, state::AbstractVector)
```

# ODEProblemParSetter

```@docs
ODEProblemParSetter
```
# ODEProblemParSetterConcrete

To support [Concrete ProblemUpdater](@ref) the following variant of 
[`ODEProblemParSetter`](@ref) is provided.
```@docs
ODEProblemParSetterConcrete
```



