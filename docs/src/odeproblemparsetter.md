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

## Setting entire state and parameter vectors
```@docs
get_u_map
```

# ODEProblemParSetter

```@docs
ODEProblemParSetter
```



