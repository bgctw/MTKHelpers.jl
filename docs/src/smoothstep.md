```@meta
CurrentModule = MTKHelpers
```

# Smooth steps

ODE solvers have a hard time with step changes, where the derivative
changes discontinuously.
The following functin help to avoid associated problems by approximating
the step by a smoother function. Argument `dx` controls the 
smoothness.

```@docs
smoothstep
```

