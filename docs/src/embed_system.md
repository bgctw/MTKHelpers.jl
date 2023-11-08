```@meta
CurrentModule = MTKHelpers
```

# Representing a system as a component

```@docs
embed_system
```

# Overriding equations of an existing system

For debugging bigger systems, it is helful to set some equations
to zero or modify/simplify the system in other ways.

Function `override_system` takes a set of equations and matches the
right-hand site to the equations of the original system and replaces
those equations.

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

# simplify the system by setting RHS ~ RHS_0 * x
function samplesystem_const(RHS0; name) 
    # demonstrating override_system by setting the RHS to constant first order rate
    m = samplesystem(;name)
    @unpack RHS,x = m
    @parameters t 
    ps = @parameters RHS_0=RHS0
    D = Differential(t)
    eqs = [
        RHS ~ RHS_0 * x,
    ]
    sys_ext = override_system(eqs, m; name, ps) 
end  

@named mc = samplesystem_const(-0.1)
@named sys = embed_system(mc)
prob = ODEProblem(sys, [mc.x => 1.0], (0.0,10.0))
sol = solve(prob)
isapprox(sol(8)[1], exp(-0.1*8), atol = 1e-5)
```

```@docs
override_system
```

# Utilities
```@docs
symbol
strip_namespace
symbols_state(::ODESystem)
```


