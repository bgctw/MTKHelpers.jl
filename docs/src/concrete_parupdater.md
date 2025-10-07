```@meta
CurrentModule = MTKHelpers
```

# Concrete ProblemUpdater

The default-implementations of `ODEParSetter` and `ParUpdater` do not
store much information in type parameters to avoid long compilation times.

The disadvantage is, that the result of some operations is not fully type-
inferred, especially when dealing with ComponentArrays that store information
in their type signature. The results of the folloginw calls are not inferred
(`x` is stands in for either `state`, `par`, or `paropt`):
- `axis_x(pset)`
- `label_x(pset, ...)` 
- `get_paropt(pset, ...)`, and `get_paropt_labeled(pset, ...)`
- `update_statepar` and `remake`
- `count_x(pset)` is inferred but no known at compile time and 
   can not be used to create StaticArrays

Therefore, this package provides function `get_concrete`, that provides
a concrete-typed version of a type.

```@docs
get_concrete
```

## Example
The following example demonstrates, how to construct a cost function that
is based on a closure in which the types are inferred.

First, lets setup a small example problem
```@example doc
using OrdinaryDiffEq, ComponentArrays, MTKHelpers, Test, ModelingToolkit   
using ModelingToolkit: t_nounits as t, D_nounits as D

function get_sys1()
    sts = @variables L(t)
    ps = @parameters k_L, k_R, k_P
    eq = [D(L) ~ 0]
    sys1 = ODESystem(eq, t, sts, vcat(ps...); name = :sys1)
end
sys1 = mtkcompile(get_sys1())
u0 = ComponentVector(L = 10.0)
p = ComponentVector(k_L = 1.0, k_R = 1 / 20, k_P = 2.0)
prob = ODEProblem(sys1,
    get_system_symbol_dict(sys1, u0), (0.0, 1.0),
    get_system_symbol_dict(sys1, p))
nothing # hide
```

Updating a problem with the default ProblemUpdater results in a non-type inferred
problem.
```@example doc
mapping = (:k_L => :k_R, :k_L => :k_P)
pg = KeysProblemParGetter(mapping, keys(u0)) 
pu = get_ode_problemupdater(pg, get_system(prob))
prob2 = pu(prob) # not inferred
nothing # hide
```

But we can use `get_concrete` to obtain a type-inferred version. 
```@example doc
puc1 = get_concrete(pu)
#prob3 = @inferred puc1(prob) # currently not working because remake not type-stable
prob3 = puc1(prob)
nothing # hide
```

This can be used to create a closure for a cost function that uses
the type-stable variant.

```@example doc
# get a concrete-type version of the ProblemParSetter and pass it 
# through a function barrier to a closure (function within let)
get_fopt = (puc=get_concrete(pu)) -> let puc=puc
    (prob) -> begin
        #prob_upd = @inferred puc(prob) # TODO inferred
        prob_upd = puc(prob)
    end # fopt function
end # let, get_fopt
fopt = get_fopt()
#prob4 = @inferred fopt(prob) # TODO inferred
prob4 = fopt(prob)
nothing # hide
```

