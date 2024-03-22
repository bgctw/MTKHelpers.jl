```@meta
CurrentModule = MTKHelpers
```

# Setting Initial state and Parameters of a problem

Often one wants to change a subset of the initial
states,`u0`, and a subset of parameters,`p`, of an AbstractODEProblem during an optimization.

Given `u0` and `p` can be expressed as ComponentVectors, 
then `popt` can be expressed as a ComponentVector of optimized parameters, 
which may include initial states. The initial states and parameter components 
must have different names.

The following example system employs a scalar and a vector-valued parameter.
```@example doc
using ModelingToolkit, OrdinaryDiffEq, ComponentArrays
using MTKHelpers
using ModelingToolkit: t_nounits as t, D_nounits as D
function samplesystem(;name,τ = 3.0, p=[1.1, 1.2]) 
    sts = @variables x(t) RHS(t)        # RHS is observed
    ps = @parameters τ=τ p[1:2] = p 
    ODESystem([ RHS  ~ p[1] + -p[2]*x + (1 - x)/τ, D(x) ~ RHS ], t, sts, vcat(ps...); name)
end                     
@named m = samplesystem()
@named sys = embed_system(m)
prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0))
nothing # hide
```

An [`ODEProblemParSetter`](@ref) then can be used to update a subset of states
and parameters in the derived problem.
Because the state of the problem can reorder components of a symbolic array
and the parameter object of the problem is complex, use functions
`get_par(pset, prob)` and `get_state(pset, prob)` or their labeled versions to 
inspect state and parameters of a problem.

```@example doc
# setup position matching, note τ is not in parameters optimized
popt = ComponentVector(state=(m₊x=0.1,), par=(m₊p=[2.1,2.2],)) 
pset = ODEProblemParSetter(sys, popt) # pass through function barrier to use type inference

# extract optimized 
get_paropt(pset, prob)          # plain vector
get_paropt_labeled(pset, prob)  # ComponentVector
name_paropt(pset, prob)         # NamedVector 

# update states and parameters
prob2 = remake(prob, popt, pset)
get_par_labeled(pset, prob2).m₊p == popt.par.m₊p # attach labels and access properties
get_state_labeled(pset, prob2).m₊x == popt.state.m₊x # attach labels and access properties
get_paropt_labeled(pset, prob2) == popt
nothing # hide
```

Note that produced labeled CompoenentArrays are not fully type-inferred, unless
a concrete versions of the ParameterSetter and function barriers are used as described 
in [Concrete ProblemUpdater](@ref).

## ProblemUpdater
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
    sys1 = ODESystem(eq, t, sts, vcat(ps...); name = :sys1)
end
sys1 = structural_simplify(get_sys1())
u0 = ComponentVector(L = 10.0)
p = ComponentVector(k_L = 1.0, k_R = 1 / 20, k_P = 2.0)
prob = ODEProblem(sys1,
    get_system_symbol_dict(sys1, u0), (0.0, 1.0),
    get_system_symbol_dict(sys1, p))
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
