```@meta
CurrentModule = MTKHelpers
```

# Setting Initial state and Parameters of a problem

Often one wants to change a subset of the initial
states,`u0`, and a subset of parameters,`p`, of an AbstractODEProblem during an optimization.

Given `u0` and `p` can be expressed as ComponentVectors, 
and `popt` can be expressed as a ComponentVector of optimized parameters, 
which may include initial states. The initial states and parameter components 
must have different names.

The following example system employs a scalar and a vector-valued parameter.
```@example doc
using ModelingToolkit, OrdinaryDiffEq, ComponentArrays
using MTKHelpers
function samplesystem(;name,τ = 3.0, p=[1.1, 1.2]) 
    @variables t 
    D = Differential(t) 
    sts = @variables x(t) RHS(t)        # RHS is observed
    ps = @parameters τ=τ p[1:2] = p 
    ODESystem([ RHS  ~ p[1] + -p[2]*x + (1 - x)/τ, D(x) ~ RHS ], t, sts, vcat(ps...); name)
end                     
@named m = samplesystem()
@named sys = embed_system(m)
prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0))
nothing # hide
```

An [`ODEProblemParSetterConcrete`](@ref) then can be used to update a subset of states
and parametes in the derived problem.

```@example doc
# setup position matching, note τ is not in parameters optimized
popt = ComponentVector(state=(m₊x=0.1,), par=(m₊p=[2.1,2.2],)) 
pset = ODEProblemParSetterConcrete(sys, popt) # pass through function barrier to use type inference

# extract optimized 
get_paropt(pset, prob)          # plain vector
get_paropt_labeled(pset, prob)  # ComponentVector
name_paropt(pset, prob)         # NamedVector 

# update states and parameters
prob2 = remake(prob, popt, pset)
prob2.p # p is still a plain vector
label_par(pset, prob2.p).m₊p == popt.par.m₊p # attach labels and access properties
label_state(pset, prob2.u0).m₊x == popt.state.m₊x # attach labels and access properties
get_paropt_labeled(pset, prob2) == popt
nothing # hide
```

Note that constructing and ODEProblemParSetterConcrete, `pset`, is only fully type-inferred 
when passing three ComponentArrays.Axis objects. This propagates to all ComponentVectors 
constructed by it, e.g. with `label_state`.
Hence, its recommended to pass `pset` across a function barrier for code
where performance matters.

## ProblemUpdater
A [`ODEProblemParSetterConcrete`](@ref) can be combined with a [`KeysProblemParGetter`](@ref)
or other specific implementations of [`AbstractProblemParGetter`](@ref) to 
update an AbstractODEProblem based on information already present in the AbstractODEProblem.

The following example updates parameters `k_R` and `k_P` in the AbstractODEProblem
to the value of `k_L`. This can be useful to ensure that these parameters
are also changed when optimizing parameter `k_L`.

An implementations of `AbstractProblemParGetter` can use any computation of
the source keys to provide its destination keys. It should implement the keys method,
so that when constructing the ProblemUpdater, consistent keys are used,
as in the example below.

```@example doc
using ModelingToolkit, OrdinaryDiffEq, ComponentArrays
using MTKHelpers
f = (u,p,t) -> p[1]*u
u0 = ComponentVector(L=1/2)
p = ComponentVector(k_L = 1.0, k_R = [2.0,3.0], k_P = [4.0,5.0], k_L2 = 6.0)
tspan = (0.,1.)
prob = ODEProblem(f,getdata(u0),tspan,getdata(p))
#
mapping = (:k_L => :k_L2, :k_R => :k_P)
pu = get_ode_problemupdater(KeysProblemParGetter(mapping, keys(u0)), u0, p)
#axis_par(par_setter(pu))
prob2 = pu(prob)
p2 = label_par(par_setter(pu), prob2.p)
p2.k_P == p.k_R
p2.k_L2 == p.k_L
nothing # hide
```
