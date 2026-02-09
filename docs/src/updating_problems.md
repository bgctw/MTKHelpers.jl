```@meta
CurrentModule = MTKHelpers
```

# Setting parts of initial state and Parameters of a problem

Often one wants to change a subset of the initial
states,`u0`, and a subset of parameters,`p`, of an AbstractODEProblem during an optimization.

Given `u0` and `p` can be expressed as ComponentVectors, 
then `popt` can be expressed as a ComponentVector of optimized parameters, 
which may include initial states. The initial states and parameter components 
must have different names.

```@example doc_depr
#The following example system employs a scalar and a vector-valued parameter.
# using ModelingToolkit, OrdinaryDiffEq, ComponentArrays
# using MTKHelpers
# using ModelingToolkit: t_nounits as t, D_nounits as D
# function samplesystem(;name,τ = 3.0, p=[1.1, 1.2]) 
#     sts = @variables x(t) RHS(t)        # RHS is observed
#     ps = @parameters τ=τ p[1:2] = p 
#     System([ RHS  ~ p[1] + -p[2]*x + (1 - x)/τ, D(x) ~ RHS ], t, sts, vcat(ps...); name)
# end                     
# @named m = samplesystem()
# @named sys = embed_system(m)
# prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0))
nothing # hide
```

```@example doc
#The following example system employs a scalar-valued parameter.
using ModelingToolkit, OrdinaryDiffEq, ComponentArrays
using ModelingToolkit: ModelingToolkit as MTK
using MTKHelpers
using ModelingToolkit: t_nounits as t, D_nounits as D
function samplesystem(;name,τ = 3.0, p1=1.1, p2=1.2) 
    sts = @variables x(t) RHS(t)        # RHS is observed
    ps = @parameters τ=τ p1 = p1 p2 = p2
    System([ RHS  ~ p1 + -p2*x + (1 - x)/τ, D(x) ~ RHS ], t, sts, vcat(ps...); name)
end                     
@named m = samplesystem()
@named sys = embed_system(m)
prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0))
nothing # hide
```


An [`ODEProblemParSetter`](@ref) then can be used to update a subset of states
and parameters in the derived problem.
Because the state of the problem can reorder components of a symbolic array
the parameter object of the problem is complex.

```@example doc
# setup position matching, note τ is not in parameters optimized
popt = ComponentVector(state=(m₊x=0.1,), par=(m₊p1=2.1,m₊p2=2.2)) 
pset = ODEProblemParSetter(sys, popt) 

# extract optimized state and parameters
get_paropt(pset, prob)          # plain vector
get_paropt_labeled(pset, prob)  # ComponentVector
name_paropt(pset, prob)         # NamedVector 

# update states and parameters
prob2 = remake(prob, popt, pset)

# check that parameters have been updated
get_paropt_labeled(pset, prob2) == popt
# alternatively MTK access
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
prob2.u0[SII.variable_index(sys, m.x)] == popt.state.m₊x
prob2.ps[m.p1] == popt.par.m₊p1
prob2.ps[:m₊p2] == popt.par.m₊p2

# Fruther, check that other parameters did not change
prob2.ps[m.τ] == initial_conditions(get_system(prob2))[m.τ]
nothing # hide
```

`MTKHelper` offers some convenience to acces paramters, states, and optimized parts
as a ComponentVector. `pset` stores the MTK indices so that they do not need to be recreated.

```@example doc
get_state_labeled(pset, prob2).m₊x == popt.state.m₊x 
get_par_labeled(pset, prob2).m₊p2 == popt.par.m₊p2 
get_paropt_labeled(pset, prob2) == popt
nothing # hide
```
There are three  [suggested ways since MTK10](https://docs.sciml.ai/ModelingToolkit/stable/basics/FAQ/#Why-are-my-parameters-some-obscure-object?) is to
- using an index object from [SymbolicIndexingInterface.jl](https://docs.sciml.ai/SymbolicIndexingInterface/stable/)
- using a setter object 
- using a SciMLStructures.jl to replace all tunable parameters
- using the Dictionary interface

Currently only the 4th variant works without problems with AD systems, but it is
not efficient.
The third variant requires changing the system definition to deterine which parameters are optimized.
The first two need quite complex [integration with PreallocationTools](https://docs.sciml.ai/ModelingToolkit/stable/basics/FAQ/#Using-ModelingToolkit-with-Optimization-/-Automatic-Differentiation) 
to be used with AD-systems.

```@example doc
#prob2 = remake(prob, [m.p2 => 3.2]) # would be nice, but not supported
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
ip2 = SII.parameter_index(sys, MTK.parse_variable(sys, "m₊p2"))
setindex!(prob.ps, 3.2, ip2)
prob.ps[m.p2] == 3.2
nothing # hide
```

```@example doc
probc = remake(prob)
setter! = SII.setp(sys, [m.p2])
setter!(probc, popt.par.m₊p2)
probc.ps[m.p2] == popt.par.m₊p2
nothing # hide
```

This package currently in the backround relies on the integration of setter objects with `PreallocationTools.jl`
for supporting `ForwardDiff` and falls back on the Dictionary approach for other AD systems.

```@example doc
using ForwardDiff
gr = ForwardDiff.gradient(
    popt -> remake(prob, popt, pset).ps[m.p1], 
    getdata(popt))
gr == [0.0, 1.0, 0.0]
nothing # hide
```
