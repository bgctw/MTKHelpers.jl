"""
    embed_system(m;name)

Embeds system `m` as the single component of a larger system.
This helps to match the naming of the states, parameters, observables to 
the namespace of the system. 

```jldocstest; output=false
using ModelingToolkit, DifferentialEquations
using MTKHelpers
function samplesystem(;name) 
    @variables t 
    D = Differential(t) 
    sts = @variables x(t) RHS(t)  # RHS is observed
    ps = @parameters τ       # parameters
    ODESystem([ RHS  ~ (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)
end                     
@named m = samplesystem()
@named sys = embed_system(m)

# note that keys are `m.x`,`m.τ` or `m.RHS`.
# Hence only m needs to be defined rather than all the states, parameters, 
# and observables of the system.
prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0), [m.τ => 3.0])
sol = solve(prob);
dx1 = sol[m.RHS][1:3] 
dx2 = sol[getproperty(m,:RHS)][1:3]  # access by symbols

# using Plots
# plot(sol, vars=[m.x,m.RHS])    
# plot(sol, vars=getproperty.(Ref(m),[:x, :RHS])) # access by symbols
length(dx2) == 3
# output
true
```
"""
function embed_system(m;name, simplify=true)
    @named _sys_embed = ODESystem(Equation[], ModelingToolkit.get_iv(m))               
    sys = compose(_sys_embed, [m]; name)
    if simplify
        sys = structural_simplify(sys)
    end
    sys
end

"""
    symbol(t)

Extract the inner symbol from a Term, Num, or BasicSymbolic object.
"""
function symbol(t::Term); symbol(t.f); end
symbol(s::SymbolicUtils.BasicSymbolic) = istree(s) ? Symbol(operation(s)) : Symbol(s)
symbol(num::Num) = symbol(num.val)
symbol(s) = Symbol(s)

"""
    strip_namespace(s)
    
Omit the part before the first dot.
"""    
function strip_namespace(s::AbstractString); match(r"[^.₊]+$",s).match; end
function strip_namespace(s::Symbol); Symbol(strip_namespace(string(s))); end

"""
    symbols_state(sys::ODESystem)
    symbols_par(sys::ODESystem)

Extract the basic symbols without namespace of system states and system parameters.
"""
function symbols_state(sys::ODESystem); symbol.(states(sys)); end
symbols_par(sys::ODESystem) = symbol.(parameters(sys))

# "apply fun to x until fun(x) == x"
# function fixpoint(fun, x, nrecur_max=12; fmap=identity)
#     nrecur_max == 0 && error("cound not find fixpoint for $f and $x.")
#     px = fun(x)
#     #@show nrecur_max, x, fmap(x), fmap(px)
#     fmap(px) == fmap(x) && return(px)
#     fixpoint(fun, px, nrecur_max-1)
# end

using ModelingToolkit: AbstractSystem, get_eqs, get_states, get_ps, get_observed, get_continuous_events, get_defaults, get_systems

"""
    override_system(eqs, basesys::AbstractSystem; 
        name::Symbol=Symbol(string(nameof(basesys))*"_ext"), 
        ps=Term[], 
        obs=Equation[], 
        evs=ModelingToolkit.SymbolicContinuousCallback[], 
        defs=Dict()
    )

Modify `basesys` by replacing some equations matched by their left-hand-side.
The keyword argument correspond to ODESystem.
"""
function override_system(eqs, basesys::AbstractSystem; 
    name::Symbol=Symbol(string(nameof(basesys))*"_ext"), 
    ps=Term[], 
    obs=Equation[], 
    evs=ModelingToolkit.SymbolicContinuousCallback[], 
    defs=Dict()
    )
    T = SciMLBase.parameterless_type(basesys)
    ivs = independent_variables(basesys)
    length(ivs) > 1 && error("Extending multivariate systems is not supported")
    eqs_base_dict = Dict(eq.lhs => eq for eq in get_eqs(basesys))
    eqs_new_keys = [eq.lhs for eq in eqs]
    is_key_present = eqs_new_keys .∈ Ref(keys(eqs_base_dict))
    !all(is_key_present) && error(
        "Expected all lhs of new equations to be present in basesys. " *
        "But following keys were not present: $(string.(eqs_new_keys[.!is_key_present]))")
    eqs_base_keys = setdiff(keys(eqs_base_dict), eqs_new_keys)
    eqs_base_no_overwrite = get.(Ref(eqs_base_dict), eqs_base_keys, missing)
    eqs_ext = union(eqs_base_no_overwrite, eqs)
    sts = get_states(basesys)
    ps_ext = union(get_ps(basesys), ps)
    obs_ext = union(get_observed(basesys), obs)
    evs_ext = union(get_continuous_events(basesys), evs)
    defs_ext = merge(get_defaults(basesys), defs) # prefer new defs 
    syss = get_systems(basesys)
    #
    if length(ivs) == 0
        T(eqs_ext, sts, ps_ext, observed = obs_ext, defaults = defs_ext, name=name, systems = syss, continuous_events=evs_ext)
    elseif length(ivs) == 1
        T(eqs_ext, ivs[1], sts, ps_ext, observed = obs_ext, defaults = defs_ext, name = name, systems = syss, continuous_events=evs_ext)
    end
end
