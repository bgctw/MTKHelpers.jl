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

Extract the inner symbol from a Term or Num object.
"""
function symbol(t::Term); symbol(t.f); end
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
(TYPEDSIGNATURES)

entend the `basesys` with `sys`, the resulting system would inherit `sys`'s name
by default. Contrary to extend, it allows replacing equations.

Experimental: subject to change
"""
function extend_overide(sys::AbstractSystem, basesys::AbstractSystem; name::Symbol=nameof(sys))
    # https://github.com/SciML/ModelingToolkit.jl/issues/1518
    T = SciMLBase.parameterless_type(basesys)
    ivs = independent_variables(basesys)
    if !(sys isa T)
        if length(ivs) == 0
            sys = convert_system(T, sys)
        elseif length(ivs) == 1
            sys = convert_system(T, sys, ivs[1])
        else
            throw("Extending multivariate systems is not supported")
        end
    end

    eqs_base_dict = Dict(eq.lhs => eq for eq in get_eqs(basesys))
    eqs_new_keys = [eq.lhs for eq in get_eqs(sys)]
    eqs_base_keys = setdiff(keys(eqs_base_dict), eqs_new_keys)
    eqs_base_no_overwrite = get.(Ref(eqs_base_dict), eqs_base_keys, missing)

    eqs = union(eqs_base_no_overwrite, get_eqs(sys))
    sts = union(get_states(basesys), get_states(sys))
    ps = union(get_ps(basesys), get_ps(sys))
    obs = union(get_observed(basesys), get_observed(sys))
    evs = union(get_continuous_events(basesys), get_continuous_events(sys))
    defs = merge(get_defaults(basesys), get_defaults(sys)) # prefer `sys`
    syss = union(get_systems(basesys), get_systems(sys))

    if length(ivs) == 0
        T(eqs, sts, ps, observed = obs, defaults = defs, name=name, systems = syss, continuous_events=evs)
    elseif length(ivs) == 1
        T(eqs, ivs[1], sts, ps, observed = obs, defaults = defs, name = name, systems = syss, continuous_events=evs)
    end
end



