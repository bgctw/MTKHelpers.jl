

"""
Support translation between parameter vectors and 
`AbstractODEProblem`.

Takes care of mapping optimized parameters to subset of
- state `u0`, and
- parameters `p`

In addition to the [`AbstractProblemParSetter`](@ref) functions,
which accessing and labelling parameter vectors,
it provides functions that access and label problem state and 
parameters:
[`axis_state`](@ref), [`count_state`](@ref), 
[`keys_state`](@ref), [`symbols_state`](@ref), [`label_state`](@ref), 
[`name_state`](@ref)
"""
abstract type AbstractODEProblemParSetter <: AbstractProblemParSetter end
struct DummyProblemParSetter <: AbstractODEProblemParSetter end # for testing error message


"""
    axis_state(pset::AbstractODEProblemParSetter)
    axis_par(pset::AbstractODEProblemParSetter)

Report the Axis, i.e. nested component symbols of problem states, and problem parameters
respectively.
Returns an AbstractAxis.
"""
function axis_state(::AbstractODEProblemParSetter) end,
function axis_par(::AbstractODEProblemParSetter) end
# need to implement in concrete types

#abstract type AbstractVectorCreator end

"""
    count_state(::AbstractODEProblemParSetter) 
    count_par(::AbstractODEProblemParSetter) 

Report the number of problem states, problem parameters and optimized parameters
respectively.    
"""
function count_state(pset::AbstractODEProblemParSetter) 
    axis_length(axis_state(pset))::Int
end,
function count_par(pset::AbstractODEProblemParSetter) 
    axis_length(axis_par(pset))::Int
end
# TODO change to length(ax) when this becomes available in ComponentArrays
# moved to AbstractProblemParSetter
# function count_paropt(pset::AbstractODEProblemParSetter) 
#     axis_length(axis_paropt(pset))
# end

@deprecate count_states(pset::AbstractProblemParSetter) count_state(pset)

"""
    keys_state(::AbstractODEProblemParSetter) 
    keys_par(::AbstractODEProblemParSetter) 

Report the keys problem states, problem parameters. 
This usually correspdonds to `keys(axis)`, but if there is a
classification, similar to [`keys_paropt`](@ref), report the keys below
this classification.
"""
function keys_state(ps::AbstractODEProblemParSetter) 
    keys(axis_state(ps))
end,
function keys_par(ps::AbstractODEProblemParSetter) 
    keys(axis_par(ps))
end

"""
    symbols_state(pset::AbstractODEProblemParSetter)
    symbols_par(pset::AbstractODEProblemParSetter)
    symbols_paropt(pset::AbstractODEProblemParSetter)

Report the names, i.e. symbols of problem states, problem parameters
respectively, i.e. the concatenation of components.
Similar to `ComponentArrays.label`, but inferred from Axis object.   
"""
function symbols_state(pset::AbstractODEProblemParSetter) 
    _ax_symbols_tuple(axis_state(pset))
end,
function symbols_par(pset::AbstractODEProblemParSetter) 
    _ax_symbols_tuple(axis_par(pset))
end

@deprecate statesyms(pset::AbstractODEProblemParSetter) symbols_state(pset)
@deprecate parsyms(pset::AbstractODEProblemParSetter) symbols_par(pset)


# dispatch to get_paropt_labeled(pset, u0, p)
function get_paropt(pset::AbstractODEProblemParSetter, prob::AbstractODEProblem; kwargs...)
    get_paropt(pset, prob.u0, prob.p; kwargs...)
end,
function get_paropt_labeled(pset::AbstractODEProblemParSetter, prob::AbstractODEProblem; kwargs...)
    get_paropt_labeled(pset, prob.u0, prob.p; kwargs...)
end,
function get_paropt(pset::AbstractODEProblemParSetter, u0, p)
    # depite pset may not be fully inferred, we can determine return type
    T = promote_type(eltype(u0), eltype(p))
    getdata(get_paropt_labeled(pset, u0, p))::Vector{T}
end
# need to implement in concrete types: get_paropt_labeled -> ComponentVector


"""
    label_state(pset::AbstractODEProblemParSetter, u::AbstractVector) 
    label_par(pset::AbstractODEProblemParSetter, par::AbstractVector) 

Produce a labeled version, i.e. a ComponentVector of initial states, parameters
 respectively.
"""
function label_state(pset::AbstractODEProblemParSetter, u)
    attach_axis(u, axis_state(pset))
end,
function label_par(pset::AbstractODEProblemParSetter, p)
    attach_axis(p, axis_par(pset))
end
# label_state(pset::AbstractODEProblemParSetter, u::SVector) = SLVector(label_state(pset, Tuple(u)))
# label_state(pset::AbstractODEProblemParSetter, u::NTuple) = NamedTuple{symbols_state(pset)}(u)


"""
    name_state(pset, u::AbstractVector) 
    name_par(pset, par::AbstractVector) 

Produce a `NamedVector` of given state, parameters, or optimized vars.
"""
function name_state(pset::AbstractODEProblemParSetter, state::AbstractVector)
    NamedArrays.NamedArray(state, (collect(symbols_state(pset))::Vector{Symbol},))
end,
function name_par(pset::AbstractODEProblemParSetter, par::AbstractVector)
    NamedArrays.NamedArray(par, (collect(symbols_par(pset))::Vector{Symbol},))
end

"""
    prob_new = remake(prob::AbstractODEProblem, popt, ps::AbstractODEProblemParSetter)
    u0new, pnew = update_statepar(ps::AbstractODEProblemParSetter, popt, u0, p) 
    deprecated: prob_new = update_statepar(ps::AbstractODEProblemParSetter, popt, prob::AbstractODEProblem) 

Return an updated problem by updating states and parameters where
values corresponding to positions in `popt` hold the values of popt.
The type is changed to the promotion type of popt, to allow working with dual numbers.
"""
function remake_pset(prob::AbstractODEProblem, popt, pset::AbstractODEProblemParSetter)
    u0, p = update_statepar(pset, popt, prob.u0, prob.p)
    remake(prob; u0, p)
end
# need to implement update_statepar in concrete subtypes
@deprecate(update_statepar(pset::AbstractODEProblemParSetter, popt, prob::AbstractODEProblem),
    remake(prob, popt, pset))


"""
    get_u_map(u_new, pset::AbstractODEProblemParSetter)
    get_p_map(p_new, pset::AbstractODEProblemParSetter)

Map each state and parameter the `AbstractODEProblemParSetter` `pset` to a position in names.

When construction an AbstractODEProblem from a ODESystem, the order of states and 
parameters may have changed compared with a previous construction.

In order to set entire state or parameter vectors, a mapping from current
to previous positions, i.e. integer indices, is required, 
so that one can get a vectors in the new format by 
- `u0_old[u_map] .= u_new`
- `p_old[p_map] .= p_new`

`u_new` can be anything for which an axis can be extracted, whose keys are used.
Specifically it can be the ComponentVector of new states itself, or a vector of symbols.

## Keyword arguments
- `is_warn_missing`: set to true to issue warnings if some ODESystem state or 
  parameter names are not found in the old names. This may give false warnings
  for System parameters that have defaults and do not need to be part
  of the parameter vector.
"""
function get_u_map(u_new, pset::AbstractODEProblemParSetter; is_warn_missing = false)
    ax_new = MTKHelpers._get_axis(u_new)
    pos_old = attach_axis(1:count_state(pset), axis_state(pset))
    if is_warn_missing
        missing_keys = setdiff(keys_state(pset), keys(ax_new))
        !isempty(missing_keys) && 
        @warn("problem states $missing_keys are missing in u_new.")
    end
    getdata(pos_old[keys(ax_new)]) 
end,
function get_p_map(p_new, pset::AbstractODEProblemParSetter; is_warn_missing = false)
    ax_new = MTKHelpers._get_axis(p_new)
    pos_old = attach_axis(1:count_par(pset), axis_par(pset))
    if is_warn_missing
        missing_keys = setdiff(keys_par(pset), keys(ax_new))
        !isempty(missing_keys) && 
        @warn("problem parameters $missing_keys are missing in u_new.")
    end
    getdata(pos_old[keys(ax_new)]) 
end

"""
    get_system(prob::AbstractODEProblem) 

Get the System associated to a problem.
"""
function get_system(prob::AbstractODEProblem)
    :sys ∈ propertynames(prob.f) || error(
        "Cannot get System for a Problem without associated system")
   prob.f.sys
end


