"""
    ODEProblemParSetterConcrete(state_template,par_template,popt_template) 
    ODEProblemParSetterConcrete(sys::ODESystem, popt_template) 

Helps keeping track of a subset of initial states and parameters to be optimized.
Similar to [`ODEProblemParSetter`](@ref), but with axis and length information
as type parameters.
`
"""
struct ODEProblemParSetterConcrete{NS, NP, POPTA <: AbstractAxis,
    SA <: AbstractAxis,
    PA <: AbstractAxis,
    POPTAS <: AbstractAxis,
} <: AbstractODEProblemParSetter
    ax_paropt::POPTA
    ax_state::SA
    ax_par::PA
    is_updated_state_i::StaticVector{NS, Bool}
    is_updated_par_i::StaticVector{NP, Bool}
    ax_paropt_scalar::POPTAS
    function ODEProblemParSetterConcrete(ax_state::AbstractAxis,
            ax_par::AbstractAxis, ax_paropt::AbstractAxis, ax_paropt_scalar::AbstractAxis,
            ::Val{isval}) where {isval}
        if isval
            is_valid, msg = validate_keys_state_par(ax_paropt_scalar, ax_state, ax_par)
            !is_valid && error(msg)
        end
        keys_paropt_state = keys(CA.indexmap(ax_paropt_scalar)[:state])
        keys_paropt_par = keys(CA.indexmap(ax_paropt)[:par])
        is_updated_state_i = isempty(keys_paropt_state) ?
                             SVector{0, Bool}() :
                             SVector((k ∈ keys_paropt_state for k in keys(ax_state))...)
        is_updated_par_i = isempty(keys_paropt_par) ?
                           SVector{0, Bool}() :
                           SVector((k ∈ keys_paropt_par for k in keys(ax_par))...)
        new{length(is_updated_state_i), length(is_updated_par_i),
            typeof(ax_paropt), typeof(ax_state), typeof(ax_par), typeof(ax_paropt_scalar)}(ax_paropt,
            ax_state, ax_par, is_updated_state_i, is_updated_par_i, ax_paropt_scalar)
    end
end

function ODEProblemParSetterConcrete(state_template, par_template, popt_template, 
    system::Union{AbstractODESystem,Nothing} = nothing;
        is_validating = Val{true}())
    ax_paropt = _get_axis(popt_template)
    ax_state = _get_axis(state_template)
    ax_par = _get_axis(par_template)
    if !(:state ∈ keys(ax_paropt) || :par ∈ keys(ax_paropt))
        ax_paropt = assign_state_par(ax_state, ax_par, ax_paropt)
    end
    (ax_state_scalar, ax_paropt_scalar) = isnothing(system) ?
                                          (ax_state, ax_paropt) :
                                          scalarize_par_and_paroptstate(ax_state,
        ax_paropt, system)
    ODEProblemParSetterConcrete(ax_state_scalar,
        ax_par,
        ax_paropt,
        ax_paropt_scalar,
        is_validating)
end

function ODEProblemParSetterConcrete(state_template,
        par_template, popt_template::Union{NTuple{N, Symbol}, AbstractVector{Symbol}},
        system::Union{AbstractODESystem,Nothing} = nothing;
        is_validating = Val{true}()) where {N}
    ax_par = _get_axis(par_template)
    ax_state = _get_axis(state_template)
    # construct a template by extracting the components of u0 and p
    u0 = attach_axis(1:axis_length(ax_state), ax_state)
    p = attach_axis(1:axis_length(ax_par), ax_par)
    popt_template_new = vcat(u0, p)[popt_template]
    ODEProblemParSetterConcrete(ax_state, ax_par, popt_template_new, system; is_validating)
end

function ODEProblemParSetterConcrete(sys::ODESystem, paropt)
    # simplify X(t) to X but keep (Y(t))[i] intact
    ODEProblemParSetterConcrete(Axis(symbol_op_scalar.(states(sys))), axis_of_nums(parameters(sys)), paropt, sys)
end
