# type alias to save typing
const VN = AbstractVector{<:SymbolicUtils.BasicSymbolic}

"""
    ODEProblemParSetterConcrete

Helps keeping track of a subset of initial states and parameters to be optimized.
Similar to [`ODEProblemParSetter`](@ref), but with axis and length information
as type parameters.

It is constructed by `get_concrete(ODEProblemParSetter(...))`.
`
"""
struct ODEProblemParSetterConcrete{NS, NP, NPS, POPTA <: AbstractAxis,
    SA <: AbstractAxis,
    SSA <: AbstractAxis,
    PA <: AbstractAxis,
    POPTAS <: AbstractAxis,
    POPTAF <: AbstractAxis,
    ES <:SymbolicUtils.Symbolic,
    NSO, EIS, 
    NPO, EIP
} <: AbstractODEProblemParSetter
    ax_paropt::POPTA
    ax_state::SA
    ax_state_scalar::SSA
    ax_par::PA
    is_updated_state_i::SVector{NS, Bool}
    is_updated_par_i::SVector{NP, Bool}
    ax_paropt_scalar::POPTAS
    ax_paropt_flat1::POPTAF
    par_nums::SVector{NPS, ES} # not isbits because BasicSymbolic{Real} is not isbits
    opt_state_nums::SVector{NSO, ES}
    opt_par_nums::SVector{NPO, ES}
    par_ind::SVector{NPS, EIP} # propb.ps -> vector
    stateopt_ind::SVector{NSO, EIS}
    popt_ind::SVector{NPO, EIP}
    function ODEProblemParSetterConcrete(
            ax_state::AbstractAxis, ax_state_scalar::AbstractAxis,
            ax_par::AbstractAxis, ax_paropt::AbstractAxis, ax_paropt_scalar::AbstractAxis,
            ax_paropt_flat1::AbstractAxis,
            par_nums::VN, opt_state_nums::VN, opt_par_nums::VN,
            par_ind, stateopt_ind::AbstractVector, popt_ind::AbstractVector
    )
        keys_paropt_state = keys(CA.indexmap(ax_paropt_scalar)[:state])
        keys_paropt_par = keys(CA.indexmap(ax_paropt)[:par])
        is_updated_state_i = isempty(keys_paropt_state) ?
                             SVector{0, Bool}() :
                             SVector((k ∈ keys_paropt_state for k in keys(ax_state))...)
        is_updated_par_i = isempty(keys_paropt_par) ?
                           SVector{0, Bool}() :
                           SVector((k ∈ keys_paropt_par for k in keys(ax_par))...)
        #
        nstate = length(is_updated_state_i) 
        npar = length(is_updated_par_i) # entries may be vectors
        npar_scalar = axis_length(ax_par)
        s_par_nums = CA.SVector{npar_scalar, eltype(par_nums)}(par_nums)
        s_opt_state_nums = CA.SVector{
            axis_length(CA.indexmap(ax_paropt)[:state]), eltype(opt_state_nums)}(opt_state_nums)
        s_opt_par_nums = CA.SVector{
            axis_length(CA.indexmap(ax_paropt)[:par]), eltype(opt_par_nums)}(opt_par_nums)
        s_par_ind = CA.SVector{npar_scalar, eltype(popt_ind)}(par_ind)
        s_stateopt_ind = CA.SVector{
            axis_length(CA.indexmap(ax_paropt)[:state]), eltype(stateopt_ind)}(stateopt_ind)
        s_popt_ind = CA.SVector{
            axis_length(CA.indexmap(ax_paropt)[:par]), eltype(popt_ind)}(popt_ind)
        #Main.@infiltrate_main
        new{nstate, npar, npar_scalar,
            typeof(ax_paropt), typeof(ax_state), typeof(ax_state_scalar),
            typeof(ax_par), typeof(ax_paropt_scalar),
            typeof(ax_paropt_flat1), eltype(par_nums),
            length(s_stateopt_ind), eltype(s_stateopt_ind),
            length(s_popt_ind), eltype(s_popt_ind)
        }(ax_paropt,
            ax_state, ax_state_scalar, ax_par, is_updated_state_i,
            is_updated_par_i, ax_paropt_scalar, ax_paropt_flat1,
            s_par_nums, s_opt_state_nums, s_opt_par_nums,
            s_par_ind, s_stateopt_ind, s_popt_ind
        )
    end
end

isconcrete(::ODEProblemParSetterConcrete) = true

# function ODEProblemParSetterConcrete(state_template, par_template, popt_template,
#         system::Union{AbstractODESystem, Nothing} = nothing;
#         is_validating = Val{true}())
#     ax_paropt = _get_axis(popt_template)
#     ax_state = _get_axis(state_template)
#     ax_par = _get_axis(par_template)
#     ax_state_array = isnothing(system) ? ax_state : axis_of_nums(unknowns(system))
#     if !(:state ∈ keys(ax_paropt) || :par ∈ keys(ax_paropt))
#         ax_paropt = assign_state_par(ax_state_array, ax_par, ax_paropt)
#     end
#     (ax_state_scalar, ax_paropt_scalar) = isnothing(system) ?
#                                           (ax_state, ax_paropt) :
#                                           scalarize_par_and_paroptstate(ax_state,
#         ax_paropt, system)
#     ODEProblemParSetterConcrete(ax_state_scalar,
#         ax_par,
#         ax_paropt,
#         ax_paropt_scalar,
#         is_validating)
# end

# function ODEProblemParSetterConcrete(state_template,
#         par_template, popt_template::Union{NTuple{N, Symbol}, AbstractVector{Symbol}},
#         system::Union{AbstractODESystem, Nothing} = nothing;
#         is_validating = Val{true}()) where {N}
#     ax_par = _get_axis(par_template)
#     #ax_state = _get_axis(state_template)
#     ax_state_array = isnothing(system) ? _get_axis(state_template) : axis_of_nums(unknowns(system))
#     #construct a template by extracting the components of u0 and p
#     u0 = attach_axis(1:axis_length(ax_state_array), ax_state_array)
#     p = attach_axis(1:axis_length(ax_par), ax_par)
#     popt_template_new = vcat(u0, p)[popt_template]
#     #popt_template_new = _get_axis(popt_template) # not the correct length of arrays
#     get_concrete(ODEProblemParSetter(state_template, ax_par, popt_template_new, system; is_validating))
# end

# function get_concrete(ODEProblemParSetter(sys::ODESystem, paropt))
#     # simplify X(t) to X but keep (Y(t))[i] intact
#     get_concrete(ODEProblemParSetter(Axis(symbol_op_scalar.(unknowns(sys)))),
#         axis_of_nums(parameters(sys)),
#         paropt,
#         sys)
# end
