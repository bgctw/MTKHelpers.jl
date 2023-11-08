"""
    ProblemParSetter(state_names,par_names,popt_names) 
    ProblemParSetter(sys::ODESystem, popt_names; strip=false) 

Helps keeping track of a subset of initial states and paraemters to be optimized.

# Arguments
- `state_names`: ComponentVector or Axis of all the initial states of the problem
- `par_names`: all the parameters of the problem
- `popt_names`: the parameter/initial states to be optimized.

If all of `state_names`, `par_names`, and `popt_names` are type-inferred Axes,
then also the constructed ProblemParSetter is type-inferred.

The states and parameters can be extracted from an `ModelingToolkit.ODESystem`.
If `strip=true`, then namespaces of parameters of a composed system are removed, 
e.g. `subcomp₊p` becomes `p`.
"""
struct ProblemParSetter{NOPT,POPTA<:AbstractAxis,SA<:AbstractAxis,PA<:AbstractAxis} <:
       AbstractProblemParSetter
    # u0_opt::NTuple{NS, Symbol}
    # p_opt::NTuple{NP, Symbol}
    ax_paropt::POPTA
    ax_state::SA
    ax_par::PA
    is_state::NTuple{NOPT,Bool}
    is_p::NTuple{NOPT,Bool}
end

function ProblemParSetter(
    ax_state::AbstractAxis,
    ax_par::AbstractAxis,
    ax_paropt::AbstractAxis,
)
    popt_names = keys(CA.indexmap(ax_paropt))
    state_names = keys(CA.indexmap(ax_state))
    par_names = keys(CA.indexmap(ax_par))
    NOPT = length(CA.indexmap(ax_paropt))
    is_state = ntuple(i -> popt_names[i] ∈ state_names, NOPT)
    is_p = ntuple(i -> popt_names[i] ∈ par_names, NOPT)
    is_u0_or_p = is_state .| is_p
    all(is_u0_or_p) || @warn(
        "missing optimization parameters in system: " *
        join(popt_names[collect(.!is_u0_or_p)], ", ")
    )
    is_u0_and_p = is_state .& is_p
    any(is_u0_and_p) && @warn(
        "expted parameter names and state names to be distinct, but occurred in both: " *
        join(popt_names[collect(is_u0_and_p)], ", ")
    )
    ProblemParSetter{NOPT,typeof(ax_paropt),typeof(ax_state),typeof(ax_par)}(
        ax_paropt,
        ax_state,
        ax_par,
        is_state,
        is_p,
    )
end

# function ProblemParSetter(state_names::NTuple{NS, Symbol}, par_names::NTuple{NP, Symbol}, popt_names::NTuple{NOPT, Symbol}) where {NS, NP, NOPT} 
#     # not type stable
#     ProblemParSetter(Axis(state_names), Axis(par_names), Axis(popt_names))
# end

function ProblemParSetter(state_template, par_template, popt_template)
    # not type-stable
    ProblemParSetter(
        _get_axis(state_template),
        _get_axis(par_template),
        _get_axis(popt_template),
    )
end

# function _get_axis(x::AbstractArray) 
#     @info("Providing Parameters as Array was deprecated for performance?")
#     # depr?: need a full-fledged axis
#     Axis(Tuple(i for i in symbol.(x)))
# end
function _get_axis(x::Tuple)
    Axis(Tuple(i for i in symbol.(x)))
end
_get_axis(x::ComponentVector) = first(getaxes(x))
_get_axis(x::AbstractAxis) = x


function ProblemParSetter(sys::ODESystem, popt_names; strip = false)
    ft = strip ? strip_namespace : identity
    state_names = ft.(symbol.(states(sys)))
    par_names = ft.(symbol.(parameters(sys)))
    ProblemParSetter(CA.Axis(state_names), CA.Axis(par_names), _get_axis(popt_names))
end

# count_state(::ProblemParSetter{N, POPTA, SA, PA}) where {N, POPTA, SA, PA} = length(CA.indexmap(SA))
# count_par(::ProblemParSetter{N, POPTA, SA, PA}) where {N, POPTA, SA, PA} = length(CA.indexmap(PA))
# count_paropt(::ProblemParSetter{N, POPTA, SA, PA}) where {N, POPTA, SA, PA} = N

# TODO change to length(ax) when this becomes available in ComponentArrays
count_state(pset::ProblemParSetter) = axis_length(pset.ax_state)
count_par(pset::ProblemParSetter) = axis_length(pset.ax_par)
count_paropt(pset::ProblemParSetter) = axis_length(pset.ax_paropt)

axis_length(ax::AbstractAxis) = lastindex(ax) - firstindex(ax) + 1
axis_length(::FlatAxis) = 0


# axis_state(::ProblemParSetter{N, POPTA, SA, PA}) where {N, POPTA, SA, PA} = SA
# axis_par(::ProblemParSetter{N, POPTA, SA, PA}) where {N, POPTA, SA, PA} = PA
# axis_paropt(::ProblemParSetter{N, POPTA, SA, PA}) where {N, POPTA, SA, PA} = POPTA

axis_state(ps::ProblemParSetter) = ps.ax_state
axis_par(ps::ProblemParSetter) = ps.ax_par
axis_paropt(ps::ProblemParSetter) = ps.ax_paropt


# function _ax_symbols(ax::Union{AbstractAxis, CA.CombinedAxis}; prefix="₊") 
#     # strip the first prefix, convert to symbol and return generator
#     (i for i in _ax_string_prefixed(ax; prefix) .|> (x -> x[(sizeof(prefix)+1):end]) .|> Symbol)
# end
# function _ax_symbols_tuple(ax::Union{AbstractAxis, CA.CombinedAxis}; kwargs...) 
#     Tuple(_ax_symbols(ax; kwargs...))::NTuple{CA.last_index(ax), Symbol}
# end
# function _ax_symbols_vector(ax::Union{AbstractAxis, CA.CombinedAxis}; kwargs...) 
#     # strip the first prefix, convert to symbol and collect into tuple
#     collect(_ax_symbols(ax; kwargs...))::Vector{Symbol}
# end

symbols_state(pset::ProblemParSetter) = _ax_symbols_tuple(axis_state(pset))
symbols_par(pset::ProblemParSetter) = _ax_symbols_tuple(axis_par(pset))
symbols_paropt(pset::ProblemParSetter) = _ax_symbols_tuple(axis_paropt(pset))

# # Using unexported interface of ComponentArrays.axis, one place to change
# "Accessor function for index from ComponentIndex"
# idx(ci::CA.ComponentIndex) = ci.idx

"""
    prob_new = update_statepar(pset::ProblemParSetter, popt, prob::ODEProblem) 
    u0new, pnew = update_statepar(pset::ProblemParSetter, popt, u0, p) 

Return an updated problem or updates states and parameters where
values corresponding to positions in `popt` are set.
"""
# function update_statepar(pset::ProblemParSetter, popt, prob::ODEProblem) 
#     u0,p = update_statepar(pset, popt, prob.u0, prob.p)
#     remake(prob; u0, p)
# end


function update_statepar(pset::ProblemParSetter, popt::TO, u0::TU, p::TP) where {TO,TU,TP}
    popt_state, popt_p = _separate_state_p(pset, popt)
    u0new =
        length(popt_state) == 0 ? u0 :
        begin
            u0_l = label_state(pset, u0)
            fdata = (u0 isa ComponentVector) ? identity : getdata # strip labels?
            TUP = typeof(similar(u0_l, promote_type(eltype(TU), eltype(TO))))
            u0new = fdata(_update_cv_top(u0_l, popt_state)::TUP)
        end
    pnew =
        length(popt_p) == 0 ? p :
        begin
            p_l = label_par(pset, p)
            fdata = (p isa ComponentVector) ? identity : getdata # strip labels?
            TPP = typeof(similar(p_l, promote_type(eltype(TP), eltype(TO))))
            pnew = fdata(_update_cv_top(p_l, popt_p)::TPP)
        end
    u0new, pnew
end

# extract p and state components of popt into separate ComponentVectors
function _separate_state_p(pset, popt)
    popt_l = label_paropt(pset, popt)
    syms_p = Tuple(k for (i, k) in enumerate(keys(popt_l)) if pset.is_p[i])
    syms_s = Tuple(k for (i, k) in enumerate(keys(popt_l)) if pset.is_state[i])
    # popt_state = _get_index_axis(popt_l, _get_axis_of_lengths_for_syms(popt_l,
    #     Tuple(k for (i,k) in enumerate(keys(popt_l)) if pset.is_state[i])))
    # #popt_p = popt_l[Axis( # WAIT use proper indexing when supported with ComponentArrays
    # popt_p = _get_index_axis(popt_l, _get_axis_of_lengths_for_syms(popt_l,
    #         Tuple(k for (i,k) in enumerate(keys(popt_l)) if pset.is_p[i])))
    #popt_state, popt_p
    _indexof(popt_l, syms_s), _indexof(popt_l, syms_p)
end

# TODO after cv[()] works with empty and mono-Tuple, replace by indexof(cv, syms)
get_empty(cv::ComponentVector) =
    ComponentVector(similar(getdata(cv), 0), (ComponentArrays.NullAxis(),))
function _indexof(cv::ComponentVector{T}, syms::NTuple{N,Symbol}) where {T,N}
    N == 1 && return (cv[KeepIndex(syms[1])])::ComponentVector{T}
    N == 0 && return (get_empty(cv))::ComponentVector{T}
    res = cv[syms]::ComponentVector{T}
end




# function typed_from_generator(type::Type, v0, vgen) where T 
#     if type <: AbstractArray
#     end
#     typeof(v0)(vgen)
# end
# # AbstractVector{T} is not more specific than ::Type, need to support all used concrete
# #typed_from_generator(::Type{AbstractVector{T}}, v0, vgen) where T = convert(typeof(v0), v0 .+ vgen)::typeof(v0)
# typed_from_generator(::Type{Vector{T}}, v0, vgen) where T = convert(typeof(v0), v0 .+ vgen)::typeof(v0)


# typed_from_generator(v0, vgen) = typeof(v0)(vgen)
# function typed_from_generator(v0::AbstractVector, vgen) 
#     # convert to std vector, because typeof(v0) does not contain all entries
#     T = Vector{eltype(v0)}  
#     convert(T, v0 .+ vgen)::T
# end

# type piracy - try to get into CompponentArrays
# Base.getindex(cv::ComponentVector, ax::AbstractAxis) = _get_index_axis(cv,ax)
# Base.getindex(cv::ComponentVector, cv_template::ComponentVector) = _get_index_axis(
#     cv,first(getaxes(cv_template)))


# function get_paropt_labeled(pset::ProblemParSetter, u0, p) 
#     ax = axis_paropt(pset)
#     keys_ax = keys(ax)
#     u0l = label_state(pset, u0)
#     pl = label_par(pset, p)
#     (i,k) = first(enumerate(keys_ax))
#     #(i,k) = (2, keys_ax[2])
#     fik = (i,k) -> begin
#         cvs = pset.is_state[i] ? getproperty(u0l,k) : (pset.is_p[i] ? getproperty(pl,k) : missing) 
#         ismissing(cvs) && return(
#             ComponentVector(NamedTuple{(k,)}((fill(missing,length(ax[k].idx)),))))
#         !(cvs isa ComponentVector) && return(
#             cvs = pset.is_state[i] ? u0l[KeepIndex(k)] : pl[KeepIndex(k)])
#         axs = ax[k].ax
#         #@show cvs, axs, typeof(cvs)
#         # TODO replace by cvs[axs] when _get_index was merged
#         _get_index_axis(cvs, axs)
#     end
#     tmp = (fik(i,k) for (i,k) in enumerate(keys_ax))
#     T = promote_type(eltype(u0), eltype(p))
#     res = reduce(vcat, tmp)::ComponentVector{T, Vector{T}}
#     # @infiltrate
#     label_paropt(pset, res) # reattach axis for type inference
# end

function get_paropt_labeled(pset::ProblemParSetter, u0, p)
    let u0_l = label_state(pset, u0), p_l = label_par(pset, p)
        fik =
            (i, k) ->
                pset.is_state[i] ? u0_l[KeepIndex(k)] :
                (pset.is_p[i] ? p_l[KeepIndex(k)] : missing)
        tmp = (fik(ik[1], ik[2]) for ik in enumerate(keys(axis_paropt(pset))))
        T = promote_type(eltype(u0), eltype(p))
        res = (
            isempty(tmp) ? ComponentVector{T}() : reduce(vcat, tmp)
        )::ComponentVector{T,Vector{T}}
        label_paropt(pset, res) # reattach axis for type inference
    end
end


# attach type in
# TODO type piracy I - until get this into ComponentArrays
# @inline CA.getdata(x::ComponentArray{T,N,A}) where {T,N,A} = getfield(x, :data)::A
# @inline CA.getdata(x::ComponentVector{T,A}) where {T,A} = getfield(x, :data)::A
#@inline CA.getdata(x::ComponentArray) = getfield(x, :data)
#@inline CA.getdata(x::ComponentVector) = getfield(x, :data)


# # extends Base.merge to work on SVector
# # ?type piracy
# merge(x::T, y::NamedTuple) where T<:SLArray = T(merge(NamedTuple(x),y)...)

# # and on Labelled Arrays
# function merge(x::T, y::NamedTuple) where T<:LArray
#     xnew = deepcopy(x)
#     for (key,val) in pairs(y)
#         xnew[key] = val
#     end
#     xnew
# end

"""
    get_u_map(names_u, pset::AbstractProblemParSetter)
    get_p_map(names_p, pset::AbstractProblemParSetter)

Map each state and parameter the `ProblemParSetter` `pset` to a position in names.

When construction an ODEProblem from a ODESystem, the order of states and 
parameters may have changed compared with a previous construction.

In order to set entire state or parameter vectors, a mapping from current
to previous positions, i.e. integer indices, is required, 
so that one can get a vectors in the new format by 
- `u0_old[u_map]`
- `p_old[p_map]`

The mapping is constructed by supplying the names of u0_old and p_old to
a ProblemParameterSetter constructed with the current ODESystem.

## Keyword arguments
- `do_warn_missing`: set to true to issue warnings if some ODESystem state or 
  parameter names are not found in the old names. This may give false warnings
  for System parameters that have defaults and do not need to be part
  of the parameter vector.
"""
function get_u_map(names_u, pset::AbstractProblemParSetter; do_warn_missing = false)
    names_uprob = symbols_state(pset)
    u_map = map(name_uprob -> findfirst(isequal(name_uprob), names_u), names_uprob)
    do_warn_missing &&
        any(isnothing.(u_map)) &&
        warning("problem states $(names_pprob[findall(isnothing.(u_map))]) not in names_u.")
    SVector(u_map)
end,
function get_p_map(names_p, pset::AbstractProblemParSetter; do_warn_missing = false)
    names_pprob = symbols_par(pset)
    p_map = map(name_pprob -> findfirst(isequal(name_pprob), names_p), names_pprob)
    # usually the default parameters, such as u_PlantPmax ~ i_L0 / β_Pi0 - imbalance_P
    # are not part of names_p -> false warning
    do_warn_missing &&
        any(isnothing.(p_map)) &&
        warning(
            "problem parameters $(names_pprob[findall(isnothing.(p_map))]) not in names_p.",
        )
    SVector(p_map)
end
